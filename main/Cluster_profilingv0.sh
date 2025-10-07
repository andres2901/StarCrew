#!/bin/bash

# Function to print help message

function print_help() {
   echo -e "Script to run syntenet pipeline using DIAMOND as sequence similarity search software, and summarize its results.
   This script perform five steps:
   1. Run preprocessing data from Syntenet.
   2. Run diamond on the preprocess data.
   3. Run interspecies synteny of syntenet that identify regions with collinearity genes.
   4. Summarize the results of syntenet on three possible modes:
     a. Raw: Return pairs that have a minimum of 8% of shared collinear genes. WARNING: False positive could be return with this mode.
     b. SSP: Return only pairs with strong synteny.
     c. Filter: Return pairs that were filter with a blast approach at the nucleotide level.
   5. Define Clusters of elements and perform a spectral clustering process to identify possible subclusters.
   
   The working directory should have the following structure:
   WorkingDiretory/
   ├── *_gff/
   │   ├── element01.gff
   │   └── element02.gff
   │   ︙
   ├── *_nucleotide/
   │   ├── element01.fa
   │   └── element02.fa
   │   ︙
   ├── *_protein/
   │   ├── element01.fa
   │   └── element02.fa
   │   ︙
   ├── captainPhylogeny/
   │   ├── element01.fa
   │   └── element02.fa
   │   ︙
   └── metadata_updated.csv"
   echo
   echo "Syntax: $0 [ -h ] -w <working_directory> -p <python_function_path> [ -m <mode> -a <anchor_points> -g <gaps> -n <min_nodes> -s <min_size> -t <threshold> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-p --python: Path to the folder where the scripts 'Blast_CleanUp.py' and 'SimpleCluster.py' are stored (required)."
   echo "-m, --mode: Specified the mode to run the summarizing process of synteny results (Available mode: Raw, SSP, Filter) (Default = Raw)."
   echo "-a, --anchors: Number of minimum anchor points for syntenet to call a collinear reagions (Integer between 5 and 25) (Default = 5)."
   echo "-g, --gaps: Number of maximum allowed gaps between anchor points for syntenet to call a collinear reagions (Integer between 5 and 25) (Default = 5)."
   echo "-s, --searchThreads: Number of threads for searching software (DIAMOND and blasd) (Default: 8)"
   echo "-h, --help: Display this help message."
}

# Initialize variables

workingDirectory_path=""
mode="Raw"
pythonFunction_path=""
anchorPoints="5"
gaps="5"
threads="8"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workingDirectory)
            shift
            workingDirectory_path="$1"
            ;;
        -m|--mode)
            shift
            mode="$1"
            ;;
        -p|--python)
            shift
            pythonFunction_path="$1"
            ;;
        -a|--anchors)
            shift
            anchorPoints="$1"
            ;;
        -g|--gaps)
            shift
            gaps="$1"
            ;;
        -s|--searchThreads)
            shift
            threads="$1"
            ;;
        -h|--help)
            help_flag=true
            ;;
        *)
            echo "Invalid option: $1"
            print_help
            exit 1
            ;;
    esac
    shift
done

# Print help if requested
if $help_flag; then
    print_help
    exit 0
fi

# Check for mandatory argument and define the path as absolute
if [[ -z "$workingDirectory_path" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if working directory exists
if [[ ! -d "$workingDirectory_path" ]]; then
    echo "Error: Directory '$workingDirectory_path' does not exist."
    exit 1
else
    workingDirectory_path=$(realpath $workingDirectory_path)
fi

# Check if mode parameter is correct
if [[ "$mode" != "Raw" && "$mode" != "SSP" && "$mode" != "Filter" ]]; then
    echo "Error: provided mode '$mode' is not accepted."
    print_help
    exit 1
fi

# Check if directory where the python function are stored exists
if [[ ! -d "$pythonFunction_path" ]]; then
    echo "Error: Directory '$pythonFunction_path' does not exist."
    exit 1
else
    pythonFunction_path=$(realpath $pythonFunction_path)
fi

# Check if python 'Blast_CleanUp.py' and 'SimpleCluster.py' scripts exists in the given path 
if [[ "${mode}" == "Filter" ]]; then
    if [[ ! -f "${pythonFunction_path}/Blast_CleanUp.py" ]]; then
        echo "Error: File '${pythonFunction_path}/Blast_CleanUp.py' does not exist."
        exit 1
    fi
fi

if [[ ! -f "${pythonFunction_path}/SimpleCluster.py" ]]; then
    echo "Error: File '${pythonFunction_path}/SimpleCluster.py' does not exist."
    exit 1
fi

# Check anchor points and maximum gaps are within allowed range
if [[ "$anchorPoints" =~ ^[0-9]+$ ]]; then
    if (( anchorPoints < 3 || anchorPoints > 25 )); then
        echo "Error: '$anchorPoints' anchor points is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$anchorPoints' is not a positive integer."
    print_help
    exit 1
fi

if [[ "$gaps" =~ ^[0-9]+$ ]]; then
    if (( gaps < 5 || gaps > 25 )); then
        echo "Error: '$gaps' anchor points is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$gaps' is not a positive integer."
    print_help
    exit 1
fi

# Check thread parameter
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

check_directory_structure() {
    local base_dir="$1"
    
    # Step 1: Locate required subdirectories and file
    local gff_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local nucleotide_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local protein_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    local CaptainPhylogeny_dir=$(find "$base_dir" -maxdepth 1 -type d -name "captainPhylogeny" 2>/dev/null)

    if [[ -z "$gff_dir" ]]; then
        echo "Error: GFF subdirectory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$nucleotide_dir" ]]; then
        echo "Error: Nucleoide subdirectory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$protein_dir" ]]; then
        echo "Error: Protein subdirectory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$CaptainPhylogeny_dir" ]]; then
        echo "Error: Captain phylogeny subdirectory not found in '$base_dir'." >&2
        exit 1
    else
        if [[ ! -f "${CaptainPhylogeny_dir}/Captain_tree.treefile" ]]; then
            echo "Error: '${CaptainPhylogeny_dir}/Captain_tree.treefile' file not found in '$base_dir'." >&2
            exit 1
        fi
    fi

    # Step 2: Check for consistent filenames across subdirectories

    # Define temporary file paths in the working directory
    local gff_files="$base_dir/gff_files.txt"
    local nucleotide_files="$base_dir/nucleotide_files.txt"
    local protein_files="$base_dir/protein_files.txt"
    
    # Get sorted list of base filenames from the GFF directory
    find "$gff_dir" -maxdepth 1 -type f -name "*.gff" | xargs -n 1 basename -s .gff | sort > "$gff_files"
    
    # Get sorted list of base filenames from the nucleotide directory
    find "$nucleotide_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$nucleotide_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$protein_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$protein_files"

    # Compare the lists. If diff finds a difference, it returns a non-zero exit code.
    if ! diff -q "$gff_files" "$nucleotide_files" >/dev/null || \
       ! diff -q "$gff_files" "$protein_files" >/dev/null; then
        echo "Error: File lists in subdirectories do not match." >&2
        echo "Details:" >&2
        echo "GFF vs. Nucleotide:" >&2
        diff "$gff_files" "$nucleotide_files" >&2
        echo "GFF vs. Protein:" >&2
        diff "$gff_files" "$protein_files" >&2
        rm "$gff_files" "$nucleotide_files" "$protein_files"
        exit 1
    fi

    # Cleanup temporary files
    rm "$gff_files" "$nucleotide_files" "$protein_files"
}

process_diamond() {
    local working_dir="$1"

    # Locate the protein directory based on the pattern *_ProcessedProtein
    local protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_ProcessedProtein" 2>/dev/null)
    
    if [[ -z "$protein_path" ]]; then
        echo "Error: Protein directory (*_ProcessedProtein) not found in '$working_dir'." >&2
        exit 1
    fi

    # Define paths for the database, temporary, and final output directories
    local NAME=$(basename "$working_dir")
    local dbs_dir="${working_dir}/dbs_diamond"
    local results_temp_dir="${working_dir}/${NAME}_results_temp"
    local diamond_results_dir="${working_dir}/${NAME}_DiamondResults"
    local all_faa="${working_dir}/all.faa"
    local temp_file_prefix="${working_dir}/temp_"

    # Create the necessary directories
    mkdir -p "$dbs_dir"
    mkdir -p "$results_temp_dir"
    mkdir -p "$diamond_results_dir"

    # Step 1: Combine all protein files into a single database file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating a single combined protein database file..."
    cat ${protein_path}/*.fasta > "$all_faa"

    # Step 2: Make the Diamond database
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Building the Diamond database..."
    local db_file="${dbs_dir}/all"
    diamond makedb --in "$all_faa" -d "$db_file" --quiet

    # Step 3: Perform all-vs-all pairwise similarity searches
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing all-vs-all similarity searches..."
    
    # Get a list of all species filenames without the .fa extension
    local fa_files=( $(find "$protein_path" -maxdepth 1 -type f -name "*.fasta") )
    
    # Create a list of all species files to check against
    find "$protein_path" -maxdepth 1 -type f -name "*.fasta" | xargs -n 1 basename -s .fasta > "${temp_file_prefix}all"
    
    for fa_file in "${fa_files[@]}"
    do
        local species_name=$(basename "$fa_file" .fasta)
        local query="$fa_file"
        local db="$db_file"
        local outfile="${results_temp_dir}/${species_name}.tsv"
        
        # Perform the blastp search using the full query path
        diamond blastp -q "$query" -d "$db" -o "$outfile" --fast -k0 --max-hsps 1 --evalue 1e-5 --matrix PAM30 --query-cover 90 -p "$threads" --quiet
        
        # The rest of the loop is adjusted to use the new species_name variable
        awk '{print $2}' "$outfile" | sed 's/_.*//g' | sort | uniq > "${temp_file_prefix}hit"
        
        # Grep the hits against the list of all species to process each pair
        grep -f "${temp_file_prefix}hit" "${temp_file_prefix}all" | while read -r line
        do
            local TARGET="$line"
            if [[ "$species_name" != "$TARGET" ]]
            then
                grep "$TARGET" "$outfile" > "${diamond_results_dir}/${species_name}_${TARGET}.tsv"
            else
                # Special case for self-comparison
                grep ".*_${line}\.[0-9]*.*_${line}\.[0-9]*" "$outfile" > "${diamond_results_dir}/${species_name}_${TARGET}.tsv"
            fi
        done
        echo "    [$(date "+%Y-%m-%d %H:%M:%S")] ${species_name} finished Diamond analysis."
    done
    
    rm -r ${dbs_dir} ${results_temp_dir} ${all_faa} ${temp_file_prefix}*
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Analysis complete. Results stored in '$diamond_results_dir'."
    return 0
}

blastn_all_vs_all() {
    local fa_dir="$1"
    local output_file="$2"
    
    # --- Error Handling ---
    if [[ ! -d "$fa_dir" ]]; then
        echo "Error: Input directory '$fa_dir' not found." >&2
        exit 1
    fi

    # Define temporary file paths using a unique ID to prevent conflicts
    local temp_id=$(date +%s%N)
    local temp_fa="${fa_dir}/temp_fa_${temp_id}.fa"
    local temp_db="${fa_dir}/temp_fa_db_${temp_id}"
    local temp_exclude_1="${fa_dir}/temp_exclude_1_${temp_id}.txt"
    local temp_exclude_2="${fa_dir}/temp_exclude_2_${temp_id}.txt"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Starting all-vs-all BLASTN analysis..."
    
    # Step 1: Merge all FASTA files into a single temporary file
    cat "${fa_dir}"/*.fa > "$temp_fa"
    
    # Step 2: Create a BLAST database from the merged FASTA file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating the database for blast filtering stage..."
    makeblastdb -dbtype nucl -parse_seqids -in "$temp_fa" -out "$temp_db"
    
    # Step 3: Run all-against-all search and append results to output file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running blastn..."
    # Get all sequence headers from the temporary fa file
    grep ">" "$temp_fa" | sed 's/>//g' | while read -r line
    do
        # Create a list of the current query sequence ID
        blastdbcmd -entry "$line" -db "$temp_db" -dbtype nucl -outfmt %i >> "$temp_exclude_1"
        
        # Create an alias list with the single ID to be excluded
        blastdb_aliastool -seqid_file_in "$temp_exclude_1" -seqid_file_out "$temp_exclude_2"
        
        # Run BLASTN, excluding the query sequence from the search results
        # We assume the query file is named "$line.fa" in the input directory.
        blastn -query "${fa_dir}/${line}.fa" -db "$temp_db" \
               -negative_seqidlist "$temp_exclude_2" -task blastn \
               -gapopen 8 -gapextend 6 -reward 5 -penalty -4 \
               -evalue 1e-60 -num_threads "$threads" \
               -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" \
               >> "$output_file"
        echo "    [$(date "+%Y-%m-%d %H:%M:%S")] ${line} finished blast analysis."
    done
    
    echo "   [$(date "+%Y-%m-%d %H:%M:%S")] BLASTN analysis finished. Results are in '$output_file'."
    
    # The trap command will automatically remove temporary files upon function exit.
    rm -f '$temp_fa' '${temp_db}.*' '$temp_exclude_1' '$temp_exclude_2'
    return 0
}

process_collinearity() {
    local working_dir="$1"
    local mode="$2"

    # Locate required subdirectories and define output path
    local collinearity_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_collinearity" 2>/dev/null)
    local gff3_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local out_path="$working_dir"
    local metadata_file="$working_dir/metadata_updated.csv"

    # Check that the subdirectories were found
    if [[ -z "$collinearity_path" || -z "$gff3_path" ]]; then
        echo "Error: Required subdirectories (*_collinearity or *_gff) not found in '$working_dir'." >&2
        exit 1
    fi

    # Generate a unique ID for temporary files to prevent conflicts.
    local temp_id=$(date +%s)
    local temp_prefix="${out_path}/temp_${temp_id}"

    # Determine number of genes per element for percentage estimation
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Counting genes per element from GFF files..."
    find "${gff3_path}" -maxdepth 1 -type f -name "*.gff" -exec grep -c "gene" {} + | awk -F'/' '{ gsub(".gff:", "\t", $NF); print $NF }' | sort > "${temp_prefix}_genes_per_element.txt"

    # Determine the number of collinear genes per element
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Counting collinear genes from collinearity files..."
    find "${collinearity_path}" -name '*.collinearity' -exec grep -E "[0-9]*-.*[0-9]*:" {} + > "${temp_prefix}_raw_collinear_lines.txt"
        awk -F '/' '{print $NF}' "${temp_prefix}_raw_collinear_lines.txt" | \
    sed -e 's/ //g; s/\.collinearity//g; s/:/\t/g' | \
    awk '
    BEGIN { OFS="\t" }
    {
        # Use a composite key to count unique values for column 3 per column 1
        key3 = $1 SUBSEP $3
        if (! (key3 in count3)) {
            count3[key3] = 1
            elements3[$1]++
        }
        
        # Use a composite key to count unique values for column 4 per column 1
        key4 = $1 SUBSEP $4
        if (! (key4 in count4)) {
            count4[key4] = 1
            elements4[$1]++
        }
    }
    END {
        for (col1_val in elements3) {
            # Print the element and the counts for both columns
            print col1_val, elements3[col1_val], elements4[col1_val]
        }
    }' | \
    sed -e 's/_/\t/g' | sort > "${temp_prefix}_collinear_genes.txt"

    # Join the two previous files to generate a file for percentage estimation
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Joining gene count data..."
    join "${temp_prefix}_collinear_genes.txt" "${temp_prefix}_genes_per_element.txt" | sort -k2 > "${temp_prefix}_join_1.txt"
    join -1 2 -2 1 "${temp_prefix}_join_1.txt" "${temp_prefix}_genes_per_element.txt" | \
        awk '{swap=$1;$1=$2;$2=swap;print $0}' | sort > "${temp_prefix}_join_2.txt"

    # Estimate percentage per element
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Estimating percentage of collinearity per element..."
    awk '{print $0"\t"100*($3/$5)"\t"100*($4/$6)}' "${temp_prefix}_join_2.txt" > "${temp_prefix}_percentage_pairwise.txt"

    # Obtain the percentage returned by the software
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Extracting software-reported general percentages..."
    # Here is an error
    find "${collinearity_path}" -name '*.collinearity' -exec grep -w -E -o "Percentage: [0-9]*\.[0-9]*" {} + > "${temp_prefix}_raw_percentage_general.txt"
        awk -F '/' '{print $NF}' "${temp_prefix}_raw_percentage_general.txt" | sed -e 's/ //g; s/.collinearity//g; s/:/\t/g; s/_/\t/g' > "${temp_prefix}_percentage_general.txt"

    # Generate the final result
    awk '{if($4>0){print}}' "${temp_prefix}_percentage_general.txt" > "${temp_prefix}_percentage_general_filter.txt"

    echo -e "element01""\t""element02""\t""General_percentage""\t""element01_percentage""\t""element02_percentage" > "${out_path}/Collinearity_percentage.txt"

    if [[ "${mode}" == "Raw" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        # Here is the error.
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';' -r -g -k3,3 | sed -e 's/;/\t/g' >> "${out_path}/Collinearity_percentage.txt"
        
    elif [[ "${mode}" == "SSP" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | \
            awk '{if(($3 >= 41) || (($4 >= 45 || $5 >= 45) && ($4/$5 >= 1.8 || $4/$5 <= 0.55 ))) print}' | \
            sed -e 's/ /;/g' | sort -t ';' -r -g -k3,3 | sed -e 's/;/\t/g' >> "${out_path}/Collinearity_percentage.txt"
        
    elif [[ "${mode}" == "Filter" ]]
    then
        # Generate files that will be filter.
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';'  >> "${temp_prefix}_Collinearity_percentage.txt"
        awk 'BEGIN{FS=";";OFS=";"}{if(($3 >= 41) || (($4 >= 45 || $5 >= 45) && ($4/$5 >= 1.8 || $4/$5 <= 0.55 ))) print $0}' "${temp_prefix}_Collinearity_percentage.txt" > "${temp_prefix}_Collinearity_percentage_high.txt"
        awk 'BEGIN{FS=";";OFS=";"}{if($3 < 41) print $0}' "${temp_prefix}_Collinearity_percentage.txt" > "${temp_prefix}_Collinearity_percentage_low.txt"
        
        if [ -s "${temp_prefix}_Collinearity_percentage_low.txt" ]; then
        
            # Select elements to compare
            mkdir "${temp_prefix}_selected_nucleotides"
            awk -F ';' '{print $1; print $2}' "${temp_prefix}_Collinearity_percentage_low.txt" | sort -u > "${temp_prefix}_elements_low.txt"
            cat "${temp_prefix}_elements_low.txt" | while read -r line
            do
                cp "${nucleotide_dir}/${line}.fa" "${temp_prefix}_selected_nucleotides/"
            done
        
            # Run blast analysis
            blastn_all_vs_all "${temp_prefix}_selected_nucleotides" "${temp_prefix}_BlastnResults.out"
        
            # Filter blastn results
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Filtering Blastn results..."
            python3 ${pythonFunction_path} -f "${temp_prefix}_BlastnResults.out" -o "${temp_prefix}_BlastnClean.out"
        
            # Filter low syntenic pairs
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Filtering False Positive pairs..."
            awk '{OFS=";"} {print $1 OFS $2}' "${temp_prefix}_BlastnClean.out" > "${temp_prefix}_BlastnPairs.out"
            join -t ';' <(sed -e 's/;/-/' "${temp_prefix}_Collinearity_percentage_low.txt" | sort) <(sed -e 's/;/-/' "${temp_prefix}_BlastnPairs.out" | sort) | sed 's/-/;/g' > "${temp_prefix}_Collinearity_percentage_lowSelected.txt"
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
            cat "${temp_prefix}_Collinearity_percentage_high.txt" "${temp_prefix}_Collinearity_percentage_lowSelected.txt" | sort -u | sort -t ';' -r -g -k3,3 | sed -e 's/;/\t/g' > "${out_path}/Collinearity_percentage.txt"
        else
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Blast filtering is being skipped due to the absence of low syntenic pairs to evaluate..."
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
            sort -t ';' -r -g -k3,3 "${temp_prefix}_Collinearity_percentage_high.txt" | sed -e 's/;/\t/g' > "${out_path}/Collinearity_percentage.txt"
        fi    

    fi

    # Remove temporary data
    rm -rf ${temp_prefix}*
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Processing complete. Final report saved to ${out_path}/Collinearity_percentage.txt"
}

process_cluster_file() {
    local working_dir="$1"
    
    local gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    local exon_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_exon" 2>/dev/null)
    local cluster_path="$working_dir/Clusters.txt"

    echo "Processing file: $cluster_path"
    
    mkdir -p "${working_dir}/Clusters"
    
    # Set the Internal Field Separator to a tab for the outer loop.
    local OLD_IFS="$IFS"
    IFS=$'\t'

    while read -r cluster_id values_string; do

        IFS=' ' read -r -a values_array <<< "$values_string"
        
        mkdir -p "${working_dir}/Clusters/${cluster_id}"
        mkdir -p "${working_dir}/Clusters/${cluster_id}/${cluster_id}_gff"
        mkdir -p "${working_dir}/Clusters/${cluster_id}/${cluster_id}_nucleotide"
        mkdir -p "${working_dir}/Clusters/${cluster_id}/${cluster_id}_protein"
        mkdir -p "${working_dir}/Clusters/${cluster_id}/${cluster_id}_exon"
        
        for value in "${values_array[@]}"; do
            cp ${gff_dir}/${value}.gff ${working_dir}/Clusters/${cluster_id}/${cluster_id}_gff/
            cp ${protein_dir}/${value}.fa ${working_dir}/Clusters/${cluster_id}/${cluster_id}_protein/
            cp ${nucleotide_dir}/${value}.fa ${working_dir}/Clusters/${cluster_id}/${cluster_id}_nucleotide/
            cp ${exon_dir}/${value}.fa ${working_dir}/Clusters/${cluster_id}/${cluster_id}_exon/
        done 
    done < "$cluster_path"

    # Restore the original IFS at the end of the function.
    IFS="$OLD_IFS"
    return 0
}

Process_treefile() {
    local working_dir="$1"
    
    local Captain="${working_dir}/CaptainPhylogeny"
    local Microsynteny="${working_dir}/MicrosyntenyPhylogeny"
    local temp_file_prefix="${working_dir}/temp_"
    local CaptainPhylogeny_dir=$(find "$working_dir" -maxdepth 1 -type d -name "captainPhylogeny" 2>/dev/null)
    local Remove_elements="${working_dir}/Remove_outsider_elements.txt"
    

    #Create politomies based on branch length and node support
    gotree collapse length -l 0.00001 -i ${Microsynteny}/*.phy.treefile -o ${temp_file_prefix}cargolength.nw
    gotree collapse length -l 0.00001 -i ${CaptainPhylogeny_dir}/Captain_tree.treefile -o ${temp_file_prefix}captainlength.nw
    gotree collapse support -s 80 -i ${temp_file_prefix}cargolength.nw -o ${temp_file_prefix}cargosupport.nw 
    gotree collapse support -s 80 -i ${temp_file_prefix}captainlength.nw -o ${temp_file_prefix}captainsupport.nw

    sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_file_prefix}cargosupport.nw 
    sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_file_prefix}captainsupport.nw

    gotree collapse support -s 80 -i ${temp_file_prefix}cargosupport.nw  -o ${temp_file_prefix}cargosupport2.nw 
    gotree collapse support -s 80 -i ${temp_file_prefix}captainsupport.nw -o ${temp_file_prefix}captainsupport2.nw

    sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_file_prefix}cargosupport2.nw 
    sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_file_prefix}captainsupport2.nw

    mv ${temp_file_prefix}cargosupport2.nw ${temp_file_prefix}cargosupport.nw
    mv ${temp_file_prefix}captainsupport2.nw ${temp_file_prefix}captainsupport.nw

    # Check any issues with the tips
    gotree labels --input ${temp_file_prefix}cargosupport.nw | sort > ${temp_file_prefix}CargoTips
    gotree labels --input ${temp_file_prefix}captainsupport.nw | sort > ${temp_file_prefix}CaptainTips

    if ! diff -q "${temp_file_prefix}CargoTips" "${temp_file_prefix}CaptainTips" >/dev/null; then

        local NumberCargo=$(wc -l ${temp_file_prefix}CargoTips | awk '{print $1}')
        local NumberCaptain=$(wc -l ${temp_file_prefix}CaptainTips | awk '{print $1}')

        if [[ ! "$NumberCargo" -eq "$NumberCaptain" ]]; then
            echo -e "
             there's difference in the number of tips between captain and Cargo tree. Comparison of trees is going to be skipped..."
            grep -v -f ${temp_file_prefix}CargoTips ${temp_file_prefix}CaptainTips > ${Remove_elements}
        else
            paste ${temp_file_prefix}CargoTips ${temp_file_prefix}CaptainTips > ${temp_file_prefix}mapfile
            gotree rename -m ${temp_file_prefix}mapfile -i ${temp_file_prefix}cargosupport.nw -o ${temp_file_prefix}cargosupport2.nw 
            mv ${temp_file_prefix}cargosupport2.nw ${temp_file_prefix}cargosupport.nw

            #Root
            gotree reroot midpoint -i ${temp_file_prefix}cargosupport.nw -o ${working_dir}/CargoMicrosyntenyPhylogeny.nw
            gotree reroot midpoint -i ${temp_file_prefix}captainsupport.nw -o ${working_dir}/CaptainPhylogeny.nw

            #Compare trees
            gotree compare trees -i ${working_dir}/CaptainPhylogeny.nw -c ${working_dir}/CargoMicrosyntenyPhylogeny.nw > ${working_dir}/Phylogenies_comparison.txt

            sed -i 's/reference/Captain/g' ${working_dir}/Phylogenies_comparison.txt
            sed -i 's/compared/Cargo/g' ${working_dir}/Phylogenies_comparison.txt
        fi  
    else
        #Root
        gotree reroot midpoint -i ${temp_file_prefix}cargosupport.nw -o ${working_dir}/CargoMicrosyntenyPhylogeny.nw
        gotree reroot midpoint -i ${temp_file_prefix}captainsupport.nw -o ${working_dir}/CaptainPhylogeny.nw

        #Compare trees
        gotree compare trees -i ${working_dir}/CaptainPhylogeny.nw -c ${working_dir}/CargoMicrosyntenyPhylogeny.nw > ${working_dir}/Phylogenies_comparison.txt

        sed -i 's/reference/Captain/g' ${working_dir}/Phylogenies_comparison.txt
        sed -i 's/compared/Cargo/g' ${working_dir}/Phylogenies_comparison.txt

    fi

    rm ${temp_file_prefix}*
}

Removed_outsider_elements() {

    local working_dir="$1"

    local gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local exon_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_exon" 2>/dev/null)

    local Remove_elements="${working_dir}/Remove_outsider_elements.txt"
    local output="${working_dir}/Remove_outsider_elements/"

    mkdir -p "${output}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Removing elements that are outsiders of this cluster based on new cargo evidence."

    cat "${Remove_elements}" | while read line
    do
        mv ${gff_dir}/${line}.gff ${output}
        mv ${protein_dir}/${line}.fa ${output}${line}_protein.fa
        mv ${nucleotide_dir}/${line}.fa ${output}${line}_nucleotide.fa
        mv ${exon_dir}/${line}.fa ${output}${line}_exon.fa
    done
}

# ==============================================================================
# Show parameters selected
# ==============================================================================
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running Syntenet module analysis with '${mode}' mode."

echo "For Step 4 (Syntenet analysis), using the following parameter:"
echo "  Minimum anchor points: $anchorPoints."
echo "  Maximum allowed gaps: $gaps."

# ==============================================================================
# Checking Working directory structure
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Checking Working directory '${workingDirectory_path}' structure."
check_directory_structure "${workingDirectory_path}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure in '${workingDirectory_path}' is valid. Proceeding."

# ==============================================================================
# Preprocessing data
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Preprocessing data for Diamond."

Rscript -e '
  # Check and Call directory argument
  args <- commandArgs(trailingOnly = TRUE)
  dir <- as.character(args[1])
  
  # Check software installation
  suppressPackageStartupMessages(library(syntenet))
  
  if (!requireNamespace("syntenet", quietly = TRUE)) {
  stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
  }

  # Call directories in the Working directory
  internal_dir <- list.dirs(dir, recursive = F)
  
  # Import data
  gff_dir <- as.character(grep("gff", internal_dir, value = TRUE))
  fa_dir <- as.character(grep("protein", internal_dir, value = TRUE))
  
  proteomes <- fasta2AAStringSetlist(fa_dir)
  annotation <- gff2GRangesList(gff_dir)
  
  # Preprocess files
  pdata <- process_input(proteomes, annotation, gene_field = "ID")
  
  #export processed sequences
  temp <- export_sequences(pdata$seq,outdir = paste(dir,"/",head(unlist(strsplit(tail(unlist(strsplit(dir,"/")),1),"_")),1),"_ProcessedProtein", sep = ""))

' "${workingDirectory_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Preprocessing finished. Proceeding."

# ==============================================================================
# Running Diamond analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running Diamond analysis."

process_diamond "${workingDirectory_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Diamond analysis finished. Proceeding."

# ==============================================================================
# Running Diamond analysis
# ==============================================================================


echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Running syntenet analysis."

mkdir ${workingDirectory_path}/MicrosyntenyPhylogeny

Rscript -e '
  # Check and Call directory argument
  args <- commandArgs(trailingOnly = TRUE)
  dir <- as.character(args[1])
  ap <- as.integer(args[2])
  mg <- as.integer(args[3])
  
  # Check software installation
  suppressPackageStartupMessages(library(syntenet))
  suppressPackageStartupMessages(library(labdsv))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library("svglite"))
  suppressPackageStartupMessages(library("ape"))
  
  
  if (!requireNamespace("syntenet", quietly = TRUE)) {
  stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
  }
  
  # Call directories in the Working directory
  internal_dir <- list.dirs(dir, recursive = F)

  # Import data
  gff_dir <- as.character(grep("gff", internal_dir, value = TRUE))
  fa_dir <- as.character(grep("protein", internal_dir, value = TRUE))
  
  proteomes <- fasta2AAStringSetlist(fa_dir)
  annotation <- gff2GRangesList(gff_dir)
  
  #preprocess files
  pdata <- process_input(proteomes, annotation, gene_field = "ID")
  
  #import blast results
  blast_dir <- as.character(grep("DiamondResults", internal_dir, value = TRUE))
  blast_list <- read_diamond(blast_dir)
  
  #perform synteny block analysis
  intersyn <- interspecies_synteny(blast_list, pdata$annotation, inter_dir = paste(dir,"/",head(unlist(strsplit(tail(unlist(strsplit(dir,"/")),1),"_")),1),"_collinearity", sep = ""), anchors = ap, max_gaps = mg)
  
  net <- parse_collinearity(intersyn)
  
  id_table <- create_species_id_table(names(proteomes))

  clusters <- cluster_network(net)

  profiles <- phylogenomic_profile(clusters)

  ProfileImage <- plot_profiles(
    profiles, 
    dist_function = labdsv::dsvdis,
    dist_params = list(index = "ruzicka")
  )
  
  ggsave(path = dir,filename = "Syntenet_profiling.svg", plot = ProfileImage, width = 16, height = 16, unit = "in")
  
  bt_mat <- binarize_and_transpose(profiles)

  if(iqtree_is_installed()) {
    phylo <- infer_microsynteny_phylogeny(bt_mat, outdir = paste(dir,"/MicrosyntenyPhylogeny/", sep = ""), threads = 1, model = "MFP")
  }

' "${workingDirectory_path}" "${anchorPoints}" "${gaps}"

rm Rplots.pdf

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Syntenet analysis finished. Proceeding."

# ==============================================================================
# Running Diamond analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Summarizing syntenet results in '${mode}' mode."

process_collinearity "${workingDirectory_path}" "${mode}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Summary finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished synteny and profiling analysis. Data is stored in ${workingDirectory_path}"

# ==============================================================================
# Running Clustering analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 6: Checking if there's a break in the cluster."

python3 ${pythonFunction_path}/SimpleCluster.py -i "${workingDirectory_path}/Collinearity_percentage.txt" -o "${workingDirectory_path}" 

ClusterNumber=$(wc -l ${workingDirectory_path}/Clusters.txt | awk '{print $1}')

if [[ $ClusterNumber -ge 2 ]]; then
    echo -e "\033[01;31mWARNING\033[m: There are more than one cluster. Breaking the files into the new clusters."
    process_cluster_file "${workingDirectory_path}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Comparison between captain tree and cargo mycrosynteny phylogenetic tree is skipped ..."
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished analysis. Data is stored in ${workingDirectory_path}. \033[01;31mWARNING\033[m: Analysis could not be completely reliable, please run again Captain phylogeny and this script for the individual clusters generated here."
else
    rm ${workingDirectory_path}/Clusters.txt
    rm ${workingDirectory_path}/network_edges.txt
    rm ${workingDirectory_path}/node_attributes.txt
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 7: Phylogenetic tree"
    Process_treefile "${workingDirectory_path}"

    if [ -s "${workingDirectory_path}/Remove_outsider_elements.txt" ]; then
        Removed_outsider_elements "${workingDirectory_path}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${workingDirectory_path}."
        echo -e "\033[01;31mWARNING\033[m: Elements have been removed, please check file '${workingDirectory_path}/Remove_elements.txt' and folder '${workingDirectory_path}/RemovedElements'."
    else
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished analysis. Data is stored in ${workingDirectory_path}"
    fi 
fi
