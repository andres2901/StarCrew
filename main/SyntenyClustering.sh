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
   │   ├── element01.fasta
   │   └── element02.fasta
   │   ︙
   └── *_protein/
   │   ├── element01.fa
   │   └── element02.fa
   │   ︙
   └── metadata_updated.csv

   NOTE: the metadata_updated.csv file is optional."
   echo
   echo "Syntax: SAT SyntenyClustering [ -h ] -w <working_directory> [ -m <mode> -a <anchor_points> -g <gaps> -n <min_nodes> -s <min_size> -t <threshold> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-m, --mode: Specified the mode to run the summarizing process of synteny results (Available mode: Raw, SSP, Filter) (Default = Raw)."
   echo "-a, --anchors: Number of minimum anchor points for syntenet to call a collinear region (Integer between 5 and 25) (Default = 8)."
   echo "-g, --gaps: Number of maximum allowed gaps between anchor points for syntenet to call a collinear region (Integer between 5 and 25) (Default = 8)."
   echo "-n, --minNodes: Minimum number of nodes in a cluster for spectral clustering to be attempted (Default: 4)."
   echo "-s, --minSize: The minimum desired size for any final sub-cluster (Default: 2)."
   echo "-t, --threshold: The minimum modularity score for a split to be accepted. Range [-0.5, 1.0] (Default: 0.05)."
   echo "-st, --searchThreads: Number of threads for searching software (DIAMOND and blast) (Default: 8)"
   echo "-p, --predictionMode: Mode of the gene prediction performed for the data (Available mode: Quick, Robust) (Default = Quick)"
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
mode="Raw"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
anchorPoints="8"
gaps="8"
minSize="2"
minNodes="4"
threshold="0.05"
threads="8"
predictionMode="Quick"
help_flag=false
metadata_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workingDirectory)
            shift
            Working_directory="$1"
            ;;
        -m|--mode)
            shift
            mode="$1"
            ;;
        -a|--anchors)
            shift
            anchorPoints="$1"
            ;;
        -g|--gaps)
            shift
            gaps="$1"
            ;;
        -n|--minNodes)
            shift
            minNodes="$1"
            ;;
        -s|--minSize)
            shift
            minSize="$1"
            ;;
        -t|--threshold)
            shift
            threshold="$1"
            ;;
        -st|--searchThreads)
            shift
            threads="$1"
            ;;
        -p|--predictionMode)
            shift
            predictionMode="$1"
            ;;
        -help)
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
if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if working directory exists
if [[ ! -d "$Working_directory" ]]; then
    echo "Error: Directory '$Working_directory' does not exist."
    exit 1
else
    Working_directory=$(realpath $Working_directory)
fi

# Check if mode parameter is correct
if [[ "$mode" != "Raw" && "$mode" != "SSP" && "$mode" != "Filter" ]]; then
    echo "Error: provided mode '$mode' is not accepted."
    print_help
    exit 1
fi

# Check if directory where the python function are stored exists
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: Directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
fi

# Check if python 'Blast_CleanUp.py' and 'Clustering.py' scripts exists in the given path 
if [[ "${mode}" == "Filter" ]]; then
    if [[ ! -f "${auxiliary_path}/Blast_CleanUp.py" ]]; then
        echo "Error: File '${auxiliary_path}/Blast_CleanUp.py' does not exist."
        exit 1
    fi
fi

if [[ ! -f "${auxiliary_path}/Clustering.py" ]]; then
    echo "Error: File '${auxiliary_path}/Clustering.py' does not exist."
    exit 1
fi

# Check anchor points and maximum gaps are within allowed range
if [[ "$anchorPoints" =~ ^[0-9]+$ ]]; then
    if (( anchorPoints < 5 || anchorPoints > 25 )); then
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

# Check parameters for subclustering

if [[ ! "$minNodes" =~ ^[0-9]+$ ]]; then
    echo "Error: '$minNodes' is not a positive integer."
    print_help
    exit 1
fi

if [[ ! "$minSize" =~ ^[0-9]+$ ]]; then
    echo "Error: '$minSize' is not a positive integer."
    print_help
    exit 1
fi

if [[ "$threshold" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    if (( $(echo "$threshold < -0.5" | bc -l) )) && (( $(echo "$threshold > 1.0" | bc -l) )); then
        echo "Error: '$threshold' is not an accepted value for clustering modularity score."
        print_help
        exit 1
    fi
else
    echo "Error: '$threshold' is not a float."
    print_help
    exit 1
fi

# Check thread parameter
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
    exit 1
fi

# Check if prediction mode parameter is correct
if [[ "$predictionMode" != "Quick" && "$predictionMode" != "Robust" ]]; then
    echo "Error: provided mode '$predictionMode' is not accepted."
    print_help
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

check_directory_structure() {
    local base_dir="$1"
    
    # Step 1: Locate required subdirectories and file
    local workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    fi

    local data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

    if [[ -z "$data_dir" ]]; then
        echo "Error: Data directory not found in '$base_dir'." >&2
        exit 1
    fi

    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    local metadata_file="${data_dir}/metadata.csv"

    if [[ -f "$metadata_file" ]]; then
        metadata_flag=true
    fi

    if [[ -z "$gff_dir" ]]; then
        echo "Error: GFF subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    if [[ -z "$protein_dir" ]]; then
        echo "Error: Protein subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    if [[ -z "$nucleotide_dir" ]]; then
        echo "Error: nucleotide subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    if [[ -z "$exon_dir" ]]; then
        echo "Error: exon subdirectory not found in '$data_dir'." >&2
        exit 1
    fi
    
    # Step 2: Check for consistent filenames across subdirectories

    # Define temporary file paths in the working directory
    local gff_files="$base_dir/gff_files.txt"
    local protein_files="$base_dir/protein_files.txt"
    local nucleotide_files="$base_dir/nucleotide_files.txt"
    local exon_files="$base_dir/exon_files.txt"
    
    # Get sorted list of base filenames from the GFF directory
    find "$gff_dir" -maxdepth 1 -type f -name "*.gff" | xargs -n 1 basename -s .gff | sort > "$gff_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$protein_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$protein_files"

    # Get sorted list of base filenames from the GFF directory
    find "$nucleotide_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$nucleotide_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$exon_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$exon_files"

    # Compare the lists. If diff finds a difference, it returns a non-zero exit code.
    if ! diff -q "$gff_files" "$protein_files" >/dev/null || \
        ! diff -q "$gff_files" "$nucleotide_files" >/dev/null || \
        ! diff -q "$gff_files" "$exon_files" >/dev/null || \
        ! diff -q "$protein_files" "$nucleotide_files" >/dev/null || \
        ! diff -q "$protein_files" "$exon_files" >/dev/null || \
        ! diff -q "$nucleotide_files" "$exon_files" >/dev/null; then
        echo "Error: File lists in subdirectories do not match." >&2
        echo "Details:" >&2
        echo "GFF vs. Protein:" >&2
        diff "$gff_files" "$protein_files" >&2
        echo "GFF vs. nucleotide:" >&2
        diff "$gff_files" "$nucleotide_files" >&2
        echo "GFF vs. exon:" >&2
        diff "$gff_files" "$exon_files" >&2
        echo "protein vs. nucleotide:" >&2
        diff "$protein_files" "$nucleotide_files" >&2
        echo "protein vs. exon:" >&2
        diff "$protein_files" "$exon_files" >&2
        echo "nucleotide vs. exon:" >&2
        diff "$nucleotide_files" "$exon_files" >&2
        rm "$gff_files" "$nucleotide_files" "$protein_files" "$exon_files"
        exit 1
    fi

    # Cleanup temporary files
    rm "$gff_files" "$nucleotide_files" "$protein_files" "$exon_files"
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    local cluster_dir="${base_dir}/Clusters/"
    local working_dir="${base_dir}/Workspace/SyntenyClustering/"
    local temp_dir="${base_dir}/Workspace/SyntenyClustering/temp/"

    # Create required subdirectories
    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}
    mkdir -p ${cluster_dir}

    # Copy require files
    cp -r ${gff_dir} ${working_dir}
    cp -r ${nucleotide_dir} ${working_dir}
    cp -r ${protein_dir} ${working_dir}
    cp -r ${exon_dir} ${working_dir}
    if $metadata_flag; then
        cp "${data_dir}/metadata.csv" ${working_dir}
    fi
}

process_diamond() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/SyntenyClustering/"
    local temp_dir="${working_dir}/temp/"

    # Locate the preprocess protein directory 
    local protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "PreprocessData" 2>/dev/null)
    
    if [[ -z "$protein_path" ]]; then
        echo "Error: PreprocessData directory not found in '$working_dir'." >&2
        exit 1
    fi

    # Define paths for the database, temporary, and final output directories
    local dbs_dir="${working_dir}/dbs_diamond"
    local results_temp_dir="${temp_dir}/results_temp"
    local diamond_results_dir="${working_dir}/DiamondResults"
    local all_faa="${working_dir}/all.faa"

    # Create the necessary directories
    mkdir -p "$dbs_dir"
    mkdir -p "$results_temp_dir"
    mkdir -p "$diamond_results_dir"

    # Combine all protein files into a single database file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating a single combined protein database file..."
    cat ${protein_path}/*.fasta > "$all_faa"

    # Make the Diamond database
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Building the Diamond database..."
    local db_file="${dbs_dir}/all"
    diamond makedb --in "$all_faa" -d "$db_file" --quiet

    # Perform all-vs-all pairwise similarity searches
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing all-vs-all similarity searches..."
    
    # Get a list of all species filenames without the .fasta extension
    local fasta_files=( $(find "$protein_path" -maxdepth 1 -type f -name "*.fasta") )
    
    # Create a list of all species files to check against
    find "$protein_path" -maxdepth 1 -type f -name "*.fasta" | xargs -n 1 basename -s .fasta > "${temp_dir}all"
    
    for fasta_file in "${fasta_files[@]}"
    do
        local species_name=$(basename "$fasta_file" .fasta)
        local query="$fasta_file"
        local db="$db_file"
        local outfile="${results_temp_dir}/${species_name}.tsv"
        
        # Perform the blastp search using the full query path
        diamond blastp -q "$query" -d "$db" -o "$outfile" --fast -k0 --max-hsps 1 --evalue 1e-5 --matrix PAM30 --query-cover 90 -p "$threads" --quiet
        
        # The rest of the loop is adjusted to use the new species_name variable
        awk '{print $2}' "$outfile" | sed 's/_.*//g' | sort | uniq > "${temp_dir}hit"
        
        # Grep the hits against the list of all species to process each pair
        grep -f "${temp_dir}hit" "${temp_dir}all" | while read -r line
        do
            local TARGET="$line"
            if [[ "$species_name" != "$TARGET" ]]
            then
                grep "${TARGET}_" "$outfile" > "${diamond_results_dir}/${species_name}_${TARGET}.tsv"
            else
                # Special case for self-comparison
                grep ".*_${line}\.[0-9]*.*_${line}\.[0-9]*" "$outfile" > "${diamond_results_dir}/${species_name}_${TARGET}.tsv"
            fi
        done

        rm "${temp_dir}hit"
        echo "    [$(date "+%Y-%m-%d %H:%M:%S")] ${species_name} finished Diamond analysis."
    done
    
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Analysis complete. Results stored in '$diamond_results_dir'."
}

blastn_all_vs_all() {
    local fasta_dir="$1"
    local output_file="$2"

    local temp_dir="${fasta_dir}/../"

    if [[ ! -d "$fasta_dir" ]]; then
        echo "Error: Directory '$fasta_dir' not found." >&2
        exit 1
    fi

    # Define temporary file paths using a unique ID to prevent conflicts
    local temp_id=$(date +%s%N)
    local temp_fasta="${temp_dir}/temp_fasta_${temp_id}.fa"
    local temp_db="${temp_dir}/temp_fasta_db_${temp_id}"
    local temp_exclude_1="${temp_dir}/temp_exclude_1_${temp_id}.txt"
    local temp_exclude_2="${temp_dir}/temp_exclude_2_${temp_id}.txt"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Starting all-vs-all BLASTN analysis..."
    
    # Step 1: Merge all FASTA files into a single temporary file
    cat "${fasta_dir}"/*.fa > "$temp_fasta"
    
    # Step 2: Create a BLAST database from the merged FASTA file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating the database for blast filtering stage..."
    makeblastdb -dbtype nucl -parse_seqids -in "$temp_fasta" -out "$temp_db"
    
    # Step 3: Run all-against-all search and append results to output file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running blastn..."
    # Get all sequence headers from the temporary fasta file
    grep ">" "$temp_fasta" | sed 's/>//g' | while read -r line
    do
        # Create a list of the current query sequence ID
        blastdbcmd -entry "$line" -db "$temp_db" -dbtype nucl -outfmt %i >> "$temp_exclude_1"
        
        # Create an alias list with the single ID to be excluded
        blastdb_aliastool -seqid_file_in "$temp_exclude_1" -seqid_file_out "$temp_exclude_2"
        
        # Run BLASTN, excluding the query sequence from the search results
        # We assume the query file is named "$line.fa" in the input directory.
        blastn -query "${fasta_dir}/${line}.fa" -db "$temp_db" \
               -negative_seqidlist "$temp_exclude_2" -task blastn \
               -gapopen 8 -gapextend 6 -reward 5 -penalty -4 \
               -evalue 1e-60 -num_threads "$threads" \
               -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" \
               >> "$output_file"
        echo "    [$(date "+%Y-%m-%d %H:%M:%S")] ${line} finished blast analysis."
    done
    
    echo "   [$(date "+%Y-%m-%d %H:%M:%S")] BLASTN analysis finished. Results are in '$output_file'."
    
    # The trap command will automatically remove temporary files upon function exit.
    rm -f '$temp_fasta' '${temp_db}.*' '$temp_exclude_1' '$temp_exclude_2'
    return 0
}

process_collinearity() {
    local base_dir="$1"
    local mode="$2"

    local working_dir="${base_dir}/Workspace/SyntenyClustering/"

    # Locate required subdirectories and define output path
    local collinearity_path=$(find "$working_dir" -maxdepth 1 -type d -name "Collinearity" 2>/dev/null)
    local gff_path=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)

    local temp_dir="${working_dir}/temp/"
    local out_path="$working_dir"

    if $metadata_flag; then
        local metadata_file="${working_dir}/metadata.csv"
    fi

    # Check that the subdirectories were found
    if [[ -z "$collinearity_path" || -z "$gff_path" || -z $nucleotide_dir ]]; then
        echo "Error: Required subdirectories (Collinearity, Gff or Nucleotide) not found in '$working_dir'." >&2
        exit 1
    fi

    # Generate a unique ID for temporary files to prevent conflicts.
    local temp_id=$(date +%s)
    local temp_prefix="${temp_dir}/temp_${temp_id}"

    # Determine number of genes per element for percentage estimation
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Counting genes per element from GFF files..."
    find "${gff_path}" -maxdepth 1 -type f -name "*.gff" -exec grep -c "gene" {} + | awk -F'/' '{ gsub(".gff:", "\t", $NF); print $NF }' | sort > "${temp_prefix}_genes_per_element.txt"

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
    if $metadata_flag; then
        cut -d ';' -f 2- "${metadata_file}" | sort -t ';' > "${temp_prefix}_metadata.txt"
        echo -e "element01""\t""element02""\t""General_percentage""\t""element01_percentage""\t""element02_percentage""\t""element01_species""\t""element01_length""\t""element02_species""\t""element02_length" > "${base_dir}/Collinearity_percentage.txt"
    else
        echo -e "element01""\t""element02""\t""General_percentage""\t""element01_percentage""\t""element02_percentage" > "${base_dir}/Collinearity_percentage.txt"
    fi

    awk '{if($4>0){print}}' "${temp_prefix}_percentage_general.txt" > "${temp_prefix}_percentage_general_filter.txt"

    if [[ "${mode}" == "Raw" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        # Here is the error.
        if $metadata_flag; then
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';'  >> "${temp_prefix}_Collinearity_percentage.txt"
        join -t ';' -j 1 -a1 "${temp_prefix}_Collinearity_percentage.txt" "${temp_prefix}_metadata.txt" > "${temp_prefix}_join.txt"
        join -t ';' -1 2 -2 1 -a1 <( sort -t ';' -k2,2 "${temp_prefix}_join.txt") "${temp_prefix}_metadata.txt" | awk 'BEGIN{FS=";";OFS=";"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -t ';' -r -g -k3,3 | sed -e 's/;/\t/g' >> "${working_dir}/Collinearity_percentage.txt"
        else
            awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';' | sed -e 's/;/\t/g' >> "${working_dir}/Collinearity_percentage.txt"
        fi
    elif [[ "${mode}" == "SSP" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        if $metadata_flag; then
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | \
            awk '{if(($3 >= 41) || (($4 >= 45 || $5 >= 45) && ($4/$5 >= 1.8 || $4/$5 <= 0.55 ))) print}' | \
            sed -e 's/ /;/g' | sort -t ';' > "${temp_prefix}_Collinearity_percentage.txt"
        join -t ';' -a1 "${temp_prefix}_Collinearity_percentage.txt" "${temp_prefix}_metadata.txt" > "${temp_prefix}_join.txt"
        join -t ';' -1 2 -a1 <( sort -t ';' -k2,2 "${temp_prefix}_join.txt") "${temp_prefix}_metadata.txt" | awk 'BEGIN{FS=";";OFS=";"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -t ';' -r -g -k3,3 | sed -e 's/;/\t/g' >> "${working_dir}/Collinearity_percentage.txt"
        else
            awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | \
            awk '{if(($3 >= 41) || (($4 >= 45 || $5 >= 45) && ($4/$5 >= 1.8 || $4/$5 <= 0.55 ))) print}' | \
            sed -e 's/ /;/g' | sort -t ';' | sed -e 's/;/\t/g' >> "${working_dir}/Collinearity_percentage.txt"
        fi
    elif [[ "${mode}" == "Filter" ]]
    then
        # Generate files that will be filter.
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';'  >> "${temp_prefix}_Collinearity_percentage.txt"
        awk 'BEGIN{FS=";";OFS=";"}{if(($3 >= 41) || (($4 >= 45 || $5 >= 45) && ($4/$5 >= 1.8 || $4/$5 <= 0.55 ))) print $0}' "${temp_prefix}_Collinearity_percentage.txt" > "${temp_prefix}_Collinearity_percentage_high.txt"
        awk 'BEGIN{FS=";";OFS=";"}{if($3 < 41) print $0}' "${temp_prefix}_Collinearity_percentage.txt" > "${temp_prefix}_Collinearity_percentage_low.txt"
        
        if [ -s "${temp_prefix}_Collinearity_percentage_low.txt" ]
        then
            # Select elements to compare
            mkdir "${temp_prefix}_selected_nucleotides"
            awk -F ';' '{print $1; print $2}' "${temp_prefix}_Collinearity_percentage_low.txt" | sort -u > "${temp_prefix}_elements_low.txt"
            cat "${temp_prefix}_elements_low.txt" | while read -r line
            do
                cp "${nucleotide_dir}/${line}.fa" "${temp_prefix}_selected_nucleotides/"
            done
        
            # Run blast analysis
            blastn_all_vs_all "${temp_prefix}_selected_nucleotides" "${working_dir}/BlastnResults.out"
        
            # Filter blastn results
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Filtering Blastn results..."
            python3 ${auxiliary_path}/Blast_CleanUp.py -f "${working_dir}/BlastnResults.out" -o "${working_dir}BlastnClean.out"
        
            # Filter low syntenic pairs
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Filtering False Positive pairs..."
            awk '{OFS=";"} {print $1 OFS $2}' "${working_dir}/BlastnClean.out" > "${working_dir}/BlastnPairs.out"
            join -t ';' <(sed -e 's/;/-/' "${temp_prefix}_Collinearity_percentage_low.txt" | sort) <(sed -e 's/;/-/' "${working_dir}/BlastnPairs.out" | sort) | sed 's/-/;/g' > "${temp_prefix}_Collinearity_percentage_lowSelected.txt"
            cat "${temp_prefix}_Collinearity_percentage_high.txt" "${temp_prefix}_Collinearity_percentage_lowSelected.txt" | sort -u > "${temp_prefix}_Collinearity_percentage_filter.txt"
        else
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Blast filtering is being skipped due to the absence of low syntenic pairs to evaluate..."
            mv "${temp_prefix}_Collinearity_percentage_high.txt" "${temp_prefix}_Collinearity_percentage_filter.txt"
        fi   
        
        # Final stage
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        if $metadata_flag; then
            join -t ';' -a1 "${temp_prefix}_Collinearity_percentage_filter.txt" "${temp_prefix}_metadata.txt" > "${temp_prefix}_join.txt"
            join -t ';' -1 2 -a1 <( sort -t ';' -k2,2 "${temp_prefix}_join.txt") "${temp_prefix}_metadata.txt" | awk 'BEGIN{FS=";";OFS=";"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -t ';' -r -g -k3,3 | sed -e 's/;/\t/g' >> "${working_dir}/Collinearity_percentage.txt"
        else
            cat "${temp_prefix}_Collinearity_percentage_filter.txt" | sed -e 's/;/\t/g' >> "${working_dir}/Collinearity_percentage.txt"
        fi
    fi

    # Remove temporary data
    rm -rf ${temp_prefix}*
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Processing complete. Final report saved to ${working_dir}/Collinearity_percentage.txt"
}

process_cluster_file() {
    local base_dir="$1"

    local cluster_dir="${base_dir}/Clusters/"
    local working_dir="${base_dir}/Workspace/SyntenyClustering/"
    
    local gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local exon_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    if $metadata_flag; then
        local metadata_file="${working_dir}/metadata.csv"
    fi

    local cluster_path="${cluster_dir}/main_clusters.txt"

    # --- Error Handling ---
    if [[ ! -f "$cluster_path" ]]; then
        echo "Error: File '$cluster_path' not found." >&2
        return 1
    fi
    
    # Set the Internal Field Separator to a tab for the outer loop.
    local OLD_IFS="$IFS"
    IFS=$'\t'

    if [[ "$predictionMode" == "Quick" ]]
    then
        while read -r cluster_id values_string; do

            IFS=' ' read -r -a values_array <<< "$values_string"
        
            echo "Processing: $cluster_id with ${#values_array[@]} elements."
        
            mkdir -p "${cluster_dir}/${cluster_id}"
            mkdir -p "${cluster_dir}/${cluster_id}/Workspace"
            mkdir -p "${cluster_dir}/${cluster_id}/Data"

            if $metadata_flag; then
                local updated_metadata="${cluster_dir}/${cluster_id}/Data/metadata.csv"
            fi
        
            for value in "${values_array[@]}"; do
                cat ${nucleotide_dir}/${value}.fa >> ${cluster_dir}/${cluster_id}/sequences.fa      
                if $metadata_flag; then
                    grep -w $value $metadata_file >> $updated_metadata
                fi   
            done
        done < "$cluster_path"
    elif [[ "$predictionMode" == "Robust" ]] 
    then
        while read -r cluster_id values_string; do

            IFS=' ' read -r -a values_array <<< "$values_string"

            echo "Processing: $cluster_id with ${#values_array[@]} elements."
        
            mkdir -p "${cluster_dir}/${cluster_id}"
            mkdir -p "${cluster_dir}/${cluster_id}/Workspace"
            mkdir -p "${cluster_dir}/${cluster_id}/Data"
            mkdir -p "${cluster_dir}/${cluster_id}/Data/Exon"
            mkdir -p "${cluster_dir}/${cluster_id}/Data/Nucleotide"
            mkdir -p "${cluster_dir}/${cluster_id}/Data/Protein"
            mkdir -p "${cluster_dir}/${cluster_id}/Data/Gff"

            if $metadata_flag; then
                local updated_metadata="${cluster_dir}/${cluster_id}/Data/metadata.csv"
            fi
        
            for value in "${values_array[@]}"; do
                cp ${gff_dir}/${value}.gff ${cluster_dir}/C${cluster_id}/Gff/
                cp ${protein_dir}/${value}.fa ${cluster_dir}/${cluster_id}/Protein/
                cp ${nucleotide_dir}/${value}.fa ${cluster_dir}/${cluster_id}/Nucleotide/    
                cp ${exon_dir}/${value}.fa ${cluster_dir}/${cluster_id}/Exon/
                if $metadata_flag; then
                    grep -w $value $metadata_file >> $updated_metadata
                fi        
            done
        done < "$cluster_path"
    fi

    # Restore the original IFS at the end of the function.
    IFS="$OLD_IFS"
    echo "Processing complete."
}

# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running Syntenet module analysis with '${mode}' mode."

echo "For Syntenet analysis, using the following parameter:"
echo "  Minimum anchor points: $anchorPoints."
echo "  Maximum allowed gaps: $gaps."

echo "For Subclustering, using the following parameter:"
echo "  Minimum number of elements in a cluster for subclustering to be attempted: $minSize."
echo "  Minimum desired size for any final sub-cluster: $minNodes."
echo "  Minimum modularity score for a split to be accepted: $threshold."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${Working_directory}' structure."
check_directory_structure "${Working_directory}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"

organize_working_directory "${Working_directory}"

# ==============================================================================
# Preprocessing data
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Preprocessing data for Diamond."

Rscript ${auxiliary_path}/syntenetPreprocess.R "${Working_directory}/Workspace/SyntenyClustering/"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

# ==============================================================================
# Running Diamond analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running Diamond analysis."

process_diamond "${Working_directory}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

# ==============================================================================
# Running Diamond analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running syntenet analysis."
Rscript ${auxiliary_path}/syntenetAnalysis.R  -d "${Working_directory}/Workspace/SyntenyClustering/" -a "${anchorPoints}" -g "${gaps}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

# ==============================================================================
# Running Diamond analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Summarizing syntenet results in '${mode}' mode."

process_collinearity "${Working_directory}" "${mode}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

# ==============================================================================
# Running Clustering analysis
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Generating element clusters."

python3 ${auxiliary_path}/Clustering.py -i "${Working_directory}/Workspace/SyntenyClustering/Collinearity_percentage.txt" -o "${Working_directory}/Clusters/" -m "${minSize}" -n "${minNodes}" -t "${threshold}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished. Proceeding."

# ==============================================================================
# Running Clusters division
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 6: Sorting elements from clusters."

process_cluster_file "${Working_directory}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished synteny and clustering analysis. Data is stored in ${Working_directory}"
