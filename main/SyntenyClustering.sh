#!/usr/bin/env bash

# ==============================================================================
# Software check block
# ==============================================================================

source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Utils.sh"
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
fi

# ==============================================================================
# Function block
# ==============================================================================

# Function to print help message
function print_help() {
   echo -e "Script to run the syntenet pipeline using DIAMOND for sequence similarity search, summarizes the results, and organizes the output data.
   It executes six main steps:
   1. Executes the initial data preprocessing step required by the syntenet pipeline.
   2. Runs DIAMOND using the preprocessed data.
   3. Runs the interspecies synteny command of syntenet to identify regions with gene collinearity.
   4. Summarizes the results of syntenet on four possible modes:
     a. Raw: Return pairs that have a minimum of 8% of shared collinear genes. \033[01;31mWARNING\033[m: This mode may yield a high rate of false positives.
     b. SSP: Return only pairs with strong synteny.
     c. FilterBlast: Returns pairs that have been filtered using a BLAST-based approach at the nucleotide level.
     d. FilterMetric: Filter and update collinearity based on a metric system.
   5. Defines initial clusters of elements and performs a spectral clustering process to identify potential subclusters.
   6. Organizes the final data output for each identified cluster.
   "
   echo ""
   echo "Syntax: SAT $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -m <string> -a <integer> -g <integer> -n <integer> -s <integer> -t <integer> -th <float> -p <string> ]"
   echo ""
   echo "Required args:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored."
   echo ""
   echo "Required args with Default:"
   echo "-m, --mode: Specified the mode to run the summarizing process of synteny results (Default = FilterMetric) [Available mode: Raw, SSP, FilterBlast, FilterMetric]."
   echo "-a, --anchors: Number of minimum anchor points for syntenet to call a collinear region (Default = 8) [range: 5 - 25]."
   echo "-g, --gaps: Number of maximum allowed gaps between anchor points for syntenet to call a collinear region (Default = 8) [range: 5 - 25]."
   echo "-n, --minNodes: Minimum number of nodes in a cluster for spectral clustering to be attempted (Default: 4)."
   echo "-s, --minSize: The minimum desired size for any final sub-cluster (Default: 2)."
   echo "-t, --threads: Number of threads for searching software (DIAMOND and blast) (Default: 8)"
   echo "-th, --threshold: The minimum modularity score for a split to be accepted (Default: 0.05) [range: -0.5 - 1.0]."
   echo ""
   echo "Required args with Default in 'FilterBlast' mode:"
   echo "-fs, --fragmentSize: The minimum fragment size of a blast alignment to be used for blastn filter (Default: 2000) [range: 1000, 5000]."
   echo "-ms, --mergeSize: The minimum merge fragment size to be used for blastn filter (Default: 5000) [range: 2000, 10000]."
   echo "-i, --identity: The minimum percentage of identity of a blast alignment to be used for blastn filter (Default: 70.0) [range: 60.0, 90.0]"
   echo "-c, --coverage: The minimum coverage of the filter merge fragments for a pair to pass the filter (Default: 20.0) [range: 10.0, 50.0]"
   echo ""
   echo "Optional args:"
   echo "-help: Display this help message."
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


    if [[ -f "${data_dir}/metadata.csv" ]]; then
        cp "${data_dir}/metadata.csv" ${working_dir}
        metadata_flag=true
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

    Total_states=${#fasta_files[@]}
    State=0
    
    for fasta_file in "${fasta_files[@]}"
    do
        State=$(($State + 1))

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
        ProgressBar $State $Total_states
    done

    echo ""
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
    makeblastdb -dbtype nucl -parse_seqids -in "$temp_fasta" -out "$temp_db" &>/dev/null
    
    # Step 3: Run all-against-all search and append results to output file
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running blastn..."
    Total_states=$(grep -c ">" "$temp_fasta")
    State=0
    # Get all sequence headers from the temporary fasta file
    grep ">" "$temp_fasta" | sed 's/>//g' | while read -r line
    do
        State=$(($State + 1))
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
               >> "$output_file" 2>/dev/null
        ProgressBar $State $Total_states
    done

    echo ""
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
    local diamond_results_dir=$(find "$working_dir" -maxdepth 1 -type d -name "DiamondResults" 2>/dev/null)

    local temp_dir="${working_dir}/temp/"
    local out_path="$working_dir"

    if $metadata_flag; then
        local metadata_file="${working_dir}/metadata.csv"
    fi

    # Check that the subdirectories were found
    if [[ -z "$collinearity_path" || -z "$gff_path" || -z $nucleotide_dir || -z $diamond_results_dir ]]; then
        echo "Error: Required subdirectories (Collinearity, DiamondResults, Gff or Nucleotide) not found in '$working_dir'." >&2
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
    
    echo -e "element01;element02;General_percentage;element01_percentage;element02_percentage" > "${temp_dir}/Collinearity_percentage.txt"
    awk '{if($4>0){print}}' "${temp_prefix}_percentage_general.txt" > "${temp_prefix}_percentage_general_filter.txt"

    if [[ "${mode}" == "Raw" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';' >> "${temp_dir}/Collinearity_percentage.txt"
        # Here is the error.
        if $metadata_flag; then
            python ${auxiliary_path}/merge_metadata.py -d "${temp_dir}/Collinearity_percentage.txt" -m "$metadata_file" -o "${working_dir}/Collinearity_percentage.txt"
        else
            sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" >> "${working_dir}/Collinearity_percentage.txt"
        fi
    elif [[ "${mode}" == "SSP" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | \
            awk '{if(($3 >= 41) || (($4 >= 45 || $5 >= 45) && ($4/$5 >= 1.8 || $4/$5 <= 0.55 ))) print}' | \
            sed -e 's/ /;/g' | sort -t ';' >> "${temp_dir}/Collinearity_percentage.txt"
        if $metadata_flag; then
            python ${auxiliary_path}/merge_metadata.py -d "${temp_dir}/Collinearity_percentage.txt" -m "$metadata_file" -o "${working_dir}/Collinearity_percentage.txt"
        else
            sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" >> "${working_dir}/Collinearity_percentage.txt" 
        fi
    elif [[ "${mode}" == "FilterBlast" ]]
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
            python ${auxiliary_path}/Blast_CleanUp.py -f "${working_dir}/BlastnResults.out" -o "${working_dir}BlastnClean.out" -fs "${fragment_size}" -i "${identity}" -ms "${merge_size}" -c "$coverage"
        
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
        cat "${temp_prefix}_Collinearity_percentage_filter.txt" >> "${temp_dir}/Collinearity_percentage.txt"
        if $metadata_flag; then
            python ${auxiliary_path}/merge_metadata.py -d "${temp_dir}/Collinearity_percentage.txt"  -m "$metadata_file" -o "${working_dir}/Collinearity_percentage.txt"
        else
            sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" >> "${working_dir}/Collinearity_percentage.txt"
        fi
    elif [[ "${mode}" == "FilterMetric" ]]
    then
        # Generate files that will be filter.
        awk 'NR == FNR {f1[$1,$2] = $0; next} $1 SUBSEP $2 in f1 {print f1[$1,$2],"\t"$7,"\t"$8}' \
            "${temp_prefix}_percentage_general_filter.txt" "${temp_prefix}_percentage_pairwise.txt" | \
            awk '{print $1 FS $2 FS $4 FS $5 FS $6}' | awk '{if($3 >= 8) print}' | \
            sed -e 's/ /;/g' | sort -t ';'  >> "${temp_prefix}_Collinearity_percentage.txt"

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Filtering results..."

        Total_states=$(($(wc -l "${temp_prefix}_Collinearity_percentage.txt" | awk '{print $1}') - 1))
        State=0

        awk 'BEGIN{FS=";";OFS=" "}NR>1{file1=$1"_"$2;file2=$2"_"$1; print file1 OFS file2}' "${temp_prefix}_Collinearity_percentage.txt" | while read line
        do
            State=$(($State + 1))
            grep -E "[0-9]*-.*[0-9]*:" ${collinearity_path}/$(echo $line | awk '{print $1}').collinearity | awk 'BEGIN{FS=OFS="\t"}{print $2 OFS $3}' | sort -u > "${temp_prefix}_collinear1.txt"
            local Max_points=$(($(wc -l "${temp_prefix}_collinear1.txt" | awk '{print $1}') * 2))
            awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' "${temp_prefix}_collinear1.txt" > "${temp_prefix}_collinear2.txt"
            grep -f "${temp_prefix}_collinear1.txt" ${diamond_results_dir}/$(echo $line | awk '{print $1}').tsv | awk 'BEGIN{FS=OFS="\t"}{if($4>100 && $3>=60){print}}' > "${temp_prefix}_hits1.txt"
            grep -f "${temp_prefix}_collinear2.txt" ${diamond_results_dir}/$(echo $line | awk '{print $2}').tsv | awk 'BEGIN{FS=OFS="\t"}{if($4>100){print}}' > "${temp_prefix}_hits2.txt"
            local hits1=$(wc -l "${temp_prefix}_hits1.txt" | awk '{print $1}') 
            local hits2=$(wc -l "${temp_prefix}_hits2.txt" | awk '{print $1}')
            local bonus1=$(awk -F '\t' '{if($3>95){sum+= 0.1*($4/200)}}END{if(sum!=""){print sum}else{print 0}}' "${temp_prefix}_hits1.txt")
            local bonus2=$(awk -F '\t' '{if($3>95){sum+= 0.1*($4/200)}}END{if(sum!=""){print sum}else{print 0}}' "${temp_prefix}_hits2.txt")
            local Metric=$(echo "$hits1 + $hits2 + $bonus1 + $bonus2" | bc)
            if (( $(bc <<< "$Metric >= $anchorPoints") )); then
                local Metric_index=$(echo "print(min(round(${Metric}/${Max_points},2),1))" | python)
                echo $(echo $line | awk '{split($1,array,"_");print array[1]";"array[2]}')";"${Metric_index} >> "${working_dir}/Metrics_selected.out"
            fi
            ProgressBar $State $Total_states
        done
        echo ""

        join -t ';' <(sed -e 's/;/-/' "${temp_prefix}_Collinearity_percentage.txt" | sort) <(sed -e 's/;/-/' "${working_dir}/Metrics_selected.out" | sort) | awk 'BEGIN{FS=OFS=";"}{$2=$2*$5;$3=$3*$5;$4=$4*$5;print $1 OFS $2 OFS $3 OFS $4}' | sed 's/-/;/g' >> "${temp_dir}/Collinearity_percentage.txt"

        # Final stage
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating final report..."
        if $metadata_flag; then
            python ${auxiliary_path}/merge_metadata.py -d "${temp_dir}/Collinearity_percentage.txt"  -m "$metadata_file" -o "${working_dir}/Collinearity_percentage.txt"
        else
            sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" >> "${working_dir}/Collinearity_percentage.txt"
        fi
    fi

    # Remove temporary data
    rm -rf ${temp_prefix}*
}

process_cluster_file() {
    local base_dir="$1"
    local data_dir="${base_dir}/Data/"

    local cluster_dir="${base_dir}/Clusters/"
    local working_dir="${base_dir}/Workspace/SyntenyClustering/"
    
    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

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

    Total_states=$(wc -l "${cluster_path}" | awk '{print $1}')
    State=0

    while read -r cluster_id values_string; do
        State=$(($State + 1))

        IFS=' ' read -r -a values_array <<< "$values_string"

        #echo "Processing: $cluster_id with ${#values_array[@]} elements."
        
        mkdir -p "${cluster_dir}/${cluster_id}"
        mkdir -p "${cluster_dir}/${cluster_id}/Workspace"
        mkdir -p "${cluster_dir}/${cluster_id}/Data"
        mkdir -p "${cluster_dir}/${cluster_id}/Data/Exon"
        mkdir -p "${cluster_dir}/${cluster_id}/Data/Nucleotide"
        mkdir -p "${cluster_dir}/${cluster_id}/Data/Protein"
        mkdir -p "${cluster_dir}/${cluster_id}/Data/Gff"

        if $metadata_flag; then
            local updated_metadata="${cluster_dir}/${cluster_id}/Data/metadata.csv"
            head -n1 $metadata_file > $updated_metadata
        fi
        
        for value in "${values_array[@]}"; do
            cp ${gff_dir}/${value}.gff ${cluster_dir}/${cluster_id}/Data/Gff/
            cp ${protein_dir}/${value}.fa ${cluster_dir}/${cluster_id}/Data/Protein/
            cp ${nucleotide_dir}/${value}.fa ${cluster_dir}/${cluster_id}/Data/Nucleotide/    
            cp ${exon_dir}/${value}.fa ${cluster_dir}/${cluster_id}/Data/Exon/
            if $metadata_flag; then
                grep -w $value $metadata_file >> $updated_metadata
                fi        
        done

        ProgressBar $State $Total_states
    done < "$cluster_path"

    # Restore the original IFS at the end of the function.
    IFS="$OLD_IFS"
    echo ""
    echo "Processing complete."
}

# ==============================================================================
# Variables block
# ==============================================================================


# Initialize variables
Working_directory=""
mode="FilterMetric"
anchorPoints="8"
gaps="8"
minNodes="4"
minSize="1"
threshold="0.05"
fragment_size="2000"
merge_size="5000"
identity="70.0"
coverage="20.0"
threads="8"
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
        -th|--threshold)
            shift
            threshold="$1"
            ;;
        -fs|--fragmentSize)
            shift
            fragment_size="$1"
            ;;
        -ms|--mergeSize)
            shift
            merge_size="$1"
            ;;
        -i|--identity)
            shift
            identity="$1"
            ;;
        -c|--coverage)
            shift
            coverage="$1"
            ;;
        -t|--threads)
            shift
            threads="$1"
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

# ==============================================================================
# Script start block
# ==============================================================================

echo "Running $(basename -s .sh "$0" ) command under the following parameters:"
echo "  Working directory: " "$Working_directory"
echo "  Mode: " "$mode"
echo "  Minimum number of anchor points: " "$anchorPoints"
echo "  Maximum number of gaps: " "$gaps"
echo "  Minimum number of nodes for Spectral clustering: " "$minNodes"
echo "  Minimum size of sub-cluster: " "$minSize"
echo "  Modularity score threshold for Spectral clustering: " "$threshold"
echo "  Number of threads: " "$threads"
echo ""

# ==============================================================================
# Check variables block
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking arguments and input files..."

# Check for mandatory argument and define the path as absolute
if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check working directory exists
if [[ ! -d "$Working_directory" ]]; then
    echo "Error: Directory '$Working_directory' does not exist."
    exit 1
else
    Working_directory=$(realpath $Working_directory)
fi
check_directory_structure "${Working_directory}"

# Check if mode parameter is correct
check_mode_parameter "$mode" "$(basename -s .sh "$0" )"

# Check for auxiliary scripts
check_auxiliary_scripts "$auxiliary_path" "$(basename -s .sh "$0" )"

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

# Check parameters for Blastn filtering
if [[ $mode == "FilterBlast" ]]; then
    if [[ "$fragment_size" =~ ^[0-9]+$ ]]; then
        if (( fragment_size < 1000 || fragment_size > 5000 )); then
            echo "Error: '$fragment_size' is not an accepted value for fragment size."
            print_help
            exit 1
        fi
    else
        echo "Error: '$fragment_size' is not a positive integer."
        print_help
        exit 1
    fi
    if [[ "$merge_size" =~ ^[0-9]+$ ]]; then
        if (( merge_size < 2000 || gaps > 10000 )); then
            echo "Error: '$merge_size' is not an accepted value for merge size."
            print_help
            exit 1
        fi
    else
        echo "Error: '$merge_size' is not a positive integer."
        print_help
        exit 1
    fi
    if [[ "$fragment_size" -gt "$merge_size" ]]; then
        echo "Error: fragment size cannot be greater than merge size."
        print_help
        exit 1
    fi
    if [[ "$identity" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
        if (( $(echo "$identity < 60.0" | bc -l) )) && (( $(echo "$identity > 90.0" | bc -l) )); then
            echo "Error: '$identity' is not an accepted value for clustering modularity score."
            print_help
            exit 1
        fi
    else
        echo "Error: '$identity' is not a float."
        print_help
        exit 1
    fi
    if [[ "$coverage" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
        if (( $(echo "$coverage < 10.0" | bc -l) )) && (( $(echo "$coverage > 50.0" | bc -l) )); then
            echo "Error: '$coverage' is not an accepted value for clustering modularity score."
            print_help
            exit 1
        fi
    else
        echo "Error: '$coverage' is not a float."
        print_help
        exit 1
    fi
fi

# Check thread parameter
check_threads "$threads"

# Check for required software
check_required_software "$(basename -s .sh "$0" )"

# ==============================================================================
# Main Block
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"
organize_working_directory "${Working_directory}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Preprocessing data for Diamond."
Rscript ${auxiliary_path}/syntenetPreprocess.R "${Working_directory}/Workspace/SyntenyClustering/"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running Diamond analysis."
process_diamond "${Working_directory}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running syntenet analysis."
Rscript ${auxiliary_path}/syntenetAnalysis.R  -d "${Working_directory}/Workspace/SyntenyClustering/" -a "${anchorPoints}" -g "${gaps}" -t "${threads}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Summarizing syntenet results in '${mode}' mode."
process_collinearity "${Working_directory}" "${mode}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Generating element clusters."
python ${auxiliary_path}/Clustering.py -i "${Working_directory}/Workspace/SyntenyClustering/Collinearity_percentage.txt" -o "${Working_directory}/Clusters/" -m "${minSize}" -n "${minNodes}" -t "${threshold}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 6: Sorting elements from clusters."
process_cluster_file "${Working_directory}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished synteny and clustering analysis."
