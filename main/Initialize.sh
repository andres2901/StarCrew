#!/usr/bin/env bash

# ==============================================================================
# SOFTWARE CHECK AND ENVIRONMENT SETUP
# ==============================================================================

# Define library and auxiliary paths
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Utils.sh"
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
    check_auxiliary_scripts "$auxiliary_path" "$(basename -s .sh "$0" )"
fi

# Validate required software for the current script
check_required_software "$(basename -s .sh "$0" )"

# Validate agat config file
agat_config="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../agat_config.yaml"
if [[ ! -f "$agat_config" ]]; then
    echo "Error: file '$agat_config' does not exist."
    exit 1
else
    agat_config=$(realpath $agat_config)
fi

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

# Function to print help message
function print_help() {
   echo "Script to organize the working directory to run the subsequent commands in the workflow."
   echo ""
   echo "Syntax: StarCrew $(basename -s .sh "$0" ) [ -help ] -f <filte_path> -g <file_path> [ -m <string> -gc <integer> -r <integer> -mg <integer> -o <string>  { -b <file_path> -s <character> -c <file_path> } -M <file_path> -t <integer> --overwrite ]"
   echo ""
   echo "Required args:"
   echo "-f, --fasta:  multifasta file wih the elements to study."
   echo ""
   echo "Required args with Default:"
   echo "-m, --mode: Mode of the input to initialize (Defaul = Simple) [Available mode: Simple, Starfish]."
   echo "-gc, --gc: integer value of gc content to filter out elements with too low gc content (Default = 0) [range: 20 - 45]."
   echo "-r, --rip: integer value of the minimum coverage of the element to be possibly affected by RIP to be filter out (Default = 0) [range: 30 - 80]"
   echo "-mg, --minGene: Minimum number of genes in an element to be include in the dataset (Default: 8) [range: 5 - 100]"
   echo "-o, --outDirectory: Specify working directory  name (Default = WorkingDirectory)."
   echo ""
   echo "Required args in 'Starfish' mode:"
   echo "-g, --gff: 2 column tsv: genome code, path to GFF. The path should be to the original gff files and not the ones formatted to run starfish."
   echo "-b, --boundaries: *.elements.feat file output of 'starfish summary' command."
   echo "-c, --captains: *_tyr.filt_intersect.fas file output of 'starfish annotate' command."
   echo ""
   echo "Required args in 'Starfish' mode with Default:"
   echo "-s, --separator: character separating genomeID from featureID that was used for Starfish run (Default = '_')."
   echo ""
   echo "Required args in 'Simple' mode:"
   echo "-g, --gff: Path to the GFF file containing gene predictions with element-relative coordinates for all elements in the fasta file."
   echo ""
   echo "Optional args:"
   echo "-M, --Metadata: csv file delimited by semicolon with the metadata information (Check wrapper documentation for more information)."
   echo "-t, --threads: threads for MetaEuk in 'Starfish' mode (Default: 24)"
   echo "--overwrite: Flag to overwrite in case there is already a previous run (Default: off)"
   echo "-help: Display this help message."
}

check_duplicates() {
    local fasta_path="$1"
    local working_dir="$2"

    seqkit rmdup -n "$fasta_path" -o /dev/null -d ${working_dir}/temp/temp_duplicated_headers &> /dev/null

    if [ -s ${working_dir}/temp/temp_duplicated_headers ]; then
        echo "Error: Duplicated headers found. Script terminated." >&2
        echo "Duplicate headers and their counts:" >&2
        echo "$duplicated_headers" >&2
        exit 1
    fi
}

# Function to filter input based on RIP-like signal and gc content
Filter_input() {
    local fasta_path="$1"
    local working_dir="$2"
    local filter_option="$3"

    local temp_dir="${working_dir}/temp/"

    if [[ "$filter_option" == 1 ]]
    then
       echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on RIP-like signal"
       python ${auxiliary_path}/rip_calculator.py "$fasta_path" -tc 0.01 -tp 1 -ts 1 -w 500 -s 100 | awk -F '\t' -v min="$rip" 'NR>1{if($4 > min){print $1}}' > ${working_dir}/Elements_filterRIPlike.txt
       removed_elements=$(wc -l ${working_dir}/Elements_filterRIPlike.txt | awk '{print $1}')
       echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove based on RIP-like signal"
       seqkit grep --quiet -n -v -f ${working_dir}/Elements_filterRIPlike.txt "$fasta_path" > ${temp_dir}/Sequences.fa
    elif [[ "$filter_option" == 2 ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on gc content"
        seqkit fx2tab -g -n "$fasta_path" | awk -v min="$filter" '{if($NF < min){$NF=""; print $0}}' | sed -e 's/ $//g' > ${working_dir}/Elements_filterGC.txt
        removed_elements=$(wc -l ${working_dir}/Elements_filterGC.txt | awk '{print $1}')
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove based on gc content"
        seqkit grep --quiet -n -v -f ${working_dir}/Elements_filterGC.txt "$fasta_path" > ${temp_dir}/Sequences-filter1.fa

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on RIP-like signal"
        python ${auxiliary_path}/rip_calculator.py "${temp_dir}/Sequences-filter1.fa" -tc 0.01 -tp 1 -ts 1 -w 500 -s 100 | awk -F '\t' -v min="$rip" 'NR>1{if($4 > min){print $1}}' > ${working_dir}/Elements_filterRIPlike.txt
        removed_elements=$(wc -l ${working_dir}/Elements_filterRIPlike.txt | awk '{print $1}')
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove based on RIP-like signal"
        seqkit grep --quiet -n -v -f ${working_dir}/Elements_filterRIPlike.txt "${temp_dir}/Sequences-filter1.fa" > ${temp_dir}/Sequences.fa
    elif [[ "$filter_option == 3" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on gc content"
        seqkit fx2tab -g -n "$fasta_path" | awk -v min="$filter" '{if($NF < min){$NF=""; print $0}}' | sed -e 's/ $//g' > ${working_dir}/Elements_filterGC.txt
        removed_elements=$(wc -l ${working_dir}/Elements_filterGC.txt | awk '{print $1}')
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove base on gc content."
        seqkit grep --quiet -n -v -f ${working_dir}/Elements_filterGC.txt "$fasta_path" > ${temp_dir}/Sequences.fa
    fi
}

# Function to convert a number to a base-62 string for unique IDs
base62() {
    local n=$1
    local base62_chars="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    local result=""
    if [ "$n" -eq 0 ]; then
        echo "0"
        return
    fi
    while [ "$n" -gt 0 ]; do
        local rem=$((n % 62))
        local char=${base62_chars:$rem:1}
        result="$char$result"
        n=$((n / 62))
    done
    printf '%s\n' "$result"
}

# Function to process headers to avoid issues during syntenet analysis
Header_processing() {
    local fasta_path="$1"
    local working_dir="$2"

    local output_fasta="${working_dir}/Sequences.fa"
    local association_csv="${working_dir}/metadata_files/sequence_head.csv"

    local max_initial_count=25
    local ID_COUNTER=0

    # Create a temporary file with initials and original headers
    seqkit fx2tab --name "$fasta_path" | awk -F'\t' '
    {
        original_header = $1
        cleaned_header = original_header
        gsub(/[^a-zA-Z0-9]/, "", cleaned_header)
        initials = substr(cleaned_header, 1, 5)
        print initials, original_header
    }
    ' > ${working_dir}/temp/temp_initials_and_headers.tsv

    # Count the occurrences of each 5-char initial group
    awk '{ print $1 }' ${working_dir}/temp/temp_initials_and_headers.tsv | sort | uniq -c > ${working_dir}/temp/temp_initial_counts.tsv

    # Reset output files to ensure a clean start
    > "$output_fasta"
    echo "new_header;original_header" > "$association_csv"

    # Create a temporary file to track used IDs to prevent collisions
    temp_used_ids="temp_used_ids.txt"
    > "$temp_used_ids"

    Total_states=$(wc -l ${working_dir}/temp/temp_initials_and_headers.tsv | awk '{print $1}')
    State=0

    # Loop through all headers and apply the renaming logic
    while read -r initials original_header; do
        State=$(($State + 1))
 
        # Get the count for the current initial group
        group_count=$(grep -w "^[[:space:]]*[0-9]*[[:space:]]*$initials$" ${working_dir}/temp/temp_initial_counts.tsv | awk '{print $1}')
 
        # Clean the header for processing
        cleaned_header=$(echo "$original_header" | sed 's/[^a-zA-Z0-9]//g')
 
        new_header=""
 
        # Check if this header belongs to an over-threshold group
        if [ "$group_count" -gt "$max_initial_count" ]; then
            # Check if the header's own last 5 characters can be used
            current_last5=$(echo "$cleaned_header" | tail -c 6)
        
            # Check if the last 5 chars are a valid length and not already used
            if [ "${#current_last5}" -eq 5 ] && ! grep -q "^$current_last5$" "$temp_used_ids"; then
                new_header="$current_last5"
                echo "$new_header" >> "$temp_used_ids" # <-- TRACK THE LAST5 ID IMMEDIATELY
            else
                # Fallback to generating a unique ID and ENSURE it is new
                while true; do
                    UNIQUE_ID=$(base62 "$ID_COUNTER")
                    # Pad to 5 characters with '0'
                    PADDED_ID=$(printf "%05s" "$UNIQUE_ID" | sed 's/ /0/g')

                    # **CRITICAL FIX: Check if the generated ID is already used**
                    if ! grep -q "^$PADDED_ID$" "$temp_used_ids"; then
                        new_header="$PADDED_ID"
                        echo "$new_header" >> "$temp_used_ids" # <-- TRACK THE BASE62 ID IMMEDIATELY
                        ID_COUNTER=$((ID_COUNTER + 1))
                        break # Exit the while true loop
                    fi
                    ID_COUNTER=$((ID_COUNTER + 1)) # Increment if collision found and try again
                done
            fi
        else
            if [ "${#cleaned_header}" -gt 30 ]; then
                while true; do
                    UNIQUE_ID=$(base62 "$ID_COUNTER")
                    # Pad to 5 characters with '0'
                    PADDED_ID=$(printf "%05s" "$UNIQUE_ID" | sed 's/ /0/g')

                    # **CRITICAL FIX: Check if the generated ID is already used**
                    if ! grep -q "^$PADDED_ID$" "$temp_used_ids"; then
                        new_header="$PADDED_ID"
                        echo "$new_header" >> "$temp_used_ids" # <-- TRACK THE BASE62 ID IMMEDIATELY
                        ID_COUNTER=$((ID_COUNTER + 1))
                        break # Exit the while true loop
                    fi
                    ID_COUNTER=$((ID_COUNTER + 1)) # Increment if collision found and try again
                done
            else
                # The group is within the limit, use the full cleaned header
                new_header="$cleaned_header"
            fi
            
            # Only track if it's a 5-char header, to prevent collisions with
            # the last-5-char logic later on.
            if [ "${#new_header}" -eq 5 ]; then
                echo "$new_header" >> "$temp_used_ids"
            fi
        fi
 
        # Print to the temporary file use for renamed
        echo -e "$new_header\t$original_header" >> ${working_dir}/temp/temp_association.tsv

        ProgressBar $State $Total_states
    done < ${working_dir}/temp/temp_initials_and_headers.tsv
    echo ""

    # Use seqkit to renames headers
    echo "    [$(date "+%Y-%m-%d %H:%M:%S")] Updating headers..."
    awk 'BEGIN {OFS="\t"} {print $2, $1}' ${working_dir}/temp/temp_association.tsv > ${working_dir}/temp/temp_association2.tsv
    seqkit replace --kv-file ${working_dir}/temp/temp_association2.tsv -p "(.*)" -r "{kv}" "$fasta_path" | seqkit seq --quiet -u > "$output_fasta"

    # Create the final CSV association file
    sed 's/\t/;/g' ${working_dir}/temp/temp_association.tsv | sed '1s/^/new_header;original_header\n/' > "$association_csv"

    # Create the updated metadata file
    if $metadata; then
        local metadata_csv="${working_dir}/metadata_files/metadata.csv"
        head -n1 $metadata_path > $metadata_csv
        sed -i -e 's/ElementID/ElementID;ElementID_updated/' $metadata_csv
        awk 'BEGIN{FS=OFS=";"}{for (i = 1; i <= NF; i++) {if($i=="")$i="NA"}; print}' $metadata_path > ${working_dir}/temp/temp_metadata.csv
        join -1 2 -2 1 -t ';' <( sort -t ";" -k2,2 $association_csv) <(sort -t ";" -k1,1 ${working_dir}/temp/temp_metadata.csv) >> $metadata_csv
    else
        echo "    [$(date "+%Y-%m-%d %H:%M:%S")] Skipping metadata file update..."
    fi
}

# Function to calculate statastics of gene content in elements
gene_stats() {
    local working_dir="$1"

    echo -e "Starship""\t""Number_genes""\t""Avg_gene_length""\t""Avg_intergenic_length" > ${working_dir}/Gene_stats.txt
    ls ${working_dir}/Data/Gff/ | xargs -n 1 basename -s .gff | while read line
    do 
        echo -e ${line}"\t"$(grep -c -P "\tgene\t" ${working_dir}/Data/Gff/${line}.gff)"\t"$(grep -P "\tgene\t" ${working_dir}/Data/Gff/${line}.gff | awk '{ sum  += $5 - $4 } END { if(NR > 0) {print sum / NR} else {print $0}}')"\t"$(grep -P "\tgene\t" ${working_dir}/Data/Gff/${line}.gff | sort -k4 -n | awk 'NR==1 {prev_col2 = $5; next} {diff = $4 - prev_col2; prev_col2 = $5; if (diff > 0) total_sum += diff} END {if((NR - 1) > 0) {print total_sum / (NR - 1)} else {print 0}}') >> ${working_dir}/Gene_stats.txt 
    done
}

organize_info() {
    local working_dir="$1"
    local CDS_flag="$2"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing gene statistic and removing elements with low gene models..."
    gene_stats "${working_dir}"

    awk -v min="$minimum_gene_content" 'NR>1{if($2>=min){print $1}}' ${working_dir}/Gene_stats.txt | sed $'s/[^[:print:]\t]//g' > ${working_dir}/temp/Good_elements.txt
    ls ${working_dir}/Data/Gff/ | xargs -n 1 basename -s .gff > ${working_dir}/temp/All_elements.txt

    grep -v -f ${working_dir}/temp/Good_elements.txt ${working_dir}/temp/All_elements.txt | while read line
    do
        rm ${working_dir}/Data/Gff/${line}.gff
    done

    echo -e "  [$(date "+%Y-%m-%d %H:%M:%S")] Organizing info..."

    Total_states=$(ls ${working_dir}/Data/Gff/ | xargs -n 1 basename -s .gff | wc -l)
    echo -e "  [$(date "+%Y-%m-%d %H:%M:%S")] Processing '$Total_states' elements that have at least '${minimum_gene_content}' genes."
    State=0

    ls ${working_dir}/Data/Gff/ | xargs -n 1 basename -s .gff | while read line 
    do
        State=$(($State + 1))
        # Dividing nucleotide sequence of elements
        echo $line > ${working_dir}/temp/temp_element.txt
        seqkit grep -n -f ${working_dir}/temp/temp_element.txt ${working_dir}/Sequences.fa -o ${working_dir}/Data/Nucleotide/${line}.fa &> /dev/null

        # Creating CDS files
        if $CDS_flag; then
            agat_sp_extract_sequences.pl --config ${agat_config} --gff ${working_dir}/Data/Gff/${line}.gff --fasta ${working_dir}/Data/Nucleotide/${line}.fa -t cds -o ${working_dir}/temp/CDS/${line}.fa &> /dev/null
        else
            agat_sp_extract_sequences.pl --config ${agat_config} --gff ${working_dir}/Data/Gff/${line}.gff --fasta ${working_dir}/Data/Nucleotide/${line}.fa -t exon --merge -o ${working_dir}/temp/CDS/${line}.fa &> /dev/null
        fi

        awk '{if($2){$1=">"$2} print $1}' ${working_dir}/temp/CDS/${line}.fa | sed 's/gene=//g' > ${working_dir}/Data/CDS/${line}.fa

        seqkit translate -f 1 ${working_dir}/Data/CDS/${line}.fa > ${working_dir}/Data/Protein/${line}.fa

        rm ${working_dir}/Data/Nucleotide/${line}.fa.index* &> /dev/null
        ProgressBar $State $Total_states
    done
    echo ""
}

# Function to process input in 'Simple' mode
Process_simple() {
    local working_dir="$1"
    local gff_file="$2"
    local CDS_flag=false

    mkdir -p ${working_dir}/Data/Protein/ ${working_dir}/Data/Gff/ ${working_dir}/Data/Nucleotide/ ${working_dir}/Data/CDS/ ${working_dir}/temp/CDS/ ${working_dir}/temp/gff/

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Preparing gff file..."

    if [[ $(grep -w -i -c "cds" ${gff_file}) -ge 10 ]]; then
        CDS_flag=true
    fi

    local Total_states=$(wc -l ${working_dir}/temp/temp_association.tsv | awk '{print $1}')
    local State=0

    cat ${working_dir}/temp/temp_association.tsv | while read line
    do
        State=$(($State + 1))     
        A=$(echo $line | awk '{print $2}')
        B=$(echo $line | awk '{print $1}')
        grep -w "^$A" ${gff_file} | sed s/$A/$B/ >> ${working_dir}/temp/gff_file.gff
        ProgressBar $State $Total_states
    done

    gff_file="${working_dir}/temp/gff_file.gff"

    echo -e "\n  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing gff file per element..."
    local Total_states=$(grep -v "^#" ${gff_file} | awk '{print $1}' | sort | uniq | wc -l)
    local State=0

    grep -v "^#" ${gff_file} | awk '{print $1}' | sort | uniq | while read line
    do
        State=$(($State + 1))
        grep "^#" ${gff_file} > ${working_dir}/temp/gff/${line}.gff
        grep -w "^$line" ${gff_file} >> ${working_dir}/temp/gff/${line}.gff
        agat_sp_keep_longest_isoform.pl --config ${agat_config} --gff ${working_dir}/temp/gff/${line}.gff -o ${working_dir}/Data/Gff/${line}.gff &> /dev/null
        sed -i -e "s/ID=/ID=${line}\./g" -e "s/Parent=/Parent=${line}\./g" ${working_dir}/Data/Gff/${line}.gff &> /dev/null
        ProgressBar $State $Total_states
    done
    echo ""

    organize_info "$working_dir" "$CDS_flag"
}

# Function to process input in 'Starfish' mode
Process_starfish() {
    local working_dir="$1"
    local boundaries_path="$2"
    local gff_path="$3"
    local CDS_flag=false

    mkdir -p ${working_dir}/Data/Protein/ ${working_dir}/Data/Gff/ ${working_dir}/Data/Nucleotide/ ${working_dir}/Data/CDS/ ${working_dir}/temp/CDS/ ${working_dir}/temp/gff/ ${working_dir}/temp/gff2/

    Total_states=$(wc -l ${working_dir}/temp/temp_association.tsv | awk '{print $1}')
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] processing '$Total_states' elements."

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating coordinate file..."
    State=0

    awk -v s="$separator" 'BEGIN{FS=OFS="\t"}NR>1{split($1,array,s);print $2 OFS array[2] OFS $4 OFS $5 OFS $7 OFS array[1]}' ${boundaries_path} > ${working_dir}/temp/coordinate_file.txt
    cat ${working_dir}/temp/temp_association.tsv | while read line
    do
        State=$(($State + 1))     
        A=$(echo $line | awk '{print $2}' | awk -F '|' '{print $1}')
        B=$(echo $line | awk '{print $1}')
        grep "^$A" ${working_dir}/temp/coordinate_file.txt | sed s/$A/$B/ >> ${working_dir}/Coordinate_file.txt
        ProgressBar $State $Total_states
    done
    
    echo -e "\n  [$(date "+%Y-%m-%d %H:%M:%S")] Preparing gff files..." 
    Total_states=$(wc -l ${gff_path} | awk '{print $1}')
    State=0

    cat ${gff_path} | while read line
    do
        State=$(($State + 1))
        GFF_FILE=$(grep -w "${line}" ${gff_path} | awk '{print $2}')
        check_gff_file "$(echo $line | awk '{print $2}')"
        sed -e s/${separator}//g $(echo $line | awk '{print $2}') > ${working_dir}/temp/gff/$(echo $line | awk '{print $1}').gff
        ProgressBar $State $Total_states
    done

    if [[ $(grep -w -i -c "cds" ${working_dir}/temp/gff/$(ls ${working_dir}/temp/gff/ | head -n1)) -ge 1 ]]; then
        CDS_flag=true
    fi

    echo -e "\n  [$(date "+%Y-%m-%d %H:%M:%S")] Extracting gene model information per element..."
    Total_states=$(awk '{print $NF}' ${working_dir}/Coordinate_file.txt | sort -u | wc -l)
    State=0
    
    awk '{print $NF}' ${working_dir}/Coordinate_file.txt | sort -u | while read line 
    do
        State=$(($State + 1))
        grep -w $line ${working_dir}/Coordinate_file.txt | cut -d$'\t' -f 1-5 > ${working_dir}/temp/temp_coordinate_file.txt
        python ${auxiliary_path}/gff_slicer.py -c ${working_dir}/temp/temp_coordinate_file.txt -i ${working_dir}/temp/gff/${line}.gff -o ${working_dir}/temp/gff2/ &> /dev/null
        awk '{print $1}' ${working_dir}/temp/temp_coordinate_file.txt | while read line; do
            agat_sp_keep_longest_isoform.pl --config ${agat_config} --gff ${working_dir}/temp/gff2/${line}.gff -o ${working_dir}/Data/Gff/${line}.gff &> /dev/null
        done
        ProgressBar $State $Total_states
    done

    rm ${working_dir}/temp/gff2/*.gff

    echo -e "\n  [$(date "+%Y-%m-%d %H:%M:%S")] Updating elements gff file..."
    Total_states=$(ls ${working_dir}/Data/Gff/ | xargs -n 1 basename -s .gff | wc -l)
    State=0
    
    ls ${working_dir}/Data/Gff/ | xargs -n 1 basename -s .gff | while read line 
    do
        State=$(($State + 1))
        sed -i -e "s/ID=/ID=${line}\./g" -e "s/Parent=/Parent=${line}\./g" ${working_dir}/Data/Gff/${line}.gff
        
        ProgressBar $State $Total_states
    done

    echo -e "\n  [$(date "+%Y-%m-%d %H:%M:%S")] Performing Captain identification based on Starfish output captain database using metaeuk..."

    metaeuk createdb ${working_dir}/Sequences.fa ${working_dir}/temp/ContigsDB --dbtype 2 -v 0 &> /dev/null
    metaeuk createdb $captains_path ${working_dir}/temp/ProteinDB --dbtype 1 -v 0 &> /dev/null

    metaeuk predictexons ${working_dir}/temp/ContigsDB ${working_dir}/temp/ProteinDB  ${working_dir}/temp/metaeukResults ${working_dir}/temp/tempFolder -s 7.5 --exhaustive-search 1 --orf-start-mode 0 --min-seq-id 0.95 --metaeuk-tcov 0.75 --min-length 200 --remove-tmp-files 1 --use-all-table-starts 1 --threads ${threads} --disk-space-limit 100G &> /dev/null
    metaeuk reduceredundancy ${working_dir}/temp/metaeukResults ${working_dir}/temp/metaeukpred ${working_dir}/temp/metaeukgroups --threads ${threads} -v 0 &> /dev/null
    metaeuk unitesetstofasta ${working_dir}/temp/ContigsDB ${working_dir}/temp/ProteinDB ${working_dir}/temp/metaeukpred ${working_dir}/temp/metaeukFinal --threads ${threads} -v 0 &> /dev/null

    sed -i -e 's/Target_ID=.*;TCS_//g' ${working_dir}/temp/metaeukFinal.gff
    sed -i -e 's/exon/CDS/g' ${working_dir}/temp/metaeukFinal.gff
    agat_convert_sp_gxf2gxf.pl --config ${agat_config} --gff ${working_dir}/temp/metaeukFinal.gff -o ${working_dir}/temp/metaeuk.gff &> /dev/null

    agat_sp_filter_by_ORF_size.pl --config ${agat_config} --gff ${working_dir}/temp/metaeuk.gff -s 200 -o ${working_dir}/temp/metaeuk_ORF.gff &> /dev/null


    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Merging gff files..."
    Total_states=$(grep -v "#" ${working_dir}/temp/metaeuk_ORF_sup200.gff | awk '{print $1}' | sort -u | wc -l)
    State=0

    rm ${working_dir}/temp/gff/*
    mkdir ${working_dir}/temp/modelsKeep/

    grep -v "#" ${working_dir}/temp/metaeuk_ORF_sup200.gff | awk '{print $1}' | sort -u | while read line 
    do
        State=$(($State + 1))

        echo "#gff version-3" > ${working_dir}/temp/gff/temp_captain_${line}.gff

        grep -w "^${line}" ${working_dir}/temp/metaeuk_ORF_sup200.gff >> ${working_dir}/temp/gff/temp_captain_${line}.gff
        agat_sp_manage_IDs.pl --config ${agat_config} --gff ${working_dir}/temp/gff/temp_captain_${line}.gff --prefix ${line}.metaeuk -o ${working_dir}/temp/gff/temp_captain2_${line}.gff &> /dev/null
        cp ${working_dir}/Data/Gff/${line}.gff ${working_dir}/temp/gff/${line}_merge.gff
        grep -v "^#" ${working_dir}/temp/gff/temp_captain2_${line}.gff >> ${working_dir}/temp/gff/${line}_merge.gff 
        python ${auxiliary_path}/merge.py ${working_dir}/temp/gff/${line}_merge.gff ${working_dir}/temp/modelsKeep/${line}.txt &> /dev/null
        agat_sp_filter_feature_from_keep_list.pl --config ${agat_config} --gff ${working_dir}/temp/gff/${line}_merge.gff  --keep_list ${working_dir}/temp/modelsKeep/${line}.txt --output ${working_dir}/temp/gff/${line}_mergeKeep.gff &> /dev/null
        rm ${working_dir}/Data/Gff/${line}.gff
        agat_sp_keep_longest_isoform.pl --config ${agat_config} --gff ${working_dir}/temp/gff/${line}_mergeKeep.gff -o ${working_dir}/Data/Gff/${line}.gff &> /dev/null

        ProgressBar $State $Total_states
    done
    echo ""

    organize_info "$working_dir" "$CDS_flag"
}

# ==============================================================================
# VARIABLES AND ARGUMENT PARSING
# ==============================================================================

fasta_path=""
mode="Simple"
metadata_path=""
filter="0"
rip="100"
minimum_gene_content="8"
out_directory="WorkingDirectory"
gff_path=""
boundaries_path=""
captains_path=""
separator="_"
threads="24"
overwrite=false
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -f|--fasta) 
	        shift
	        fasta_path="$1"
	        ;;
        -m|--mode)
            shift
            mode="$1"
            ;;
        -M|--Metadata)
            shift
            metadata_path="$1"
            ;;
        -gc|--gc)
            shift
            filter="$1"
            ;;
        -r|--rip)
            shift
            rip="$1"
            ;;
        -mg|--minGene)
            shift
            minimum_gene_content="$1"
            ;;
        -o|--outDirectory)
            shift
            out_directory="$1"
            ;;
        -g|--gff)
            shift
            gff_path="$1"
            ;;
        -b|--boundaries)
            shift
            boundaries_path="$1"
            ;;
        -c|--captains)
            shift
            captains_path="$1"
            ;;
        -s|--separator)
            shift
            separator="$1"
            ;; 
        -t|--threads)
            shift
            threads="$1"
            ;;
        --overwrite)
            overwrite=true
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

if $help_flag; then
    print_help
    exit 0
fi

echo "Running $(basename -s .sh "$0" ) command under the following parameters:"
echo "  Fasta file: " "$fasta_path"
echo "  Mode: " "$mode"
echo "  Minimum %gc: " "$filter"
echo "  Maximum %RIP-like signal: " "$rip"
echo "  Minumum gene content: " "$minimum_gene_content"
echo "  Output directory name: " "$out_directory"
echo "  Gff file: " "$gff_path"
echo "  Starfish boundary file: " "$boundaries_path"
echo "  Starfish captain file: " "$captains_path"
echo "  Separator used during starfish run: " "$separator"
echo "  Metadata file: " "$metadata_path" 
echo "  Threads: " "$threads"
echo "  Overwrite previous run: " "$overwrite"
echo ""

# ==============================================================================
# ARGUMENTS AND INPUT CHECK
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking arguments and input files..."

check_mode_parameter "$mode" "$(basename -s .sh "$0" )"

if [[ "$mode" == "Simple" ]]; then
    if [[ -z "$fasta_path" || -z "$gff_path" ]]; then
        echo "Error: Missing required argument(s)."
        print_help
        exit 1
    else
        check_fasta_dna "$fasta_path"
        fasta_path=$(realpath $fasta_path)

        check_gff_file "$gff_path"
        gff_path=$(realpath $gff_path)
    fi
elif [[ "$mode" == "Starfish" ]]; then
    if [[ -z "$fasta_path" || -z "$gff_path" || -z "$boundaries_path" || -z "$captains_path" || -z "$separator" ]]; then
        echo "Error: Missing required argument(s)."
        print_help
        exit 1
    else
        check_fasta_dna "$fasta_path"
        fasta_path=$(realpath $fasta_path)

        check_gff_paths "$gff_path"
        gff_path=$(realpath $gff_path)

        check_boundaries_file "$boundaries_path"
        boundaries_path=$(realpath $boundaries_path)

        check_fasta_protein "$captains_path"
        captains_path=$(realpath $captains_path)

        impossible_separator=":;|"
        impossible_separator_pattern="[${impossible_separator}]"
        if [[ "${separator}" == "${impossible_separator_pattern}" ]]; then
            echo "Error: '$separator' is not accepted."
            exit 1
        fi
    fi
fi

if [[ -z "$metadata_path" ]]; then
    metadata=false
else
    check_metadata_file "$metadata_path" "$fasta_path"
    metadata=true
    metadata_path=$(realpath $metadata_path)
fi

if [[ "$minimum_gene_content" =~ ^[0-9]+$ ]]; then
    if (( $minimum_gene_content < 5 || $minimum_gene_content > 100 )); then
        echo "Error: '$minimum_gene_content' minimum gene content is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$minimum_gene_content' is not a positive integer."
    print_help
    exit 1
fi

if [[ "$filter" =~ ^[0-9]+$ ]]; then
    if (( $filter != 0 && ($filter < 20 || $filter > 45) )); then
        echo "Error: '$filter' gc content is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$filter' is not a positive integer."
    print_help
    exit 1
fi

if [[ "$rip" =~ ^[0-9]+$ ]]; then
    if (( $rip != 100 && ($rip < 30 || $rip > 80) )); then
        echo "Error: '$rip' rip coverage is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$rip' is not a positive integer."
    print_help
    exit 1
fi

check_threads "${threads}"

if [[ -d "$out_directory" ]]; then
    if $overwrite; then
        rm -r $out_directory
        mkdir $out_directory
        mkdir -p ${out_directory}/Data ${out_directory}/Workspace ${out_directory}/metadata_files ${out_directory}/temp/ 
    else
        echo "Error: directory '$out_directory' already exist."
        echo "If you want to overwrite this previous run, add the '--overwrite' flag to the command line."
        exit 1
    fi
else
    mkdir $out_directory
    mkdir -p ${out_directory}/Data ${out_directory}/Workspace ${out_directory}/metadata_files ${out_directory}/temp/ 
fi

# ==============================================================================
# MAIN SCRIPT
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Checking for duplicate headers."
check_duplicates "$fasta_path" "$out_directory"
echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> No duplicate headers found. Proceeding."

if [[ "$rip" -lt 100 ]]
then
    if [[ "$filter" == 0 ]]
    then
       echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: filtering elements with more than '${rip}' percent of sequence with RIP-like signal."
       Filter_input "$fasta_path" "$out_directory" "1"
       echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding."
    else
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: filtering elements with gc content below '${filter}' percent or with more than '${rip}' percent of sequence with RIP-like signal."
        Filter_input "$fasta_path" "$out_directory" "2"
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding..."
    fi
else
    if [[ "$filter" == 0 ]]
    then
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] skipping step 2 of gc content and rip filtering..."
        cp $fasta_path ${out_directory}/temp/Sequences.fa
    else
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: filtering elements with gc content below '${filter}' percent."
        Filter_input "$fasta_path" "$out_directory" "3"
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding..."
    fi
fi

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Processing headers, creating sequence header association file, and updating metadata if available."
Header_processing "${out_directory}/temp/Sequences.fa" "$out_directory"
echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding..."

if [[ "${mode}" == "Starfish" ]]; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Processing Starfish input files."
    Process_starfish "$out_directory" "$boundaries_path" "$gff_path"
elif [[ "${mode}" == "Simple" ]]; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Processing Simple input files."
    Process_simple "$out_directory" "$gff_path"
fi

rm -r ${out_directory}/temp/
echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> Process Complete."
