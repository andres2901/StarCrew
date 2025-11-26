#!/bin/bash

# Function to print help message
function print_help() {
   echo "Script to organize the working directory to run the subsequent commands in the workflow."
   echo
   echo "Syntax: SAT Initialize [ -help ] -f <filte_path> [ -m <file_path> -o <string> -g <integer> -r ]"
   echo "options:"
   echo "-f, --fasta:  multifasta file wih the elements to study (required)."
   echo "-m, --metadata: csv file delimited by semicolon without headers with the information of the elements as follows: seqID;species;seqlength (optional)."
   echo "-g, --gff: 2 column tsv: genome code, path to GFF. The path should be to the original gff files and not the ones formatted to run starfish (Optional)."
   echo "-b, --boundaries: *.elements.feat file output of 'starfish summary' command (Optional, Mandatory if used the -g/--gff parameter)."
   echo "-s, --separator: character separating genomeID from featureID that was used for starfish run (Optional, Mandatory if used the -g/--gff parameter)."
   echo "-c, --captains: *_tyr.filt_intersect.fas file output of 'starfish annotate' command (Optional, Mandatory if used the -g/--gff parameter)."
   echo "-gc, --gc: integer value of gc content to filter out elements with too low gc content (Default = 0) [range: 20 - 45]"
   echo "-o, --outDirectory: Specify working directory  name (Default = WorkingDirectory)."
   echo "-r, --rip: integer value of the minimum coverage of the element to be possibly affected by RIP to be filter out (Default = 0) [range: 30 - 80]"
   echo "-help: Display this help message."
}

# Initialize variables
out_directory="WorkingDirectory"
filter="0"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
rip="0"
fasta_path=""
gff_path=""
boundaries_path=""
separator=""
captains_path=""
metadata_path=""
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -f|--fasta) 
	        shift
	        fasta_path="$1"
	        ;;
         -m|--metadata)
            shift
            metadata_path="$1"
            ;;
        -g|--gff)
            shift
            gff_path="$1"
            ;;
        -b|--boundaries)
            shift
            boundaries_path="$1"
            ;;
        -s|--separator)
            shift
            separator="$1"
            ;;
        -c|--captains)
            shift
            captains_path="$1"
            ;;
        -gc|--gc)
            shift
            filter="$1"
            ;;
        -r|--rip)
            shift
            rip="$1"
            ;;
	    -o|--outDirectory)
            shift
	        out_directory="$1"
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

# Check for mandatory arguments
if [[ -z "$fasta_path" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if fasta file exist
if [[ ! -f "$fasta_path" ]]; then
    echo "Error: file '$fasta_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${fasta_path:0:1}" == "/" ]]; then
        fasta_path=$(realpath $fasta_path)
    fi
fi

# Check if metadata file exist
if [[ -z "$metadata_path" ]]; then
    metadata=false
else
    if [[ ! -f "$metadata_path" ]]; then
        echo "Error: file '$metadata_path' does not exist."
        exit 1
    else
        if [[ ! "${metadata_path:0:1}" == "/" ]]; then
            metadata_path=$(realpath $metadata_path)
        fi
        metadata=true
    fi
fi

# Check for starfish output
if [[ ! -z "$gff_path" ]]; then
    if [[ ! -f "$gff_path" ]]; then
        echo "Error: file '$gff_path' does not exist."
        exit 1
    else
        if [[ ! "${gff_path:0:1}" == "/" ]]; then
            gff_path=$(realpath $gff_path)
        fi
    fi

    if [[ -z "${boundaries_path}" || -z "${separator}" || -z "${captains_path}" ]]; then
        echo "Error: missing argument for starfish related initialization."
        print_help
        exit 1
    else
        if [[ ! -f "$boundaries_path" ]]; then
            echo "Error: file '$boundaries_path' does not exist."
            exit 1
        else
        # Check if the path is absolute
            if [[ ! "${boundaries_path:0:1}" == "/" ]]; then
                boundaries_path=$(realpath $boundaries_path)
            fi
        fi

#        if [[ ! -f "$captains_path" ]]; then
#            echo "Error: file '$captains_path' does not exist."
#            exit 1
#        else
        # Check if the path is absolute
#            if [[ ! "${captains_path:0:1}" == "/" ]]; then
#                captains_path=$(realpath $captains_path)
#            fi
#        fi

        impossible_separator=":;|"
        impossible_separator_pattern="[${impossible_separator}]"
        if [[ "s${separator}" == "${impossible_separator_pattern}" ]]; then
            echo "Error: '$separator' is not accepted."
            exit 1
        fi
    fi
fi

# Check filter parameter
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
    if (( $rip != 0 && ($rip < 30 || $rip > 80) )); then
        echo "Error: '$rip' rip coverage is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$rip' is not a positive integer."
    print_help
    exit 1
fi

# Check for rip calculator
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
    # Check the presence of the specific scripts
    if [[ ! -f "${auxiliary_path}/rip_calculator.py" ]]; then
        echo "Error: file '${auxiliary_path}/rip_calculator.py' does not exist."
        exit 1
    fi

    if [[ ! -f "${auxiliary_path}/gff_slicer.py" ]]; then
        echo "Error: file '${auxiliary_path}/gff_slicer.py' does not exist."
        exit 1
    fi
fi

# Check for required software
if [[ -z "$(which seqkit)" ]]; then
    echo "Error: Missing seqkit function."
    exit 1
fi

# check if output directory exist
if [[ -d "$out_directory" ]]; then
    echo "Error: directory '$out_directory' already exist."
    exit 1
else
    mkdir $out_directory
    mkdir ${out_directory}/temp/
fi

# ==============================================================================
# Check and prepare fasta file
# ==============================================================================

# Character set for generating unique 5-character IDs
# Using base62 for a wide range of short IDs
base62_chars="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
id_counter=0
max_initial_count=25 # Threshold for assigning a new ID

# Function to convert a number to a base-62 string for unique IDs
base62() {
    local n=$1
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

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Checking for duplicate headers."

# Verify if ther is any header duplicated.
seqkit rmdup -n "$fasta_path" -o /dev/null -d ${out_directory}/temp/temp_duplicated_headers &> /dev/null

if [ -s ${out_directory}/temp/temp_duplicated_headers ]; then
    echo "Error: Duplicated headers found. Script terminated." >&2
    echo "Duplicate headers and their counts:" >&2
    echo "$duplicated_headers" >&2
    exit 1
fi

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> No duplicate headers found. Proceeding."

# Filter stage
if [[ $rip > 0 ]]
then
    if [[ $filter == 0 ]]
    then
       echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: filtering elements with more than '${rip}' percent of sequence with RIP-like signal."
       echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on RIP-like signal"
       python ${auxiliary_path}/rip_calculator.py "$fasta_path" -tc 0.01 -tp 1 -ts 1 -w 500 -s 100 | awk -F '\t' -v min="$rip" 'NR>1{if($4 > min){print $1}}' > ${out_directory}/Elements_filterRIPlike.txt
       removed_elements=$(wc -l ${out_directory}/Elements_filterRIPlike.txt | awk '{print $1}')
       echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove based on RIP-like signal"
       seqkit grep --quiet -n -v -f ${out_directory}/Elements_filterRIPlike.txt "$fasta_path" > ${fasta_path}.filtered.fa
       fasta_path=$(realpath ${fasta_path}.filtered.fa)
       echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding."
    else
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: filtering elements with gc content below '${filter}' percent or with more than '${rip}' percent of sequence with RIP-like signal."
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on gc content"
        seqkit fx2tab -g -n "$fasta_path" | awk -v min="$filter" '{if($NF < min){$NF=""; print $0}}' | sed -e 's/ $//g' > ${out_directory}/Elements_filterGC.txt
        removed_elements=$(wc -l ${out_directory}/Elements_filterGC.txt | awk '{print $1}')
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove based on gc content"
        seqkit grep --quiet -n -v -f ${out_directory}/Elements_filterGC.txt "$fasta_path" > ${fasta_path}.filtered.fa

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on RIP-like signal"
        python ${auxiliary_path}/rip_calculator.py "${fasta_path}.filtered.fa" -tc 0.01 -tp 1 -ts 1 -w 500 -s 100 | awk -F '\t' -v min="$rip" 'NR>1{if($4 > min){print $1}}' > ${out_directory}/Elements_filterRIPlike.txt
        removed_elements=$(wc -l ${out_directory}/Elements_filterRIPlike.txt | awk '{print $1}')
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove based on RIP-like signal"
        seqkit grep --quiet -n -v -f ${out_directory}/Elements_filterRIPlike.txt "${fasta_path}.filtered.fa" > ${fasta_path}.filtered2.fa

        fasta_path=$(realpath ${fasta_path}.filtered2.fa)
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding..."
    fi
else
    if [[ $filter == 0 ]]
    then
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] skipping step 2 of gc content and rip filtering..."
    else
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: filtering elements with gc content below '${filter}' percent."
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] filtering elements based on gc content"
        seqkit fx2tab -g -n "$fasta_path" | awk -v min="$filter" '{if($NF < min){$NF=""; print $0}}' | sed -e 's/ $//g' > ${out_directory}/Elements_filterGC.txt
        removed_elements=$(wc -l ${out_directory}/Elements_filterGC.txt | awk '{print $1}')
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] '${removed_elements}' elements were remove base on gc content."
        seqkit grep --quiet -n -v -f ${out_directory}/Elements_filterGC.txt "$fasta_path" > ${fasta_path}.filtered.fa
        fasta_path=$(realpath ${fasta_path}.filtered.fa)
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Proceeding..."
    fi
fi

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Processing headers, creating sequence header association file, and updating metadata if available."

mkdir -p ${out_directory}/Data ${out_directory}/Workspace ${out_directory}/metadata_files

# Define output file names
output_fasta="${out_directory}/Sequences.fa"
association_csv="${out_directory}/metadata_files/sequence_head.csv"

# Create a temporary file with initials and original headers
seqkit fx2tab --name "$fasta_path" | awk -F'\t' '
    {
        original_header = $1
        cleaned_header = original_header
        gsub(/[^a-zA-Z0-9_]/, "", cleaned_header)
        initials = substr(cleaned_header, 1, 5)
        print initials, original_header
    }
' > ${out_directory}/temp/temp_initials_and_headers.tsv

# Count the occurrences of each 5-char initial group
awk '{ print $1 }' ${out_directory}/temp/temp_initials_and_headers.tsv | sort | uniq -c > ${out_directory}/temp/temp_initial_counts.tsv

# Reset output files to ensure a clean start
> "$output_fasta"
echo "new_header;original_header" > "$association_csv"

# Create a temporary file to track used IDs to prevent collisions
temp_used_ids="temp_used_ids.txt"
> "$temp_used_ids"

# Loop through all headers and apply the renaming logic
while read -r initials original_header; do
    
    # Get the count for the current initial group
    group_count=$(grep -w "^[[:space:]]*[0-9]*[[:space:]]*$initials$" ${out_directory}/temp/temp_initial_counts.tsv | awk '{print $1}')
    
    # Clean the header for processing
    cleaned_header=$(echo "$original_header" | sed 's/[^a-zA-Z0-9_]//g')
    
    new_header=""
    
    # Check if this header belongs to an over-threshold group
    if [ "$group_count" -gt "$max_initial_count" ]; then
        # Check if the header's own last 5 characters can be used
        current_last5=$(echo "$cleaned_header" | tail -c 6)
        
        # Check if the last 5 chars are a valid length and not already used
        if [ "${#current_last5}" -eq 5 ] && ! grep -q "^$current_last5$" "$temp_used_ids"; then
            new_header="$current_last5"
        else
            # Fallback to generating a unique ID
            UNIQUE_ID=$(base62 "$ID_COUNTER")
            PADDED_ID=$(printf "%05s" "$UNIQUE_ID" | sed 's/ /0/g')
            new_header="$PADDED_ID"
            ID_COUNTER=$((ID_COUNTER + 1))
        fi
    else
        # The group is within the limit, use the cleaned header
        new_header="$cleaned_header"
    fi
    
    # Add the new header to the list of used IDs to prevent future collisions
    echo "$new_header" >> "$temp_used_ids"
    
    # Print to the temporary file use for renamed
    echo -e "$new_header\t$original_header" >> ${out_directory}/temp/temp_association.tsv
done < ${out_directory}/temp/temp_initials_and_headers.tsv

# Use seqkit to renames headers
awk 'BEGIN {OFS="\t"} {print $2, $1}' ${out_directory}/temp/temp_association.tsv > ${out_directory}/temp/temp_association2.tsv
seqkit replace --kv-file ${out_directory}/temp/temp_association2.tsv -p "(.*)" -r "{kv}" "$fasta_path" | seqkit seq -u > "$output_fasta"

# Create the final CSV association file
sed 's/\t/;/g' ${out_directory}/temp/temp_association.tsv | sed '1s/^/new_header;original_header\n/' > "$association_csv"

# Create the updated metadata file
if $metadata; then
    metadata_csv="${out_directory}/metadata_files/metadata.csv"
    awk 'BEGIN{FS=OFS=";"}{if($2=="")$2="NA"; else if($3=="")$3="NA"; print}' $metadata_path > ${out_directory}/temp/temp_metadata.csv
    join -1 2 -2 1 -t ';' <( sort -t ";" -k2,2 $association_csv) <(sort -t ";" -k1,1 ${out_directory}/temp/temp_metadata.csv) > $metadata_csv
else
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Skipping metadata file update..."
fi

if [[ ! -z "${gff_path}" ]]; then
    mkdir -p ${out_directory}/Data/Protein/ ${out_directory}/Data/Gff/ ${out_directory}/Data/Nucleotide/ ${out_directory}/Data/Exon/ ${out_directory}/temp/exon/
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating coordinate file..."
    awk -v s="$separator" 'BEGIN{FS=OFS="\t"}NR>1{split($1,array,s);print $2 OFS array[2] OFS $4 OFS $5 OFS $7 OFS array[1]}' ${boundaries_path} > ${out_directory}/temp/coordinate_file.txt

    cat ${out_directory}/temp/temp_association.tsv | while read line; do     A=$(echo $line | awk '{print $2}' | awk -F '|' '{print $1}'); B=$(echo $line | awk '{print $1}'); grep $A ${out_directory}/temp/coordinate_file.txt | sed s/$A/$B/ >> ${out_directory}/Coordinate_file.txt; done
    
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Preparing gff files..."
    mkdir ${out_directory}/temp/gff/
    cat ${gff_path} | while read line
    do
        GFF_FILE=$(grep -w "${line}" ${gff_path} | awk '{print $2}')
        if [[ -f ${GFF_FILE} ]]; then
            sed -e s/${separator}//g $(echo $line | awk '{print $2}') > ${out_directory}/temp/gff/$(grep -w "${line}" ${gff_path} | awk '{print $1}').gff
        else
            echo "ERROR: Do not find GFF file for $(grep "${line}" ${gff_path} | awk '{print $1}') genome"
            exit 1
        fi
    done
        
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing gff file per element..."
    awk '{print $NF}' ${out_directory}/Coordinate_file.txt | sort -u | while read line 
    do
        grep -w $line ${out_directory}/Coordinate_file.txt | cut -d$'\t' -f 1-5 > ${out_directory}/temp/temp_coordinate_file.txt
        python ${auxiliary_path}/gff_slicer.py -c ${out_directory}/temp/temp_coordinate_file.txt -i ${out_directory}/temp/gff/${line}.gff -o ${out_directory}/Data/Gff/
    done

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Organizing info..."

    grep ">" $output_fasta | awk -F '>' '{print $2}' | while read line 
    do
        # Dividing nucleotide sequence of elements
        echo $line > ${out_directory}/temp/temp_element.txt
        seqkit grep -n -f ${out_directory}/temp/temp_element.txt $output_fasta -o ${out_directory}/Data/Nucleotide/${line}.fa &> /dev/null

        # Creating Exome
        agat_sp_extract_sequences.pl --gff ${out_directory}/Data/Gff/${line}.gff --fasta $output_fasta -t cds --merge -o ${out_directory}/temp/exon/${line}.fa &> /dev/null
        awk '{if($2){$1=">"$2} print $1}' ${out_directory}/temp/exon/${line}.fa| sed 's/gene=//g' > ${out_directory}/Data/Exon/${line}.fa

        # Creating proteome
        seqkit translate ${out_directory}/Data/Exon/${line}.fa --trim > ${out_directory}/Data/Protein/${line}.fa
    done

fi

# Cleanup temporary file
#rm -r ${out_directory}/temp/

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> Processing complete. Output files created:"
echo "     - New FASTA file: $output_fasta"
echo "     - Association file: $association_csv"
if $metadata; then
echo "     - Update metadata file: $metadata_csv"
fi
