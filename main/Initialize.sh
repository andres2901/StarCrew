#!/bin/bash

# Function to print help message

function print_help() {
   echo "Script to perform a quick and raw prediction of genes for starships. It determines the statistics of each element prediction, filter the elements based on a minimum gene content (8) and organized it based on family association from the metadata."
   echo
   echo "Syntax: SAT Initialize [ -help ] -f <genome_file> [ -m <metadata_file> -o <name> ]"
   echo "options:"
   echo "-f, --fasta:  multifasta file wih the elements to study (required)."
   echo "-m, --metadata: csv file delimited by semicolon without headers with the information of the elements as follows: seqID;species;seqlength (optional)."
   echo "-n, --name: Specify working directory  name (Default = WorkingDirectory)."
   echo "-help: Display this help message."
}

# Initialize variables

out_directory="WorkingDirectory"
fasta_path=""
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
	    -n|--name)
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

# check if fasta file exist

if [[ ! -f "$fasta_path" ]]; then
    echo "Error: file '$fasta_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${fasta_path:0:1}" == "/" ]]; then
        fasta_path=$(realpath $fasta_path)
    fi
fi

# check if metadata file exist

if [[ -z "$metadata_path" ]]; then
    metadata=false
else
    if [[ ! -f "$metadata_path" ]]; then
        echo "Error: file '$metadata_path' does not exist."
        exit 1
    else
    metadata=true
    fi
fi

# check if output directory exist

if [[ -d "$out_directory" ]]; then
    echo "Error: directory '$out_directory' already exist."
    exit 1
else
    mkdir $out_directory
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

seqkit rmdup -n "$fasta_path" -o /dev/null -d temp_duplicated_headers &> /dev/null

if [ -s temp_duplicated_headers ]; then
    echo "Error: Duplicated headers found. Script terminated." >&2
    echo "Duplicate headers and their counts:" >&2
    echo "$duplicated_headers" >&2
    exit 1
fi

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> No duplicate headers found. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Processing headers, creating sequence header association file, and updating metadata if available."

mkdir -p ${out_directory}/Data ${out_directory}/Workspace ${out_directory}/metadata_files

# Define output file names
output_fasta="${out_directory}/sequences.fa"
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
' > temp_initials_and_headers.tsv

# Count the occurrences of each 5-char initial group
awk '{ print $1 }' temp_initials_and_headers.tsv | sort | uniq -c > temp_initial_counts.tsv

# Reset output files to ensure a clean start
> "$output_fasta"
echo "new_header;original_header" > "$association_csv"

# Create a temporary file to track used IDs to prevent collisions
temp_used_ids="temp_used_ids.txt"
> "$temp_used_ids"

# Loop through all headers and apply the renaming logic
while read -r initials original_header; do
    
    # Get the count for the current initial group
    group_count=$(grep -w "^[[:space:]]*[0-9]*[[:space:]]*$initials$" temp_initial_counts.tsv | awk '{print $1}')
    
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
    
    # Print to the temporary file for `seqkit`
    echo -e "$new_header\t$original_header" >> temp_association.tsv
done < temp_initials_and_headers.tsv

# Use seqkit to renames headers
awk 'BEGIN {OFS="\t"} {print $2, $1}' temp_association.tsv > temp_association2.tsv
seqkit replace --kv-file temp_association2.tsv -p "(.*)" -r "{kv}" "$fasta_path" > "$output_fasta"

# Create the final CSV association file
sed 's/\t/;/g' temp_association.tsv | sed '1s/^/new_header;original_header\n/' > "$association_csv"

# Create the updated metadata file

if $metadata; then
    metadata_csv="${out_directory}/metadata_files/metadata.csv"
    awk 'BEGIN{FS=OFS=";"}{if($2=="")$2="NA"; else if($3=="")$3="NA"; print}' $metadata_path > temp_metadata.csv
    #mv  $metadata_path
    join -1 2 -2 1 -t ';' <( sort -t ";" -k2,2 $association_csv) <(sort -t ";" -k1,1 temp_metadata.csv) > $metadata_csv
    #awk 'BEGIN{FS=OFS=";"}{print $1 OFS $2 OFS $3}'
    #| awk 'BEGIN {FS=OFS=";"} {for (i=3; i<=5; i++) if ($i == "" || $i == " ") $i = "NA"; print}'
else
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Skipping metadata file update..."
fi

# Cleanup temporary file
rm temp*

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> Processing complete. Output files created:"
echo "     - New FASTA file: $output_fasta"
echo "     - Association file: $association_csv"
if $metadata; then
echo "     - Update metadata file: $metadata_csv"
fi