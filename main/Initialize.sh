#!/usr/bin/env bash
# ==============================================================================
# NAME:        Initialize.sh
# DESCRIPTION: Organizes the working directory and prepares input data to run
#              the subsequent commands in the StarCrew workflow.
# USAGE:       StarCrew Initialize [options]
#              StarCrew Initialize -help
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        12/Mar/2026
# VERSION:     1.0.0
# ==============================================================================

set -uo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

readonly SCRIPT_DIR="$( dirname -- "$( readlink -f -- "$0" )" )"
readonly LIB_DIR="${SCRIPT_DIR}/../lib"
readonly AUXILIARY_DIR="${SCRIPT_DIR}/../aux"

readonly AGAT_CONFIG="${SCRIPT_DIR}/../agat_config.yaml"

# Source required libraries
if [[ ! -f "${LIB_DIR}/Utils.sh" ]]; then
  echo "Error: file '${LIB_DIR}/Utils.sh' does not exist." >&2
  exit 1
fi
source "${LIB_DIR}/Utils.sh"

if [[ ! -f "${LIB_DIR}/Check.sh" ]]; then
  echo "Error: file '${LIB_DIR}/Check.sh' does not exist." >&2
  exit 1
fi
source "${LIB_DIR}/Check.sh"

# Validate auxiliary directory and scripts
if [[ ! -d "$AUXILIARY_DIR" ]]; then
  echo "Error: directory '${AUXILIARY_DIR}' does not exist." >&2
  exit 1
fi
check_auxiliary_scripts "$(realpath "${AUXILIARY_DIR}")" "$(basename -s .sh "$0")"

# Validate required software
check_required_software "$(basename -s .sh "$0")"

# Validate agat config file
if [[ ! -f "$AGAT_CONFIG" ]]; then
  echo "Error: file '${AGAT_CONFIG}' does not exist." >&2
  exit 1
fi
readonly AGAT_CONFIG_PATH="$(realpath "${AGAT_CONFIG}")"

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays the help message with all available options.
# Arguments:
#   None
# Returns:
#   0 always
print_help() {
  echo "Script to organize the working directory to run the subsequent commands in the workflow."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -f <file_path> -g <file_path>"
  echo "       [ -m <string> -gc <integer> -r <integer> -mg <integer> -o <string>"
  echo "         { -b <file_path> -s <character> -c <file_path> } -M <file_path>"
  echo "         -t <integer> --overwrite ]"
  echo
  echo "Required args:"
  echo "  -f, --fasta         Multifasta file with the elements to study."
  echo "  -g, --gff           In 'Simple' mode: path to GFF with element-relative coordinates."
  echo "                      In 'Starfish' mode: 2-column TSV (genome code, path to GFF)."
  echo
  echo "Required args with defaults:"
  echo "  -m, --mode          Input mode (Default: Simple) [Available: Simple, Starfish]."
  echo "  -gc, --gc           Minimum GC content threshold; elements below are filtered out"
  echo "                      (Default: 0) [range: 20-45]."
  echo "  -r, --rip           Maximum RIP-like signal coverage; elements above are filtered out"
  echo "                      (Default: 0) [range: 30-80]."
  echo "  -mg, --minGene      Minimum number of genes per element (Default: 8) [range: 5-100]."
  echo "  -o, --outDirectory  Working directory name (Default: WorkingDirectory)."
  echo
  echo "Required args in 'Starfish' mode:"
  echo "  -b, --boundaries    *.elements.feat file from 'starfish summary'."
  echo "  -c, --captains      *_tyr.filt_intersect.fas file from 'starfish annotate'."
  echo "  -s, --separator     Character separating genomeID from featureID (Default: '_')."
  echo
  echo "Optional args:"
  echo "  -M, --Metadata      CSV file (semicolon-delimited) with metadata information."
  echo "  -t, --threads       Threads for MetaEuk in 'Starfish' mode (Default: 24)."
  echo "  --overwrite         Overwrite an existing previous run (Default: off)."
  echo "  -help               Display this help message."
}

# Checks for duplicate sequence headers in a FASTA file.
# Arguments:
#   $1 - path to the FASTA file
#   $2 - working directory path
# Returns:
#   0 if no duplicates found, exits with 1 if duplicates are detected
check_duplicates() {
  local fasta_path="$1"
  local working_dir="$2"

  seqkit rmdup -n "$fasta_path" -o /dev/null \
    -d "${working_dir}/temp/temp_duplicated_headers" &> /dev/null

  if [[ -s "${working_dir}/temp/temp_duplicated_headers" ]]; then
    echo "Error: duplicated headers found. Script terminated." >&2
    echo "Duplicate headers and their counts:" >&2
    cat "${working_dir}/temp/temp_duplicated_headers" >&2
    exit 1
  fi
}

# Filters input sequences based on RIP-like signal and/or GC content.
# Arguments:
#   $1 - path to the FASTA file
#   $2 - working directory path
#   $3 - filter option: 1=RIP only, 2=GC+RIP, 3=GC only
# Returns:
#   0 on success
filter_input() {
  local fasta_path="$1"
  local working_dir="$2"
  local filter_option="$3"
  local temp_dir="${working_dir}/temp"
  local removed_elements

  if [[ "$filter_option" == "1" || "$filter_option" == "2" ]]; then
    if [[ "$filter_option" == "2" ]]; then
      log_step "Filtering elements based on GC content..."
      seqkit fx2tab -g -n "$fasta_path" \
        | awk -v min="$gc_filter" '{if($NF < min){$NF=""; print $0}}' \
        | sed -e 's/ $//g' > "${working_dir}/Elements_filterGC.txt"
      removed_elements=$(wc -l < "${working_dir}/Elements_filterGC.txt")
      log_step "'${removed_elements}' elements removed based on GC content."
      seqkit grep --quiet -n -v -f "${working_dir}/Elements_filterGC.txt" \
        "$fasta_path" > "${temp_dir}/Sequences-filter1.fa"
      fasta_path="${temp_dir}/Sequences-filter1.fa"
    fi

    log_step "Filtering elements based on RIP-like signal..."
    python "${AUXILIARY_DIR}/rip_calculator.py" "$fasta_path" \
      -tc 0.01 -tp 1 -ts 1 -w 500 -s 100 \
      | awk -F '\t' -v min="$rip" 'NR>1{if($4 > min){print $1}}' \
      > "${working_dir}/Elements_filterRIPlike.txt"
    removed_elements=$(wc -l < "${working_dir}/Elements_filterRIPlike.txt")
    log_step "'${removed_elements}' elements removed based on RIP-like signal."
    seqkit grep --quiet -n -v -f "${working_dir}/Elements_filterRIPlike.txt" \
      "$fasta_path" > "${temp_dir}/Sequences.fa"

  elif [[ "$filter_option" == "3" ]]; then
    log_step "Filtering elements based on GC content..."
    seqkit fx2tab -g -n "$fasta_path" \
      | awk -v min="$gc_filter" '{if($NF < min){$NF=""; print $0}}' \
      | sed -e 's/ $//g' > "${working_dir}/Elements_filterGC.txt"
    removed_elements=$(wc -l < "${working_dir}/Elements_filterGC.txt")
    log_step "'${removed_elements}' elements removed based on GC content."
    seqkit grep --quiet -n -v -f "${working_dir}/Elements_filterGC.txt" \
      "$fasta_path" > "${temp_dir}/Sequences.fa"
  fi
}

# Converts a non-negative integer to a base-62 string for unique ID generation.
# Arguments:
#   $1 - non-negative integer to convert
# Returns:
#   Prints the base-62 string to stdout
base62() {
  local n="$1"
  local base62_chars="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  local result=""

  if [[ "$n" -eq 0 ]]; then
    echo "0"
    return
  fi

  while [[ "$n" -gt 0 ]]; do
    local rem=$(( n % 62 ))
    result="${base62_chars:$rem:1}${result}"
    n=$(( n / 62 ))
  done

  printf '%s\n' "$result"
}

# Generates a unique padded base-62 ID not already present in the used-IDs file.
# Writes the result and updated counter to a temp file to avoid subshell issues.
# Arguments:
#   $1 - current counter value (integer)
#   $2 - path to the file tracking used IDs
#   $3 - path to temp file to write "new_id counter" result line
# Returns:
#   0 always
generate_unique_id() {
  local counter="$1"
  local used_ids_file="$2"
  local out_file="$3"
  local padded_id

  while true; do
    padded_id=$(printf "%05s" "$(base62 "$counter")" | sed 's/ /0/g')
    if ! grep -q "^${padded_id}$" "$used_ids_file"; then
      echo "$padded_id" >> "$used_ids_file"
      counter=$(( counter + 1 ))
      echo "${padded_id} ${counter}" > "$out_file"
      return
    fi
    counter=$(( counter + 1 ))
  done
}

# Renames sequence headers to avoid issues during syntenet analysis,
# and creates a header-association CSV and an updated metadata file.
# Arguments:
#   $1 - path to the input FASTA file
#   $2 - working directory path
# Returns:
#   0 on success
process_headers() {
  local fasta_path="$1"
  local working_dir="$2"
  local output_fasta="${working_dir}/Sequences.fa"
  local association_csv="${working_dir}/metadata_files/sequence_head.csv"
  local temp_dir="${working_dir}/temp"
  local max_initial_count=25
  local id_counter=0
  local temp_used_ids="${temp_dir}/temp_used_ids.txt"
  local temp_id_result="${temp_dir}/temp_id_result.txt"

  # Build initial-to-header mapping
  seqkit fx2tab --name "$fasta_path" | awk -F'\t' '{
    original_header = $1
    cleaned_header  = original_header
    gsub(/[^a-zA-Z0-9]/, "", cleaned_header)
    initials = substr(cleaned_header, 1, 5)
    print initials, original_header
  }' > "${temp_dir}/temp_initials_and_headers.tsv"

  # Count occurrences per 5-char initial group
  awk '{print $1}' "${temp_dir}/temp_initials_and_headers.tsv" \
    | sort | uniq -c | sed 's/^ //g' > "${temp_dir}/temp_initial_counts.tsv"

  # Initialize output files
  > "$output_fasta"
  > "$temp_used_ids"
  > "${temp_dir}/temp_association.tsv"
  echo "new_header;original_header" > "$association_csv"

  local total_states
  total_states=$(wc -l < "${temp_dir}/temp_initials_and_headers.tsv")
  local state=0

  while IFS=$' ' read -r initials original_header; do
    state=$(( state + 1 ))

    local group_count
    group_count=$(awk -F' ' -v key="$initials" '$2==key{print $1}' \
      "${temp_dir}/temp_initial_counts.tsv")
    group_count="${group_count:-0}"

    local cleaned_header
    cleaned_header=$(echo "$original_header" | sed 's/[^a-zA-Z0-9]//g')

    local new_header=""

    if [[ "$group_count" -gt "$max_initial_count" ]]; then
      local current_last5
      current_last5=$(echo "$cleaned_header" | tail -c 6)

      if [[ "${#current_last5}" -eq 5 ]] && ! grep -q "^${current_last5}$" "$temp_used_ids"; then
        new_header="$current_last5"
        echo "$new_header" >> "$temp_used_ids"
      else
        generate_unique_id "$id_counter" "$temp_used_ids" "$temp_id_result"
        read -r new_header id_counter < "$temp_id_result"
      fi
    else
      if [[ "${#cleaned_header}" -gt 15 ]]; then
        generate_unique_id "$id_counter" "$temp_used_ids" "$temp_id_result"
        read -r new_header id_counter < "$temp_id_result"
      else
        new_header="$cleaned_header"
      fi

      if [[ "${#new_header}" -eq 5 ]]; then
        echo "$new_header" >> "$temp_used_ids"
      fi
    fi

    echo -e "${new_header}\t${original_header}" >> "${temp_dir}/temp_association.tsv"
    ProgressBar "$state" "$total_states"
  done < "${temp_dir}/temp_initials_and_headers.tsv"
  echo ""


  # Rename headers using seqkit
  log_step "Updating headers..."
  awk 'BEGIN{OFS="\t"}{print $2, $1}' "${temp_dir}/temp_association.tsv" \
    > "${temp_dir}/temp_association2.tsv"
  seqkit replace --kv-file "${temp_dir}/temp_association2.tsv" \
    -p "(.*)" -r "{kv}" "$fasta_path" \
    | seqkit seq --quiet -u > "$output_fasta"

  # Build final CSV association file
  sed 's/\t/;/g' "${temp_dir}/temp_association.tsv" \
    | sed '1s/^/new_header;original_header\n/' > "$association_csv"

  # Update metadata file if provided
  if $metadata; then
    local metadata_csv="${working_dir}/metadata_files/metadata.csv"
    head -n1 "$metadata_path" > "$metadata_csv"
    sed -i -e 's/ElementID/ElementID;ElementID_updated/' "$metadata_csv"
    awk 'BEGIN{FS=OFS=";"}{for(i=1;i<=NF;i++){if($i=="")$i="NA"}; print}' \
      "$metadata_path" > "${temp_dir}/temp_metadata.csv"
    join -1 2 -2 1 -t ';' \
      <(sort -t ";" -k2,2 "$association_csv") \
      <(sort -t ";" -k1,1 "${temp_dir}/temp_metadata.csv") \
      >> "$metadata_csv"
  else
    log_step "Skipping metadata file update..."
  fi
}

# Computes gene count and average gene/intergenic length for each element.
# Arguments:
#   $1 - working directory path
# Returns:
#   0 on success; writes results to Gene_stats.txt
compute_gene_stats() {
  local working_dir="$1"

  echo -e "Starship\tNumber_genes\tAvg_gene_length\tAvg_intergenic_length" \
    > "${working_dir}/Gene_stats.txt"

  ls "${working_dir}/Data/Gff/" | xargs -n1 basename -s .gff | while read -r element; do
    local gff="${working_dir}/Data/Gff/${element}.gff"
    echo -e "${element}\t$(grep -c -P "\tgene\t" "$gff")\t$(
      grep -P "\tgene\t" "$gff" \
        | awk '{sum += $5 - $4} END {if(NR>0) print sum/NR; else print $0}'
    )\t$(
      grep -P "\tgene\t" "$gff" | sort -k4 -n \
        | awk 'NR==1{prev=$5; next}{diff=$4-prev; prev=$5; if(diff>0) total+=diff}
               END{if((NR-1)>0) print total/(NR-1); else print 0}'
    )" >> "${working_dir}/Gene_stats.txt"
  done
}

# Filters elements by minimum gene content and generates per-element sequence,
# CDS, and protein files.
# Arguments:
#   $1 - working directory path
#   $2 - CDS flag (true/false): use CDS features if true, exon+merge if false
# Returns:
#   0 on success
organize_info() {
  local working_dir="$1"
  local cds_flag="$2"

  log_step "Computing gene statistics and removing low gene-model elements..."
  compute_gene_stats "$working_dir"

  awk -v min="$minimum_gene_content" 'NR>1{if($2>=min){print $1}}' \
    "${working_dir}/Gene_stats.txt" \
    | sed $'s/[^[:print:]\t]//g' > "${working_dir}/temp/Good_elements.txt"

  ls "${working_dir}/Data/Gff/" | xargs -n1 basename -s .gff \
    > "${working_dir}/temp/All_elements.txt"

  grep -v -f "${working_dir}/temp/Good_elements.txt" \
    "${working_dir}/temp/All_elements.txt" \
    | while read -r element; do
      rm "${working_dir}/Data/Gff/${element}.gff"
    done

  log_step "Organizing sequence and annotation data..."

  local total_states
  total_states=$(ls "${working_dir}/Data/Gff/" | wc -l)
  log_step "Processing '${total_states}' elements with at least '${minimum_gene_content}' genes."
  local state=0

  ls "${working_dir}/Data/Gff/" | xargs -n1 basename -s .gff | while read -r element; do
    state=$(( state + 1 ))

    # Extract nucleotide sequence for this element
    echo "$element" > "${working_dir}/temp/temp_element.txt"
    seqkit grep -n -f "${working_dir}/temp/temp_element.txt" \
      "${working_dir}/Sequences.fa" \
      -o "${working_dir}/Data/Nucleotide/${element}.fa" &> /dev/null

    # Extract CDS or exon sequences
    if $cds_flag; then
      agat_sp_extract_sequences.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${working_dir}/Data/Gff/${element}.gff" \
        --fasta "${working_dir}/Data/Nucleotide/${element}.fa" \
        -t cds -o "${working_dir}/temp/CDS/${element}.fa" &> /dev/null
    else
      agat_sp_extract_sequences.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${working_dir}/Data/Gff/${element}.gff" \
        --fasta "${working_dir}/Data/Nucleotide/${element}.fa" \
        -t exon --merge -o "${working_dir}/temp/CDS/${element}.fa" &> /dev/null
    fi

    awk '{if($2){$1=">"$2} print $1}' "${working_dir}/temp/CDS/${element}.fa" \
      | sed 's/gene=//g' > "${working_dir}/Data/CDS/${element}.fa"

    seqkit translate -f 1 "${working_dir}/Data/CDS/${element}.fa" \
      > "${working_dir}/Data/Protein/${element}.fa"

    rm -f "${working_dir}/Data/Nucleotide/${element}.fa.index"*
    ProgressBar "$state" "$total_states"
  done
  echo ""
}

# Processes input in 'Simple' mode: renames headers in GFF, splits by element,
# and calls organize_info.
# Arguments:
#   $1 - working directory path
#   $2 - path to the GFF file
# Returns:
#   0 on success
process_simple() {
  local working_dir="$1"
  local gff_file="$2"
  local cds_flag=false

  mkdir -p "${working_dir}/Data/Protein" "${working_dir}/Data/Gff" \
           "${working_dir}/Data/Nucleotide" "${working_dir}/Data/CDS" \
           "${working_dir}/temp/CDS" "${working_dir}/temp/gff"

  log_step "Preparing GFF file..."

  if [[ $(grep -w -i -c "cds" "$gff_file") -ge 10 ]]; then
    cds_flag=true
  fi

  local total_states
  total_states=$(wc -l < "${working_dir}/temp/temp_association.tsv")
  local state=0

  while read -r new_id original_id; do
    state=$(( state + 1 ))
    grep -w "^${original_id}" "$gff_file" \
      | sed "s/${original_id}/${new_id}/" \
      >> "${working_dir}/temp/gff_file.gff"
    ProgressBar "$state" "$total_states"
  done < "${working_dir}/temp/temp_association.tsv"

  gff_file="${working_dir}/temp/gff_file.gff"

  echo ""
  log_step "Splitting GFF file per element..."
  local unique_elements
  unique_elements=$(grep -v "^#" "$gff_file" | awk '{print $1}' | sort | uniq)
  total_states=$(echo "$unique_elements" | wc -l)
  state=0

  echo "$unique_elements" | while read -r element; do
    state=$(( state + 1 ))
    {
      grep "^#" "$gff_file"
      grep -w "^${element}" "$gff_file"
    } > "${working_dir}/temp/gff/${element}.gff"

    agat_sp_keep_longest_isoform.pl --config "${AGAT_CONFIG_PATH}" \
      --gff "${working_dir}/temp/gff/${element}.gff" \
      -o "${working_dir}/Data/Gff/${element}.gff" &> /dev/null

    sed -i -e "s/ID=/ID=${element}./g" \
           -e "s/Parent=/Parent=${element}./g" \
      "${working_dir}/Data/Gff/${element}.gff" &> /dev/null

    ProgressBar "$state" "$total_states"
  done
  echo ""

  organize_info "$working_dir" "$cds_flag"
}

# Processes input in 'Starfish' mode: builds coordinate file, slices GFF files,
# runs MetaEuk for captain identification, merges GFF models, and calls organize_info.
# Arguments:
#   $1 - working directory path
#   $2 - path to the boundaries file
#   $3 - path to the GFF paths TSV file
# Returns:
#   0 on success
process_starfish() {
  local working_dir="$1"
  local boundaries_path="$2"
  local gff_path="$3"
  local cds_flag=false

  mkdir -p "${working_dir}/Data/Protein" "${working_dir}/Data/Gff" \
           "${working_dir}/Data/Nucleotide" "${working_dir}/Data/CDS" \
           "${working_dir}/temp/CDS" "${working_dir}/temp/gff" \
           "${working_dir}/temp/gff2"

  local total_states
  total_states=$(wc -l < "${working_dir}/temp/temp_association.tsv")
  log_step "Processing '${total_states}' elements."

  # Build coordinate file
  log_step "Creating coordinate file..."
  awk -v s="$separator" 'BEGIN{FS=OFS="\t"}NR>1{
    split($1,array,s); print $2 OFS array[2] OFS $4 OFS $5 OFS $7 OFS array[1]
  }' "$boundaries_path" > "${working_dir}/temp/coordinate_file.txt"

  local state=0
  while read -r new_id original_id; do
    state=$(( state + 1 ))
    grep "^${original_id}" "${working_dir}/temp/coordinate_file.txt" \
      | sed "s/${original_id}/${new_id}/" \
      >> "${working_dir}/Coordinate_file.txt"
    ProgressBar "$state" "$total_states"
  done < "${working_dir}/temp/temp_association.tsv"

  # Prepare per-genome GFF files
  echo ""
  log_step "Preparing GFF files..."
  total_states=$(wc -l < "$gff_path")
  state=0

  while read -r genome_code gff_file; do
    state=$(( state + 1 ))
    check_gff_file "$gff_file"
    sed -e "s/${separator}//g" "$gff_file" \
      > "${working_dir}/temp/gff/${genome_code}.gff"
    ProgressBar "$state" "$total_states"
  done < "$gff_path"

  local first_gff
  first_gff=$(ls "${working_dir}/temp/gff/" | head -n1)
  if [[ $(grep -w -i -c "cds" "${working_dir}/temp/gff/${first_gff}") -ge 1 ]]; then
    cds_flag=true
  fi

  # Extract gene models per element
  echo ""
  log_step "Extracting gene model information per element..."
  total_states=$(awk '{print $NF}' "${working_dir}/Coordinate_file.txt" | sort -u | wc -l)
  state=0

  awk '{print $NF}' "${working_dir}/Coordinate_file.txt" | sort -u | while read -r element; do
    state=$(( state + 1 ))
    grep -w "$element" "${working_dir}/Coordinate_file.txt" \
      | cut -d$'\t' -f1-5 > "${working_dir}/temp/temp_coordinate_file.txt"

    python "${AUXILIARY_DIR}/gff_slicer.py" \
      -c "${working_dir}/temp/temp_coordinate_file.txt" \
      -i "${working_dir}/temp/gff/${element}.gff" \
      -o "${working_dir}/temp/gff2/" &> /dev/null

    awk '{print $1}' "${working_dir}/temp/temp_coordinate_file.txt" \
      | while read -r sub_element; do
        agat_sp_keep_longest_isoform.pl --config "${AGAT_CONFIG_PATH}" \
          --gff "${working_dir}/temp/gff2/${sub_element}.gff" \
          -o "${working_dir}/Data/Gff/${sub_element}.gff" &> /dev/null
      done
    ProgressBar "$state" "$total_states"
  done

  rm "${working_dir}/temp/gff2/"*.gff

  # Update GFF ID/Parent prefixes
  echo ""
  log_step "Updating element GFF files..."
  total_states=$(ls "${working_dir}/Data/Gff/" | wc -l)
  state=0

  ls "${working_dir}/Data/Gff/" | xargs -n1 basename -s .gff | while read -r element; do
    state=$(( state + 1 ))
    sed -i -e "s/ID=/ID=${element}./g" \
           -e "s/Parent=/Parent=${element}./g" \
      "${working_dir}/Data/Gff/${element}.gff"
    ProgressBar "$state" "$total_states"
  done

  # Captain identification via MetaEuk
  echo ""
  log_step "Identifying captains using MetaEuk based on Starfish output..."
  metaeuk createdb "${working_dir}/Sequences.fa" \
    "${working_dir}/temp/ContigsDB" --dbtype 2 -v 0 &> /dev/null
  metaeuk createdb "$captains_path" \
    "${working_dir}/temp/ProteinDB" --dbtype 1 -v 0 &> /dev/null

  metaeuk predictexons \
    "${working_dir}/temp/ContigsDB" "${working_dir}/temp/ProteinDB" \
    "${working_dir}/temp/metaeukResults" "${working_dir}/temp/tempFolder" \
    -s 7.5 --exhaustive-search 1 --orf-start-mode 0 \
    --min-seq-id 0.95 --metaeuk-tcov 0.75 --min-length 200 \
    --remove-tmp-files 1 --use-all-table-starts 1 \
    --threads "${threads}" --disk-space-limit 100G &> /dev/null

  metaeuk reduceredundancy \
    "${working_dir}/temp/metaeukResults" \
    "${working_dir}/temp/metaeukpred" \
    "${working_dir}/temp/metaeukgroups" \
    --threads "${threads}" -v 0 &> /dev/null

  metaeuk unitesetstofasta \
    "${working_dir}/temp/ContigsDB" "${working_dir}/temp/ProteinDB" \
    "${working_dir}/temp/metaeukpred" "${working_dir}/temp/metaeukFinal" \
    --threads "${threads}" -v 0 &> /dev/null

  sed -i -e 's/Target_ID=.*;TCS_//g' \
         -e 's/exon/CDS/g' \
    "${working_dir}/temp/metaeukFinal.gff"

  agat_convert_sp_gxf2gxf.pl --config "${AGAT_CONFIG_PATH}" \
    --gff "${working_dir}/temp/metaeukFinal.gff" \
    -o "${working_dir}/temp/metaeuk.gff" &> /dev/null

  agat_sp_filter_by_ORF_size.pl --config "${AGAT_CONFIG_PATH}" \
    --gff "${working_dir}/temp/metaeuk.gff" \
    -s 200 -o "${working_dir}/temp/metaeuk_ORF.gff" &> /dev/null

  # Merge GFF models per element
  log_step "Merging GFF files..."
  total_states=$(grep -v "#" "${working_dir}/temp/metaeuk_ORF_sup200.gff" \
    | awk '{print $1}' | sort -u | wc -l)
  state=0

  rm "${working_dir}/temp/gff/"*
  mkdir -p "${working_dir}/temp/modelsKeep"

  grep -v "#" "${working_dir}/temp/metaeuk_ORF_sup200.gff" \
    | awk '{print $1}' | sort -u | while read -r element; do
      state=$(( state + 1 ))

      echo "#gff version-3" \
        > "${working_dir}/temp/gff/temp_captain_${element}.gff"
      grep -w "^${element}" "${working_dir}/temp/metaeuk_ORF_sup200.gff" \
        >> "${working_dir}/temp/gff/temp_captain_${element}.gff"

      agat_sp_manage_IDs.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${working_dir}/temp/gff/temp_captain_${element}.gff" \
        --prefix "${element}.metaeuk" \
        -o "${working_dir}/temp/gff/temp_captain2_${element}.gff" &> /dev/null

      cp "${working_dir}/Data/Gff/${element}.gff" \
         "${working_dir}/temp/gff/${element}_merge.gff"
      grep -v "^#" "${working_dir}/temp/gff/temp_captain2_${element}.gff" \
        >> "${working_dir}/temp/gff/${element}_merge.gff"

      python "${AUXILIARY_DIR}/merge.py" \
        "${working_dir}/temp/gff/${element}_merge.gff" \
        "${working_dir}/temp/modelsKeep/${element}.txt" &> /dev/null

      agat_sp_filter_feature_from_keep_list.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${working_dir}/temp/gff/${element}_merge.gff" \
        --keep_list "${working_dir}/temp/modelsKeep/${element}.txt" \
        --output "${working_dir}/temp/gff/${element}_mergeKeep.gff" &> /dev/null

      rm "${working_dir}/Data/Gff/${element}.gff"
      agat_sp_keep_longest_isoform.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${working_dir}/temp/gff/${element}_mergeKeep.gff" \
        -o "${working_dir}/Data/Gff/${element}.gff" &> /dev/null

      ProgressBar "$state" "$total_states"
    done
  echo ""

  organize_info "$working_dir" "$cds_flag"
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

fasta_path=""
mode="Simple"
metadata_path=""
gc_filter="0"
rip="100"
minimum_gene_content="8"
out_directory="WorkingDirectory"
gff_path=""
boundaries_path=""
captains_path=""
separator="_"
threads="24"
overwrite=false
metadata=false

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--fasta)         shift; fasta_path="$1" ;;
    -m|--mode)          shift; mode="$1" ;;
    -M|--Metadata)      shift; metadata_path="$1" ;;
    -gc|--gc)           shift; gc_filter="$1" ;;
    -r|--rip)           shift; rip="$1" ;;
    -mg|--minGene)      shift; minimum_gene_content="$1" ;;
    -o|--outDirectory)  shift; out_directory="$1" ;;
    -g|--gff)           shift; gff_path="$1" ;;
    -b|--boundaries)    shift; boundaries_path="$1" ;;
    -c|--captains)      shift; captains_path="$1" ;;
    -s|--separator)     shift; separator="$1" ;;
    -t|--threads)       shift; threads="$1" ;;
    --overwrite)        overwrite=true ;;
    -help)
      print_help
      exit 0
      ;;
    *)
      echo "Error: invalid option '$1'." >&2
      print_help
      exit 1
      ;;
  esac
  shift
done

# ==============================================================================
# PARAMETER SUMMARY
# ==============================================================================

echo "Running $(basename -s .sh "$0") command with the following parameters:"
echo "  Fasta file:                      ${fasta_path}"
echo "  Mode:                            ${mode}"
echo "  Minimum GC (%):                  ${gc_filter}"
echo "  Maximum RIP-like signal (%):     ${rip}"
echo "  Minimum gene content:            ${minimum_gene_content}"
echo "  Output directory:                ${out_directory}"
echo "  GFF file/path:                   ${gff_path}"
echo "  Starfish boundaries file:        ${boundaries_path}"
echo "  Starfish captains file:          ${captains_path}"
echo "  Starfish separator:              ${separator}"
echo "  Metadata file:                   ${metadata_path}"
echo "  Threads:                         ${threads}"
echo "  Overwrite previous run:          ${overwrite}"
echo ""

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Checking arguments and input files..."

check_mode_parameter "$mode" "$(basename -s .sh "$0")"

if [[ "$mode" == "Simple" ]]; then
  if [[ -z "$fasta_path" || -z "$gff_path" ]]; then
    echo "Error: missing required argument(s) for 'Simple' mode." >&2
    print_help
    exit 1
  fi
  check_fasta_dna "$fasta_path"
  fasta_path=$(realpath "$fasta_path")
  check_gff_file "$gff_path"
  gff_path=$(realpath "$gff_path")

elif [[ "$mode" == "Starfish" ]]; then
  if [[ -z "$fasta_path" || -z "$gff_path" || -z "$boundaries_path" \
     || -z "$captains_path" || -z "$separator" ]]; then
    echo "Error: missing required argument(s) for 'Starfish' mode." >&2
    print_help
    exit 1
  fi
  check_fasta_dna "$fasta_path"
  fasta_path=$(realpath "$fasta_path")
  check_gff_paths "$gff_path"
  gff_path=$(realpath "$gff_path")
  check_boundaries_file "$boundaries_path"
  boundaries_path=$(realpath "$boundaries_path")
  check_fasta_protein "$captains_path"
  captains_path=$(realpath "$captains_path")

  impossible_separator=":;|"
  impossible_separator_pattern="[${impossible_separator}]"
  if [[ "${separator}" =~ "${impossible_separator_pattern}" ]]; then
    echo "Error: separator '${separator}' is not accepted." >&2
    exit 1
  fi
fi

if [[ -z "$metadata_path" ]]; then
  metadata=false
else
  check_metadata_file "$metadata_path" "$fasta_path"
  metadata=true
  metadata_path=$(realpath "$metadata_path")
fi

if [[ "$minimum_gene_content" =~ ^[0-9]+$ ]]; then
  if (( minimum_gene_content < 5 || minimum_gene_content > 100 )); then
    echo "Error: minimum gene content '${minimum_gene_content}' is out of range [5-100]." >&2
    print_help; exit 1
  fi
else
  echo "Error: minimum gene content '${minimum_gene_content}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$gc_filter" =~ ^[0-9]+$ ]]; then
  if (( gc_filter != 0 && (gc_filter < 20 || gc_filter > 45) )); then
    echo "Error: GC content value '${gc_filter}' is out of range [20-45]." >&2
    print_help; exit 1
  fi
else
  echo "Error: GC content value '${gc_filter}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$rip" =~ ^[0-9]+$ ]]; then
  if (( rip != 100 && (rip < 30 || rip > 80) )); then
    echo "Error: RIP coverage value '${rip}' is out of range [30-80]." >&2
    print_help; exit 1
  fi
else
  echo "Error: RIP coverage value '${rip}' is not a positive integer." >&2
  print_help; exit 1
fi

check_threads "${threads}"

if [[ -d "$out_directory" ]]; then
  if $overwrite; then
    rm -r "$out_directory"
  else
    echo "Error: directory '${out_directory}' already exists." >&2
    echo "Use '--overwrite' to overwrite a previous run." >&2
    exit 1
  fi
fi

mkdir -p "${out_directory}/Data" "${out_directory}/Workspace" \
         "${out_directory}/metadata_files" "${out_directory}/temp"

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

log_info "Step 1: Checking for duplicate headers."
check_duplicates "$fasta_path" "$out_directory"
log_step "-> No duplicate headers found. Proceeding."

# Step 2: Sequence filtering
if (( rip < 100 )) && (( gc_filter == 0 )); then
  log_info "Step 2: Filtering elements with more than '${rip}%' RIP-like signal."
  filter_input "$fasta_path" "$out_directory" "1"
elif (( rip < 100 )) && (( gc_filter > 0 )); then
  log_info "Step 2: Filtering elements with GC < '${gc_filter}%' or RIP-like signal > '${rip}%'."
  filter_input "$fasta_path" "$out_directory" "2"
elif (( rip == 100 )) && (( gc_filter == 0 )); then
  log_info "Step 2: Skipping GC and RIP filtering."
  cp "$fasta_path" "${out_directory}/temp/Sequences.fa"
else
  log_info "Step 2: Filtering elements with GC content below '${gc_filter}%'."
  filter_input "$fasta_path" "$out_directory" "3"
fi
log_step "Proceeding."

log_info "Step 3: Processing headers and creating sequence association file."
process_headers "${out_directory}/temp/Sequences.fa" "$out_directory"
log_step "Proceeding."

log_info "Step 4: Processing input files in '${mode}' mode."
if [[ "$mode" == "Starfish" ]]; then
  process_starfish "$out_directory" "$boundaries_path" "$gff_path"
elif [[ "$mode" == "Simple" ]]; then
  process_simple "$out_directory" "$gff_path"
fi

rm -r "${out_directory}/temp/"
log_step "-> Process complete."
