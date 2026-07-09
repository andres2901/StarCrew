#!/usr/bin/env bash
# ==============================================================================
# NAME:        GeneHomogenization.sh
# DESCRIPTION: Perform an homogenization of the gene predictions in a given 
#              dataset
#              Executes 3 main steps:
#              1. Look-up for elements that share at least 50% of coverage using BLASTN.
#              2. Lift and merge gene models.
#              3. Updated the dataset.
# USAGE:       StarCrew GeneHomogenization [options]
#              StarCrew GeneHomogenization -help
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        9/Jul/2026
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

if [[ ! -d "$AUXILIARY_DIR" ]]; then
  echo "Error: directory '${AUXILIARY_DIR}' does not exist." >&2
  exit 1
fi
check_auxiliary_scripts "$(realpath "${AUXILIARY_DIR}")" "$(basename -s .sh "$0")"

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
  echo "Command to Perform an homogenization of the gene predictions in a given dataset."
  echo "Performs 3 main steps:"
  echo "  1. Look-up for elements that share at least 50% of coverage using BLASTN."
  echo "  2. Lift and merge gene models."
  echo "  3. Updated the dataset."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -w <directory_path>"
  echo "       [ -t <integer> --overwrite ]"
  echo
  echo "Required args:"
  echo "  -w, --workingDirectory  Working directory where all data are stored."
  echo
  echo "Optional args:"
  echo "  -t, --threads  Number of threads (Default: 8)."
  echo "  --overwrite  Overwrite a previous run (Default: off)."
  echo "  -help        Display this help message."
}

# Sets up the workspace directory structure for this command.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if a previous workspace run exists
organize_working_directory() {
  local base_dir="$1"
  local data_dir="${base_dir}/Data"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  # Locate required input subdirectories
  local gff_dir nucleotide_dir protein_dir cds_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  if [[ -d "$working_dir" ]]; then
    echo "Error: a previous run exists in '${working_dir}'." >&2
    echo "Use '--overwrite' to overwrite it." >&2
    exit 1
  fi

  mkdir -p "$working_dir" "$temp_dir"

  # Copy required data into workspace
  cp -r "$gff_dir" "$working_dir"
  cp -r "$nucleotide_dir" "$working_dir"
  cp -r "$protein_dir" "$working_dir"
  cp -r "$cds_dir" "$working_dir"
}

# Runs BLASTN of each element against the full set and filter the elements to which share 50% coverage
# Arguments:
#   $1 - - base working directory path
# Returns:
#   0 on success
nucleotide_match() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local match_dir="${working_dir}/Match/"

  mkdir "$match_dir"

  ls ${working_dir}/Nucleotide/ | sed -e 's/\.fa//' > "${temp_dir}/queries.txt"

  local query_list="${temp_dir}/queries.txt"

  cat ${working_dir}/Nucleotide/*.fa > "${working_dir}/Sequences.fa"

  makeblastdb -dbtype nucl \
      -in "${working_dir}/Sequences.fa" &> /dev/null

  local total_states
  total_states=$(wc -l < "$query_list")
  local state=0
 
  log_step "Running BLASTN for ${total_states} elements..."
 
  while read -r query; do
    state=$(( state + 1 ))
       
    blastn \
      -query "${working_dir}/Nucleotide/${query}.fa" \
      -db "${working_dir}/Sequences.fa" \
      -task blastn \
      -gapopen 8 -gapextend 6 -reward 5 -penalty -4 \
      -evalue 1e-60 -num_threads "$threads" \
      -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" \
      > "$temp_dir/${query}-blastraw.txt" 2> /dev/null

    python "${AUXILIARY_DIR}/Blast_CleanUp.py" \
      -f "$temp_dir/${query}-blastraw.txt" \
      -o "$temp_dir/${query}-blastclean.txt" \
      -fs "2000" -i "70" -m "Classification" \
      -ms "5000" -c "50" &> ${temp_dir}/blast_cleanup_log.txt

    cut -f2 "$temp_dir/${query}-blastclean.txt" | sort -u > "${match_dir}/${query}.out"
 
    ProgressBar "$state" "$total_states"
  done < "$query_list"
  echo ""

  find ${match_dir}/ -size 0 -delete

  rm -r ${temp_dir}/* 2> /dev/null
}

# Runs liftoff per element using as reference a pool of similar elements, filter and merge the results.
# Arguments:
#   $1 - - base working directory path
# Returns:
#   0 on success
run_liftoff() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local match_dir="${working_dir}/Match/"
  local liftoff_dir="${working_dir}/liftoff/"

  gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)

  mkdir "${liftoff_dir}"

  ls ${match_dir}/ | sed -e 's/\.out//' > "${temp_dir}/queries.txt"

  local query_list="${temp_dir}/queries.txt"

  local total_states
  total_states=$(wc -l < "$query_list")
  local state=0
 
  log_step "Running liftoff for ${total_states} elements..."
 
  while read -r query; do
    state=$(( state + 1 ))

    mkdir "${temp_dir}/${query}-liftoff"

    echo "#gff3" > "${temp_dir}/${query}-liftoff/Subject.gff"

    cat ${match_dir}/${query}.out | while read subject;
    do
      cat ${nucleotide_dir}/${subject}.fa >> "${temp_dir}/${query}-liftoff/Subject.fa"
      grep -v "^#" ${gff_dir}/${subject}.gff >> "${temp_dir}/${query}-liftoff/Subject.gff"
    done

    liftoff -g "${temp_dir}/${query}-liftoff/Subject.gff" -dir "${temp_dir}/${query}-liftoff/intermediate_files" -o ${temp_dir}/${query}-liftoff/${query}-liftoff.gff \
    ${nucleotide_dir}/${query}.fa "${temp_dir}/${query}-liftoff/Subject.fa" -polish -exclude_partial &>/dev/null

    grep "valid_ORFs=0" ${temp_dir}/${query}-liftoff/${query}-liftoff.gff_polished | cut -f9 | cut -d";" -f1 | sed -e 's/ID=//' > ${temp_dir}/${query}-liftoff/liftoffremove.txt
    if [ -s ${temp_dir}/${query}-liftoff/liftoffremove.txt ]; then
      agat_sp_filter_feature_from_kill_list.pl --config "${AGAT_CONFIG_PATH}" \
        --gff ${temp_dir}/${query}-liftoff/${query}-liftoff.gff_polished --kill_list ${temp_dir}/${query}-liftoff/liftoffremove.txt \
        --output ${temp_dir}/${query}-liftoff/liftoff-final.gff &>/dev/null
      agat_sp_merge_annotations.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${gff_dir}/${query}.gff" --gff "${temp_dir}/${query}-liftoff/liftoff-final.gff" \
        --out "${temp_dir}/${query}-liftoff/Merge.gff" &>/dev/null
    else
      agat_sp_merge_annotations.pl --config "${AGAT_CONFIG_PATH}" \
        --gff "${gff_dir}/${query}.gff" --gff ${temp_dir}/${query}-liftoff/${query}-liftoff.gff_polished \
        --out "${temp_dir}/${query}-liftoff/Merge.gff" &>/dev/null
    fi

    python ${AUXILIARY_DIR}/merge.py "${temp_dir}/${query}-liftoff/Merge.gff" "${temp_dir}/${query}-liftoff/ModelsKeep.txt" &>/dev/null
    agat_sp_filter_feature_from_keep_list.pl --config "${AGAT_CONFIG_PATH}" \
      --gff "${temp_dir}/${query}-liftoff/Merge.gff" --keep_list "${temp_dir}/${query}-liftoff/ModelsKeep.txt" \
      --output "${temp_dir}/${query}-liftoff/Merge-Keep.gff" &>/dev/null
    agat_sp_keep_longest_isoform.pl --config "${AGAT_CONFIG_PATH}" \
      --gff "${temp_dir}/${query}-liftoff/Merge-Keep.gff" -o "${liftoff_dir}/${query}.gff" &>/dev/null
    sed -i -e "s/ID=/ID=${query}./g" \
           -e "s/Parent=/Parent=${query}./g" "${liftoff_dir}/${query}.gff"

    ProgressBar "$state" "$total_states"
  done < "$query_list"
  echo ""

  grep -v "^#" ${liftoff_dir}/* | grep -w "Liftoff" | cut -f1 | cut -d":" -f2 | sort -u > ${working_dir}/Updated_elements.txt
  local updated_elements=$(wc -l < "${working_dir}/Updated_elements.txt")

  rm -r ${temp_dir}/*

  log_step "${updated_elements} elements have updated gene models"
}

# Moves elements information with updated gene models to a save space
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
store_old_gene_data() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local data_dir
  data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

  local gff_dir protein_dir cds_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  local store_elements="${working_dir}/Updated_elements.txt"
  local output="${data_dir}/Old_gene_models"
  mkdir -p "$output" "${output}/Gff" "${output}/Protein" "${output}/CDS"

  log_step "storing old gene data from elements with updated gene models..."

  while read -r element; do
    mv "${gff_dir}/${element}.gff"        "${output}/Gff/"
    mv "${protein_dir}/${element}.fa"     "${output}/Protein/"
    mv "${cds_dir}/${element}.fa"         "${output}/CDS/"
  done < "$store_elements"
}

# Organize and update GFF, CDS, and protein files per element.
# Arguments:
#   $1 - working directory path
# Returns:
#   0 on success
organize_information() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local updated_elements="${working_dir}/Updated_elements.txt"

  local total_states
  total_states=$(wc -l < "$updated_elements")
  log_step "Organizing annotation data for '${total_states}' elements."
  local state=0

  cat ${updated_elements} | while read -r element; 
  do
    state=$(( state + 1 ))
    cp "${working_dir}/liftoff/${element}.gff" "${base_dir}/Data/Gff/"
    agat_sp_extract_sequences.pl --config "${AGAT_CONFIG_PATH}" \
      --gff "${base_dir}/Data/Gff/${element}.gff" \
      --fasta "${base_dir}/Data/Nucleotide/${element}.fa" \
      -t cds -o "${temp_dir}/${element}.fa" &> /dev/null

    awk '{if($2){$1=">"$2} print $1}' "${temp_dir}/${element}.fa" \
    | sed 's/gene=//g' > "${base_dir}/Data/CDS/${element}.fa"

    seqkit translate -f 1 "${base_dir}/Data/CDS/${element}.fa" \
    > "${base_dir}/Data/Protein/${element}.fa"

    rm -f ${base_dir}/Data/Nucleotide/${element}.fa.index*
    ProgressBar "$state" "$total_states"
  done

  echo ""
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

working_directory=""
threads="8"
overwrite=false

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--workingDirectory)  shift; working_directory="$1" ;;
    -t|--threads)           shift; threads="$1" ;;
    --overwrite)            overwrite=true ;;
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
echo "  Working directory:               ${working_directory}"
echo "  Threads:                         ${threads}"
echo "  Overwrite previous run:          ${overwrite}"
echo ""

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Checking arguments and input files..."

if [[ -z "$working_directory" ]]; then
  echo "Error: missing required argument '-w / --workingDirectory'." >&2
  print_help; exit 1
fi
if [[ ! -d "$working_directory" ]]; then
  echo "Error: directory '${working_directory}' does not exist." >&2
  exit 1
fi
working_directory=$(realpath "$working_directory")
check_directory_structure "${working_directory}"

check_threads "$threads"

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

if $overwrite; then
  log_info "Removing previous run if exists..."
  overwrite "${working_directory}" "$(basename -s .sh "$0")"
fi

log_info "Organizing workspace..."
organize_working_directory "${working_directory}"

log_info "Step 1: All-vs-All BLASTN..."
nucleotide_match "${working_directory}"
log_info "-> Step 1 finished. Proceeding."

log_info "Step 2: Lifting gene models..."
run_liftoff "${working_directory}"
log_info "-> Step 2 finished. Proceeding."

log_info "Step 3: organizing results..."
store_old_gene_data "${working_directory}"
organize_information "${working_directory}"
log_info "-> Step 3 finished. Proceeding."

rm -r "${working_directory}/Workspace/$(basename -s .sh "$0")/temp" > /dev/null
log_info "$(basename -s .sh "$0") finished."
