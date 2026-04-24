#!/usr/bin/env bash
# ==============================================================================
# NAME:        OrthogroupsOverrepresentation.sh
# DESCRIPTION: Identifies orthogroups overrepresented in a dataset.
#              Executes three main steps:
#              1. OrthoFinder run with DIAMOND ultra-sensitive mode.
#              2. Removal of orthogroups associated with captain genes.
#              3. Overrepresentation analysis in one of two modes:
#                 Outliers   - Flags orthogroups with abnormal copy numbers
#                              using IQR-based upper fence (Standard or Skew).
#                 Enrichment - Identifies orthogroups enriched in a metadata
#                              group using one-sided Fisher's exact test.
# USAGE:       StarCrew OrthogroupsOverrepresentation [options]
#              StarCrew OrthogroupsOverrepresentation -help
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        24/Apr/2026
# VERSION:     1.0.0
# ==============================================================================

set -uo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

readonly SCRIPT_DIR="$( dirname -- "$( readlink -f -- "$0" )" )"
readonly LIB_DIR="${SCRIPT_DIR}/../lib"
readonly AUXILIARY_DIR="${SCRIPT_DIR}/../aux"

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

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays the help message with all available options.
# Arguments:
#   None
# Returns:
#   0 always
print_help() {
  echo "Command to identify overrepresented orthogroups in a dataset."
  echo "Performs three main steps:"
  echo "  1. OrthoFinder with DIAMOND ultra-sensitive mode."
  echo "  2. Removal of captain-associated orthogroups."
  echo "  3. Overrepresentation analysis:"
  echo "     Outliers   - Flags orthogroups with abnormal copy numbers via"
  echo "                  skew IQR upper fence: Q3 + n*e^(4MC)*IQR)."
  echo "     Enrichment - Identifies orthogroups enriched in a metadata group"
  echo "                  using one-sided Fisher's exact test."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -w <directory_path>"
  echo "       [ -m <string> { -a <string> -c <float> -cm <string>"
  echo "         | -n <string> -v <string> -p <float> -pa <string> }"
  echo "         -t <integer> { --overwrite | --skip-orthofinder } ]"
  echo
  echo "Required args:"
  echo "  -w, --workingDirectory  Working directory where all data are stored."
  echo
  echo "Required args with defaults:"
  echo "  -m, --mode         Analysis mode (Default: Outliers) [Available: Outliers, Enrichment]."
  echo "  -s, --scoreMatrix  DIAMOND score matrix for OrthoFinder (Default: BLOSUM62)"
  echo "                     [Available: BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80,"
  echo "                      BLOSUM90, PAM250, PAM70, PAM30]."
  echo
  echo "Required args with defaults in 'Outliers' mode:"
  echo "  -c, --coefficient    Fence coefficient (Default: 1.5) [range: 1.5-3.0]."
  echo "  -cm, --countMode     Count type for outliers (Default: Gene) [Available: Gene, Ship]."
  echo
  echo "Required args in 'Enrichment' mode:"
  echo "  -n, --name   Metadata column name to use for grouping."
  echo "  -v, --value  Target value within the column to compare against the rest."
  echo
  echo "Required args with defaults in 'Enrichment' mode:"
  echo "  -p, --pValue  Significance threshold (Default: 0.05) [range: 0.001-0.1]."
  echo "  -pa, --pAdj   Multiple testing correction (Default: bonferroni)"
  echo "                [Available: BH, BY, bonferroni, fdr, hochberg, holm, hommel]."
  echo
  echo "Optional args:"
  echo "  -t, --threads          Threads (Default: 8)."
  echo "  --overwrite            Overwrite a previous run. Not compatible with '--skip-orthofinder'."
  echo "  --skip-orthofinder     Skip OrthoFinder and re-run only the analysis step."
  echo "                         Not compatible with '--overwrite'."
  echo "  -help                  Display this help message."
}

# Sets up the workspace by copying protein, GFF, and captain data.
# Sets global: captain_phylogeny
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 on missing directories/files or previous run
organize_working_directory() {
  local base_dir="$1"
  local data_dir="${base_dir}/Data"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  if [[ -d "$working_dir" ]]; then
    echo "Error: a previous run exists in '${working_dir}'." >&2
    echo "Use '--overwrite' or '--skip-orthofinder' to proceed." >&2
    exit 1
  fi

  local protein_dir gff_dir metadata_dir captain_dir
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  metadata_dir=$(find "$base_dir" -maxdepth 1 -type d -name "metadata_files" 2>/dev/null)
  captain_dir=$(find "$base_dir" -maxdepth 1 -type d -name "CaptainIdentification" 2>/dev/null)

  if [[ ! -d "$captain_dir" ]]; then
    echo "Error: CaptainIdentification directory not found in '${base_dir}'." >&2
    exit 1
  fi
  if [[ ! -f "${captain_dir}/CaptainsID.txt" ]]; then
    echo "Error: CaptainsID.txt not found in '${captain_dir}'." >&2
    exit 1
  fi

  mkdir -p "$working_dir" "$temp_dir"
  cp -r "$protein_dir" "$working_dir"
  cp -r "$gff_dir" "$working_dir"
  cp "${captain_dir}/CaptainsID.txt" "$working_dir"

  if [[ -f "${captain_dir}/CaptainPhylogeny.nw" ]]; then
    cp "${captain_dir}/CaptainPhylogeny.nw" "$working_dir"
    captain_phylogeny=true
  else
    captain_phylogeny=false
  fi

  if [[ "$mode" == "Enrichment" ]]; then
    if [[ -f "${metadata_dir}/metadata.csv" ]]; then
      cp "${metadata_dir}/metadata.csv" "$working_dir"
    else
      echo "Error: metadata.csv not found in '${metadata_dir}'." >&2
      exit 1
    fi
  fi
}

# Runs OrthoFinder with DIAMOND ultra-sensitive mode, optionally guided by the
# captain phylogeny. Sets global: orthofinder_flag
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 always; sets orthofinder_flag=false on failure
run_orthofinder() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local output_dir="${working_dir}/Orthofinder"
  local captain_phylogeny_file="${working_dir}/CaptainPhylogeny.nw"

  local protein_dir
  protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)

  local element_number
  element_number=$(ls "$protein_dir" | wc -l)

  if (( $(ulimit -Sn) < element_number * 2 )); then
    echo "Error: system file descriptor limit ($(ulimit -Sn)) is too low." >&2
    echo "Please raise it to at least '$((element_number * 2))' using 'ulimit -n'." >&2
    exit 1
  fi

  if $captain_phylogeny; then
    orthofinder \
      -a "$(( threads / 2 ))" -t "${threads}" \
      -f "${protein_dir}" -A mafft \
      -S diamond_ultra_sens --matrix "${score_matrix}" \
      -s "${captain_phylogeny_file}" --scores-v2 \
      -o "${output_dir}" -n characterization \
      &> "${working_dir}/orthofinder.log"
  else
    orthofinder \
      -a "$(( threads / 2 ))" -t "${threads}" \
      -f "${protein_dir}" -A mafft \
      -S diamond_ultra_sens --matrix "${score_matrix}" \
      --scores-v2 -o "${output_dir}" -n characterization \
      &> "${working_dir}/orthofinder.log"
  fi

  local results_path="${output_dir}/Results_characterization/Orthogroups/Orthogroups.GeneCount.tsv"

  if [[ -f "$results_path" ]]; then
    orthofinder_flag=true
    cp "$results_path" "$working_dir"

    awk 'BEGIN {FS=OFS="\t"}
      NR==1 { TOTAL_COLUMNS=NF; print $0; next }
      {
        ROW_TOTAL=0
        sub(/[[:space:]]+$/, "", $0)
        printf "%s", $1
        for (i=2; i<=TOTAL_COLUMNS; i++) {
          is_present = (i<=NF) ? (($i!="") ? 1 : 0) : 0
          ROW_TOTAL += is_present
          printf "%s%d", OFS, is_present
        }
        printf "%s%d\n", OFS, ROW_TOTAL
      }' \
      "${output_dir}/Results_characterization/Orthogroups/Orthogroups_UnassignedGenes.tsv" \
      | sed '1d' >> "${working_dir}/Orthogroups.GeneCount.tsv"

    cp "${output_dir}/Results_characterization/Orthogroups/Orthogroups.txt" \
       "$working_dir"
  else
    orthofinder_flag=false
  fi
}

# Filters out orthogroups containing captain genes from the OrthoFinder output.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if required files are missing
remove_captain_orthogroups() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local count_file="${working_dir}/Orthogroups.GeneCount.tsv"
  local gene_file="${working_dir}/Orthogroups.txt"
  local captain_file="${working_dir}/CaptainsID.txt"

  if [[ ! -f "$count_file" || ! -f "$gene_file" ]]; then
    echo "Error: required OrthoFinder output files are missing." >&2
    exit 1
  fi
  if [[ ! -f "$captain_file" ]]; then
    echo "Error: CaptainsID.txt is missing." >&2
    exit 1
  fi

  grep -w -v -f "$captain_file" "$gene_file" \
    > "${working_dir}/Orthogroups-captainless.txt"
  awk -F ':' '{print $1}' "${working_dir}/Orthogroups-captainless.txt" \
    > "${temp_dir}/Orthogroups-captainlessID.txt"
  head -n1 "$count_file" > "${working_dir}/Orthogroups.GeneCount-captainless.tsv"
  grep -w -f "${temp_dir}/Orthogroups-captainlessID.txt" "$count_file" \
    >> "${working_dir}/Orthogroups.GeneCount-captainless.tsv"
}

# Copies analysis results from the workspace to the final output directory.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if workspace directory is missing
organize_information() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local orthogroups_dir="${working_dir}/Orthofinder/Results_characterization/Orthogroup_Sequences"
  local output_dir="${base_dir}/$(basename -s .sh "$0")"

  if [[ ! -d "$working_dir" ]]; then
    echo "Error: workspace directory '${working_dir}' not found." >&2
    exit 1
  fi

  mkdir -p "$output_dir"

  [[ -f "${working_dir}/OrthogroupsSizeHistogram.svg" ]] && \
    cp "${working_dir}/OrthogroupsSizeHistogram.svg" "${output_dir}/"

  if [[ -f "${working_dir}/Overrepresented_orthogroups.txt" ]]; then
    cp "${working_dir}/Overrepresented_orthogroups.txt" "${output_dir}/"
    mkdir -p "${output_dir}/Orthogroups"
    while read -r orthogroup; do
      cp "${orthogroups_dir}/${orthogroup}.fa" "${output_dir}/Orthogroups/"
    done < <(awk '{print $1}' "${working_dir}/Overrepresented_orthogroups.txt")
  fi

  [[ -f "${working_dir}/Enrichment_results.txt" ]] && \
    cp "${working_dir}/Enrichment_results.txt" "${output_dir}/"

  if [[ -f "${working_dir}/Enrich_orthogroups.txt" ]]; then
    cp "${working_dir}/Enrich_orthogroups.txt" "${output_dir}/"
    mkdir -p "${output_dir}/Orthogroups"
    while read -r orthogroup; do
      cp "${orthogroups_dir}/${orthogroup}.fa" "${output_dir}/Orthogroups/"
    done < <(awk '{print $1}' "${working_dir}/Enrich_orthogroups.txt")
  fi
}

# Removes previous analysis outputs while preserving the OrthoFinder workspace.
# Used by --skip-orthofinder to allow re-running only the analysis step.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if metadata is missing in Enrichment mode
clean_previous_run() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local output_dir="${base_dir}/$(basename -s .sh "$0")"
  local metadata_dir
  metadata_dir=$(find "$base_dir" -maxdepth 1 -type d -name "metadata_files" 2>/dev/null)

  rm -rf "$output_dir"
  rm -f "${working_dir}/OrthogroupsSizeHistogram.svg" \
        "${working_dir}/Overrepresented_orthogroups.txt" \
        "${working_dir}/Enrichment_results.txt" \
        "${working_dir}/Enrich_orthogroups.txt"

  if [[ "$mode" == "Enrichment" ]]; then
    if [[ -f "${metadata_dir}/metadata.csv" ]]; then
      cp "${metadata_dir}/metadata.csv" "$working_dir"
    else
      echo "Error: metadata.csv not found in '${metadata_dir}'." >&2
      exit 1
    fi
  fi
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

working_directory=""
mode="Outliers"
score_matrix="BLOSUM62"
coefficient="1.5"
count_mode="Gene"
column=""
value=""
pvalue="0.05"
padjust="bonferroni"
threads="8"
overwrite=false
skip_orthofinder=false

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--workingDirectory)  shift; working_directory="$1" ;;
    -m|--mode)              shift; mode="$1" ;;
    -s|--scoreMatrix)       shift; score_matrix="$1" ;;
    -c|--coefficient)       shift; coefficient="$1" ;;
    -cm|--countMode)        shift; count_mode="$1" ;;
    -n|--name)              shift; column="$1" ;;
    -v|--value)             shift; value="$1" ;;
    -p|--pValue)            shift; pvalue="$1" ;;
    -pa|--pAdj)             shift; padjust="$1" ;;
    -t|--threads)           shift; threads="$1" ;;
    --overwrite)            overwrite=true ;;
    --skip-orthofinder)     skip_orthofinder=true ;;
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
echo "  Mode:                            ${mode}"
echo "  DIAMOND score matrix:            ${score_matrix}"
echo "  Fence coefficient (Outliers):    ${coefficient}"
echo "  Count mode (Outliers):           ${count_mode}"
echo "  Metadata column (Enrichment):    ${column}"
echo "  Target value (Enrichment):       ${value}"
echo "  p-value (Enrichment):            ${pvalue}"
echo "  Multiple testing correction:     ${padjust}"
echo "  Threads:                         ${threads}"
echo "  Overwrite previous run:          ${overwrite}"
echo "  Skip OrthoFinder:                ${skip_orthofinder}"
echo ""

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Checking arguments and input files..."

check_mode_parameter "$mode" "$(basename -s .sh "$0")"

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

valid_matrices="BLOSUM45 BLOSUM50 BLOSUM62 BLOSUM80 BLOSUM90 PAM250 PAM70 PAM30"
if [[ ! " ${valid_matrices} " =~ " ${score_matrix} " ]]; then
  echo "Error: score matrix '${score_matrix}' is not a valid DIAMOND matrix." >&2
  print_help; exit 1
fi

if [[ "$mode" == "Enrichment" ]]; then
  if [[ -z "$column" || -z "$value" ]]; then
    echo "Error: '-n / --name' and '-v / --value' are required in Enrichment mode." >&2
    print_help; exit 1
  fi

  if [[ "$pvalue" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    if (( $(echo "$pvalue < 0.001" | bc -l) )) || \
       (( $(echo "$pvalue > 0.1" | bc -l) )); then
      echo "Error: p-value '${pvalue}' is out of range [0.001-0.1]." >&2
      print_help; exit 1
    fi
  else
    echo "Error: p-value '${pvalue}' is not a valid float." >&2
    print_help; exit 1
  fi

  valid_padjust="BH BY bonferroni fdr hochberg holm hommel"
  if [[ ! " ${valid_padjust} " =~ " ${padjust} " ]]; then
    echo "Error: correction method '${padjust}' is not valid." >&2
    print_help; exit 1
  fi

  if [[ ! -f "${working_directory}/metadata_files/metadata.csv" ]]; then
    echo "Error: metadata.csv not found in '${working_directory}/metadata_files/'." >&2
    exit 1
  fi

  column_number=$(sed -n "1 s/${column}.*//p" \
    "${working_directory}/metadata_files/metadata.csv" \
    | sed 's/[^;*]//g' | wc -c)

  if (( column_number == 0 )); then
    echo "Error: column '${column}' does not exist in metadata." >&2
    exit 1
  elif (( column_number <= 2 )); then
    echo "Error: the target column cannot be the ID column." >&2
    exit 1
  fi

  if [[ "$value" == "NA" ]]; then
    echo "Error: 'Enrichment' mode does not accept 'NA' as a target value." >&2
    exit 1
  fi

  local_temp="${working_directory}/temp_dataset"
  find "${working_directory}/Data/Protein" -maxdepth 1 -type f -name "*.fa" \
    | xargs -n1 basename -s .fa > "$local_temp"

  value_number=$(grep -f "$local_temp" \
    "${working_directory}/metadata_files/metadata.csv" \
    | cut -d ";" -f "$column_number" \
    | grep -c "$value")
  rm -f "$local_temp"

  if (( value_number == 0 )); then
    echo "Error: value '${value}' not found in column '${column}'." >&2
    exit 1
  elif (( value_number < 5 )); then
    echo "Error: value '${value}' has only ${value_number} elements — too few for enrichment." >&2
    exit 1
  fi
  log_info "'${value_number}' elements are in the ingroup for enrichment analysis."
fi

if [[ "$mode" == "Outliers" ]]; then
  if [[ "$count_mode" != "Gene" && "$count_mode" != "Ship" ]]; then
    echo "Error: count mode '${count_mode}' is not valid." >&2
    print_help; exit 1
  fi

  if [[ "$coefficient" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    if (( $(echo "$coefficient < 1.5" | bc -l) )) || \
       (( $(echo "$coefficient > 3.0" | bc -l) )); then
      echo "Error: coefficient '${coefficient}' is out of range [1.5-3.0]." >&2
      print_help; exit 1
    fi
  else
    echo "Error: coefficient '${coefficient}' is not a valid float." >&2
    print_help; exit 1
  fi
fi

if $overwrite && $skip_orthofinder; then
  echo "Error: '--overwrite' and '--skip-orthofinder' are mutually exclusive." >&2
  exit 1
fi

check_threads "$threads"

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

if $overwrite; then
  log_info "Removing previous run if exists..."
  overwrite "${working_directory}" "$(basename -s .sh "$0")" "${mode}"
fi

check_directory_structure "$working_directory"

if ! $skip_orthofinder; then
  log_info "Organizing workspace..."
  organize_working_directory "${working_directory}"

  log_info "Step 1: Running OrthoFinder..."
  run_orthofinder "${working_directory}"
  log_info "-> Step 1 finished. Proceeding."
else
  log_info "Skipping workspace organization and OrthoFinder run."
  clean_previous_run "${working_directory}"
fi

log_info "Step 2: Removing captain-associated orthogroups..."
remove_captain_orthogroups "${working_directory}"
log_info "-> Step 2 finished. Proceeding."

log_info "Step 3: Running overrepresentation analysis in '${mode}' mode..."
Rscript "${AUXILIARY_DIR}/OverrepresentationAnalysis.R" \
  -d "${working_directory}/Workspace/$(basename -s .sh "$0")" \
  -m "${mode}" \
  -c "${coefficient}" \
  -n "${column}" \
  -v "${value}" \
  -p "${pvalue}" \
  --padjust "${padjust}" \
  --countmode "${count_mode}"
log_info "-> Step 3 finished. Proceeding."

log_info "Organizing output..."
organize_information "${working_directory}"
log_info "OrthogroupsOverrepresentation finished."
