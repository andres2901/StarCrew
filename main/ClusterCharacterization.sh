#!/usr/bin/env bash
# ==============================================================================
# NAME:        ClusterCharacterization.sh
# DESCRIPTION: Characterizes each element cluster through eight main steps:
#              1. Orthogroup identification via OrthoFinder.
#              2. All-vs-all BLASTN for nucleotide synteny visualization.
#              3. Hierarchical clustering based on orthogroup gene counts.
#              4. Full cluster connectivity, cargo heatmap, and synteny figure.
#              5. Detection of individual nesting events within the cluster.
#              6. Core gene identification (general and subcluster-specific).
#              7. Putative cargo movement event detection across subclusters.
#              8. Discordance detection between cargo clustering and captain tree.
# USAGE:       StarCrew ClusterCharacterization [options]
#              StarCrew ClusterCharacterization -help
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
  echo "Command to characterize each element cluster."
  echo "Performs eight steps per cluster:"
  echo "  1. Orthogroup identification via OrthoFinder."
  echo "  2. All-vs-all BLASTN for nucleotide synteny visualization."
  echo "  3. Hierarchical clustering based on orthogroup gene counts."
  echo "  4. Full cluster connectivity, cargo heatmap, and synteny figure."
  echo "  5. Detection of individual nesting events."
  echo "  6. Core gene identification:"
  echo "     6.1. General: orthogroups present in >= 80% of elements."
  echo "     6.2. Specific: orthogroups present in >= 80% of elements"
  echo "          within subclusters (>= 5 elements)."
  echo "  7. Putative cargo movement event detection across subclusters."
  echo "  8. Discordance detection between cargo clustering and captain tree."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -w <directory_path>"
  echo "       [ -l <integer> -i <float> -t <integer> --overwrite ]"
  echo
  echo "Required args:"
  echo "  -w, --workingDirectory  Working directory where all data are stored."
  echo
  echo "Required args with defaults:"
  echo "  -l, --length    Minimum BLAST alignment length for synteny visualization"
  echo "                  (Default: 1000) [range: 200-5000]."
  echo "  -i, --identity  Minimum BLAST identity % for synteny visualization"
  echo "                  (Default: 70.0) [range: 50.0-95.0]."
  echo
  echo "Optional args:"
  echo "  -t, --threads   Threads for OrthoFinder and BLAST (Default: 8)."
  echo "  --overwrite     Overwrite a previous run (Default: off)."
  echo "  -help           Display this help message."
}

# Validates that the ClustersAnalyzed.txt file exists and reports cluster count.
# Sets global: cluster_number
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if the file is missing
check_clusters() {
  local base_dir="$1"
  local cluster_information="${base_dir}/Clusters/ClustersAnalyzed.txt"

  if [[ ! -f "$cluster_information" ]]; then
    echo "Error: 'ClustersAnalyzed.txt' not found in '${base_dir}/Clusters/'." >&2
    exit 1
  fi

  cluster_number=$(wc -l < "$cluster_information")
  log_info "Analyzing '${cluster_number}' clusters.\n"
}

# Validates that a captain phylogeny file exists for the given cluster directory.
# Sets global: captainremoval_number
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 on success, exits with 1 if phylogeny file is missing
check_captain_information() {
  local base_dir="$1"
  local captain_phylogeny="${base_dir}/CaptainIdentification/CaptainPhylogeny.nw"

  if [[ ! -f "$captain_phylogeny" ]]; then
    echo "Error: captain phylogeny not found in '${base_dir}/CaptainIdentification'." >&2
    exit 1
  fi

  local captainremoval_dir
  captainremoval_dir=$(find "${base_dir}" -maxdepth 1 \
    -type d -name "Captainless_elements" 2>/dev/null)

  if [[ -z "$captainremoval_dir" ]]; then
    captainremoval_number="0"
  else
    captainremoval_number=$(find "$captainremoval_dir" \
      -maxdepth 1 -type f -name "*.gff" | wc -l)
  fi
}

# Sets up the workspace for a cluster by copying data directories and building
# a merged GFF file.
# Arguments:
#   $1 - cluster base directory path
#   $2 - path to the sub_clusters.txt file
#   $3 - cluster ID used to filter subcluster entries
# Returns:
#   0 on success, exits with 1 if a previous workspace run exists
organize_working_directory() {
  local base_dir="$1"
  local subcluster_file="$2"
  local subcluster_id="$3"

  local data_dir="${base_dir}/Data"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local full_gff="${working_dir}/Final_model.gff"
  local captain_phylogeny="${base_dir}/CaptainIdentification/CaptainPhylogeny.nw"

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
  cp -r "$gff_dir" "$working_dir"
  cp -r "$nucleotide_dir" "$working_dir"
  cp -r "$protein_dir" "$working_dir"
  cp -r "$cds_dir" "$working_dir"
  cp "$captain_phylogeny" "$working_dir"

  grep "${subcluster_id}" "$subcluster_file" \
    > "${working_dir}/Subclusters_MCL.txt"

  echo "##gff-version 3" > "$full_gff"
  grep -v "#" "${gff_dir}"/* >> "$full_gff"
}

# Runs OrthoFinder on the cluster protein directory guided by the captain
# phylogeny. Sets global: orthofinder_flag
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 always; sets orthofinder_flag=false on failure
run_orthofinder() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local output_dir="${working_dir}/Orthofinder"
  local captain_phylogeny="${working_dir}/CaptainPhylogeny.nw"

  local protein_dir
  protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)

  local element_number
  element_number=$(ls "$protein_dir" | wc -l)

  if (( $(ulimit -Sn) < element_number * 3 )); then
    log_warning "System open file limit per process may be too low for OrthoFinder." \
      "Consider raising it to at least '$((element_number * 3))'" \
      "using 'ulimit -n'."
  fi

  orthofinder \
    -a "$(( threads / 2 ))" -t "${threads}" \
    -f "${protein_dir}" -A mafft -S diamond \
    --matrix PAM30 -s "${captain_phylogeny}" \
    --scores-v2 -o "${output_dir}" -n characterization \
    &> "${working_dir}/orthofinder.log"

  local results_path="${output_dir}/Results_characterization/Orthogroups/Orthogroups.GeneCount.tsv"

  if [[ -f "$results_path" ]]; then
    orthofinder_flag=true
    cp "$results_path" "$working_dir"

    # Append singleton genes as binary presence/absence rows
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

# Runs all-vs-all BLASTN within the cluster and filters results by length.
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 on success
run_blast() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local nucleotide_dir
  nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)

  cat "${nucleotide_dir}"/*.fa > "${temp_dir}/sequence.fasta"

  makeblastdb -dbtype nucl \
    -in "${temp_dir}/sequence.fasta" \
    -out "${temp_dir}/Cluster" &> /dev/null

  blastn \
    -query "${temp_dir}/sequence.fasta" \
    -db "${temp_dir}/Cluster" \
    -evalue 1e-60 \
    -num_threads "${threads}" \
    -outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen" \
    -task blastn \
    -gapopen 8 -gapextend 6 -reward 5 -penalty -4 \
    -out "${temp_dir}/blastresults.txt"

  awk -v min="$blast_length" \
    'BEGIN{FS=OFS="\t"} {if($8>=min) print}' \
    "${temp_dir}/blastresults.txt" \
    > "${working_dir}/Blast_CleanResults.txt"
}

# Checks for identified core gene files and copies orthogroup sequences.
# Appends the cluster ID to the global ClusterCore.txt if core genes are found.
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 always
check_core() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local orthogroups_dir="${working_dir}/Orthofinder/Results_characterization/Orthogroup_Sequences"
  local temp_dir="${working_dir}/temp"

  ls "${working_dir}/Core_genes"*.txt > "${temp_dir}/core_files.txt" 2>/dev/null

  local file_number
  file_number=$(wc -l < "${temp_dir}/core_files.txt")

  if (( file_number > 0 )); then
    log_step "Core genes identified. Processing..."

    grep -w "${cluster_id}" "${base_dir}/../ClustersAnalyzed.txt" \
      >> "${base_dir}/../ClusterCore.txt"

    mkdir -p "${working_dir}/Core_genes"
    cat "${working_dir}/Core_genes"*.txt | sort -u > "${temp_dir}/Full_core.txt"

    while read -r orthogroup; do
      cp "${orthogroups_dir}/${orthogroup}.fa" "${working_dir}/Core_genes/"
    done < "${temp_dir}/Full_core.txt"
  fi
}

# Checks for cargo movement files and copies orthogroup sequences per subcluster.
# Appends the cluster ID to the global ClusterMovement.txt if movement is found.
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 always
check_movement() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local orthogroups_dir="${working_dir}/Orthofinder/Results_characterization/Orthogroup_Sequences"
  local temp_dir="${working_dir}/temp"

  ls "${working_dir}/"*_moveOrthologs.txt \
    > "${temp_dir}/movement_files.txt" 2>/dev/null

  local file_number
  file_number=$(wc -l < "${temp_dir}/movement_files.txt")

  if (( file_number > 0 )); then
    log_step "Cargo movement genes identified. Processing..."

    grep -w "${cluster_id}" "${base_dir}/../ClustersAnalyzed.txt" \
      >> "${base_dir}/../ClusterMovement.txt"

    xargs -n1 basename -s .txt < "${temp_dir}/movement_files.txt" \
      | while read -r subcluster; do
          mkdir -p "${working_dir}/${subcluster}"
          while read -r orthogroup; do
            cp "${orthogroups_dir}/${orthogroup}.fa" \
               "${working_dir}/${subcluster}/"
          done < "${working_dir}/${subcluster}.txt"
        done
  fi
}

# Copies final characterization outputs from the workspace to the cluster's
# output directory.
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 always
organize_information() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local output_dir="${base_dir}/$(basename -s .sh "$0")"

  mkdir -p "$output_dir"

  if [[ ! -d "$working_dir" ]]; then
    echo "Error: workspace directory '${working_dir}' not found." >&2
    return 1
  fi

  local tree_count
  tree_count=$(find "$working_dir" -name "CargoHierarchicalTree*.nwk" \
    2>/dev/null | wc -l)

  if (( tree_count >= 1 )); then
    cp "${working_dir}/CargoHierarchicalTree"*.nwk "${output_dir}/"
  else
    log_warning "No cargo hierarchical tree found for this cluster."
  fi

  [[ -d "${working_dir}/Figures" ]] && \
    cp -r "${working_dir}/Figures" "${output_dir}/"

  if [[ -d "${working_dir}/Core_genes" ]]; then
    mkdir -p "${output_dir}/Core"
    cp -r "${working_dir}/Core_genes" "${output_dir}/Core/Orthogroups"
    cp "${working_dir}/Core_genes-"* "${output_dir}/Core/"
  fi

  local movement_count
  movement_count=$(ls "${working_dir}/"*_moveOrthologs.txt 2>/dev/null | wc -l)

  if (( movement_count >= 1 )); then
    mkdir -p "${output_dir}/Movement_genes"
    ls "${working_dir}/"*_moveOrthologs.txt \
      | xargs -n1 basename -s .txt \
      | while read -r line; do
          local subcluster
          subcluster=$(echo "$line" | awk -F '_' '{print $1}')
          mkdir -p "${output_dir}/Movement_genes/${subcluster}"
          cp -r "${working_dir}/${subcluster}_moveOrthologs" \
                "${output_dir}/Movement_genes/${subcluster}/Orthogroups"
          cp "${working_dir}/${subcluster}_moveOrthologs.txt" \
             "${output_dir}/Movement_genes/${subcluster}/OrthogroupsID.txt"
          cp "${working_dir}/${subcluster}_moveOrthologsTable.csv" \
             "${output_dir}/Movement_genes/${subcluster}/Matrix.txt"
          cp "${working_dir}/${subcluster}"*_moveOrthologsMatrix.csv \
             "${output_dir}/Movement_genes/${subcluster}/" 2>/dev/null || true
        done
  fi

  [[ -f "${working_dir}/Discordant_elements.txt" ]] && \
    cp "${working_dir}/Discordant_elements.txt" "${output_dir}/"
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

working_directory=""
threads="8"
blast_length="1000"
identity="70"
overwrite=false

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--workingDirectory)  shift; working_directory="$1" ;;
    -t|--threads)           shift; threads="$1" ;;
    -l|--length)            shift; blast_length="$1" ;;
    -i|--identity)          shift; identity="$1" ;;
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
echo "  Minimum alignment length:        ${blast_length}"
echo "  Minimum identity (%):            ${identity}"
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

if [[ "$blast_length" =~ ^[0-9]+$ ]]; then
  if (( blast_length < 200 || blast_length > 5000 )); then
    echo "Error: length '${blast_length}' is out of range [200-5000]." >&2
    print_help; exit 1
  fi
else
  echo "Error: length '${blast_length}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$identity" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
  if (( $(echo "$identity < 50.0" | bc -l) )) || \
     (( $(echo "$identity > 95.0" | bc -l) )); then
    echo "Error: identity '${identity}' is out of range [50.0-95.0]." >&2
    print_help; exit 1
  fi
else
  echo "Error: identity '${identity}' is not a valid float." >&2
  print_help; exit 1
fi

check_threads "$threads"

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

check_clusters "${working_directory}"

# Remove per-run control files when overwriting
if $overwrite; then
  rm -f "${working_directory}/Clusters/ClusterOrthogroups.txt" \
        "${working_directory}/Clusters/ClusterCore.txt" \
        "${working_directory}/Clusters/ClusterMovement.txt"
fi

# Iterate over each cluster
while read -r cluster_id; do
  internal_dir="${working_directory}/Clusters/${cluster_id}"
  subcluster_number=$(
    grep -w "${cluster_id}" \
      "${working_directory}/Clusters/ClustersAnalyzed.txt" \
      | awk '{print $3}' \
      | sed $'s/[^[:print:]\t]//g'
  )

  log_info "Analyzing cluster '${cluster_id}'."
  check_directory_structure "${internal_dir}"
  check_captain_information "${internal_dir}"

  if $overwrite; then
    log_info "Removing previous run for cluster '${cluster_id}' if exists..."
    overwrite "${internal_dir}" "$(basename -s .sh "$0")"
  fi

  log_info "Organizing workspace..."
  organize_working_directory \
    "${internal_dir}" \
    "${working_directory}/Clusters/sub_clusters.txt" \
    "${cluster_id}"

  log_info "Step 1: Running OrthoFinder..."
  run_orthofinder "${internal_dir}"
  log_step "Checking OrthoFinder results..."
  if ! $orthofinder_flag; then
    log_warning "OrthoFinder failed for cluster '${cluster_id}'. Skipping."
    continue
  fi
  grep -w "${cluster_id}" \
    "${working_directory}/Clusters/ClustersAnalyzed.txt" \
    >> "${working_directory}/Clusters/ClusterOrthogroups.txt"
  log_info "-> Step 1 finished. Proceeding."

  log_info "Step 2: Running BLASTN for nucleotide synteny..."
  run_blast "${internal_dir}"
  log_info "-> Step 2 finished. Proceeding."

  log_info "Step 3: Running cluster characterization..."
  Rscript "${AUXILIARY_DIR}/ClusterAnalysis.R" \
    -d "${internal_dir}/Workspace/$(basename -s .sh "$0")/" \
    -s "${subcluster_number}" \
    -c "${captainremoval_number}" \
    -p "${identity}"

  mkdir -p "${internal_dir}/Workspace/$(basename -s .sh "$0")/Figures"
  mv "${internal_dir}/Workspace/$(basename -s .sh "$0")"/*.svg \
     "${internal_dir}/Workspace/$(basename -s .sh "$0")/Figures/" 2>/dev/null || true

  log_step "Checking for core genes..."
  check_core "${internal_dir}"

  log_step "Checking for cargo movement genes..."
  check_movement "${internal_dir}"
  log_info "-> Step 3 finished. Proceeding."

  log_info "Organizing output for cluster '${cluster_id}'..."
  organize_information "${internal_dir}"
  rm -r "${internal_dir}/Workspace/$(basename -s .sh "$0")/temp" > /dev/null
  log_info "Cluster '${cluster_id}' complete.\n"

done < <(
  awk '{print $1}' "${working_directory}/Clusters/ClustersAnalyzed.txt" \
    | sed $'s/[^[:print:]\t]//g'
)

log_info "All clusters have been analyzed."
