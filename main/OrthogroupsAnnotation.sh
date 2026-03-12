#!/usr/bin/env bash
# ==============================================================================
# NAME:        OrthogroupsAnnotation.sh
# DESCRIPTION: Functional annotation of orthogroups selected by mode.
#              Executes four main steps:
#              1. Workspace organization and orthogroup selection by mode:
#                 Core, MoveAssociated, All, or Overrepresented.
#              2. Protein characterization using three approaches:
#                 2.1. InterProScan (CDD, Gene3D, HAMAP, PANTHER, Pfam, etc.)
#                 2.2. Foldseek (3D-structure-based homology search).
#                 2.3. hhblits (domain search against PfamA).
#              3. Per-orthogroup internal summary (CSV per orthogroup).
#              4. General summary across orthogroups retaining annotations
#                 shared by >= 50% of proteins within each orthogroup.
# USAGE:       StarCrew OrthogroupsAnnotation [options]
#              StarCrew OrthogroupsAnnotation -help
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
readonly DATABASE_DIR="${SCRIPT_DIR}/../databases"
readonly INTERPRO_DIR="${SCRIPT_DIR}/../interproscan"

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

check_required_software "$(basename -s .sh "$0")"

if [[ ! -d "$DATABASE_DIR" ]]; then
  echo "Error: database directory '${DATABASE_DIR}' does not exist." >&2
  exit 1
fi
check_databases "$DATABASE_DIR" "$(basename -s .sh "$0")"

readonly FOLDSEEK_PATH="${DATABASE_DIR}/Foldseek"
readonly HHSUITE_PATH="${DATABASE_DIR}/hhsuite"

check_interpro_software "$INTERPRO_DIR"
readonly INTERPRO_PATH="$(realpath "$INTERPRO_DIR")"

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays the help message with all available options.
# Arguments:
#   None
# Returns:
#   0 always
print_help() {
  echo "Command to run functional annotation for orthogroups."
  echo "Performs four steps:"
  echo "  1. Orthogroup selection by mode:"
  echo "     Core           - Core orthogroups from ClusterCharacterization."
  echo "     MoveAssociated - Orthogroups linked to cargo movement events."
  echo "     All            - All orthogroups from ClusterCharacterization."
  echo "     Overrepresented- Orthogroups from OrthogroupsOverrepresentation."
  echo "  2. Protein characterization:"
  echo "     2.1. InterProScan (CDD, Gene3D, HAMAP, PANTHER, Pfam, PIRSF,"
  echo "          PRINTS, PROSITEPATTERNS, PROSITEPROFILES, SFLD, SMART,"
  echo "          SUPERFAMILY, TIGRFAM)."
  echo "     2.2. Foldseek (3D-structure-based homology search)."
  echo "     2.3. hhblits (PfamA domain search)."
  echo "  3. Per-orthogroup internal summary CSV."
  echo "  4. General summary retaining annotations shared by >= 50% of proteins."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -w <directory_path>"
  echo "       [ -m <string> -f <string> -t <integer> --overwrite ]"
  echo
  echo "Required args:"
  echo "  -w, --workingDirectory  Working directory where all data are stored."
  echo
  echo "Required args with defaults:"
  echo "  -m, --mode       Orthogroups to annotate (Default: All)"
  echo "                   [Available: Core, MoveAssociated, All, Overrepresented]."
  echo "  -f, --foldseekdb Foldseek database (Default: afdb_swissprot)"
  echo "                   [Available: pdb, afdb_swissprot]."
  echo "  -t, --threads    Threads for all analyses (Default: 8)."
  echo
  echo "Optional args:"
  echo "  --overwrite  Overwrite a previous run (Default: off)."
  echo "  -help        Display this help message."
}

# Validates the presence of the cluster control file for the selected mode and
# sets the global clusters_file path.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if the required file is missing
check_clusters() {
  local base_dir="$1"
  local cluster_dir
  cluster_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Clusters" 2>/dev/null)

  case "$mode" in
    All)
      clusters_file="${cluster_dir}/ClusterOrthogroups.txt"
      ;;
    MoveAssociated)
      clusters_file="${cluster_dir}/ClusterMovement.txt"
      ;;
    Core)
      clusters_file="${cluster_dir}/ClusterCore.txt"
      ;;
  esac

  if [[ ! -f "$clusters_file" ]]; then
    echo "Error: '$(basename "$clusters_file")' not found in '${cluster_dir}'." >&2
    exit 1
  fi

  local cluster_number
  cluster_number=$(wc -l < "$clusters_file")
  log_info "Analyzing '${cluster_number}' clusters.\n"
}

# Validates that OrthogroupsOverrepresentation output exists for Overrepresented mode.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 0 (no work to do) or 1 on error
check_overrepresentation() {
  local base_dir="$1"
  local overrepresentation_dir
  overrepresentation_dir=$(find "$base_dir" -maxdepth 1 \
    -type d -name "OrthogroupsOverrepresentation" 2>/dev/null)

  if [[ ! -d "$overrepresentation_dir" ]]; then
    echo "Error: 'OrthogroupsOverrepresentation' directory not found." >&2
    echo "Run 'OrthogroupsOverrepresentation' before this command." >&2
    exit 1
  fi

  if [[ ! -d "${overrepresentation_dir}/Orthogroups" ]]; then
    log_info "No overrepresented orthogroups found in this dataset. Nothing to do."
    exit 0
  fi
}

# Checks that the required ClusterCharacterization workspace subdirectories
# exist for the selected mode. Sets global: directory_flag
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 always; sets directory_flag=false if required directories are missing
check_internal_directory_structure() {
  local base_dir="$1"

  local workspace_dir
  workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)
  if [[ -z "$workspace_dir" ]]; then
    echo "Error: Workspace directory not found in '${base_dir}'." >&2
    exit 1
  fi

  local annotation_dir
  annotation_dir=$(find "$workspace_dir" -maxdepth 1 \
    -type d -name "ClusterCharacterization" 2>/dev/null)
  if [[ -z "$annotation_dir" ]]; then
    echo "Error: ClusterCharacterization directory not found in '${workspace_dir}'." >&2
    echo "Run 'ClusterCharacterization' before this command." >&2
    exit 1
  fi

  directory_flag=true

  case "$mode" in
    All)
      local orthogroups_dir
      orthogroups_dir=$(find "$annotation_dir" -maxdepth 3 \
        -type d -name "Orthogroup_Sequences" 2>/dev/null)
      if [[ -z "$orthogroups_dir" ]]; then
        echo "Error: Orthogroup_Sequences not found in '${annotation_dir}'." >&2
        directory_flag=false
      fi
      ;;
    MoveAssociated)
      local move_dir
      move_dir=$(find "$annotation_dir" -maxdepth 1 \
        -type d -name "*_moveOrthologs" 2>/dev/null)
      if [[ -z "$move_dir" ]]; then
        echo "Warning: No movement orthogroup directories found in '${annotation_dir}'." >&2
        directory_flag=false
      fi
      ;;
    Core)
      local core_dir
      core_dir=$(find "$annotation_dir" -maxdepth 1 \
        -type d -name "Core_genes" 2>/dev/null)
      if [[ -z "$core_dir" ]]; then
        echo "Error: Core_genes directory not found in '${annotation_dir}'." >&2
        directory_flag=false
      fi
      ;;
  esac
}

# Copies orthogroup FASTA files into the workspace based on the selected mode.
# Removes stop-codon asterisks from sequences after copying.
# Arguments:
#   $1 - cluster base directory path
# Returns:
#   0 on success, exits with 1 if a previous workspace run exists
organize_working_directory() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local annotation_dir="${base_dir}/Workspace/ClusterCharacterization"
  local orthogroups_dir="${working_dir}/Orthogroups"
  local temp_dir="${working_dir}/temp"

  if [[ -d "$working_dir" ]]; then
    echo "Error: a previous run exists in '${working_dir}'." >&2
    echo "Use '--overwrite' to overwrite it." >&2
    exit 1
  fi

  mkdir -p "$orthogroups_dir" "$temp_dir"

  case "$mode" in
    All)
      local data_dir
      data_dir=$(find "$annotation_dir" -maxdepth 3 \
        -type d -name "Orthogroup_Sequences" 2>/dev/null)
      cp "${data_dir}"/* "$orthogroups_dir/"
      ;;
    MoveAssociated)
      find "$annotation_dir" -maxdepth 1 -type d -name "*_moveOrthologs" \
        | while read -r move_dir; do
            cp "${move_dir}"/* "$orthogroups_dir/"
          done
      ;;
    Core)
      local core_dir
      core_dir=$(find "$annotation_dir" -maxdepth 1 \
        -type d -name "Core_genes" 2>/dev/null)
      cp "${core_dir}"/* "$orthogroups_dir/"
      ;;
  esac

  sed -i -e 's/\*//g' "${orthogroups_dir}"/*
}

# Copies overrepresented orthogroup FASTA files into the workspace.
# Removes stop-codon asterisks from sequences after copying.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if a previous workspace run exists
organize_working_directory_overrepresentation() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local orthogroups_dir="${working_dir}/Orthogroups"
  local temp_dir="${working_dir}/temp"

  if [[ -d "$working_dir" ]]; then
    echo "Error: a previous run exists in '${working_dir}'." >&2
    echo "Use '--overwrite' to overwrite it." >&2
    exit 1
  fi

  mkdir -p "$orthogroups_dir" "$temp_dir"

  local overrepresented_dir data_dir
  overrepresented_dir=$(find "$base_dir" -maxdepth 1 \
    -type d -name "OrthogroupsOverrepresentation" 2>/dev/null)
  data_dir=$(find "$overrepresented_dir" -maxdepth 3 \
    -type d -name "Orthogroups" 2>/dev/null)

  cp "${data_dir}"/* "$orthogroups_dir/"
  sed -i -e 's/\*//g' "${orthogroups_dir}"/*
}

# Runs Foldseek easy-search for each orthogroup against the selected database.
# Joins results with annotation index files and removes empty output files.
# Arguments:
#   $1 - cluster or working directory base path
# Returns:
#   0 on success
run_foldseek() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local orthogroups_dir="${working_dir}/Orthogroups"
  local temp_dir="${working_dir}/temp"
  local foldseek_results="${working_dir}/Foldseek"

  mkdir -p "$foldseek_results"

  local total_states state=0
  total_states=$(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" | wc -l)

  while read -r og_id; do
    state=$(( state + 1 ))

    foldseek easy-search \
      "${orthogroups_dir}/${og_id}.fa" \
      "${FOLDSEEK_PATH}/${foldseekdb}" \
      "${temp_dir}/${og_id}.m8" \
      "${temp_dir}/tmp" \
      --prostt5-model "${FOLDSEEK_PATH}/weights" \
      -e 0.001 -c 0.5 --cov-mode 0 -v 0 \
      --threads "${threads}" &> /dev/null

    if [[ "$foldseekdb" == "pdb" ]]; then
      awk 'BEGIN{FS=OFS="\t"}{split($2,a,"-");$2=a[1];print}' \
        "${temp_dir}/${og_id}.m8" > "${temp_dir}/${og_id}-2.m8"
      join -t $'\t' -i -1 2 -2 1 \
        <(sort -k2,2 "${temp_dir}/${og_id}-2.m8") \
        "${FOLDSEEK_PATH}/entries_update.idx" \
        | awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' \
        | sort -k1,1 | sed -e 's/;/,/g' \
        > "${foldseek_results}/${og_id}.m8"
    elif [[ "$foldseekdb" == "afdb_swissprot" ]]; then
      awk 'BEGIN{FS=OFS="\t"}{split($2,a,"-");$2=a[2];print}' \
        "${temp_dir}/${og_id}.m8" > "${temp_dir}/${og_id}-2.m8"
      join -t $'\t' -i -1 2 -2 1 \
        <(sort -k2,2 "${temp_dir}/${og_id}-2.m8") \
        "${FOLDSEEK_PATH}/Accession_swissprot.txt" \
        | awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' \
        | sort -k1,1 | sed -e 's/;/,/g' \
        > "${foldseek_results}/${og_id}.m8"
    fi

    ProgressBar "$state" "$total_states"
  done < <(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" \
    | xargs -n1 basename -s .fa)
  echo ""

  find "$foldseek_results" -size 0 -delete
}

# Aligns each orthogroup with MAFFT and runs hhblits against PfamA.
# Filters results to e-value <= 0.001 and removes empty output files.
# Arguments:
#   $1 - cluster or working directory base path
# Returns:
#   0 on success
run_hhblits() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local orthogroups_dir="${working_dir}/Orthogroups"
  local temp_dir="${working_dir}/temp"
  local hhblits_results="${working_dir}/hhblits"

  mkdir -p "$hhblits_results"

  local total_states state=0
  total_states=$(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" | wc -l)

  while read -r og_id; do
    state=$(( state + 1 ))

    mafft --maxiterate 1000 --genafpair --reorder \
      --thread "${threads}" \
      "${orthogroups_dir}/${og_id}.fa" \
      > "${temp_dir}/${og_id}_aligned.fa" 2>/dev/null

    hhblits \
      -i "${temp_dir}/${og_id}_aligned.fa" \
      -o "${hhblits_results}/${og_id}.hhr" \
      -blasttab "${temp_dir}/${og_id}.txt" \
      -d "${HHSUITE_PATH}/pfam" \
      -e 0.001 -n 6 -M 50 -z 2 -Z 10 \
      -realign_old_hits -cov 50 \
      -cpu "${threads}" &>/dev/null

    if [[ -f "${temp_dir}/${og_id}.txt" ]]; then
      awk 'BEGIN{FS=OFS="\t"} {if($11<=0.001) print}' \
        "${temp_dir}/${og_id}.txt" \
        > "${hhblits_results}/${og_id}.txt"
    fi

    ProgressBar "$state" "$total_states"
  done < <(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" \
    | xargs -n1 basename -s .fa)
  echo ""

  find "$hhblits_results" -size 0 -delete
}

# Runs InterProScan for each orthogroup using a fixed set of applications.
# Removes empty output files after completion.
# Arguments:
#   $1 - cluster or working directory base path
# Returns:
#   0 on success
run_interproscan() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local orthogroups_dir="${working_dir}/Orthogroups"
  local interpro_results="${working_dir}/InterProScan"

  mkdir -p "$interpro_results"

  local total_states state=0
  total_states=$(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" | wc -l)

  while read -r og_id; do
    state=$(( state + 1 ))

    "${INTERPRO_PATH}/interproscan.sh" \
      -i "${orthogroups_dir}/${og_id}.fa" \
      -f tsv \
      -appl CDD,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PRINTS,\
PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM \
      --goterms \
      -o "${interpro_results}/${og_id}.tsv" &>/dev/null

    ProgressBar "$state" "$total_states"
  done < <(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" \
    | xargs -n1 basename -s .fa)
  echo ""

  find "$interpro_results" -size 0 -delete
}

# Builds per-orthogroup and general annotation summary tables from InterProScan,
# Foldseek, and hhblits results.
# Arguments:
#   $1 - cluster or working directory base path
# Returns:
#   0 on success
create_summary_table() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local orthogroups_dir="${working_dir}/Orthogroups"
  local temp_dir="${working_dir}/temp"
  local interpro_results="${working_dir}/InterProScan"
  local hhblits_results="${working_dir}/hhblits"
  local foldseek_results="${working_dir}/Foldseek"
  local summary_results="${working_dir}/Summary"

  mkdir -p "$summary_results"
  echo "OrthogroupID;hhblits_domains;InterProScan;InterProScan_GO;Foldseek" \
    > "${working_dir}/General_summary.csv"

  while read -r og_id; do
    local protein_number
    protein_number=$(grep -c ">" "${orthogroups_dir}/${og_id}.fa")

    echo "ProteinID;InterProScan;InterProScan_GO;Foldseek" \
      > "${summary_results}/${og_id}.csv"
    grep ">" "${orthogroups_dir}/${og_id}.fa" \
      | awk -F '>' '{print $2}' | sort -k1,1 \
      > "${temp_dir}/${og_id}.csv"

    local interpro_general="" interprogo_general=""

    if [[ -f "${interpro_results}/${og_id}.tsv" ]]; then
      join -t ";" -a1 "${temp_dir}/${og_id}.csv" \
        <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $5"-\""$6"\""}' \
            "${interpro_results}/${og_id}.tsv" \
          | sed -e 's/[^[:print:]]$//' -e 's/;/,/g' \
          | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' \
          | sort -k1,1 | sed 's/\t/;/g') \
        | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' \
        > "${temp_dir}/${og_id}-2.csv"

      interpro_general=$(
        awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $5"-\""$6"\""}' \
          "${interpro_results}/${og_id}.tsv" \
          | sed -e 's/[^[:print:]]$//' -e 's/;/,/g' \
          | sort -k1,2 -u \
          | awk -F '\t' '{print $2}' \
          | sort | uniq -c | sed 's/^ *//g' \
          | awk -v n="$protein_number" '{if($1>=(n*0.5)) print}' \
          | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//'
      )

      awk 'BEGIN{FS=OFS="\t"}{if($14!="-") print $1 OFS $14}' \
        "${interpro_results}/${og_id}.tsv" \
        > "${temp_dir}/${og_id}-IPS.tsv"

      if [[ -s "${temp_dir}/${og_id}-IPS.tsv" ]]; then
        join -t ";" -a1 "${temp_dir}/${og_id}-2.csv" \
          <(sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' \
              "${temp_dir}/${og_id}-IPS.tsv" \
            | sort -k1,1 | sed 's/\t/;/g') \
          | awk 'BEGIN{FS=OFS=";"}{if($3==""){print $0";"}else{print $0}}' \
          > "${temp_dir}/${og_id}-3.csv"

        interprogo_general=$(
          sort -k1,2 -u "${temp_dir}/${og_id}-IPS.tsv" \
            | awk -F '\t' '{print $2}' \
            | sort | uniq -c | sed 's/^ *//g' \
            | awk -v n="$protein_number" '{if($1>=(n*0.5)) print}' \
            | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//'
        )
      else
        sed -e 's/$/;/g' "${temp_dir}/${og_id}-2.csv" \
          > "${temp_dir}/${og_id}-3.csv"
      fi
    else
      sed -e 's/$/;;/g' "${temp_dir}/${og_id}.csv" \
        > "${temp_dir}/${og_id}-3.csv"
    fi

    local foldseek_general=""

    if [[ -f "${foldseek_results}/${og_id}.m8" ]]; then
      local fs_col_awk
      if [[ "$foldseekdb" == "pdb" ]]; then
        fs_col_awk='{print $1 OFS $13}'
      else
        fs_col_awk='{print $1 OFS "\""$13"\""}'
      fi

      join -t ";" -a1 "${temp_dir}/${og_id}-3.csv" \
        <(awk "BEGIN{FS=OFS=\"\t\"} ${fs_col_awk}" \
            "${foldseek_results}/${og_id}.m8" \
          | sort -k1,2 -u \
          | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' \
          | sort -k1,1 | sed 's/\t/;/g') \
        | awk 'BEGIN{FS=OFS=";"}{if($8==""){print $0";"}else{print $0}}' \
        > "${temp_dir}/${og_id}-4.csv"

      foldseek_general=$(
        awk "BEGIN{FS=OFS=\"\t\"} ${fs_col_awk}" \
          "${foldseek_results}/${og_id}.m8" \
          | sort -k1,2 -u \
          | awk -F '\t' '{print $2}' \
          | sort | uniq -c | sed 's/^ *//g' \
          | awk -v n="$protein_number" '{if($1>=(n*0.5)) print}' \
          | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//'
      )
    else
      sed -e 's/$/;/g' "${temp_dir}/${og_id}-3.csv" \
        > "${temp_dir}/${og_id}-4.csv"
    fi

    cat "${temp_dir}/${og_id}-4.csv" >> "${summary_results}/${og_id}.csv"

    local hhblits_domains=""
    if [[ -f "${hhblits_results}/${og_id}.txt" ]]; then
      awk -F '\t' '{print $2}' "${hhblits_results}/${og_id}.txt" \
        | sort -u > "${temp_dir}/${og_id}-hhblits.txt"

      hhblits_domains=$(
        grep -f "${temp_dir}/${og_id}-hhblits.txt" \
          "${hhblits_results}/${og_id}.hhr" \
          | grep ">" \
          | awk -F '>' '{print $2}' \
          | sort -u \
          | sed -e 's/ ; /|/1' -e 's/ ; /-/' \
                -e 's/|/ \"/' -e 's/$/\"/' \
          | tr -s '\n' ',' | sed 's/,$//g'
      )
    fi

    echo "${og_id};${hhblits_domains};${interpro_general};${interprogo_general};${foldseek_general}" \
      >> "${working_dir}/General_summary.csv"
  done < <(find "$orthogroups_dir" -maxdepth 1 -name "*.fa" \
    | xargs -n1 basename -s .fa)
}

# Copies annotation outputs from the workspace to the final output directory.
# Arguments:
#   $1 - cluster or working directory base path
# Returns:
#   0 always
organize_information() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")-${mode}"
  local output_dir="${base_dir}/$(basename -s .sh "$0")-${mode}"

  mkdir -p "$output_dir"

  if [[ ! -d "$working_dir" ]]; then
    echo "Error: workspace directory '${working_dir}' not found." >&2
    return 1
  fi

  [[ -f "${working_dir}/General_summary.csv" ]] && \
    cp "${working_dir}/General_summary.csv" "${output_dir}/"

  [[ -d "${working_dir}/Summary" ]] && \
    cp -r "${working_dir}/Summary" "${output_dir}/Orthogroups_summary"
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

working_directory=""
mode="All"
foldseekdb="afdb_swissprot"
threads="8"
overwrite=false
clusters_file=""

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--workingDirectory)  shift; working_directory="$1" ;;
    -m|--mode)              shift; mode="$1" ;;
    -t|--threads)           shift; threads="$1" ;;
    -f|--foldseekdb)        shift; foldseekdb="$1" ;;
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
echo "  Mode:                            ${mode}"
echo "  Foldseek database:               ${foldseekdb}"
echo "  Threads:                         ${threads}"
echo "  Overwrite previous run:          ${overwrite}"
echo ""

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Checking arguments and input files..."

check_mode_parameter "${mode}" "$(basename -s .sh "$0")"

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

if [[ "$foldseekdb" != "pdb" && "$foldseekdb" != "afdb_swissprot" ]]; then
  echo "Error: foldseek database '${foldseekdb}' is not valid." >&2
  print_help; exit 1
fi
check_foldseek_databases "${FOLDSEEK_PATH}" "${foldseekdb}"

check_threads "${threads}"

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

# Helper to run annotation steps 1-4 on a given directory
run_annotation_steps() {
  local target_dir="$1"

  log_info "Step 1: Running InterProScan..."
  run_interproscan "${target_dir}"
  log_info "-> Step 1 finished. Proceeding."

  log_info "Step 2: Running Foldseek..."
  run_foldseek "${target_dir}"
  log_info "-> Step 2 finished. Proceeding."

  log_info "Step 3: Running hhblits..."
  run_hhblits "${target_dir}"
  log_info "-> Step 3 finished. Proceeding."

  log_info "Step 4: Creating summary table..."
  create_summary_table "${target_dir}"
  rm -r "${target_dir}/Workspace/$(basename -s .sh "$0")-${mode}/temp" \
    2>/dev/null || true
  log_info "Organizing output..."
  organize_information "${target_dir}"
  log_info "-> Step 4 finished. Proceeding."
}

if [[ "$mode" == "All" || "$mode" == "MoveAssociated" || "$mode" == "Core" ]]; then
  check_clusters "${working_directory}"

  while read -r cluster_id; do
    internal_dir="${working_directory}/Clusters/${cluster_id}"
    log_info "Analyzing cluster '${cluster_id}'."
    log_info "Checking workspace directory structure..."
    check_internal_directory_structure "${internal_dir}"

    if $directory_flag; then
      log_info "Directory structure valid. Proceeding."

      if $overwrite; then
        log_info "Removing previous run for cluster '${cluster_id}' if exists..."
        overwrite "${internal_dir}" "$(basename -s .sh "$0")" "${mode}"
      fi

      log_info "Organizing workspace..."
      organize_working_directory "${internal_dir}"

      run_annotation_steps "${internal_dir}"
      log_info "Cluster '${cluster_id}' complete.\n"
    else
      log_warning "Cluster '${cluster_id}' missing required directories for '${mode}' mode. Skipping.\n"
    fi
  done < <(
    awk '{print $1}' "$clusters_file" | sed $'s/[^[:print:]\t]//g'
  )

  log_info "All clusters have been analyzed."

elif [[ "$mode" == "Overrepresented" ]]; then
  check_overrepresentation "${working_directory}"

  if $overwrite; then
    log_info "Removing previous run if exists..."
    overwrite "${working_directory}" "$(basename -s .sh "$0")" "${mode}"
  fi

  log_info "Organizing workspace..."
  organize_working_directory_overrepresentation "${working_directory}"

  run_annotation_steps "${working_directory}"

  log_info "Annotation finished."
fi
