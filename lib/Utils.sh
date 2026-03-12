#!/usr/bin/env bash
# ==============================================================================
# NAME:        Utils.sh
# DESCRIPTION: Shared utility library for StarCrew pipeline scripts.
#              Provides logging, progress display, and workspace cleanup
#              functions used across all pipeline commands.
# USAGE:       source Utils.sh
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        12/Mar/2026
# VERSION:     1.0.0
# ==============================================================================

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays an inline progress bar with timestamp.
# Arguments:
#   $1 - current step (integer)
#   $2 - total steps (integer)
# Returns:
#   0 always
ProgressBar() {
  local _progress _done _left _fill _empty
  (( _progress = ($1 * 100 / $2 * 100) / 100 )) || true
  (( _done     = (_progress * 4) / 10 ))         || true
  (( _left     = 40 - _done ))                   || true
  _fill=$(printf "%${_done}s")
  _empty=$(printf "%${_left}s")
  printf "\rProgress $(date "+%Y-%m-%d %H:%M:%S") : [${_fill// /#}${_empty// /-}] ${_progress}%%"
}

# Removes previous run output and workspace directories for a given command.
# For CaptainIdentification, also restores captainless elements to Data/.
# For OrthogroupsAnnotation, also removes the mode-suffixed directories.
# For SyntenyClustering, also removes the Clusters/ directory.
# Arguments:
#   $1 - base working directory path
#   $2 - command name
#   $3 - mode suffix (optional, used by OrthogroupsAnnotation)
# Returns:
#   0 always
overwrite() {
  local working_dir="$1"
  local command="$2"
  local mode="${3:-}"
  local workspace="${working_dir}/Workspace"

  # Remove top-level output directory
  if [[ -d "${working_dir}/${command}" ]]; then
    log_info "Removing '${command}' output directory..."
    rm -r "${working_dir}/${command}"
  fi

  # Remove workspace directory
  if [[ -d "${workspace}/${command}" ]]; then
    log_info "Removing '${command}' workspace directory..."
    rm -r "${workspace}/${command}"
  fi

  case "$command" in

    SyntenyClustering)
      if [[ -d "${working_dir}/Clusters" ]]; then
        log_info "Removing 'Clusters' directory..."
        rm -r "${working_dir}/Clusters"
      fi
      ;;

    CaptainIdentification)
      local captainless_dir="${working_dir}/Data/Captainless_elements"
      if [[ -d "$captainless_dir" ]]; then
        log_info "Restoring captainless elements to Data/..."

        local data_dir gff_dir protein_dir nucleotide_dir cds_dir
        data_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)
        gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
        protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
        nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
        cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

        while read -r element; do
          mv "${captainless_dir}/${element}.gff"            "${gff_dir}/"
          mv "${captainless_dir}/${element}_protein.fa"     "${protein_dir}/${element}.fa"
          mv "${captainless_dir}/${element}_nucleotide.fa"  "${nucleotide_dir}/${element}.fa"
          mv "${captainless_dir}/${element}_CDS.fa"         "${cds_dir}/${element}.fa"
        done < <(find "$captainless_dir" -maxdepth 1 -name "*.gff" \
          | xargs -n1 basename -s .gff)

        rm -r "$captainless_dir"
      fi
      ;;

    OrthogroupsAnnotation)
      if [[ -n "$mode" ]]; then
        if [[ -d "${working_dir}/${command}-${mode}" ]]; then
          log_info "Removing '${command}-${mode}' output directory..."
          rm -r "${working_dir}/${command}-${mode}"
        fi
        if [[ -d "${workspace}/${command}-${mode}" ]]; then
          log_info "Removing '${command}-${mode}' workspace directory..."
          rm -r "${workspace}/${command}-${mode}"
        fi
      fi
      ;;

  esac
}

# Prints a timestamped informational message to stdout.
# Arguments:
#   $1 - message to display
# Returns:
#   0 always
log_info() {
  echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] $1"
}

# Prints a timestamped warning message in red to stdout.
# Arguments:
#   $1 - message to display
# Returns:
#   0 always
log_warning() {
  echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] \033[01;31mWARNING\033[m: $1"
}

# Prints a timestamped indented step-level message to stdout.
# Arguments:
#   $1 - message to display
# Returns:
#   0 always
log_step() {
  echo -e "  [$(date "+%Y-%m-%d %H:%M:%S")] $1"
}