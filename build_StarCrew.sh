#!/usr/bin/env bash
# ==============================================================================
# NAME:        build_StarCrew.sh
# DESCRIPTION: Sets up the StarCrew environment, including conda environment
#              creation, R package installation, and database downloads.
# USAGE:       ./build_StarCrew.sh [options]
#              ./build_StarCrew.sh -help
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        11/Mar/2026
# VERSION:     1.0.0
# ==============================================================================

set -o pipefail

# ==============================================================================
# CONSTANTS
# ==============================================================================

readonly SCRIPT_DIR="$( dirname -- "$( readlink -f -- "$0" )" )"
readonly MAIN_DIR="$(realpath "${SCRIPT_DIR}")"
readonly LIB_DIR="${MAIN_DIR}/lib"
readonly DATABASES_DIR="${MAIN_DIR}/databases"

readonly INTERPROSCAN_VERSION="5.77-108.0"
readonly INTERPROSCAN_ARCHIVE="interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz"
readonly INTERPROSCAN_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPROSCAN_VERSION}/${INTERPROSCAN_ARCHIVE}"

# ==============================================================================
# LOGGING
# ==============================================================================

# Prints a timestamped informational message to stdout.
# Arguments:
#   $1 - message to display
# Returns:
#   0 always
log_info() {
  echo "[$(date "+%Y-%m-%d %H:%M:%S")] $1"
}

# Prints a timestamped warning message to stdout.
# Arguments:
#   $1 - message to display
# Returns:
#   0 always
log_warning() {
  echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] \033[01;31mWARNING\033[m: $1"
}

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays the help message with all available options.
# Arguments:
#   None
# Returns:
#   0 always
print_help() {
  echo "Script to set up the StarCrew environment."
  echo
  echo "Usage: bash $(basename "$0") [options]"
  echo
  echo "Options:"
  echo "  --skipDatabases   Skip the installation of databases for OrthogroupsAnnotation."
  echo -e "                    \033[01;31mWARNING\033[m:: Leaves the OrthogroupsAnnotation command unusable."
  echo "                    Recommended when disk space is limited (Default: off)."
  echo "  --onlyDatabases   Only install databases; skip conda environment setup."
  echo -e "                    \033[01;31mWARNING\033[m:: Expects a previous conda environment setup."
  echo "                    Not compatible with '--skipDatabases' (Default: off)."
  echo "  -help             Display this help message."
}

# Verifies that required library files exist and sources Check.sh.
# Arguments:
#   None
# Returns:
#   0 if all libraries are present, exits with 1 otherwise
check_libraries() {
  log_info "Verifying library files..."

  if [[ ! -d "$LIB_DIR" ]]; then
    echo "Error: folder '${LIB_DIR}' does not exist." >&2
    exit 1
  fi

  if [[ ! -f "${LIB_DIR}/Check.sh" ]]; then
    echo "Error: file '${LIB_DIR}/Check.sh' does not exist." >&2
    exit 1
  fi
  source "${LIB_DIR}/Check.sh"

  if [[ ! -f "${LIB_DIR}/Utils.sh" ]]; then
    echo "Error: file '${LIB_DIR}/Utils.sh' does not exist." >&2
    exit 1
  fi
}

# Installs and verifies required R packages.
# Arguments:
#   None
# Returns:
#   0 on success, exits with 1 if any package is missing after installation
install_r_packages() {
  log_info "Verifying required R packages..."

  Rscript -e '
    # Install packages
    install.packages("gggenomes", repos = "https://cloud.r-project.org/")
      
    required_pkgs <- c(
      "optparse", "ggplot2", "ape", "ggtree", "reshape2",
      "viridis", "dplyr", "gggenomes", "scales", "dendextend",
      "NbClust", "syntenet", "ggnewscale", "bc3net", "mrfDepth"
    )

    suppressPackageStartupMessages({
      for (pkg in required_pkgs) { library(pkg, character.only = TRUE) }
    })

    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste0("Package \"", pkg, "\" not installed."), call. = FALSE)
      }
    }

    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] R packages installed successfully.\n"))
  '
}

# Updates the OrthoFinder config.json using the json_updater.py helper script.
# Arguments:
#   None
# Returns:
#   0 on success
configure_orthofinder() {
  log_info "Modifying OrthoFinder config.json file..."

  local orthofinder_config
  orthofinder_config="$(dirname "$(which orthofinder)")/src/orthofinder/run/config.json"
  python ${MAIN_DIR}/json_updater.py "${orthofinder_config}" --cmd-update
}

# Creates or validates the StarCrew conda environment and installs R packages.
# Arguments:
#   None
# Returns:
#   0 on success, exits with 1 on error
setup_conda_environment() {
  local conda_env_count
  conda_env_count=$(conda info --envs | grep -w -c "StarCrew")

  if [[ "$conda_env_count" -eq 0 ]]; then
    log_info "Creating conda environment with required libraries..."
    conda env create -f "${MAIN_DIR}/StarCrew_environment.yml"
  else
    log_info "Existing 'StarCrew' environment found. Checking installation..."
  fi

  source activate StarCrew
  check_required_software "All"

  log_info "Checking main directory..."
  check_main_directory "${MAIN_DIR}"
  chmod +x "${MAIN_DIR}/bin/StarCrew"
  ln -sf "${MAIN_DIR}/bin/StarCrew" "$(which seqkit | sed -e 's/seqkit//')"

  install_r_packages
  configure_orthofinder

  log_info "Setting up autocompletion..."
  cp "${MAIN_DIR}/starcrew_completion.sh" "$(which StarCrew | sed -e 's/bin\/StarCrew//')etc/conda/activate.d/"
  cp "${MAIN_DIR}/starcrew_deactivate_completion.sh" "$(which StarCrew | sed -e 's/bin\/StarCrew//')etc/conda/deactivate.d/"
}

# Downloads, verifies, and installs InterProScan and annotation databases.
# Arguments:
#   None
# Returns:
#   0 on success, exits with 1 on checksum failure
install_databases() {
  log_info "Installing and configuring InterProScan..."

  mkdir -p "${MAIN_DIR}/interproscan"
  wget "${INTERPROSCAN_URL}"
  wget "${INTERPROSCAN_URL}.md5"

  md5sum -c "${INTERPROSCAN_ARCHIVE}.md5"
  tar -pxvzf "${INTERPROSCAN_ARCHIVE}" -C "${MAIN_DIR}/interproscan" --strip-components=1
  rm "${INTERPROSCAN_ARCHIVE}"*

  cd "${MAIN_DIR}/interproscan"
  python3 setup.py -f interproscan.properties
  cd "${MAIN_DIR}"

  # Foldseek databases
  log_info "Creating directories for protein annotation databases..."
  mkdir -p "${DATABASES_DIR}/Foldseek" "${DATABASES_DIR}/hhsuite"

  cd "${DATABASES_DIR}/Foldseek"
  log_info "Downloading Foldseek databases..."
  foldseek databases ProstT5 weights tmp
  foldseek databases PDB pdb tmp
  foldseek databases Alphafold/Swiss-Prot afdb_swissprot tmp

  log_info "Generating index file for PDB database..."
  wget "https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
  awk '{FS=OFS="\t"}NR>2{if($2!="DNA" && $2!="DNA-RNA HYBRID"){print}}' entries.idx \
    | sort -k1,1 > entries_update.idx
  rm entries.idx

  log_info "Generating index file for AlphaFold/Swiss-Prot database..."
  wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
  gunzip uniprot_sprot.fasta.gz
  grep ">" uniprot_sprot.fasta \
    | awk -F '|' '{print $2 OFS $3}' \
    | sed -e 's/[A-Z]*=.*//g' -e 's/ $//' \
    | awk '{$2="\t"} 1' \
    | sed 's/ \t /\t/' \
    | sort -k1,1 > Accession_swissprot.txt
  rm uniprot_sprot.fasta

  # HH-suite databases
  cd "${DATABASES_DIR}/hhsuite"
  log_info "Downloading HH-suite pfam database..."
  wget "https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz"
  tar -zxvf pfamA_35.0.tar.gz
  rm pfamA_35.0.tar.gz

  cd "${MAIN_DIR}"
}

# Parses arguments and orchestrates the setup workflow.
# Arguments:
#   $@ - original script arguments
# Returns:
#   0 on success, 1 on invalid argument combination
main() {
  local only_databases=false
  local skip_databases=false

  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --skipDatabases)
        skip_databases=true
        ;;
      --onlyDatabases)
        only_databases=true
        ;;
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

  # Validate incompatible flags
  if $only_databases && $skip_databases; then
    echo "Error: '--onlyDatabases' and '--skipDatabases' are mutually exclusive." >&2
    exit 1
  fi

  check_libraries

  # Conda environment setup
  if $only_databases; then
    log_info "Skipping conda environment setup..."
    local conda_env_count
    conda_env_count=$(conda info --envs | grep -w -c "StarCrew")
    if [[ "$conda_env_count" -eq 0 ]]; then
      echo "Error: no 'StarCrew' environment found. Run without '--onlyDatabases' first." >&2
      exit 1
    fi
  else
    setup_conda_environment
  fi

  # Database installation
  if $skip_databases; then
    log_info "Skipping database installation."
    log_warning "The 'OrthogroupsAnnotation' command will not be available."
  else
    install_databases
    log_info "Checking full installation..."
    check_installation "${MAIN_DIR}" "${skip_databases}"
  fi

  log_info "Installation successful."
}

main "$@"
