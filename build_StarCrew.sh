#!/usr/bin/env bash

# ==============================================================================
# STARCREW ENVIRONMENT SETUP
# ==============================================================================

# ==============================================================================
# DIRECTORY CHECK
# ==============================================================================

echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Verifying library files...\n"

# Setting main directory and paths
Main_dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
Main_dir="$(realpath $Main_dir)"
databases_dir="${Main_dir}/databases/"
LIB_DIR="${Main_dir}/lib/"

# Library Verification
if [[ ! -d "$LIB_DIR" ]]; then
    echo "Error: folder '$LIB_DIR' does not exist."
    exit 1
else
    LIB_DIR=$(realpath $LIB_DIR)  
    if [[ ! -f "${LIB_DIR}/Check.sh" ]]; then
        echo "Error: file '${LIB_DIR}/Check.sh' does not exist."
        exit 1
    else
        source "${LIB_DIR}/Check.sh"
    fi

    if [[ ! -f "${LIB_DIR}/Utils.sh" ]]; then
        echo "Error: file '${LIB_DIR}/Utils.sh' does not exist."
        exit 1
    fi
fi

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

function print_help() {
    echo -e "Script to setup StarCrew environment."
    echo ""
    echo "Syntax: bash $(basename "$0" ) [options...]"
    echo ""
    echo "Optional args:"
    echo "--skipDatabases: Flag to skip the installation of the databases for the OrthogroupsAnnotation command. \033[01;31mWARNING\033[m: This setup left the OrthogroupsAnnotation command unusable,
          recommended when there is not enough disk space available for the databases (Default: off)."
    echo "--onlyDatabases: Flag to only install the databases and skip the conda environment setup. \033[01;31mWARNING\033[m: This flag expect a previous setup of the conda environment.
                              Not compatible with '--skipDatabases' flag (Default: off). "
    echo "-help: Display this help message."
}

# ==============================================================================
# VARIABLES AND ARGUMENT PARSING
# ==============================================================================

# Initialize variables
only_databases=false
skip_databases=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skipDatabases)
            skip_databases=true
            ;;
        --onlyDatabases)
            only_databases=true
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

if $only_databases && $skip_databases; then
    echo "Error: '--onlyDatabases' and '--skip-skipDatabases' flags are both 'on' and those are incompatible."
    exit 1
fi

# Conda Environment Setup
Conda_environment_presence=$(conda info --envs | grep -w -c "StarCrew")

if $only_databases; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] skipping conda environment setup..."
    if [[ $Conda_environment_presence -eq 0 ]]; then
        echo "Error: There's no 'Starcrew' environment. Please run this script without the '--only_databases' flag."
        print_help
        exit 1
    fi
else
    if [[ $Conda_environment_presence -eq 0 ]]; then
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Creating conda environment with required libraries...\n"

        conda env create -f ${Main_dir}/StarCrew_environment.yml
        source activate StarCrew

        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Checking environment installation...\n"
        check_required_software "All"
    
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Checking main directory...\n"
        check_main_directory "${Main_dir}"
        chmod +x ${Main_dir}/bin/StarCrew
        ln -sf ${Main_dir}/bin/StarCrew $(which seqkit | sed -e 's/seqkit//')
    else
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] There's an existing 'StarCrew' environment, checking if everything is installed...\n"
        source activate StarCrew
        check_required_software "All"

        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Checking main directory...\n"
        check_main_directory "${Main_dir}"
        chmod +x ${Main_dir}/bin/StarCrew
        ln -sf ${Main_dir}/bin/StarCrew $(which seqkit | sed -e 's/seqkit//')
    fi

    # R Package Installation & Verification
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Intalling gggenomes and Checking the correct installation of R packages...\n"

    Rscript -e '
      # Install packages
      install.packages("gggenomes", repos = "https://cloud.r-project.org/")
  
      # Checking the installation
      required_pkgs <- c(
        "optparse", "ggplot2", "ape", "ggtree", "reshape2", 
        "viridis", "dplyr", "gggenomes", "scales", "dendextend", 
        "NbClust", "syntenet", "ggnewscale", "bc3net", "mrfDepth"
      )

      suppressPackageStartupMessages({
        library(optparse)
        library(ggplot2)
        library(ape)
        library(ggtree)
        library(reshape2)
        library(viridis)
        library(dplyr)
        library(gggenomes)
        library(scales)
        library(dendextend)
        library(NbClust)
        library(syntenet)
        library(ggnewscale)
        library(bc3net)
        library(mrfDepth)
      })

      for (pkg in required_pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(paste0("Package \"", pkg, "\" not installed."), call. = FALSE)
        }
      }

      cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Successful instalation of R packages"," \n", sep=""))
    '

    # Orthofinder configuration
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Modifying config.json file of Orthofinder...\n"

    python json_updater.py "$(echo $(which orthofinder | awk 'BEGIN{FS=OFS="/"}NF--')"/src/orthofinder/run/config.json")" --cmd-update
fi

if $skip_databases; then
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Skipping database installation"
    echo -e "\033[01;31mWARNING\033[m: you're not going to be able to use the 'OrthogroupsAnnotation' command"
else
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Installing and configuring Interproscan...\n"
    mkdir -p ${Main_dir}/interproscan
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz.md5

    md5sum -c interproscan-5.76-107.0-64-bit.tar.gz.md5
    tar -pxvzf interproscan-5.76-107.0-64-bit.tar.gz -C interproscan --strip-components=1

    cd ${Main_dir}/interproscan/
    python3 setup.py -f interproscan.properties
    cd ${Main_dir}
    rm interproscan-5.76-107.0-64-bit.tar.gz*

    # Database Download (Foldseek & HH-suite)
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Creating directories for protein annotation databases..\n"
    cd ${databases_dir}
    mkdir -p Foldseek hhsuite

    cd ${databases_dir}/Foldseek/
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Downloading Foldseek databases...\n"
    foldseek databases ProstT5 weights tmp
    foldseek databases PDB pdb tmp
    foldseek databases Alphafold/Swiss-Prot afdb_swissprot tmp

    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Generating index file for PDB database...\n"
    wget "https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
    awk '{FS=OFS="\t"}NR>2{if($2!="DNA" && $2!="DNA-RNA HYBRID"){print}}' entries.idx | sort -k1,1 > entries_update.idx
    rm entries.idx

    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Generating index file for alphaphold database...\n"
    wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip uniprot_sprot.fasta.gz
    grep ">" uniprot_sprot.fasta | awk -F '|' '{print $2 OFS $3}' | sed -e 's/[A-Z]*=.*//g' -e 's/ $//' | awk '{$2="\t"} 1' | sed 's/ \t /\t/' | sort -k1,1 > Accession_swissprot.txt
    rm uniprot_sprot.fasta

    cd ${databases_dir}/hhsuite/
    echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Downloading hhsuite pfam database...\n"
    wget "https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz"
    tar -zxvf pfamA_35.0.tar.gz
    rm pfamA_35.0.tar.gz

    cd ${Main_dir}
fi

echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Checking full instalation..."
check_installation "${Main_dir}" "${skip_databases}"
echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] installation Successful."
