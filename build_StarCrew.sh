#!/usr/bin/env bash

# Exit with any non-zero status from a command.
echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Verifying library files...\n"

# Setting amain directory.
Main_dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
databases_dir="${Main_dir}/databases/"

# Check library directory
LIB_DIR="${Main_dir}/lib/"

# Check that the libraries are there
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

# Check conda environment
Conda_environment_presence=$(conda info --envs | grep -w -c "StarCrew")

if [[ $Conda_environment_presence -eq 0 ]]; then
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Creating conda environment with required libraries...\n"

    conda env create -f ${Main_dir}/StarCrew_environment.yml
    source activate StarCrew

    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Checking environment installation...\n"
    check_required_software "All"
    
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Checking main directory...\n"
    check_main_directory "${Main_dir}"
    chmod +x ${Main_dir}/bin/StarCrew
else
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] There's an existing 'StarCrew' environment, checking if everything is installed...\n"
    source activate StarCrew
    check_required_software "All"

    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Checking main directory...\n"
    check_main_directory "${Main_dir}"
    chmod +x ${Main_dir}/bin/StarCrew
fi

echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Intalling gggenomes and Checking the correct installation of R packages...\n"

Rscript -e '
  # Install packages
  install.packages("gggenomes", repos = "https://cloud.r-project.org/")
  
  # Checking the installation
  suppressPackageStartupMessages(library(optparse))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ape))
  suppressPackageStartupMessages(library(ggtree))
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(viridis))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(gggenomes))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(dendextend))
  suppressPackageStartupMessages(library(NbClust))
  suppressPackageStartupMessages(library(syntenet))
  suppressPackageStartupMessages(library(ggnewscale))
  suppressPackageStartupMessages(library(bc3net))
  suppressPackageStartupMessages(library(mrfDepth))

  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
   stop("Package \"ggplot2\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
   stop("Package \"ape\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("ggtree", quietly = TRUE)) {
   stop("Package \"ggtree\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
   stop("Package \"reshape2\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
   stop("Package \"viridis\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
   stop("Package \"dplyr\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("gggenomes", quietly = TRUE)) {
   stop("Package \"gggenomes\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
   stop("Package \"scales\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("dendextend", quietly = TRUE)) {
   stop("Package \"dendextend\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("NbClust", quietly = TRUE)) {
   stop("Package \"NbClust\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("syntenet", quietly = TRUE)) {
   stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
   stop("Package \"ggnewscale\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("bc3net", quietly = TRUE)) {
   stop("Package \"bc3net\" not installed. Please install it to run this script.", call. = FALSE)
  }
  if (!requireNamespace("mrfDepth", quietly = TRUE)) {
   stop("Package \"mrfDepth\" not installed. Please install it to run this script.", call. = FALSE)
  }

  cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Successful instalation of R packages"," \n", sep=""))
'

echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Modifying config.json file of Orthofinder...\n"

python json_updater.py "$(echo $(which orthofinder | awk 'BEGIN{FS=OFS="/"}NF--')"/src/orthofinder/run/config.json")" --cmd-update

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

cd ${databases_dir}

echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Creating directories for protein annotation databases..\n"

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

cd ${databases_dir}

echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] Checking full instalation..."

check_installation "${Main_dir}"

echo -e "\n[$(date "+%Y-%m-%d %H:%M:%S")] installation Successful."