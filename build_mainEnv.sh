#!/bin/bash

# Setting the database directory
Main_dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
aux_dir="${Main_dir}/aux/"
databases_dir="${Main_dir}/databases/"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Creating conda environment with required libraries..."

conda env create -f ${Main_dir}/SAT_environment.yml
source activate SAT

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Installing R packages..."
# There may be an issue for the installation
Rscript -e '
  # Install packages
  install.packages(c("ape","reshape2","viridis","dplyr","gggenomes","scales","dendextend","NbClust","svglite","optparse","future"), repos = "https://cloud.r-project.org/")
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

  cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Successful instalation of R packages"," \n", sep=""))
'

echo "Donwloading MACSE..."
wget -O ${aux_dir}/macse.jar "https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Modifying config.json file of Orthofinder..."

python python ${Main_dir}/json_updater.py "$(echo $(which orthofinder | awk 'BEGIN{FS=OFS="/"}NF--')"/src/orthofinder/run/config.json")" --cmd-update

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Downloading Mycomobilome database..."

cd ${databases_dir}
curl -O "https://zenodo.org/records/17037469/files/MycoMobilome_v1.0.tar.gz"
tar -zxvf MycoMobilome_v1.0.tar.gz
mv MycoMobilome_v1.0/ MycoMobilome_db
rm MycoMobilome_v1.0.tar.gz
