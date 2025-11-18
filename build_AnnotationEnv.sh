z#!/bin/bash

# Setting the database directory
Main_dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
databases_dir="${Main_dir}/databases/"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Creating conda environment with required libraries..."

conda env create -f ${Main_dir}/Annotation_environment.yml

conda init

source activate Annotationtest2

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Installing DeepFRI and related database..."

git clone "https://github.com/flatironinstitute/DeepFRI.git"
cd DeepFRI/
pip install .

wget -O trained_models.tar.gz "https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz"
tar -xvzf trained_models.tar.gz -C ./
rm trained_models.tar.gz

python ${Main_dir}/json_updater.py ${Main_dir}/DeepFRI/trained_models/model_config.json --path-replace ${Main_dir}/DeepFRI/trained_models

cd ${databases_dir}

echo "Creating directories.."

mkdir -p Foldseek hhsuite

cd Foldseek/
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Downloading Foldseek databases..."

foldseek databases ProstT5 weights tmp
foldseek databases PDB pdb tmp
foldseek databases Alphafold/Swiss-Prot afdb_swissprot tmp

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Generating index file for PDB database..."

wget "https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
awk '{FS=OFS="\t"}NR>2{if($2!="DNA" && $2!="DNA-RNA HYBRID"){print}}' entries.idx | sort -k1,1 > entries_update.idx
rm entries.idx

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Generating index file for alphaphold database..."

wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
gunzip uniprot_sprot.fasta.gz
grep ">" uniprot_sprot.fasta | awk -F '|' '{print $2 OFS $3}' | sed -e 's/[A-Z]*=.*//g' -e 's/ $//' | awk '{$2="\t"} 1' | sed 's/ \t /\t/' | sort -k1,1 > Accession_swissprot.txt
rm uniprot_sprot.fasta

cd ../hhsuite/
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Downloading hhsuite pfam database..."

wget "https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz"
tar -zxvf pfamA_35.0.tar.gz
rm pfamA_35.0.tar.gz

cd ../
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Downloading Mycomobilome database..."

curl -O "https://zenodo.org/records/17037469/files/MycoMobilome_v1.0.tar.gz"
tar -zxvf MycoMobilome_v1.0.tar.gz
mv MycoMobilome_v1.0/ MycoMobilome_db
rm MycoMobilome_v1.0.tar.gz

