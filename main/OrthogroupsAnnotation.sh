#!/bin/bash

# Function to print help message

function print_help() {
   echo -e "Script to run a functional annotation for orthogroups.
   This script perform xxxxx steps:
   1. Identify orthogroups through OrthoFinder software.
   2. Perform all-vs-all Blastn.
   3. Perform a hierarchical clustering of the elements based on Orthogroup gene count including singletons.
   4. Determine the full conection of the cluster and create a cargo orthogroups heatmap and synteny image for the cluster.
   5. Identify possible individual nesting events inside the cluster.
   6. Identify core genes in the cluster in two ways:
     6.1. General core: orthogroups that are present in at least 80% of the elements in the cluster.
     6.2. Specific core:
       6.2.1. Divide the Cluster in subclusters of a height above 0.8 in the hierarchical clustering.
       6.2.2. If subslusters are generated identify core genes in each one that have at least 5 elements using the same logic of general core.
   7. If subclusters are present it try to identify putative cargo movement events including the specific orthogroups involve.
   8. Determine if there are discordances at 'Clade' lavel between CArgo hierarchical clustering and Captain phylogenetic tree.
   "
   echo
   echo "Syntax: SAT ClusterCharacterization [ -help ] -w <directory_path> -c <file_path> [ -t <integer> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-m, --mode: Define the orthogroups to be analyzed (Default = Core) [Available mode: MoveAssociated, Core, All]."
   echo "-c, --clusters: file with a list of clusters to be analyzed, each line correspond to a single cluster ID (required)."
   echo "-f, --foldseekdb: Name of the Foldseek database to use (Default = pdb) [Available: pdb, afdb_swissprot]."
   echo "-t, --threads: Number of threads for all analysis (Default: 8)"
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
mode="All"
clusters_file=""
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
database_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../databases/"
DeepFRI_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../DeepFRI/"
foldseekdb="pdb"
threads="8"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workingDirectory)
            shift
            Working_directory="$1"
            ;;
        -m|--mode)
            shift
            mode="$1"
            ;;
        -c|--clusters)
            shift
            clusters_file="$1"
            ;;
        -t|--threads)
            shift
            threads="$1"
            ;;
        -f|--foldseekdb)
            shift
            foldseekdb="$1"
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

# Print help if requested
if $help_flag; then
    print_help
    exit 0
fi

# Check for mandatory argument and define the path as absolute
if [[ -z "$Working_directory" || -z "$clusters_file" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if working directory exists
if [[ ! -d "$Working_directory" ]]; then
    echo "Error: Directory '$Working_directory' does not exist."
    exit 1
else
    Working_directory=$(realpath $Working_directory)
fi

# Check if cluster file exist
if [[ ! -f "$clusters_file" ]]; then
    echo "Error: File '$clusters_file' does not exist."
    exit 1
else
    clusters_file=$(realpath $clusters_file)
fi

# Check if mode parameter is correct
if [[ "$mode" != "All" && "$mode" != "MoveAssociated" && "$mode" != "Core" ]]; then
    echo "Error: provided mode '$mode' is not accepted."
    print_help
    exit 1
fi

# Check thread parameter
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
    exit 1
fi

# Check if database directory exists
if [[ ! -d "$database_path" ]]; then
    echo "Error: Directory '$database_path' does not exist."
    exit 1
else
    database_path=$(realpath $database_path)
    if [[ ! -d "${database_path}/Foldseek/" ]]; then
        echo "Error: Directory 'Foldseek' does not exist in '$database_path'."
        exit 1
    else
        foldseek_path="${database_path}/Foldseek/"
        if [[ "$foldseekdb" != "pdb" && "$foldseekdb" != "afdb_swissprot" ]]; then
            echo "Error: '$foldseekdb' is not accepted as a foldseek database."
            print_help
            exit 1
        fi

        if [[ ! -d "$foldseek_path/weights/" ]]; then
            echo "Error: Directory 'weights' does not exist in '${foldseek_path}'."
            exit 1
        fi

        if [[ ! -f "$foldseek_path/${foldseekdb}" ]]; then
            echo "Error: 'pdb' database does not exist in '${foldseek_path}'."
            exit 1
        fi

        if [[ "$foldseekdb" == "pdb" ]]; then
            if [[ ! -f "$foldseek_path/entries_update.idx" ]]; then
                echo "Error: 'entries_update.idx' file does not exist in '${foldseek_path}'."
                exit 1
            fi
        fi

        if [[ "$foldseekdb" == "afdb_swissprot" ]]; then
            if [[ ! -f "$foldseek_path/Accession_swissprot.txt" ]]; then
                echo "Error: 'Accession_swissprot.txt' file does not exist in '${foldseek_path}'."
                exit 1
            fi
        fi
    fi

    if [[ ! -d "${database_path}/hhsuite/" ]]; then
        echo "Error: Directory 'hhsuite' does not exist in '$database_path'."
        exit 1
    else
        hhsuite_path="${database_path}/hhsuite/"
        if [[ ! -f "$hhsuite_path/pfam.md5sum" ]]; then
            echo "Error: 'pfam' database does not exist in '${hhsuite_path}'."
            exit 1
        fi
    fi
fi

# Check for required software
if [[ -z "$(which mafft)" ]]; then
    echo "Error: Missing mafft function."
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

check_clusters() {
    # Modify based on the input file to compare with the current set of Clusters ID 
    local base_dir="$1"
    local file="$2"

    local cluster_original="$base_dir/cluster_file.txt"
    ls -d ${base_dir}/Clusters/*/ | awk -F '/' '{print $(NF - 1)}' > ${cluster_original}

    local diff=$(comm -13 <(sort ${cluster_original}) <(sort ${clusters_file}))

    if [[ $diff != "" ]]; then
        rm ${cluster_original} 
        echo "ERROR: there are additional lines no compatible to current ClusterID in ${clusters_file}."
        echo "Check for this lines: ${diff}"
        exit 1
    else
        rm ${cluster_original}
    fi
}

check_directory_structure() {
    local base_dir="$1"
    
    # Locate required subdirectories and file
    local workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    fi

    local ClusterCharacterization_dir=$(find "$workspace_dir" -maxdepth 1 -type d -name "ClusterCharacterization" 2>/dev/null)

    if [[ -z "$ClusterCharacterization_dir" ]]; then
        echo "Error: ClusterCharacterization directory not found in '$workspace_dir'." >&2
        echo "The ClusterCharacterization command should be run before this analysis" >&2
        exit 1
    fi
    
    local moveOrthologs_dir=$(find "$ClusterCharacterization_dir" -maxdepth 1 -type d -name "*_moveOrthologs" 2>/dev/null)
    local CoreGenes_dir=$(find "$ClusterCharacterization_dir" -maxdepth 1 -type d -name "Core_genes" 2>/dev/null)
    local Orthogroups_dir=$(find "$ClusterCharacterization_dir" -maxdepth 3 -type d -name "Orthogroup_Sequences" 2>/dev/null)

    directory_flag=true

    if [[ $mode == "All" ]]; then
        if [[ -z "$Orthogroups_dir" ]]; then
            echo "Error: GFF subdirectory not found in '$ClusterCharacterization_dir'." >&2
            directory_flag=false
        fi
    elif [[ $mode == "MoveAssociated" ]]; then
        if [[ -z "$moveOrthologs_dir" ]]; then
            echo "Orthogroups subdirectory for moving genes not found in '$ClusterCharacterization_dir'." >&2
            directory_flag=false
        fi
    elif [[ $mode == "Core" ]]; then
        if [[ -z "$CoreGenes_dir" ]]; then
            echo "Error: GFF subdirectory not found in '$ClusterCharacterization_dir'." >&2
            directory_flag=false
        fi
    fi
}

organize_working_directory() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/OrthogroupsAnnotation/"
    local ClusterCharacterization_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"

    mkdir -p ${working_dir}
    mkdir -p ${Orthogroups_dir}
    mkdir -p ${temp_dir}

    if [[ $mode == "All" ]]; then
        local data_dir=$(find "$ClusterCharacterization_dir" -maxdepth 3 -type d -name "Orthogroup_Sequences" 2>/dev/null)
        cp ${data_dir}/* ${Orthogroups_dir}/
    elif [[ $mode == "MoveAssociated" ]]; then
        local moveOrthologs_dir=$(find "$ClusterCharacterization_dir" -maxdepth 1 -type d -name "*_moveOrthologs" 2>/dev/null)
        echo $moveOrthologs_dir | awk '{OFS=RS;$1=$1}1' | while read line; 
        do 
            cp ${line}/* ${Orthogroups_dir}/
        done
    elif [[ $mode == "Core" ]]; then
        local CoreGenes_dir=$(find "$ClusterCharacterization_dir" -maxdepth 1 -type d -name "Core_genes" 2>/dev/null)
        cp ${CoreGenes_dir}/* ${Orthogroups_dir}/
    fi
}

run_deepfri() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/OrthogroupsAnnotation/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local DeepFRI_results="${working_dir}/DeepFRI/"

    mkdir -p "${DeepFRI_results}"

    ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
    do
        python ${DeepFRI_path}/predict.py --fasta_fn ${Orthogroups_dir}/${OrthogroupID}.fa --model_config ${DeepFRI_path}/trained_models/model_config.json -ont mf -o ${temp_dir}/${OrthogroupID} &>/dev/null
        grep -v "^#" ${temp_dir}/${OrthogroupID}_MF_predictions.csv | awk '{FS=OFS=","}{if($3 >= 0.5){print $0}}' | sed 's/,/\t/;s/,/\t/;s/,/\t/' > ${DeepFRI_results}/${OrthogroupID}_MF_predictions.txt
        python ${DeepFRI_path}/predict.py --fasta_fn ${Orthogroups_dir}/${OrthogroupID}.fa --model_config ${DeepFRI_path}/trained_models/model_config.json -ont bp -o ${temp_dir}/${OrthogroupID} &>/dev/null
        grep -v "^#" ${temp_dir}/${OrthogroupID}_BP_predictions.csv | awk '{FS=OFS=","}{if($3 >= 0.5){print $0}}' | sed 's/,/\t/;s/,/\t/;s/,/\t/' > ${DeepFRI_results}/${OrthogroupID}_BP_predictions.txt
        python ${DeepFRI_path}/predict.py --fasta_fn ${Orthogroups_dir}/${OrthogroupID}.fa --model_config ${DeepFRI_path}/trained_models/model_config.json -ont cc -o ${temp_dir}/${OrthogroupID} &>/dev/null
        grep -v "^#" ${temp_dir}/${OrthogroupID}_CC_predictions.csv | awk '{FS=OFS=","}{if($3 >= 0.5){print $0}}' | sed 's/,/\t/;s/,/\t/;s/,/\t/' > ${DeepFRI_results}/${OrthogroupID}_CC_predictions.txt
        python ${DeepFRI_path}/predict.py --fasta_fn ${Orthogroups_dir}/${OrthogroupID}.fa --model_config ${DeepFRI_path}/trained_models/model_config.json -ont ec -o ${temp_dir}/${OrthogroupID} &>/dev/null
        grep -v "^#" ${temp_dir}/${OrthogroupID}_EC_predictions.csv | awk '{FS=OFS=","}{if($3 >= 0.5){print $0}}' | sed 's/,/\t/;s/,/\t/;s/,/\t/' > ${DeepFRI_results}/${OrthogroupID}_EC_predictions.txt
    done

    find ${DeepFRI_results} -size 0 -delete
}

run_foldseek() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/OrthogroupsAnnotation/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local foldseek_results="${working_dir}/Foldseek/"

    mkdir -p "${foldseek_results}"

    if [[ $foldseekdb == "pdb" ]]; then
        ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
        do
            foldseek easy-search ${Orthogroups_dir}/${OrthogroupID}.fa ${foldseek_path}/${foldseekdb} ${temp_dir}/${OrthogroupID}.m8 ${temp_dir}/tmp --prostt5-model ${foldseek_path}/weights -e 0.001 -c 0.5 --cov-mode 0 -v 0 --threads ${threads} &>/dev/null
            awk '{FS=OFS="\t"}{split($2,array,"-");$2=array[1];print}' ${temp_dir}/${OrthogroupID}.m8 > ${temp_dir}/${OrthogroupID}-2.m8
            join -t $'\t' -i -1 2 -2 1 <(sort -k2,2 ${temp_dir}/${OrthogroupID}-2.m8) ${foldseek_path}/entries_update.idx | awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -k1,1 > ${foldseek_results}/${OrthogroupID}.m8
        done
    elif [[ $foldseekdb == "afdb_swissprot" ]]; then
        ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
        do
            foldseek easy-search ${Orthogroups_dir}/${OrthogroupID}.fa ${foldseek_path}/${foldseekdb} ${temp_dir}/${OrthogroupID}.m8 ${temp_dir}/tmp --prostt5-model ${foldseek_path}/weights -e 0.001 -c 0.5 --cov-mode 0 -v 0 --threads ${threads} &>/dev/null
            awk '{FS=OFS="\t"}{split($2,array,"-");$2=array[2];print}' ${temp_dir}/${OrthogroupID}.m8 > ${temp_dir}/${OrthogroupID}-2.m8
            join -t $'\t' -i -1 2 -2 1 <(sort -k2,2 ${temp_dir}/${OrthogroupID}-2.m8) ${foldseek_path}/Accession_swissprot.txt | awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -k1,1 > ${foldseek_results}/${OrthogroupID}.m8
        done
    fi

    find ${foldseek_results} -size 0 -delete
}

run_hhblits() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/OrthogroupsAnnotation/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local hhblits_results="${working_dir}/hhblist/"

    mkdir -p "${hhblits_results}"

    ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
    do
        mafft --maxiterate 1000 --genafpair --thread ${threads} ${Orthogroups_dir}/${OrthogroupID}.fa > ${temp_dir}/${OrthogroupID}_aligned.fa 2>/dev/null
        hhblits -i ${temp_dir}/${OrthogroupID}_aligned.fa -o ${hhblits_results}/${OrthogroupID}.hhr -blasttab ${temp_dir}/${OrthogroupID}.txt -d ${hhsuite_path}/pfam -e 0.001 -n 6 -M 50 -z 2 -Z 10 -noprefilt -cpu ${threads} &>/dev/null
        awk '{FS=OFS="\t"}{if($11<=0.001){print}}' ${temp_dir}/${OrthogroupID}.txt > ${hhblits_results}/${OrthogroupID}.txt
    done

    find ${hhblits_results} -size 0 -delete
}

create_summary_table(){
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/OrthogroupsAnnotation/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local hhblits_results="${working_dir}/hhblist/"
    local foldseek_results="${working_dir}/Foldseek/"
    local DeepFRI_results="${working_dir}/DeepFRI/"
    local Summary_results="${working_dir}/Summary/"

    mkdir -p "${Summary_results}"
    echo "OrthogroupID;hhblits_domains;DeepFRI_MF;DeepFRI_BP;DeepFRI_CC;DeepFRI_EC;Foldseek" > ${working_dir}/General_summary.csv

    #Create summary for each one

    ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
    do
        local ProteinNumber=$(grep -c ">" ${Orthogroups_dir}/${OrthogroupID}.fa)

        echo "ProteinID;DeepFRI_MF;DeepFRI_BP;DeepFRI_CC;DeepFRI_EC;Foldseek" > ${Summary_results}/${OrthogroupID}.csv
        grep ">" ${Orthogroups_dir}/${OrthogroupID}.fa | awk -F '>' '{print $2}' | sort -n > ${temp_dir}/${OrthogroupID}.csv

        if [[ -f ${DeepFRI_results}/${OrthogroupID}_MF_predictions.txt ]]; then
            join -t ";" -a1 ${temp_dir}/${OrthogroupID}.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $2" "$4}' ${DeepFRI_results}/${OrthogroupID}_MF_predictions.txt | sed 's/[^[:print:]]$//' | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 -n | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-2.csv
            MF_general=$(awk -F '\t' '{print $2" "$4}' ${DeepFRI_results}/${OrthogroupID}_MF_predictions.txt | sed 's/[^[:print:]]$//' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
        else
            sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}.csv > ${temp_dir}/${OrthogroupID}-2.csv
            MF_general=""
        fi

        if [[ -f ${DeepFRI_results}/${OrthogroupID}_BP_predictions.txt ]]; then
            join -t ";" -a1 ${temp_dir}/${OrthogroupID}-2.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $2" "$4}' ${DeepFRI_results}/${OrthogroupID}_BP_predictions.txt | sed 's/[^[:print:]]$//' | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 -n | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-3.csv
            BP_general=$(awk 'BEGIN{FS="\t"}{print $2" "$4}' ${DeepFRI_results}/${OrthogroupID}_BP_predictions.txt | sed 's/[^[:print:]]$//' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
        else
            sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}-2.csv > ${temp_dir}/${OrthogroupID}-3.csv
            BP_general=""
        fi

        if [[ -f ${DeepFRI_results}/${OrthogroupID}_CC_predictions.txt ]]; then
            join -t ";" -a1 ${temp_dir}/${OrthogroupID}-3.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $2" "$4}' ${DeepFRI_results}/${OrthogroupID}_CC_predictions.txt | sed 's/[^[:print:]]$//' | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 -n | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-4.csv
            CC_general=$(awk 'BEGIN{FS="\t"}{print $2" "$4}' ${DeepFRI_results}/${OrthogroupID}_CC_predictions.txt | sed 's/[^[:print:]]$//' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
        else
            sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}-3.csv > ${temp_dir}/${OrthogroupID}-4.csv
            CC_general=""
        fi

        if [[ -f ${DeepFRI_results}/${OrthogroupID}_EC_predictions.txt ]]; then
            join -t ";" -a1 ${temp_dir}/${OrthogroupID}-4.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $2}' ${DeepFRI_results}/${OrthogroupID}_EC_predictions.txt | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 -n | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-5.csv
            EC_general=$(awk 'BEGIN{FS=OFS="\t"}{print $2}' ${DeepFRI_results}/${OrthogroupID}_EC_predictions.txt | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print $2}}' | tr -s '\n' ',' | sed 's/,$//')
        else
            sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}-4.csv > ${temp_dir}/${OrthogroupID}-5.csv
            EC_general=""
        fi

        if [[ -f ${foldseek_results}/${OrthogroupID}.m8 ]]; then
            if [[ $foldseekdb == "pdb" ]]; then
                join -t ";" -a1 ${temp_dir}/${OrthogroupID}-5.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $13}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u |sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 -n | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-6.csv
                Foldseek_general=$(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $13}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u | awk -F '\t' '{print $2}' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
            elif [[ $foldseekdb == "afdb_swissprot" ]]; then
                join -t ";" -a1 ${temp_dir}/${OrthogroupID}-5.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS "\""$13"\""}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 -n | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-6.csv
                Foldseek_general=$(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS "\""$13"\""}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u | awk -F '\t' '{print $2}' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
                #| sed -e 's/^/\"/' -e 's/$/\"/' 
            fi
        else
            sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}-5.csv > ${temp_dir}/${OrthogroupID}-6.csv
            Foldseek_general=""
        fi

        cat ${temp_dir}/${OrthogroupID}-6.csv >> ${Summary_results}/${OrthogroupID}.csv


        if [[ -f ${hhblits_results}/${OrthogroupID}.txt ]]; then
            awk -F '\t' '{print $2}' ${hhblits_results}/${OrthogroupID}.txt | sort -u > ${temp_dir}/${OrthogroupID}-hhblits.txt

            echo $OrthogroupID";"$(grep -f ${temp_dir}/${OrthogroupID}-hhblits.txt ${hhblits_results}/${OrthogroupID}.hhr | grep ">" | awk -F '>' '{print $2}' | sed -e 's/ ; /|/1' -e 's/ ; /-/' -e 's/|/ \"/' -e 's/$/\"/' | tr -s '\n' ',' | sed 's/,$//g')";"$MF_general";"$BP_general";"$CC_general";"$EC_general";"$Foldseek_general >> ${working_dir}/General_summary.csv
        else
            echo $OrthogroupID";"";"$MF_general";"$BP_general";"$CC_general";"$EC_general";"$Foldseek_general >> ${working_dir}/General_summary.csv
        fi
    done

}

# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running ClusterCharacterization Module."
check_clusters "${Working_directory}" "${clusters_file}"

cat ${clusters_file} | sed $'s/[^[:print:]\t]//g' | while read ClusterId
do
    internal_dir="${Working_directory}/Clusters/${ClusterId}/"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
    check_directory_structure "${internal_dir}"

    if $directory_flag; then

        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."
    
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
        organize_working_directory "${internal_dir}"

        # ==============================================================================
        # Running DeepFRI
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running DeepFRI..."
        run_deepfri "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

        # ==============================================================================
        # Running Foldseek
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running foldseek..."
        run_foldseek "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

        # ==============================================================================
        # Running hhblits
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running hhblits..."
        run_hhblits "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

        # ==============================================================================
        # Creating Summary table
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Creating summary table..."
        create_summary_table "${internal_dir}"
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding. \n"
    else
        echo -e "Cluster $ClusterId do not have the required directory for '$mode' mode. \n"
    fi
done
echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"