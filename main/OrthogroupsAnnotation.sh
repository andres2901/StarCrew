#!/usr/bin/env bash

# ==============================================================================
# SOFTWARE CHECK AND ENVIRONMENT SETUP
# ==============================================================================

# Define library
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Utils.sh"
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

# Validate required software for the current script
check_required_software "$(basename -s .sh "$0" )"

# Database Verification
database_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../databases/"
check_databases "${database_path}" "$(basename -s .sh "$0" )"
foldseek_path="${database_path}/Foldseek/"
hhsuite_path="${database_path}/hhsuite/"

Interpro_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../interproscan/"
check_interpro_software "${Interpro_path}"
Interpro_path=$(realpath $Interpro_path)

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

# Function to print help message
function print_help() {
   echo -e "Script to run a functional annotation for orthogroups.
   This script perform three steps:
   1. Organize the Orthogroups that are selected base on the mode.
     1.1 Core: Orthogoups that were identify as core by the ClusterCharacteriation command.
     1.2 MoveAssociated: Orthogoups that were identify as part of a putative movement event between subclusters by the ClusterCharacterization command.
     2.3 All: All orthogroups identify by the ClusterCharacterization command.
     2.4 Overrepresented: Orthogroups that were identified as Overrepresented by the OrthogroupsOverrepresentation command.
   2. Perform the characterization of the Orthogroup proteins with four approaches:
     2.1 InterProScan: Using all default applications except COILS and MOBIDB.
     2.3 Foldseek: Search for homologs proteins against a database based on the 3D structure.
     2.4 hhblits: Search domains against the PfamA database.
   3. Summarize the results of the previous step:
     3.1 Internal summary: For each Orthogroups summarize the results per protein in a csv
     3.2 General summary: Return a summary for the Orthogroup under the assumption that all proteins in each Orthogroups have the same function.
     It return only those 'chracteristics' that are shared for at least 50% of the proteins in the Orthogroup.
   "
   echo
   echo "Syntax: StarCrew $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -m <string> -f <string> -t <integer> --overwrite ]"
   echo ""
   echo "Required args:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored."
   echo ""
   echo "Required args with Default:"
   echo "-m, --mode: Define the orthogroups to be analyzed (Default = All) [Available mode: MoveAssociated, Core, All, Overrepresented]."
   echo "-f, --foldseekdb: Name of the Foldseek database to use (Default = afdb_swissprot) [Available: pdb, afdb_swissprot]."
   echo "-t, --threads: Number of threads for all analysis (Default: 8)."
   echo ""
   echo "Optional args:"
   echo "--overwrite: Flag to overwrite in case there is already a previous run of $(basename -s .sh "$0" ) (Default: off)."
   echo "-help: Display this help message."
}

check_clusters() {
    local base_dir="$1"
    local Cluster_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Clusters" 2>/dev/null)

    if [[ $mode == "All" ]]; then
        if [[ ! -f "${Cluster_dir}/ClusterOrthogroups.txt" ]]; then
            echo "Error: ClusterOrthogroups.txt file is not find in '$Cluster_dir'." >&2
            exit 1
        else
            clusters_file="${Cluster_dir}/ClusterOrthogroups.txt"
        fi
    elif [[ $mode == "MoveAssociated" ]]; then
        if [[ ! -f "${Cluster_dir}/ClusterMovement.txt" ]]; then
            echo "Error: ClusterMovement.txt file is not find in '$Cluster_dir'." >&2
            exit 1
        else
            clusters_file="${Cluster_dir}/ClusterMovement.txt"
        fi
    elif [[ $mode == "Core" ]]; then
        if [[ ! -f "${Cluster_dir}/ClusterCore.txt" ]]; then
            echo "Error: ClusterCore.txt file is not find in '$Cluster_dir'." >&2
            exit 1
        else
            clusters_file="${Cluster_dir}/ClusterCore.txt"
        fi
    fi

    Cluster_number=$(wc -l "$clusters_file" | awk '{print $1}')
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing '${Cluster_number}' clusters.\n"
}

check_overrepresentation() {
    local base_dir="$1"
    local overrepresentation_dir=$(find "$base_dir" -maxdepth 1 -type d -name "OrthogroupsOverrepresentation" 2>/dev/null)

    if [[ ! -d "${overrepresentation_dir}" ]]; then
        echo "Error: 'OrthogroupsOverrepresentation' folder do not exist. Please run 'OrthogroupsOverrepresentation' command." >&2
        exit 1
    else
        if [[ ! -d "${overrepresentation_dir}/Orthogroups/" ]]; then
            echo "[$(date "+%Y-%m-%d %H:%M:%S")] There's no 'Orthogroups' folder in 'overrepresentation_dir' folder. Check if there is any Overrepresented Orthogroup in this dataset."
            exit 0
        fi
    fi
}

check_internal_directory_structure() {
    local base_dir="$1"
    
    # Locate required subdirectories and file
    local workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    fi

    local ClusterAnnotation_dir=$(find "$workspace_dir" -maxdepth 1 -type d -name "ClusterCharacterization" 2>/dev/null)

    if [[ -z "$ClusterAnnotation_dir" ]]; then
        echo "Error: ClusterCharacterization directory not found in '$workspace_dir'." >&2
        echo "The ClusterCharacterization command should be run before this analysis" >&2
        exit 1
    fi
    
    local moveOrthologs_dir=$(find "$ClusterAnnotation_dir" -maxdepth 1 -type d -name "*_moveOrthologs" 2>/dev/null)
    local CoreGenes_dir=$(find "$ClusterAnnotation_dir" -maxdepth 1 -type d -name "Core_genes" 2>/dev/null)
    local Orthogroups_dir=$(find "$ClusterAnnotation_dir" -maxdepth 3 -type d -name "Orthogroup_Sequences" 2>/dev/null)

    directory_flag=true

    if [[ $mode == "All" ]]; then
        if [[ -z "$Orthogroups_dir" ]]; then
            echo "Error: Orthogroups subdirectory not found in '$ClusterAnnotation_dir'." >&2
            directory_flag=false
        fi
    elif [[ $mode == "MoveAssociated" ]]; then
        if [[ -z "$moveOrthologs_dir" ]]; then
            echo "Orthogroups subdirectory for moving genes not found in '$ClusterAnnotation_dir'." >&2
            directory_flag=false
        fi
    elif [[ $mode == "Core" ]]; then
        if [[ -z "$CoreGenes_dir" ]]; then
            echo "Error: Orthogroups subdirectory not found in '$ClusterAnnotation_dir'." >&2
            directory_flag=false
        fi
    fi
}

organize_working_directory() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"
    local ClusterAnnotation_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"

    if [[ -d "$working_dir" ]]; then
        echo "Error: There's a previous run in the Workspace."
        echo "If you want to overwrite this previous run, add the '--overwrite' flag to the command line."
        exit 1
    fi

    mkdir -p ${working_dir}
    mkdir -p ${Orthogroups_dir}
    mkdir -p ${temp_dir}

    if [[ $mode == "All" ]]; then
        local data_dir=$(find "$ClusterAnnotation_dir" -maxdepth 3 -type d -name "Orthogroup_Sequences" 2>/dev/null)
        cp ${data_dir}/* ${Orthogroups_dir}/
    elif [[ $mode == "MoveAssociated" ]]; then
        local moveOrthologs_dir=$(find "$ClusterAnnotation_dir" -maxdepth 1 -type d -name "*_moveOrthologs" 2>/dev/null)
        echo $moveOrthologs_dir | awk '{OFS=RS;$1=$1}1' | while read line; 
        do 
            cp ${line}/* ${Orthogroups_dir}/
        done
    elif [[ $mode == "Core" ]]; then
        local CoreGenes_dir=$(find "$ClusterAnnotation_dir" -maxdepth 1 -type d -name "Core_genes" 2>/dev/null)
        cp ${CoreGenes_dir}/* ${Orthogroups_dir}/
    fi
}

organize_working_directory_overrepresentation() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"
    local overrepresented_dir=$(find "$base_dir" -maxdepth 1 -type d -name "OrthogroupsOverrepresentation" 2>/dev/null)
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"

    if [[ -d "$working_dir" ]]; then
        echo "Error: There's a previous run in the Workspace."
        echo "If you want to overwrite this previous run, add the '--overwrite' flag to the command line."
        exit 1
    fi

    mkdir -p ${working_dir}
    mkdir -p ${Orthogroups_dir}
    mkdir -p ${temp_dir}

    local data_dir=$(find "$overrepresented_dir" -maxdepth 3 -type d -name "Orthogroups" 2>/dev/null)
    cp ${data_dir}/* ${Orthogroups_dir}/
}

run_foldseek() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local foldseek_results="${working_dir}/Foldseek/"

    mkdir -p "${foldseek_results}"

    Total_states=$(ls ${Orthogroups_dir} | wc -l)
    State=0

    if [[ $foldseekdb == "pdb" ]]; then
        ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
        do
            State=$(($State + 1))
            foldseek easy-search ${Orthogroups_dir}/${OrthogroupID}.fa ${foldseek_path}/${foldseekdb} ${temp_dir}/${OrthogroupID}.m8 ${temp_dir}/tmp --prostt5-model ${foldseek_path}/weights -e 0.001 -c 0.5 --cov-mode 0 -v 0 --threads ${threads} &>/dev/null
            awk '{FS=OFS="\t"}{split($2,array,"-");$2=array[1];print}' ${temp_dir}/${OrthogroupID}.m8 > ${temp_dir}/${OrthogroupID}-2.m8
            join -t $'\t' -i -1 2 -2 1 <(sort -k2,2 ${temp_dir}/${OrthogroupID}-2.m8) ${foldseek_path}/entries_update.idx | awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -k1,1 | sed -e 's/;/,/g' > ${foldseek_results}/${OrthogroupID}.m8
            ProgressBar $State $Total_states
        done
    elif [[ $foldseekdb == "afdb_swissprot" ]]; then
        ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
        do
            State=$(($State + 1))
            foldseek easy-search ${Orthogroups_dir}/${OrthogroupID}.fa ${foldseek_path}/${foldseekdb} ${temp_dir}/${OrthogroupID}.m8 ${temp_dir}/tmp --prostt5-model ${foldseek_path}/weights -e 0.001 -c 0.5 --cov-mode 0 -v 0 --threads ${threads} &>/dev/null
            awk '{FS=OFS="\t"}{split($2,array,"-");$2=array[2];print}' ${temp_dir}/${OrthogroupID}.m8 > ${temp_dir}/${OrthogroupID}-2.m8
            join -t $'\t' -i -1 2 -2 1 <(sort -k2,2 ${temp_dir}/${OrthogroupID}-2.m8) ${foldseek_path}/Accession_swissprot.txt | awk 'BEGIN{FS=OFS="\t"}{swap=$1;$1=$2;$2=swap;print $0}' | sort -k1,1 | sed -e 's/;/,/g' > ${foldseek_results}/${OrthogroupID}.m8
            ProgressBar $State $Total_states
        done
    fi

    echo ""

    find ${foldseek_results} -size 0 -delete
}

run_hhblits() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local hhblits_results="${working_dir}/hhblist/"

    mkdir -p "${hhblits_results}"

    Total_states=$(ls ${Orthogroups_dir} | wc -l)
    State=0

    ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
    do
        State=$(($State + 1))
        mafft --maxiterate 1000 --genafpair --reorder --thread ${threads} ${Orthogroups_dir}/${OrthogroupID}.fa > ${temp_dir}/${OrthogroupID}_aligned.fa 2>/dev/null
        hhblits -i ${temp_dir}/${OrthogroupID}_aligned.fa -o ${hhblits_results}/${OrthogroupID}.hhr -blasttab ${temp_dir}/${OrthogroupID}.txt -d ${hhsuite_path}/pfam -e 0.001 -n 6 -M 50 -z 2 -Z 10 -realign_old_hits -cov 50 -cpu ${threads} &>/dev/null
        if [[ -f ${temp_dir}/${OrthogroupID}.txt ]]; then
            awk '{FS=OFS="\t"}{if($11<=0.001){print}}' ${temp_dir}/${OrthogroupID}.txt > ${hhblits_results}/${OrthogroupID}.txt
        fi
        ProgressBar $State $Total_states
    done

    echo ""

    find ${hhblits_results} -size 0 -delete
}

run_interproScan() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local interpro_results="${working_dir}/InterProScan/"

    mkdir -p "${interpro_results}"

    Total_states=$(ls ${Orthogroups_dir} | wc -l)
    State=0

    ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
    do
        State=$(($State + 1))
        ${Interpro_path}/interproscan.sh -i ${Orthogroups_dir}/${OrthogroupID}.fa -f tsv -appl CDD,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PRINTS,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM --goterms -o ${interpro_results}/${OrthogroupID}.tsv &>/dev/null
        ProgressBar $State $Total_states
    done

    echo ""

    find ${interpro_results} -size 0 -delete
}

create_summary_table(){
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"
    local Orthogroups_dir="${working_dir}/Orthogroups/"
    local temp_dir="${working_dir}/temp/"
    local interpro_results="${working_dir}/InterProScan/"
    local hhblits_results="${working_dir}/hhblist/"
    local foldseek_results="${working_dir}/Foldseek/"
    local DeepFRI_results="${working_dir}/DeepFRI/"
    local Summary_results="${working_dir}/Summary/"

    mkdir -p "${Summary_results}"
    echo "OrthogroupID;hhblits_domains;InterProScan;InterProScan_GO;Foldseek" > ${working_dir}/General_summary.csv

    #Create summary for each one

    ls ${Orthogroups_dir} | xargs -n 1 basename -s .fa | while read OrthogroupID 
    do
        local ProteinNumber=$(grep -c ">" ${Orthogroups_dir}/${OrthogroupID}.fa)

        echo "ProteinID;InterProScan;InterProScan_GO;Foldseek" > ${Summary_results}/${OrthogroupID}.csv
        grep ">" ${Orthogroups_dir}/${OrthogroupID}.fa | awk -F '>' '{print $2}' | sort -k1,1 > ${temp_dir}/${OrthogroupID}.csv

        if [[ -f ${interpro_results}/${OrthogroupID}.tsv ]]; then
            join -t ";" -a1 ${temp_dir}/${OrthogroupID}.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $5"-\""$6"\""}' ${interpro_results}/${OrthogroupID}.tsv | sed -e 's/[^[:print:]]$//' -e 's/;/,/g' | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($2==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-2.csv
            InterPro_General=$(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $5"-\""$6"\""}' ${interpro_results}/${OrthogroupID}.tsv | sed -e 's/[^[:print:]]$//' -e 's/;/,/g' | sort -k1,2 -u | awk -F '\t' '{print $2}'  | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
            
            awk 'BEGIN{FS=OFS="\t"}{if($14!="-"){print $1 OFS $14}}' ${interpro_results}/${OrthogroupID}.tsv > ${temp_dir}/${OrthogroupID}-IPS.tsv
            if [[ -s ${temp_dir}/${OrthogroupID}-IPS.tsv ]]; then
                join -t ";" -a1 ${temp_dir}/${OrthogroupID}-2.csv <(sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' ${temp_dir}/${OrthogroupID}-IPS.tsv | sort -k1,1 | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($3==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-3.csv
                InterProGO_General=$(sort -k1,2 -u ${temp_dir}/${OrthogroupID}-IPS.tsv | awk -F '\t' '{print $2}'  | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
            else
                sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}-2.csv > ${temp_dir}/${OrthogroupID}-3.csv
                InterProGO_General=""
            fi
        else
            sed -e 's/$/;;/g' ${temp_dir}/${OrthogroupID}.csv > ${temp_dir}/${OrthogroupID}-3.csv
            InterPro_General=""
            InterProGO_General=""
        fi

        if [[ -f ${foldseek_results}/${OrthogroupID}.m8 ]]; then
            if [[ $foldseekdb == "pdb" ]]; then
                join -t ";" -a1 ${temp_dir}/${OrthogroupID}-3.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $13}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u |sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($8==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-4.csv
                Foldseek_general=$(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $13}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u | awk -F '\t' '{print $2}' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
            elif [[ $foldseekdb == "afdb_swissprot" ]]; then
                join -t ";" -a1 ${temp_dir}/${OrthogroupID}-3.csv <(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS "\""$13"\""}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1,/;t1;P;D' | sort -k1,1 | sed 's/\t/;/g') | awk 'BEGIN{FS=OFS=";"}{if($8==""){print $0";"}else{print $0}}' > ${temp_dir}/${OrthogroupID}-4.csv
                Foldseek_general=$(awk 'BEGIN{FS=OFS="\t"}{print $1 OFS "\""$13"\""}' ${foldseek_results}/${OrthogroupID}.m8 | sort -k1,2 -u | awk -F '\t' '{print $2}' | sort | uniq -c | sed 's/^ *//g' | awk -v Num=$ProteinNumber '{if($1>=(Num*0.5)){print}}' | cut -d ' ' -f2- | tr -s '\n' ',' | sed 's/,$//')
            fi
        else
            sed -e 's/$/;/g' ${temp_dir}/${OrthogroupID}-3.csv > ${temp_dir}/${OrthogroupID}-4.csv
            Foldseek_general=""
        fi

        cat ${temp_dir}/${OrthogroupID}-4.csv >> ${Summary_results}/${OrthogroupID}.csv


        if [[ -f ${hhblits_results}/${OrthogroupID}.txt ]]; then
            awk -F '\t' '{print $2}' ${hhblits_results}/${OrthogroupID}.txt | sort -u > ${temp_dir}/${OrthogroupID}-hhblits.txt

            echo $OrthogroupID";"$(grep -f ${temp_dir}/${OrthogroupID}-hhblits.txt ${hhblits_results}/${OrthogroupID}.hhr | grep ">" | awk -F '>' '{print $2}' | sort -u | sed -e 's/ ; /|/1' -e 's/ ; /-/' -e 's/|/ \"/' -e 's/$/\"/' | tr -s '\n' ',' | sed 's/,$//g')";"$InterPro_General";"$InterProGO_General";"$Foldseek_general >> ${working_dir}/General_summary.csv
        else
            echo $OrthogroupID";"";"$InterPro_General";"$InterProGO_General";"$Foldseek_general >> ${working_dir}/General_summary.csv
        fi
    done
}

organize_information() {
    local base_dir="$1"
    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/"

    local Annotation_dir="${base_dir}/$(basename -s .sh "$0" )-${mode}/"

    mkdir -p ${Annotation_dir}

    if [[ ! -d "${working_dir}" ]]; then
        echo "Error in working directory"
    fi

    if [[ -f "${working_dir}/General_summary.csv" ]]; then
        cp ${working_dir}/General_summary.csv ${Annotation_dir}/
    fi

    if [[ -d "${working_dir}/Summary/" ]]; then
        cp -r ${working_dir}/Summary/ ${Annotation_dir}/Orthogroups_summary
    fi
}

# ==============================================================================
# VARIABLES AND ARGUMENT PARSING
# ==============================================================================

# Initialize variables
Working_directory=""
mode="All"
foldseekdb="afdb_swissprot"
threads="8"
overwrite=false
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
        -t|--threads)
            shift
            threads="$1"
            ;;
        -f|--foldseekdb)
            shift
            foldseekdb="$1"
            ;;
        --overwrite)
            overwrite=true
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

if $help_flag; then
    print_help
    exit 0
fi

echo "Running $(basename -s .sh "$0" ) command under the following parameters:"
echo "  Working directory: " "$Working_directory"
echo "  Cluster file: " "$clusters_file"
echo "  Mode: " "$mode"
echo "  Foldseek database: " "$foldseekdb"
echo "  Threads: " "$threads"
echo "  Overwrite previous run: " "$overwrite"
echo ""

# ==============================================================================
# ARGUMENTS AND INPUT CHECK
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking arguments and input files..."

check_mode_parameter "${mode}" "$(basename -s .sh "$0" )"

if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
else
    if [[ ! -d "$Working_directory" ]]; then
        echo "Error: Directory '$Working_directory' does not exist."
        exit 1
    else
        Working_directory=$(realpath $Working_directory)
        check_directory_structure "${Working_directory}"
    fi
fi

check_clusters "${Working_directory}"

if [[ "$foldseekdb" != "pdb" && "$foldseekdb" != "afdb_swissprot" ]]; then
    echo "Error: '$foldseekdb' is not accepted as a foldseek database."
    print_help
    exit 1
else
    check_foldseek_databases "${foldseek_path}" "${foldseekdb}"
fi

check_threads "${threads}"

# ==============================================================================
# MAIN SCRIPT
# ==============================================================================

if [[ "$mode" == "All" || "$mode" == "MoveAssociated" || "$mode" == "Core" ]]; then
    awk '{print $1}' ${clusters_file} | sed $'s/[^[:print:]\t]//g' | while read ClusterId
    do
        internal_dir="${Working_directory}/Clusters/${ClusterId}/"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
        check_internal_directory_structure "${internal_dir}"

        if $directory_flag; then
            echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."
            if $overwrite; then
                echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking and removing previous run if exists..."
                overwrite "${internal_dir}" "$(basename -s .sh "$0" )" "${mode}"
            fi
    
            echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
            organize_working_directory "${internal_dir}"

            echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running InterProScan..."
            run_interproScan "${internal_dir}"
            echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

            echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running foldseek..."
            run_foldseek "${internal_dir}"
            echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

            echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running hhblits..."
            run_hhblits "${internal_dir}"
            echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

            echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Creating summary table..."
            create_summary_table "${internal_dir}"
            rm -r "${internal_dir}/Workspace/$(basename -s .sh "$0" )-${mode}/temp/" 2> /dev/null
            echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing information to the main directory."
            organize_information "${internal_dir}"
            echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding. \n"
        else
            echo -e "Cluster $ClusterId do not have the required directory for '$mode' mode. \n"
        fi
    done
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"

elif [[ ${mode} == "Overrepresented" ]]; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${Working_directory}' structure."
    check_overrepresentation "${Working_directory}"

    if $overwrite; then
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking and removing previous run if exists..."
        overwrite "${Working_directory}" "$(basename -s .sh "$0" )" "${mode}"
    fi

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
    organize_working_directory_overrepresentation "${Working_directory}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running InterProScan..."
    run_interproScan "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running foldseek..."
    run_foldseek "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running hhblits..."
    run_hhblits "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Creating summary table..."
    create_summary_table "${Working_directory}"
    rm -r "${Working_directory}/Workspace/$(basename -s .sh "$0" )-${mode}/temp/" 2> /dev/null
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing information to the main directory."
    organize_information "${Working_directory}"
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding. \n"
fi
