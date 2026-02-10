#!/usr/bin/env bash

# ==============================================================================
# SOFTWARE CHECK AND ENVIRONMENT SETUP
# ==============================================================================

# Define library and auxiliary paths
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Utils.sh"
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
    check_auxiliary_scripts "${auxiliary_path}" "$(basename -s .sh "$0" )"
fi

# Validate required software for the current script
check_required_software "$(basename -s .sh "$0" )"

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

function print_help() {
    echo -e "Script to identify orthogroups that are overrepresented in a specific dataset.
    This script perform three main steps:
    1. Run Orthofinder with DIAMOND ultra-sensitive mode.
    2. Remove orthogroups associated with captains.
    3. Perform the analysis depending on the selected mode:
      3.1. Outliers: Identify orthogroups that have an abnormal number of representative in the dataset using interquartile (IQR) upper fence [IQR = Q3 - Q1], depending on two approches for this kind of outlier identification:
        3.1.1. Standard: Identified orthogroups as outliers using as fence the following value: Q3 + \e[3mn\e[0m * IQR.
        3.1.2. Skew: Identify orthogroups as outliers using as fence the following value: Q3 + \e[3mn\e[0me^4MC * IQR.
      3.2. Enrichment: Identify orthogroups enriched in a group of elements based on qualitative variables in the metadata using the one-sided Fisher's exact test."
    echo ""
    echo "Syntax: StarCrew $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -m <string> { -a <string> -c <integer> -cm <string> | -n <string> -v <string> -p <float> } -t <integer> { --overwrite | --skip-orthofinder } ]"
    echo ""
    echo "Required args:"
    echo "-w, --workingDirectory: Specify the working directory where all data are stored."
    echo ""
    echo "Required args with Default:"
    echo "-m, --mode: Define the mode that will be used to flag orthogroups (Default: Outliers) [Available mode: Outliers, Enrichment]."
    echo ""
    echo "Required args with Default in 'Outliers' mode:"
    echo "-a, --approximation: Define the IQR approximation that is going to be used to defined outliers (Default: Skew) [Available mode: Standard, Skew]."
    echo "-c,--coefficient: In the case of 'IQR' rule, determine the coefficient for the fence definition (Default: 1.5) [range: 1 - 3]."
    echo "-cm, --countMode: Counts to be used for outliers (Default: Gene) [Available mode: Gene, Ship]."
    echo ""
    echo "Required args in 'Enrichment' mode:"
    echo "-n, --name: column name of the variable in the metadata file to be used."
    echo "-v, --value: value from the variable to be compare against the rest."
    echo ""
    echo "Required args with Default in 'Enrichment' mode:"
    echo "-p, --pValue: p-value to determined if an orthogroup is significantly enriched (Default: 0.05) [range: 0.001 - 0.1]."
    echo ""
    echo "Optional args:"
    echo "-t, --threads: Number of threads (Default: 8)"
    echo "--overwrite: Flag to overwrite in case there is already a previous run of $(basename -s .sh "$0" ). Not compatible wit '--skip-orthofinder' flag (Default: off)"
    echo "--skip-orthofinder: Flag to skip orthofinder in case a previous run was done and only want to change the mode or other value of the analysis. 
                              Not compatible with '--overwrite' flag (Default: off) "
    echo "-help: Display this help message."
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"
    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/"
    local temp_dir="${working_dir}/temp/"

    if [[ -d "$working_dir" ]]; then
        echo "Error: There's a previous run in the Workspace."
        echo "If you want to overwrite this previous run, add the '--overwrite' or '--skip-orthofinder' flag to the command line."
        exit 1
    fi

    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local metadata_dir=$(find "$base_dir" -maxdepth 1 -type d -name "metadata_files" 2>/dev/null)
    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)


    local Captain_dir=$(find "$base_dir" -maxdepth 1 -type d -name "CaptainIdentification" 2>/dev/null)
    if [[ ! -d "$Captain_dir" ]]; then
        echo "Error: CaptainIdentification results folder was not found in the working directory."
        exit 1
    else
        if [[ ! -f "${Captain_dir}/CaptainsID.txt" ]]; then
            echo "Error: No captain ID file was found."
            exit 1
        fi
    fi

    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    cp -r ${protein_dir} ${working_dir}
    cp -r ${gff_dir} ${working_dir}
    cp ${Captain_dir}/CaptainsID.txt ${working_dir}
    if [[ -f "${Captain_dir}/CaptainPhylogeny.nw" ]]; then
        cp ${Captain_dir}/CaptainPhylogeny.nw ${working_dir}
        captain_phylogeny=true
    else
        captain_phylogeny=false
    fi

    if [[ "${mode}" = "Enrichment" ]]; then
        if [[ -f "${metadata_dir}/metadata.csv" ]]; then
            cp "${metadata_dir}/metadata.csv" ${working_dir}
        else
            echo "Error: metadata file was not found."
            exit 1
        fi
    fi
}

run_orthofinder() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/"
    local protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local captainPhylogeny="${working_dir}/CaptainPhylogeny.nw"
    local output_dir="${working_dir}/Orthofinder"
    local temp_dir="${working_dir}/temp/"

    local element_number=$(ls ${protein_dir} | wc -l)

    echo ""

    if [[ $(ulimit -Sn) -lt "$((element_number + 124))" ]]; then
        echo -e "Error: The system limits on the number of files a process can open is probably too low (Current: '$(ulimit -Sn)'). Please increase it at least to '$((element_number + 124))' or more."
        exit 1
    fi  

    if $captain_phylogeny; then
        orthofinder -a "$(( ${threads} / 2 ))" -t "${threads}" -f ${protein_dir} -A mafft -S diamond_ultra_sens --matrix PAM30 -s ${captainPhylogeny} --scores-v2 -o ${output_dir} -n characterization
    else
        orthofinder -a "$(( ${threads} / 2 ))" -t "${threads}" -f ${protein_dir} -A mafft -S diamond_ultra_sens --matrix PAM30 --scores-v2 -o ${output_dir} -n characterization
        #&> ${working_dir}/orthofinder.log
    fi

    local results_path="${output_dir}/Results_characterization/Orthogroups/Orthogroups.GeneCount.tsv"

    if [ -f "${results_path}" ]; then
        orthofinder_flag=true
        cp ${results_path} ${working_dir}
        awk 'BEGIN {FS=OFS="\t"} 
         NR==1 { TOTAL_COLUMNS = NF; 
         print $0; 
         next}
         {ROW_TOTAL = 0;
         sub(/[[:space:]]+$/, "", $0); 
         printf "%s", $1; 
         for (i=2; i<=TOTAL_COLUMNS; i++) {
         is_present = (i <= NF) ? (($i != "") ? 1 : 0) : 0;
         ROW_TOTAL += is_present;
         printf "%s%d", OFS, is_present;
         }
         printf "%s%d\n", OFS, ROW_TOTAL;
        }' ${output_dir}/Results_characterization/Orthogroups/Orthogroups_UnassignedGenes.tsv | sed '1d' >> ${working_dir}/Orthogroups.GeneCount.tsv
        cp ${output_dir}/Results_characterization/Orthogroups/Orthogroups.txt ${working_dir}
    else
        orthofinder_flag=false
    fi
}

remove_captain_orthogroups() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/"
    local temp_dir="${working_dir}/temp/"
    local count_file="${working_dir}/Orthogroups.GeneCount.tsv"
    local gene_file="${working_dir}/Orthogroups.txt"
    local captain_file="${working_dir}/CaptainsID.txt"

    if [[ ! -f "${count_file}" || ! -f "${gene_file}" ]]; then
        echo "Error: require files from orthofinder are missing."
        exit 1
    fi

    if [[ ! -f "${captain_file}" ]]; then
        echo "Error: require Captain ID file is missing."
        exit 1
    fi

    grep -w -v -f ${captain_file} ${gene_file} > ${working_dir}/Orthogroups-captainless.txt
    awk -F ':' '{print $1}' ${working_dir}/Orthogroups-captainless.txt > ${temp_dir}/Orthogroups-captainlessID.txt
    head -n1 ${count_file} > ${working_dir}/Orthogroups.GeneCount-captainless.tsv
    grep -w -f ${temp_dir}/Orthogroups-captainlessID.txt ${count_file} >> ${working_dir}/Orthogroups.GeneCount-captainless.tsv
}

organize_information() {
    local base_dir="$1"
    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )"
    local orthogroups_dir="${working_dir}/Orthofinder/Results_characterization/Orthogroup_Sequences/"

    local Overrepresented_dir="${base_dir}/$(basename -s .sh "$0" )/"
    mkdir -p ${Overrepresented_dir}

    if [[ ! -d "${working_dir}" ]]; then
        echo "Error in working directory"
        exit 1
    fi

    if [[ -f "${working_dir}/OrthogroupsSizeHistogram.svg" ]]; then
        cp ${working_dir}/OrthogroupsSizeHistogram.svg ${Overrepresented_dir}
    fi

    if [[ -f "${working_dir}/Overrepresented_orthogroups.txt" ]]; then
        cp ${working_dir}/Overrepresented_orthogroups.txt ${Overrepresented_dir}
        mkdir -p ${Overrepresented_dir}/Orthogroups/
        awk '{print $1}' ${working_dir}/Overrepresented_orthogroups.txt | while read Orthogroup
        do 
            cp ${orthogroups_dir}/${Orthogroup}.fa ${Overrepresented_dir}/Orthogroups/
        done
    fi

    if [[ -f "${working_dir}/Enrichment_results.txt" ]]; then
        cp ${working_dir}/Enrichment_results.txt ${Overrepresented_dir}
    fi

    if [[ -f "${working_dir}/Enrich_orthogroups.txt" ]]; then
        cp ${working_dir}/Enrich_orthogroups.txt ${Overrepresented_dir}
        mkdir -p ${Overrepresented_dir}/Orthogroups/
        awk '{print $1}' ${working_dir}/Enrich_orthogroups.txt | while read Orthogroup
        do 
            cp ${orthogroups_dir}/${Orthogroup}.fa ${Overrepresented_dir}/Orthogroups/
        done
    fi
}

# ==============================================================================
# VARIABLES AND ARGUMENT PARSING
# ==============================================================================

# Initialize variables
Working_directory=""
mode="Outliers"
rule="Skew"
coefficient="1.5"
countmode="Gene"
column=""
value=""
pvalue="0.05"
threads="8"
overwrite=false
skip_orthofinder=false
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
        -a|--approximation)
            shift
            rule="$1"
            ;;
        -c|--coefficient)
            shift
            coefficient="$1"
            ;;
        -cm|--countMode)
            shift
            countmode="$1"
            ;;
        -n|--name)
            shift
            column="$1"
            ;;
        -v|--value)
            shift
            value="$1"
            ;;
        -p|--pValue)
            shift
            pvalue="$1"
            ;;
        -t|--threads)
            shift
            threads="$1"
            ;;
        --overwrite)
            overwrite=true
            ;;
        --skip-orthofinder)
            skip_orthofinder=true
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
echo "  Mode: " "$mode"
echo "  approximation for 'Outliers' mode: " "$rule"
echo "  Coeeficient for fence definition: " "$coefficient"
echo "  Count to use for 'Outliers' mode: " "$countmode"
echo "  Column in metadata for 'Enrichment' mode: " "$column"
echo "  Value in metadata for 'Enrichment' mode: " "$value"
echo "  p-Value in 'Enrichment' mode:" "$pvalue"
echo "  Overwrite: " "$overwrite"
echo "  Skip orthofinder: " "$skip_orthofinder"
echo "  Number of threads: " "$threads"
echo ""

# ==============================================================================
# ARGUMENTS AND INPUT CHECK
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking arguments and input files..."

check_mode_parameter "$mode" "$(basename -s .sh "$0" )"

if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
else
    if [[ ! -d "$Working_directory" ]]; then
        echo "Error: folder '$Working_directory' does not exist."
        exit 1
    else
        if [[ ! "${Working_directory:0:1}" == "/" ]]; then
            Working_directory=$(realpath $Working_directory)
            check_directory_structure "${Working_directory}"
        fi
    fi
fi

if [[ "$mode" == "Enrichment" ]]; then
    if [[ -z "$column" || -z "$value" ]]; then
        echo "Error: Missing required arguments for 'Enrichment' mode."
        print_help
        exit 1
    fi

    if [[ "$pvalue" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
        if (( $(echo "$pvalue < 0.001" | bc -l) )) || (( $(echo "$pvalue > 0.1" | bc -l) )); then
            echo "Error: '$pvalue' is not an accepted value for p-value."
            print_help
            exit 1
        fi
    else
        echo "Error: '$pvalue' is not a float."
        print_help
        exit 1
    fi

    if [[ ! -f "${Working_directory}/metadata_files/metadata.csv" ]]; then
        echo "Error: metadata file wasn't found in the working folder"
        exit 1
    fi

    column_number=$(sed -n "1 s/${column}.*//p" ${Working_directory}/metadata_files/metadata.csv | sed 's/[^;*]//g' | wc -c)
    if [[ $column_number -eq 0 ]]; then
        echo "Error: column '$column' do not exist in metadata."
        exit 1
    elif [[ $column_number -le 2 ]]; then
        echo "Error: the column to analyzed cannot be the ID column."
        exit 1
    elif [[ $column_number -gt 2 ]]; then
        if [[ "$value" == "NA" ]]; then
            echo "Error: 'Enrichment' mode do not accept to analyze 'NA' as target value"
        else
            value_number=$(cut -d ";" -f $column_number ${Working_directory}/metadata_files/metadata.csv | grep -c "$value")

            if [[ $value_number -eq 0 ]]; then
                echo "Error: provide value '$value' do not exist in the column '$column' of metadata."
                exit 1
            elif [[ $value_number -lt 10 ]]; then
                echo "Error: provide value '$value' have an n of '$value_number' and it's too low for enrichment analysis."
                exit 1
            fi
        fi
    fi
fi

if [[ "$mode" == "Outliers" ]]; then
    if [[ "$rule" != "Standard" && "$rule" != "Skew" ]]; then
        echo "Error: '$rule' approximation is not a valid value"
        print_help
        exit 1
    fi

    if [[ "$countmode" != "Gene" && "$countmode" != "Ship" ]]; then
        echo -e "Error: '$countmode' count mode is not a valid value.\n"
        print_help
        exit 1
    fi

    if [[ "$coefficient" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
        if (( $(echo "$coefficient < 1.5" | bc -l) )) || (( $(echo "$coefficient > 3.0" | bc -l) )); then
            echo "Error: '$coefficient' is not an accepted value for the fence coefficient."
            print_help
            exit 1
        fi
    else
        echo "Error: '$coefficient' is not a float."
        print_help
        exit 1
    fi
fi

if $overwrite && $skip_orthofinder; then
    echo "Error: '--overwrite' and '--skip-orthofinder' flags are both 'on' and those are incompatible."
    exit 1
fi

check_threads "$threads"

# ==============================================================================
# MAIN SCRIPT
# ==============================================================================

if $overwrite; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking and removing previous run if exists..."
    overwrite "${Working_directory}" "$(basename -s .sh "$0" )" "${mode}"
fi

check_directory_structure "$Working_directory"

if ! $skip_orthofinder; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
    organize_working_directory "${Working_directory}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running orthofinder..."
    run_orthofinder "${Working_directory}"
else
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Skipping the working directory organization and orthofinder run."
fi

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Removing orthogroups associated with Captains..."
remove_captain_orthogroups "${Working_directory}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Performing analysis..."
Rscript ${auxiliary_path}/OverrepresentationAnalysis.R -d "${Working_directory}/Workspace/$(basename -s .sh "$0" )" -m "${mode}" -a "${rule}" -c "${coefficient}" -n "${column}" -v "${value}" -p "${pvalue}" -cm "${countmode}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing information to the main directory."
organize_information "${Working_directory}"