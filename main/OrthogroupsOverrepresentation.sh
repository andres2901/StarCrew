#!/usr/bin/env bash

# ==============================================================================
# Software check block
# ==============================================================================

source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Utils.sh"
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

# ==============================================================================
# Function block
# ==============================================================================

# Function to print help message
function print_help() {
    echo "Script to identify orthogroups that are overrepresented in the a specific dataset.
    This script perform three steps:
    1. Run earlgrey TE prediction.
    2. Organize the results."
    echo
    echo "Syntax: StarClust $(basename -s .sh "$0" ) [ -help ] -w <directory_path> -d <file_path> [ -m <string> -t <integer> ]"
    echo ""
    echo "Required args:"
    echo "-w, --workingDirectory: Specify the working directory where all data are stored."
    echo ""
    echo "Required args with Default:"
    echo "-m, --mode: Define the mode that will be used to flag orthogroups (Default = Outliers) [Available mode: Outliers, Enrichment]."
    echo ""
    echo "Required args in 'Outliers' mode with Default:"
    echo "-r, --rule: Define the IQR rule that is going to be used to defined outliers (Default = Standard) [Available mode: Standard, Skew, Percentile]."
    echo "-p, --percentile: In the case of percentile rule, the percentile that will be used (Default = 95) [Range: 90 - 99]"
    echo ""
    echo "Required args in 'Enrichment' mode:"
    echo "-c, --column: column of the metadata file to be used."
    echo "-v, --value: value from the column to be compare against the rest."
    echo ""
    echo "Optional args:"
    echo "-t, --threads: Number of threads (Default = 8)"
    echo "-help: Display this help message."
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/"
    local temp_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/temp/"

    # Create required subdirectories
    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local metadata_dir=$(find "$base_dir" -maxdepth 1 -type d -name "metadata_files" 2>/dev/null)

    local Captain_dir=$(find "$base_dir" -maxdepth 1 -type d -name "CaptainIdentification" 2>/dev/null)
    if [[ ! -d "$Captain_dir" ]]; then
        echo "Error: CaptainIdentification results folder was not found in the working directory."
        exit 1
    else
        if [[ ! -f "${Captain_dir}/CaptainsID.txt" ]]; then
            echo "Error: No captain ID file was found."
            exit 1
        fi
        if [[ ! -f "${Captain_dir}/CaptainPhylogeny.nw" ]]; then
            echo "Error: No captain phylogeny file was found."
            exit 1
        fi
    fi

    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    cp -r ${protein_dir} ${working_dir}
    cp ${Captain_dir}/CaptainsID.txt ${working_dir}
    cp ${Captain_dir}/CaptainPhylogeny.nw ${working_dir}

    if [[ -f "${metadata_dir}/metadata.csv" ]]; then
        cp "${metadata_dir}/metadata.csv" ${working_dir}
        metadata_flag=true
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

    if [[ ${element_number} -gt 500 ]]; then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Cluster is too big. Dividing into two for Orthofinder run.."
        ulimit -n $((element_number + 128))
        gotree prune -i ${captainPhylogeny} --random "$(( $element_number / 2 ))" | gotree reroot midpoint -o ${temp_dir}/Captain_subtree.nw
        gotree stats tips -i ${temp_dir}/Captain_subtree.nw | awk 'NR>1{print $4}' > ${temp_dir}/Tips.txt
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Organizing information.."
        mkdir ${temp_dir}/protein1 ${temp_dir}/protein2
        ls ${protein_dir} | grep -f ${temp_dir}/Tips.txt | while read line
        do
            cp ${protein_dir}/${line} ${temp_dir}/protein1
        done
        ls ${protein_dir} | grep -v -f ${temp_dir}/Tips.txt | while read line
        do
            cp ${protein_dir}/${line} ${temp_dir}/protein2
        done

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] First Orthofinder run.."
        orthofinder -a "$(( ${threads} / 2 ))" -t "${threads}" -f ${temp_dir}/protein1 -A mafft -S diamond -I 4 -T iqtree3 --matrix PAM30 -s ${temp_dir}/Captain_subtree.nw --scores-v2 -o ${output_dir} -n Initial &> ${working_dir}/orthofinder1.log

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Second Orthofinder run.."
        orthofinder -a "$(( ${threads} / 2 ))" -t "${threads}" -A mafft -S diamond -I 4 -T iqtree3 --matrix PAM30 -s ${captainPhylogeny} --scores-v2 --assign ${temp_dir}/protein2 --core ${output_dir}/Results_Initial -n characterization &> ${working_dir}/orthofinder2.log
    else
        orthofinder -a "$(( ${threads} / 2 ))" -t "${threads}" -f ${protein_dir} -A mafft -S diamond -I 4 -T iqtree3 --matrix PAM30 -s ${captainPhylogeny} --scores-v2 -o ${output_dir} -n characterization &> ${working_dir}/orthofinder.log
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
        cp ${output_dir}/Results_characterization/Orthogroups/Orthogroups.txt ${working_dir}/
    else
        orthofinder_flag=false
    fi
}

remove_captain_orthogroups() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/"
    local temp_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/temp/"
    local count_file="${working_dir}/Orthogroups.GeneCount.tsv"
    local gene_file="${working_dir}/Orthogroups.txt"
    local captain_file="${working_dir}/CaptainsID.txt"

    grep -w -v -f ${captain_file} ${gene_file} > ${working_dir}/Orthogroups-captainless.txt
    awk -F ':' '{print $1}' ${working_dir}/Orthogroups-captainless.txt > ${temp_dir}/Orthogroups-captainlessID.txt
    head -n1 ${count_file} > ${working_dir}/Orthogroups.GeneCount-captainless.tsv
    grep -w -f ${temp_dir}/Orthogroups-captainlessID.txt ${count_file} >> ${working_dir}/Orthogroups.GeneCount-captainless.tsv
}

# ==============================================================================
# Variables block
# ==============================================================================

# Initialize variables
Working_directory=""
mode="Outliers"
rule="Standard"
percentile="95"
column=""
value=""
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
        -r|--rule)
            shift
            rule="$1"
            ;;
        -p|--percentile)
            shift
            percentile="$1"
            ;;
        -c|--column)
            shift
            column="$1"
            ;;
        -v|--value)
            shift
            value="$1"
            ;;
        -t|--threads)
            shift
            threads="$1"
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

# ==============================================================================
# Script start block
# ==============================================================================

echo "Running $(basename -s .sh "$0" ) command under the following parameters:"
echo "  Working directory: " "$Working_directory"
echo "  Mode: " "$mode"
echo "  Rule for Outliers mode: " "$rule"
echo "  Percentile: " "$percentile"
echo "  Column in metadata for 'Enrichment' mode: " "$column"
echo "  Value in metadata for 'Enrichment' mode: " "$value"
echo "  Number of threads: " "$threads"
echo ""

# ==============================================================================
# Check variables block
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking arguments and input files..."

# Check for mandatory arguments
if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# check if working directory exist
if [[ ! -d "$Working_directory" ]]; then
    echo "Error: folder '$Working_directory' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${Working_directory:0:1}" == "/" ]]; then
        Working_directory=$(realpath $Working_directory)
    fi
fi
check_directory_structure "${Working_directory}"

# Check mode parameter
check_mode_parameter "$mode" "$(basename -s .sh "$0" )"

# Check parameters in 'Enrichment' mode
if [[ "$mode" == "Enrichment" ]]; then
    if [[ -z "$column" || -z "$value" ]]; then
        echo "Error: Missing required arguments for 'Enrichment' mode."
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
            value_number=$(cut -d ";" -f $column_number metadata.csv | grep -c)

            if [[ $value_number -eq 0 ]]; then
                echo "Error: provide value '$value' do not exist in the column '$column' of metadata."
                exit 1
            elif [[ $value_number -le 10 ]]; then
                echo "Error: provide value '$value' have an n of '$value_number' and it's too low for enrichment analysis."
                exit 1
            fi
        fi
    fi
fi

check_threads "$threads"

# Check for software presence
check_required_software "$(basename -s .sh "$0" )"

# ==============================================================================
# Main Block
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running RobustGenePrediction Module with the following parameters:"
echo "  database: ${database}"
echo -e "  Mode: ${mode}\n"

if [[ "${mode}" == "All" ]]
then
    # ==============================================================================
    # Checking Working directory structure
    # ==============================================================================
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${Working_directory}' structure."
    check_directory_structure "$Working_directory"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"
    organize_working_directory "${Working_directory}"

    # ==============================================================================
    # Run earlgrey
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running earlGreyAnnotationOnly command..."
    process_earlgrey "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    # ==============================================================================
    # Organize results
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Organizing results..."
    organize_files "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished."

elif [[ "${mode}" == "Cluster" ]]
then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running in '${mode}' mode."
    check_clusters "${Working_directory}" "${clusters_file}"
    awk '{print $1}' ${Working_directory}/Clusters/SelectedClusters.txt | sed $'s/[^[:print:]\t]//g' | while read ClusterId
    do
        internal_dir="${Working_directory}/Clusters/${ClusterId}/"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
        check_directory_structure "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"
        organize_working_directory "${internal_dir}"

        # ==============================================================================
        # process earlgrey
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: processing braker results.."
        process_earlgrey "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

        # ==============================================================================
        # Running metaeuk for captain identification
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running metaeuk for specialized captain prediction.."
        organize_files "${internal_dir}"
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding. \n"
    done
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"
fi
