#!/bin/bash

# Function to print help message

function print_help() {
   echo "Script to predict TEs in the sequences based on earlgrey approach and Mycomobilome database.
   This script perform five steps:
   1. Run earlgrey TE prediction.
   2. Organize the data.
   3. Merge and filter both previous predictions.
   4. Determines the statistics of each element prediction.
   - (All mode) Filter the elements based on a minimum gene content
   5. Organized the files in the working directory for future modules."
   echo
   echo "Syntax: SAT TEPrediction [ -help ] -w <directory_path> -d <file_path> [ -m <string> -t <integer> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-d, --database: path to Mycomobylome database (required)."
   echo "-m, --mode: Define the data that will be use for the TE prediction. This can be perform for all the data or for each cluster (Available mode: Cluster, All) (Default = Cluster)."
   echo "-t, --threads: Number of threads for earlgrey (Default = 8)"
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
mode="Cluster"
database=""
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
        -d|--database)
            shift
            database="$1"
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

# Check for mandatory arguments

if [[ -z "$Working_directory" || -z "$database" ]]; then
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

# Check if database file exists
if [[ ! -f "$database" ]]; then
    echo "Error: file '$database' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${database:0:1}" == "/" ]]; then
        database=$(realpath $database)
    fi
fi

# Check if mode parameter is correct
if [[ "$mode" != "All" && "$mode" != "Cluster" ]]; then
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

# Check for software presence
if [[ -z "$(which earlGreyAnnotationOnly)" ]]; then
    echo "Error: Missing earlGreyAnnotationOnly function."
    exit 1
fi

if [[ -z "$(which seqkit)" ]]; then
    echo "Error: Missing seqkit function."
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

check_directory_structure() {
    local base_dir="$1"
    
    # Locate required subdirectories and file
    local workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)
    local data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

    local fasta_file="${base_dir}/Sequences.fa"

    if [[ ! -f $fasta_file ]]; then
        echo "Error: sequence fasta file does not exist in '$base_dir'."
        exit 1
    else
        fasta_path=$(realpath $fasta_file)
        local input_size=$(seqkit stats ${fasta_path} | grep "FASTA" | awk '{print $4}')
        echo -e "  Number of input elements: ${input_size}\n"
    fi

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$data_dir" ]]; then
        echo "Error: Data directory not found in '$base_dir'." >&2
        exit 1
    fi
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local working_dir="${base_dir}/Workspace/TEPrediction/"

    # Create required subdirectories
    mkdir -p ${working_dir}

    # Copy sequence data making sure that is not softmasked
    seqkit seq -u $fasta_path > ${working_dir}/Sequences.fa
}

process_earlgrey() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/TEPrediction/"

    earlGreyAnnotationOnly -g "${working_dir}/Sequences.fa" -s TE -o "${working_dir}" -l ${database} -t 16 -m yes &> /dev/null
}

organize_files() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/TEPrediction/"
    local input_file="${working_dir}/TE_EarlGrey/TE_summaryFiles/TE.filteredRepeats.bed"
    local Output_dir="${base_dir}/Data/"

    if [[ ! -f $input_file ]]; then
        echo "Error: Earlgrey results are not found."
        exit 1
    fi

    mkdir -p ${Output_dir}/TEs/

    # Divide the bed result for each element
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing bed for individual elements."
    awk 'NR>1{print $1}' ${input_file} | sort | uniq | while read line
    do 
        grep -w ^"$line" ${input_file} > ${Output_dir}/TEs/${line}.bed
    done
}

check_clusters() {
    local base_dir="$1"

    local cluster_information="${base_dir}/Clusters/SelectedClusters.txt"

    if [[ ! -f $cluster_information ]]; then
        echo "Error: SelectedClusters.txt file does not exist in '${base_dir}/Clusters/'."
        exit 1
    else
        Cluster_number=$(wc -l "$cluster_information" | awk '{print $1}')
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing '${Cluster_number}' clusters."
    fi
}

# ==============================================================================
# Start the process
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
    check_clusters "${Working_directory}"
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
