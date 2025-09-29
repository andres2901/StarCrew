#!/bin/bash

# Function to print help message

function print_help() {
   echo "Script to perform a quick and raw prediction of genes for starships. It determines the statistics of each element prediction, filter the elements based on a minimum gene content (8) and organized it based on family association from the metadata."
   echo
   echo "Syntax: SAT RobustGenePrediction [ -h ] -w <genome_file> [ -m <mode> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-p, --proteinDB: protein database fasta file (required)."
   echo "-h, --headers: headers of the protein database (required)."
   echo "-m, --mode: Define the data that will be use for the gene prediction. This can be perform for all the data or for each cluster (Available mode: Cluster, All) (Default = Cluster)."
   echo "-ms, --minSize: Minimum size of a Cluster to be include in the analyzis when running the 'Cluster' mode (Default = 5) [range: 5 - 10]"
   echo "-t, --threads: Number of threads for Braker (Default = 8)"
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
protein_path=""
mode="Cluster"
minimum_size="5"
threads="8"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -f|--fasta) 
	        shift
	        fasta_path="$1"
	        ;; 
        -p|--proteinDB)
            shift
            protein_path="$1"
            ;;
        -m|--mode)
            shift
            mode="$1"
            ;;
        -ms|--minSize)
            shift
            minimum_size="$1"
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

if [[ -z "$Working_directory" || -z "$protein_path" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# check if working directory exist
if [[ ! -f "$Working_directory" ]]; then
    echo "Error: folder '$Working_directory' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${Working_directory:0:1}" == "/" ]]; then
        Working_directory=$(realpath $Working_directory)
    fi
fi

# Check if protein file exists
if [[ ! -f "$protein_path" ]]; then
    echo "Error: file '$protein_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${protein_path:0:1}" == "/" ]]; then
        protein_path=$(realpath $protein_path)
    fi
fi

# Check if mode parameter is correct
if [[ "$mode" != "All" && "$mode" != "Cluster" ]]; then
    echo "Error: provided mode '$mode' is not accepted."
    print_help
    exit 1
fi

# check if python folder exist
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: folder '$auxiliary_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${auxiliary_path:0:1}" == "/" ]]; then
        auxiliary_path=$(realpath $auxiliary_path)
    fi
    # Check the presence of the specific scripts
    if [[ ! -f "${auxiliary_path}/merge.py" ]]; then
        echo "Error: file '${auxiliary_path}/merge.py' does not exist."
        exit 1
    fi

    if [[ ! -f "${auxiliary_path}/gff_filter.py" ]]; then
        echo "Error: file '${auxiliary_path}/gff_filter.py' does not exist."
        exit 1
    fi
fi

# Check minimum size is within allowed range
if [[ "$minimum_size" =~ ^[0-9]+$ ]]; then
    if (( $minimum_size < 5 || $minimum_size > 10 )); then
        echo "Error: '$minimum_size' minimum size is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$minimum_size' is not a positive integer."
    print_help
    exit 1
fi

# Check thread parameter
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
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

    local fasta_file="${base_dir}/sequences.fa"

    if [[ ! -f $fasta_file ]]; then
        echo "Error: sequence fasta file does not exist in '$base_dir'."
        exit 1
    else
        fasta_path=$(realpath $fasta_file)
        local input_size=$(seqkit stats ${fasta_path} | grep "FASTA" | awk '{print $4}')
        echo -e "Number of input elements: ${input_size}\n"
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

generate_database() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"

    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database..."
    metaeuk createdb  $fasta_path ${working_dir}/ContigsDB --dbtype 2 -v 0

    metaeuk createdb $protein_path ${working_dir}/ProteinDB --dbtype 1 -v 0

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk..."
    metaeuk predictexons ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --max-seqs 500 --use-all-table-starts 1 --start-sens 7.5 -v 0

    metaeuk unitesetstofasta ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/metaeukFinal

    awk '{print $6}' ${working_dir}/metaeukFinal.headersMap.tsv | awk -F '|' '{print $1}' | sort -u > ${working_dir}/selected_headers.txt

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating database for Braker..."

    # Return headers from database
    seqkit seq -n $protein_path > ${working_dir}/database_headers.txt

    # Select proteins in the database
    seqkit grep -n -f ${working_dir}/database_headers.txt $protein_path > ${working_dir}/Selected_database.fa

    rm ${working_dir}/metaeuk*
    rm ${working_dir}/ProteinDB*
}

run_braker() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local temp_dir="${working_dir}/temp/"
    local output="${working_dir}/braker_filter.gff"

    braker --genome ${base_dir}/Sequences.fa --softmasking_off --downsampling_lambda=0 --prot_seq ${working_dir}/Selected_database.fa --gff3 --fungus --alternatives-from-evidence=false --augustus_args "--genemodel=complete --noInFrameStop=true" --threads=8 --workingdir ${working_dir}/braker --useexisting
}



# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running RobustGenePrediction Module with the following parameters:"
echo "  Protein database: ${protein_path}"
echo -e "  Mode: ${mode}\n"

if [[ "${mode}" == "All" ]]
then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${Working_directory}' structure."
    check_directory_structure "$Working_directory"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Creating required subdirectories"
    mkdir -p ${Working_directory}/Workspace/RobustGenePrediction/
    mkdir -p ${Working_directory}/Workspace/RobustGenePrediction/temp/

    # ==============================================================================
    # Creating protein database for braker run
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running metaeuk prediction to generate specialized protein database..."
    generate_database "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    # ==============================================================================
    # Running braker
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running braker gene prediction.."
    run_braker "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."
    
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${out_directory}."

elif [[ "${mode}" == "Cluster" ]]
then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running in '${mode}' mode. Analyzing clusters with a minimum size of ${minimum_size}."
    awk -v min="$minimum_size" 'NR>1{if($2>=min){print $1}}' ${Working_directory}/Clusters/cluster_stats.txt | sed $'s/[^[:print:]\t]//g' | | while read ClusterId
    do
        internal_dir="${Working_directory}/Clusters/${ClusterId}/"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
        check_directory_structure "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."
        
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Creating required subdirectories"
        mkdir -p ${internal_dir}/Workspace/RobustGenePrediction/
        mkdir -p ${internal_dir}/Workspace/RobustGenePrediction/temp/

        # ==============================================================================
        # Creating protein database for braker run
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running metaeuk prediction to generate specialized protein database..."
        generate_database "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

        # ==============================================================================
        # Running braker
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running braker gene prediction.."
        run_braker "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."
    done
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"
fi

