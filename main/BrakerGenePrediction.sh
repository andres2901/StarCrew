#!/bin/bash

# Function to print help message
function print_help() {
   echo "Script to perform BRAKER gene prediction for use with the subsequent RobustGenePrediction command.
   It executes two main steps:
   1. Identifies potential proteins using MetaEuk gene prediction and creates a smaller, focused protein database from a user-customized input database.Identify proteins that can be produce by the elements trough metaeuk geneprediction and create a smaller protein database from a user-custom protein database.
   2. Executes BRAKER gene prediction using the newly created, smaller protein database."
   echo
   echo "Syntax: SAT RobustGenePrediction [ -help ] -w <directory_path> -p <file_path> [ -m <string> -ms <integer> -t <integer> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-p, --proteinDB: protein database fasta file (required)."
   echo "-m, --mode: Define the data that will be use for the gene prediction. This can be perform for all the data or for each cluster (Default = Cluster) [Available mode: Cluster, All]."
   echo "-ms, --minSize: Minimum size of a Cluster to be include in the analyzis when running the 'Cluster' mode (Default = 4) [range: 4 - 10]"
   echo "-t, --threads: Number of threads for Braker (Default = 8)"
   echo "-help: Display this help message."
}

# Initialize variables
Working_directory=""
protein_path=""
mode="Cluster"
minimum_size="4"
threads="8"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
help_flag=false
metadata_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -w|--workingDirectory) 
	        shift
	        Working_directory="$1"
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

# Check if working directory exist
if [[ ! -d "$Working_directory" ]]; then
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

# Check if python folder exist
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
    if (( $minimum_size < 4 || $minimum_size > 10 )); then
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

# Check for required software
if [[ -z "$(which braker)" ]]; then
    echo "Error: Missing braker function."
    echo "Make sure that the function name is braker and no braker3, braker.pl or any other variation."
    exit 1
fi

if [[ -z "$(which seqkit)" ]]; then
    echo "Error: Missing seqkit function."
    exit 1
fi

if [[ -z "$(which metaeuk)" ]]; then
    echo "Error: Missing metaeuk function."
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

    # Check for the presence of required files/folders
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

    if [[ "$mode" == "All" ]]; then
        local metadata_file="${base_dir}/metadata_files/metadata.csv"
        if [[ -f $metadata_file ]]; then
            metadata_flag=true
        fi
    elif [[ "$mode" == "Cluster" ]]; then
        local metadata_file="${data_dir}/metadata.csv"
        if [[ -f $metadata_file ]]; then
            metadata_flag=true
        fi
    fi
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local temp_dir="${base_dir}/Workspace/RobustGenePrediction/temp/"

    # Create required subdirectories
    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    # Copy sequence data making sure that is not softmasked
    seqkit seq -u $fasta_path > ${temp_dir}/Sequences.fa
    mv ${temp_dir}/Sequences.fa $fasta_path

    if $metadata_flag; then
        if [[ "$mode" == "All" ]]; then
            local metadata_file="${base_dir}/metadata_files/metadata.csv"
        elif [[ "$mode" == "Cluster" ]]; then
            local metadata_file="${data_dir}/metadata.csv"
        fi
        cp "${metadata_file}" ${working_dir}
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

    metaeuk predictexons ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --max-seqs 500 --use-all-table-starts 1 --start-sens 7.5 -v 0 &> /dev/null
    metaeuk unitesetstofasta ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/metaeukFinal -v 0

    awk '{print $6}' ${working_dir}/metaeukFinal.headersMap.tsv | awk -F '|' '{print $1}' | sort -u > ${working_dir}/selected_headers.txt

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating database for Braker..."

    # Return headers from database
    seqkit seq -n $protein_path > ${working_dir}/database_headers.txt

    grep -f ${working_dir}/selected_headers.txt ${working_dir}/database_headers.txt > ${working_dir}/selected_full_headers.txt

    # Select proteins in the database
    seqkit grep --quiet -n -f ${working_dir}/selected_full_headers.txt $protein_path > ${working_dir}/Selected_database.fa

    sed -i -e 's/\///g' -e 's/(//g' -e 's/)//g' -e 's/|//g' -e 's/ //g' ${working_dir}/Selected_database.fa

    #Remove unnecesary files
    rm ${working_dir}/metaeuk*
    rm ${working_dir}/ProteinDB*
    rm ${working_dir}/ContigsDB*
    rm ${working_dir}/database_headers.txt
}

run_braker() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"

    braker --genome ${fasta_path} --softmasking_off --downsampling_lambda=0 --prot_seq ${working_dir}/Selected_database.fa --gff3 --fungus --alternatives-from-evidence=false --augustus_args "--genemodel=complete --noInFrameStop=true" --threads=8 --workingdir ${working_dir}/braker --useexisting &> /dev/null
}

run_braker_second() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"

    braker --genome ${fasta_path} --softmasking_off --downsampling_lambda=0 --prot_seq ${working_dir}/Selected_database.fa --gff3 --fungus --threads=8 --workingdir ${working_dir}/braker --useexisting &> /dev/null
}

check_braker() {
    local base_dir="$1"

    local braker_file="${base_dir}/Workspace/RobustGenePrediction/braker/braker.gff3"
    local output="${base_dir}/Workspace/RobustGenePrediction/PreliminarGenePrediction.gff"

    if [[ -f $braker_file ]]; then
        braker_flag=false
        cp ${braker_file} ${output}
    else
        braker_flag=true 
    fi
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

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"

    organize_working_directory "${Working_directory}"

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
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
    check_braker "${Working_directory}"
    if $braker_flag; then
        echo -e "\033[01;31mWARNING\033[m: Braker fails, this cluster needs to be manually checked."
        rm -r "${Working_directory}/Workspace/RobustGenePrediction/braker/"
    else
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  Successfull run of Braker."
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."
    
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${out_directory}."

elif [[ "${mode}" == "Cluster" ]]
then
    Cluster_number=$(awk -v min="$minimum_size" 'NR>1{if($2>=min){print $1}}' ${Working_directory}/Clusters/cluster_stats.txt | wc -l)
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running in '${mode}' mode. Analyzing '${Cluster_number}' clusters with a minimum size of ${minimum_size}."
    awk -v min="$minimum_size" 'NR>1{if($2>=min){print $1}}' ${Working_directory}/Clusters/cluster_stats.txt | sed $'s/[^[:print:]\t]//g' | while read ClusterId
    do
        internal_dir="${Working_directory}/Clusters/${ClusterId}/"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
        check_directory_structure "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."
        
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
        organize_working_directory "${internal_dir}"

        # ==============================================================================
        # Creating protein database for braker run
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running metaeuk prediction to generate specialized protein database..."
        generate_database "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

        # ==============================================================================
        # Running braker
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running braker gene prediction..."
        run_braker "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
        check_braker "${internal_dir}"
        if $braker_flag; then
            echo -e "\033[01;31mWARNING\033[m: Braker fails, trying a second time..."
            rm -r "${internal_dir}/Workspace/RobustGenePrediction/braker/"
            run_braker_second "${internal_dir}"
            check_braker "${internal_dir}"
            if $braker_flag; then
                echo -e "\033[01;31mWARNING\033[m: Braker fails for a second time, this cluster needs to be manually checked."
                rm -r "${internal_dir}/Workspace/RobustGenePrediction/braker/"
            else
                echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  Successfull run of Braker. Storing this cluster for further analysis."
                grep -w ${ClusterId} ${Working_directory}/Clusters/cluster_stats.txt >> ${Working_directory}/Clusters/SelectedClusters.txt
            fi
        else
            echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  Successfull run of Braker. Storing this cluster for further analysis."
            grep -w ${ClusterId} ${Working_directory}/Clusters/cluster_stats.txt >> ${Working_directory}/Clusters/SelectedClusters.txt
        fi
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding.\n"
    done
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"
fi
