#!/usr/bin/env bash

# Function to print help message
function print_help() {
   echo "Script to perform BRAKER gene prediction for use with the subsequent RobustGenePrediction command.
   It executes two main steps:
   1. Identifies potential proteins using MetaEuk gene prediction and creates a smaller, focused protein database from a user-customized input database.Identify proteins that can be produce by the elements trough metaeuk geneprediction and create a smaller protein database from a user-custom protein database.
   2. Executes BRAKER gene prediction using the newly created, smaller protein database."
   echo
   echo "Syntax: SAT RobustGenePrediction [ -help ] -s <file_path> -p <file_path> [ -m <string> -ms <integer> -t <integer> ]"
   echo "options:"
   echo "-s, --Sequence: Specify the path to the fasta file with ships (required)."
   echo "-p, --proteinDB: protein database fasta file (required)."
   echo "-t, --threads: Number of threads for Braker (Default = 8)"
   echo "-help: Display this help message."
}

# Initialize variables
fasta_path=""
protein_path=""
threads="8"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -s|--Sequence) 
	        shift
	        fasta_path="$1"
	        ;; 
        -p|--proteinDB)
            shift
            protein_path="$1"
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
if [[ -z "$fasta_path" || -z "$protein_path" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if working directory exist
if [[ ! -f "$fasta_path" ]]; then
    echo "Error: file '$fasta_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${fasta_path:0:1}" == "/" ]]; then
        fasta_path=$(realpath $fasta_path)
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

if [[ -z "$(which cd-hit)" ]]; then
    echo "Error: Missing metaeuk function."
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

generate_database() {
    local nucleotide="$1"
    local protein="$2"

    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database..."

    metaeuk createdb  $nucleotide ContigsDB --dbtype 2
    metaeuk createdb $protein ProteinDB --dbtype 1

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk..."

    metaeuk predictexons ContigsDB ProteinDB metaeukResults tempFolder -s 7.5 --orf-start-mode 0 --min-length 50 --min-seq-id 0.3 --cov-mode 1 -c 0.50 --remove-tmp-files 1 --max-seqs 1000 --use-all-table-starts 1 --min-aln-len 50
    metaeuk reduceredundancy metaeukResults metaeukpred metaeukgroups
    metaeuk unitesetstofasta ContigsDB ProteinDB metaeukpred metaeukFinal

    awk '{print $6}' metaeukFinal.headersMap.tsv | awk -F '|' '{print $1}' | sort -u > selected_headers.txt

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating database for Braker..."

    # Return headers from database
    seqkit seq -n $protein > database_headers.txt

    grep -f selected_headers.txt database_headers.txt > selected_full_headers.txt

    # Select proteins in the database
    seqkit grep --quiet -n -f selected_full_headers.txt $protein > Selected_database.fa

    sed -i -e 's/\///g' -e 's/(//g' -e 's/)//g' -e 's/|//g' -e 's/ //g' Selected_database.fa

    #Remove unnecesary files
    rm metaeuk*
    rm ProteinDB*
    rm ContigsDB*
    rm database_headers.txt
}

run_braker() {
    local nucleotide="$1"
    local protein="$2"

    braker --genome ${nucleotide} --softmasking_off --downsampling_lambda=0 --prot_seq ${protein} --gff3 --fungus --augustus_args "--genemodel=complete --noInFrameStop=true" --threads=${threads} --useexisting
}

# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running RobustGenePrediction Module with the following parameters:"


# ==============================================================================
# Creating protein database for braker run
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running metaeuk prediction to generate specialized protein database..."
generate_database "${fasta_path}" "${protein_path}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

# ==============================================================================
# Running braker
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running braker gene prediction.."
run_braker "${fasta_path}" "Selected_database.fa"
