#!/bin/bash

# Function to print help message

function print_help() {
   echo "Script to perform a quick and raw prediction of genes for starships. It determines the statistics of each element prediction, filter the elements based on a minimum gene content (8) and organized it based on family association from the metadata."
   echo
   echo "Syntax: SAT QuickGenePrediction [ -help ] -w <working_directory> -p <protein_db> [ -m <gene_number> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-p, --proteinDB: protein database fasta file (required)."
   echo "-m, --minGene: Minimum number of genes in an element to be include in the dataset (Default: 8) [range: 5 - 100]"
   echo "-help: Display this help message."
}

# Initialize variables
protein_path=""
Working_directory=""
minimum_gene_content="8"
metadata_flag=false
help_flag=false

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
        -m|--minGene)
            shift
            minimum_gene_content="$1"
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
if [[ -z "$protein_path" || -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
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

# Check minimum gene content is within allowed range
if [[ "$minimum_gene_content" =~ ^[0-9]+$ ]]; then
    if (( $minimum_gene_content < 5 || $minimum_gene_content > 100 )); then
        echo "Error: '$minimum_gene_content' minimum gene content is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$minimum_gene_content' is not a positive integer."
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

    local fasta_file="${base_dir}/Sequences.fa"
    local metadata_file="${base_dir}/metadata_files/metadata.csv"

    if [[ ! -f $fasta_file ]]; then
        echo "Error: sequence fasta file does not exist in '$base_dir'."
        exit 1
    else
        fasta_path=$(realpath $fasta_file)
        local input_size=$(seqkit stats ${fasta_path} | grep "FASTA" | awk '{print $4}')
        echo -e "Number of input elements: ${input_size}\n"
    fi

    if [[ -f $metadata_file ]]; then
        metadata_flag=true
    fi

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$data_dir" ]]; then
        echo "Error: Organized_data directory not found in '$base_dir'." >&2
        exit 1
    fi
}

gene_prediction() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/QuickGenePrediction/"

    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database..."
    metaeuk createdb  $fasta_path ${working_dir}/ContigsDB --dbtype 2 -v 0

    metaeuk createdb $protein_path ${working_dir}/ProteinDB --dbtype 1 -v 0

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk predictexons..."
    metaeuk predictexons ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --max-seqs 500 --use-all-table-starts 1 --start-sens 7.5 -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing redundancy from metaeuk..."
    metaeuk reduceredundancy ${working_dir}/metaeukResults ${working_dir}/metaeukpred ${working_dir}/metaeukgroups -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating gff file from metaeuk results..."
    metaeuk unitesetstofasta ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukpred ${working_dir}/metaeukFinal -v 0

    # Modify gff, so Agat can recognized it structure
    sed -e 's/Target_ID=.*;TCS_//g' ${working_dir}/metaeukFinal.gff > metaeuk.gff

    # Remove elements with incomplete gene coding models
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing gene models without start and stop codon..."
    agat_sp_filter_incomplete_gene_coding_models.pl --gff metaeuk.gff --fasta $fasta_path -o ${working_dir}/metaeuk_fix.gff &> /dev/null

    find ./ -name "*.agat.log" -delete
}

gene_stats() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/QuickGenePrediction/"
    local temp_directory="${base_dir}/Workspace/QuickGenePrediction/temp/"

    awk 'NR>1{print $1}' ${working_dir}/metaeuk_fix.gff | sort | uniq | while read line
    do 
        echo -e ${line}"\t"$(grep ${line} ${working_dir}/metaeuk_fix.gff | grep -c "gene")"\t"$(grep ${line} ${working_dir}/metaeuk_fix.gff | grep "gene" | awk '{ sum  += $5 - $4 } END { print sum / NR }') >> ${temp_directory}stats_gff3.txt 
    done

    awk 'NR>1{print $1}' ${working_dir}/metaeuk_fix.gff | sort | uniq | while read line
    do 
        grep -E "${line}|gene" ${working_dir}/metaeuk_fix.gff | sort -k4 -n | awk 'NR==1 {prev_col2 = $5; next} {diff = prev_col2 - $4; if (diff > 0) total_sum += diff; prev_col2 = $5} END {print total_sum / (NR - 1)}' >> ${temp_directory}intergenic.txt
    done

    echo -e "Starship""\t""Number_genes""\t""Avg_gene_length""\t""Avg_intergenic_length" > ${base_dir}/Gene_stats.txt

    paste ${temp_directory}stats_gff3.txt ${temp_directory}intergenic.txt >> ${base_dir}/Gene_stats.txt
}

organize_files() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/QuickGenePrediction/"
    local Output_dir="${base_dir}/Data/"
    local temp_directory="${base_dir}/Workspace/QuickGenePrediction/temp/"
    local fasta_sequence=""
    if $metadata_flag; then
        local metadata_file="${base_dir}/metadata_files/metadata.csv"
        local updated_metadata="${Output_dir}metadata.csv"
    fi

    mkdir -p ${temp_directory}multiple ${temp_directory}gff/ ${temp_directory}protein/ ${temp_directory}exon/ ${Output_dir}/Protein/ ${Output_dir}/Gff/ ${Output_dir}/Nucleotide/ ${Output_dir}/Exon/

    # Divide the gff result for each element
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing gff for individual elements."
    awk 'NR>1{print $1}' ${working_dir}/metaeuk_fix.gff | sort | uniq | while read line
    do 
        head -n1 ${working_dir}/metaeuk_fix.gff > ${temp_directory}multiple/${line}.gff
        grep $line ${working_dir}/metaeuk_fix.gff >> ${temp_directory}multiple/${line}.gff
    done

    awk 'NR>1{print $1}' ${working_dir}/metaeuk_fix.gff | sort | uniq | while read line
    do 
        agat_sp_manage_IDs.pl --gff ${temp_directory}multiple/${line}.gff --prefix ${line}. -o ${temp_directory}gff/${line}.gff &> /dev/null
    done

    find ./ -name "*.agat.log" -delete

    # Filter gff files to select elements above a threshold
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Filtering elements based on a minimum of '$minimum_gene_content' predicted genes and storing related files..."
    awk -v min="$minimum_gene_content" 'NR>1{if($2>=min){print $1}}' ${base_dir}/Gene_stats.txt | sed $'s/[^[:print:]\t]//g' | while read line
    do
        # Filtering gff files
        cp ${temp_directory}gff/${line}.gff ${Output_dir}/Gff/
        
        # Dividing nucleotide of elements
        echo $line > ${temp_directory}/temp_element.txt
        seqkit grep -n -f ${temp_directory}/temp_element.txt $fasta_path -o ${Output_dir}/Nucleotide/${line}.fa &> /dev/null

        # Creating Exome
        agat_sp_extract_sequences.pl --gff ${Output_dir}/Gff/${line}.gff --fasta $fasta_path -t exon --merge -o ${temp_directory}exon/${line}.fa &> /dev/null
        awk '{if($2){$1=">"$2} print $1}' ${temp_directory}exon/${line}.fa | sed 's/gene=//g' > ${Output_dir}/Exon/${line}.fa

        # Creating proteome
        seqkit translate ${Output_dir}/Exon/${line}.fa --trim > ${Output_dir}/Protein/${line}.fa

        # Updating metadata
        if $metadata_flag; then
            grep -w $line $metadata_file >> $updated_metadata
        fi

        find ./ -name "*.agat.log" -delete
    done

    local element_number=$(ls ${Output_dir}/Gff/ | wc -l)

    echo -e "  Elements that pass the filter stage: ${element_number}"

    rm ${base_dir}/*index*
    rm ${temp_directory}temp_element.txt
}

# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running QuickGenePrediction Module with the following parameters:"
echo "  Protein database: ${protein_path}"
echo -e "  Minimum gene content for an element to be mantain: ${minimum_gene_content}\n"

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${Working_directory}' structure."
check_directory_structure "$Working_directory"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Creating required subdirectories"

mkdir -p ${Working_directory}/Workspace/QuickGenePrediction/
mkdir -p ${Working_directory}/Workspace/QuickGenePrediction/temp/

# ==============================================================================
# Generate gene prediction
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running metaeuk for gene prediction."

gene_prediction "$Working_directory"

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> Gene prediction finished. Proceeding."

# ==============================================================================
# Generate gene prediction statistics
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Generating statistics of gene prediction."

gene_stats "$Working_directory"

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> Gene prediction statistics finished. Proceeding."

# ==============================================================================
# Dividing and organizing results
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Organizing and filtering results."

organize_files "$Working_directory"

# Clean temporary directory from workspace
rm -r ${Working_directory}/Workspace/QuickGenePrediction/temp/

echo "  [$(date "+%Y-%m-%d %H:%M:%S")] -> Organization and filtering of the results is finished."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${Working_directory}."

