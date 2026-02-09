#!/usr/bin/env bash

# Function to print help message

function print_help() {
   echo "Script to compact information from a gene prediction.
   This script perform five steps:
   1. Filter results from a previous gene prediction analysis of the elements.
   2. Identify putative captain genes through metaEuk.
   3. Merge and filter both previous predictions.
   4. Determines the statistics of each element prediction.
   5. Organized the files in the working directory for future modules."
   echo
   echo "Syntax: bash RobustGenePrediction.sh [ -help ] -w <directory_path> [ -m <string> -mg <integer> ]"
   echo "options:"
   echo "-s, --Sequence: Specify the path to the fasta file with ships (required)."
   echo "-g; --gff: Specify the path to the gff file from braker (required)"
   echo "-a; --auxiliary: Specify the path where the 'merge.py' and 'gff_filter.py' scripts are stored (required)."
   echo "-h; --hmm: Specify the path to the 'Captain' hmm profile (required)."
   echo "-c; --captain: Specify the path to Captain protein database (required)."
   echo "-help: Display this help message."
}

# Initialize variables

fasta_path=""
gff_path=""
auxiliary_path=""
hmm_path=""
captain_path=""
help_flag=false
metadata_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -s|--Sequence) 
            shift
            fasta_path="$1"
            ;; 
        -g|--gff) 
            shift
            gff_path="$1"
            ;; 
        -a|--auxiliary)
            shift
            auxiliary_path="$1"
            ;;
        -h|--hmm)
            shift
            hmm_path="$1"
            ;;
        -c|--captain)
            shift
            captain_path="$1"
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

if [[ -z "$fasta_path" || -z "$auxiliary_path" || -z "$hmm_path" || -z "$captain_path" || -z "$gff_path" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# check if working directory exist
if [[ ! -f "$fasta_path" ]]; then
    echo "Error: file '$fasta_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${fasta_path:0:1}" == "/" ]]; then
        fasta_path=$(realpath $fasta_path)
    fi
fi

if [[ ! -f "$gff_path" ]]; then
    echo "Error: file '$gff_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${gff_path:0:1}" == "/" ]]; then
        gff_path=$(realpath $gff_path)
    fi
fi

if [[ ! -f "$captain_path" ]]; then
    echo "Error: file '$captain_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${captain_path:0:1}" == "/" ]]; then
        captain_path=$(realpath $captain_path)
    fi
fi

if [[ ! -f "$hmm_path" ]]; then
    echo "Error: file '$hmm_path' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${hmm_path:0:1}" == "/" ]]; then
        hmm_path=$(realpath $hmm_path)
    fi
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

# Check for software presence
if [[ -z "$(which python)" ]]; then
    echo "Error: Missing python function."
    exit 1
fi

if [[ -z "$(which hmmsearch)" ]]; then
    echo "Error: Missing hmmsearch function."
    exit 1
fi

if [[ -z "$(which agat_sp_filter_incomplete_gene_coding_models.pl)" ]]; then
    echo "Error: Missing agat functions."
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

process_braker() {
    local temp_dir="temp/"
    local output="braker_filter.gff"

    agat_sp_filter_incomplete_gene_coding_models.pl --gff ${gff_path} --fasta ${fasta_path} -o ${temp_dir}/braker.gff

    python ${auxiliary_path}/gff_filter.py -i ${temp_dir}/braker.gff -o ${temp_dir}/selectedgenes.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${temp_dir}/braker.gff --keep_list ${temp_dir}/selectedgenes.txt --output ${output}
}

run_captain_metaeuk() {
    local temp_dir="temp/"
    local hmmer_results="metaeuk_hmmer.txt"
    local output="metaeuk_filter.gff"

    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database."
    metaeuk createdb $fasta_path ContigsDB --dbtype 2 -v 0

    metaeuk createdb $captain_path ProteinDB --dbtype 1 -v 0

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk predictexons."
    metaeuk predictexons ContigsDB ProteinDB metaeukResults tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --use-all-table-starts 1 --start-sens 7.5 --orf-start-mode 0 -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing redundancy from metaeuk."
    metaeuk reduceredundancy metaeukResults metaeukpred metaeukgroups -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating gff file from metaeuk results."
    metaeuk unitesetstofasta ContigsDB ProteinDB metaeukpred metaeukFinal -v 0

    # Modify gff that agat can recognize
    sed -e 's/Target_ID=.*;TCS_//g' metaeukFinal.gff > metaeuk.gff

    #Remove elements with incomplete gene coding models
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing gene models without start and stop codon."

    agat_sp_filter_incomplete_gene_coding_models.pl --gff metaeuk.gff --fasta ${fasta_path} -o metaeuk_fix.gff

    agat_sp_manage_IDs.pl --gff metaeuk_fix.gff --prefix metaeuk. -o metaeuk_fix2.gff

    agat_sp_extract_sequences.pl --gff metaeuk_fix2.gff --fasta $fasta_path -t exon --merge -p -o metaeuk_protein.fa

    sed -i 's/_mRNA//g' metaeuk_protein.fa

    hmmsearch --max --noali --domE 1e-3 --domtblout ${hmmer_results} ${hmm_path} metaeuk_protein.fa

    grep -v "#" ${hmmer_results} | awk '{print $1}' | sort -u > ${temp_dir}/selected_captains.txt

    agat_sp_filter_feature_from_keep_list.pl --gff metaeuk_fix2.gff --keep_list ${temp_dir}/selected_captains.txt --output ${output}

    # Cleanup temporary files and intermediate results
    rm metaeukResults*
    rm metaeukpred*
    rm metaeukgroups*
    rm metaeukFinal*
    rm ContigsDB*
    rm ProteinDB*
    rm *_incomplete.gff
}

merge_models() {
    local braker_results="braker_filter.gff"
    local metaeuk_results="metaeuk_filter.gff"
    local temp_dir="temp/"
    local out_gff="Final_model.gff"

    if [[ ! -f "$braker_results" || ! -f "$metaeuk_results" ]]; then
        echo "Error: result from one of the gene predictors cannot be found." >&2
        exit 1
    fi

    agat_sp_keep_longest_isoform.pl --gff ${braker_results} -o ${temp_dir}/LongIso.gff

    agat_sp_merge_annotations.pl --gff ${temp_dir}/LongIso.gff --gff ${metaeuk_results} --out ${temp_dir}/merge.gff

    python ${auxiliary_path}/merge.py ${temp_dir}/merge.gff ${temp_dir}/modelsKeep.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${temp_dir}/merge.gff --keep_list ${temp_dir}/modelsKeep.txt --output ${temp_dir}/merge_keep.gff

    agat_sp_keep_longest_isoform.pl --gff ${temp_dir}/merge_keep.gff -o ${out_gff}
}

gene_stats() {
    local temp_directory="temp/"

    awk 'NR>1{print $1}' Final_model.gff | sort | uniq | while read line
    do 
        echo -e ${line}"\t"$(grep -w ${line} Final_model.gff | grep -c "gene")"\t"$(grep -w ${line} Final_model.gff | grep "gene" | awk '{ sum  += $5 - $4 } END { print sum / NR }') >> ${temp_directory}stats_gff3.txt 
    done

    awk 'NR>1{print $1}' Final_model.gff | sort | uniq | while read line
    do 
        grep -E "${line}|gene" Final_model.gff | sort -k4 -n | awk 'NR==1 {prev_col2 = $5; next} {diff = prev_col2 - $4; if (diff > 0) total_sum += diff; prev_col2 = $5} END {print total_sum / (NR - 1)}' >> ${temp_directory}intergenic.txt
    done

    echo -e "Starship""\t""Number_genes""\t""Avg_gene_length""\t""Avg_intergenic_length" > Gene_stats.txt

    paste ${temp_directory}stats_gff3.txt ${temp_directory}intergenic.txt >> Gene_stats.txt
}


# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running RobustGenePrediction Module with the following parameters:"

mkdir temp/
# ==============================================================================
# Process braker results
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: processing braker results..."
process_braker
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

# ==============================================================================
# Running metaeuk for captain identification
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running metaeuk for specialized captain prediction.."
run_captain_metaeuk
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

# ==============================================================================
# merging models to select a putative captain
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: merging Braker and metaeuk results.."
merge_models
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

# ==============================================================================
# Generate gene prediction statistics
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Creating gene statistics file.."
gene_stats
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished."
