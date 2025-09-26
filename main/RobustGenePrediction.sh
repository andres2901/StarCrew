#!/bin/bash

# Function to print help message

function print_help() {
   echo "Script to perform a quick and raw prediction of genes for starships. It determines the statistics of each element prediction, filter the elements based on a minimum gene content (8) and organized it based on family association from the metadata."
   echo
   echo "Syntax: SAT RobustGenePrediction [ -h ] -w <genome_file> [ -m <mode> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-m, --mode: Define the data that will be use for the gene prediction. This can be perform for the whole dataset or specifics clusters (Available mode: Cluster, Whole) (Default = Cluster)."
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
mode="Cluster"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
hmm_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../hmm/"
captain_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../databases/captains.fa"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
	    -f|--fasta) 
	        shift
	        fasta_path="$1"
	        ;; 
        -m|--mode)
            shift
            mode="$1"
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

if [[ -z "$Working_directory" ]]; then
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

# ==============================================================================
# Bash function block
# ==============================================================================

check_directory_structure() {
    local base_dir="$1"
    
    # Locate required subdirectories and file
    local workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)
    local data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

    local fasta_file="${base_dir}/sequences.fa"
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

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)

    local cluster_dir="${base_dir}/Clusters/"
    local working_dir="${base_dir}/Workspace/SyntenyClustering/"
    local temp_dir="${base_dir}/Workspace/SyntenyClustering/temp/"

    # Create required subdirectories
    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}
    mkdir -p ${cluster_dir}

    # Copy require files
    cp -r ${gff_dir} ${working_dir}
    cp -r ${nucleotide_dir} ${working_dir}
    cp -r ${protein_dir} ${working_dir}
    if $metadata_flag; then
        cp "${data_dir}/metadata.csv" ${working_dir}
    fi
}

generate_database() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/QuickGenePrediction/"

    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database..."
    metaeuk createdb  $fasta_path ${working_dir}/ContigsDB --dbtype 2 -v 0

    metaeuk createdb $protein_path ${working_dir}/ProteinDB --dbtype 1 -v 0

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk predictexons..."
    metaeuk easy-predict ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --max-seqs 500 --use-all-table-starts 1 --start-sens 7.5 -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing redundancy from metaeuk..."
    metaeuk reduceredundancy ${working_dir}/metaeukResults ${working_dir}/metaeukpred ${working_dir}/metaeukgroups -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating gff file from metaeuk results..."
    metaeuk unitesetstofasta ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukpred ${working_dir}/metaeukFinal -v 0

    # Modify gff, so Agat can recognized it structure
    sed -e 's/Target_ID=.*;TCS_//g' metaeukFinal.gff > metaeuk.gff

    # Remove elements with incomplete gene coding models
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing gene models without start and stop codon..."
    agat_sp_filter_incomplete_gene_coding_models.pl --gff metaeuk.gff --fasta $fasta_path -o ${working_dir}/metaeuk_fix.gff &> /dev/null

    find ./ -name "*.agat.log" -delete
}

run_braker() {
    local work_dir="$1"

    local temp_prefix="${work_dir}/temp_"
    local output="${work_dir}/braker_filter.gff"

    braker --genome <sequences> --softmasking_off --downsampling_lambda=0 --prot_seq <protein_database> --gff3 --fungus --alternatives-from-evidence=false --augustus_args "--genemodel=complete --noInFrameStop=true" --threads=8 --useexisting

    agat_sp_filter_incomplete_gene_coding_models.pl --gff "$gff_path" --fasta "$fasta_path" -o ${temp_prefix}braker.gff &> /dev/null

    python ${auxiliary_path}/gff_filter.py -i ${temp_prefix}braker.gff -o ${temp_prefix}selectedgenes.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${temp_prefix}braker.gff --keep_list ${temp_prefix}selectedgenes.txt --output ${output} &> /dev/null

    rm ${temp_prefix}*

}

run_metaeuk() {
    local work_dir="$1"
    local genome="$2"
    local protein="$3"


    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database."
    metaeuk createdb  $genome ContigsDB --dbtype 2 -v 0

    metaeuk createdb $protein ProteinDB --dbtype 1 -v 0

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk predictexons."
    metaeuk predictexons ContigsDB ProteinDB metaeukResults tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --use-all-table-starts 1 --start-sens 7.5 --orf-start-mode 0 -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing redundancy from metaeuk."
    metaeuk reduceredundancy metaeukResults metaeukpred metaeukgroups -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating gff file from metaeuk results."
    metaeuk unitesetstofasta ContigsDB ProteinDB metaeukpred metaeukFinal -v 0

    # Modify gff that agat recognized

    sed -e 's/Target_ID=.*;TCS_//g' metaeukFinal.gff > metaeuk.gff

    #Remove elements with incomplete gene coding models
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing gene models without start and stop codon."
    agat_sp_filter_incomplete_gene_coding_models.pl --gff metaeuk.gff --fasta $genome -o ${work_dir}/metaeuk_fix.gff &> /dev/null

    agat_sp_extract_sequences.pl --gff ${work_dir}/metaeuk_fix.gff --fasta $genome -t exon --merge -p -o ${work_dir}/metaeuk_protein.fa &> /dev/null

    sed -i 's/_mRNA//g' ${work_dir}/metaeuk_protein.fa

    # Cleanup temporary files and intermediate results
    rm metaeuk*
    rm -r temp*
    rm ContigsDB*
    rm ProteinDB*
    rm -r ${work_dir}/*_incomplete.gff

}

select_models() {
    local work_dir="$1"
    local hmm_profile="$2"

    local protein="${work_dir}/metaeuk_protein.fa"
    local gff_path="${work_dir}/metaeuk_fix.gff"
    local hmmer_results="${work_dir}/metaeuk_hmmer.txt"
    local temp_prefix="${work_dir}/temp_"
    local out_gff="${work_dir}/metaeuk_selected.gff"

    # Check presence of require files 

    if [[ ! -f "$protein" ]]; then
        echo "Error: metaeuk extracted protein not found in '$work_dir'." >&2
        exit 1
    fi

    if [[ ! -f "$gff_path" ]]; then
        echo "Error: metaeuk gff results not found in '$work_dir'." >&2
        exit 1
    fi

    hmmsearch --max --noali --domE 10e-6 --domtblout ${hmmer_results} ${hmm_profile} ${protein} &> /dev/null

    grep -v "#" ${hmmer_results} | awk '{print $1}' | sort -u > ${temp_prefix}selected_captains.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${gff_path} --keep_list ${temp_prefix}selected_captains.txt --output ${out_gff} &> /dev/null

    rm ${temp_prefix}*
    rm ${gff_path}
    rm ${protein}
    rm ${hmmer_results}
}

merge_models() {
    local work_dir="$1"
    local braker_gff="$2" 
    local temp_prefix="${work_dir}/temp_"
    local gff_path="${work_dir}/metaeuk_selected.gff"
    local out_gff="${work_dir}/Final_model.gff"


    if [[ ! -f "$gff_path" ]]; then
        echo "Error: metaeuk gff with selected results not found in '$work_dir'." >&2
        exit 1
    fi

    agat_sp_keep_longest_isoform.pl -g ${braker_gff} -o ${temp_prefix}LongIso.gff &> /dev/null

    agat_sp_merge_annotations.pl --gff ${temp_prefix}LongIso.gff --gff ${gff_path} --out ${temp_prefix}merge.gff &> /dev/null

    python3 ${auxiliary_path}/merge.py ${temp_prefix}merge.gff ${temp_prefix}modelsKeep.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${temp_prefix}merge.gff --keep_list ${temp_prefix}modelsKeep.txt --output ${out_gff} &> /dev/null

    rm ${temp_prefix}*

}

# ==============================================================================
# Run metaeuk gene prediction
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running metaeuk prediction only with captain proteins."

run_metaeuk "${out_directory}" "${fasta_path}" "${captain_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

# ==============================================================================
# Select gene models
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Selecting gene models based on hmm profile."

select_models "${out_directory}" "${hmm_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

# ==============================================================================
# Merge and generate final gene model
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Filtering braker results."

filter_braker "${out_directory}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Merging Braker results and metaeuk captain results."

merge_models "${out_directory}" "${out_directory}/braker_filter.gff"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

# ==============================================================================
# Generate gene prediction statistics
# ==============================================================================

gff_path="${out_directory}/Final_model.gff"

if [[ ! -f "$gff_path" ]]; then
    echo "Error: File $gff_path with the results of the merging process not found." >&2
    exit 1
fi

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Generating statistics of gene prediction."

# Generate statistics

mkdir temp_directory

awk 'NR>1{print $1}' $gff_path | sort -u | while read line
do 
  echo -e ${line}"\t"$(grep ${line} ${gff_path} | grep -c "gene")"\t"$(grep ${line} ${gff_path} | grep "gene" | awk '{ sum  += $5 - $4 } END { print sum / NR }') >> temp_directory/stats_gff3.txt 
done

awk 'NR>1{print $1}' ${gff_path} | sort -u | while read line
do 
  grep -E "${line}|gene" ${gff_path} | sort -k4 -n | awk 'NR==1 {prev_col2 = $5; next} {diff = prev_col2 - $4; if (diff > 0) total_sum += diff; prev_col2 = $5} END {print total_sum / (NR - 1)}' >> temp_directory/intergenic.txt
done

echo -e "Starship""\t""Number_genes""\t""Avg_gene_length""\t""Avg_intergenic_length" > ${out_directory}/gene_stats.txt

paste temp_directory/stats_gff3.txt temp_directory/intergenic.txt >> ${out_directory}/gene_stats.txt

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished. Proceeding."

# ==============================================================================
# Dividing and organizing results
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 6: Organizing results."

# Create directories
mkdir temp_directory/multiple ${out_directory}/updated_gff/ ${out_directory}/updated_protein/ ${out_directory}/updated_exon/ ${out_directory}/updated_nucleotide/ ${out_directory}/temp_gff/ ${out_directory}/temp_exon/ ${out_directory}/temp_protein/

# Divide the gff result for each element
echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing gff for individual elements."
awk 'NR>1{print $1}' ${gff_path} | sort -u | while read line
do 
  head -n1 ${gff_path} > temp_directory/multiple/${line}.gff
  grep -w $line ${gff_path} >> temp_directory/multiple/${line}.gff
done

awk 'NR>1{print $1}' ${gff_path} | sort -u | while read line
do 
  agat_sp_keep_longest_isoform.pl --gff temp_directory/multiple/${line}.gff -o ${out_directory}/temp_gff/${line}.gff &> /dev/null
done

awk 'NR>1{print $1}' ${gff_path} | sort -u | while read line
do 
  agat_sp_manage_IDs.pl --gff ${out_directory}/temp_gff/${line}.gff --prefix ${line}. -o ${out_directory}/updated_gff/${line}.gff &> /dev/null
done

# Cleanup temporary directory
rm -r temp*
find ./ -name "*.agat.log" -delete

# Generate proteome fasta for each element
echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generate proteome and exome fasta for each element."
ls ${out_directory}/updated_gff/ | awk -F '/' '{print $NF}' | sed 's/\.gff//g' | while read line
do 
  agat_sp_extract_sequences.pl --gff ${out_directory}/updated_gff/${line}.gff --fasta $fasta_path -t exon --merge -p -o ${out_directory}/temp_protein/${line}.fa &> /dev/null
  agat_sp_extract_sequences.pl --gff ${out_directory}/updated_gff/${line}.gff --fasta $fasta_path -t exon --merge -o ${out_directory}/temp_exon/${line}.fa &> /dev/null
  awk '{if($2){$1=">"$2} print $1}' ${out_directory}/temp_protein/${line}.fa | sed 's/gene=//g' > ${out_directory}/updated_protein/${line}.fa
  awk '{if($2){$1=">"$2} print $1}' ${out_directory}/temp_exon/${line}.fa | sed 's/gene=//g' > ${out_directory}/updated_exon/${line}.fa
done

# Generate nucleotide individual files of filter elements
echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing nucleotide fasta for individual elements."
awk -v min="$minimum_gene_content" 'NR>1{if($2>=min){print $1}}' ${out_directory}/gene_stats.txt | sed $'s/[^[:print:]\t]//g' | while read line
do
  echo $line > ${out_directory}/temp_element.txt
  seqkit grep -n -f ${out_directory}/temp_element.txt $fasta_path -o ${out_directory}/updated_nucleotide/${line}.fa &> /dev/null
done

find ./ -name "*.agat.log" -delete

# Cleanup temporary directories and files in output path
rm -r ${out_directory}/temp*

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 6 finished."
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${out_directory}."

