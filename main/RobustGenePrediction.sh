#!/bin/bash

# Function to print help message

function print_help() {
   echo "Script to compact information from a gene prediction.
   This script perform five steps:
   1. Filter results from a previous gene prediction analysis of the elements.
   2. Identify putative captain genes through metaEuk.
   3. Merge and filter both previous predictions.
   4. Determines the statistics of each element prediction.
   - (All mode) Filter the elements based on a minimum gene content
   5. Organized the files in the working directory for future modules."
   echo
   echo "Syntax: SAT RobustGenePrediction [ -help ] -w <directory_path> [ -m <string> -mg <integer> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-m, --mode: Define the data that will be use for the gene prediction. This can be perform for all the data or for each cluster (Available mode: Cluster, All) (Default = Cluster)."
   echo "-mg, --minGene: Minimum number of genes in an element to be include in the dataset when running the 'All' mode (Default: 8) [range: 5 - 100]"
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
mode="Cluster"
minimum_gene_content="8"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
hmm_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../hmm/"
captain_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../databases/Captains.fa"
help_flag=false
metadata_flag=false

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
        -mg|--minGene)
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
if [[ ! -d "$Working_directory" ]]; then
    echo "Error: folder '$Working_directory' does not exist."
    exit 1
else
    # Check if the path is absolute
    if [[ ! "${Working_directory:0:1}" == "/" ]]; then
        Working_directory=$(realpath $Working_directory)
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
        echo -e "Number of input elements: ${input_size}\n"
    fi

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    else
        local Working_dir=$(find "$workspace_dir" -maxdepth 1 -type d -name "RobustGenePrediction" 2>/dev/null)
        if [[ -z "$Working_dir" ]]; then
            echo "Error: RobustGenePrediction directory not found in '$workspace_dir'." >&2
            echo "Run the module BrakerGenePrediction before this module." >&2
            exit 1
        else
            braker_file="${Working_dir}/PreliminarGenePrediction.gff"
            if [[ ! -f "$braker_file" ]]; then
                echo "Error: preliminary gene prediction not found in '$Working_dir'." >&2
                echo "Run the module BrakerGenePrediction before this module." >&2
                exit 1
            else
                braker_path=$(realpath $braker_file)
            fi

            local metadata_file="${Working_dir}/metadata.csv"
            if [[ -f $metadata_file ]]; then
                metadata_flag=true
            fi
        fi
    fi

    if [[ -z "$data_dir" ]]; then
        echo "Error: Data directory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -f $metadata_file ]]; then
        metadata_flag=true
    fi
}

process_braker() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local temp_dir="${working_dir}/temp/"
    local output="${working_dir}/braker_filter.gff"

    agat_sp_filter_incomplete_gene_coding_models.pl --gff ${braker_path} --fasta ${base_dir}/Sequences.fa -o ${temp_dir}/braker.gff &> /dev/null

    python ${auxiliary_path}/gff_filter.py -i ${temp_dir}/braker.gff -o ${temp_dir}/selectedgenes.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${temp_dir}/braker.gff --keep_list ${temp_dir}/selectedgenes.txt --output ${output} &> /dev/null

}

run_captain_metaeuk() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local temp_dir="${working_dir}/temp/"
    local hmmer_results="${working_dir}/metaeuk_hmmer.txt"
    local output="${working_dir}/metaeuk_filter.gff"

    # Create databases for metaeuk
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Creating metaeuk database."
    metaeuk createdb  $fasta_path ${working_dir}/ContigsDB --dbtype 2 -v 0

    metaeuk createdb $captain_path ${working_dir}/ProteinDB --dbtype 1 -v 0

    # Run metaeuk gene prediction
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Running metaeuk predictexons."
    metaeuk predictexons ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukResults ${working_dir}/tempFolder -s 7.5 --exhaustive-search-filter 1 --filter-msa 1 --chain-alignments 1  --remove-tmp-files 1 --use-all-table-starts 1 --start-sens 7.5 --orf-start-mode 0 -v 0 &> /dev/null

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing redundancy from metaeuk."
    metaeuk reduceredundancy ${working_dir}/metaeukResults ${working_dir}/metaeukpred ${working_dir}/metaeukgroups -v 0

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Generating gff file from metaeuk results."
    metaeuk unitesetstofasta ${working_dir}/ContigsDB ${working_dir}/ProteinDB ${working_dir}/metaeukpred ${working_dir}/metaeukFinal -v 0

    # Modify gff that agat can recognize
    sed -e 's/Target_ID=.*;TCS_//g' ${working_dir}/metaeukFinal.gff > ${working_dir}/metaeuk.gff

    #Remove elements with incomplete gene coding models
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing gene models without start and stop codon."

    agat_sp_filter_incomplete_gene_coding_models.pl --gff ${working_dir}/metaeuk.gff --fasta ${base_dir}/Sequences.fa -o ${working_dir}/metaeuk_fix.gff &> /dev/null

    agat_sp_extract_sequences.pl --gff ${working_dir}/metaeuk_fix.gff --fasta ${base_dir}/Sequences.fa -t exon --merge -p -o ${working_dir}/metaeuk_protein.fa &> /dev/null

    sed -i 's/_mRNA//g' ${working_dir}/metaeuk_protein.fa

    hmmsearch --max --noali --domE 10e-6 --domtblout ${hmmer_results} ${hmm_path}/Captain.hmm ${working_dir}/metaeuk_protein.fa &> /dev/null

    grep -v "#" ${hmmer_results} | awk '{print $1}' | sort -u > ${temp_dir}/selected_captains.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${working_dir}/metaeuk_fix.gff --keep_list ${temp_dir}/selected_captains.txt --output ${output} &> /dev/null

    # Cleanup temporary files and intermediate results
    rm ${working_dir}/metaeukResults*
    rm ${working_dir}/metaeukpred*
    rm ${working_dir}/metaeukgroups*
    rm ${working_dir}/metaeukFinal*
    rm ${working_dir}/ContigsDB*
    rm ${working_dir}/ProteinDB*
    rm ${working_dir}/*_incomplete.gff
}

merge_models() {

    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local braker_results="${working_dir}/braker_filter.gff"
    local metaeuk_results="${working_dir}/metaeuk_filter.gff"
    local temp_dir="${working_dir}/temp/"
    local out_gff="${working_dir}/Final_model.gff"

    if [[ ! -f "$braker_results" || ! -f "$metaeuk_results" ]]; then
        echo "Error: result from one of the gene predictors cannot be found." >&2
        exit 1
    fi

    agat_sp_keep_longest_isoform.pl --gff ${braker_results} -o ${temp_dir}/LongIso.gff &> /dev/null

    agat_sp_merge_annotations.pl --gff ${temp_dir}/LongIso.gff --gff ${metaeuk_results} --out ${temp_dir}/merge.gff &> /dev/null

    python ${auxiliary_path}/merge.py ${temp_dir}/merge.gff ${temp_dir}/modelsKeep.txt

    agat_sp_filter_feature_from_keep_list.pl --gff ${temp_dir}/merge.gff --keep_list ${temp_dir}/modelsKeep.txt --output ${temp_dir}/merge_keep.gff &> /dev/null

    agat_sp_keep_longest_isoform.pl --gff ${temp_dir}/merge_keep.gff -o ${out_gff} &> /dev/null
}

gene_stats() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local temp_directory="${base_dir}/Workspace/RobustGenePrediction/temp/"

    awk 'NR>1{print $1}' ${working_dir}/Final_model.gff | sort | uniq | while read line
    do 
        echo -e ${line}"\t"$(grep -w ${line} ${working_dir}/Final_model.gff | grep -c "gene")"\t"$(grep -w ${line} ${working_dir}/Final_model.gff | grep "gene" | awk '{ sum  += $5 - $4 } END { print sum / NR }') >> ${temp_directory}stats_gff3.txt 
    done

    awk 'NR>1{print $1}' ${working_dir}/Final_model.gff | sort | uniq | while read line
    do 
        grep -E "${line}|gene" ${working_dir}/Final_model.gff | sort -k4 -n | awk 'NR==1 {prev_col2 = $5; next} {diff = prev_col2 - $4; if (diff > 0) total_sum += diff; prev_col2 = $5} END {print total_sum / (NR - 1)}' >> ${temp_directory}intergenic.txt
    done

    echo -e "Starship""\t""Number_genes""\t""Avg_gene_length""\t""Avg_intergenic_length" > ${base_dir}/Gene_stats.txt

    paste ${temp_directory}stats_gff3.txt ${temp_directory}intergenic.txt >> ${base_dir}/Gene_stats.txt
}

organize_files() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/RobustGenePrediction/"
    local Output_dir="${base_dir}/Data/"
    local temp_directory="${base_dir}/Workspace/RobustGenePrediction/temp/"

    mkdir -p ${temp_directory}multiple ${temp_directory}gff/ ${temp_directory}protein/ ${temp_directory}exon/ ${Output_dir}/Protein/ ${Output_dir}/Gff/ ${Output_dir}/Nucleotide/ ${Output_dir}/Exon/

    # Divide the gff result for each element
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Dividing gff for individual elements."
    awk 'NR>1{print $1}' ${working_dir}/Final_model.gff | sort | uniq | while read line
    do 
        head -n1 ${working_dir}/Final_model.gff > ${temp_directory}multiple/${line}.gff
        grep -w $line ${working_dir}/Final_model.gff >> ${temp_directory}multiple/${line}.gff
    done

    awk 'NR>1{print $1}' ${working_dir}/Final_model.gff | sort | uniq | while read line
    do 
        agat_sp_manage_IDs.pl --gff ${temp_directory}multiple/${line}.gff --prefix ${line}. -o ${temp_directory}gff/${line}.gff &> /dev/null
    done

    find ./ -name "*.agat.log" -delete

    # Filter gff files to select elements above a threshold
    if [[ "${mode}" == "All" ]]
    then
        if $metadata_flag; then
            local metadata_file="${working_dir}/metadata.csv"
            local updated_metadata="${Output_dir}metadata.csv"
        fi
        
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
    elif [[ "${mode}" == "Cluster" ]]
    then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Storing related files for each element gene prediction..."
        awk 'NR>1{print $1}' ${base_dir}/Gene_stats.txt | sed $'s/[^[:print:]\t]//g' | while read line
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

            find ./ -name "*.agat.log" -delete
        done
    fi

    rm ${temp_directory}temp_element.txt
    rm ${base_dir}/*index*
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
if [[ "${mode}" == "All" ]]
then
    echo "  Mode: ${mode}"
    echo -e "  Minimum gene content: ${minimum_gene_content}\n"
elif [[ "${mode}" == "Cluster" ]]
then
    echo -e "  Mode: ${mode}\n"
fi

if [[ "${mode}" == "All" ]]
then
    # ==============================================================================
    # Checking Working directory structure
    # ==============================================================================
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${Working_directory}' structure."
    check_directory_structure "$Working_directory"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."

    # ==============================================================================
    # Process braker results
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: processing braker results..."
    process_braker "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    # ==============================================================================
    # Running metaeuk for captain identification
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running metaeuk for specialized captain prediction.."
    run_captain_metaeuk "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

    # ==============================================================================
    # merging models to select a putative captain
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: merging Braker and metaeuk results.."
    merge_models "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

    # ==============================================================================
    # Generate gene prediction statistics
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Creating gene statistics file.."
    gene_stats "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

    # ==============================================================================
    # Dividing and organizing results
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Organizing files.."
    organize_files "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${out_directory}."

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

        # ==============================================================================
        # process braker results
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: processing braker results.."
        process_braker "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

        # ==============================================================================
        # Running metaeuk for captain identification
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running metaeuk for specialized captain prediction.."
        run_captain_metaeuk "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

        # ==============================================================================
        # merging models to select a putative captain
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: merging Braker and metaeuk results.."
        merge_models "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

        # ==============================================================================
        # Generate gene prediction statistics
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Creating gene statistics file.."
        gene_stats "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

        # ==============================================================================
        # Dividing and organizing results
        # ==============================================================================

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Organizing files.."
        organize_files "${internal_dir}"
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished. Proceeding.\n"
    done
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"
fi
