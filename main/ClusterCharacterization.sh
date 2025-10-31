#!/bin/bash

# Function to print help message

function print_help() {
   echo -e "Script to run the characterization of each selected cluster.
   This script perform eight steps:
   1. Identify orthogroups through OrthoFinder software.
   2. Perform all-vs-all Blastn.
   3. Perform a hierarchical clustering of the elements based on Orthogroup gene count including singletons.
   4. Determine the full conection of the cluster and create a cargo orthogroups heatmap and synteny image for the cluster.
   5. Identify possible individual nesting events inside the cluster.
   6. Identify core genes in the cluster in two ways:
     6.1. General core: orthogroups that are present in at least 80% of the elements in the cluster.
     6.2. Specific core:
       6.2.1. Divide the Cluster in subclusters of a height above 0.8 in the hierarchical clustering.
       6.2.2. If subslusters are generated identify core genes in each one that have at least 5 elements using the same logic of general core.
   7. If subclusters are present it try to identify putative cargo movement events including the specific orthogroups involve.
   8. Determine if there are discordances at 'Clade' lavel between CArgo hierarchical clustering and Captain phylogenetic tree.
   "
   echo
   echo "Syntax: SAT ClusterCharacterization [ -help ] -w <directory_path> [ -t <integer> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-t, --threads: Number of threads for orthofinder and blast (Default: 8)"
   echo "-help: Display this help message."
}

# Initialize variables

Working_directory=""
mode="Raw"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
threads="8"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workingDirectory)
            shift
            Working_directory="$1"
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

# Check for mandatory argument and define the path as absolute
if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if working directory exists
if [[ ! -d "$Working_directory" ]]; then
    echo "Error: Directory '$Working_directory' does not exist."
    exit 1
else
    Working_directory=$(realpath $Working_directory)
fi

# Check if directory where the R function are stored exists
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: Directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
fi

if [[ ! -f "${auxiliary_path}/ClusterAnalysis.R" ]]; then
    echo "Error: File '${auxiliary_path}/ClusterAnalysis.R' does not exist."
    exit 1
fi

# Check thread parameter
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
    exit 1
fi

# Check for software presence
if [[ -z "$(which Rscript)" ]]; then
    echo "Error: Missing Rscript function."
    exit 1
fi

if [[ -z "$(which blastn)" ]]; then
    echo "Error: Missing blastn function."
    exit 1
fi

if [[ -z "$(which makeblastdb)" ]]; then
    echo "Error: Missing makeblastdb function."
    exit 1
fi

if [[ -z "$(which orthofinder)" ]]; then
    echo "Error: Missing orthofinder function."
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

check_clusters() {
    local base_dir="$1"

    local cluster_information="${base_dir}/Clusters/ClustersAnalyzed.txt"

    if [[ ! -f $cluster_information ]]; then
        echo "Error: ClustersAnalyzed.txt file does not exist in '${base_dir}/Clusters/'."
        exit 1
    else
        Cluster_number=$(wc -l "$cluster_information" | awk '{print $1}')
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing '${Cluster_number}' clusters.\n"
    fi
}

check_directory_structure() {
    local base_dir="$1"
    
    # Locate required subdirectories and file
    local workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)

    if [[ -z "$workspace_dir" ]]; then
        echo "Error: Workspace directory not found in '$base_dir'." >&2
        exit 1
    fi

    local data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

    if [[ -z "$data_dir" ]]; then
        echo "Error: Data directory not found in '$base_dir'." >&2
        exit 1
    fi
    
    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    if [[ -z "$gff_dir" ]]; then
        echo "Error: GFF subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    if [[ -z "$protein_dir" ]]; then
        echo "Error: Protein subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    if [[ -z "$nucleotide_dir" ]]; then
        echo "Error: nucleotide subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    if [[ -z "$exon_dir" ]]; then
        echo "Error: exon subdirectory not found in '$data_dir'." >&2
        exit 1
    fi

    # Step 2: Check for consistent filenames across subdirectories

    # Define temporary file paths in the working directory
    local gff_files="$base_dir/gff_files.txt"
    local protein_files="$base_dir/protein_files.txt"
    local nucleotide_files="$base_dir/nucleotide_files.txt"
    local exon_files="$base_dir/exon_files.txt"
    
    # Get sorted list of base filenames from the GFF directory
    find "$gff_dir" -maxdepth 1 -type f -name "*.gff" | xargs -n 1 basename -s .gff | sort > "$gff_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$protein_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$protein_files"

    # Get sorted list of base filenames from the GFF directory
    find "$nucleotide_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$nucleotide_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$exon_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$exon_files"

    # Compare the lists. If diff finds a difference, it returns a non-zero exit code.
    if ! diff -q "$gff_files" "$protein_files" >/dev/null || \
        ! diff -q "$gff_files" "$nucleotide_files" >/dev/null || \
        ! diff -q "$gff_files" "$exon_files" >/dev/null || \
        ! diff -q "$protein_files" "$nucleotide_files" >/dev/null || \
        ! diff -q "$protein_files" "$exon_files" >/dev/null || \
        ! diff -q "$nucleotide_files" "$exon_files" >/dev/null; then
        echo "Error: File lists in subdirectories do not match." >&2
        echo "Details:" >&2
        echo "GFF vs. Protein:" >&2
        diff "$gff_files" "$protein_files" >&2
        echo "GFF vs. nucleotide:" >&2
        diff "$gff_files" "$nucleotide_files" >&2
        echo "GFF vs. exon:" >&2
        diff "$gff_files" "$exon_files" >&2
        echo "protein vs. nucleotide:" >&2
        diff "$protein_files" "$nucleotide_files" >&2
        echo "protein vs. exon:" >&2
        diff "$protein_files" "$exon_files" >&2
        echo "nucleotide vs. exon:" >&2
        diff "$nucleotide_files" "$exon_files" >&2
        rm "$gff_files" "$nucleotide_files" "$protein_files" "$exon_files"
        exit 1
    fi

    # Cleanup temporary files
    rm "$gff_files" "$nucleotide_files" "$protein_files" "$exon_files"   

    local captainPhylogeny="${base_dir}/CaptainPhylogeny.nw"

    if [[ ! -f $captainPhylogeny ]]; then
        echo "Error: captain phylogeny file does not exist in '$base_dir'."
        exit 1
    fi

    local captainremoval_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Captainless_elements" 2>/dev/null)

    if [[ -z "$captainremoval_dir" ]]; then
        captainremoval_number="0"
    else
        captainremoval_number=$(find "$captainremoval_dir" -maxdepth 1 -type f -name "*.gff" | sort -u | wc -l)
    fi
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    local working_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local temp_dir="${base_dir}/Workspace/ClusterCharacterization/temp/"

    local full_gff="${working_dir}/Final_model.gff"

    local captainPhylogeny="${base_dir}/CaptainPhylogeny.nw"

    # Create required subdirectories
    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    # Copy require files
    cp -r ${gff_dir} ${working_dir}
    cp -r ${nucleotide_dir} ${working_dir}
    cp -r ${protein_dir} ${working_dir}
    cp -r ${exon_dir} ${working_dir}
    cp ${captainPhylogeny} ${working_dir}/

    echo "##gff-version 3" > ${full_gff}

    cat ${gff_dir}/* | grep -v "#" >> ${full_gff}
}

run_orthofinder() {
    local base_dir="$1"

    local captainPhylogeny="${base_dir}/CaptainPhylogeny.nw"

    local working_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local output_dir="${working_dir}/Orthofinder"

    orthofinder -t "${threads}" -f ${protein_dir} -A mafft -S diamond -I 4 -T iqtree3 --matrix PAM30 -s ${captainPhylogeny} -o ${output_dir} -n characterization &> ${working_dir}/orthofinder.log

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
    else
        orthofinder_flag=false
    fi
}

run_blast() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local temp_dir="${working_dir}/temp/"
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local output_dir="${working_dir}/Orthofinder"

    cat ${nucleotide_dir}/*.fa > ${temp_dir}/sequence.fasta

    makeblastdb -dbtype nucl -parse_seqids -in ${temp_dir}/sequence.fasta -out ${temp_dir}/Cluster &>/dev/null

    blastn -query ${temp_dir}/sequence.fasta -db ${temp_dir}/Cluster -evalue 1e-60 -num_threads "${threads}" -outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen" -task blastn -gapopen 8 -gapextend 6 -reward 5 -penalty -4 -out ${temp_dir}/blastresults.txt

    awk 'begin{fs=ofs="\t"}{if($8>=2000) {print}}' ${temp_dir}/blastresults.txt > ${working_dir}/Blast_CleanResults.txt
}

check_core() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local orthogroups_dir="${working_dir}/Orthofinder/Results_characterization/Orthogroup_Sequences/"
    local temp_dir="${working_dir}/temp/"

    ls ${working_dir}/Core_genes*.txt > ${temp_dir}/core_files.txt 2>/dev/null

    File_number=$(wc -l "${temp_dir}/core_files.txt" | awk '{print $1}')

    if [[ $File_number -gt 0 ]]; then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Core genes have been identified. Processing..."

        grep -w ${ClusterId} ${base_dir}/../ClustersAnalyzed.txt >> ${base_dir}/../ClusterCore.txt
        
        mkdir -p "${working_dir}/Core_genes"

        cat ${working_dir}/Core_genes*.txt | sort -u > ${temp_dir}/Full_core.txt

        while read Orthogroup; 
        do 
            cp ${orthogroups_dir}/${Orthogroup}.fa ${working_dir}/Core_genes/
        done < ${temp_dir}/Full_core.txt
    fi
}

check_movement() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/ClusterCharacterization/"
    local orthogroups_dir="${working_dir}/Orthofinder/Results_characterization/Orthogroup_Sequences/"
    local temp_dir="${working_dir}/temp/"

    ls ${working_dir}/*_moveOrthologs.txt > ${temp_dir}/movement_files.txt 2>/dev/null

    File_number=$(wc -l "${temp_dir}/movement_files.txt" | awk '{print $1}')

    if [[ $File_number -gt 0 ]]; then
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Genes in movement have been identified. Processing..."

        grep -w ${ClusterId} ${base_dir}/../ClustersAnalyzed.txt >> ${base_dir}/../ClusterMovement.txt

        cat ${temp_dir}/movement_files.txt | xargs -n 1 basename -s .txt | while read SubCluster
        do
            mkdir -p "${working_dir}/${SubCluster}"
            while read Orthogroup; 
            do 
                cp ${orthogroups_dir}/${Orthogroup}.fa ${working_dir}/${SubCluster}/
            done < ${working_dir}/${SubCluster}.txt
        done
    fi
}

# ==============================================================================
# Start the process
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running ClusterCharacterization Module."
check_clusters "${Working_directory}"

awk '{print $1}' ${Working_directory}/Clusters/ClustersAnalyzed.txt | sed $'s/[^[:print:]\t]//g' | while read ClusterId
do
    internal_dir="${Working_directory}/Clusters/${ClusterId}/"
    subcluster_number=$(grep -w ${ClusterId} ${Working_directory}/Clusters/ClustersAnalyzed.txt | awk '{print $3}' | sed $'s/[^[:print:]\t]//g')
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
    check_directory_structure "${internal_dir}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure is valid. Proceeding."
    
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
    organize_working_directory "${internal_dir}"

    # ==============================================================================
    # Running orthofinder
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running orthofinder..."
    run_orthofinder "${internal_dir}"
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
    if ! $orthofinder_flag; then
        echo -e "  \033[01;31mERROR\033[m: There was an error with orthofinder in this cluster.\n"
        rm -r "${internal_dir}/Workspace/ClusterCharacterization/"
        continue
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    # ==============================================================================
    # Running blast
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running blast for nucleotide synteny..."
    run_blast "${internal_dir}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

    # ==============================================================================
    # Running characterization
    # ==============================================================================

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running characterization of the cluster..."
    Rscript ${auxiliary_path}/ClusterAnalysis.R -d "${internal_dir}/Workspace/ClusterCharacterization/" -s "${subcluster_number}" -c $captainremoval_number

    mkdir -p "${internal_dir}/Workspace/ClusterCharacterization/Images"
    mv ${internal_dir}/Workspace/ClusterCharacterization/*.svg "${internal_dir}/Workspace/ClusterCharacterization/Images/"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking for core genes..."
    check_core "${internal_dir}"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking for identifiable genes that participate in movement..."
    check_movement "${internal_dir}"
    
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding.\n"
done
echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"