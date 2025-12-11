#!/usr/bin/env bash

# ==============================================================================
# Software check block
# ==============================================================================

source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
    check_auxiliary_scripts "${auxiliary_path}" "$(basename -s .sh "$0" )"
fi

# ==============================================================================
# Function block
# ==============================================================================

# Function to print help message
function print_help() {
   echo -e "Script to run the characterization of each selected cluster.
   This script perform eight steps per cluster to analyzed:
   1. Identify orthogroups through OrthoFinder.
   2. Perform all-vs-all Blastn for synteny visualization.
   3. Perform a hierarchical clustering of the elements based on Orthogroup gene count including singletons.
   4. Determine the full conection of the cluster and create a cargo orthogroups heatmap and synteny image for the cluster.
   5. Identify possible individual nesting events inside the cluster.
   6. Identify core genes in the cluster in two ways:
     6.1. General core: orthogroups that are present in at least 80% of the elements in the cluster.
     6.2. Specific core: Orthogroups that are present in at least 80% of the elements for subclusters generated at a 0.8 height of the hierarchical tree of cargo content.
       6.2.1. Divide the Cluster in subclusters of a height above 0.8 in the hierarchical clustering.
       6.2.2. If subslusters are presen, identify core genes in each one that have at least 5 elements using the same logic of general core.
   7. If subclusters are present it try to identify putative cargo movement events and try to avoid 'General core' genes.
   8. Determine if there are discordances at 'Clade' lavel between Cargo hierarchical clustering and Captain phylogenetic tree.
   "
   echo
   echo "Syntax: StarClust $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -t <integer> ]"
   echo ""
   echo "Required args:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored."
   echo ""
   echo "Required args with Default:"
   echo "-t, --threads: Number of threads for orthofinder and blast (Default: 8)"
   echo ""
   echo "Optional args:"
   echo "--overwrite: Flag to overwrite in case there is already a previous run of $(basename -s .sh "$0" ) (Default: off)"
   echo "-help: Display this help message."
}

# Modify this part when I have finish the SyntenyClustering fix
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

check_captain_information() {
    local base_dir="$1"   

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

    #if [[ -d "$working_dir" ]]; then
    #    echo "Error: There's a previous run in the Workspace."
    #    echo "If you want to overwrite this previous run, add the '--overwrite' flag to the command line."
    #    exit 1
    #fi

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

    awk 'begin{fs=ofs="\t"}{if($8>=1000) {print}}' ${temp_dir}/blastresults.txt > ${working_dir}/Blast_CleanResults.txt
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

organize_information() {
    local base_dir="$1"
    local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0" )/"

    local Characterization_dir="${base_dir}/$(basename -s .sh "$0" )/"

    mkdir -p ${Characterization_dir}

    if [[ ! -d "${working_dir}" ]]; then
        echo "Error in working directory"
    fi

    if [[ -f "${working_dir}/CargoHC.nwk" ]]; then
        cp ${working_dir}/CargoHC.nwk ${Characterization_dir}/CargoHierarchicalTree.nwk
    fi

    if [[ -d "${working_dir}/Images/" ]]; then
        cp -r ${working_dir}/Images/ ${Characterization_dir}/
    fi

    if [[ -d "${working_dir}/Core_genes/" ]]; then
        mkdir -p ${Characterization_dir}/Core/
        cp -r ${working_dir}/Core_genes/ ${Characterization_dir}/Core/Orthogroups/
        cp ${working_dir}/Core_genes-* ${Characterization_dir}/Core/
    fi

    local movingFolders=$(ls "${working_dir}/*_moveOrthologs.txt" | wc -l | awk '{print $1}')

    if [[ $movingFolders -ge 1 ]]; then
        mkdir -p ${Characterization_dir}/Movement_genes/
        ls ${working_dir}/*_moveOrthologs.txt | xargs -n1 basename -s .txt | while read line
        do
            local SubCluster=$(echo $line | awk -F '_' '{print $1}')
            mkdir -p ${Characterization_dir}/Movement_genes/${SubCluster}/
            cp -r ${working_dir}/${SubCluster}_moveOrthologs/ ${Characterization_dir}/Movement_genes/${SubCluster}/Orthogroups
            cp ${SubCluster}_moveOrthologs.txt ${Characterization_dir}/Movement_genes/${SubCluster}/OrthogroupsID.txt
            cp ${SubCluster}_moveOrthologsTable.csv ${Characterization_dir}/Movement_genes/${SubCluster}/Matrix.txt
        done
    fi
}

# ==============================================================================
# Variables block
# ==============================================================================

# Initialize variables
Working_directory=""
threads="8"
overwrite=false
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
        --overwrite)
            overwrite=true
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
echo "  Number of threads: " "$threads"
echo "  Overwrite previous run: " "$overwrite"
echo ""

# ==============================================================================
# Check variables block
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking arguments and input files..."

# Check for mandatory argument and define the path as absolute
if [[ -z "$Working_directory" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

if [[ ! -d "$Working_directory" ]]; then
    echo "Error: Directory '$Working_directory' does not exist."
    exit 1
else
    Working_directory=$(realpath $Working_directory)
fi
check_directory_structure "${Working_directory}"

# Check thread parameter
check_threads "$threads"

# Check for software presence
check_required_software "$(basename -s .sh "$0" )"

# ==============================================================================
# Main Block
# ==============================================================================

check_clusters "${Working_directory}"

# Modify this part when I have finish the SyntenyClustering fix
awk '{print $1}' ${Working_directory}/Clusters/ClustersAnalyzed.txt | sed $'s/[^[:print:]\t]//g' | while read ClusterId
do
    internal_dir="${Working_directory}/Clusters/${ClusterId}/"
    subcluster_number=$(grep -w ${ClusterId} ${Working_directory}/Clusters/ClustersAnalyzed.txt | awk '{print $3}' | sed $'s/[^[:print:]\t]//g')
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
    check_directory_structure "${internal_dir}"
    check_captain_information "${internal_dir}"

    if $overwrite; then
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking and removing previous run if exists..."
        overwrite "${Working_directory}" "$(basename -s .sh "$0" )"
    fi
    
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace..."
    organize_working_directory "${internal_dir}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Running orthofinder..."
    run_orthofinder "${internal_dir}"
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
    if ! $orthofinder_flag; then
        echo -e "  \033[01;31mERROR\033[m: There was an error with orthofinder in this cluster.\n"
        rm -r "${internal_dir}/Workspace/ClusterCharacterization/"
        continue
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Running blast for nucleotide synteny..."
    run_blast "${internal_dir}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Running characterization of the cluster..."
    Rscript ${auxiliary_path}/ClusterAnalysis.R -d "${internal_dir}/Workspace/ClusterCharacterization/" -s "${subcluster_number}" -c $captainremoval_number

    mkdir -p "${internal_dir}/Workspace/ClusterCharacterization/Images"
    mv ${internal_dir}/Workspace/ClusterCharacterization/*.svg "${internal_dir}/Workspace/ClusterCharacterization/Images/"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking for core genes..."
    check_core "${internal_dir}"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking for identifiable genes that participate in movement..."
    check_movement "${internal_dir}"
    
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing information to the main directory."
    organize_information "${internal_dir}"
    rm -r "${internal_dir}/Workspace/$(basename -s .sh "$0" )/temp/" > /dev/null
    echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Proceeding.\n"
done
echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"