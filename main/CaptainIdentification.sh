#!/usr/bin/env bash

# ==============================================================================
# Software check block
# ==============================================================================

source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Utils.sh"
source "$( dirname -- "$( readlink -f -- "$0"; )"; )""/../lib/Check.sh"

auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
    check_auxiliary_scripts "${auxiliary_path}" "$(basename -s .sh "$0" )"
fi

hmmprofile_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../hmm/"
if [[ ! -d "$hmmprofile_path" ]]; then
    echo "Error: directory '$hmmprofile_path' does not exist."
    exit 1
else
    hmmprofile_path=$(realpath $hmmprofile_path)
    if [[ ! -f "${hmmprofile_path}/CAT_domain.hmm" || ! -f "${hmmprofile_path}/CAT_domain.hmm.h3f" || ! -f "${hmmprofile_path}/CAT_domain.hmm.h3i" || ! -f "${hmmprofile_path}/CAT_domain.hmm.h3m" || ! -f "${hmmprofile_path}/CAT_domain.hmm.h3p" ]]; then
        echo "Error: binary compressed datafiles for hmmscan of 'CAT_domain' does not exist."
        exit 1
    else
        CAT_hmm="${hmmprofile_path}/CAT_domain.hmm"
    fi
    if [[ ! -f "${hmmprofile_path}/DUF3435.hmm" || ! -f "${hmmprofile_path}/DUF3435.hmm.h3f" || ! -f "${hmmprofile_path}/DUF3435.hmm.h3i" || ! -f "${hmmprofile_path}/DUF3435.hmm.h3m" || ! -f "${hmmprofile_path}/DUF3435.hmm.h3p" ]]; then
        echo "Error: binary compressed datafiles for hmmscan of 'DUF3435' does not exist."
        exit 1
    else
        DUF_hmm="${hmmprofile_path}/DUF3435.hmm"
    fi
    if [[ ! -f "${hmmprofile_path}/Captain.hmm" || ! -f "${hmmprofile_path}/Captain.hmm.h3f" || ! -f "${hmmprofile_path}/Captain.hmm.h3i" || ! -f "${hmmprofile_path}/Captain.hmm.h3m" || ! -f "${hmmprofile_path}/Captain.hmm.h3p" ]]; then
        echo "Error: binary compressed datafiles for hmmscan of 'Captain' does not exist."
        exit 1
    else
        CAPTAIN_hmm="${hmmprofile_path}/Captain.hmm"
    fi
fi

database_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../databases/"
if [[ ! -d "$database_path" ]]; then
    echo "Error: directory '$database_path' does not exist."
    exit 1
else
    database_path=$(realpath $database_path)
    # Check the presence of the specific scripts
    if [[ ! -f "${database_path}/Captains_exon.fa" ]]; then
        echo "Error: file '${database_path}/Captains_exon.fa' does not exist."
        exit 1
    else
        check_fasta_dna "${database_path}/Captains_exon.fa"
    fi
fi

agat_config="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../agat_config.yaml"
if [[ ! -f "$agat_config" ]]; then
    echo "Error: file '$agat_config' does not exist."
    exit 1
else
    agat_config=$(realpath $agat_config)
fi

# ==============================================================================
# Function block
# ==============================================================================

# Function to print help message
function print_help() {
   echo -e "Script to identify captain genes within each element and construct a phylogenetic tree based on these captains.
   It executes five main steps:
   1. Run hmmsearch using hmm profiles of captain domains against the proteome of each element.
   2. Processes the data to identify Captains and regions suitable for phylogenetic analysi. Three minimum confidence levels can be used for Captain identification:
      2.1 Only a match with the Captain HMM profile from Starfish (\033[01;31mWARNING\033[m: This may lead to false positive identifications and result in an unreliable phylogenetic analysis).
      2.2 A match with the Captain HMM profile plus a match with the DUF3435 HMM profile.
      2.3 A match with the Captain HMM profile and DUF3435 HMM profile plus a match with the Integrase catalytic core HMM profile.
   3. For elements lacking a confident Captain gene, the script searches for a putative Captain pseudogene at the beginning and end of the element.
   4. Aligns exonic sequences using MACSE (with amino acid output) and preprocesses the alignment with Clipkit.
   5. Runs maximum-likelihood phylogenetic tree inference, when there at least two unique captain sequences:
      5.1 Run IQ-TREE with 1000 UFBotstrap and 1000 sh-aLRT if there at least 4 unique sequence, in other case run it without support.
      5.2 Collapsed near-zero and low-support (when available) branches.
      5.3 mid-root the tree.

   There are three available mode:
   -Cluster: Analyzes and performs all five steps per cluster, and remove elements without a suitable Captain gene or pseudogene from the main dataset.
   -FullAll: Analyzes and performs all five steps on the whole dataset, and remove elements without a suitable Captain gene or pseudogene from the main dataset..
   -AllID: Analyzes and performs only the first three steps (Identification) on the whole dataset, and remove elements without a suitable Captain gene or pseudogene from the main dataset.
   "
   echo
   echo "Syntax: SAT $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -l <integer> -c <integer> -m <string> -t <integer> ]"
   echo ""
   echo "Required args:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo ""
   echo "Required args with Default:"
   echo "-m, --mode: specified the mode (Default = Cluster) [Available mode: Cluster, FullAll, AllID]."
   echo "-l, --length: Minimum length of the protein to be identify as captain (Default: 250) [range: 200 - 800]."
   echo "-c, --confidenceLevel: Minimum confidence level to call a captain. Note: the script is always going to try to return the captain with the highest level of confidence (Default: 2) [range: 1 - 3]."
   echo "-r, --rangeKb:The distance (as a number of kilobases) from the beginning or end of the element within which a gene must fall to be considered a captain (Default: 10) [range: 3 - 20]"
   echo "-t, --threads: Number of threads to use for phylogenetic tree inference (Default: 1)."
   echo ""
   echo "Optional args:"
   echo "-help: Display this help message."
}

organize_working_directory() {
    local base_dir="$1"

    local data_dir="${base_dir}/Data/"

    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"
    local temp_dir="${base_dir}/Workspace/CaptainIdentification/temp/"

    # Create required subdirectories
    mkdir -p ${working_dir}
    mkdir -p ${temp_dir}

    # Copy require files
    cp -r ${gff_dir} ${working_dir}
    cp -r ${nucleotide_dir} ${working_dir}
    cp -r ${protein_dir} ${working_dir}
    cp -r ${exon_dir} ${working_dir}
}

process_hmmsearch() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"

    # Locate the protein directory
    local protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    
    # Define path for the output directory
    local hmmer_results_prefix="${working_dir}/Hmmsearch_"
    local hmmer_results_CAT="${hmmer_results_prefix}CAT"
    local hmmer_results_DUF="${hmmer_results_prefix}DUF"
    local hmmer_results_CAPTAIN="${hmmer_results_prefix}Captain"

    # Create the necessary output directory
    mkdir -p "${hmmer_results_CAT}"
    mkdir -p "${hmmer_results_DUF}"
    mkdir -p "${hmmer_results_CAPTAIN}"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing profile searches..."
    
    # Get a list of all species filenames without the .fa extension
    local fa_files=( $(find "$protein_path" -maxdepth 1 -type f -name "*.fa") )
    
    for fa_file in "${fa_files[@]}"
    do
        local species_name=$(basename "$fa_file" .fa)
        local query="$fa_file"
        local outfileCAPTAIN="${hmmer_results_CAPTAIN}/${species_name}.txt"
        local outfileDUF="${hmmer_results_DUF}/${species_name}.txt"
        local outfileCAT="${hmmer_results_CAT}/${species_name}.txt"

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Looking at element '${species_name}'"
        
        # Perform profile search
        hmmscan --max --noali --cpu ${threads} --domE 0.001 --domtblout ${outfileCAPTAIN} ${CAPTAIN_hmm} ${query} >/dev/null
        hmmscan --max --noali --cpu ${threads} --domE 0.001 --domtblout ${outfileDUF} ${DUF_hmm} ${query} >/dev/null
        hmmscan --max --noali --cpu ${threads} --domE 0.001 --domtblout ${outfileCAT} ${CAT_hmm} ${query} >/dev/null
    done

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Analysis complete. Results stored in '$working_dir'."
    return 0
}

Captain_identification() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"

    local gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local CAPTAIN_path=$(find "$working_dir" -maxdepth 1 -type d -name "Hmmsearch_Captain" 2>/dev/null)
    local CAT_path=$(find "$working_dir" -maxdepth 1 -type d -name "Hmmsearch_CAT" 2>/dev/null)
    local DUF_path=$(find "$working_dir" -maxdepth 1 -type d -name "Hmmsearch_DUF" 2>/dev/null)
    local results_path="${working_dir}/CaptainsID.txt"
    local empty_elements="${working_dir}/EmptyElements.txt"

    python ${auxiliary_path}/hmmscan_process.py --hmm1 "${CAPTAIN_path}" --hmm2 "${DUF_path}" --hmm3 "${CAT_path}" --gff "${gff_dir}" --fasta "${nucleotide_dir}" --output "${results_path}" --empty "${empty_elements}" --min_common "${level}" --min_length "${length}" --range_kb "${range}"
}

Captain_pseudogene() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"

    # Locate required subdirectories and define output path
    local nucleotide_path=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local exon_path=$(find "$working_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)
    local captain_file="$working_dir/CaptainsID.txt"
    
    local empty_elements="${working_dir}/EmptyElements.txt"
    local pseudoExons="${working_dir}/Captains_pseudo.fa"
    local Remove_elements="${working_dir}/Captainless_elements.txt"
    local temp_dir="${working_dir}/temp/"


    if [ ! -s "${captain_file}" ]; then
        echo -e "\033[01;31mWARNING\033[m: there is no captain gene identify in this set of data. Looking for pseudogenes only."
        cp ${database_path}/Captains_exon.fa ${temp_dir}/Captains_exon.fa
    else
        # Select captains that were correctly identify 
        cat ${exon_path}/*.fa > ${temp_dir}/exon.fa
        seqkit grep --quiet -f ${captain_file} ${temp_dir}/exon.fa -o ${working_dir}/Captains_exon.fa
        cat ${database_path}/Captains_exon.fa ${working_dir}/Captains_exon.fa > ${temp_dir}/Captains_exon.fa
    fi

    if [ -s "${empty_elements}" ]; then
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] \033[01;31mWARNING\033[m: Not all elements have an identifible captain. Trying to identify a region possibly associated to a captain pseudogene..."

        #Search for possible pseudogenes in the empty elements
        mkdir -p ${temp_dir}/blast
        cat ${empty_elements} | while read line 
        do 
            # Search for captain pseudogene at the beginning in the positive strand
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Looking for captain pseudogene in element '${line}'"

            seqkit subseq --quiet -r 1:20000 ${nucleotide_path}/${line}.fa > ${temp_dir}/blast/${line}_start.fa

            makeblastdb -in ${temp_dir}/blast/${line}_start.fa -dbtype nucl -out ${temp_dir}/blast/${line}_start >/dev/null

            blastn -query ${temp_dir}/Captains_exon.fa -db ${temp_dir}/blast/${line}_start -outfmt "6 sseqid sstart send" | awk 'BEGIN{FS=OFS="\t"}{if($2 < $3){print}}' | sort -k2 -n -u | awk -F'\t' '
            BEGIN {
            last_start = -1;
            last_end = -1;
            }
            {
                current_start = ($2 < $3) ? $2 : $3;
                current_end = ($2 < $3) ? $3 : $2;

                if (last_start == -1) {
                    last_start = current_start;
                    last_end = current_end;
                } else if (current_start <= last_end) {
                if (current_end > last_end) {
                    last_end = current_end;
                }
            } else {
                print $1 "\t" last_start "\t" last_end;
                last_start = current_start;
                last_end = current_end;
            }
            }
            END {
                if (last_start != -1) {
                    print $1 "\t" last_start "\t" last_end;
                }
            }' > ${temp_dir}/${line}_start.bed

            local exonNumber=$(wc -l ${temp_dir}/${line}_start.bed | awk '{print $1}')
            local exon_length=$(awk '{$4 = $3 - $2} {sum += $4} END {print sum}' ${temp_dir}/${line}_start.bed)

            if [[ ${exonNumber} -ge 2 && ${exon_length} -ge 600 ]]; then
                seqkit subseq --quiet --bed ${temp_dir}/${line}_start.bed ${nucleotide_path}/${line}.fa  | grep -v ">" | sed -z  's/\n//g' | sed "1i >${line}" | sed -e '$a\' >> ${pseudoExons}
            else
                # Search for captain pseudogene at the end in the negative strand
                seqkit subseq --quiet -r -20000:-1 ${nucleotide_path}/${line}.fa | seqkit seq --quiet --reverse --complement -v --seq-type dna > ${temp_dir}/blast/${line}_end.fa
                makeblastdb -in ${temp_dir}/blast/${line}_end.fa -dbtype nucl -out ${temp_dir}/blast/${line}_end >/dev/null

                blastn -query ${temp_dir}/Captains_exon.fa -db ${temp_dir}/blast/${line}_end -outfmt "6 sseqid sstart send" | awk 'BEGIN{FS=OFS="\t"}{if($2 < $3){print}}' | sort -k2 -n -u | awk -F'\t' '
                BEGIN {
                last_start = -1;
                last_end = -1;
                }
                {
                current_start = ($2 < $3) ? $2 : $3;
                current_end = ($2 < $3) ? $3 : $2;

                if (last_start == -1) {
                    last_start = current_start;
                    last_end = current_end;
                } else if (current_start <= last_end) {
                if (current_end > last_end) {
                    last_end = current_end;
                }
                } else {
                print $1 "\t" last_start "\t" last_end;
                last_start = current_start;
                last_end = current_end;
                }
                }
                END {
                if (last_start != -1) {
                    print $1 "\t" last_start "\t" last_end;
                }
                }' > ${temp_dir}/${line}_end.bed

                local exonNumber=$(wc -l ${temp_dir}/${line}_end.bed | awk '{print $1}')
                local exon_length=$(awk '{$4 = $3 - $2} {sum += $4} END {print sum}' ${temp_dir}/${line}_end.bed)

                if [[ ${exonNumber} -ge 2 && ${exon_length} -ge 600 ]]; then
                    seqkit subseq --quiet --bed ${temp_dir}/${line}_end.bed ${temp_dir}/blast/${line}_end.fa  | grep -v ">" | sed -z  's/\n//g' | sed "1i >${line}" | sed -e '$a\' >> ${pseudoExons}
                else 
                    echo -e "  \033[01;31mWARNING\033[m: Element \033[1m'${line}'\033[m do not have an identifiable confident pseudogene. It will be removed from the final alignment."
                    echo "${line}" >> ${Remove_elements}
                fi
            fi 
        done

        rm ${nucleotide_path}/*seqkit*
    else
        rm ${empty_elements}
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] All elements have a confident captain gene model"
    fi
}

Alignment() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"

    # Locate required subdirectories and define output path
    local pseudoExons="${working_dir}/Captains_pseudo.fa"
    local captain_file="${working_dir}/Captains_exon.fa"
    local temp_dir="${working_dir}/temp/"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing captain alignment..."

    if [[ ! -s "${pseudoExons}"  && -s "${captain_file}"  ]]; then
        java -jar ${auxiliary_path}/macse.jar -prog alignSequences -seq ${working_dir}/Captains_exon.fa -out_AA ${working_dir}/Captain_proteins_aligned.fa >/dev/null
        captainless_flag=false
    elif [[ -s "${pseudoExons}"  && -s "${captain_file}" ]]; then
        java -jar ${auxiliary_path}/macse.jar -prog alignSequences -seq ${working_dir}/Captains_exon.fa -seq_lr ${pseudoExons} -out_AA ${working_dir}/Captain_proteins_aligned.fa >/dev/null
        captainless_flag=false
    elif [[ -s "${pseudoExons}"  && ! -s "${captain_file}" ]]; then
        #Select a subset of captain exons from the big database that looks similar to our pseudogenes
        makeblastdb -in ${database_path}/Captains_exon.fa -dbtype nucl -out ${temp_dir}/blast/Captains_exon_database >/dev/null
        blastn -query ${pseudoExons} -db ${temp_dir}/blast/Captains_exon_database -task blastn -perc_identity 60 -qcov_hsp_perc 60 -evalue 0 -max_hsps 1 -outfmt "6 sseqid" | sort -u > ${temp_dir}/Selected_captain_exons.txt
        seqkit grep --quiet -f ${temp_dir}/Selected_captain_exons.txt ${database_path}/Captains_exon.fa -o ${temp_dir}/Selected_captain_exons.fa

        #Perform alignment with this selected sequences
        java -jar ${auxiliary_path}/macse.jar -prog alignSequences -seq ${temp_dir}/Selected_captain_exons.fa -seq_lr ${pseudoExons} -out_AA ${temp_dir}/Captain_proteins_aligned_pre.fa >/dev/null

        #Remove sequences from the database
        seqkit grep --quiet -v -f ${temp_dir}/Selected_captain_exons.txt ${temp_dir}/Captain_proteins_aligned_pre.fa -o ${working_dir}/Captain_proteins_aligned.fa
        captainless_flag=false
    elif [[ ! -s "${pseudoExons}" && ! -s "${captain_file}" ]]; then
        echo -e "\033[01;31mERROR\033[m:: there is no captain gene or pseudogene identify in this set of data"
        captainless_flag=true
        return 1
    fi

    clipkit ${working_dir}/Captain_proteins_aligned.fa -m gappy -g 0.90 -l -q

    awk -F '.' '{print $1}' ${working_dir}/Captain_proteins_aligned.fa.clipkit > ${working_dir}/Captain_proteins_aligned_trimmed.fa

    rm ${working_dir}/Captain_proteins_aligned.fa.clipkit
}

Tree_inference() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"

    local protein_file="${working_dir}/Captain_proteins_aligned_trimmed.fa"
    local outdir="${working_dir}/captainPhylogeny"
    local temp_dir="${working_dir}/temp/"

    if [[ ! -f "$protein_file" ]]; then
        echo "Error: Input file '$protein_file' not found." >&2
        exit 1
    fi

    # Check that all captains are not the same
    seqkit rmdup --quiet -s ${protein_file} -o ${temp_dir}/unique &> /dev/null

    local unique_sequences=$(grep -c ">" ${temp_dir}/unique)

    if [[ "$unique_sequences" -ge 4 ]]; then
        mkdir -p ${outdir}
        iqtree3 -T ${threads} -m MFP --prefix ${outdir}/Captain_tree -B 1000 --alrt 1000 -s ${protein_file} -quiet --polytomy

        # Collapse near-zero branches and split with low support
        gotree collapse length -l 0.00001 -i ${outdir}/Captain_tree.treefile -o ${temp_dir}/captainlength.nw
        gotree collapse support -s 80 -i ${temp_dir}/captainlength.nw -o ${temp_dir}/captainsupport.nw

        sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_dir}/captainsupport.nw
        gotree collapse support -s 95 -i ${temp_dir}/captainsupport.nw -o ${temp_dir}/captainsupport2.nw
        sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_dir}/captainsupport2.nw
        mv ${temp_dir}/captainsupport2.nw ${temp_dir}/captainsupport.nw

        # Root tree
        gotree reroot midpoint -i ${temp_dir}/captainsupport.nw -o ${base_dir}/CaptainPhylogeny.nw
    elif [[ "$unique_sequences" -ge 2 ]]; then
        echo -e "  [$(date "+%Y-%m-%d %H:%M:%S")] \033[01;31mWARNING\033[m: There's only '${unique_sequences}' unique sequences in the current dataset. Tree inference will be perform without support values."

        mkdir -p ${outdir}
        iqtree3 -T ${threads} -m MFP --prefix ${outdir}/Captain_tree -s ${protein_file} -quiet --polytomy

        # Collapse near-zero branches
        gotree collapse length -l 0.00001 -i ${outdir}/Captain_tree.treefile -o ${temp_dir}/captainlength.nw

        # Root tree
        gotree reroot midpoint -i ${temp_dir}/captainlength.nw -o ${base_dir}/CaptainPhylogeny.nw
    else
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] There's only one unique sequence in the current dataset. Skipping tree inference."
    fi
}

Removed_empty_elements() {
    local base_dir="$1"

    local working_dir="${base_dir}/Workspace/CaptainIdentification/"

    local data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

    local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
    local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
    local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
    local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

    local Remove_elements="${working_dir}/Captainless_elements.txt"
    local output="${data_dir}/Captainless_elements/"

    mkdir -p "${output}"

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Removing elements without suitable gene or pseudogene model of captain from the final dataset."

    cat "${Remove_elements}" | while read line
    do
        mv ${gff_dir}/${line}.gff ${output}
        mv ${protein_dir}/${line}.fa ${output}${line}_protein.fa
        mv ${nucleotide_dir}/${line}.fa ${output}${line}_nucleotide.fa
        mv ${exon_dir}/${line}.fa ${output}${line}_exon.fa
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
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing '${Cluster_number}' clusters.\n"
    fi
}

check_phylogeny() {
    local base_dir="$1"

    local phylogeny_file="${base_dir}/CaptainPhylogeny.nw"

    if [[ -f $phylogeny_file ]]; then
        phylogeny_flag=false
    else
        phylogeny_flag=true 
    fi
}

# ==============================================================================
# Variables block
# ==============================================================================

# Initialize variables
Working_directory=""
mode="Cluster"
length="250"
level="2"
range="10"
threads="1"
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
        -l|--length)
            shift
            length="$1"
            ;;
        -c|--confidenceLevel)
            shift
            level="$1"
            ;;
        -r|--rangeKb)
            shift
            range="$1"
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

# ==============================================================================
# Script start block
# ==============================================================================

echo "Running $(basename -s .sh "$0" ) command under the following parameters:"
echo "  Working directory: " "$Working_directory"
echo "  Mode: " "$mode"
echo "  Minimum number of anchor points: " "$anchorPoints"
echo "  Maximum number of gaps: " "$gaps"
echo "  Minimum number of nodes for Spectral clustering: " "$minNodes"
echo "  Minimum size of sub-cluster: " "$minSize"
echo "  Modularity score threshold for Spectral clustering: " "$threshold"
echo "  Number of threads: " "$threads"
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

# Check if mode parameter is correct
check_mode_parameter "$mode" "$(basename -s .sh "$0" )"

# Check that arguments are within allowed range
if [[ "$level" =~ ^[0-9]+$ ]]; then
    if (( $level < 1 || $level > 3 )); then
        echo "Error: '$level' confidence level is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$level' is not a positive integer."
    print_help
    exit 1
fi

if [[ "$length" =~ ^[0-9]+$ ]]; then
    if (( $length < 200 || $length > 800 )); then
        echo "Error: '$length' length is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$length' is not a positive integer."
    print_help
    exit 1
fi

if [[ "$range" =~ ^[0-9]+$ ]]; then
    if (( $range < 3 || $length > 20 )); then
        echo "Error: '$range' range is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$range' is not a positive integer."
    print_help
    exit 1
fi

if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
    exit 1
fi

# Check for required software
check_required_software "$(basename -s .sh "$0" )"

# ==============================================================================
# Main Block
# ==============================================================================

if [[ "${mode}" == "FullAll" ]]
then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"
    organize_working_directory "${Working_directory}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Perform hmmscan."
    process_hmmsearch "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Identifying captains from hmmscan results."
    Captain_identification "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Identifying if there are pseudogenes..."
    Captain_pseudogene "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding..."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Group captains and performed alignment."
    Alignment "${Working_directory}"
    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
    if $captainless_flag; then
        echo -e "  \033[01;31mERROR\033[m: There is no captain identify in this set of data."
        exit 1
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Perform phylogenetic tree inference of captains."
    Tree_inference "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished."

    if [ -s "${Working_directory}/Workspace/CaptainIdentification/Captainless_elements.txt" ]; then
        Removed_empty_elements "${Working_directory}"
        echo -e "  \033[01;31mWARNING\033[m: Elements have been removed, check this dataset."
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished.\n"
    else
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished.\n"
    fi
elif [[ "${mode}" == "AllID" ]]
then
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"
    organize_working_directory "${Working_directory}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Perform hmmscan profile."
    process_hmmsearch "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Identifying captains from hmmscan results."
    Captain_identification "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Identifying if there are pseudogenes..."
    Captain_pseudogene "${Working_directory}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding..."

    if [ -s "${Working_directory}/Workspace/CaptainIdentification/Captainless_elements.txt" ]; then
        Removed_empty_elements "${Working_directory}"
        echo -e "  \033[01;31mWARNING\033[m: Elements have been removed, check this dataset."
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished.\n"
    else
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished.\n"
    fi
elif [[ "${mode}" == "Cluster" ]]
then
    check_clusters "${Working_directory}"
    awk '{print $1}' ${Working_directory}/Clusters/SelectedClusters.txt | sed $'s/[^[:print:]\t]//g' | while read ClusterId
    do
        internal_dir="${Working_directory}/Clusters/${ClusterId}/"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Analyzing Cluster '$ClusterId'."
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Working directory '${internal_dir}' structure."
        check_directory_structure "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure in '${internal_dir}' is valid. Proceeding."

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Organizing workspace"
        organize_working_directory "${internal_dir}"

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Perform hmmscan."
        process_hmmsearch "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 1 finished. Proceeding."

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Identifying captains from hmmscan results."
        Captain_identification "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Identifying if there are pseudogenes..."
        Captain_pseudogene "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding..."

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Group captains and performed alignment."
        Alignment "${internal_dir}"
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
        if $captainless_flag; then
            echo -e "  \033[01;31mWARNING\033[m: There is no captain identify in this set of data.\n"
            continue
        fi
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Perform phylogenetic tree inference of captains."
        Tree_inference "${internal_dir}"
        echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished."

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Checking results.."
        check_phylogeny "${internal_dir}"
        if $phylogeny_flag; then
            echo -e "\033[01;31mWARNING\033[m: There's no captain phylogeny file, this cluster needs to be manually checked."
        else
            echo -e "[$(date "+%Y-%m-%d %H:%M:%S")]  Successfull run. Storing this cluster for further analysis."
            grep -w ${ClusterId} ${Working_directory}/Clusters/SelectedClusters.txt >> ${Working_directory}/Clusters/ClustersAnalyzed.txt
        fi

        if [ -s "${internal_dir}/Workspace/CaptainIdentification/Captainless_elements.txt" ]; then
            Removed_empty_elements "${internal_dir}"
            echo -e "  \033[01;31mWARNING\033[m: Elements have been removed, check this cluster."
            echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished.\n"
        else
            echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] Finished.\n"
        fi
    done
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] All clusters have been analyze"
fi
