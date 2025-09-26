#!/bin/bash

# Function to print help message

function print_help() {
   echo -e "Script to run captain identification in each element and phylogenetic tree of this captains.
   This script perform five steps:
   1. Run hmmsearch against protein database for each element.
   2. Process the data and identify captains and regions used for phylogenetic analysis. For the captain identification, there are tree level of minimum confidence that can be used:
      2.1 Only a match with the catain hmm profile from starfish (\033[01;31mWARNING\033[m: This could lead to false positive identification leading to a bad phylogenetic analysis).
      2.2 Plus a match with the DUF3435 hmm profile.
      2.3 Plus a match with the integrase catalitic core hmm profile.
   3. Align exonic sequence with MACSE with aminoacid output and preprocess alignment with Clipkit.
   4. Run maximum-likelihood phylogenetic tree inference.
   
   The working directory should have the following structure:
   WorkingDiretory/
   ├── *_gff/
   │   ├── element01.gff
   │   └── element02.gff
   │   ︙
   └── *_protein/
       ├── element01.fa
       └── element02.fa
       ︙"
   echo
   echo "Syntax: $0 [ -h ] -w <working_directory> [ -c <confidence_level> -t <num_threads> ]"
   echo "options:"
   echo "-w, --workingDirectory: Specify the working directory where all data are stored (required)."
   echo "-l, --length: minimum length of the protein to be identify as captain [range: 300 - 800] (Default: 500)."
   echo "-c, --confidenceLevel: Minimum confidence level to call a captain. Note: the script is always going to try to return the captain with the highest level of confidence [range: 1 - 3] (Default: 2)"
   echo "-t, --threads: Tnumber of threads to use for alignment and phylogenetic tree inference (Default: 1)."
   echo "-help: Display this help message."
}

# Initialize variables

workingDirectory_path=""
hmmprofile_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../hmm/"
auxiliary_path="$( dirname -- "$( readlink -f -- "$0"; )"; )""/../aux/"
level="2"
length="500"
threads="1"
help_flag=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workingDirectory)
            shift
            workingDirectory_path="$1"
            ;;
        -h|--hmm)
            shift
            hmmprofile_path="$1"
            ;;
        -c|--confidenceLevel)
            shift
            level="$1"
            ;;
        -l|--length)
            shift
            length="$1"
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
if [[ -z "$workingDirectory_path" ]]; then
    echo "Error: Missing required arguments."
    print_help
    exit 1
fi

# Check if working directory exists
if [[ ! -d "$workingDirectory_path" ]]; then
    echo "Error: Directory '$workingDirectory_path' does not exist."
    exit 1
else
    workingDirectory_path=$(realpath $workingDirectory_path)
fi

# Check if hmmprofile directory exists
if [[ ! -d "$hmmprofile_path" ]]; then
    echo "Error: Directory '$hmmprofile_path' does not exist."
    exit 1
else
    hmmprofile_path=$(realpath $hmmprofile_path)
    if [[ ! -f "${hmmprofile_path}/CAT_domain.hmm" ]]; then
        echo "Error: File '${hmmprofile_path}/CAT_domain.hmm' does not exist."
        exit 1
    else
        CAT_hmm="${hmmprofile_path}/CAT_domain.hmm"
    fi
    if [[ ! -f "${hmmprofile_path}/DUF3435.hmm" ]]; then
        echo "Error: File '${hmmprofile_path}/DUF3435.hmm' does not exist."
        exit 1
    else
        DUF_hmm="${hmmprofile_path}/DUF3435.hmm"
    fi
    if [[ ! -f "${hmmprofile_path}/Captain.hmm" ]]; then
        echo "Error: File '${hmmprofile_path}/Captain.hmm' does not exist."
        exit 1
    else
        CAPTAIN_hmm="${hmmprofile_path}/Captain.hmm"
    fi
fi


# Check if the python function file exists
if [[ ! -d "$auxiliary_path" ]]; then
    echo "Error: directory '$auxiliary_path' does not exist."
    exit 1
else
    auxiliary_path=$(realpath $auxiliary_path)
fi

# Check confidence level is within allowed range
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

# Check length is within allowed range
if [[ "$length" =~ ^[0-9]+$ ]]; then
    if (( $length < 300 || $length > 800 )); then
        echo "Error: '$length' length is not an accepted value."
        print_help
        exit 1
    fi
else
    echo "Error: '$length' is not a positive integer."
    print_help
    exit 1
fi

# Check thread parameter
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: '$threads' is not a positive integer."
    print_help
    exit 1
fi

# Check for software presence
if [[ -z "$(which hmmsearch)" ]]; then
    echo "Error: Missing hmmsearch function."
    exit 1
fi

if [[ -z "$(which seqkit)" ]]; then
    echo "Error: Missing seqkit function."
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

if [[ -z "$(which clipkit)" ]]; then
    echo "Error: Missing clipkit function."
    exit 1
fi

if [[ -z "$(which gotree)" ]]; then
    echo "Error: Missing gotree function."
    exit 1
fi

# ==============================================================================
# Bash function block
# ==============================================================================

check_directory_structure() {
    local base_dir="$1"
    
    # Step 1: Locate required subdirectories and file
    local gff_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local protein_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    local nucleotide_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local exon_dir=$(find "$base_dir" -maxdepth 1 -type d -name "*_exon" 2>/dev/null)

    if [[ -z "$gff_dir" ]]; then
        echo "Error: GFF subdirectory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$protein_dir" ]]; then
        echo "Error: Protein subdirectory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$nucleotide_dir" ]]; then
        echo "Error: nucleotide subdirectory not found in '$base_dir'." >&2
        exit 1
    fi

    if [[ -z "$exon_dir" ]]; then
        echo "Error: exon subdirectory not found in '$base_dir'." >&2
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
}

process_hmmsearch() {
    local working_dir="$1"

    # Locate the protein directory based on the pattern *_protein
    local protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    
    # Define path for the output directory
    local NAME=$(basename "$working_dir")
    local hmmer_results_prefix="${working_dir}/${NAME}_Hmmsearch"
    local hmmer_results_CAT="${hmmer_results_prefix}CAT"
    local hmmer_results_DUF="${hmmer_results_prefix}DUF"
    local hmmer_results_CAPTAIN="${hmmer_results_prefix}Captain"

    # Create the necessary directoriy
    mkdir -p "${hmmer_results_CAT}"
    mkdir -p "${hmmer_results_DUF}"
    mkdir -p "${hmmer_results_CAPTAIN}"

    # Step 1: Perform hmmsearch
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
        hmmsearch --max --noali --cpu ${threads} --domE 10e-6 --domtblout ${outfileCAPTAIN} ${CAPTAIN_hmm} ${query} >/dev/null
        hmmsearch --max --noali --cpu ${threads} --domE 10e-6 --domtblout ${outfileDUF} ${DUF_hmm} ${query} >/dev/null
        hmmsearch --max --noali --cpu ${threads} --domE 10e-6 --domtblout ${outfileCAT} ${CAT_hmm} ${query} >/dev/null
    done

    echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Analysis complete. Results stored in '$working_dir'."
    return 0
}

Captain_identification() {
    local working_dir="$1"

    local gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local CAPTAIN_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_HmmsearchCaptain" 2>/dev/null)
    local CAT_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_HmmsearchCAT" 2>/dev/null)
    local DUF_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_HmmsearchDUF" 2>/dev/null)
    local results_path="${working_dir}/Captains.txt"
    local empty_elements="${working_dir}/EmptyElements.txt"

    python ${auxiliary_path}/hmmer_process.py --hmm1 "${CAPTAIN_path}" --hmm2 "${DUF_path}" --hmm3 "${CAT_path}" --gff "${gff_dir}" --fasta "${nucleotide_dir}" --output "${results_path}" --empty "${empty_elements}" --min_common "${level}" --min_length "${length}"

}

Alignment() {
    local working_dir="$1"

    # Locate required subdirectories and define output path
    local protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    local nucleotide_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local exon_path=$(find "$working_dir" -maxdepth 1 -type d -name "*_exon" 2>/dev/null)
    local captain_file="$working_dir/Captains.txt"
    local temp_prefix="${working_dir}/temp_"
    local empty_elements="${working_dir}/EmptyElements.txt"
    local pseudoExons="${working_dir}/Captains_pseudo.fa"
    local Remove_elements="${working_dir}/Remove_elements.txt"
 

    if [ ! -s "${captain_file}" ]; then
        echo "ERROR: there is no captain identify in this set of data"
        exit 1
    fi

    # Select captains that were correctly identify 
    cat ${exon_path}/*.fa > ${temp_prefix}exon.fa
    awk '{print $1}' ${captain_file} > ${temp_prefix}CaptainIDs.txt
    seqkit grep --quiet -f ${temp_prefix}CaptainIDs.txt ${temp_prefix}exon.fa -o ${working_dir}/Captains_exon.fa

    if [ -s "${empty_elements}" ]; then
        echo -e "[$(date "+%Y-%m-%d %H:%M:%S")] \033[01;31mWARNING\033[m: Not all elements have an identifible captain. Trying to identify a region possibly associated to a captain pseudogene..."

        #Search for possible pseudogenes in the empty elements
        mkdir -p ${temp_prefix}blast
        cat ${empty_elements} | while read line 
        do 
            echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Looking for captain pseudogene in element '${line}'"

            seqkit subseq --quiet -r 1:20000 ${nucleotide_path}/${line}.fa > ${temp_prefix}blast/${line}_start.fa

            makeblastdb -in ${temp_prefix}blast/${line}_start.fa -dbtype nucl -out ${temp_prefix}blast/${line}_start >/dev/null

            blastn -query ${working_dir}/Captains_exon.fa -db ${temp_prefix}blast/${line}_start -outfmt "6 sseqid sstart send" | sort -k2 -n -u | awk -F'\t' '
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
            }' > ${temp_prefix}${line}.bed

            local exonNumber=$(wc -l ${temp_prefix}${line}.bed | awk '{print $1}')

            if [[ ${exonNumber} -lt 3 ]]; then

                seqkit subseq --quiet -r -20000:-1 ${nucleotide_path}/${line}.fa | seqkit seq --quiet --reverse --complement -v --seq-type dna > ${temp_prefix}blast/${line}_end.fa
                makeblastdb -in ${temp_prefix}blast/${line}_end.fa -dbtype nucl -out ${temp_prefix}blast/${line}_end >/dev/null

                blastn -query ${working_dir}/Captains_exon.fa -db ${temp_prefix}blast/${line}_end -outfmt "6 sseqid sstart send" | sort -k2 -n -u | awk -F'\t' '
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
                }' > ${temp_prefix}${line}-2.bed

                local exonNumber=$(wc -l ${temp_prefix}${line}-2.bed | awk '{print $1}')

                if [[ ${exonNumber} -ge 3 ]]; then
                    seqkit subseq --quiet --bed ${temp_prefix}${line}-2.bed ${nucleotide_path}/${line}.fa  | grep -v ">" | sed -z  's/\n//g' | sed "1i >${line}" | sed -e '$a\' >> ${pseudoExons}
                else 
                    echo -e "  \033[01;31mWARNING\033[m: Element \033[1m'${line}'\033[m do not have an identifiable confident pseudogene. It will be removed from the final alignment."
                    echo "${line}" >> ${Remove_elements}
                fi
            else
                local exonNumber=$(wc -l ${temp_prefix}${line}.bed | awk '{print $1}')

                if [[ ${exonNumber} -ge 3 ]]; then
                    seqkit subseq --quiet --bed ${temp_prefix}${line}.bed ${nucleotide_path}/${line}.fa  | grep -v ">" | sed -z  's/\n//g' | sed "1i >${line}" | sed -e '$a\' >> ${pseudoExons}
                else 
                    echo "${exonNumber}"
                    echo -e "  \033[01;31mWARNING\033[m: Element '${line}' do not have an identifiable confident pseudogene. It will be removed from the final alignment."
                    echo "${line}" >> ${Remove_elements}

                fi
            fi 
        done

        rm ${nucleotide_path}/*seqkit*

        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing captain alignment..."

        if [ ! -s "${pseudoExons}" ]; then
            java -jar ${auxiliary_path}/macse.jar -prog alignSequences -seq ${working_dir}/Captains_exon.fa -out_AA ${working_dir}/Captain_proteins_aligned.fa >/dev/null
        else
            java -jar ${auxiliary_path}/macse.jar -prog alignSequences -seq ${working_dir}/Captains_exon.fa -seq_lr ${pseudoExons} -out_AA ${working_dir}/Captain_proteins_aligned.fa >/dev/null
        fi

    else
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Performing captain alignment..."

        java -jar ${auxiliary_path} -prog alignSequences -seq ${working_dir}/Captains_exon.fa -out_AA ${working_dir}/Captain_proteins_aligned.fa >/dev/null

        rm ${empty_elements}
    fi

    clipkit ${working_dir}/Captain_proteins_aligned.fa -m gappy -g 0.90 -l -q

    awk -F '.' '{print $1}' ${working_dir}/Captain_proteins_aligned.fa.clipkit > ${working_dir}/Captain_proteins_aligned_trimmed.fa

    rm -r ${temp_prefix}*
    rm ${working_dir}/Captain_proteins_aligned.fa.clipkit
}

Tree_inference() {
    local working_dir="$1"
    local protein_file="${working_dir}/Captain_proteins_aligned_trimmed.fa"
    local outdir="${working_dir}/captainPhylogeny"
    local temp_prefix="${working_dir}/temp_"

    # --- Error Handling ---
    if [[ ! -f "$protein_file" ]]; then
        echo "Error: Input file '$protein_file' not found." >&2
        exit 1
    fi

    seqkit rmdup --quiet -s ${protein_file} -o ${temp_prefix}unique &> /dev/null

    local unique_sequences=$(grep -c ">" ${temp_prefix}unique)

    if [[ "$unique_sequences" -ge 4 ]]; then
        mkdir -p ${outdir}
        iqtree3 -T ${threads} -m MFP --prefix ${outdir}/Captain_tree -B 1000 --alrt 1000 -s ${protein_file} -quiet --polytomy

        gotree collapse length -l 0.00001 -i ${outdir}/Captain_tree.treefile -o ${temp_prefix}captainlength.nw
        gotree collapse support -s 80 -i ${temp_prefix}captainlength.nw -o ${temp_prefix}captainsupport.nw

        sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_prefix}captainsupport.nw
        gotree collapse support -s 95 -i ${temp_prefix}captainsupport.nw -o ${temp_prefix}captainsupport2.nw
        sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' ${temp_prefix}captainsupport2.nw
        mv ${temp_prefix}captainsupport2.nw ${temp_prefix}captainsupport.nw

        gotree reroot midpoint -i ${temp_prefix}captainsupport.nw -o ${working_dir}/CaptainPhylogeny.nw
    else
        echo "  [$(date "+%Y-%m-%d %H:%M:%S")] Skipping tree inference due to low number of unique captain sequences. There's only '$unique_sequences' unique sequence(s) in the current dataset"
    fi

    rm ${temp_prefix}*
}

Removed_empty_elements() {

    local working_dir="$1"

    local gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_gff" 2>/dev/null)
    local protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_protein" 2>/dev/null)
    local nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_nucleotide" 2>/dev/null)
    local exon_dir=$(find "$working_dir" -maxdepth 1 -type d -name "*_exon" 2>/dev/null)

    local Remove_elements="${working_dir}/Remove_elements.txt"
    local output="${working_dir}/Remove_elements/"

    mkdir -p "${output}"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Removing elements without suitable gene or pseudogene model of captain from the final dataset."

    cat "${Remove_elements}" | while read line
    do
        mv ${gff_dir}/${line}.gff ${output}
        mv ${protein_dir}/${line}.fa ${output}${line}_protein.fa
        mv ${nucleotide_dir}/${line}.fa ${output}${line}_nucleotide.fa
        mv ${exon_dir}/${line}.fa ${output}${line}_exon.fa
    done
}

# ==============================================================================
# Checking Working directory structure
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Running captain identification and phylogenetic tree reconstruction of elements."
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 1: Checking Working directory '${workingDirectory_path}' structure."
check_directory_structure "${workingDirectory_path}"
echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> The directory structure in '${workingDirectory_path}' is valid. Proceeding."

# ==============================================================================
# Preprocessing data
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 2: Perform hmmsearch profile."

process_hmmsearch "${workingDirectory_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 2 finished. Proceeding."

# ==============================================================================
# Identify captains
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 3: Identifying captains from hmmsearch results."

Captain_identification "${workingDirectory_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 3 finished. Proceeding."

# ==============================================================================
# Alignment
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 4: Group captains and performed alignment."

Alignment "${workingDirectory_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 4 finished. Proceeding."


# ==============================================================================
# Alignment
# ==============================================================================

echo "[$(date "+%Y-%m-%d %H:%M:%S")] Step 5: Perform phylogenetic tree inference of captains."

Tree_inference "${workingDirectory_path}"

echo "[$(date "+%Y-%m-%d %H:%M:%S")]  -> Step 5 finished."

if [ -s "${workingDirectory_path}/Remove_elements.txt" ]; then
    Removed_empty_elements "${workingDirectory_path}"
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${workingDirectory_path}."
    echo -e "\033[01;31mWARNING\033[m: Elements have been removed, please check file '${workingDirectory_path}/Remove_elements.txt' and folder '${workingDirectory_path}/RemovedElements'."
else
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished. All results are store in ${workingDirectory_path}."
fi

