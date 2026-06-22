#!/usr/bin/env bash
# ==============================================================================
# NAME:        CaptainIdentification.sh
# DESCRIPTION: Identifies captain genes/pseudogenes within each element and
#              constructs a phylogenetic tree based on those sequences.
#              Executes five main steps:
#              1. HMM profile search against each element proteome.
#              2. Captain identification from HMM results by position/confidence.
#              3. Pseudogene search in elements lacking a confident captain.
#              4. Sequence alignment with MACSE and trimming with ClipKit.
#              5. Maximum-likelihood phylogenetic tree inference with IQ-TREE.
# USAGE:       StarCrew CaptainIdentification [options]
#              StarCrew CaptainIdentification -help
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        04/May/2026
# VERSION:     1.0.0
# ==============================================================================

set -uo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

readonly SCRIPT_DIR="$( dirname -- "$( readlink -f -- "$0" )" )"
readonly LIB_DIR="${SCRIPT_DIR}/../lib"
readonly AUXILIARY_DIR="${SCRIPT_DIR}/../aux"
readonly HMM_DIR="${SCRIPT_DIR}/../hmm"
readonly DATABASE_DIR="${SCRIPT_DIR}/../databases"

# Source required libraries
if [[ ! -f "${LIB_DIR}/Utils.sh" ]]; then
  echo "Error: file '${LIB_DIR}/Utils.sh' does not exist." >&2
  exit 1
fi
source "${LIB_DIR}/Utils.sh"

if [[ ! -f "${LIB_DIR}/Check.sh" ]]; then
  echo "Error: file '${LIB_DIR}/Check.sh' does not exist." >&2
  exit 1
fi
source "${LIB_DIR}/Check.sh"

# Validate auxiliary directory and scripts
if [[ ! -d "$AUXILIARY_DIR" ]]; then
  echo "Error: directory '${AUXILIARY_DIR}' does not exist." >&2
  exit 1
fi
check_auxiliary_scripts "$(realpath "${AUXILIARY_DIR}")" "$(basename -s .sh "$0")"

# Validate required software
check_required_software "$(basename -s .sh "$0")"

# Validate HMM profile directory and required files
if [[ ! -d "$HMM_DIR" ]]; then
  echo "Error: HMM profile directory '${HMM_DIR}' not found." >&2
  exit 1
fi

readonly HMM_PATH="$(realpath "${HMM_DIR}")"
readonly HMM_DOMAINS=("CAT_domain" "DUF3435" "Captain")
readonly HMM_EXTENSIONS=("" ".h3f" ".h3i" ".h3m" ".h3p")

for dom in "${HMM_DOMAINS[@]}"; do
  for ext in "${HMM_EXTENSIONS[@]}"; do
    if [[ ! -f "${HMM_PATH}/${dom}.hmm${ext}" ]]; then
      echo "Error: missing HMM data for '${dom}' (${dom}.hmm${ext})." >&2
      exit 1
    fi
  done
done

readonly CAT_HMM="${HMM_PATH}/CAT_domain.hmm"
readonly DUF_HMM="${HMM_PATH}/DUF3435.hmm"
readonly CAPTAIN_HMM="${HMM_PATH}/Captain.hmm"

# Validate databases
check_databases "${DATABASE_DIR}" "$(basename -s .sh "$0")"
readonly DATABASE_PATH="$(realpath "${DATABASE_DIR}")"

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays the help message with all available options.
# Arguments:
#   None
# Returns:
#   0 always
print_help() {
  echo "Command to identify captain genes/pseudogenes and build a phylogenetic tree."
  echo "Executes five main steps:"
  echo "  1. HMM profile search against each element proteome."
  echo "  2. Captain identification from HMM results. Three confidence levels:"
  echo "     1 - Captain HMM match only."
  echo "         WARNING: May produce false positives and unreliable phylogenies."
  echo "     2 - Captain HMM + DUF3435 HMM match."
  echo "     3 - Captain HMM + DUF3435 HMM + YR Recombinase Active Site HMM match."
  echo "  3. Pseudogene search at element boundaries for captainless elements."
  echo "  4. MACSE alignment (AA output) trimmed with ClipKit."
  echo "  5. IQ-TREE phylogenetic inference (>=4 unique sequences: 1000 UFBootstrap"
  echo "     + sh-aLRT; 2-3 unique sequences: no support; 1 sequence: skipped)."
  echo
  echo "Available modes:"
  echo "  Cluster  - Runs all five steps per cluster; removes captainless elements."
  echo "  FullAll  - Runs all five steps on the full dataset; removes captainless elements."
  echo "             WARNING: MACSE may fail silently on large datasets."
  echo "  AllID    - Runs only steps 1-3 (identification) on the full dataset."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -w <directory_path>"
  echo "       [ -m <string> -l <integer> -c <integer> -r <integer>"
  echo "         -ms <integer> -t <integer> --overwrite ]"
  echo
  echo "Required args:"
  echo "  -w, --workingDirectory  Working directory where all data are stored."
  echo
  echo "Required args with defaults:"
  echo "  -m, --mode              Run mode (Default: AllID) [Available: Cluster, FullAll, AllID]."
  echo "  -l, --length            Minimum captain protein length in aa (Default: 250) [range: 200-800]."
  echo "  -c, --confidenceLevel   Minimum confidence level for captain calling (Default: 2) [range: 1-3]."
  echo "  -r, --rangeKb           Distance in kb from element boundary for captain search (Default: 10) [range: 3-20]."
  echo
  echo "Required args with defaults in 'Cluster' mode:"
  echo "  -ms, --minSize          Minimum cluster size to include in analysis (Default: 4) [range: 4-10]."
  echo
  echo "Optional args:"
  echo "  -t, --threads           Threads for hmmscan and IQ-TREE (Default: 1)."
  echo "  --overwrite             Overwrite a previous run (Default: off)."
  echo "  -help                   Display this help message."
}

# Checks for and loads captain data from a previous run if available.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 always; sets previous_captain_run=true if prior data exists
check_previous_information() {
  local base_dir="$1"

  captain_dir=$(find "$base_dir" -maxdepth 1 -type d -name "$(basename -s .sh "$0")" 2>/dev/null)

  if [[ ! -d "$captain_dir" ]]; then
    previous_captain_run=false
    return 0
  fi

  if [[ -f "${captain_dir}/Captains_CDS.fa" ]]; then
    check_fasta_dna "${captain_dir}/Captains_CDS.fa"
  fi

  if [[ -f "${captain_dir}/Captains_pseudo.fa" ]]; then
    check_fasta_dna "${captain_dir}/Captains_pseudo.fa"
  fi

  previous_captain_run=true
}

# Sets up the workspace by copying Data subdirectories and any prior captain files.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if a previous run exists without --overwrite
organize_working_directory() {
  local base_dir="$1"
  local data_dir="${base_dir}/Data"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local gff_dir nucleotide_dir protein_dir cds_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  if [[ -d "$working_dir" && ! $previous_captain_run ]]; then
    echo "Error: a previous run exists in '${working_dir}'." >&2
    echo "Use '--overwrite' to overwrite it." >&2
    exit 1
  elif [[ -d "$working_dir" && $previous_captain_run ]]; then
    rm -r "$working_dir"
  fi

  mkdir -p "$working_dir" "$temp_dir"
  cp -r "$gff_dir" "$working_dir"
  cp -r "$nucleotide_dir" "$working_dir"
  cp -r "$protein_dir" "$working_dir"
  cp -r "$cds_dir" "$working_dir"

  if $previous_captain_run; then
    local captain_dir="${base_dir}/$(basename -s .sh "$0")"
    [[ -f "${captain_dir}/Captains_CDS.fa" ]] && \
      cp "${captain_dir}/Captains_CDS.fa" "$working_dir"
    [[ -f "${captain_dir}/Captains_pseudo.fa" ]] && \
      cp "${captain_dir}/Captains_pseudo.fa" "$working_dir"
  fi
}

# Runs hmmscan for all three HMM profiles against each element proteome.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
process_hmmscan() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local protein_path
  protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)

  local hmmer_prefix="${working_dir}/Hmmsearch_"
  local hmmer_cat="${hmmer_prefix}CAT"
  local hmmer_duf="${hmmer_prefix}DUF"
  local hmmer_captain="${hmmer_prefix}Captain"

  mkdir -p "$hmmer_cat" "$hmmer_duf" "$hmmer_captain"

  log_step "Performing HMM profile searches..."

  local fa_files=()
  while IFS= read -r -d '' f; do
    fa_files+=("$f")
  done < <(find "$protein_path" -maxdepth 1 -type f -name "*.fa" -print0)

  local total_states=${#fa_files[@]}
  local state=0

  for fa_file in "${fa_files[@]}"; do
    state=$(( state + 1 ))
    local starship_name
    starship_name=$(basename "$fa_file" .fa)

    hmmscan --max --noali --cpu "${threads}" --domE 0.001 \
      --domtblout "${hmmer_captain}/${starship_name}.txt" \
      "${CAPTAIN_HMM}" "${fa_file}" > /dev/null
    hmmscan --max --noali --cpu "${threads}" --domE 0.001 \
      --domtblout "${hmmer_duf}/${starship_name}.txt" \
      "${DUF_HMM}" "${fa_file}" > /dev/null
    hmmscan --max --noali --cpu "${threads}" --domE 0.001 \
      --domtblout "${hmmer_cat}/${starship_name}.txt" \
      "${CAT_HMM}" "${fa_file}" > /dev/null

    ProgressBar "$state" "$total_states"
  done
  echo ""
  log_step "HMM profile searches complete."
}

# Identifies captain genes from hmmscan results using position and confidence level.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
captain_identification() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"

  local gff_dir captain_path cat_path duf_path
  gff_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  captain_path=$(find "$working_dir" -maxdepth 1 -type d -name "Hmmsearch_Captain" 2>/dev/null)
  cat_path=$(find "$working_dir" -maxdepth 1 -type d -name "Hmmsearch_CAT" 2>/dev/null)
  duf_path=$(find "$working_dir" -maxdepth 1 -type d -name "Hmmsearch_DUF" 2>/dev/null)

  python "${AUXILIARY_DIR}/hmmscan_process.py" \
    --hmm1 "${captain_path}" \
    --hmm2 "${duf_path}" \
    --hmm3 "${cat_path}" \
    --gff "${gff_dir}" \
    --fasta "${nucleotide_dir}" \
    --output "${working_dir}/CaptainsID.txt" \
    --empty "${working_dir}/EmptyElements.txt" \
    --not_passed "${working_dir}/EmptyElements-reason.txt" \
    --min_common "${level}" \
    --min_length "${length}" \
    --range_kb "${range}" > ${working_dir}/hmmscan_process_log.txt
}

# Searches for captain pseudogenes at element boundaries for captainless elements.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
captain_pseudogene() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local nucleotide_path cds_path
  nucleotide_path=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  cds_path=$(find "$working_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  local captain_file="${working_dir}/CaptainsID.txt"
  local empty_elements="${working_dir}/EmptyElements.txt"
  local pseudo_exons="${working_dir}/Captains_pseudo.fa"
  local remove_elements="${working_dir}/Captainless_elements.txt"
  local full_range=$(( range * 1000 ))

  if [[ ! -s "$captain_file" ]]; then
    log_warning "No captain gene identified in this dataset. Looking for pseudogenes only."
    cp "${DATABASE_PATH}/Captains_CDS.fa" "${temp_dir}/Captains_CDS.fa"
  else
    cat "${cds_path}"/*.fa > "${temp_dir}/CDS.fa"
    seqkit grep --quiet -f "$captain_file" "${temp_dir}/CDS.fa" \
      -o "${working_dir}/Captains_CDS.fa"
    cat "${DATABASE_PATH}/Captains_CDS.fa" "${working_dir}/Captains_CDS.fa" \
      > "${temp_dir}/Captains_CDS.fa"
  fi

  if [[ ! -s "$empty_elements" ]]; then
    rm "$empty_elements"
    log_step "All elements have a confident captain gene model."
    return 0
  fi

  log_warning "Not all elements have an identifiable captain. Searching for pseudogenes..."

  local total_states
  total_states=$(wc -l < "$empty_elements")
  local state=0

  mkdir -p "${temp_dir}/blast"

  while read -r element; do
    state=$(( state + 1 ))

    # Adjust range if element is shorter than 2x full range
    local element_length
    element_length=$(seqkit stats "${nucleotide_path}/${element}.fa" \
      | grep "DNA" | awk '{print $5}' | sed 's/,//')
    if (( element_length < full_range * 2 )); then
      full_range=$(( element_length / 2 ))
    fi

    # Search at the start of element (positive strand)
    seqkit subseq --quiet -r "1:${full_range}" \
      "${nucleotide_path}/${element}.fa" > "${temp_dir}/blast/${element}_start.fa"
    makeblastdb -in "${temp_dir}/blast/${element}_start.fa" \
      -dbtype nucl -out "${temp_dir}/blast/${element}_start" > /dev/null

    blastn -query "${temp_dir}/Captains_CDS.fa" \
      -db "${temp_dir}/blast/${element}_start" \
      -outfmt "6 sseqid sstart send" \
      | awk 'BEGIN{FS=OFS="\t"}{if($2 < $3){print}}' \
      | sort -k2 -n -u \
      | awk -F'\t' '
        BEGIN { last_start=-1; last_end=-1 }
        {
          current_start = ($2 < $3) ? $2 : $3
          current_end   = ($2 < $3) ? $3 : $2
          if (last_start == -1) {
            last_start = current_start; last_end = current_end
          } else if (current_start <= last_end) {
            if (current_end > last_end) last_end = current_end
          } else {
            print $1 "\t" last_start "\t" last_end
            last_start = current_start; last_end = current_end
          }
        }
        END { if (last_start != -1) print $1 "\t" last_start "\t" last_end }
      ' > "${temp_dir}/${element}_start.bed"

    local exon_number exon_length
    exon_number=$(wc -l < "${temp_dir}/${element}_start.bed")
    exon_length=$(awk '{sum += $3 - $2} END {print sum+0}' \
      "${temp_dir}/${element}_start.bed")

    if (( exon_number >= 2 && exon_length >= 600 )); then
      seqkit subseq --quiet --bed "${temp_dir}/${element}_start.bed" \
        "${nucleotide_path}/${element}.fa" \
        | grep -v ">" | sed -z 's/\n//g' \
        | sed "1i >${element}" | sed -e '$a\' >> "$pseudo_exons"
    else
      # Search at the end of element (negative strand)
      seqkit subseq --quiet -r "-${full_range}:-1" \
        "${nucleotide_path}/${element}.fa" \
        | seqkit seq --quiet --reverse --complement -v --seq-type dna \
        > "${temp_dir}/blast/${element}_end.fa"
      makeblastdb -in "${temp_dir}/blast/${element}_end.fa" \
        -dbtype nucl -out "${temp_dir}/blast/${element}_end" > /dev/null

      blastn -query "${temp_dir}/Captains_CDS.fa" \
        -db "${temp_dir}/blast/${element}_end" \
        -outfmt "6 sseqid sstart send" \
        | awk 'BEGIN{FS=OFS="\t"}{if($2 < $3){print}}' \
        | sort -k2 -n -u \
        | awk -F'\t' '
          BEGIN { last_start=-1; last_end=-1 }
          {
            current_start = ($2 < $3) ? $2 : $3
            current_end   = ($2 < $3) ? $3 : $2
            if (last_start == -1) {
              last_start = current_start; last_end = current_end
            } else if (current_start <= last_end) {
              if (current_end > last_end) last_end = current_end
            } else {
              print $1 "\t" last_start "\t" last_end
              last_start = current_start; last_end = current_end
            }
          }
          END { if (last_start != -1) print $1 "\t" last_start "\t" last_end }
        ' > "${temp_dir}/${element}_end.bed"

      exon_number=$(wc -l < "${temp_dir}/${element}_end.bed")
      exon_length=$(awk '{sum += $3 - $2} END {print sum+0}' \
        "${temp_dir}/${element}_end.bed")

      if (( exon_number >= 2 && exon_length >= 600 )); then
        seqkit subseq --quiet --bed "${temp_dir}/${element}_end.bed" \
          "${temp_dir}/blast/${element}_end.fa" \
          | grep -v ">" | sed -z 's/\n//g' \
          | sed "1i >${element}" | sed -e '$a\' >> "$pseudo_exons"
      else
        log_warning "Element '${element}' has no identifiable captain pseudogene. It will be removed."
        echo "${element}" >> "$remove_elements"
      fi
    fi

    full_range=$(( range * 1000 ))
    ProgressBar "$state" "$total_states"
  done < "$empty_elements"
  echo ""

  rm -f "${nucleotide_path}/"*seqkit* 2> /dev/null
}

# Aligns captain CDS and pseudogene sequences with MACSE and trims with ClipKit.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, 1 if no sequences are available or MACSE fails
run_alignment() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local pseudo_exons="${working_dir}/Captains_pseudo.fa"
  local captain_file="${working_dir}/Captains_CDS.fa"

  log_step "Performing captain alignment..."

  if [[ ! -s "$pseudo_exons" && -s "$captain_file" ]]; then
    macse -prog alignSequences \
      -seq "${captain_file}" \
      -out_AA "${working_dir}/Captain_proteins_aligned.fa" > /dev/null

  elif [[ -s "$pseudo_exons" && -s "$captain_file" ]]; then
    macse -prog alignSequences \
      -seq "${captain_file}" \
      -seq_lr "${pseudo_exons}" \
      -out_AA "${working_dir}/Captain_proteins_aligned.fa" > /dev/null

  elif [[ -s "$pseudo_exons" && ! -s "$captain_file" ]]; then
    # Select a subset of database captains similar to the pseudogenes
    makeblastdb -in "${DATABASE_PATH}/Captains_CDS.fa" \
      -dbtype nucl -out "${temp_dir}/blast/Captains_CDS_database" > /dev/null
    blastn -query "${pseudo_exons}" \
      -db "${temp_dir}/blast/Captains_CDS_database" \
      -task blastn -perc_identity 60 -qcov_hsp_perc 60 \
      -evalue 1e-5 -max_hsps 1 -outfmt "6 sseqid" \
      | sort -u > "${temp_dir}/Selected_captain_CDSs.txt"
    seqkit grep --quiet \
      -f "${temp_dir}/Selected_captain_CDSs.txt" \
      "${DATABASE_PATH}/Captains_CDS.fa" \
      -o "${temp_dir}/Selected_captain_CDSs.fa"

      if [[ ! -s "${temp_dir}/Selected_captain_CDSs.fa" ]]; then
        rm -r ${temp_dir}/blast ${temp_dir}/Selected_captain_CDSs*
        if [[ $mode == "Cluster" ]]; then
          captain_dir=$(find "${base_dir}/../../" -maxdepth 1 -type d -name "$(basename -s .sh "$0")" 2>/dev/null)
          if [[ ! -z $captain_dir ]]; then
            makeblastdb -in "${captain_dir}/Captains_CDS.fa" \
              -dbtype nucl -out "${temp_dir}/blast/Captains_CDS_database" > /dev/null
            blastn -query "${pseudo_exons}" \
              -db "${temp_dir}/blast/Captains_CDS_database" \
              -task blastn -perc_identity 60 -qcov_hsp_perc 60 \
              -evalue 1e-5 -max_hsps 1 -outfmt "6 sseqid" \
              | sort -u > "${temp_dir}/Selected_captain_CDSs.txt"
            seqkit grep --quiet \
              -f "${temp_dir}/Selected_captain_CDSs.txt" \
              "${captain_dir}/Captains_CDS.fa" \
              -o "${temp_dir}/Selected_captain_CDSs.fa"
          else
            log_warning "No captain gene available for this dataset alignment." >&2
            captainless_flag=true
            return 1     
          fi
        else
          echo "Error: no captain gene or pseudogene identified in this dataset." >&2
          captainless_flag=true
          return 1
        fi
      fi

    macse -prog alignSequences \
      -seq "${temp_dir}/Selected_captain_CDSs.fa" \
      -seq_lr "${pseudo_exons}" \
      -out_AA "${temp_dir}/Captain_proteins_aligned_pre.fa" > /dev/null

    # Remove database sequences from final alignment
    seqkit grep --quiet -v \
      -f "${temp_dir}/Selected_captain_CDSs.txt" \
      "${temp_dir}/Captain_proteins_aligned_pre.fa" \
      -o "${working_dir}/Captain_proteins_aligned.fa"

  else
    log_warning "No captain gene available for this dataset alignment." >&2
    captainless_flag=true
    return 1
  fi

  if [[ ! -f "${working_dir}/Captain_proteins_aligned.fa" ]]; then
    echo "Error: no captain gene or pseudogene identified in this dataset." >&2
    captainless_flag=true
    return 1
  fi

  clipkit "${working_dir}/Captain_proteins_aligned.fa" -m gappy -g 0.90 -l -q
  awk -F '.' '{print $1}' "${working_dir}/Captain_proteins_aligned.fa.clipkit" \
    > "${working_dir}/Captain_proteins_aligned_trimmed.fa"
  rm "${working_dir}/Captain_proteins_aligned.fa.clipkit"
}

# Runs IQ-TREE phylogenetic inference and post-processes the resulting tree.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
run_tree_inference() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local protein_file="${working_dir}/Captain_proteins_aligned_trimmed.fa"
  local outdir="${working_dir}/captainPhylogeny"

  if [[ ! -f "$protein_file" ]]; then
    echo "Error: input file '${protein_file}' not found." >&2
    return 1
  fi

  seqkit rmdup --quiet -s "$protein_file" -o "${temp_dir}/unique" &> /dev/null

  local unique_sequences
  unique_sequences=$(grep -c ">" "${temp_dir}/unique")

  if (( unique_sequences >= 4 )); then
    mkdir -p "$outdir"
    iqtree3 -T "${threads}" -m MFP \
      --prefix "${outdir}/Captain_tree" \
      -B 1000 --alrt 1000 \
      -s "$protein_file" -quiet --polytomy

    gotree collapse length -l 0.00001 \
      -i "${outdir}/Captain_tree.treefile" \
      -o "${temp_dir}/captainlength.nw"

    gotree collapse support -s 80 \
      -i "${temp_dir}/captainlength.nw" \
      -o "${temp_dir}/captainsupport.nw"
    sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' \
      "${temp_dir}/captainsupport.nw"

    gotree collapse support -s 95 \
      -i "${temp_dir}/captainsupport.nw" \
      -o "${temp_dir}/captainsupport2.nw"
    sed -i 's/)\([0-9.]*\)\/\([0-9.]*\):/)\2\/\1:/g' \
      "${temp_dir}/captainsupport2.nw"
    mv "${temp_dir}/captainsupport2.nw" "${temp_dir}/captainsupport.nw"

    gotree reroot midpoint \
      -i "${temp_dir}/captainsupport.nw" \
      -o "${working_dir}/CaptainPhylogeny.nw"

  elif (( unique_sequences >= 2 )); then
    log_warning "Only '${unique_sequences}' unique sequences found. Tree will be inferred without support values."
    mkdir -p "$outdir"
    iqtree3 -T "${threads}" -m MFP \
      --prefix "${outdir}/Captain_tree" \
      -s "$protein_file" -quiet --polytomy

    gotree collapse length -l 0.00001 \
      -i "${outdir}/Captain_tree.treefile" \
      -o "${temp_dir}/captainlength.nw"
    gotree reroot midpoint \
      -i "${temp_dir}/captainlength.nw" \
      -o "${working_dir}/CaptainPhylogeny.nw"
  else
    log_step "Only one unique sequence in dataset. Skipping tree inference."
  fi
}

# Moves captainless elements out of the active Data directories.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
remove_captainless_elements() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local data_dir
  data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

  local gff_dir nucleotide_dir protein_dir cds_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  local remove_elements="${working_dir}/Captainless_elements.txt"
  local output="${data_dir}/Captainless_elements"
  mkdir -p "$output"

  log_step "Removing captainless elements from the final dataset..."

  while read -r element; do
    mv "${gff_dir}/${element}.gff"        "${output}/"
    mv "${protein_dir}/${element}.fa"     "${output}/${element}_protein.fa"
    mv "${nucleotide_dir}/${element}.fa"  "${output}/${element}_nucleotide.fa"
    mv "${cds_dir}/${element}.fa"         "${output}/${element}_CDS.fa"
  done < "$remove_elements"
}

# Checks that cluster_stats.txt exists and counts qualifying clusters.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if file is missing
check_clusters() {
  local base_dir="$1"
  local cluster_stats="${base_dir}/Clusters/cluster_stats.txt"

  if [[ ! -f "$cluster_stats" ]]; then
    echo "Error: 'cluster_stats.txt' not found in '${base_dir}/Clusters/'." >&2
    exit 1
  fi

  cluster_number=$(awk -v min="$minimum_size" \
    'NR>1{if($2>=min){print $1}}' "$cluster_stats" | wc -l)
}

# Sets phylogeny_flag based on whether a phylogeny file was produced.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 always
check_phylogeny() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"

  if [[ -f "${working_dir}/CaptainPhylogeny.nw" ]]; then
    phylogeny_flag=true
  else
    phylogeny_flag=false
  fi
}

# Copies key result files from the workspace to the main output directory.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
organize_information() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local captain_dir="${base_dir}/$(basename -s .sh "$0")"

  mkdir -p "$captain_dir"

  if [[ ! -d "$working_dir" ]]; then
    echo "Error: working directory '${working_dir}' not found." >&2
    return 1
  fi

  [[ -f "${working_dir}/CaptainsID.txt" ]] && \
    cp "${working_dir}/CaptainsID.txt" "$captain_dir"
  [[ -f "${working_dir}/Captains_CDS.fa" ]] && \
    cp "${working_dir}/Captains_CDS.fa" "$captain_dir"
  [[ -f "${working_dir}/Captains_pseudo.fa" ]] && \
    cp "${working_dir}/Captains_pseudo.fa" "$captain_dir"

  if $phylogeny_flag; then
    cp "${working_dir}/CaptainPhylogeny.nw" "$captain_dir"
  fi
}

# Runs steps 1-3 (HMM search, captain identification, pseudogene search) for
# a given directory, and removes captainless elements if any are found.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
run_identification_steps() {
  local base_dir="$1"

  log_info "Step 1: Running hmmscan profile search."
  process_hmmscan "${base_dir}"
  log_info "-> Step 1 finished. Proceeding."

  log_info "Step 2: Identifying captains from hmmscan results."
  captain_identification "${base_dir}"
  log_info "-> Step 2 finished. Proceeding."

  log_info "Step 3: Searching for captain pseudogenes..."
  captain_pseudogene "${base_dir}"

  if [[ -s "${base_dir}/Workspace/$(basename -s .sh "$0")/Captainless_elements.txt" ]]; then
    remove_captainless_elements "${base_dir}"
    log_warning "Elements have been removed from the dataset. Please check."
  fi
  log_info "-> Step 3 finished. Proceeding."
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

working_directory=""
mode="AllID"
length="250"
level="2"
range="10"
threads="1"
minimum_size="4"
overwrite=false
previous_captain_run=false
captainless_flag=false
alignmentless_flag=false
phylogeny_flag=false
cluster_number=0

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--workingDirectory)  shift; working_directory="$1" ;;
    -m|--mode)              shift; mode="$1" ;;
    -l|--length)            shift; length="$1" ;;
    -c|--confidenceLevel)   shift; level="$1" ;;
    -r|--rangeKb)           shift; range="$1" ;;
    -t|--threads)           shift; threads="$1" ;;
    -ms|--minSize)          shift; minimum_size="$1" ;;
    --overwrite)            overwrite=true ;;
    -help)
      print_help
      exit 0
      ;;
    *)
      echo "Error: invalid option '$1'." >&2
      print_help
      exit 1
      ;;
  esac
  shift
done

# ==============================================================================
# PARAMETER SUMMARY
# ==============================================================================

echo "Running $(basename -s .sh "$0") command with the following parameters:"
echo "  Working directory:               ${working_directory}"
echo "  Mode:                            ${mode}"
echo "  Minimum captain length (aa):     ${length}"
echo "  Confidence level:                ${level}"
echo "  Range (kb):                      ${range}"
echo "  Minimum cluster size:            ${minimum_size}"
echo "  Threads:                         ${threads}"
echo "  Overwrite previous run:          ${overwrite}"
echo ""

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Checking arguments and input files..."

check_mode_parameter "$mode" "$(basename -s .sh "$0")"

if [[ -z "$working_directory" ]]; then
  echo "Error: missing required argument '-w / --workingDirectory'." >&2
  print_help
  exit 1
fi

if [[ ! -d "$working_directory" ]]; then
  echo "Error: directory '${working_directory}' does not exist." >&2
  exit 1
fi
working_directory=$(realpath "$working_directory")
check_directory_structure "${working_directory}"

if [[ "$level" =~ ^[0-9]+$ ]]; then
  if (( level < 1 || level > 3 )); then
    echo "Error: confidence level '${level}' is out of range [1-3]." >&2
    print_help; exit 1
  fi
else
  echo "Error: confidence level '${level}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$length" =~ ^[0-9]+$ ]]; then
  if (( length < 200 || length > 800 )); then
    echo "Error: length '${length}' is out of range [200-800]." >&2
    print_help; exit 1
  fi
else
  echo "Error: length '${length}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ ! "$minimum_size" =~ ^[0-9]+$ ]] || (( minimum_size < 4 || minimum_size > 10 )); then
  echo "Error: minimum size '${minimum_size}' is out of range [4-10]." >&2
  print_help; exit 1
fi

if [[ "$range" =~ ^[0-9]+$ ]]; then
  if (( range < 3 || range > 20 )); then
    echo "Error: range '${range}' is out of range [3-20]." >&2
    print_help; exit 1
  fi
else
  echo "Error: range '${range}' is not a positive integer." >&2
  print_help; exit 1
fi

check_threads "${threads}"

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

if [[ "${mode}" == "FullAll" ]]; then
  if $overwrite; then
    log_info "Removing previous run if exists..."
    overwrite "${working_directory}" "$(basename -s .sh "$0")"
    previous_captain_run=false
  else
    log_info "Checking for data from a previous run..."
    check_previous_information "${working_directory}"
  fi

  log_info "Organizing workspace..."
  organize_working_directory "${working_directory}"

  if ! $previous_captain_run; then
    run_identification_steps "${working_directory}"
  else
    log_info "Previous run data found. Steps 1, 2, and 3 will be skipped."
  fi

  log_info "Step 4: Grouping captains and performing alignment."
  run_alignment "${working_directory}"
  if $captainless_flag || $alignmentless_flag; then exit 1; fi
  log_info "-> Step 4 finished. Proceeding."

  log_info "Step 5: Running captain phylogenetic tree inference."
  run_tree_inference "${working_directory}"
  check_phylogeny "${working_directory}"
  log_info "-> Step 5 finished."

  log_info "Organizing results to main directory."
  organize_information "${working_directory}"
  rm -r "${working_directory}/Workspace/$(basename -s .sh "$0")/temp" > /dev/null
  log_info "Finished."

elif [[ "${mode}" == "AllID" ]]; then
  if $overwrite; then
    log_info "Removing previous run if exists..."
    overwrite "${working_directory}" "$(basename -s .sh "$0")"
  else
    previous_captain_run=false
  fi

  log_info "Organizing workspace..."
  organize_working_directory "${working_directory}"
  run_identification_steps "${working_directory}"

  log_info "Organizing results to main directory."
  organize_information "${working_directory}"
  rm -r "${working_directory}/Workspace/$(basename -s .sh "$0")/temp" > /dev/null
  log_info "Finished."

elif [[ "${mode}" == "Cluster" ]]; then
  check_clusters "${working_directory}"

  if $overwrite; then
    rm -f "${working_directory}/Clusters/ClustersAnalyzed.txt"
  fi

  log_info "Analyzing '${cluster_number}' clusters.\n"

  awk -v min="$minimum_size" 'NR>1{if($2>=min){print $1}}' \
    "${working_directory}/Clusters/cluster_stats.txt" \
    | sed $'s/[^[:print:]\t]//g' \
    | while read -r cluster_id; do

    local_dir="${working_directory}/Clusters/${cluster_id}"

    log_info "Analyzing cluster '${cluster_id}'."
    check_directory_structure "${local_dir}"

    if $overwrite; then
      log_info "Removing previous run if exists..."
      overwrite "${local_dir}" "$(basename -s .sh "$0")"
    else
      log_info "Checking for data from a previous run..."
      check_previous_information "${local_dir}"
    fi

    log_info "Organizing workspace..."
    organize_working_directory "${local_dir}"

    if ! $previous_captain_run; then
      run_identification_steps "${local_dir}"
    else
      log_info "Previous run data found. Steps 1, 2, and 3 will be skipped."
    fi

    log_info "Step 4: Grouping captains and performing alignment."
    captainless_flag=false
    alignmentless_flag=false
    run_alignment "${local_dir}"

    if $captainless_flag || $alignmentless_flag; then
      log_warning "There's no alignment available for this cluster."
      continue
    fi
    log_info "-> Step 4 finished. Proceeding."

    log_info "Step 5: Running captain phylogenetic tree inference."
    run_tree_inference "${local_dir}"
    check_phylogeny "${local_dir}"
    log_info "-> Step 5 finished."

    log_info "Organizing results to main directory."
    organize_information "${local_dir}"

    if ! $phylogeny_flag; then
      log_warning "No captain phylogeny file produced. This cluster requires manual inspection."
    else
      log_info "Successful run. Storing cluster '${cluster_id}' for further analysis."
      grep -w "${cluster_id}" \
        "${working_directory}/Clusters/cluster_stats.txt" \
        >> "${working_directory}/Clusters/ClustersAnalyzed.txt"
    fi

    rm -r "${local_dir}/Workspace/$(basename -s .sh "$0")/temp" > /dev/null
    log_info "Cluster '${cluster_id}' finished.\n"
  done

  log_info "All clusters have been analyzed."
fi
