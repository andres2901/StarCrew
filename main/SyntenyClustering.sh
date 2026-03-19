#!/usr/bin/env bash
# ==============================================================================
# NAME:        SyntenyClustering.sh
# DESCRIPTION: Clusters elements based on collinear regions using the syntenet
#              pipeline with DIAMOND for sequence similarity search.
#              Executes six main steps:
#              1. Data preprocessing for the syntenet pipeline.
#              2. All-vs-all DIAMOND similarity search.
#              3. Interspecies synteny detection with syntenet.
#              4. Collinearity summarization in one of four modes:
#                 Raw, SSP, FilterBlast, or FilterMetric.
#              5. Spectral clustering to define element clusters.
#              6. Organization of per-cluster output data.
# USAGE:       StarCrew SyntenyClustering [options]
#              StarCrew SyntenyClustering -help
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        12/Mar/2026
# VERSION:     1.0.0
# ==============================================================================

set -uo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

readonly SCRIPT_DIR="$( dirname -- "$( readlink -f -- "$0" )" )"
readonly LIB_DIR="${SCRIPT_DIR}/../lib"
readonly AUXILIARY_DIR="${SCRIPT_DIR}/../aux"

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

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Displays the help message with all available options.
# Arguments:
#   None
# Returns:
#   0 always
print_help() {
  echo "Command to cluster elements based on collinear regions using the syntenet pipeline."
  echo "Executes six main steps:"
  echo "  1. Data preprocessing for syntenet."
  echo "  2. All-vs-all DIAMOND similarity search."
  echo "  3. Interspecies synteny detection with syntenet."
  echo "  4. Collinearity summarization. Available modes:"
  echo "     Raw         - Pairs with >= 8% shared collinear genes."
  echo "                   WARNING: High false positive rate."
  echo "     SSP         - Strong synteny pairs only."
  echo "     FilterBlast - Raw pairs filtered by nucleotide-level BLAST."
  echo "     FilterMetric- Pairs filtered and updated by a metric system."
  echo "  5. Spectral clustering to define element clusters."
  echo "  6. Per-cluster data organization."
  echo
  echo "Usage: StarCrew $(basename -s .sh "$0") [-help] -w <directory_path>"
  echo "       [ -m <string> -a <integer> -g <integer> -e <float> -n <integer>"
  echo "         -s <integer> -th <float> { -fs <integer> -ms <integer>"
  echo "         -i <float> -c <float> } -t <integer>"
  echo "         --preCluster --captainInfo --overwrite --skip-syntenet ]"
  echo
  echo "Required args:"
  echo "  -w, --workingDirectory  Working directory where all data are stored."
  echo
  echo "Required args with defaults:"
  echo "  -m, --mode              Summarization mode (Default: FilterMetric)"
  echo "                          [Available: Raw, SSP, FilterBlast, FilterMetric]."
  echo "  -a, --anchors           Minimum anchor points for syntenet (Default: 8) [range: 3-25]."
  echo "  -g, --gaps              Maximum gaps between anchors for syntenet (Default: 8) [range: 5-25]."
  echo "  -e, --evalue            E-value threshold for syntenet (Default: 0.00001) [range: 0.00001-0.01]."
  echo "  -n, --minNodes          Minimum nodes for spectral clustering (Default: 4)."
  echo "  -s, --minSize           Minimum final sub-cluster size (Default: 1)."
  echo "  -th, --threshold        Minimum modularity score for sub-clusters (Default: 0) [range: -0.5-1.0]."
  echo
  echo "Required args with defaults in 'FilterBlast' mode:"
  echo "  -fs, --fragmentSize     Minimum BLAST fragment size (Default: 2000) [range: 1000-5000]."
  echo "  -ms, --mergeSize        Minimum merge fragment size (Default: 5000) [range: 2000-10000]."
  echo "  -i, --identity          Minimum BLAST identity % (Default: 70.0) [range: 60.0-90.0]."
  echo "  -c, --coverage          Minimum coverage of merged fragments (Default: 20.0) [range: 10.0-50.0]."
  echo
  echo "Optional args:"
  echo "  -t, --threads           Threads for DIAMOND and BLAST (Default: 8)."
  echo "  --preCluster            Pre-cluster before syntenet analysis; recommended for large datasets (Default: off)."
  echo "  --captainInfo           Include captain CDS/pseudogene info in each cluster (Default: off)."
  echo "  --overwrite             Overwrite a previous run (Default: off)."
  echo "  --skip-syntenet         Skip collinearity detection to re-run mode or clustering parameters."
  echo "                          Not compatible with '--overwrite' (Default: off)."
  echo "  -help                   Display this help message."
}

# Sets up the workspace by copying Data subdirectories and metadata if present.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if a previous workspace run exists
organize_working_directory() {
  local base_dir="$1"
  local data_dir="${base_dir}/Data"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local cluster_dir="${base_dir}/Clusters"

  local gff_dir nucleotide_dir protein_dir cds_dir metadata_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)
  metadata_dir=$(find "$base_dir" -maxdepth 1 -type d -name "metadata_files" 2>/dev/null)

  if [[ -d "$working_dir" ]]; then
    echo "Error: a previous run exists in '${working_dir}'." >&2
    echo "Use '--overwrite' to overwrite it." >&2
    exit 1
  fi

  mkdir -p "$working_dir" "$temp_dir" "$cluster_dir"
  cp -r "$gff_dir" "$working_dir"
  cp -r "$nucleotide_dir" "$working_dir"
  cp -r "$protein_dir" "$working_dir"
  cp -r "$cds_dir" "$working_dir"

  if [[ -f "${metadata_dir}/metadata.csv" ]]; then
    cp "${metadata_dir}/metadata.csv" "$working_dir"
    metadata_flag=true
  fi
}

# Runs all-vs-all DIAMOND blastp searches across all preprocessed protein files.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if PreprocessData directory is missing
process_diamond() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local protein_path
  protein_path=$(find "$working_dir" -maxdepth 1 -type d -name "PreprocessData" 2>/dev/null)
  if [[ -z "$protein_path" ]]; then
    echo "Error: PreprocessData directory not found in '${working_dir}'." >&2
    exit 1
  fi

  local dbs_dir="${working_dir}/dbs_diamond"
  local results_temp_dir="${temp_dir}/results_temp"
  local diamond_results_dir="${working_dir}/DiamondResults"
  local all_faa="${working_dir}/all.faa"

  mkdir -p "$dbs_dir" "$results_temp_dir" "$diamond_results_dir"

  log_step "Creating combined protein database file..."
  cat "${protein_path}"/*.fasta > "$all_faa"

  log_step "Building DIAMOND database..."
  local db_file="${dbs_dir}/all"
  diamond makedb --in "$all_faa" -d "$db_file" --quiet

  log_step "Performing all-vs-all similarity searches..."

  local fa_files=()
  while IFS= read -r -d '' f; do
    fa_files+=("$f")
  done < <(find "$protein_path" -maxdepth 1 -type f -name "*.fasta" -print0)

  find "$protein_path" -maxdepth 1 -type f -name "*.fasta" \
    | xargs -n1 basename -s .fasta > "${temp_dir}/all"

  local total_states=${#fa_files[@]}
  local state=0

  for fasta_file in "${fa_files[@]}"; do
    state=$(( state + 1 ))
    local starship_name
    starship_name=$(basename "$fasta_file" .fasta)
    local outfile="${results_temp_dir}/${starship_name}.tsv"

    diamond blastp -q "$fasta_file" -d "$db_file" -o "$outfile" \
      --fast -k0 --max-hsps 1 --evalue 1e-5 --matrix PAM30 \
      --query-cover 70 --subject-cover 70 -p "$threads" \
      --hit-membuf --quiet

    awk '{print $2}' "$outfile" \
      | awk -F '_' '{print $2}' \
      | awk -F '.' '{print $1}' \
      | sed 's/[^[:print:]]//g' \
      | sort | uniq > "${temp_dir}/hit"

    grep -w -f "${temp_dir}/hit" "${temp_dir}/all" | while read -r target; do
      if [[ "$starship_name" != "$target" ]]; then
        grep "_${target}\." "$outfile" \
          > "${diamond_results_dir}/${starship_name}_${target}.tsv"
        echo -e "${starship_name}\t${target}" >> "${temp_dir}/Diamond_files.txt"
      else
        grep "^.*_${target}\..*_${target}\." "$outfile" \
          > "${diamond_results_dir}/${starship_name}_${target}.tsv"
      fi
    done

    rm "${temp_dir}/hit"
    ProgressBar "$state" "$total_states"
  done
  echo ""
}

# Pre-clusters elements using DIAMOND results before running syntenet per cluster.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success
process_precluster() {
  local base_dir="$1"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"
  local collinearity_path="${working_dir}/Collinearity"
  local precluster_dir="${working_dir}/precluster"

  mkdir -p "$precluster_dir" "$collinearity_path"

  local gff_path protein_dir diamond_results_dir
  gff_path=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  protein_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  diamond_results_dir=$(find "$working_dir" -maxdepth 1 -type d -name "DiamondResults" 2>/dev/null)

  log_step "Starting pre-cluster process..."
  python "${AUXILIARY_DIR}/PreCluster.py" \
    -i "${temp_dir}/Diamond_files.txt" \
    -o "${temp_dir}" &> ${working_dir}/pre_cluster_log.txt

  local total_states
  total_states=$(wc -l < "${temp_dir}/Clusters.txt")
  local state=0

  if (( total_states >= 2 )); then
    log_step "Preparing information for '${total_states}' clusters..."

    while IFS=$'\t' read -r cluster_id values_string; do
      state=$(( state + 1 ))
      local values_array
      read -r -a values_array <<< "$values_string"

      mkdir -p "${precluster_dir}/${cluster_id}/Protein" \
               "${precluster_dir}/${cluster_id}/Gff" \
               "${precluster_dir}/${cluster_id}/DiamondResults"

      for value in "${values_array[@]}"; do
        cp "${gff_path}/${value}.gff"          "${precluster_dir}/${cluster_id}/Gff/"
        cp "${protein_dir}/${value}.fa"        "${precluster_dir}/${cluster_id}/Protein/"
        cp "${diamond_results_dir}/${value}"*  "${precluster_dir}/${cluster_id}/DiamondResults/"
      done

      ProgressBar "$state" "$total_states"
    done < "${temp_dir}/Clusters.txt"
    echo ""

    state=0
    ls -d "${precluster_dir}/"*/ | while read -r cluster_path; do
      state=$(( state + 1 ))
      log_step "Running syntenet for pre-cluster '${state}'."
      Rscript "${AUXILIARY_DIR}/syntenetAnalysis.R" \
        -d "${cluster_path}" \
        -a "${anchor_points}" -g "${gaps}" \
        -t "${threads}" -e "${evalue}"
      find "${cluster_path}/Collinearity/" -type f -name "*.collinearity" \
        -exec cp {} "${collinearity_path}/" \;
      log_step "Pre-cluster '${state}' finished."
    done
  else
    log_step "No sub-clusters found. Running syntenet on the full dataset..."
    Rscript "${AUXILIARY_DIR}/syntenetAnalysis.R" \
      -d "${working_dir}" \
      -a "${anchor_points}" -g "${gaps}" \
      -t "${threads}" -e "${evalue}"
  fi
}


# Runs BLASTN only for the element pairs listed in a collinearity file,
# grouping all subjects per query into a single per-query database to
# avoid redundant searches.
# Arguments:
#   $1 - path to the low-confidence collinearity pairs file (;-separated,
#        col1=query element, col2=subject element)
#   $2 - directory containing per-element nucleotide FASTA files
#   $3 - output file path for BLAST results
# Returns:
#   0 on success, exits with 1 if inputs are not found
run_blastn_pairwise() {
  local pairs_file="$1"
  local fasta_dir="$2"
  local output_file="$3"
  local temp_dir
  temp_dir=$(dirname "$output_file")
 
  if [[ ! -f "$pairs_file" ]]; then
    echo "Error: pairs file '${pairs_file}' not found." >&2
    exit 1
  fi
  if [[ ! -d "$fasta_dir" ]]; then
    echo "Error: FASTA directory '${fasta_dir}' not found." >&2
    exit 1
  fi
 
  local temp_id
  temp_id=$(date +%s%N)
 
  # Extract sorted unique query list (column 1)
  local query_list="${temp_dir}/temp_queries_${temp_id}.txt"
  awk -F ';' '{print $1}' "$pairs_file" | sort -u > "$query_list"
 
  local total_states
  total_states=$(wc -l < "$query_list")
  local state=0
 
  log_step "Running pairwise BLASTN for ${total_states} query elements..."
 
  while read -r query; do
    state=$(( state + 1 ))
 
    # Collect all subjects for this query from the pairs file
    local subject_list="${temp_dir}/temp_subjects_${temp_id}_${query}.txt"
    awk -F ';' -v q="$query" '$1==q {print $2}' "$pairs_file" \
      | sort -u > "$subject_list"
 
    # Build a database with only the relevant subjects
    local subject_fasta="${temp_dir}/temp_subjects_${temp_id}_${query}.fa"
    local subject_db="${temp_dir}/temp_db_${temp_id}_${query}"
 
    while read -r subject; do
      if [[ -f "${fasta_dir}/${subject}.fa" ]]; then
        cat "${fasta_dir}/${subject}.fa" >> "$subject_fasta"
      else
        log_warning "FASTA file for subject '${subject}' not found. Skipping."
      fi
    done < "$subject_list"
 
    if [[ ! -s "$subject_fasta" ]]; then
      log_warning "No valid subjects found for query '${query}'. Skipping."
      rm -f "$subject_list" "$subject_fasta"
      continue
    fi
 
    makeblastdb -dbtype nucl \
      -in "$subject_fasta" -out "$subject_db" &> /dev/null
 
    blastn \
      -query "${fasta_dir}/${query}.fa" \
      -db "$subject_db" \
      -task blastn \
      -gapopen 8 -gapextend 6 -reward 5 -penalty -4 \
      -evalue 1e-60 -num_threads "$threads" \
      -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" \
      >> "$output_file" 2> /dev/null
 
    rm -f "$subject_list" "$subject_fasta" "${subject_db}."*
    ProgressBar "$state" "$total_states"
  done < "$query_list"
  echo ""
 
  rm -f "$query_list"
}

# Summarizes syntenet collinearity results into a percentage-based report.
# Applies one of four filtering modes: Raw, SSP, FilterBlast, or FilterMetric.
# Arguments:
#   $1 - base working directory path
#   $2 - summarization mode (Raw | SSP | FilterBlast | FilterMetric)
# Returns:
#   0 on success, exits with 1 if required directories are missing
process_collinearity() {
  local base_dir="$1"
  local run_mode="$2"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"
  local temp_dir="${working_dir}/temp"

  local collinearity_path gff_path nucleotide_dir diamond_results_dir
  collinearity_path=$(find "$working_dir" -maxdepth 1 -type d -name "Collinearity" 2>/dev/null)
  gff_path=$(find "$working_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  diamond_results_dir=$(find "$working_dir" -maxdepth 1 -type d -name "DiamondResults" 2>/dev/null)

  if [[ -z "$collinearity_path" || -z "$gff_path" \
     || -z "$nucleotide_dir" || -z "$diamond_results_dir" ]]; then
    echo "Error: required subdirectories not found in '${working_dir}'." >&2
    exit 1
  fi

  mkdir -p "$temp_dir"

  local metadata_file=""
  if $metadata_flag; then
    metadata_file="${working_dir}/metadata.csv"
  fi

  local temp_id
  temp_id=$(date +%s)
  local temp_prefix="${temp_dir}/temp_${temp_id}"

  # Count genes per element
  log_step "Counting genes per element from GFF files..."
  find "${gff_path}" -maxdepth 1 -type f -name "*.gff" \
    -exec grep -c "gene" {} + \
    | awk -F'/' '{ gsub(".gff:", "\t", $NF); print $NF }' \
    | sort -k1,1 > "${temp_prefix}_genes_per_element.txt"

  # Count collinear genes per element pair
  log_step "Counting collinear genes from collinearity files..."
  find "${collinearity_path}" -name '*.collinearity' \
    -exec grep -E "[0-9]*-.*[0-9]*:" {} + \
    > "${temp_prefix}_raw_collinear_lines.txt"

  awk -F '/' '{print $NF}' "${temp_prefix}_raw_collinear_lines.txt" \
    | sed -e 's/ //g; s/\.collinearity//g; s/:/\t/g' \
    | awk '
      BEGIN { OFS="\t" }
      {
        key3 = $1 SUBSEP $3
        if (!(key3 in count3)) { count3[key3]=1; elements3[$1]++ }
        key4 = $1 SUBSEP $4
        if (!(key4 in count4)) { count4[key4]=1; elements4[$1]++ }
      }
      END {
        for (k in elements3) print k, elements3[k], elements4[k]
      }' \
    | sed -e 's/_/\t/g' | sort -k1,1 > "${temp_prefix}_collinear_genes.txt"

  # Join counts with gene totals
  log_step "Joining gene count data..."
  join "${temp_prefix}_collinear_genes.txt" \
       "${temp_prefix}_genes_per_element.txt" \
    | sort -k2,2 > "${temp_prefix}_join_1.txt"
  join -1 2 -2 1 "${temp_prefix}_join_1.txt" \
       "${temp_prefix}_genes_per_element.txt" \
    | awk '{swap=$1;$1=$2;$2=swap;print $0}' \
    | sort -k1,1 > "${temp_prefix}_join_2.txt"

  # Estimate per-element collinearity percentage
  log_step "Estimating collinearity percentages per element..."
  awk '{print $0"\t"100*($3/$5)"\t"100*($4/$6)}' \
    "${temp_prefix}_join_2.txt" > "${temp_prefix}_percentage_pairwise.txt"

  # Extract software-reported general percentages
  log_step "Extracting software-reported general percentages..."
  find "${collinearity_path}" -name '*.collinearity' \
    -exec grep -w -E -o "Percentage: [0-9]*\.[0-9]*" {} + \
    > "${temp_prefix}_raw_percentage_general.txt"
  awk -F '/' '{print $NF}' "${temp_prefix}_raw_percentage_general.txt" \
    | sed -e 's/ //g; s/.collinearity//g; s/:/\t/g; s/_/\t/g' \
    > "${temp_prefix}_percentage_general.txt"

  # Filter and initialize final output
  echo -e "element01;element02;General_percentage;element01_percentage;element02_percentage" \
    > "${temp_dir}/Collinearity_percentage.txt"
  awk '{if($4>0){print}}' "${temp_prefix}_percentage_general.txt" \
    > "${temp_prefix}_percentage_general_filter.txt"

  # Shared awk join pattern used across modes
  local join_awk='NR==FNR{f1[$1,$2]=$0;next} $1 SUBSEP $2 in f1{print f1[$1,$2],"\t"$7,"\t"$8}'
  local format_awk='{print $1 FS $2 FS $4 FS $5 FS $6}'

  log_step "Generating final report in '${run_mode}' mode..."

  if [[ "$run_mode" == "Raw" ]]; then
    awk "$join_awk" \
        "${temp_prefix}_percentage_general_filter.txt" \
        "${temp_prefix}_percentage_pairwise.txt" \
      | awk "$format_awk" \
      | awk '{if($3 >= 8) print}' \
      | sed -e 's/ /;/g' | sort -t ';' \
      >> "${temp_dir}/Collinearity_percentage.txt"

    if $metadata_flag; then
      python "${AUXILIARY_DIR}/merge_metadata.py" \
        -d "${temp_dir}/Collinearity_percentage.txt" \
        -m "$metadata_file" \
        -o "${working_dir}/Collinearity_percentage.txt" &> ${working_dir}/merge_metadata_log.txt
    else
      sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" \
        >> "${working_dir}/Collinearity_percentage.txt"
    fi

  elif [[ "$run_mode" == "SSP" ]]; then
    awk "$join_awk" \
        "${temp_prefix}_percentage_general_filter.txt" \
        "${temp_prefix}_percentage_pairwise.txt" \
      | awk "$format_awk" \
      | awk '{if(($3>=41)||(($4>=45||$5>=45)&&($4/$5>=1.8||$4/$5<=0.55))) print}' \
      | sed -e 's/ /;/g' | sort -t ';' \
      >> "${temp_dir}/Collinearity_percentage.txt"

    if $metadata_flag; then
      python "${AUXILIARY_DIR}/merge_metadata.py" \
        -d "${temp_dir}/Collinearity_percentage.txt" \
        -m "$metadata_file" \
        -o "${working_dir}/Collinearity_percentage.txt" &> ${working_dir}/merge_metadata_log.txt
    else
      sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" \
        >> "${working_dir}/Collinearity_percentage.txt"
    fi

  elif [[ "$run_mode" == "FilterBlast" ]]; then
    awk "$join_awk" \
        "${temp_prefix}_percentage_general_filter.txt" \
        "${temp_prefix}_percentage_pairwise.txt" \
      | awk "$format_awk" \
      | awk '{if($3 >= 8) print}' \
      | sed -e 's/ /;/g' | sort -t ';' \
      >> "${temp_prefix}_Collinearity_percentage.txt"

    awk 'BEGIN{FS=";";OFS=";"}{
      if(($3>=41)||(($4>=45||$5>=45)&&($4/$5>=1.8||$4/$5<=0.55))) print
    }' "${temp_prefix}_Collinearity_percentage.txt" \
      > "${temp_prefix}_Collinearity_percentage_high.txt"

    awk 'BEGIN{FS=";";OFS=";"}{if($3 < 41) print}' \
      "${temp_prefix}_Collinearity_percentage.txt" \
      > "${temp_prefix}_Collinearity_percentage_low.txt"

    if [[ -s "${temp_prefix}_Collinearity_percentage_low.txt" ]]; then
      run_blastn_pairwise \
        "${temp_prefix}_Collinearity_percentage_low.txt" \
        "${nucleotide_dir}" \
        "${working_dir}/BlastnResults.out"

      log_step "Filtering BLASTN results..."
      python "${AUXILIARY_DIR}/Blast_CleanUp.py" \
        -f "${working_dir}/BlastnResults.out" \
        -o "${working_dir}/BlastnClean.out" \
        -fs "${fragment_size}" -i "${identity}" \
        -ms "${merge_size}" -c "${coverage}" &> ${working_dir}/blast_cleanup_log.txt

      log_step "Filtering false positive pairs..."
      awk '{OFS=";"} {print $1 OFS $2}' "${working_dir}/BlastnClean.out" \
        > "${working_dir}/BlastnPairs.out"
      grep -f "${working_dir}/BlastnPairs.out" \
        "${temp_prefix}_Collinearity_percentage_low.txt" > "${temp_prefix}_Collinearity_percentage_lowSelected.txt"
      cat "${temp_prefix}_Collinearity_percentage_high.txt" \
          "${temp_prefix}_Collinearity_percentage_lowSelected.txt" \
        | sort -u > "${temp_prefix}_Collinearity_percentage_filter.txt"
    else
      log_step "No low-synteny pairs to evaluate. Skipping BLAST filtering."
      mv "${temp_prefix}_Collinearity_percentage_high.txt" \
         "${temp_prefix}_Collinearity_percentage_filter.txt"
    fi

    cat "${temp_prefix}_Collinearity_percentage_filter.txt" \
      >> "${temp_dir}/Collinearity_percentage.txt"
    if $metadata_flag; then
      python "${AUXILIARY_DIR}/merge_metadata.py" \
        -d "${temp_dir}/Collinearity_percentage.txt" \
        -m "$metadata_file" \
        -o "${working_dir}/Collinearity_percentage.txt" &> ${working_dir}/merge_metadata_log.txt
    else
      sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" \
        >> "${working_dir}/Collinearity_percentage.txt"
    fi

  elif [[ "$run_mode" == "FilterMetric" ]]; then
    local metric_threshold
    metric_threshold=$(python -c "print(max(6, ${anchor_points}))")

    awk "$join_awk" \
        "${temp_prefix}_percentage_general_filter.txt" \
        "${temp_prefix}_percentage_pairwise.txt" \
      | awk "$format_awk" \
      | awk '{if($3 >= 8) print}' \
      | sed -e 's/ /;/g' | sort -t ';' \
      >> "${temp_prefix}_Collinearity_percentage.txt"

    log_step "Filtering results using metric threshold '${metric_threshold}'..."
    python "${AUXILIARY_DIR}/FilterMetric.py" \
      --collinearity "${temp_prefix}_Collinearity_percentage.txt" \
      --syntenet "${collinearity_path}/" \
      --diamond "${diamond_results_dir}/" \
      --threshold "${metric_threshold}" \
      --output "${temp_dir}/Collinearity_percentage.txt" &> ${working_dir}/filter_metric_log.txt

    if $metadata_flag; then
      python "${AUXILIARY_DIR}/merge_metadata.py" \
        -d "${temp_dir}/Collinearity_percentage.txt" \
        -m "$metadata_file" \
        -o "${working_dir}/Collinearity_percentage.txt" &> ${working_dir}/merge_metadata_log.txt
    else
      sed -e 's/;/\t/g' "${temp_dir}/Collinearity_percentage.txt" \
        >> "${working_dir}/Collinearity_percentage.txt"
    fi
  fi

  rm -rf "${temp_prefix}"*
}

# Organizes per-cluster output directories and copies relevant data files.
# Optionally includes captain CDS/pseudogene information per element.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 if cluster file is missing
process_cluster_file() {
  local base_dir="$1"
  local data_dir="${base_dir}/Data"
  local cluster_dir="${base_dir}/Clusters"
  local working_dir="${base_dir}/Workspace/$(basename -s .sh "$0")"

  local gff_dir nucleotide_dir protein_dir cds_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  local metadata_file=""
  if $metadata_flag; then
    metadata_file="${working_dir}/metadata.csv"
  fi

  local cds_file="" id_file="" pseudo_file=""
  if $captain_info; then
    [[ -f "${captain_dir}/Captains_CDS.fa" ]] && cds_file="${captain_dir}/Captains_CDS.fa"
    [[ -f "${captain_dir}/CaptainsID.txt" ]]  && id_file="${captain_dir}/CaptainsID.txt"
    [[ -f "${captain_dir}/Captains_pseudo.fa" ]] && pseudo_file="${captain_dir}/Captains_pseudo.fa"
  fi

  local cluster_path="${cluster_dir}/main_clusters.txt"
  if [[ ! -f "$cluster_path" ]]; then
    echo "Error: file '${cluster_path}' not found." >&2
    return 1
  fi

  local total_states
  total_states=$(wc -l < "$cluster_path")
  local state=0

  while IFS=$'\t' read -r cluster_id values_string; do
    state=$(( state + 1 ))
    local values_array
    read -r -a values_array <<< "$values_string"

    mkdir -p "${cluster_dir}/${cluster_id}/Workspace" \
             "${cluster_dir}/${cluster_id}/Data/CDS" \
             "${cluster_dir}/${cluster_id}/Data/Nucleotide" \
             "${cluster_dir}/${cluster_id}/Data/Protein" \
             "${cluster_dir}/${cluster_id}/Data/Gff"

    if $captain_info; then
      mkdir -p "${cluster_dir}/${cluster_id}/CaptainIdentification"
    fi

    local updated_metadata=""
    if $metadata_flag; then
      mkdir -p "${cluster_dir}/${cluster_id}/metadata_files"
      updated_metadata="${cluster_dir}/${cluster_id}/metadata_files/metadata.csv"
      head -n1 "$metadata_file" > "$updated_metadata"
    fi

    for value in "${values_array[@]}"; do
      cp "${gff_dir}/${value}.gff"       "${cluster_dir}/${cluster_id}/Data/Gff/"
      cp "${protein_dir}/${value}.fa"    "${cluster_dir}/${cluster_id}/Data/Protein/"
      cp "${nucleotide_dir}/${value}.fa" "${cluster_dir}/${cluster_id}/Data/Nucleotide/"
      cp "${cds_dir}/${value}.fa"        "${cluster_dir}/${cluster_id}/Data/CDS/"

      if $metadata_flag; then
        grep -w "$value" "$metadata_file" >> "$updated_metadata"
      fi

      if $captain_info; then
        if [[ -n "$cds_file" ]]; then
          local cds_value
          cds_value=$(grep "${value}\." "$cds_file" | awk -F '>' '{print $2}')
          if [[ -n "$cds_value" ]]; then
            seqkit grep --quiet -p "$cds_value" \
              "${cluster_dir}/${cluster_id}/Data/CDS/${value}.fa" \
              >> "${cluster_dir}/${cluster_id}/CaptainIdentification/Captains_CDS.fa"
            grep "^${value}\." "$id_file" \
              >> "${cluster_dir}/${cluster_id}/CaptainIdentification/CaptainsID.txt"
          else
            seqkit grep --quiet -p "$value" "$pseudo_file" \
              >> "${cluster_dir}/${cluster_id}/CaptainIdentification/Captains_pseudo.fa"
          fi
        else
          seqkit grep --quiet -p "$value" "$pseudo_file" \
            >> "${cluster_dir}/${cluster_id}/CaptainIdentification/Captains_pseudo.fa"
        fi
      fi
    done

    ProgressBar "$state" "$total_states"
  done < "$cluster_path"
  echo ""
  log_step "Cluster file processing complete."
}

# ==============================================================================
# DEFAULT VARIABLE VALUES
# ==============================================================================

working_directory=""
mode="FilterMetric"
anchor_points="8"
gaps="8"
evalue="0.00001"
min_nodes="4"
min_size="2"
threshold="0.02"
fragment_size="2000"
merge_size="5000"
identity="70.0"
coverage="20.0"
threads="8"
precluster=false
captain_info=false
overwrite=false
skip_syntenet=false
metadata_flag=false
captain_dir=""

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--workingDirectory)  shift; working_directory="$1" ;;
    -m|--mode)              shift; mode="$1" ;;
    -a|--anchors)           shift; anchor_points="$1" ;;
    -g|--gaps)              shift; gaps="$1" ;;
    -e|--evalue)            shift; evalue="$1" ;;
    -n|--minNodes)          shift; min_nodes="$1" ;;
    -s|--minSize)           shift; min_size="$1" ;;
    -th|--threshold)        shift; threshold="$1" ;;
    -fs|--fragmentSize)     shift; fragment_size="$1" ;;
    -ms|--mergeSize)        shift; merge_size="$1" ;;
    -i|--identity)          shift; identity="$1" ;;
    -c|--coverage)          shift; coverage="$1" ;;
    -t|--threads)           shift; threads="$1" ;;
    --preCluster)           precluster=true ;;
    --captainInfo)          captain_info=true ;;
    --overwrite)            overwrite=true ;;
    --skip-syntenet)        skip_syntenet=true ;;
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
echo "  Minimum anchor points:           ${anchor_points}"
echo "  Maximum gaps:                    ${gaps}"
echo "  E-value threshold:               ${evalue}"
echo "  Minimum nodes (spectral):        ${min_nodes}"
echo "  Minimum sub-cluster size:        ${min_size}"
echo "  Modularity threshold:            ${threshold}"
echo "  Threads:                         ${threads}"
echo "  Pre-clustering:                  ${precluster}"
echo "  Captain information:             ${captain_info}"
echo "  Overwrite previous run:          ${overwrite}"
echo "  Skip syntenet:                   ${skip_syntenet}"
echo ""

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Checking arguments and input files..."

if [[ -z "$working_directory" ]]; then
  echo "Error: missing required argument '-w / --workingDirectory'." >&2
  print_help; exit 1
fi

if [[ ! -d "$working_directory" ]]; then
  echo "Error: directory '${working_directory}' does not exist." >&2
  exit 1
fi
working_directory=$(realpath "$working_directory")
check_directory_structure "${working_directory}"

check_mode_parameter "$mode" "$(basename -s .sh "$0")"

if [[ "$anchor_points" =~ ^[0-9]+$ ]]; then
  if (( anchor_points < 3 || anchor_points > 25 )); then
    echo "Error: anchor points '${anchor_points}' is out of range [3-25]." >&2
    print_help; exit 1
  fi
else
  echo "Error: anchor points '${anchor_points}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$gaps" =~ ^[0-9]+$ ]]; then
  if (( gaps < 5 || gaps > 25 )); then
    echo "Error: gaps '${gaps}' is out of range [5-25]." >&2
    print_help; exit 1
  fi
else
  echo "Error: gaps '${gaps}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$evalue" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
  if (( $(echo "$evalue < 0.00001" | bc -l) )) || \
     (( $(echo "$evalue > 0.01" | bc -l) )); then
    echo "Error: e-value '${evalue}' is out of range [0.00001-0.01]." >&2
    print_help; exit 1
  fi
else
  echo "Error: e-value '${evalue}' is not a valid float." >&2
  print_help; exit 1
fi

if [[ ! "$min_nodes" =~ ^[0-9]+$ ]]; then
  echo "Error: minimum nodes '${min_nodes}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ ! "$min_size" =~ ^[0-9]+$ ]]; then
  echo "Error: minimum size '${min_size}' is not a positive integer." >&2
  print_help; exit 1
fi

if [[ "$threshold" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
  if (( $(echo "$threshold < -0.5" | bc -l) )) || \
     (( $(echo "$threshold > 1.0" | bc -l) )); then
    echo "Error: threshold '${threshold}' is out of range [-0.5-1.0]." >&2
    print_help; exit 1
  fi
else
  echo "Error: threshold '${threshold}' is not a valid float." >&2
  print_help; exit 1
fi

if [[ "$mode" == "FilterBlast" ]]; then
  if [[ "$fragment_size" =~ ^[0-9]+$ ]]; then
    if (( fragment_size < 1000 || fragment_size > 5000 )); then
      echo "Error: fragment size '${fragment_size}' is out of range [1000-5000]." >&2
      print_help; exit 1
    fi
  else
    echo "Error: fragment size '${fragment_size}' is not a positive integer." >&2
    print_help; exit 1
  fi

  if [[ "$merge_size" =~ ^[0-9]+$ ]]; then
    if (( merge_size < 2000 || merge_size > 10000 )); then
      echo "Error: merge size '${merge_size}' is out of range [2000-10000]." >&2
      print_help; exit 1
    fi
  else
    echo "Error: merge size '${merge_size}' is not a positive integer." >&2
    print_help; exit 1
  fi

  if (( fragment_size > merge_size )); then
    echo "Error: fragment size cannot be greater than merge size." >&2
    print_help; exit 1
  fi

  if [[ "$identity" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    if (( $(echo "$identity < 60.0" | bc -l) )) || \
       (( $(echo "$identity > 90.0" | bc -l) )); then
      echo "Error: identity '${identity}' is out of range [60.0-90.0]." >&2
      print_help; exit 1
    fi
  else
    echo "Error: identity '${identity}' is not a valid float." >&2
    print_help; exit 1
  fi

  if [[ "$coverage" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    if (( $(echo "$coverage < 10.0" | bc -l) )) || \
       (( $(echo "$coverage > 50.0" | bc -l) )); then
      echo "Error: coverage '${coverage}' is out of range [10.0-50.0]." >&2
      print_help; exit 1
    fi
  else
    echo "Error: coverage '${coverage}' is not a valid float." >&2
    print_help; exit 1
  fi
fi

check_threads "$threads"

if $captain_info; then
  if [[ ! -d "${working_directory}/CaptainIdentification" ]]; then
    log_warning "No CaptainIdentification folder found. Proceeding without captain information."
    captain_info=false
  else
    captain_dir=$(realpath "${working_directory}/CaptainIdentification")
    local_captain_count=$(ls "${captain_dir}/Captain"* 2>/dev/null | wc -l)
    if (( local_captain_count == 0 )); then
      log_warning "'${captain_dir}' is empty. Proceeding without captain information."
      captain_info=false
    else
      while read -r captain_fa; do
        check_fasta_dna "$captain_fa"
      done < <(ls "${captain_dir}/Captain"*.fa 2>/dev/null)
    fi
  fi
fi

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

if $overwrite; then
  log_info "Removing previous run if exists..."
  overwrite "${working_directory}" "$(basename -s .sh "$0")"
fi

if ! $skip_syntenet; then
  log_info "Organizing workspace..."
  organize_working_directory "${working_directory}"
  log_step "Metadata: ${metadata_flag}"

  log_info "Step 1: Preprocessing data for DIAMOND."
  Rscript "${AUXILIARY_DIR}/syntenetPreprocess.R" \
    "${working_directory}/Workspace/$(basename -s .sh "$0")/"
  log_info "-> Step 1 finished. Proceeding."

  log_info "Step 2: Running DIAMOND analysis."
  process_diamond "${working_directory}"
  log_info "-> Step 2 finished. Proceeding."

  log_info "Step 3: Running syntenet analysis."
  if $precluster; then
    process_precluster "${working_directory}"
  else
    Rscript "${AUXILIARY_DIR}/syntenetAnalysis.R" \
      -d "${working_directory}/Workspace/$(basename -s .sh "$0")/" \
      -a "${anchor_points}" -g "${gaps}" \
      -t "${threads}" -e "${evalue}"
  fi
  log_info "-> Step 3 finished. Proceeding."

else
  log_info "Skipping syntenet analysis."
  log_info "Removing previous cluster results if available..."
  rm -r "${working_directory}/Clusters/"* 2> /dev/null
  if [[ -f "${working_directory}/metadata_files/metadata.csv" ]]; then
    metadata_flag=true
  fi
  log_step "Metadata: ${metadata_flag}"
fi

log_info "Step 4: Summarizing syntenet results in '${mode}' mode."
process_collinearity "${working_directory}" "${mode}"
log_info "-> Step 4 finished. Proceeding."

log_info "Step 5: Generating element clusters."
python "${AUXILIARY_DIR}/Clustering.py" \
  -i "${working_directory}/Workspace/$(basename -s .sh "$0")/Collinearity_percentage.txt" \
  -o "${working_directory}/Clusters/" \
  -m "${min_size}" -n "${min_nodes}" -t "${threshold}" &> ${working_directory}/Workspace/$(basename -s .sh "$0")/clustering_log.txt
log_info "-> Step 5 finished. Proceeding."

log_info "Step 6: Sorting elements into clusters."
process_cluster_file "${working_directory}"
rm -r "${working_directory}/Workspace/$(basename -s .sh "$0")/temp" 2> /dev/null
log_info "Synteny and clustering analysis finished."
