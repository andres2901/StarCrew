#!/usr/bin/env bash
# ==============================================================================
# NAME:        Check.sh
# DESCRIPTION: Shared validation library for StarCrew pipeline scripts.
#              Provides functions to check software availability, directory
#              structure, file integrity, and database completeness.
# USAGE:       source Check.sh
# AUTHOR:      Andres F. Lizcano Salas
# DATE:        12/Mar/2026
# VERSION:     1.0.0
# ==============================================================================

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Verifies that the StarCrew main directory contains all required subdirectories,
# scripts, HMM profiles, and database files.
# Arguments:
#   $1 - path to the StarCrew main installation directory
# Returns:
#   0 on success, exits with 1 on any missing component
check_main_directory() {
  local main_dir="$1"

  # bin/
  if [[ ! -d "${main_dir}/bin" ]]; then
    echo "Error: directory '${main_dir}/bin' does not exist." >&2
    exit 1
  fi
  if [[ ! -f "${main_dir}/bin/StarCrew" ]]; then
    echo "Error: StarCrew script not found in '${main_dir}/bin'." >&2
    exit 1
  fi

  # main/
  if [[ ! -d "${main_dir}/main" ]]; then
    echo "Error: directory '${main_dir}/main' does not exist." >&2
    exit 1
  fi
  local main_scripts=(
    CaptainIdentification
    ClusterCharacterization
    Initialize
    OrthogroupsAnnotation
    OrthogroupsOverrepresentation
    SyntenyClustering
  )
  for script in "${main_scripts[@]}"; do
    if [[ ! -f "${main_dir}/main/${script}.sh" ]]; then
      echo "Error: missing command script '${script}.sh' in '${main_dir}/main'." >&2
      exit 1
    fi
  done

  # aux/
  if [[ ! -d "${main_dir}/aux" ]]; then
    echo "Error: directory '${main_dir}/aux' does not exist." >&2
    exit 1
  fi
  check_auxiliary_scripts "${main_dir}/aux" "All"

  # hmm/
  if [[ ! -d "${main_dir}/hmm" ]]; then
    echo "Error: directory '${main_dir}/hmm' does not exist." >&2
    exit 1
  fi
  local hmm_path="${main_dir}/hmm"
  local domains=("CAT_domain" "DUF3435" "Captain")
  local extensions=("" ".h3f" ".h3i" ".h3m" ".h3p")
  for dom in "${domains[@]}"; do
    for ext in "${extensions[@]}"; do
      if [[ ! -f "${hmm_path}/${dom}.hmm${ext}" ]]; then
        echo "Error: missing HMM file '${dom}.hmm${ext}' in '${hmm_path}'." >&2
        exit 1
      fi
    done
  done

  # databases/
  if [[ ! -d "${main_dir}/databases" ]]; then
    echo "Error: directory '${main_dir}/databases' does not exist." >&2
    exit 1
  fi
  local db_path="${main_dir}/databases"
  if [[ ! -f "${db_path}/Captains.fa" ]]; then
    echo "Error: 'Captains.fa' not found in '${db_path}'." >&2
    exit 1
  fi
  check_fasta_protein "${db_path}/Captains.fa"

  if [[ ! -f "${db_path}/Captains_CDS.fa" ]]; then
    echo "Error: 'Captains_CDS.fa' not found in '${db_path}'." >&2
    exit 1
  fi
  check_fasta_dna "${db_path}/Captains_CDS.fa"

  # Required root files
  local root_files=("agat_config.yaml" "json_updater.py" "StarCrew_environment.yml")
  for f in "${root_files[@]}"; do
    if [[ ! -f "${main_dir}/${f}" ]]; then
      echo "Error: '${f}' not found in '${main_dir}'." >&2
      exit 1
    fi
  done
}

# Validates that the number of threads is a positive integer.
# Arguments:
#   $1 - threads value to validate
# Returns:
#   0 on success, exits with 1 if invalid
check_threads() {
  local threads="$1"
  if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: threads '${threads}' is not a positive integer." >&2
    print_help
    exit 1
  fi
}

# Validates that a file exists and is a valid DNA FASTA file.
# Arguments:
#   $1 - path to the FASTA file
# Returns:
#   0 on success, exits with 1 if missing or invalid
check_fasta_dna() {
  local fasta_path="$1"
  if [[ ! -f "$fasta_path" ]]; then
    echo "Error: file '${fasta_path}' does not exist." >&2
    exit 1
  fi
  local fasta_check
  fasta_check=$(seqkit seq --quiet --seq-type dna -v "$fasta_path")
  if [[ -z "$fasta_check" ]]; then
    echo "Error: '${fasta_path}' is not a valid DNA FASTA file." >&2
    exit 1
  fi
}

# Validates that a file exists and is a valid protein FASTA file.
# Arguments:
#   $1 - path to the FASTA file
# Returns:
#   0 on success, exits with 1 if missing or invalid
check_fasta_protein() {
  local fasta_path="$1"
  if [[ ! -f "$fasta_path" ]]; then
    echo "Error: file '${fasta_path}' does not exist." >&2
    exit 1
  fi
  local fasta_check
  fasta_check=$(seqkit seq --quiet -t protein -v "$fasta_path")
  if [[ -z "$fasta_check" ]]; then
    echo "Error: '${fasta_path}' is not a valid protein FASTA file." >&2
    exit 1
  fi
}

# Validates that a GFF3 file exists and has consistent 9-column structure,
# valid strand values, and ID attributes on all features.
# Arguments:
#   $1 - path to the GFF file
# Returns:
#   0 on success, exits with 1 if missing or malformed
check_gff_file() {
  local gff_path="$1"
  if [[ ! -f "$gff_path" ]]; then
    echo "Error: file '${gff_path}' does not exist." >&2
    exit 1
  fi

  local col_check
  col_check=$(grep -v "^#" "$gff_path" | awk -F '\t' '{print NF}' | sort -u)
  if [[ "$col_check" -ne 9 ]]; then
    echo "Error: '${gff_path}' has inconsistent column count (expected 9)." >&2
    exit 1
  fi

  local strand_check
  strand_check=$(grep -v "^#" "$gff_path" \
    | awk -F '\t' '{print $7}' | grep -Evc "\+|-|\.")
  if [[ "$strand_check" -ne 0 ]]; then
    echo "Error: '${gff_path}' has invalid values in the strand column." >&2
    exit 1
  fi

  local id_check
  id_check=$(grep -v "^#" "$gff_path" | grep -vc "ID=")
  if [[ "$id_check" -ne 0 ]]; then
    echo "Error: '${gff_path}' has lines missing the ID attribute." >&2
    exit 1
  fi
}

# Validates a two-column file mapping element IDs to GFF file paths, and checks
# that all referenced GFF files exist on disk.
# Arguments:
#   $1 - path to the ome2gff mapping file
# Returns:
#   0 on success, exits with 1 if missing, malformed, or files not found
check_gff_paths() {
  local gff_path="$1"
  if [[ ! -f "$gff_path" ]]; then
    echo "Error: file '${gff_path}' does not exist." >&2
    exit 1
  fi

  local col_check
  col_check=$(awk -F '\t' '{print NF}' "$gff_path" | sort -u)
  if [[ "$col_check" -ne 2 ]]; then
    echo "Error: '${gff_path}' does not have exactly 2 columns." >&2
    exit 1
  fi

  while read -r element gff_file; do
    if [[ ! -f "$gff_file" ]]; then
      echo "Error: GFF file '${gff_file}' for element '${element}' not found." >&2
      exit 1
    fi
  done < "$gff_path"
}

# Validates that a boundaries file exists and has 21 columns with valid strand
# and coordinate values.
# Arguments:
#   $1 - path to the boundaries file
# Returns:
#   0 on success, exits with 1 if missing or malformed
check_boundaries_file() {
  local boundaries_path="$1"
  if [[ ! -f "$boundaries_path" ]]; then
    echo "Error: file '${boundaries_path}' does not exist." >&2
    exit 1
  fi

  local col_check
  col_check=$(grep -v "^#" "$boundaries_path" \
    | awk -F '\t' '{print NF}' | sort -u)
  if [[ "$col_check" -ne 21 ]]; then
    echo "Error: '${boundaries_path}' does not have 21 columns." >&2
    exit 1
  fi

  local strand_check start_check end_check
  strand_check=$(grep -v "^#" "$boundaries_path" \
    | awk -F '\t' '{print $7}' | grep -Evc "\+|-")
  start_check=$(grep -v "^#" "$boundaries_path" \
    | awk -F '\t' '{print $4}' | grep -Ewvc "[0-9]*")
  end_check=$(grep -v "^#" "$boundaries_path" \
    | awk -F '\t' '{print $5}' | grep -Ewvc "[0-9]*")

  if [[ "$strand_check" -ne 0 || "$start_check" -ne 0 || "$end_check" -ne 0 ]]; then
    echo "Error: '${boundaries_path}' has inconsistent column values." >&2
    exit 1
  fi
}

# Validates a metadata CSV file: checks existence, column consistency, unique
# headers, and that at least some elements match the provided FASTA headers.
# Arguments:
#   $1 - path to the metadata CSV file
#   $2 - path to a FASTA file whose headers are used for cross-validation
# Returns:
#   0 on success, exits with 1 if missing or malformed
check_metadata_file() {
  local metadata_path="$1"
  local fasta_path="$2"

  if [[ ! -f "$metadata_path" ]]; then
    echo "Error: file '${metadata_path}' does not exist." >&2
    exit 1
  fi

  if [[ ! -s "$metadata_path" ]]; then
    return 0
  fi

  if [[ $(head -n1 "$metadata_path" | grep -c "ElementID") -ne 1 ]]; then
    return 0
  fi

  local col_counts col_unique
  col_counts=$(awk -F ';' '{print NF}' "$metadata_path" | sort -u)
  col_unique=$(echo "$col_counts" | wc -l)

  if (( col_unique >= 2 )); then
    echo "Error: inconsistent column count across rows in '${metadata_path}'." >&2
    exit 1
  fi

  local col_number
  col_number=$(echo "$col_counts")
  if (( col_number == 1 )); then
    echo "Error: only one column found in '${metadata_path}'." >&2
    exit 1
  fi

  local dup_check
  dup_check=$(head -n1 "$metadata_path" \
    | tr ';' '\n' | sort | uniq -c | awk '{print $1}' | sort -u)
  if [[ "$dup_check" -ne 1 ]]; then
    echo "Error: duplicate header columns in '${metadata_path}'." >&2
    exit 1
  fi

  local headers_pattern metadata_count
  headers_pattern=$(grep "^>" "$fasta_path" \
    | awk -F '>' '{print $2}' | tr '\n' '|' | sed 's/|$//')
  metadata_count=$(grep -Ec "$headers_pattern" "$metadata_path")
  if (( metadata_count == 0 )); then
    echo "Error: metadata does not contain any elements from the FASTA file." >&2
    exit 1
  fi
}

# Validates that all required auxiliary Python and R scripts exist for the
# given command.
# Arguments:
#   $1 - path to the auxiliary scripts directory
#   $2 - command name (Initialize | SyntenyClustering | CaptainIdentification |
#        ClusterCharacterization | OrthogroupsOverrepresentation | All)
# Returns:
#   0 on success, exits with 1 if any script is missing
check_auxiliary_scripts() {
  local auxiliary_path="$1"
  local command="$2"

  local scripts=()
  case "$command" in
    Initialize)
      scripts=("rip_calculator.py" "gff_slicer.py" "gff_split.py" "merge.py")
      ;;
    SyntenyClustering)
      scripts=(
        "PreCluster.py" "Blast_CleanUp.py" "Clustering.py"
        "merge_metadata.py" "FilterMetric.py"
        "syntenetAnalysis.R" "syntenetPreprocess.R"
      )
      ;;
    CaptainIdentification)
      scripts=("hmmscan_process.py")
      ;;
    ClusterCharacterization)
      scripts=("ClusterAnalysis.R")
      ;;
    OrthogroupsOverrepresentation)
      scripts=("OverrepresentationAnalysis.R")
      ;;
    OrthogroupsAnnotation)
      scripts=()
      ;;
    All)
      scripts=(
        "rip_calculator.py" "gff_slicer.py" "gff_split.py" "merge.py"
        "PreCluster.py" "Blast_CleanUp.py" "Clustering.py"
        "merge_metadata.py" "FilterMetric.py"
        "syntenetAnalysis.R" "syntenetPreprocess.R"
        "hmmscan_process.py" "ClusterAnalysis.R"
        "OverrepresentationAnalysis.R"
      )
      ;;
    *)
      echo "Error: unknown command '${command}' in check_auxiliary_scripts." >&2
      exit 1
      ;;
  esac

  for script in "${scripts[@]}"; do
    if [[ ! -f "${auxiliary_path}/${script}" ]]; then
      echo "Error: missing auxiliary script '${script}' in '${auxiliary_path}'." >&2
      exit 1
    fi
  done
}

# Validates that the provided mode is accepted for the given command.
# Arguments:
#   $1 - mode value to validate
#   $2 - command name
# Returns:
#   0 on success, exits with 1 if the mode is not accepted
check_mode_parameter() {
  local mode="$1"
  local command="$2"

  local valid_modes=""
  case "$command" in
    Initialize)
      valid_modes="Simple Starfish"
      ;;
    SyntenyClustering)
      valid_modes="Raw SSP FilterBlast FilterMetric"
      ;;
    CaptainIdentification)
      valid_modes="FullAll AllID Cluster"
      ;;
    OrthogroupsAnnotation)
      valid_modes="All MoveAssociated Core Overrepresented"
      ;;
    OrthogroupsOverrepresentation)
      valid_modes="Outliers Enrichment"
      ;;
    *)
      echo "Error: unknown command '${command}' in check_mode_parameter." >&2
      exit 1
      ;;
  esac

  if [[ ! " ${valid_modes} " =~ " ${mode} " ]]; then
    echo "Error: mode '${mode}' is not accepted for command '${command}'." >&2
    print_help
    exit 1
  fi
}

# Validates that all required software tools are available in PATH for the
# given command.
# Arguments:
#   $1 - command name (Initialize | SyntenyClustering | CaptainIdentification |
#        ClusterCharacterization | OrthogroupsOverrepresentation |
#        OrthogroupsAnnotation | All)
# Returns:
#   0 on success, exits with 1 if any tool is missing
check_required_software() {
  local command="$1"

  if [[ -z "$(which seqkit)" ]]; then
    echo "Error: missing required tool 'seqkit'." >&2
    exit 1
  fi

  local tools=()
  case "$command" in
    Initialize)
      tools=("python" "metaeuk" "agat_sp_extract_sequences.pl")
      ;;
    SyntenyClustering)
      tools=(
        "python" "Rscript" "diamond"
        "blastn" "makeblastdb" "blastdb_aliastool" "blastdbcmd"
      )
      ;;
    CaptainIdentification)
      tools=(
        "python" "macse" "hmmscan"
        "blastn" "makeblastdb" "clipkit" "iqtree3" "gotree" "mafft"
      )
      ;;
    ClusterCharacterization)
      tools=("Rscript" "blastn" "makeblastdb" "orthofinder")
      ;;
    OrthogroupsOverrepresentation)
      tools=("orthofinder" "gotree" "Rscript")
      ;;
    OrthogroupsAnnotation)
      tools=("mafft" "foldseek" "hhblits")
      ;;
    All)
      tools=(
        "python" "metaeuk" "agat_sp_extract_sequences.pl"
        "Rscript" "diamond"
        "blastn" "makeblastdb" "blastdb_aliastool" "blastdbcmd"
        "orthofinder" "gotree"
        "macse" "hmmscan" "clipkit" "iqtree3" "mafft"
        "foldseek" "hhblits"
      )
      ;;
    *)
      echo "Error: unknown command '${command}' in check_required_software." >&2
      exit 1
      ;;
  esac

  for tool in "${tools[@]}"; do
    if [[ -z "$(which "$tool")" ]]; then
      echo "Error: missing required tool '${tool}'." >&2
      exit 1
    fi
  done
}

# Validates the working directory structure, checking that Data/, Workspace/,
# and all four data subdirectories exist with consistent filenames across them.
# Arguments:
#   $1 - base working directory path
# Returns:
#   0 on success, exits with 1 on any structural inconsistency
check_directory_structure() {
  local base_dir="$1"

  mkdir -p "${base_dir}/temp"

  local workspace_dir data_dir
  workspace_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Workspace" 2>/dev/null)
  data_dir=$(find "$base_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

  if [[ -z "$workspace_dir" ]]; then
    echo "Error: Workspace directory not found in '${base_dir}'." >&2
    exit 1
  fi
  if [[ -z "$data_dir" ]]; then
    echo "Error: Data directory not found in '${base_dir}'." >&2
    exit 1
  fi

  local gff_dir protein_dir nucleotide_dir cds_dir
  gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
  protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
  nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
  cds_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

  [[ -z "$gff_dir" ]]        && { echo "Error: Gff subdirectory not found in '${data_dir}'." >&2; exit 1; }
  [[ -z "$protein_dir" ]]    && { echo "Error: Protein subdirectory not found in '${data_dir}'." >&2; exit 1; }
  [[ -z "$nucleotide_dir" ]] && { echo "Error: Nucleotide subdirectory not found in '${data_dir}'." >&2; exit 1; }
  [[ -z "$cds_dir" ]]        && { echo "Error: CDS subdirectory not found in '${data_dir}'." >&2; exit 1; }

  local gff_list="${base_dir}/temp/gff_files.txt"
  local protein_list="${base_dir}/temp/protein_files.txt"
  local nucleotide_list="${base_dir}/temp/nucleotide_files.txt"
  local cds_list="${base_dir}/temp/cds_files.txt"

  find "$gff_dir"        -maxdepth 1 -type f -name "*.gff" \
    | xargs -n1 basename -s .gff | sort > "$gff_list"
  find "$protein_dir"    -maxdepth 1 -type f -name "*.fa" \
    | xargs -n1 basename -s .fa  | sort > "$protein_list"
  find "$nucleotide_dir" -maxdepth 1 -type f -name "*.fa" \
    | xargs -n1 basename -s .fa  | sort > "$nucleotide_list"
  find "$cds_dir"        -maxdepth 1 -type f -name "*.fa" \
    | xargs -n1 basename -s .fa  | sort > "$cds_list"

  local mismatch=false
  local pairs=(
    "Gff:Protein:${gff_list}:${protein_list}"
    "Gff:Nucleotide:${gff_list}:${nucleotide_list}"
    "Gff:CDS:${gff_list}:${cds_list}"
    "Protein:Nucleotide:${protein_list}:${nucleotide_list}"
    "Protein:CDS:${protein_list}:${cds_list}"
    "Nucleotide:CDS:${nucleotide_list}:${cds_list}"
  )

  for pair in "${pairs[@]}"; do
    IFS=':' read -r label_a label_b file_a file_b <<< "$pair"
    if ! diff -q "$file_a" "$file_b" > /dev/null; then
      echo "Error: file lists do not match between ${label_a} and ${label_b}." >&2
      diff "$file_a" "$file_b" >&2
      mismatch=true
    fi
  done

  rm -r "${base_dir}/temp"

  if $mismatch; then
    exit 1
  fi
}

# Validates that an InterProScan installation directory exists and contains
# a working interproscan.sh script and a valid properties file.
# Arguments:
#   $1 - path to the InterProScan installation directory
# Returns:
#   0 on success, exits with 1 if missing or misconfigured
check_interpro_software() {
  local interpro_path="$1"

  if [[ ! -d "$interpro_path" ]]; then
    echo "Error: InterProScan directory '${interpro_path}' does not exist." >&2
    exit 1
  fi

  if [[ ! -f "${interpro_path}/interproscan.sh" ]]; then
    echo "Error: 'interproscan.sh' not found in '${interpro_path}'." >&2
    exit 1
  fi

  if [[ ! $(${interpro_path}/interproscan.sh -version \
      | grep -c "InterProScan") -eq 2 ]]; then
    echo "Error: InterProScan is not properly installed in '${interpro_path}'." >&2
    exit 1
  fi

  local props="${interpro_path}/interproscan.properties"
  if [[ ! -f "$props" ]]; then
    echo "Error: 'interproscan.properties' not found in '${interpro_path}'." >&2
    exit 1
  fi

  local expected_apps="antifam|gene3d|hamap|ncbifam|panther|pfam-a|pirsf|pirsr|sfld|superfamily"
  if [[ ! $(grep -Ec "$expected_apps" "$props") -eq 10 ]]; then
    echo "Error: 'interproscan.properties' is missing required application entries." >&2
    exit 1
  fi
  if [[ ! $(grep -E "$expected_apps" "$props" \
      | awk -F '=' '{print NF}' | sort -u) -eq 2 ]]; then
    echo "Error: 'interproscan.properties' has malformed application entries." >&2
    exit 1
  fi
}

# Validates that required annotation database files exist for the given command.
# Arguments:
#   $1 - path to the databases directory
#   $2 - command name (OrthogroupsAnnotation | CaptainIdentification | All)
# Returns:
#   0 on success, exits with 1 if any database file or directory is missing
check_databases() {
  local database_path="$1"
  local command="$2"

  if [[ ! -d "$database_path" ]]; then
    echo "Error: database directory '${database_path}' does not exist." >&2
    exit 1
  fi

  case "$command" in
    OrthogroupsAnnotation|All)
      local foldseek_path="${database_path}/Foldseek"
      local hhsuite_path="${database_path}/hhsuite"

      if [[ ! -d "$foldseek_path" ]]; then
        echo "Error: 'Foldseek' directory not found in '${database_path}'." >&2
        exit 1
      fi
      if [[ ! -d "${foldseek_path}/weights" ]]; then
        echo "Error: 'weights' directory not found in '${foldseek_path}'." >&2
        exit 1
      fi
      if [[ ! -f "${foldseek_path}/weights/prostt5-f16.gguf" ]]; then
        echo "Error: Foldseek weight file 'prostt5-f16.gguf' not found." >&2
        exit 1
      fi

      if [[ ! -d "$hhsuite_path" ]]; then
        echo "Error: 'hhsuite' directory not found in '${database_path}'." >&2
        exit 1
      fi
      if [[ ! -f "${hhsuite_path}/pfam.md5sum" ]]; then
        echo "Error: 'pfam.md5sum' not found in '${hhsuite_path}'." >&2
        exit 1
      fi
      if [[ ! $(ls "${hhsuite_path}/pfam_"* 2>/dev/null | wc -l) -eq 6 ]]; then
        echo "Error: missing hhblits PfamA database files in '${hhsuite_path}'." >&2
        exit 1
      fi

      if [[ "$command" == "All" ]]; then
        if [[ ! -f "${database_path}/Captains.fa" ]]; then
          echo "Error: 'Captains.fa' not found in '${database_path}'." >&2
          exit 1
        fi
        check_fasta_protein "${database_path}/Captains.fa"
        if [[ ! -f "${database_path}/Captains_CDS.fa" ]]; then
          echo "Error: 'Captains_CDS.fa' not found in '${database_path}'." >&2
          exit 1
        fi
        check_fasta_dna "${database_path}/Captains_CDS.fa"
      fi
      ;;

    CaptainIdentification)
      if [[ ! -f "${database_path}/Captains_CDS.fa" ]]; then
        echo "Error: 'Captains_CDS.fa' not found in '${database_path}'." >&2
        exit 1
      fi
      check_fasta_dna "${database_path}/Captains_CDS.fa"

      if [[ ! -f "${database_path}/Captains.fa" ]]; then
        echo "Error: 'Captains.fa' not found in '${database_path}'." >&2
        exit 1
      fi
      check_fasta_protein "${database_path}/Captains.fa"
      ;;

    *)
      echo "Error: unknown command '${command}' in check_databases." >&2
      exit 1
      ;;
  esac
}

# Validates that all required files for a given Foldseek database exist.
# Arguments:
#   $1 - path to the Foldseek database directory
#   $2 - database name (pdb | afdb_swissprot | All)
# Returns:
#   0 on success, exits with 1 if any database file is missing
check_foldseek_databases() {
  local foldseek_path="$1"
  local foldseekdb="$2"

  case "$foldseekdb" in
    pdb)
      if [[ ! $(ls "${foldseek_path}/${foldseekdb}"* 2>/dev/null | wc -l) -eq 40 ]]; then
        echo "Error: missing files for Foldseek '${foldseekdb}' database." >&2
        exit 1
      fi
      if [[ ! -f "${foldseek_path}/entries_update.idx" ]]; then
        echo "Error: 'entries_update.idx' not found in '${foldseek_path}'." >&2
        exit 1
      fi
      ;;
    afdb_swissprot)
      if [[ ! -f "${foldseek_path}/Accession_swissprot.txt" ]]; then
        echo "Error: 'Accession_swissprot.txt' not found in '${foldseek_path}'." >&2
        exit 1
      fi
      if [[ ! $(ls "${foldseek_path}/${foldseekdb}"* 2>/dev/null | wc -l) -eq 16 ]]; then
        echo "Error: missing files for Foldseek '${foldseekdb}' database." >&2
        exit 1
      fi
      ;;
    All)
      if [[ ! $(ls "${foldseek_path}/pdb"* 2>/dev/null | wc -l) -eq 40 ]]; then
        echo "Error: missing files for Foldseek 'pdb' database." >&2
        exit 1
      fi
      if [[ ! -f "${foldseek_path}/entries_update.idx" ]]; then
        echo "Error: 'entries_update.idx' not found in '${foldseek_path}'." >&2
        exit 1
      fi
      if [[ ! -f "${foldseek_path}/Accession_swissprot.txt" ]]; then
        echo "Error: 'Accession_swissprot.txt' not found in '${foldseek_path}'." >&2
        exit 1
      fi
      if [[ ! $(ls "${foldseek_path}/afdb_swissprot"* 2>/dev/null | wc -l) -eq 16 ]]; then
        echo "Error: missing files for Foldseek 'afdb_swissprot' database." >&2
        exit 1
      fi
      ;;
    *)
      echo "Error: unknown Foldseek database '${foldseekdb}'." >&2
      exit 1
      ;;
  esac
}

# Validates that the InterProScan installation is working and all databases
# are present for a full pipeline run.
# Arguments:
#   $1 - path to the StarCrew main installation directory
# Returns:
#   0 on success, exits with 1 on any missing component
check_installation() {
  local main_dir="$1"

  log_info "Checking InterProScan installation..."
  check_interpro_software "${main_dir}/interproscan"

  log_info "Checking databases..."
  check_databases "${main_dir}/databases" "All"
  check_foldseek_databases "${main_dir}/databases/Foldseek" "All"
}