#!/usr/bin/env bash
# ==============================================================================
# NAME:        starcrew_completion.sh
# DESCRIPTION: Bash completion for the StarCrew command-line interface.
#              Provides Tab-completion for subcommands and their flags,
#              including value suggestions for enumerated arguments.
#
# INSTALL (conda environment):
#   mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d"
#   mkdir -p "${CONDA_PREFIX}/etc/conda/deactivate.d"
#   cp starcrew_completion.sh "${CONDA_PREFIX}/etc/conda/activate.d/"
#   echo 'complete -r StarCrew 2>/dev/null; unset _starcrew_complete' \
#     > "${CONDA_PREFIX}/etc/conda/deactivate.d/starcrew_completion.sh"
#
# MANUAL LOAD (current session only):
#   source starcrew_completion.sh
# ==============================================================================

_starcrew_complete() {
  local cur prev words cword
  _init_completion || return

  # All available subcommands
  local subcommands=(
    Initialize
    GeneHomogenization
    CaptainIdentification
    SyntenyClustering
    ClusterCharacterization
    OrthogroupsOverrepresentation
    OrthogroupsAnnotation
  )

  # ============================================================
  # Position 1: complete subcommand names
  # ============================================================
  if (( cword == 1 )); then
    COMPREPLY=( $(compgen -W "${subcommands[*]}" -- "$cur") )
    return
  fi

  local subcmd="${words[1]}"

  # ============================================================
  # Flags shared across all subcommands
  # ============================================================
  local common_flags="-help --overwrite -w --workingDirectory -t --threads"

  # ============================================================
  # Per-subcommand flag definitions
  # ============================================================
  local flags=""
  case "$subcmd" in

    Initialize)
      flags="${common_flags} 
             -f --fasta
             -g --gff
             -m --mode
             -gc --gc
             -r --rip
             -mg --minGene
             -o --outDirectory
             -M --Metadata
             -b --boundaries
             -c --captains
             -s --separator"
      ;;

    CaptainIdentification)
      flags="${common_flags}
             -m --mode
             -l --length
             -c --confidenceLevel
             -r --rangeKb
             -ms --minSize"
      ;;

    SyntenyClustering)
      flags="${common_flags}
             -m --mode
             -a --anchors
             -g --gaps
             -e --evalue
             -n --minNodes
             -s --minSize
             -th --threshold
             -fs --fragmentSize
             -ms --mergeSize
             -i --identity
             -c --coverage
             --preCluster
             --captainInfo
             --skip-syntenet"
      ;;

    ClusterCharacterization)
      flags="${common_flags}
             -l --length
             -i --identity
             -c --constellation"
      ;;

    OrthogroupsOverrepresentation)
      flags="${common_flags}
             -m --mode
             -s --scoreMatrix
             -c --coefficient
             -cm --countMode
             -n --name
             -v --value
             -p --pValue
             -pa --pAdj
             --skip-orthofinder"
      ;;

    OrthogroupsAnnotation)
      flags="${common_flags}
             -m --mode
             -f --foldseekdb"
      ;;

  esac

  # ============================================================
  # Value completion for flags with fixed option sets
  # ============================================================
  case "$prev" in

    # --- Shared ---
    -m|--mode)
      case "$subcmd" in
        CaptainIdentification)
          COMPREPLY=( $(compgen -W "Cluster FullAll AllID" -- "$cur") )
          ;;
        SyntenyClustering)
          COMPREPLY=( $(compgen -W "Raw SSP FilterBlast FilterMetric Constellation" -- "$cur") )
          ;;
        OrthogroupsOverrepresentation)
          COMPREPLY=( $(compgen -W "Outliers Enrichment" -- "$cur") )
          ;;
        OrthogroupsAnnotation)
          COMPREPLY=( $(compgen -W "Core ShareAccessory SubAccessory Overrepresented" -- "$cur") )
          ;;
        Initialize)
          COMPREPLY=( $(compgen -W "Simple Starfish" -- "$cur") )
          ;;
      esac
      return
      ;;

    # --- OrthogroupsOverrepresentation ---
    -s|--scoreMatrix)
      COMPREPLY=( $(compgen -W \
        "BLOSUM45 BLOSUM50 BLOSUM62 BLOSUM80 BLOSUM90 PAM250 PAM70 PAM30" \
        -- "$cur") )
      return
      ;;
    -cm|--countMode)
      COMPREPLY=( $(compgen -W "Gene Ship" -- "$cur") )
      return
      ;;
    -pa|--pAdj)
      COMPREPLY=( $(compgen -W \
        "BH BY bonferroni fdr hochberg holm hommel" \
        -- "$cur") )
      return
      ;;

    # --- OrthogroupsAnnotation ---
    -f|--foldseekdb)
      case "$subcmd" in
        OrthogroupsAnnotation)
          COMPREPLY=( $(compgen -W "pdb afdb_swissprot" -- "$cur") )
          ;;
        Initialize)
          COMPREPLY=( $(compgen -f -- "$cur") )
          ;;
      esac
      return
      ;;

  esac

  # ============================================================
  # If cur starts with '-', suggest flags for this subcommand.
  # Otherwise fall back to filesystem completion (files + dirs).
  # ============================================================
  if [[ "$cur" == -* ]]; then
    COMPREPLY=( $(compgen -W "$flags" -- "$cur") )
  else
    COMPREPLY=( $(compgen -f -- "$cur") )
  fi
}

complete -F _starcrew_complete StarCrew
