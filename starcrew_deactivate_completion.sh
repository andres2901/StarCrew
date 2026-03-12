#!/usr/bin/env bash
# Removes StarCrew Tab-completion when the conda environment is deactivated.
complete -r StarCrew 2>/dev/null
unset -f _starcrew_complete 2>/dev/null
