#!/usr/bin/env bash

# 1. Create ProgressBar function
# 1.1 Input is currentState($1) and totalState($2)
function ProgressBar {
# Process data
    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
# Build progressbar string lengths
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

# 1.2 Build progressbar strings and print the ProgressBar line
# 1.2.1 Output example:                           
# 1.2.1.1 Progress : [########################################] 100%
printf "\rProgress $(date "+%Y-%m-%d %H:%M:%S") : [${_fill// /#}${_empty// /-}] ${_progress}%%"

}

overwrite() {
    local working_dir="$1"
    local command="$2"

    local workspace="${working_dir}/Workspace/"

    if [[ -d "${working_dir}/${command}/" ]]; then
        echo "Removing '${command}' folder from the working directory."
        rm -r ${working_dir}/${command}/
    fi

    if [[ -d "${workspace}/${command}/" ]]; then
        echo "Removing '${command}' folder from the Workspace directory."
        rm -r ${workspace}/${command}/
    fi

    if [[ ${command} == "SyntenyClustering" ]]; then
        if [[ -d "${working_dir}/Clusters/" ]]; then
            echo "Removing 'Clusters' folder from the working directory."
            rm -r ${working_dir}/Clusters/
        fi
    fi

    if [[ ${command} == "CaptainIdentification" ]]; then
        if [[ -d "${working_dir}/Data/Captainless_elements/" ]]; then
            echo "Restoring captainless elements in the 'Data' folder."

            local data_dir=$(find "$working_dir" -maxdepth 1 -type d -name "Data" 2>/dev/null)

            local gff_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Gff" 2>/dev/null)
            local protein_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Protein" 2>/dev/null)
            local nucleotide_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Nucleotide" 2>/dev/null)
            local exon_dir=$(find "$data_dir" -maxdepth 1 -type d -name "Exon" 2>/dev/null)

            ls ${working_dir}/Data/Captainless_elements/*.gff | xargs -n 1 basename -s .fa | while read line
            do
                mv ${data_dir}/Captainless_elements/${line}.gff ${gff_dir}/
                mv ${data_dir}/Captainless_elements/${line}_protein.fa ${protein_dir}/${line}.fa
                mv ${data_dir}/Captainless_elements/${line}_nucleotide.fa ${nucleotide_dir}/${line}.fa
                mv ${data_dir}/Captainless_elements/${line}_exon.fa ${exon_dir}/${line}.fa
            done
            rm -r ${working_dir}/Data/Captainless_elements/
        fi
    fi
}