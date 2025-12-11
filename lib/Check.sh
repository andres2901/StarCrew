#!/usr/bin/env bash

check_installation() {
    local Main_directory="$1"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking Interpro installation..."
    check_interpro_software "${Main_directory}/interproscan/"

    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Checking databases..."
    check_databases "${Main_directory}/databases/" "All"
    check_foldseek_databases "${Main_directory}/databases/Foldseek/" "All"
}

check_main_directory(){
    local Main_directory="$1"

    if [[ ! -d "${Main_directory}/bin" ]]; then
        echo "Error: folder '${Main_directory}/bin' does not exist."
        exit 1
    else
        if [[ ! -f "${Main_directory}/bin/StarClust" ]]; then
            echo "Error: StarClust script was not found in '${Main_directory}/bin'."
            exit 1
        else
            bash ${Main_directory}/bin/StarClust -help
        fi
    fi

    if [[ ! -d "${Main_directory}/main" ]]; then
        echo "Error: folder '${Main_directory}/main' does not exist."
        exit 1
    else
        if [[ ! -f "${Main_directory}/main/CaptainIdentification.sh" || ! -f "${Main_directory}/main/ClusterCharacterization.sh" || ! -f "${Main_directory}/main/Initialize.sh" || ! -f "${Main_directory}/main/OrthogroupsAnnotation.sh" || ! -f "${Main_directory}/main/SyntenyClustering.sh" || ! -f "${Main_directory}/main/TEPrediction.sh" ]]; then
            echo "Error: There's a missing Main script in '${Main_directory}/main' folder."
            exit 1
        fi
    fi

    if [[ ! -d "${Main_directory}/aux" ]]; then
        echo "Error: folder '${Main_directory}/aux' does not exist."
        exit 1
    else
        check_auxiliary_scripts "${Main_directory}/aux" "All"
    fi

    if [[ ! -d "${Main_directory}/hmm" ]]; then
        echo "Error: directory '$hmmprofile_path' does not exist."
        exit 1
    else
        hmmprofile_path="${Main_directory}/hmm"
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

    if [[ ! -d "${Main_directory}/databases" ]]; then
        echo "Error: folder '${Main_directory}/databases' does not exist."
        exit 1
    else
        if [[ ! -f "${Main_directory}/databases/Captains.fa" ]]; then
            echo "Error: fasta file with captains protein information was not found was not found in '${Main_directory}/databases'."
            exit 1
        else
            check_fasta_protein "${Main_directory}/databases/Captains.fa"
        fi
        if [[ ! -f "${Main_directory}/databases/Captains_exon.fa" ]]; then
            echo "Error: fasta file with captains exon information was not found was not found in '${Main_directory}/databases'."
            exit 1
        else
            check_fasta_protein "${Main_directory}/databases/Captains_exon.fa"
        fi
    fi

    # Check for basic files
    if [[ ! -f "${Main_directory}/agat_config.yaml" ]]; then
        echo "Error: 'agat_config.yaml' file is not found in the '${Main_directory}' folder."
        exit 1
    fi

    if [[ ! -f "${Main_directory}/json_updater.py" ]]; then
        echo "Error: 'json_updater.py' script is not found in the '${Main_directory}' folder."
        exit 1
    fi

    if [[ ! -f "${Main_directory}/StarClust_environment.yml" ]]; then
        echo "Error: 'StarClust_environment.yml' file is not found in the '${Main_directory}' folder."
        exit 1
    fi
}

check_working_directory() {
	local working_directoy="$1"

	if [[ ! -d "$working_directoy" ]]; then
        echo "Error: folder '$working_directoy' does not exist."
        exit 1
    else
    if [[ ! "${working_directoy:0:1}" == "/" ]]; then
        working_directoy=$(realpath $working_directoy)
    fi
    fi

    echo $working_directoy
}

check_threads() {
	local threads="$1"

	if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
        echo "Error: '$threads' is not a positive integer."
        print_help
        exit 1
    fi
}

check_fasta_dna() {
	local fasta_path="$1"

	if [[ ! -f "$fasta_path" ]]; then
        echo "Error: file '$fasta_path' does not exist."
        exit 1
    else
    	# Check if fasta file is valid
    	local Fasta_check=$(seqkit seq --quiet --seq-type dna -v $fasta_path)
    	if [[ "$Fasta_check" == "" ]]; then
    		echo "Error: '${fasta_path}' is an invalid DNA fasta file"
    		exit 1
        else
            return 0
        fi
    fi
}

check_fasta_protein() {
	local fasta_path="$1"

	if [[ ! -f "$fasta_path" ]]; then
        echo "Error: file '$fasta_path' does not exist."
        exit 1
    else
    	# Check if fasta file is valid
    	local Fasta_check=$(seqkit seq --quiet -t protein -v $fasta_path)
    	if [[ "$Fasta_check" == "" ]]; then
    		echo "Error: '${fasta_path}' is an invalid protein fasta file"
    		exit 1
        else
            return 0
    	fi
    fi
}

check_gff_file() {
	local gff_path="$1"

	if [[ ! -f "$gff_path" ]]; then
        echo "Error: file '${gff_path}' does not exist."
        exit 1
    else
    	# Check if fasta file is valid
    	local gff_check=$(grep -v "^#" $gff_path | awk -F '\t' '{print NF}' | sort -u)

    	if [[ $gff_check -eq 9 ]]; then
    		local gff_check2=$(grep -v "^#" $gff_path | awk -F '\t' '{print $7}' | egrep -v -c "\+|-|\.")

    		if [[ $gff_check2 -eq 0 ]]; then
    			local gff_check3=$(grep -v "^#" $gff_path | egrep -v -c "ID=")

    			if [[ $gff_check3 -eq 0 ]]; then
                    gff_check4=$(grep -v "^#" $gff_path | egrep -v -c "Parent=")
                    if [[ $gff_check4 == $(grep -v "^#" $gff_path | egrep -c "gene") ]]; then
                        return 0
                    else
                        echo "Error: There are lines without Parent_ID that are not gene in the gff file '${gff_path}'"
                        exit 1
                    fi
                else
                	echo "Error: There are lines without ID in the gff file '${gff_path}'"
            	    exit 1
                fi
            else
            	echo "Error: Check gff file '${gff_path}' due to inconsistencies in the strand column"
            	exit 1
            fi
        else
        	echo "Error: Check gff file '${gff_path}' due to inconsistencies in the number of columns"
        	exit 1
        fi
    fi
}

check_gff_paths() {
    local gff_path="$1"

    if [[ ! -f "$gff_path" ]]; then
        echo "Error: file '$gff_path' does not exist."
        exit 1
    else
        # Check if fasta file is valid
        local gff_check=$(awk -F '\t' '{print NF}' $gff_path | sort -u)

        if [[ $gff_check -eq 2 ]]; then
            cat ${gff_path} | while read line
            do
                GFF_FILE=$(grep -w "${line}" ${gff_path} | awk '{print $2}')
                if [[ ! -f ${GFF_FILE} ]]; then
                    echo "ERROR: Do not find GFF file for $(grep "${line}" ${gff_path} | awk '{print $1}') genome"
                    exit 1
                fi
            done
            return 0
        else
            echo "Error: Check gff file due to inconsistencies in the number of columns"
            exit 1
        fi
    fi
}

check_boundaries_file() {
    local gff_path="$1"

    if [[ ! -f "$gff_path" ]]; then
        echo "Error: file '$gff_path' does not exist."
        exit 1
    else
        # Check if fasta file is valid
        local gff_check=$(grep -v "^#" $gff_path | awk -F '\t' '{print NF}' | sort -u)

        if [[ $gff_check -eq 21 ]]; then
            local gff_check2=$(grep -v "^#" $gff_path | awk -F '\t' '{print $7}' | egrep -v -c "\+|-")
            local gff_check3=$(grep -v "^#" $gff_path | awk -F '\t' '{print $4}' | egrep -w -v -c "[0-9]]")
            local gff_check4=$(grep -v "^#" $gff_path | awk -F '\t' '{print $5}' | egrep -w -v -c "[0-9]]")

            if [[ $gff_check2 -eq 0 || $gff_check3 -eq 0 || $gff_check4 -eq 0  ]]; then
                return 0
            else
                echo "Error: Check gff file due to inconsistencies in its columns"
                exit 1
            fi
        else
            echo "Error: Check boundary file due to inconsistencies in the number of columns."
            exit 1
        fi
    fi
}

check_metadata_file () {
	local metadata_path=$1
	local fasta_path=$2

    if [[ ! -f "$metadata_path" ]]; then
        echo "Error: file '$metadata_path' does not exist."
        exit 1
    else
        if [[ -s $metadata_path ]]; then
        	if [[ $(head -n1 $metadata_path | grep -c "ElementID") -eq 1 ]]; then

        	    local Column_number=$(awk -F ';' '{print NF}' $metadata_path | sort -u)
        	    local Columns=$(awk -F ';' '{print NF}' $metadata_path | sort -u | wc -l)

            	if [[ $Columns -ge 2 ]]; then
            		echo "Error: There are discordances in the number of columns per row in the metadata file."
            		exit
        	    else
                    if [[ $Column_number -eq 1 ]]; then
        	            echo "Error: There's only one column in metadata file"
            	        exit 1
                    else
                        if [[ $(head -n1 $metadata_path | tr ';' '\n' | sort | uniq -c | awk '{print $1}' |sort -u) -eq 1 ]]; then

                            local Headers_name=$(grep "^>" $fasta_path | awk -F '>' '{print $2}' | tr '\n' '|' | sed 's/|$//')
        	                local Headers_count=$(grep -c "^>" $fasta_path)
        	                local Metadata_count=$(egrep $Headers_name $metadata_path | wc -l)
        	                if [[ $Metadata_count -eq 0 ]]; then
        	        	        echo "Error: Metadata do not contain the elements in the fasta file"
        	        	        exit 1
            	            fi
                        else
                            echo "Error: There are duplicated header columns"
                            exit 1
                        fi
        	        fi
        	    fi
        	fi
        fi
    fi
}

check_auxiliary_scripts() {
	local auxiliary_path="$1"
	local command="$2"

    if [[ $command == "Initialize" ]]; then
        if [[ ! -f "${auxiliary_path}/rip_calculator.py" ]]; then
            echo "Error: file '${auxiliary_path}/rip_calculator.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/gff_slicer.py" ]]; then
            echo "Error: file '${auxiliary_path}/gff_slicer.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/merge.py" ]]; then
            echo "Error: file '${auxiliary_path}/merge.py' does not exist."
            exit 1
        fi
    elif [[ $command == "SyntenyClustering" ]]; then
        if [[ ! -f "${auxiliary_path}/Blast_CleanUp.py" ]]; then
            echo "Error: File '${auxiliary_path}/Blast_CleanUp.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/Clustering.py" ]]; then
            echo "Error: File '${auxiliary_path}/Clustering.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/merge_metadata.py" ]]; then
            echo "Error: File '${auxiliary_path}/merge_metadata.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/syntenetAnalysis.R" ]]; then
            echo "Error: File '${auxiliary_path}/syntenetAnalysis.R' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/syntenetPreprocess.R" ]]; then
            echo "Error: File '${auxiliary_path}/syntenetPreprocess.R' does not exist."
            exit 1
        fi
    elif [[ $command == "CaptainIdentification" ]]; then
        if [[ ! -f "${auxiliary_path}/hmmscan_process.py" ]]; then
            echo "Error: file '${auxiliary_path}/hmmscan_process.py' does not exist."
            exit 1
        fi
    elif [[ $command == "ClusterCharacterization" ]]; then
    	if [[ ! -f "${auxiliary_path}/ClusterAnalysis.R" ]]; then
            echo "Error: file '${auxiliary_path}/ClusterAnalysis.R' does not exist."
            exit 1
        fi
    elif [[ $command == "All" ]]; then
        if [[ ! -f "${auxiliary_path}/rip_calculator.py" ]]; then
            echo "Error: file '${auxiliary_path}/rip_calculator.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/gff_slicer.py" ]]; then
            echo "Error: file '${auxiliary_path}/gff_slicer.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/merge.py" ]]; then
            echo "Error: file '${auxiliary_path}/merge.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/Blast_CleanUp.py" ]]; then
            echo "Error: File '${auxiliary_path}/Blast_CleanUp.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/Clustering.py" ]]; then
            echo "Error: File '${auxiliary_path}/Clustering.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/merge_metadata.py" ]]; then
            echo "Error: File '${auxiliary_path}/merge_metadata.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/syntenetAnalysis.R" ]]; then
            echo "Error: File '${auxiliary_path}/syntenetAnalysis.R' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/syntenetPreprocess.R" ]]; then
            echo "Error: File '${auxiliary_path}/syntenetPreprocess.R' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/hmmscan_process.py" ]]; then
            echo "Error: file '${auxiliary_path}/hmmscan_process.py' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/ClusterAnalysis.R" ]]; then
            echo "Error: file '${auxiliary_path}/ClusterAnalysis.R' does not exist."
            exit 1
        fi
    else
        echo "Error: command '$command' is incorrect."
        exit 1
    fi
}

check_minimum_gene_content() {
	local minimum_gene_content="$1"

	if [[ "$minimum_gene_content" =~ ^[0-9]+$ ]]; then
        if (( $minimum_gene_content < 5 || $minimum_gene_content > 100 )); then
            echo "Error: '$minimum_gene_content' minimum gene content is not an accepted value."
            print_help
            exit 1
        fi
    else
        echo "Error: '$minimum_gene_content' is not a positive integer."
        print_help
        exit 1
    fi
}

check_mode_parameter() {
	local mode="$1"
	local command="$2"

	if [[ $command == "Initialize" ]]; then
		if [[ "$mode" != "Simple" && "$mode" != "Starfish" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
	elif [[ $command == "SyntenyClustering" ]]; then
		if [[ "$mode" != "Raw" && "$mode" != "SSP" && "$mode" != "FilterBlast" && "$mode" != "FilterMetric" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
    elif [[ $command == "RobustGenePrediction" || $command == "BrakerGenePrediction" || $command == "TEPrediction" ]]; then
    	if [[ "$mode" != "All" && "$mode" != "Cluster" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
    elif [[ $command == "CaptainIdentification" ]]; then
    	if [[ "$mode" != "FullAll" && "$mode" != "AllID"  && "$mode" != "Cluster" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
    elif [[ $command == "OrthogroupsAnnotation" ]]; then
    	if [[ "$mode" != "All" && "$mode" != "MoveAssociated" && "$mode" != "Core" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
    else
    	echo "Error: command '$command' is incorrect."
    	exit 1
    fi
}

check_required_software() {
	local command="$1"

	if [[ -z "$(which seqkit)" ]]; then
        echo "Error: Missing seqkit function."
        exit 1
    fi

	if [[ $command == "Initialize" ]]; then
		if [[ -z "$(which python)" || -z "$(which metaeuk)" || -z "$(which agat_sp_extract_sequences.pl)" ]]; then
            echo "Error: Missing require function(s) for Initialize"
            exit 1
        fi
    elif [[ $command == "SyntenyClustering" ]]; then
    	if [[ -z "$(which python)" || -z "$(which Rscript)" || -z "$(which diamond)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which blastdb_aliastool)" || -z "$(which blastdbcmd)" ]]; then
            echo "Error: Missing require function(s) for SyntenyClustering"
            exit 1
        fi
    elif [[ $command == "TEPrediction" ]]; then
        if [[ -z "$(which earlGreyAnnotationOnly)" ]]; then
            echo "Error: Missing require function(s) for TEPrediction"
            exit 1
        fi
    elif [[ $command == "CaptainIdentification" ]]; then
        if [[ -z "$(which python)" || -z "$(which macse)" || -z "$(which hmmscan)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which clipkit)" || -z "$(which iqtree3)" || -z "$(which gotree)" ]]; then
            echo "Error: Missing require function(s) for CaptainIdentification"
            exit 1
        fi
    elif [[ $command == "ClusterCharacterization" ]]; then
        if [[ -z "$(which Rscript)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which orthofinder)" ]]; then
            echo "Error: Missing require function(s) for ClusterCharacterization"
            exit 1
        fi
    elif [[ $command == "OrthogroupsAnnotation" ]]; then
        if [[ -z "$(which mafft)" || -z "$(which foldseek)" || -z "$(which hhblits)" ]]; then
            echo "Error: Missing require function(s) for OrthogroupsAnnotation"
            exit 1
        fi
    elif [[ $command == "All" ]]; then
        if [[ -z "$(which python)" || -z "$(which metaeuk)" || -z "$(which agat_sp_extract_sequences.pl)" ]]; then
            echo "Error: Missing require function(s) for Initialize"
            exit 1
        fi
        if [[ -z "$(which python)" || -z "$(which Rscript)" || -z "$(which diamond)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which blastdb_aliastool)" || -z "$(which blastdbcmd)" ]]; then
            echo "Error: Missing require function(s) for SyntenyClustering"
            exit 1
        fi
        if [[ -z "$(which earlGreyAnnotationOnly)" ]]; then
            echo "Error: Missing require function(s) for TEPrediction"
            exit 1
        fi
        if [[ -z "$(which python)" || -z "$(which macse)" || -z "$(which hmmscan)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which clipkit)" || -z "$(which iqtree3)" || -z "$(which gotree)" ]]; then
            echo "Error: Missing require function(s) for CaptainIdentification"
            exit 1
        fi
        if [[ -z "$(which Rscript)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which orthofinder)" ]]; then
            echo "Error: Missing require function(s) for ClusterCharacterization"
            exit 1
        fi
        if [[ -z "$(which mafft)" || -z "$(which foldseek)" || -z "$(which hhblits)" ]]; then
            echo "Error: Missing require function(s) for OrthogroupsAnnotation"
            exit 1
        fi
    else
        echo "Error: Unknown command."
        exit 1
    fi
}

check_directory_structure() {
    local base_dir="$1"

    mkdir -p ${base_dir}/temp

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
    
    # Check for consistent filenames across subdirectories
    local gff_files="$base_dir/temp/gff_files.txt"
    local protein_files="$base_dir/temp/protein_files.txt"
    local nucleotide_files="$base_dir/temp/nucleotide_files.txt"
    local exon_files="$base_dir/temp/exon_files.txt"
    
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
    rm -r "${base_dir}/temp"
}

check_interpro_software() {
    local Interpro_path="$1"
    if [[ ! -d "$Interpro_path" ]]; then
        echo "Error: Sofware of Interpro software does not exist."
        exit 1
    else
        if [[ ! -f "${Interpro_path}/interproscan.sh" ]]; then
            echo "Error: interproscan script file is not available in the '${Interpro_path}' folder."
            exit 1
        else
            if [[ ! $(${Interpro_path}/interproscan.sh -version | grep -c "InterProScan") -eq 2 ]]; then
               echo "Error: interproscan script is not properly install."
               exit 1 
           fi
        fi
        if [[ ! -f "${Interpro_path}/interproscan.properties" ]]; then
            echo "Error: property file for interproscan is not available '$Interpro_path' folder."
            exit 1
        else
            if [[ ! $(egrep -c "antifam|gene3d|hamap|ncbifam|panther|pfam-a|pirsf|pirsr|sfld|superfamily" ${Interpro_path}/interproscan.properties) -eq 10 ]]; then
                echo "Error: property file for interproscan is corrupted."
                exit 1
            elif [[ ! $(egrep "antifam|gene3d|hamap|ncbifam|panther|pfam-a|pirsf|pirsr|sfld|superfamily" ${Interpro_path}/interproscan.properties | awk -F '=' '{print NF}' | sort -u) -eq 2 ]]; then
                echo "Error: property file for interproscan is corrupted."
                exit 1
            fi
        fi
    fi
}

check_databases() {
    local database_path="$1"
    local command="$2"

    if [[ $command == "OrthogroupsAnnotation" ]]; then
        if [[ ! -d "$database_path" ]]; then
            echo "Error: Directory '$database_path' does not exist."
            exit 1
        else
            if [[ ! -d "${database_path}/Foldseek/" ]]; then
                echo "Error: Directory 'Foldseek' does not exist in '$database_path'."
                exit 1
            else
                foldseek_path="${database_path}/Foldseek/"

                if [[ ! -d "$foldseek_path/weights/" ]]; then
                    echo "Error: Directory 'weights' does not exist in '${foldseek_path}'."
                    exit 1
                else 
                    if [[ ! -f "$foldseek_path/weights/prostt5-f16.gguf" ]]; then
                        echo "Error: weight file for Foldseek is not available at $foldseek_path/weights/"
                        exit 1
                    fi
                fi
            fi

            if [[ ! -d "${database_path}/hhsuite/" ]]; then
                echo "Error: Directory 'hhsuite' does not exist in '$database_path'."
                exit 1
            else
                hhsuite_path="${database_path}/hhsuite/"
                if [[ ! -f "$hhsuite_path/pfam.md5sum" ]]; then
                    echo "Error: 'pfam' database does not exist in '${hhsuite_path}'."
                    exit 1
                fi

                if [[ ! $(ls $hhsuite_path/pfam_* | wc -l) -eq 6 ]]; then
                    echo "Error: There are missing files of hhblits pfam database."
                    exit 1
                fi
            fi
            
        fi
    elif [[ $command == "TEPrediction" ]]; then
        if [[ ! -d "$database_path/MycoMobilome_db/" ]]; then
            echo "Error: file '$database_path/MycoMobilome_db/' does not exist."
            exit 1
        else
            if [[ ! $(ls $database_path/MycoMobilome_db/*fasta | wc -l) -eq 3 ]]; then
                echo "Error: There are missing files for MycoMobilome database."
                exit 1
            fi 
        fi
    elif [[ $command == "All" ]]; then
        if [[ ! -d "$database_path/MycoMobilome_db/" ]]; then
            echo "Error: file '$database_path/MycoMobilome_db/' does not exist."
            exit 1
        else
            if [[ ! $(ls $database_path/MycoMobilome_db/*.fasta | wc -l) -eq 3 ]]; then
                echo "Error: There are missing files for MycoMobilome database."
                exit 1
            fi
    
        fi
        if [[ ! -d "$database_path" ]]; then
            echo "Error: Directory '$database_path' does not exist."
            exit 1
        else
            if [[ ! -d "${database_path}/Foldseek/" ]]; then
                echo "Error: Directory 'Foldseek' does not exist in '$database_path'."
                exit 1
            else
                foldseek_path="${database_path}/Foldseek/"

                if [[ ! -d "$foldseek_path/weights/" ]]; then
                    echo "Error: Directory 'weights' does not exist in '${foldseek_path}'."
                    exit 1
                else 
                    if [[ ! -f "$foldseek_path/weights/prostt5-f16.gguf" ]]; then
                        echo "Error: weight file for Foldseek is not available at $foldseek_path/weights/"
                        exit 1
                    fi
                fi
            fi

            if [[ ! -d "${database_path}/hhsuite/" ]]; then
                echo "Error: Directory 'hhsuite' does not exist in '$database_path'."
                exit 1
            else
                hhsuite_path="${database_path}/hhsuite/"
                if [[ ! -f "$hhsuite_path/pfam.md5sum" ]]; then
                    echo "Error: 'pfam' database does not exist in '${hhsuite_path}'."
                    exit 1
                fi

                if [[ ! $(ls $hhsuite_path/pfam_* | wc -l) -eq 6 ]]; then
                    echo "Error: There are missing files of hhblits pfam database."
                    exit 1
                fi
            fi
            
        fi
    else
        echo "Error: command '$command' is incorrect."
        exit 1
    fi
}

check_foldseek_databases() {
    local foldseek_path="$1"
    local foldseekdb="$2"

    if [[ "$foldseekdb" == "pdb" ]]; then
        if [[ ! $(ls $foldseek_path/${foldseekdb}* | wc -l) -eq 40 ]]; then
            echo "Error: There are missing files for '${foldseekdb}' database."
            exit 1
        fi
        if [[ ! -f "$foldseek_path/entries_update.idx" ]]; then
            echo "Error: 'entries_update.idx' file does not exist in '${foldseek_path}'."
            exit 1
        fi
    elif [[ "$foldseekdb" == "afdb_swissprot" ]]; then
        if [[ ! -f "$foldseek_path/Accession_swissprot.txt" ]]; then
            echo "Error: 'Accession_swissprot.txt' file does not exist in '${foldseek_path}'."
            exit 1
        fi
        if [[ ! $(ls $foldseek_path/${foldseekdb}* | wc -l) -eq 16 ]]; then
            echo "Error: There are missing files for '${foldseekdb}' database."
            exit 1
        fi
    elif [[ "$foldseekdb" == "All" ]]; then
        if [[ ! $(ls $foldseek_path/pdb* | wc -l) -eq 40 ]]; then
            echo "Error: There are missing files for '$pdb' database."
            exit 1
        fi
        if [[ ! -f "$foldseek_path/entries_update.idx" ]]; then
            echo "Error: 'entries_update.idx' file does not exist in '${foldseek_path}'."
            exit 1
        fi
        if [[ ! -f "$foldseek_path/Accession_swissprot.txt" ]]; then
            echo "Error: 'Accession_swissprot.txt' file does not exist in '${foldseek_path}'."
            exit 1
        fi
        if [[ ! $(ls $foldseek_path/afdb_swissprot* | wc -l) -eq 16 ]]; then
            echo "Error: There are missing files for 'afdb_swissprot' database."
            exit 1
        fi
    else
        echo "Error: Database '$foldseekdb' is incorrect."
        exit 1
    fi
}
