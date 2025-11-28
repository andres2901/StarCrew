#!/usr/bin/env bash

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

    return $working_directoy
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
    	local Fasta_check=$(seqkit seq --quiet -t dna -v $fasta_path)
    	if [[ ! -s "$Fasta_check" ]]; then
    		echo "Check input fasta file"
    		exit 1
    	fi
        # Check if the path is absolute
        if [[ ! "${fasta_path:0:1}" == "/" ]]; then
            fasta_path=$(realpath $fasta_path)
        fi

    fi
    return $fasta_path
}

check_fasta_protein() {
	local fasta_path="$1"

	if [[ ! -f "$fasta_path" ]]; then
        echo "Error: file '$fasta_path' does not exist."
        exit 1
    else
    	# Check if fasta file is valid
    	local Fasta_check=$(seqkit seq --quiet -t protein -v $fasta_path)
    	if [[ ! -s "$Fasta_check" ]]; then
    		echo "Check input fasta file"
    		exit 1
    	fi
        # Check if the path is absolute
        if [[ ! "${fasta_path:0:1}" == "/" ]]; then
            fasta_path=$(realpath $fasta_path)
        fi

    fi
    return $fasta_path
}

check_gff() {
	local gff_path="$1"

	if [[ ! -f "$gff_path" ]]; then
        echo "Error: file '$gff_path' does not exist."
        exit 1
    else
    	# Check if fasta file is valid
    	local gff_check=$(grep -v "^#" $gff_path | awk -F '\t' '{print NF}' | sort -u)

    	if [[ $gff_check -eq 9 ]]; then
    		local gff_check2=$(grep -v "^#" $gff_path | awk -F '\t' '{print $7}' | egrep -v -c "+|-|\.")

    		if [[ $gff_check2 -eq 0 ]]; then
    			local gff_check3$(grep -v "^#" $gff_path | egrep -v -c "ID=")

    			if [[ $gff_check3 -eq 0 ]]; then

                    # Check if the path is absolute
                    if [[ ! "${gff_path:0:1}" == "/" ]]; then
                        gff_path=$(realpath $gff_path)
                    fi

                    return $gff_path
                else
                	echo "Error: There are lines without ID in the gff file"
            	    exit 1
                fi
            else
            	echo "Error: Check gff file due to inconsistencies in the strand column"
            	exit 1
            fi
        else
        	echo "Error: Check gff file due to inconsistencies in the number of columns"
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
        if [[ ! "${metadata_path:0:1}" == "/" ]]; then
            metadata_path=$(realpath $metadata_path)
        fi

        if [[ -s $metadata_path ]]; then
        	if [[ $(awk 'NR==1{print $0}' $metadata_path | grep -c "ElementID") -eq 1 ]]; then

        	    local Column_number=$(awk -F ';' '{print NF}' $metadata_path | sort -u)
        	    local Columns=$(awk -F ';' '{print NF}' $metadata_path | sort -u | wc -l)

            	if [[ $Columns -ge 2 ]]; then
            		echo "Error: There are discordances in the number of columns in the metadata file"
            		exit
        	    else
                    if [[ $Column_number -eq 1 ]]; then
        	            echo "Error: There's only one column in metadata file"
            	        exit 1
                    else
            	        local Headers_name=$(grep "^>" $fasta_path | awk -F '>' '{print $2}' | tr '\n' '|' | sed 's/|$//')
        	            local Headers_count=$(grep -c "^>" $fasta_path)
        	            local Metadata_count=$(egrep $Headers_name $metadata_path | wc -l)
        	            if [[ $Metadata_count -eq 0 ]]; then
        	        	    echo "Error: Metadata do not contain the elements in the fasta file"
        	        	    exit 1
            	        fi
        	        fi
        	    fi
        	fi
        fi
    fi

    return $metadata_path
}

check_auxiliary_scripts() {
	local auxiliary_path="$1"
	local command="$2"

	if [[ ! -d "$auxiliary_path" ]]; then
        echo "Error: directory '$auxiliary_path' does not exist."
        exit 1
    else
        auxiliary_path=$(realpath $auxiliary_path)
    fi

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
        if [[ ! -f "${auxiliary_path}/syntenetAnalysis.R" ]]; then
            echo "Error: File '${auxiliary_path}/syntenetAnalysis.R' does not exist."
            exit 1
        fi
        if [[ ! -f "${auxiliary_path}/syntenetPreprocess.R" ]]; then
            echo "Error: File '${auxiliary_path}/syntenetPreprocess.R' does not exist."
            exit 1
        fi
    elif [[ $command == "RobustGenePrediction" || $command == "BrakerGenePrediction" ]]; then
    	if [[ ! -f "${auxiliary_path}/merge.py" ]]; then
            echo "Error: file '${auxiliary_path}/merge.py' does not exist."
            exit 1
        fi

        if [[ ! -f "${auxiliary_path}/gff_filter.py" ]]; then
            echo "Error: file '${auxiliary_path}/gff_filter.py' does not exist."
            exit 1
        fi
    elif [[ $command == "CaptainIdentification" ]]; then
    	if [[ ! -f "${auxiliary_path}/macse.jar" ]]; then
            echo "Error: file '${auxiliary_path}/macse.jar' does not exist."
            exit 1
        fi
    elif [[ $command == "ClusterCharacterization" ]]; then
    	if [[ ! -f "${auxiliary_path}/ClusterAnalysis.R" ]]; then
            echo "Error: file '${auxiliary_path}/ClusterAnalysis.R' does not exist."
            exit 1
        fi
    fi

    return $auxiliary_path
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
		if [[ "$mode" != "Simple" && "$mode" != "Full" && "$mode" != "Starfish" ]]; then
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
		if [[ -z "$(which python)" || -z "$(which metaeuk)" ]]; then
            echo "Error: Missing require function(s) for Initialize"
            exit 1
        fi
    elif [[ $command == "QuickGenePrediction" ]]; then
    	if [[ -z "$(which agat_sp_filter_incomplete_gene_coding_models.pl)" || -z "$(which seqkit)" || -z "$(which metaeuk)" ]]; then
            echo "Error: Missing require function(s) for QuickGenePrediction"
            exit 1
        fi
    elif [[ $command == "SyntenyClustering" ]]; then
    	if [[ -z "$(which python)" || -z "$(which Rscript)" || -z "$(which diamond)" || -z "$(which blastn)" || -z "$(which makeblastdb)" || -z "$(which blastdb_aliastool)" || -z "$(which blastdbcmd)" ]]; then
            echo "Error: Missing require function(s) for SyntenyClustering"
            exit 1
        fi
    elif [[ $command == "BrakerGenePrediction" ]]; then
    	if [[ -z "$(which braker)" || -z "$(which seqkit)" || -z "$(which metaeuk)" ]]; then
            echo "Error: Missing require function(s) for BrakerGenePrediction"
            exit 1
        fi
    fi
}