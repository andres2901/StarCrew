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
        if [[ ! -f "${Main_directory}/bin/StarCrew" ]]; then
            echo "Error: StarCrew script was not found in '${Main_directory}/bin'."
            exit 1
        fi
    fi

    if [[ ! -d "${Main_directory}/main" ]]; then
        echo "Error: folder '${Main_directory}/main' does not exist."
        exit 1
    else
        scripts_file=("CaptainIdentification" "ClusterCharacterization" "Initialize" "OrthogroupsAnnotation" "OrthogroupsOverrepresentation" "SyntenyClustering")
        for file in "${scripts_file[@]}"; do
            if [[ ! -f "${Main_directory}/main/${file}.sh" ]]; then
                echo "Error: Missing HMM data for '${file}'."
                exit 1
            fi
        done
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
        # Validate binary compressed datafiles for hmmscan
        domains=("CAT_domain" "DUF3435" "Captain")
        extensions=("" ".h3f" ".h3i" ".h3m" ".h3p")

        for dom in "${domains[@]}"; do
            for ext in "${extensions[@]}"; do
                if [[ ! -f "${hmmprofile_path}/${dom}.hmm${ext}" ]]; then
                    echo "Error: Missing HMM data for '${dom}' (${ext})."
                    exit 1
                fi
            done
        done
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
        if [[ ! -f "${Main_directory}/databases/Captains_CDS.fa" ]]; then
            echo "Error: fasta file with captains CDS information was not found was not found in '${Main_directory}/databases'."
            exit 1
        else
            check_fasta_dna "${Main_directory}/databases/Captains_CDS.fa"
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
                    return 0
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
            echo "Error: Check ome2gff file due to inconsistencies in the number of columns"
            exit 1
        fi
    fi
}

check_boundaries_file() {
    local boundaries_path="$1"

    if [[ ! -f "$boundaries_path" ]]; then
        echo "Error: file '$boundaries_path' does not exist."
        exit 1
    else
        # Check if fasta file is valid
        local boundaries_check=$(grep -v "^#" $boundaries_path | awk -F '\t' '{print NF}' | sort -u)

        if [[ $boundaries_check -eq 21 ]]; then
            local boundaries_check2=$(grep -v "^#" $boundaries_path | awk -F '\t' '{print $7}' | egrep -v -c "\+|-")
            local boundaries_check3=$(grep -v "^#" $boundaries_path | awk -F '\t' '{print $4}' | egrep -w -v -c "[0-9]]")
            local boundaries_check4=$(grep -v "^#" $boundaries_path | awk -F '\t' '{print $5}' | egrep -w -v -c "[0-9]]")

            if [[ $boundaries_check2 -eq 0 || $boundaries_check3 -eq 0 || $boundaries_check4 -eq 0  ]]; then
                return 0
            else
                echo "Error: Check boundary file due to inconsistencies in its columns"
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
        scripts_file=("rip_calculator.py" "gff_slicer.py" "merge.py")
        for file in "${scripts_file[@]}"; do
            if [[ ! -f "${auxiliary_path}/${file}" ]]; then
                echo "Error: Missing auxiliary script '${file}'."
                exit 1
            fi
        done
    elif [[ $command == "SyntenyClustering" ]]; then
        scripts_file=("PreCluster.py" "Blast_CleanUp.py" "Clustering.py" "merge_metadata.py" "syntenetAnalysis.R" "syntenetPreprocess.R")
        for file in "${scripts_file[@]}"; do
            if [[ ! -f "${auxiliary_path}/${file}" ]]; then
                echo "Error: Missing auxiliary script '${file}'."
                exit 1
            fi
        done
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
    elif [[ $command == "OrthogroupsOverrepresentation" ]]; then
        if [[ ! -f "${auxiliary_path}/OverrepresentationAnalysis.R" ]]; then
            echo "Error: file '${auxiliary_path}/OverrepresentationAnalysis.R' does not exist."
            exit 1
        fi
    elif [[ $command == "All" ]]; then
        scripts_file=("rip_calculator.py" "gff_slicer.py" "merge.py" "PreCluster.py" "Blast_CleanUp.py" "Clustering.py" "merge_metadata.py" "syntenetAnalysis.R" "syntenetPreprocess.R" "hmmscan_process.py" "ClusterAnalysis.R" "OverrepresentationAnalysis.R")
        for file in "${scripts_file[@]}"; do
            if [[ ! -f "${auxiliary_path}/${file}" ]]; then
                echo "Error: Missing auxiliary script '${file}'."
                exit 1
            fi
        done
    else
        echo "Error: command '$command' is incorrect."
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
    elif [[ $command == "CaptainIdentification" ]]; then
    	if [[ "$mode" != "FullAll" && "$mode" != "AllID"  && "$mode" != "Cluster" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
    elif [[ $command == "OrthogroupsAnnotation" ]]; then
    	if [[ "$mode" != "All" && "$mode" != "MoveAssociated" && "$mode" != "Core" && "$mode" != "Overrepresented" ]]; then
            echo "Error: provided mode '$mode' is not accepted."
            print_help
            exit 1
        fi
    elif [[ $command == "OrthogroupsOverrepresentation" ]]; then
        if [[ "$mode" != "Outliers" && "$mode" != "Enrichment" ]]; then
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
        function_list=("python" "metaeuk" "agat_sp_extract_sequences.pl")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
    elif [[ $command == "SyntenyClustering" ]]; then
        function_list=("python" "Rscript" "diamond" "blastn" "makeblastdb" "blastdb_aliastool" "blastdbcmd")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
    elif [[ $command == "CaptainIdentification" ]]; then
        function_list=("python" "macse" "hmmscan" "blastn" "makeblastdb" "clipkit" "iqtree3" "gotree" "mafft")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
    elif [[ $command == "ClusterCharacterization" ]]; then
        function_list=("Rscript" "blastn" "makeblastdb" "orthofinder")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
    elif [[ $command == "OrthogroupsOverrepresentation" ]]; then
        function_list=("orthofinder" "gotree" "Rscript")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
    elif [[ $command == "OrthogroupsAnnotation" ]]; then
        function_list=("mafft" "foldseek" "hhblits")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
    elif [[ $command == "All" ]]; then
        function_list=("python" "metaeuk" "agat_sp_extract_sequences.pl" "Rscript" "diamond" "blastn" "makeblastdb" "blastdb_aliastool" "blastdbcmd" "orthofinder" "gotree" "macse" "hmmscan" "clipkit" "iqtree3" "mafft" "foldseek" "hhblits")
        for fun in "${function_list[@]}"; do
            if [[ -z "$(which $fun)" ]]; then
                echo "Error: Missing require function '${fun}'."
                exit 1
            fi
        done
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
    local CDS_dir=$(find "$data_dir" -maxdepth 1 -type d -name "CDS" 2>/dev/null)

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

    if [[ -z "$CDS_dir" ]]; then
        echo "Error: CDS subdirectory not found in '$data_dir'." >&2
        exit 1
    fi
    
    # Check for consistent filenames across subdirectories
    local gff_files="$base_dir/temp/gff_files.txt"
    local protein_files="$base_dir/temp/protein_files.txt"
    local nucleotide_files="$base_dir/temp/nucleotide_files.txt"
    local CDS_files="$base_dir/temp/b_files.txt"
    
    # Get sorted list of base filenames from the GFF directory
    find "$gff_dir" -maxdepth 1 -type f -name "*.gff" | xargs -n 1 basename -s .gff | sort > "$gff_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$protein_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$protein_files"

    # Get sorted list of base filenames from the GFF directory
    find "$nucleotide_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$nucleotide_files"
    
    # Get sorted list of base filenames from the protein directory
    find "$CDS_dir" -maxdepth 1 -type f -name "*.fa" | xargs -n 1 basename -s .fa | sort > "$CDS_files"

    # Compare the lists. If diff finds a difference, it returns a non-zero exit code.
    if ! diff -q "$gff_files" "$protein_files" >/dev/null || \
        ! diff -q "$gff_files" "$nucleotide_files" >/dev/null || \
        ! diff -q "$gff_files" "$CDS_files" >/dev/null || \
        ! diff -q "$protein_files" "$nucleotide_files" >/dev/null || \
        ! diff -q "$protein_files" "$CDS_files" >/dev/null || \
        ! diff -q "$nucleotide_files" "$CDS_files" >/dev/null; then
        echo "Error: File lists in subdirectories do not match." >&2
        echo "Details:" >&2
        echo "GFF vs. Protein:" >&2
        diff "$gff_files" "$protein_files" >&2
        echo "GFF vs. nucleotide:" >&2
        diff "$gff_files" "$nucleotide_files" >&2
        echo "GFF vs. CDS:" >&2
        diff "$gff_files" "$CDS_files" >&2
        echo "protein vs. nucleotide:" >&2
        diff "$protein_files" "$nucleotide_files" >&2
        echo "protein vs. CDS:" >&2
        diff "$protein_files" "$CDS_files" >&2
        echo "nucleotide vs. CDS:" >&2
        diff "$nucleotide_files" "$CDS_files" >&2
        rm "$gff_files" "$nucleotide_files" "$protein_files" "$CDS_files"
        exit 1
    fi

    # Cleanup temporary files
    rm -r "${base_dir}/temp"
}

check_interpro_software() {
    local Interpro_path="$1"
    if [[ ! -d "${Interpro_path}" ]]; then
        echo "Error: Folder of Interpro software does not exist."
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
    elif [[ $command = "CaptainIdentification" ]]; then
        if [[ ! -d "$database_path" ]]; then
            echo "Error: directory '$database_path' does not exist."
            exit 1
        else
            if [[ ! -f "${database_path}/Captains_CDS.fa" ]]; then
                echo "Error: file '${database_path}/Captains_CDS.fa' does not exist."
                exit 1
            else
                check_fasta_dna "${database_path}/Captains_CDS.fa"
            fi

            if [[ ! -f "${database_path}/Captains.fa" ]]; then
                echo "Error: file '${database_path}/Captains.fa' does not exist."
                exit 1
            else
                check_fasta_protein "${database_path}/Captains.fa"
            fi
        fi
    elif [[ $command == "All" ]]; then
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
