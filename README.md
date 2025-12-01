# StarClust: STARship CLUSTering tool

## Installation

This tool was specifically written to be run on Linux and requires the following software and dependencies to be installed and accessible via the system path:
- python 3.10.
- java.
- hmmer.
- python with the following packages: networkx, biopython, gffutils, pandas, numpy, scikit-learn.
- seqkit.
- makeblastdb.
- blast.
- clipkit.
- iqtree3.
- gotree.
- orthofinder.
- R v4.4.3 with the following packages: ape, reshape2, viridis, dplyr, gggenomes, scales, dendextend, NbClust, svglite, ggplot2=3.5.2, ggtree, syntenet, optparse.
- metaeuk.
- agat.
- diamond.
- macse.
- braker (the executable command must be named **braker** instead of braker.pl).

To run this toolkit, you must first clone this repository locally. Afterward, grant executable permissions to all scripts located in the `bin/` `aux/` and `main/` folders. 

### sorftware dependencies installation

A single Conda environment, sufficient for most main commands, can be set up by following this installation process:

```
# Create the basic environment with R and python version
conda create -y -n StarClust -c conda-forge r-essentials r-base=4.4.3 python=3.10

# Activate the environment
conda activate StarClust

#Install bioinformatic tools
conda install -y -c bioconda -c conda-forge metaeuk
conda install -y -c bioconda -c conda-forge seqkit agat
conda install -y -c bioconda -c conda-forge orthofinder=3.1.0
conda install -y -c bioconda -c conda-forge gotree hmmer clipkit

# Install packages required for auxiliary scripts
conda install -y conda-forge::openjdk
conda install -y -c bioconda -c conda-forge networkx biopython gffutils pandas numpy scikit-learn

# Install R packages required that can be smoothly install with conda
conda install -y -c bioconda -c conda-forge bioconductor-syntenet bioconductor-ggtree
conda install -y -c conda-forge r-ggplot2=3.5.2

# Start R environment
R

# Install the required R packages
install.packages(c("ape","reshape2","viridis","dplyr","gggenomes","scales","dendextend","NbClust","svglite","optparse"))

## Check installation of all libraries

quit()

# Modify orthofinder config.json file

vi ${Path_to_StarClust_environment}/bin/src/orthofinder/run/config.json

### Standard orthofinder version:

#    "diamond":{ 
#    "program_type": "search",
#    "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
#    "search_cmd": "diamond blastp --ignore-warnings -d DATABASE -q INPUT -o OUTPUT --more-sensitive -p 1 --quiet -e 0.001 --compress 1"
#    },

### Modified version:

#    "diamond":{ 
#    "program_type": "search",
#    "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
#    "search_cmd": "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT --query-cover 90 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND --fast -p 1 --quiet -e 1e-5 --compress 1"
#    },
```

Alternatively, the environment can be set up using the YAML file provided in this repository. After creating the environment, you will need to manually install the required R packages and make the Orthofinder modifications as follows:

```
# Create environment
conda env create -f StarClust_environment.yml

# Activate the environment
conda activate StarClust

# Start R environment
R

# Install the required R packages
install.packages(c("ape","reshape2","viridis","dplyr","gggenomes","scales","dendextend","NbClust","svglite","optparse"))

## Check installation of all libraries

quit()

# Modify orthofinder config.json file

vi ${Path_to_StarClust_environment}/bin/src/orthofinder/run/config.json

### Standard orthofinder version:

#    "diamond":{ 
#    "program_type": "search",
#    "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
#    "search_cmd": "diamond blastp --ignore-warnings -d DATABASE -q INPUT -o OUTPUT --more-sensitive -p 1 --quiet -e 0.001 --compress 1"
#    },

### Modified:

#    "diamond":{ 
#    "program_type": "search",
#    "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
#    "search_cmd": "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT --query-cover 90 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND --fast -p 1 --quiet -e 1e-5 --compress 1"
#    },
```

The aligner MACSE (Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons) must be manually downloaded and copied into the `aux/` folder, using the specific file name `macse.jar`.

For the command `StarClust BrakerGenePrediction`, a separate Conda environment must be set up due to Perl dependency conflicts between BRAKER, AGAT, and recent MetaEuk versions. To set up this environment, follow the instructions provided in the BRAKER repository (https://github.com/Gaius-Augustus/BRAKER). Additionally, ensure that MetaEuk is reachable in the system path, which may require installing an older, compatible version.

## Citing StarClust and software called by StarClust

StarClust is a toolkit that calls different bioinformatic tools, for that reason any publication of results obtained by StarClust required the citation of the tools that were called.

| Command | Dependency | Citation |
|:---:|:---:|:---|
|`Initialize`| `metaeuk` | [Karin et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32245390/) |
|`SyntenyClustering`| `` | [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Marcais et al 2018](https://pubmed.ncbi.nlm.nih.gov/29373581/) |
|`CaptainIdentification`| `` | |
|`ClusterCharacterization`| `` |  |
|`OrthogroupsAnnotation`| `` |  |
