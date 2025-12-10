# StarClust: STARship CLUSTering tool

## Description

StarClust is a tool designed to systematically analyze Starships cargo. It composeed by seven command and it create a work folder where everything is going to be properly organize to avoid any issue with downstream analysis. It's design to analyze from captain identification (as a confirmation from upstream analysis) to orthogroups functional annotation in cargo genes.

## Table of contents

- [System requirements](#System-requirements)
- [Software requirements](#Software-requirements)
- [Database requirements](#Database-requirements)
- [Installation](#Installation)
- [Input](#Input)
- [Main commands](#Main-commands)
    - [Initialize](#Initialize)
    - [CaptainIdentification](#CaptainIdentification)
    - [SyntenyClustering](#SyntenyClustering)
    - [OrthogroupsOverrepresentation](#OrthogroupsOverrepresentation)
    - [TEPrediction](#TEPrediction)
    - [ClusterCharacterization](#ClusterCharacterization)
    - [OrthogroupsAnnotation](#OrthogroupsAnnotation)
- [Pipeline modes](#Pipeline-modes)
- [](#)
- [Citing StarClust and software called by StarClust](#Citing-StarClust-and-software-called-by-StarClust)
- [License](#License)

## System requirements

Waiting for full tool development to check the final requirements.

## Software requirements

This tool was specifically written to be run on Linux and requires the following software and dependencies to be installed and accessible via the system path:

- python 3.9.
- java.
- hmmer.
- python with the following packages: networkx, biopython, gffutils, pandas, numpy, scikit-learn.
- seqkit.
- earlgrey.
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
- interproscan.
- hhsuite.
- foldseek.

## Database requirements

The commands related to functional annotation and Transposable element identification require databases that are not distributed under this repositorie and need to be properly dowwnload:

- interproscan related databases.
- pfamA from hhblits.
- foldseek: ProstT5, PDB and alphaphold.
- Mycomobilome.

## Installation

Before installation, ensure you have **Anaconda3** and **git** installed and accessible in your system path, as they are required for the full tool installation. After cloning, execute the `build_StarClust.sh` script to properly install all necessary software dependencies and databases required by the tool's commands, as shown below:

```
# clone the script
git clone https://github.com/andres2901/StarClust.git

cd StarClust/
bash build_StarClust.sh

```

## Input

## Main commands

The tool is composed of seven main commands and some of them are sequential in the way that they should run while others can be run independently for specific purposes.

### Initialize

```
Script to organize the working directory to run the subsequent commands in the workflow.

Syntax: StarClust Initialize [ -help ] -f <filte_path> [ -m <string> -gc <integer> -r <integer> -mg <integer> -o <string> -g <file_path> -b <file_path> -s <character> -c <file_path> -M <file_path> --overwrite ]"

Required args:
-f, --fasta:  multifasta file wih the elements to study.

Required args with Default:
-m, --mode: Mode of the input to initialize (Defaul = Simple) [Available mode: Simple, Starfish]
-gc, --gc: integer value of gc content to filter out elements with too low gc content (Default = 0) [range: 20 - 45]
-r, --rip: integer value of the minimum coverage of the element to be possibly affected by RIP to be filter out (Default = 0) [range: 30 - 80]
-mg, --minGene: Minimum number of genes in an element to be include in the dataset (Default: 8) [range: 5 - 100]
-o, --outDirectory: Specify working directory  name (Default = WorkingDirectory).

Required args in 'Starfish' mode:
-g, --gff: 2 column tsv: genome code, path to GFF. The path should be to the original gff files and not the ones formatted to run starfish.
-b, --boundaries: *.elements.feat file output of 'starfish summary' command.
-s, --separator: character separating genomeID from featureID that was used for starfish run.
-c, --captains: *_tyr.filt_intersect.fas file output of 'starfish annotate' command.

Required args in 'Simple' mode:
-g, --gff: Path to the GFF file containing gene predictions with element-relative coordinates for all elements in the fasta file.

Optional args:
-M, --Metadata: csv file delimited by semicolon with the metadata information: ElemenID;<data1>;<data2>;....
--overwrite: Flag to overwrite in case there is already a previous run (Default: off)
-help: Display this help message.
```

This command will filter the input data base on a series of filters and organize the resulting database in a project folder that will be the starting point for the rest of the commands. The main tree filter that you can use in this approach are three:
 
 - RIP-like signal: In this case will remove any element that have a RIP-like signal. For this purpose, we design a python script inspire in theRipper web server (https://theripper.hawk.rocks/#/home) with some modifications in the default threshold base on [Margolin et al. 1988](https://pmc.ncbi.nlm.nih.gov/articles/PMC1460257/) and [Lewis et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2661801/). In this case, the removal of elements is performed in the following way:
     - Calculation of 'rip product index', 'substrate denominator', and 'rip substrate index' in windows of 1000bp and sliding windows of 500bp.
     - Identification of RIP-like regions based on: 'rip product index' >

 - gc content: In this case will remove any element that have a gc content percentage below a specific treshold. The idea behind this filter is to remove elements that might have gone trough RIP or similar deleterious processess and do not have a RIP-like signal in it. This is to avoid any possible issue in downstream analysis. The allow range is base on previous literature of the general gc content of Pezizomycotina species.

### CaptainIdentification

### SyntenyClustering

### OrthogroupsOverrepresentation

### TEPrediction

### ClusterCharacterization

### OrthogroupsAnnotation

## Pipeline modes

As previously mentioned, this tool can be used for multiple purposes and therefore some commands are sequential. Here, we're going to mention the two main purposes that this tool can be used to and how we recommend the sequential running of the pipelines.

## Citing StarClust and software called by StarClust

Please cite our work if you use `StarClust` in your research:

StarClust is a tool that calls different bioinformatic software, for that reason any publication of results obtained by StarClust required the citation of the tools that were called.

| Command | mode | Dependency | Citation |
|:---:|:---:|:---| :---|
|`Initialize`| `Starfish` | `seqkit`, `agat`, `metaeuk` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Dainat](https://nbisweden.github.io/AGAT/how_to_cite/), [Karin et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32245390/) |
| `Initialize` | `Simple` | `seqkit`, `agat` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Dainat](https://nbisweden.github.io/AGAT/how_to_cite/) |
|`SyntenyClustering`| `Raw`, `SSP`, `FilterMetric` | `seqkit`, `syntenet`, `DIAMOND`, `scikit-learn` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Almeida-Silva et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36539202/), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Pedregosa et al. 2011](https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html) |
|`SyntenyClustering`| `FilterBlast` | `seqkit`, `syntenet`, `DIAMOND`, `scikit-learn`, `blast+` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Almeida-Silva et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36539202/), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Pedregosa et al. 2011](https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/) |
|`CaptainIdentification`| `AllID` | `hmmscan`, `seqkit`, `blast+`, `macse`  | [hmmer](http://hmmer.org/), [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Ranwez et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30165589/) |
|`CaptainIdentification`| `Cluster`, `FullAll`| `hmmscan`, `seqkit`, `blast+`, `macse`, `iqtree3`, `gotree` | [hmmer](http://hmmer.org/), [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Ranwez et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30165589/), [Wong et al. 2025](https://ecoevorxiv.org/repository/view/8916/), [Lemoine & Gascuel 2021](https://pubmed.ncbi.nlm.nih.gov/34396097/) |
|`ClusterCharacterization`| - | `orthofinder`, `DIAMOND`, `blast+`, `ape`, `ggtree`, `gggenomes` | [Emms et al. 2025](https://www.biorxiv.org/content/10.1101/2025.07.15.664860v1), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Paradis et al. 2004](https://pubmed.ncbi.nlm.nih.gov/14734327/), [Yu et al. 2016](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12628), [Hackl et al. 2024](https://arxiv.org/abs/2411.13556) |
|`OrthogroupsAnnotation`| - | `foldseek`, `mafft`, `hhblits`, `interproscan`,| [van Kempen et al 2024](https://pubmed.ncbi.nlm.nih.gov/37156916/), [Katoh & Standley 2013](https://pubmed.ncbi.nlm.nih.gov/23329690/), [Steinegger et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31521110/), [Jones et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24451626/) |
|`TEPrediction`| - | `earlgrey`, `Micomobilome` | [Baril et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38577785/), [Baril & Croll 2025](https://www.biorxiv.org/content/10.1101/2025.10.28.685023v1)  |

## License

Waiting for License decision