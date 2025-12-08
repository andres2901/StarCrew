# StarClust: STARship CLUSTering tool

## Description

StarClust is a tool designed to systematically analyze Starships cargo. It composeed by seven command and it create a work folder where everything is going to be properly organize to avoid any issue with subsequent analysis.

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

## Citing StarClust and software called by StarClust

Please cite our work if you use `StarClust` in your research:

StarClust is a tool that calls different bioinformatic software, for that reason any publication of results obtained by StarClust required the citation of the tools that were called.

| Command | mode | Dependency | Citation |
|:---:|:---:|:---| :---|
|`Initialize`| `Starfish` | `seqkit`, `agat`, `metaeuk` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Dainat](https://nbisweden.github.io/AGAT/how_to_cite/), [Karin et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32245390/) |
| ^^ | `Simple` | `seqkit`, `agat` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Dainat](https://nbisweden.github.io/AGAT/how_to_cite/) |
|`SyntenyClustering`| `Raw`, `SSP`, `FilterMetric` | `seqkit`, `syntenet`, `DIAMOND`, `scikit-learn` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Almeida-Silva et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36539202/), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Pedregosa et al. 2011] (https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html) |
|`SyntenyClustering`| `FilterBlast` | `seqkit`, `syntenet`, `DIAMOND`, `scikit-learn`, `blast+` | [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Almeida-Silva et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36539202/), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Pedregosa et al. 2011] (https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/) |
|`CaptainIdentification`| `AllID` | `hmmscan`, `seqkit`, `blast+`, `macse`  | [hmmer](http://hmmer.org/), [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Ranwez et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30165589/) |
|`CaptainIdentification`| `Cluster`, `FullAll`| `hmmscan`, `seqkit`, `blast+`, `macse`, `iqtree3`, `gotree` | [hmmer](http://hmmer.org/), [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Ranwez et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30165589/), [Wong et al. 2025](https://ecoevorxiv.org/repository/view/8916/), [Lemoine & Gascuel 2021](https://pubmed.ncbi.nlm.nih.gov/34396097/) |
|`ClusterCharacterization`| - | `orthofinder`, `DIAMOND`, `blast+`, `ape`, `ggtree`, `gggenomes` | [Emms et al. 2025](https://www.biorxiv.org/content/10.1101/2025.07.15.664860v1), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Paradis et al. 2004](https://pubmed.ncbi.nlm.nih.gov/14734327/), [Yu et al. 2016](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12628), [Hackl et al. 2024](https://arxiv.org/abs/2411.13556) |
|`OrthogroupsAnnotation`| - | `foldseek`, `mafft`, `hhblits`, `interproscan`,| [van Kempen et al 2024](https://pubmed.ncbi.nlm.nih.gov/37156916/), [Katoh & Standley 2013](https://pubmed.ncbi.nlm.nih.gov/23329690/), [Steinegger et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31521110/), [Jones et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24451626/) |
|`TEPrediction`| - | `earlgrey`, `Micomobilome` | [Baril et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38577785/), [Baril & Croll 2025](https://www.biorxiv.org/content/10.1101/2025.10.28.685023v1)  |
