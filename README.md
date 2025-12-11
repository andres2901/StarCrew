# StarClust: STARship CLUSTering tool

## Overview

**StarClust** is a tool specifically designed to systematically analyze the cargo genes of Starship elements. The tool is composed of seven distinct commands and creates a well-organized project folder to facilitate downstream analysis. StarClust provides an initial analysis workflow, ranging from captain identification (to confirm upstream analysis) to the functional annotation of orthogroups within the cargo genes. The primary goal of this tool is to offer researchers a simple, integrated workflow for an initial exploratory analysis and comparison of Starship cargo genes. This process is expected to help identify biological patterns or generate hypotheses that can be further tested in either dry-lab or wet-lab environments.

In addition to the main commands, StarClust is distributed with a diverse set of auxiliary scripts. Although these scripts are primarily used for specific tasks within the main workflow, users can utilize them independently for their own purposes in other bioinformatics settings. These auxiliary scripts cover a diverse range of tasks often encountered in a genomic analysis workflow, including: 1)A modified RIP-like signal calculator, 2)Filtering genes or isoforms from a GFF file based on intron density, 3) Extracting gene information in batch from multiple genomic regions within a GFF3 file, 4) And other specific tasks.

## Table of contents

- [System requirements](#System-requirements)
- [Software requirements](#Software-requirements)
- [Database requirements](#Database-requirements)
- [Installation](#Installation)
- [Input](#Input)
    - [Simple input files](#Simple-input-files)
    - [Starfish input files](#Starfish-input-files)
- [Main commands](#Main-commands)
    - [Initialize](#Initialize)
    - [CaptainIdentification](#CaptainIdentification)
    - [SyntenyClustering](#SyntenyClustering)
    - [OrthogroupsOverrepresentation](#OrthogroupsOverrepresentation)
    - [TEPrediction](#TEPrediction)
    - [ClusterCharacterization](#ClusterCharacterization)
    - [OrthogroupsAnnotation](#OrthogroupsAnnotation)
- [Project folder organization](#Project-folder-organization)
- [Pipeline modes](#Pipeline-modes)
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
# clone the repository
git clone https://github.com/andres2901/StarClust.git

cd StarClust/
bash build_StarClust.sh
```

## Input

The tool currently supports two distinct types of input data. Generally, both modes require a minimum of two kind of information: the element nucleotide sequences and a gene prediction in GFF3 format.

An optional metadata file can also be supplied, containing as much information as the user requires. However, the user must adhere to the following specifications:

- **Format:** It must be a `.csv` file using a semi-colon (`;`) as the separator. There should be no tab characters within the file.
- **Header:** The first row must represent the header, containing the name of each column.
- **Identifier Column:** The first column must represent the name/ID of the element, and its header must be `ElementID`.

### Simple input files

This input format is intended for use when analyzing **manually identified** or **curated** Starship elements, or when you wish to study a specific subset, or the totality, of the elements retrieved from [Starbase](https://starbase.serve.scilifelab.se/).

In this scenario, you only require two basic input files:

1.  **Elements File (Multifasta):** A multifasta file containing the nucleotide sequences of the elements you intend to study.

2.  **Gene Prediction File (GFF3):** A GFF3 file containing the gene predictions for all elements listed in the `Elements` fasta file.
    * The `seqID` in the first column **must be the header/ID of the element** itself, **not** the ID of the contig from which it was extracted.
    * Similarly, the coordinates must be **relative to the element's sequence**, not the contig from which it was extracted.

### Starfish input files

This input mode is used when you have processed your dataset using the **[Starfish toolkit](https://github.com/egluckthaler/starfish)** and wish to analyze the obtained "ships" (elements).

**It is strongly recommended** that you first perform an initial inspection and removal of false positives, as described in the Starfish **[step-by-step tutorial](https://github.com/egluckthaler/starfish/wiki/Step-by-step-tutorial)**.

To use this tool, you must have the results from both the `geneFinder` and `elementFinder` modules of Starfish. **StarClust** utilizes three resulting Starfish files and one user-constructed file:

1.  **Elements File (Multifasta):**
    * This is the `*.elements.fna` file resulting from the `starfish summarize` command.
    * This file contains the nucleotide sequences of the elements to be studied.
    * **Recommendation:** Ideally, the user should remove elements identified as false positives during their preliminary inspection.

2.  **Boundaries File (Metadata):**
    * This is the `*.elements.feat` file resulting from the `starfish summarize` command.
    * This file contains the element metadata from Starfish, which is used for delimiting the element gene content information.
    * **Recommendation:** The user should maintain the header and remove the corresponding false positive elements from this table.

3.  **Captain File (Multifasta):**
    * This is the `*_tyr.filt_intersect.fas` file resulting from the `starfish annotate` command.
    * This file contains the putative YR recombinase genes identified by Starfish.
    * **Note:** While filtering is not strictly required, users working with large datasets can reduce the file to include only Captains of **true positive** elements to decrease run time.

4.  **Gene Prediction Path File (TSV):**
    * This is a user-constructed TSV file with a two-column structure (similar to the optional Starfish input file): `genome code` and `path to GFF`.
    * **Crucial Points:**
        * The `genome code` must be **identical** to the code used in the Starfish analysis.
        * The GFF path must point to the **original GFF file**, and **not** the files returned by the `starfish format` or `starfish format-ncbi` commands.
        * This tool requires the gene, mRNA, intron, exon, and CDS information, and needs the original contig `seqID` to be maintained. The files returned by the `starfish format` or `starfish format-ncbi` commands only retain mRNA information and alter the IDs to suit the Starfish workflow.
    * **Acquisition:** The user can obtain the necessary GFF file using the same command line as detailed in the Starfish step-by-step tutorial:

```
realpath gff3/* | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > ome2gff.txt
```

## Main commands

The tool is composed of seven main commands. Most of these commands are sequential (meaning they must be run in a specific order) while others can be executed independently for specific purposes. The only command that is required for the rest of the workflow and must always be run first is the `StarClust Initialize` command, which generates the dedicated project folder.

### Initialize

```
Script to organize the working directory to run the subsequent commands in the workflow.

Syntax: StarClust Initialize [ -help ] -f <filte_path> -g <file_path> [ -m <string> -gc <integer> -r <integer> -mg <integer> -o <string> -b <file_path> -s <character> -c <file_path> -M <file_path> --overwrite ]"

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
-s, --separator: character separating genomeID from featureID that was used for Starfish run.
-c, --captains: *_tyr.filt_intersect.fas file output of 'starfish annotate' command.

Required args in 'Simple' mode:
-g, --gff: Path to the GFF file containing gene predictions with element-relative coordinates for all elements in the fasta file.

Optional args:
-M, --Metadata: csv file delimited by semicolon with the metadata information: ElementID;<data1>;<data2>;....
--overwrite: Flag to overwrite in case there is already a previous run (Default: off)
-help: Display this help message.
```

This command filters the input data based on a series of criteria and organizes the resulting database within a project folder. This folder serves as the starting point for all subsequent commands. The user should be aware that the element naming may change, with a specific, new ID designated for the Project. However, an association file will always be provided to map the original element ID to the new ID.

The command provides three main filtering approaches you can utilize:

1. RIP-like Signal Filtering: If enabled, this filter removes any element that shows a Repeat-Induced Point Mutation (RIP)-like signal over an specific threshold.
    * For this purpose, we utilize a custom Python script inspired by [theRipper web server](https://theripper.hawk.rocks/#/home), with modifications to the default thresholds based on [Margolin et al. 1988](https://pmc.ncbi.nlm.nih.gov/articles/PMC1460257/) and [Lewis et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2661801/).
    * The removal of elements is performed through the following steps:
        * **Index Calculation:** Calculation of the 'rip substrate index,' 'rip product index,' and 'rip composite index' in 1000 bp windows with a 250 bp sliding step.
        * **RIP-like Window Identification:** A window is identified as RIP-like if the 'rip substrate index' $< 1.05$, the 'rip product index' $> 0.8$, and the 'rip composite index' $> 0$.
        * **RIP-like Region Calculation:** Three sequential RIP-like windows must be identified, and only the overlapping section among these windows is designated as a RIP-like region. This strictness is implemented to avoid over-identifying RIP-like regions.
        * **Coverage Calculation:** The ratio of nucleotides in RIP-like regions to the full element length is calculated. This specific coverage value is used to perform the element filtering.

2. GC Content Filtering: This filter removes any element with a GC content percentage below a specific threshold.
    * The rationale behind this filter is to remove elements that may have undergone RIP or similar deleterious processes but do not exhibit a strong RIP-like signal.
    * The allowed range is based on previous literature concerning the general genome GC content of Pezizomycotina species.

3. Gene Content Filtering: This filter removes elements that are empty or contain a very low number of cargo genes.
    * Since the primary goal of this tool is to analyze cargo content, maintaining elements without sufficient cargo genes is inefficient. This filter ensures that only elements useful for downstream analysis are retained.

In the case of 'Starfish' mode, the input data goes trough a specific preprocess to be able to enter the pipeline:

1. The gene information from the original gff files is extracted and updated base on the positions from the metadata file from starfish run.
2. As Starfish perform a _de novo_ gene prediction of captains, we have to perform again this approach to be able to obtain the full genetic information of this genes (exon, intron, CDS...). For our purpose, we performed a modify version of that approach:
    * Based on the YR recombinase database obtain in Starfish, we run meateuk against the nucleotide file looking for a minimum sequence identity of 0.95 and minimum coverage of 0.95.
    * In the case there was a previous gene model in the same position, this tool will select the longest gene model. This present a difference with Starfish logic that always select metaeuk gene model in those cases.

When using the 'Starfish' input mode, the provided data undergoes specific preprocessing before entering the main pipeline:

1.  **GFF Information Update:** The gene information from the original GFF files is extracted and updated based on the positions provided in the metadata file resulting from the Starfish run.
2.  **Captain Gene Prediction:** As Starfish performs a *de novo* prediction of Captain genes, we must re-execute a similar approach to obtain the complete genetic information (exon, intron, CDS, etc.) for these genes. We use a modified version of the Starfish approach:
    * Using the YR recombinase database obtained from Starfish, MetaEuk is run against the nucleotide file, requiring a minimum sequence identity of 0.95 and a minimum coverage of 0.95.
    * If a previous gene model exists at the same position, this tool selects the longest gene model. This differs from the original Starfish logic, which consistently selects the MetaEuk gene model in such conflict cases.

### CaptainIdentification

```
Script to identify captain genes within each element and construct a phylogenetic tree based on these captains.
It executes five main steps:
1. Run hmmscan using hmm profiles of specific domains in captains against the proteome of each element.
2. Processes the data to identify Captains and regions suitable for phylogenetic analysis. Three minimum confidence levels can be used for Captain identification:
  2.1 Only a match with the Captain HMM profile from Starfish. WARNING: This may lead to false positive identifications and result in an unreliable phylogenetic analysis).
  2.2 A match with the Captain HMM profile plus a match with the DUF3435 HMM profile.
  2.3 A match with the Captain HMM profile and DUF3435 HMM profile plus a match with the Integrase catalytic core HMM profile.
3. For elements lacking a confident Captain gene, the script searches for a putative Captain pseudogene at the beginning and end of the element.
4. Aligns exonic sequences using MACSE (with amino acid output) and preprocesses the alignment with Clipkit.
5. Runs maximum-likelihood phylogenetic tree inference, when there at least two unique captain sequences:
  5.1 Run IQ-TREE with 1000 UFBotstrap and 1000 sh-aLRT if there at least 4 unique sequence, in other case run it without support.
  5.2 Collapsed near-zero and low-support (when available) branches.
  5.3 mid-root the tree.

There are three available mode:
-Cluster: Analyzes and performs all five steps per cluster, and remove elements without a suitable Captain gene or pseudogene from the main dataset.
-FullAll: Analyzes and performs all five steps on the whole dataset, and remove elements without a suitable Captain gene or pseudogene from the main dataset..
-AllID: Analyzes and performs only the first three steps (Identification) on the whole dataset, and remove elements without a suitable Captain gene or pseudogene from the main dataset.

Syntax: StarClust CaptainIdentification [ -help ] -w <directory_path> [ -l <integer> -c <integer> -m <string> -t <integer> -ms <integer> --overwrite ]

Required args:
-w, --workingDirectory: Specify the working directory where all data are stored.

Required args with Default:
-m, --mode: specified the mode (Default = AllID) [Available mode: Cluster, FullAll, AllID].
-l, --length: Minimum length of the protein to be identify as captain (Default: 250) [range: 200 - 800].
-c, --confidenceLevel: Minimum confidence level to call a captain. Note: the script is always going to try to return the captain with the highest level of confidence (Default: 2) [range: 1 - 3].
-r, --rangeKb:The distance (as a number of kilobases) from the beginning or end of the element within which a gene must fall to be considered a captain (Default: 10) [range: 3 - 20]

Required args with Default in 'Cluster' mode:
-ms, --minSize: Minimum size of a Cluster to be include in the analysis when running the 'Cluster' mode (Default = 4) [range: 4 - 10]

Optional args:
-t, --threads: Number of threads to use for phylogenetic tree inference (Default: 1).
--overwrite: Flag to overwrite in case there is already a previous run of CaptainIdentification (Default: off)
-help: Display this help message.
```

The identification of the Captain gene is based on its basic description: a Tyrosine Recombinase (YR recombinase) containing a DUF3435 domain ([Gluck-Thaler et al. 2022](https://academic.oup.com/mbe/article/39/5/msac109/6588634), [Gluck-Thaler & Vogan 2024](https://academic.oup.com/nar/article/52/10/5496/7660083?login=true)). For this purpose, we construct or retrieve the following three sets of HMM profiles:

1.  **Starfish Captain HMM:** The HMM profile for Captain proteins distributed with the [Starfish toolkit](https://github.com/egluckthaler/starfish).
2.  **DUF3435 HMM:** The HMM profile for the DUF3435 domain (Pfam accession: PF11917.13).
3.  **YR Recombinase Active Site HMMs:** A database of HMM profiles associated with the active domain of YR recombinases, including:
    * DNA breaking-rejoining enzymes (Superfamily accession: SSF56349)
    * Phage integrase family (Pfam accession: PF00589.27)
    * Integrase catalytic core (CATH-Gene3D accession: G3DSA:1.10.443.10)
    * TYR_RECOMBINASE (PROSITE accession: PS51898 version 3)
    * TYR_RECOMBINASE_FLP (PROSITE accession: PS51899 version 3)

To be considered a "true" Captain, the predicted gene must pass a series of strict filters. This approach differs from the original Starfish methodology, which relies on a single filter (matching the Captain HMM profile). Our increased strictness aims to remove "false positive" or highly "degraded" elements lacking a strong Captain signal.

The implemented filters are:

1.  **HMM Profile Count:** Must match a minimum number (defined by the user) of the HMM profiles listed above.
2.  **Minimum Length:** The protein must meet a minimum length of amino acids, which is modifiable by the user.
3.  **Exon Count:** The exon count must be between 2 and 11. This constraint is based on the observation that Captains in the current dataset contain mostly intron-containing gene models. Also, preliminary analysis indicated that gene models with an excessive number of introns often correspond to pseudogenes, where gene predictors introduce many introns to avoid intra-frame stop codons.
4.  **Genomic Position:** The gene must be located at the beginning of the element in the positive strand (or at the end if the element sequence is in its reverse complement). The acceptable range for this position is defined by the user.

In the scenario where no gene passes the filters to be called a 'Captain', a pseudogene will be identified if any sequence homology exists based on the exon sequence of known Captain genes.

**Note 1:** Be aware that this command will always remove elements that do not have an identifiable Captain gene or pseudogene, as these elements are not useful for the downstream analysis.

**Note 2:** These filters are still a work in progress as the community continues to gather and document information about the gene model and protein structure of Captains. Currently, available literature lacks comprehensive documentation on this specific gene model.

### SyntenyClustering

```
Script to run the syntenet pipeline using DIAMOND for sequence similarity search, summarizes the results, and organizes the output data.
It executes six main steps:
1. Executes the initial data preprocessing step required by the syntenet pipeline.
2. Runs DIAMOND using the preprocessed data.
3. Runs the interspecies synteny command of syntenet to identify regions with gene collinearity.
4. Summarizes the results of syntenet on four possible modes:
 a. Raw: Return pairs that have a minimum of 8% of shared collinear genes. WARNING: This mode may yield a high rate of false positives.
 b. SSP: Return only pairs with strong synteny.
 c. FilterBlast: Returns pairs that have been filtered using a BLAST-based approach at the nucleotide level.
 d. FilterMetric: Filter and update collinearity based on a metric system.
5. Defines initial clusters of elements and performs a spectral clustering process to identify potential subclusters.
6. Organizes the final data output for each identified cluster.

Syntax: SAT $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -m <string> -a <integer> -g <integer> -n <integer> -s <integer> -t <integer> -th <float> -p <string> --captainInfo --overwrite ]

Required args:
-w, --workingDirectory: Specify the working directory where all data are stored.

Required args with Default:
-m, --mode: Specified the mode to run the summarizing process of synteny results (Default = FilterMetric) [Available mode: Raw, SSP, FilterBlast, FilterMetric].
-a, --anchors: Number of minimum anchor points for syntenet to call a collinear region (Default = 8) [range: 3 - 25].
-g, --gaps: Number of maximum allowed gaps between anchor points for syntenet to call a collinear region (Default = 8) [range: 5 - 25].
-n, --minNodes: Minimum number of nodes in a cluster for spectral clustering to be attempted (Default: 4).
-s, --minSize: The minimum desired size for any final sub-cluster (Default: 2).
-th, --threshold: The minimum modularity score for a split to be accepted (Default: 0.05) [range: -0.5 - 1.0].
 
 Required args with Default in 'FilterBlast' mode:
-fs, --fragmentSize: The minimum fragment size of a blast alignment to be used for blastn filter (Default: 2000) [range: 1000, 5000].
-ms, --mergeSize: The minimum merge fragment size to be used for blastn filter (Default: 5000) [range: 2000, 10000].
-i, --identity: The minimum percentage of identity of a blast alignment to be used for blastn filter (Default: 70.0) [range: 60.0, 90.0].
-c, --coverage: The minimum coverage of the filter merge fragments for a pair to pass the filter (Default: 20.0) [range: 10.0, 50.0].

Optional args:
-t, --threads: Number of threads for searching software (DIAMOND and blast) (Default: 8).
--captainInfo: Flag to check and used the captain exon/pseudoexon information for the cluster (Default: off).
--overwrite: Flag to overwrite in case there is already a previous run of $(basename -s .sh "$0" ) (Default: off).
-help: Display this help message.
```

The clustering of elements is based on collinearity regions. In this case, 

### ClusterCharacterization

```
Script to run the characterization of each selected cluster.
This script perform eight steps per cluster to analyzed:
1. Identify orthogroups through OrthoFinder.
2. Perform all-vs-all Blastn for synteny visualization.
3. Perform a hierarchical clustering of the elements based on Orthogroup gene count including singletons.
4. Determine the full conection of the cluster and create a cargo orthogroups heatmap and synteny image for the cluster.
5. Identify possible individual nesting events inside the cluster.
6. Identify core genes in the cluster in two ways:
 6.1. General core: orthogroups that are present in at least 80% of the elements in the cluster.
 6.2. Specific core: Orthogroups that are present in at least 80% of the elements for subclusters generated at a 0.8 height of the hierarchical tree of cargo content.
   6.2.1. Divide the Cluster in subclusters of a height above 0.8 in the hierarchical clustering.
   6.2.2. If subslusters are presen, identify core genes in each one that have at least 5 elements using the same logic of general core.
7. If subclusters are present it try to identify putative cargo movement events and try to avoid 'General core' genes.
8. Determine if there are discordances at 'Clade' lavel between Cargo hierarchical clustering and Captain phylogenetic tree.

Syntax: StarClust $(basename -s .sh "$0" ) [ -help ] -w <directory_path> [ -t <integer> ]

Required args:
-w, --workingDirectory: Specify the working directory where all data are stored.

Required args with Default:
-t, --threads: Number of threads for orthofinder and blast (Default: 8)

Optional args:
--overwrite: Flag to overwrite in case there is already a previous run of $(basename -s .sh "$0" ) (Default: off)
-help: Display this help message.
```

### OrthogroupsOverrepresentation

```
Waiting to be develop
```

### OrthogroupsAnnotation

```
Script to run a functional annotation for orthogroups.
This script perform three steps:
1. Organize the Orthogroups that are selected base on the mode.
 1.1 Core: Orthogoups that were identify as core by the ClusterCharacteriation command.
 1.2 MoveAssociated: Orthogoups that were identify as part of a putative movement event between subclusters by the ClusterCharacterization command.
 2.3 All: All orthogroups identify by the ClusterCharacterization command.
 2.4 Overrepresented: Orthogroups that were identified as Overrepresented by the OrthogroupsOverrepresentation command.
2. Perform the characterization of the Orthogroup proteins with four approaches:
 2.1 InterProScan: Using all default applications except COILS and MOBIDB.
 2.3 Foldseek: Search for homologs proteins against a database based on the secondary structure.
 2.4 hhblits: Search domains against PfamA database.
3. Summarize the results of the previous step:
 3.1 Internal summary: For each Orthogroups summarize the results per protein in a csv
 3.2 General summary: Return a summary for the Orthogroup under the assumption that all proteins in each Orthogroups have the same function. It return only those 'chracteristics' that are shared for at least 50% of the proteins in the Orthogroup.

Syntax: StarClust OrthogroupsAnnotation [ -help ] -w <directory_path> [ -m <string> -f <string> -t <integer> --overwrite ]

Required args:
-w, --workingDirectory: Specify the working directory where all data are stored.

Required args with Default:
-m, --mode: Define the orthogroups to be analyzed (Default = Core) [Available mode: MoveAssociated, Core, All, Overrepresented].
-f, --foldseekdb: Name of the Foldseek database to use (Default = afdb_swissprot) [Available: pdb, afdb_swissprot].
-t, --threads: Number of threads for all analysis (Default: 8).

Optional args:
--overwrite: Flag to overwrite in case there is already a previous run of OrthogroupsAnnotation (Default: off).
-help: Display this help message.
```

### TEPrediction

```
Script to predict TEs in the sequences based on earlgrey approach and Mycomobilome database.
This script perform two steps:
1. Run earlgrey TE prediction.
2. Organize the results."

Syntax: StarClust $(basename -s .sh "$0" ) [ -help ] -w <directory_path> -d <file_path> [ -m <string> -t <integer> ]

Required args:
-w, --workingDirectory: Specify the working directory where all data are stored.

Required args with Default:
-d, --database: Mycomobilome database to be use (Default = allConsensus) [Available type: allConsensus, proteinEvidence, unknown].
-m, --mode: Define the data that will be use for the TE prediction. This can be perform for all the data or for each cluster (Default = Cluster) [Available mode: Cluster, All]."
-t, --threads: Number of threads for earlgrey (Default = 8)

Required args in 'Cluster' mode:
-c, --clusters: file with a list of clusters to be analyzed, each line correspond to a single cluster ID.

Optional args:
-help: Display this help message.
```

This script was specifically designed for the purpose of identifying putative Transposable Elements (TEs) located inside starships. It utilizes a recently developed, curated database of TEs sourced from fungi ([Mycomobilome](https://github.com/TobyBaril/MycoMobilome)). This command is designed to run independently and is not connected to, nor does it require, any of the other main workflow commands (except for the `Initialize` command, which is necessary to create the project folder structure).

## Project folder organization

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
|`CaptainIdentification`| `AllID` | `hmmscan`, `seqkit`, `blast+`  | [hmmer](http://hmmer.org/), [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/) |
|`CaptainIdentification`| `Cluster`, `FullAll`| `hmmscan`, `seqkit`, `blast+`, `macse`, `iqtree3`, `gotree` | [hmmer](http://hmmer.org/), [Shen et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38898985/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Ranwez et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30165589/), [Wong et al. 2025](https://ecoevorxiv.org/repository/view/8916/), [Lemoine & Gascuel 2021](https://pubmed.ncbi.nlm.nih.gov/34396097/) |
|`ClusterCharacterization`| - | `orthofinder`, `DIAMOND`, `blast+`, `ape`, `ggtree`, `gggenomes` | [Emms et al. 2025](https://www.biorxiv.org/content/10.1101/2025.07.15.664860v1), [Buchfink et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33828273/), [Camacho et al. 2009](https://pubmed.ncbi.nlm.nih.gov/20003500/), [Paradis et al. 2004](https://pubmed.ncbi.nlm.nih.gov/14734327/), [Yu et al. 2016](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12628), [Hackl et al. 2024](https://arxiv.org/abs/2411.13556) |
|`OrthogroupsAnnotation`| - | `foldseek`, `mafft`, `hhblits`, `interproscan`,| [van Kempen et al 2024](https://pubmed.ncbi.nlm.nih.gov/37156916/), [Katoh & Standley 2013](https://pubmed.ncbi.nlm.nih.gov/23329690/), [Steinegger et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31521110/), [Jones et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24451626/) |
|`TEPrediction`| - | `earlgrey`, `Micomobilome` | [Baril et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38577785/), [Baril & Croll 2025](https://www.biorxiv.org/content/10.1101/2025.10.28.685023v1)  |

## License

Waiting for License decision