# Introduction
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/3d67b9cc-48d5-4065-a350-f7d76e8a6a43)
STELLARSCOPE SETUP. (A) Alignments are filtered according to a user-provided list of passing barcodes ("whitelist"). (B) The cell barcode (CB) and unique molecular identifier (UMI) from valid fragments are stored internally. (C) Initial weight matrix with fragments as rows and candidate assigned features as columns. (D) Values for the initial weight matrix setup result from intersecting each fragment's alignment(s) with the TE features annotation and selecting the best alignment score for each fragment for each locus.
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/697f2268-d7de-4b3f-94a1-dd33e4994c5c)
MULTIMAPPER-AWARE UMI DEDUPLICATION. (E) fragments that contain the same CB+UMI combination (i.e. duplicates) and their alignment positions are identified. (F) An undirected weighted graph is built for each CB+UMI combination with fragments as nodes and shared alignments as edge weights. For each component the most informative read according to alignment quality and ambiguity criteria is selected as representative. This method identifies and corrects non obvious duplicates (e.g. f1-f2, and f1-f3).
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/112230d4-bdb3-4a15-b202-7786e7f739ca)
MODEL FITTING. Stellarscope fits a Bayesian mixture model to the deduplicated weight matrix using an expectation maximization algorithm for each cell (G), for all cells (H), and for each cell type (I) in pooling modes Individual, Pseudobulk, and Celltype, respectively.!
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/e385fdde-4331-4f6e-b228-2593da57af65)
REASSIGNMENT. Once the model is fitted and parameters are estimated, Stellarscope uses the posterior probability matrix to reassign ambiguous fragments to their final generating locus. Stellarscope provides a variety of reassignment strategies (J) including filtering based on a threshold, excluding fragments with multiple optimal alignments, and randomly selecting from multiple optimal alignments; these criteria result in a different number of excluded alignments (shaded in grey). The output from Stellarscope (K) includes an umi-tracking file with the graphs and representative reads selection; a log file with the fitted models, the number of observations and parameters estimated, and a log likelihood for the fitted model; an updated BAM file; and a sparse single-cell counts matrix compatible with all the generally used analysis tools.
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/05db4ccd-b0d2-47e2-8255-cff19705cd01)

## Brief Overview of Stellarscope

Stellarscope is a computational biology tool for quantifying
Transposable Elements (TEs) transcripts -at single-locus resolution- in
single cell RNA sequencing data (e.g. scRNA-seq, the RNA-seq portion of
CITE-seq or Immune Profiling, etc). This tutorial introduces the usage
of Stellarscope using a freely available data set of Peripheral Blood
Mononuclear Cells (PBMCs) from a healthy donor as an example.

- The input to Stellarscope is a single-cell BAM alignment file and a
  list of cell barcodes.
- Additionally, Stellarscope requires a TE transcripts annotation file
  in GTF format.
- The output of Stellarscope includes an updated BAM file and a sparse
  single-cell counts matrix compatible with various downstream analysis
  tools.

Single-cell data sets often contain reads originating from transposable element (TE) transcripts, and they are commonly discarded due to their tendency to map to multiple locations in the genome. Stellarscope addresses the challenge of determining their true origin by reassigning reads to their most probable origin using a a bayesian mixture model and an Expectation-Maximization (EM) algorithm.

Understanding how Stellarscope reassigns reads from transposable element (TE) transcripts using a Bayesian mixture model and an Expectation-Maximization (EM) algorithm involves several steps and concepts. Let's break it down:

### 1. Transposable Elements (TEs) in Single-Cell RNA Sequencing
	• TE Transcripts: TEs are sequences in the genome that can move around and often have multiple copies. When TEs are transcribed, their RNA can map to many places in the genome, making it hard to determine their exact origin.
	• Challenge: In single-cell RNA sequencing (scRNA-seq), reads from TEs are often discarded because they can map to multiple locations, complicating downstream analyses.
### 2. Bayesian Mixture Model
	• Bayesian Approach: A Bayesian mixture model is a probabilistic model that assumes data points (in this case, sequencing reads) are generated from a mixture of several distributions, each representing a different source or origin.
	• Components: For TE reads, each component of the mixture model represents a potential origin (genomic location) where the read might have come from.
	• Prior Information: Bayesian models incorporate prior knowledge or assumptions about the distribution of reads across different locations.
### 3. Expectation-Maximization (EM) Algorithm
	• Goal: The EM algorithm is used to find the most probable origins of the TE reads by iteratively refining the estimates of the parameters of the mixture model.
	• Two Steps:
		• Expectation (E) Step: Calculate the expected value of the latent variables (i.e., the probabilities that each read originated from each possible location) based on the current parameter estimates.
		• Maximization (M) Step: Update the parameters of the model (e.g., the probabilities of reads originating from each location) to maximize the likelihood given the expected values from the E step.
### How Stellarscope Works
	1. Initialization: Start with initial guesses for the parameters of the mixture model, such as the initial probabilities of reads coming from different genomic locations.
	2. Iterative Refinement:
		• E Step: For each read, calculate the probability that it came from each possible location using the current parameters. This involves considering the likelihood of the read mapping to each location and the prior probability of that location.
		• M Step: Update the parameters (e.g., location probabilities) to maximize the overall likelihood of observing the read data given the assignments from the E step.
	• Convergence: Repeat the E and M steps until the parameters converge to stable values, meaning the model has found the most probable assignments of reads to their origins.
### Example Simplified
	1. Reads: Suppose you have a read that could map to three different locations.
	2. Prior Probabilities: Assume initial probabilities for each location (e.g., 0.3, 0.3, 0.4).
	3. E Step: Calculate the likelihood of the read mapping to each location, and use this to update the probabilities (e.g., after the E step, the probabilities might be 0.2, 0.2, 0.6).
	4. M Step: Update the model parameters to maximize the likelihood based on the new probabilities.
	5. Repeat: Continue iterating until the probabilities don't change much, indicating that the most likely origin for the read has been found.
### Conclusion
Stellarscope uses a Bayesian mixture model and the EM algorithm to tackle the problem of ambiguous read mapping for TEs in scRNA-seq data. By iteratively refining the probable origins of reads, it ensures that reads are assigned to their most likely genomic locations, improving the accuracy and utility of the data for downstream analysis. This approach allows researchers to include TE reads in their analyses, providing a more complete picture of the transcriptome.

## Importance

Single-cell data sets often contain reads originating from transposable
element (TE) transcripts, and they are commonly discarded due to their
tendency to map to multiple locations in the genome.

Stellarscope addresses the challenge of determining their true origin by
reassigning reads to their most probable origin using a a bayesian
mixture model and an Expectation-Maximization (EM) algorithm.

By facilitating the quantification of TE transcripts in single-cell
analyses, Stellarscope contributes to advancing our understanding of
biological processes at the single-cell level.

[Learn more about
Stellarscope](https://www.biorxiv.org/content/10.1101/2023.12.28.573568v1.full)

# Requirements for this tutorial

- Linux (basic command line knowledge)
  - [curl](https://github.com/curl/curl)
- Internet connection

This tutorial is designed for beginners, so it includes information
about the sequencing files and instructions on obtaining basic
bioinformatics software, such as conda, which are commonly used in the
field. Experienced users may find the contents index helpful for swiftly
navigating through the tutorial and determining which steps to bypass.

# Setup

## Step 1: Install Miniconda

Download and use the
[Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/)
installer script:

``` bash
# download installer
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# run installer script
bash Miniconda3-latest-Linux-x86_64.sh
```

- Follow the prompts on the installer screens (press enter to continue,
  review and accept the license agreement, select a directory where the
  miniconda installation will live).

- You can now clean up (remove) the installer
  (`rm Miniconda3-latest-Linux-x86_64.sh`).

## Step 2: Prepare Working Directory

``` bash
# First create a directory `stellarscope_analysis` in a location of your preference
mkdir stellarscope_analysis

# then change into that directory
cd stellarscope_analysis
```

### Obtain scRNA-seq data

In this demonstration of Stellarscope’s usage, you will use single-cell
RNA-seq data obtained from 500 PBMCs from a healthy donor. The data was
generated by 10x Genomics using the Chromium X platform (“[500 Human
PBMCs, 3’ LT v3.1, Chromium
X](https://www.10xgenomics.com/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-x-3-1-low-6-1-0)”).

- Download the freely available raw reads (fastq files) from 10x
  Genomics. This should take a couple of minutes (depending on your
  connection)

``` bash
# use curl to download data
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_fastqs.tar
```

- Extract the contents of the tar file

``` bash
tar -xf 500_PBMC_3p_LT_Chromium_X_fastqs.tar
```

- After running this command, you will have a directory called
  `500_PBMC_3p_LT_Chromium_X_fastqs` in the working directory with the
  FASTQ files inside:

``` bash
# see the contents of the `500_PBMC_3p_LT_Chromium_X_fastqs` directory
tree 500_PBMC_3p_LT_Chromium_X_fastqs
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/b1b804ba-3adb-433a-9708-4127b042c731)

``` markdown
500_PBMC_3p_LT_Chromium_X_fastqs
|-- 500_PBMC_3p_LT_Chromium_X_S4_L003_I1_001.fastq.gz
|-- 500_PBMC_3p_LT_Chromium_X_S4_L003_I2_001.fastq.gz
|-- 500_PBMC_3p_LT_Chromium_X_S4_L003_R1_001.fastq.gz
|-- 500_PBMC_3p_LT_Chromium_X_S4_L003_R2_001.fastq.gz
|-- 500_PBMC_3p_LT_Chromium_X_S4_L004_I1_001.fastq.gz
|-- 500_PBMC_3p_LT_Chromium_X_S4_L004_I2_001.fastq.gz
|-- 500_PBMC_3p_LT_Chromium_X_S4_L004_R1_001.fastq.gz
`-- 500_PBMC_3p_LT_Chromium_X_S4_L004_R2_001.fastq.gz
```

In this case, the reads in the `R1` files are 28 nucleotides long, the
first 16 nt correspond to the cell barcode (CB) sequence and the rest
(10 nt) to the unique molecular identifier (UMI) sequence; these lengths
are specific to the 3’ assay v3 chemistry 10x protocol. The reads in
`R2` files correspond to the sequenced cDNA (i.e. generated from the
transcripts). The `I1` and `I2` files contain Illumina multiplexing
indexes and won’t be used in this protocol.

``` bash
(base) liumy@shpc-3412-instance-ANJXbtSK:~/software/stellarscope$ zcat 500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L003_I1_001.fastq.gz | head
@A00519:986:HFFLJDSX2:3:1101:11143:1000 1:N:0:GAGACGCACG+ATGTTCATAG
GAGACGCACG
+
FFFFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:12970:1000 1:N:0:GAGACGCACG+ATGTTCATAG
GAGACGCACG
+
FFFFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:13729:1000 1:N:0:GAGACGCACG+ATGTTCATAG
GAGACGCACG


(base) liumy@shpc-3412-instance-ANJXbtSK:~/software/stellarscope$ zcat 500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L003_R1_001.fastq.gz | head
@A00519:986:HFFLJDSX2:3:1101:11143:1000 1:N:0:GAGACGCACG+ATGTTCATAG
ANTAACCTCTCCGAGGAATAGCTAAGCA
+
F#FFFFFFFFFFFFFFFFFFFFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:12970:1000 1:N:0:GAGACGCACG+ATGTTCATAG
ANCCTTGCAAGTAGTACAAATCGGTGTC
+
F#FFFFFFFFFFFFFFFFFFFFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:13729:1000 1:N:0:GAGACGCACG+ATGTTCATAG
TNACAGGAGAAATTGCATATGCTCAGAG


(base) liumy@shpc-3412-instance-ANJXbtSK:~/software/stellarscope$ zcat 500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L003_I2_001.fastq.gz | head
@A00519:986:HFFLJDSX2:3:1101:11143:1000 2:N:0:GAGACGCACG+ATGTTCATAG
ATGTTCATAG
+
:FFFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:12970:1000 2:N:0:GAGACGCACG+ATGTTCATAG
ATGTTCATAG
+
,:FFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:13729:1000 2:N:0:GAGACGCACG+ATGTTCATAG
ATGTTCATAG


(base) liumy@shpc-3412-instance-ANJXbtSK:~/software/stellarscope$ zcat 500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L003_R2_001.fastq.gz | head
@A00519:986:HFFLJDSX2:3:1101:11143:1000 2:N:0:GAGACGCACG+ATGTTCATAG
CTGCTGGAGGGCAGGGCCACAGGAGGGGTGAACAGGAGATGCCACAAGGTCACCATGAAGCACGACCTCAAATGATGTCATTTAGCTGGG
+
FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFF
@A00519:986:HFFLJDSX2:3:1101:12970:1000 2:N:0:GAGACGCACG+ATGTTCATAG
TAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F
@A00519:986:HFFLJDSX2:3:1101:13729:1000 2:N:0:GAGACGCACG+ATGTTCATAG
AATATTTATTGTCAAATCAGTGATTCTGTAGAGGGGTAAATAGCTAAATTTCCTTTTAAAAATTGTATGATGTGCTGATTTCAAATTGAT
```

### Obtain resources

We have prepared a resources bundle in the form of a directory that
contains the annotation files and some additional files to carry out the
Stellarscope analyses in this tutorial. We’ve made them available for
download through zenodo.

Use the curl command to download `resources.tar.gz` and untar to obtain
the `resources` directory. Please be aware that the download may take
some time as the file size is approximately 25 GB.

``` bash
# download file from zenodo
# curl -OJL https://zenodo.org/records/10671822/files/resources.tar.gz?download=1
# curl: (56) Recv failure: Connection timed out
wget -c -O resources.tar.gz https://zenodo.org/records/10671822/files/resources.tar.gz?download=1

# decompress file
tar -xvf resources.tar.gz
```

You should now see the `resources` directory with the following
contents:

``` bash
# the contents of `resources/`
tree resources
```

``` markdown
resources
|-- celltypes_tsv
|   |-- pbmc500_azimuth_l1.tsv
|   `-- pbmc500_azimuth_l2.tsv
|-- retro.hg38.v1.gtf
|-- STAR_GRCh38.d1.vd1_gencode.v38
|   `-- (contents not shown)
|-- stellarscope_protocol_env.yaml
`-- whitelist_10x
    `-- 3M-february-2018.txt
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/b63f03a3-61df-49b8-8551-2240925090bf)

## Step 3: Create a Stellarscope conda environment

Conda environments provide a way to encapsulate dependencies for
different analysis projects. Each self-contained environment can have
its own set of packages and dependencies, ensuring both reproducibility
and that they do not interfere with each other.

To go through this protocol, you will need [Python
≥3.5](https://www.python.org/), some additional libraries, the
[STAR](https://github.com/alexdobin/STAR) RNA-seq aligner, and
[Stellarscope](https://github.com/nixonlab/stellarscope). All of these
will be installed into a conda environment called
`stellarscope_protocol_env`.

- Use the yaml file we provide to create the environment:

``` bash
# Activate the conda base environment
conda activate base

# create a new conda environment named 'stellarscope_protocol_env'
conda env create --yes -f resources/stellarscope_protocol_env.yaml

# activate Stellarscope environment
conda activate stellarscope_protocol_env

# another method
conda create -n stellarscope
conda activate stellarscope 

conda install python samtools 'star==2.7.10b'
pip install future pyyaml cython numpy scipy 'pysam>=0.19' htslib intervaltree pandas packaging anndata 
git clone https://github.com/nixonlab/stellarscope.git
cd stellarscope
pip install -e .
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/e6973974-78f0-416c-9db5-ad48174fddba)

# Prepare Stellarscope’s input

## Step 1: STAR alignment

The first step is to align the raw reads from PBMCs to the reference
genome. The STAR aligner is used to complete this task. It is crucial
that the alignment includes multimapping reads (see STAR arguments
below), as these reads provide essential information such as alignment
quality, that is necessary for subsequent steps, i.e. the reassignment
process with Stellarscope. The output from this step will be a BAM file
`results/star_alignment/Aligned.sortedByCoord.out.bam`.

To expedite the process, 18 CPU cores are utilized. On average, with
this configuration, the alignment completes within ~20 minutes.

``` bash
# first create directory for STAR results
mkdir -p results/star_alignment

# run STAR alignment
STAR \
  --runThreadN 15 \
  --genomeDir resources/STAR_GRCh38.d1.vd1_gencode.v38 \
  --readFilesIn 500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L003_R2_001.fastq.gz,500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L004_R2_001.fastq.gz 500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L003_R1_001.fastq.gz,500_PBMC_3p_LT_Chromium_X_fastqs/500_PBMC_3p_LT_Chromium_X_S4_L004_R1_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --soloCBwhitelist resources/whitelist_10x/3M-february-2018.txt \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS NM nM MD CR CY UR UY CB UB GX GN sS sQ sM \
  --outSAMtype BAM SortedByCoordinate \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --limitOutSJcollapsed 5000000 \
  --outFilterMultimapNmax 500 \
  --outFilterMultimapScoreRange 5 \
  --outFileNamePrefix results/star_alignment/
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/84f14cda-b5ad-4b3f-80a0-8f0334d26086)

- The STAR reference genome index within the downloaded resources was
  built using the human reference genome GRCh38 version and GENCODE 38
  annotation (`--genomeDir`).

- Notice that STAR requires that the first FASTQ file(s) contain the
  cDNA read and the second FASTQ file(s) contain the cell barcode+UMI
  read, and that they’re separated by a blank space (`--readFilesIn`).

- In 10x Genomics single-cell sequencing, a cell barcode serves as a
  molecular “tag” that allows us to distinguish from which cell came the
  transcript that originated each sequencing read in the heterogeneous
  mixture that was sequenced. A
  [whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist)
  provides the list of known cell barcodes. The whitelist files are
  usually found as part of the
  [cellranger](https://github.com/10XGenomics/cellranger) software. For
  convenience, we have provided the correct whitelist in the downloaded
  resources (`--soloCBwhitelist`), if you use data from a different 10x
  protocol you will need to provide the appropriate whitelist.

- If you use data from a different protocol, with a different [barcode
  strategy/structure](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#barcode-geometry),
  you will need to adjust the `--soloCBstart`, `--soloCBlen`,
  `--soloUMIstart`, `--soloUMIlen` parameters.

- Refer to [STAR’s
  documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#matching-cellranger-4xx-and-5xx-results)
  for information on the adapter clipping, and the strategy to match
  whitelist barcodes, the UMI filtering and deduplication.

- In this alignment multimapping reads are allowed by setting the
  argument `--outFilterMultimapNmax` to `500`, and the argument
  `--outFilterMultimapScoreRange` is set to `5` so that for each
  multimapping read, its alignments additional to the best one will
  actually be included in the output BAM file.

### output

``` bash
# explore the contents of the star alignment directory
tree -t results/star_alignment
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/c70e58fc-409b-432a-9977-ef75925df643)

``` markdown
results/star_alignment/
|-- SJ.out.tab
|-- Solo.out
|   |-- Barcodes.stats
|   `-- Gene
|       |-- Features.stats
|       |-- raw
|       |   |-- features.tsv
|       |   |-- matrix.mtx
|       |   `-- barcodes.tsv
|       |-- filtered
|       |   |-- features.tsv
|       |   |-- matrix.mtx
|       |   `-- barcodes.tsv
|       |-- Summary.csv
|       `-- UMIperCellSorted.txt
|-- Aligned.sortedByCoord.out.bam
|-- Log.final.out
|-- Log.out
|-- Log.progress.out
`-- _STARtmp/
    `-- (contents not shown)
```

## Step 2: Stellarscope Cellsort

The second step to prepare Stellarscope’s input is sorting the alignment
BAM file so that all the cell barcodes bundled together, all the
alignments for each read are together, and the cell barcodes that are
not in the whitelist are filtered out.

This is achieved using the `stellarscope cellsort [samfile] [whitelist]`
command. Additionally, the argument `--outfile` is used to indicate the
output sorted BAM file (in this case it’ll be called
`results/stellarscope/Aligned.sortedByCB.bam`), and a `--tempdir`
(e.g. `/tmp`, most modern Unix-like systems, including Linux
distributions and macOS, typically include a /tmp directory as part of
their file system hierarchy).

Using 12 CPU cores, `stellarscope cellsort` completed within 13 minutes.

``` bash
# create directory for Stellarscope results
mkdir -p results/stellarscope

# sort reads alignments by CB
stellarscope cellsort \
  --nproc 15 \
  --tempdir /tmp \
  --outfile results/stellarscope/Aligned.sortedByCB.bam \
  results/star_alignment/Aligned.sortedByCoord.out.bam \
  results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv 
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/488494e2-9e66-49e6-9519-c4c0e11529eb)
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/0c9da6b1-d31e-464b-a718-40b9d000a1dc)

##error
``` bash
pip uninstall stellarscope
pip install -e .
```

# Basic Stellarscope analysis

The basic usage of Stellarscope is
`stellarscope assign [samfile] [gtffile]`. Where the `samfile` is the
alignment containing multimapping reads (e.g. the CB-sorted
`results/stellarscope/Aligned.sortedByCB.bam` file), and `gtffile` is
the TE annotation in GTF format.

The file `resources/retro.hg38.v1.gtf` defines the transcriptional unit
of each TE to be quantified, for more details see the Stellarscope
paper.

``` bash
# you can use this command to explore the annotation
head resources/retro.hg38.v1.gtf 
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/3da4f926-800b-4cdc-846f-4b59b4a7ce13)

The `stellarscope assign` command will do the following operations:
filter alignments by the user-provided list of passing barcodes
(`--whitelist`), create an initial weight matrix with fragments and
candidate assigned features, perform our multimapper-aware UMI
deduplication, fit a bayesian model to the deduplicated weight matrix
using an EM, once the EM is complete, the posterior probability matrix
is used to reassign each ambiguous fragment to their generating locus.

Outputs from `stellarscope assign` include a sparse counts matrix with
values for the TEs in the annotation, a stats file with a report of the
run, and (if requested) an updated BAM. Additionally, at specific
checkpoints throughout Stellarscope’s execution, Python pickle files are
produced. These files mark the operations described above ensuring that
key stages of the process are preserved and enabling easy retrival and
continuation in case of interruptions, or for future reference.

In this Stellarscope analysis, the argument `--pooling_mode` is set to
`individual`, thus fitting one model for each cell barcode (i.e
resolving ambiguous alignments to TEs within each cell).

Using 12 CPU cores, `stellarscope assign` completes within ~1 hour.

``` bash
# create a directory for stellarscope individual results
mkdir -p results/stellarscope/individual

# stellarscope analysis (pooling mode individual)
nohup stellarscope assign \
  --exp_tag pbmc500_individual \
  --outdir results/stellarscope/individual \
  --nproc 15 \
  --stranded_mode F \
  --whitelist results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
  --pooling_mode individual \
  --reassign_mode best_exclude \
  --max_iter 500 --debug --seed 240626 \
  --updated_sam \
  results/stellarscope/Aligned.sortedByCB.bam \
  resources/retro.hg38.v1.gtf > results/stellarscope/log/stellarscope_assign_individual_240626.log 2>&1 &
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/3ef8fb05-4b03-4a4e-9328-b3bed1eb873e)

set seed for stellarscope assign

2h30min 

- The string provided to `--exp_tag` sets the basename for all the
  Stellarscope output files.

- The argument `--stranded_mode` is used to consider feature strand when
  assigning reads, and here it is set to F as [10x libraries are
  stranded](https://kb.10xgenomics.com/hc/en-us/articles/360004396971-Are-10x-Single-Cell-gene-expression-libraries-strand-specific)
  (you can find more information on each protocol’s sheet).

- In this case the `--reassign_mode` argument is set to
  “`best_exclude`”. This means that fragments with multiple best
  assignments (i.e. their values in the posterior probability matrix are
  equal) are excluded from the final counts.

- The argument `--max_iter` controls the maximum number of iterations
  for the EM algorithm.

- The argument `--updated_sam` is useful to obtain an updated alignment
  file (`results/stellarscope/pbmc500_individual-updated.bam`)

These are all the outputs in the `results/stellarscope/individual`
directory:

``` bash
# inspect the analysis results
tree -t results/stellarscope
```

``` markdown
results/stellarscope/
|-- Aligned.sortedByCB.bam
`-- individual
    |-- pbmc500_individual-other.bam
    |-- pbmc500_individual-tmp_tele.bam
    |-- pbmc500_individual-checkpoint.load_alignment.pickle
    |-- pbmc500_individual-umi_tracking.txt
    |-- pbmc500_individual-checkpoint.dedup_umi.pickle
    |-- pbmc500_individual-barcodes.tsv
    |-- pbmc500_individual-features.tsv
    |-- pbmc500_individual-TE_counts.mtx
    |-- pbmc500_individual-stats.final.tsv
    |-- pbmc500_individual-updated.bam
    `-- pbmc500_individual-checkpoint.final.pickle
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/a7a27142-e6de-45de-9823-c0dfb812862e)


You can see the operations that `stellarscope assign` performs mirrored
in the output files. First, `pbmc500_individual-other.bam`,
`pbmc500_individual-tmp_tele.bam`, and
`pbmc500_individual-checkpoint.load_alignment.pickle` are generated.
These BAM files contain the alignments that do and do not overlap the TE
annotation, respectively; the pickle file corresponds to Stellarscope’s
loading of alignments that intersect with the TE annotation. Second, the
files `pbmc500_individual-umi_tracking.txt` and
`pbmc500_individual-checkpoint.dedup_umi.pickle` are created. They
contain a record and checkpoint of the PCR duplicates removal process.
Next, Stellarscope outputs the sparse matrix files (`*-barcodes.tsv`,
`*-features.tsv`, and `*-TE_counts.mtx`). Finally, the stats report file
(`*-stats.final.tsv`), the updated BAM (`*-updated.bam`) and final
checkpoint (`*-checkpoint.final.pickle`) are saved.

# Advanced Stellarscope analysis: multiple pooling modes

For certain analyses, the need to conduct Stellarscope’s reassigning
using multiple pooling modes may arise. In this scenario it will be
beneficial to start by loading the data and generating the checkpoint
file corresponding to the UMI deduplication stage. This will provide the
flexibility to employ different pooling modes and arguments,
facilitating reassignment procedures.

You can achieve this by combining the commands `stellarscope assign`
with the flag `--skip_em` and `stellarscope resume`.

## Step 1: Stellarscope load

In this case (`--skip_em`), there is no pooling mode argument, no
reassign mode argument, and no iterations argument.

Stellarscope will load the alignments intersecting TE annotation and
remove PCR duplicates in around 45 minutes.

``` bash
# load the alignment data into stellarscope
nohup stellarscope assign \
  --exp_tag pbmc500_stload \
  --outdir results/stellarscope \
  --skip_em \
  --stranded_mode F \
  --whitelist results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
  --updated_sam --seed 240626 \
  results/stellarscope/Aligned.sortedByCB.bam \
  resources/retro.hg38.v1.gtf \
  --logfile results/stellarscope/pbmc500_stload.log > results/stellarscope/log/stellarscope_assign_stload_240626.log 2>&1 &
```
[18600]

Upon inspection of the `results/stellarscope` directory, you can see
that the two checkpoints
`pbmc500_stload-checkpoint.load_alignment.pickle` and
`pbmc500_stload-checkpoint.dedup_umi.pickle` have been created.

In this case, the argument `--logfile` was provided, and thus, the file
`results/stellarscope/pbmc500_stload.log` is created. This file captures
all the information that would have been displayed in the terminal
during the Stellarscope’s execution. It includes messages, warnings, and
errors (if any) generated as it ran.

``` bash
# the stellarscope files are output with the prefix `pbmc500_stload`
tree -t results/stellarscope
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/d1779b53-6a1b-4a31-bce2-5496011020db)

``` markdown
results/stellarscope/
|-- Aligned.sortedByCB.bam
|-- individual/
|   `-- (contents not shown)
|-- pbmc500_stload-other.bam
|-- pbmc500_stload-tmp_tele.bam
|-- pbmc500_stload-checkpoint.load_alignment.pickle
|-- pbmc500_stload-umi_tracking.txt
|-- pbmc500_stload-checkpoint.dedup_umi.pickle
|-- pbmc500_stload.log
`-- pbmc500_stload-stats.final.tsv
```

## Step 2: Stellarscope resume

Now, two Stellarscope analyses will demonstrate how to use the
`stellarscope_resume` command.

- In the first one, `--pooling_mode` pseudobulk is used, and matrices
  employing every reassignment mode are generated.

- In the second one, cell type labels from annotating the PBMCs using
  the reference transcriptomes from
  [azimuth](https://azimuth.hubmapconsortium.org/) are used to perform
  Stellarscope’s reassignment with `--pooling_mode` celltype, and the
  final counts matrix is generated by setting a threshold and
  considering assignments that exceed it.

### Stellarscope Pseudobulk

In the pooling mode `pseudobulk` a single model is fitted, meaning that
ambiguous alignments are resolved sharing information across all cells.

For this analysis, the `stellarscope resume` command accepts
reassignment arguments (pooling mode, iterations number, etc.) and a
`*-checkpoint.dedup_umi.pickle` file (instead of a BAM file).

The reassignment mode is the method that Stellarscope will use to
resolve which one is the “best” reassignment for fragments that have
multiple possible reassignments. In this case, the argument
`--use_every_reasign_mode` is added so that Stellarscope outputs seven
`.mtx` files in total.

One for each `best_exclude` (fragments with multiple best assignments
are excluded), `best_conf` (assignments with values above a threshold
set by `--conf_prob` are counted), `best_random` (the best assignment is
randomly chosen from the best assignments set), `best_average` (the
fragment is divided evenly among the best assignments), `initial_unique`
(only the originally uniquely aligned reads are included; EM is not
considered), `initial_random` (an assignment is randomly selected from
the initial alignments; EM is not considered), and `total_hits` (all
assignments are included; EM is not considered). The final three modes
do not perform reassignment or model fitting but are included for
comparison.

Note that the matrix corresponding to `best_exclude` will be called
`pbmc500_pseudobulk-TE_counts.mtx`, as opposed to the other matrices,
which will have the reassignment mode indicated in their name
(e.g. `pbmc500_pseudobulk-TE_counts.best_conf.mtx`). This is because
`best_exclude` is the default option.

This Stellarscope resume analysis took less than 2 minutes in total.

``` bash
# create a directory for stellarscope pseudobulk results
mkdir -p results/stellarscope/pseudobulk

# use stellarscope resume to continue from the deduplication checkpoint
stellarscope resume \
            --exp_tag pbmc500_pseudobulk \
            --outdir results/stellarscope/pseudobulk \
            --nproc 12 \
            --pooling_mode pseudobulk \
            --use_every_reassign_mode \
            --max_iter 500\
            --updated_sam \
            results/stellarscope/pbmc500_stload-checkpoint.dedup_umi.pickle \
            --logfile results/stellarscope/pseudobulk/pbmc500_pseudobulk.log
```

You can see `results/stellarscope/pseudobulk/` directory that you obtain
a matrix for each reassignment mode:

``` bash
# now the resutls from pooling mode pseudobulk are in the `pseudobulk` directory
tree -t results/stellarscope
```

``` markdown
results/stellarscope/
|-- Aligned.sortedByCB.bam
|-- individual/
|   `-- (contents not shown)
|-- pbmc500_stload-other.bam
|-- pbmc500_stload-tmp_tele.bam
|-- pbmc500_stload-checkpoint.load_alignment.pickle
|-- pbmc500_stload-umi_tracking.txt
|-- pbmc500_stload-checkpoint.dedup_umi.pickle
|-- pbmc500_stload.log
|-- pbmc500_stload-stats.final.tsv
`-- pseudobulk/
    |-- pbmc500_pseudobulk-barcodes.tsv
    |-- pbmc500_pseudobulk-features.tsv
    |-- pbmc500_pseudobulk-TE_counts.mtx
    |-- pbmc500_pseudobulk-TE_counts.best_conf.mtx
    |-- pbmc500_pseudobulk-TE_counts.best_random.mtx
    |-- pbmc500_pseudobulk-TE_counts.best_average.mtx
    |-- pbmc500_pseudobulk-TE_counts.initial_unique.mtx
    |-- pbmc500_pseudobulk-TE_counts.initial_random.mtx
    |-- pbmc500_pseudobulk-TE_counts.total_hits.mtx
    |-- pbmc500_pseudobulk-stats.final.tsv
    |-- pbmc500_pseudobulk-checkpoint.final.pickle
    `-- pbmc500_pseudobulk.log
```

### Stellarscope Celltype

This is an example of fitting one model per cell type. This strategy
depends on user-provided cell identities. Using this `celltype` pooling
mode carries the assumption that the (retro)transcriptome is shared more
closely between cells that are classified as the same cell type.
Annotating/solving cell identity is an entire research line, so users
should be aware of their goals and the strategies they use to determine
cell types in their data.

To exemplify this Stellarscope analysis, we have provided within the
resources a file `resources/celltypes_tsv/pbmc500_azimuth_l1.tsv` where
each of the PBMC’s cell barcodes has an assigned cell types obtained
from azimuth. This is provided to Stellarscope using the argument
`--celltype_tsv`.

This Stellarscope analysis took around 7 minutes in total.

``` bash
mkdir -p results/stellarscope/celltype

stellarscope resume \
            --exp_tag pbmc500_l1 \
            --outdir results/stellarscope/celltype \
            --nproc 12 \
            --pooling_mode celltype \
            --celltype_tsv resources/celltypes_tsv/pbmc500_azimuth_l1.tsv \
            --reassign_mode best_conf \
            --conf_prob 0.95 \
            --max_iter 500\
            --updated_sam \
            results/stellarscope/pbmc500_stload-checkpoint.dedup_umi.pickle \
            --logfile results/stellarscope/celltype/pbmc500_l1.log
```

Inspect the outputs in the `results/stellarscope/celltype` directory:

Since only one reassignment mode was requested (`best_conf` with a
`0.95` threshold), only one `.mtx` file is produced, and the
reassignment mode is not indicated in its name, it is simply called
`results/stellarscope/celltype/pbmc500_l1-TE_counts.mtx`.

``` bash
# stellarscope pooling mode celltype results
tree -t results/stellarscope
```

``` markdown
results/stellarscope/
|-- Aligned.sortedByCB.bam
|-- individual/
|   `-- (contents not shown)
|-- pbmc500_stload-other.bam
|-- pbmc500_stload-tmp_tele.bam
|-- pbmc500_stload-checkpoint.load_alignment.pickle
|-- pbmc500_stload-umi_tracking.txt
|-- pbmc500_stload-checkpoint.dedup_umi.pickle
|-- pbmc500_stload.log
|-- pbmc500_stload-stats.final.tsv
|-- pseudobulk/
|   `-- (contents not shown)
`-- celltype/
    |-- pbmc500_l1-barcodes.tsv
    |-- pbmc500_l1-features.tsv
    |-- pbmc500_l1-TE_counts.mtx
    |-- pbmc500_l1-stats.final.tsv
    |-- pbmc500_l1-checkpoint.final.pickle
    `-- pbmc500_l1.log
```
