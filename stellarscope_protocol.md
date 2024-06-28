# 1. Introduction
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/3d67b9cc-48d5-4065-a350-f7d76e8a6a43)
**STELLARSCOPE SETUP**. (A) Alignments are filtered according to a user-provided list of passing barcodes ("whitelist"). (B) The cell barcode (CB) and unique molecular identifier (UMI) from valid fragments are stored internally. (C) Initial weight matrix with **fragments as rows and candidate assigned features as columns**. (D) Values for the initial weight matrix setup result from **intersecting** each fragment's alignment(s) with the TE features annotation and selecting the best alignment score for each fragment for each locus.

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/697f2268-d7de-4b3f-94a1-dd33e4994c5c)
**MULTIMAPPER-AWARE UMI DEDUPLICATION**. (E) fragments that contain the same CB+UMI combination (i.e. duplicates) and their alignment positions are identified. (F) An undirected weighted graph is built for each CB+UMI combination with fragments as nodes and shared alignments as edge weights. For each component the most informative read according to alignment quality and ambiguity criteria is selected as representative. This method identifies and corrects non obvious duplicates (e.g. f1-f2, and f1-f3).

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/112230d4-bdb3-4a15-b202-7786e7f739ca)
**MODEL FITTING**. Stellarscope fits a **Bayesian mixture model to the deduplicated weight matrix** using an **expectation maximization algorithm** for each cell (G), for all cells (H), and for each cell type (I) in pooling modes Individual, Pseudobulk, and Celltype, respectively.!

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/e385fdde-4331-4f6e-b228-2593da57af65)
**REASSIGNMENT**. Once the model is fitted and parameters are estimated, Stellarscope uses the **posterior probability matrix** to reassign ambiguous fragments to their final generating locus. Stellarscope provides a variety of reassignment strategies (J) including filtering based on a threshold, excluding fragments with multiple optimal alignments, and randomly selecting from multiple optimal alignments; these criteria result in a different number of excluded alignments (shaded in grey). The output from Stellarscope (K) includes an **umi-tracking file** with the graphs and representative reads selection; **a log file** with the fitted models, the number of observations and parameters estimated, and a log likelihood for the fitted model; an **updated BAM file**; and **a sparse single-cell counts matrix** compatible with all the generally used analysis tools.

## 1.1 Brief Overview of Stellarscope

Stellarscope is a computational biology tool for **quantifying Transposable Elements (TEs) transcripts** -at single-locus resolution- in single cell RNA sequencing data (e.g. scRNA-seq, the RNA-seq portion of CITE-seq or Immune Profiling, etc). **This tutorial introduces the usage of Stellarscope using a freely available data set of Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor as an example.**

- The input to Stellarscope is **a single-cell BAM alignment file and a list of cell barcodes**.
- Additionally, Stellarscope requires a **TE transcripts annotation file in GTF format**.
- The output of Stellarscope includes **an updated BAM file and a sparse single-cell counts matrix** compatible with various downstream analysis tools.

Single-cell data sets often contain reads originating from transposable element (TE) transcripts, and they are commonly discarded due to their tendency to map to multiple locations in the genome. Stellarscope addresses the challenge of **determining their true origin by reassigning reads to their most probable origin using a a bayesian mixture model and an Expectation-Maximization (EM) algorithm**.

> ### Stellarscope reassigns reads from transposable element (TE) transcripts using a Bayesian mixture model and an Expectation-Maximization (EM) algorithm
> #### 1. Transposable Elements (TEs) in Single-Cell RNA Sequencing
>> • TE Transcripts: TEs are sequences in the genome that can move around and often have multiple copies. When TEs are transcribed, their RNA can map to many places in the genome, making it hard to determine their exact  origin.
>> 
>> • Challenge: In single-cell RNA sequencing (scRNA-seq), reads from TEs are often discarded because they can map to multiple locations, complicating downstream analyses.
> #### 2. Bayesian Mixture Model
>> • Bayesian Approach: A Bayesian mixture model is a probabilistic model that assumes data points (in this case, sequencing reads) are generated from a mixture of several distributions, each representing a different source or origin.
>> 
>> • Components: For TE reads, each component of the mixture model represents a potential origin (genomic location) where the read might have come from.
>> 
>> • Prior Information: Bayesian models incorporate prior knowledge or assumptions about the distribution of reads across different locations.
> #### 3. Expectation-Maximization (EM) Algorithm
>> • Goal: The EM algorithm is used to find the most probable origins of the TE reads by iteratively refining the estimates of the parameters of the mixture model.
>> 
>> • Two Steps: 1) Expectation (E) Step: Calculate the expected value of the latent variables (i.e., the probabilities that each read originated from each possible location) based on the current parameter estimates. 2) Maximization (M) Step: Update the parameters of the model (e.g., the probabilities of reads originating from each location) to maximize the likelihood given the expected values from the E step.

## 1.2 Importance

By facilitating the quantification of TE transcripts in single-cell analyses, Stellarscope contributes to advancing our understanding of biological processes at the single-cell level.

[Learn more about Stellarscope: A single-cell transposable element atlas of human cell identity](https://www.biorxiv.org/content/10.1101/2023.12.28.573568v1.full)

## 1.3 Requirements for this tutorial
This tutorial is designed for beginners, so it includes information about the sequencing files and instructions on obtaining basic bioinformatics software, such as conda, which are commonly used in the field. Experienced users may find the contents index helpful for swiftly navigating through the tutorial and determining which steps to bypass.
- Linux (basic command line knowledge): [curl](https://github.com/curl/curl)
- Internet connection

# 2 Setup
## 2.1 Install Miniconda

Download and use the [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/) installer script:

``` bash
# download installer
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# run installer script
bash Miniconda3-latest-Linux-x86_64.sh
```

- Follow the prompts on the installer screens (press enter to continue, review and accept the license agreement, select a directory where the miniconda installation will live).
- You can now clean up (remove) the installer (`rm Miniconda3-latest-Linux-x86_64.sh`).

## 2.2 Prepare Working Directory

``` bash
# First create a directory `stellarscope` in a location of your preference
mkdir ~/software/stellarscope

# then change into that directory
cd ~/software/stellarscope
```

### 2.2.1 Obtain scRNA-seq data

In this demonstration of Stellarscope’s usage, you will use single-cell RNA-seq data obtained from 500 PBMCs from a healthy donor. The data was generated by 10x Genomics using the Chromium X platform (“[500 Human
PBMCs, 3’ LT v3.1, Chromium X](https://www.10xgenomics.com/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-x-3-1-low-6-1-0)”). Download the freely available raw reads (fastq files) from 10x Genomics. This should take a couple of minutes (depending on your connection)

``` bash
# use curl to download data
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_fastqs.tar

# Extract the contents of the tar file
tar -xf 500_PBMC_3p_LT_Chromium_X_fastqs.tar

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

> #### file check: 500_PBMC_3p_LT_Chromium_X_fastqs
> In this case, the reads in the `R1` files are 28 nucleotides long, the first 16 nt correspond to the cell barcode (CB) sequence and the rest (12 nt) to the unique molecular identifier (UMI) sequence; these lengths
are specific to the 3’ assay v3 chemistry 10x protocol. The reads in `R2` files correspond to the sequenced cDNA (i.e. generated from the transcripts). The `I1` and `I2` files contain Illumina multiplexing indexes and won’t be used in this protocol.

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

### 2.2.2 Obtain resources

We have prepared a resources bundle in the form of a directory that contains **the annotation files and some additional files** to carry out the Stellarscope analyses in this tutorial. We’ve made them available for download through zenodo. Use the curl command to download `resources.tar.gz` and untar to obtain the `resources` directory. Please be aware that the download may take some time as the file size is approximately 25 GB.

``` bash
# download file from zenodo
# curl -OJL https://zenodo.org/records/10671822/files/resources.tar.gz?download=1
# curl: (56) Recv failure: Connection timed out
wget -c -O resources.tar.gz https://zenodo.org/records/10671822/files/resources.tar.gz?download=1

# decompress file
tar -xvf resources.tar.gz

# the contents of `resources/`
tree resources
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/b63f03a3-61df-49b8-8551-2240925090bf)

``` markdown
resources
├── celltypes_tsv
│   ├── pbmc500_azimuth_l1.tsv # celltype barcode file
│   └── pbmc500_azimuth_l2.tsv
├── retro.hg38.v1.gtf  # TE annotation files
├── src
├── STAR_GRCh38.d1.vd1_gencode.v38
│   ├── chrLength.txt
│   ├── chrNameLength.txt
│   ├── chrName.txt
│   ├── chrStart.txt
│   ├── exonGeTrInfo.tab
│   ├── exonInfo.tab
│   ├── geneInfo.tab
│   ├── Genome
│   ├── genomeParameters.txt
│   ├── SA
│   ├── SAindex
│   ├── sjdbInfo.txt
│   ├── sjdbList.fromGTF.out.tab
│   ├── sjdbList.out.tab
│   └── transcriptInfo.tab
├── stellarscope_protocol_env_1.yaml 
├── stellarscope_protocol_env.yaml # environment.yaml file
└── whitelist_10x
    └── 3M-february-2018.txt # white list data

4 directories, 21 files
```

#### check the resource files
**`head resources/STAR_GRCh38.d1.vd1_gencode.v38/geneInfo.tab`**
``` markdown
60649
ENSG00000223972.5	DDX11L1	transcribed_unprocessed_pseudogene
ENSG00000227232.5	WASH7P	unprocessed_pseudogene
ENSG00000278267.1	MIR6859-1	miRNA
ENSG00000243485.5	MIR1302-2HG	lncRNA
ENSG00000284332.1	MIR1302-2	miRNA
ENSG00000237613.2	FAM138A	lncRNA
ENSG00000268020.3	OR4G4P	unprocessed_pseudogene
ENSG00000240361.2	OR4G11P	transcribed_unprocessed_pseudogene
ENSG00000186092.7	OR4F5	protein_coding
```

**`head resources/STAR_GRCh38.d1.vd1_gencode.v38/sjdbList.out.tab`**
``` markdown
chr1	12058	12178	+
chr1	12228	12612	+
chr1	12698	12974	+
chr1	12722	13220	+
chr1	13053	13220	+
chr1	13375	13452	+
chr1	14502	15004	-
chr1	15039	15795	-
chr1	15948	16606	-
chr1	16766	16857	-
```

**`head resources/STAR_GRCh38.d1.vd1_gencode.v38/sjdbList.fromGTF.out.tab`**
``` markdown
chr1	12058	12178	+	1
chr1	12228	12612	+	1
chr1	12698	12974	+	1
chr1	12722	13220	+	1
chr1	13053	13220	+	1
chr1	13375	13452	+	1
chr1	14502	15004	-	2
chr1	15039	15795	-	2
chr1	15948	16606	-	2
chr1	16766	16857	-	2
```

**`head resources/STAR_GRCh38.d1.vd1_gencode.v38/sjdbInfo.txt`**
``` markdown
390941	100
12057	12177	0	0	0	1
12227	12611	1	0	2	1
12697	12973	1	0	0	1
12721	13219	1	2	1	1
13052	13219	3	2	3	1
13374	13451	1	0	0	1
14500	15002	0	1	0	2
15038	15794	2	0	2	2
15947	16605	2	0	3	2
```

## 2.3 Create a Stellarscope conda environment

Conda environments provide a way to encapsulate dependencies for different analysis projects. Each self-contained environment can have its own set of packages and dependencies, ensuring both reproducibility and that they do not interfere with each other. To go through this protocol, you will need [Python ≥3.5](https://www.python.org/), some additional libraries, the [STAR](https://github.com/alexdobin/STAR) RNA-seq aligner, and [Stellarscope](https://github.com/nixonlab/stellarscope). All of these will be installed into a conda environment called `stellarscope_protocol_env`.

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

> #### file check: stellarscope_protocol_env.yaml
``` markdown
name: stellarscope_protocol_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python >=3.6
  - future
  - pip
  - pyyaml
  - cython
  - numpy
  - scipy
  - pysam >=0.19
  - htslib
  - intervaltree
  - pandas
  - samtools >=1.16
  - packaging
  - anndata
  - star ==2.7.10b
  - pip:
      - -e git+ssh://git@github.com/nixonlab/stellarscope.git@dev#egg=stellarscope

# change the last one if you can't access the ssh & no dev
#  - -e git+https://github.com/nixonlab/stellarscope.git@main#egg=stellarscope
```

> ### -e git+https://github.com/nixonlab/stellarscope.git@main#egg=stellarscope
> -e: This flag tells pip to install the package in "editable" mode. This means that changes to the source code will be reflected immediately without needing to reinstall the package.
> 
> git+https://github.com/nixonlab/stellarscope.git: This specifies the URL of the Git repository to clone.
> 
> @dev or @main: This specifies the branch to check out. In this case, it’s the dev branch.
> 
> #egg=stellarscope: This is a naming convention used by pip to name the package after installing it. The egg syntax is used to ensure that the package is named correctly in the environment.

# 3 Prepare Stellarscope’s input
## 3.1 STAR alignment

The first step is to align the raw reads from PBMCs to the reference genome. The **STAR aligner** is used to complete this task. It is crucial that the alignment includes multimapping reads (see STAR arguments below), as these reads provide essential information such as **alignment quality**, that is necessary for subsequent steps, i.e. the reassignment process with Stellarscope. The output from this step will be a BAM file `results/star_alignment/Aligned.sortedByCoord.out.bam`. 

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

>#### arguments check: STAR
>>- The STAR reference genome index within the downloaded resources was built using the human reference genome GRCh38 version and GENCODE 38 annotation (`--genomeDir`).
>>- Notice that STAR requires that the first FASTQ file(s) contain the cDNA read and the second FASTQ file(s) contain the cell barcode+UMI read, and that they’re separated by a blank space (`--readFilesIn`).
>>- In 10x Genomics single-cell sequencing, a cell barcode serves as a molecular “tag” that allows us to distinguish from which cell came the transcript that originated each sequencing read in the heterogeneous mixture that was sequenced. A [whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist) provides the list of known cell barcodes. The whitelist files are usually found as part of the [cellranger](https://github.com/10XGenomics/cellranger) software. For convenience, we have provided the correct whitelist in the downloaded resources (`--soloCBwhitelist`), if you use data from a different 10x   protocol you will need to provide the appropriate whitelist.
>>> ##### whitelist with 6794880 cell barcodes
>>> download form https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist
>>> 
>>> `head resources/whitelist_10x/3M-february-2018.txt`
>>> 
>>> ![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/2307530e-603e-4412-8a71-de8df3c35e56)
>>>
>>- If you use data from a different protocol, with a different [barcode strategy/structure](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#barcode-geometry), you will need to adjust the `--soloCBstart`, `--soloCBlen`, `--soloUMIstart`, `--soloUMIlen` parameters.
>>- Refer to [STAR’s documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#matching-cellranger-4xx-and-5xx-results) for information on the adapter clipping, and the strategy to match whitelist barcodes, the UMI filtering and deduplication.
>>- In this alignment multimapping reads are allowed by setting the argument `--outFilterMultimapNmax` to `500`, and the argument `--outFilterMultimapScoreRange` is set to `5` so that for each multimapping read, its alignments additional to the best one will actually be included in the output BAM file.
>>- Unmapped reads can be output into the SAM/BAM Aligned.* le(s) with `--outSAMunmapped Within` option.

### Output

``` bash
# explore the contents of the star alignment directory
tree -t results/star_alignment
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/c70e58fc-409b-432a-9977-ef75925df643)

``` markdown
results/star_alignment
├── SJ.out.tab
├── Solo.out
│   ├── Barcodes.stats
│   └── Gene
│       ├── Features.stats
│       ├── raw
│       │   ├── features.tsv
│       │   ├── barcodes.tsv
│       │   └── matrix.mtx
│       ├── filtered
│       │   ├── barcodes.tsv
│       │   ├── features.tsv
│       │   └── matrix.mtx
│       ├── Summary.csv
│       └── UMIperCellSorted.txt
├── Aligned.sortedByCoord.out.bam
├── Log.final.out 
├── Log.out # main log le with a lot of detailed information about the run. This le is most useful for troubleshooting and debugging.
└── Log.progress.out

4 directories, 15 files
```

### check the output
**`head results/star_alignment/SJ.out.tab`**

> SJ.out.tab contains high condence collapsed splice junctions in tab-delimited format. Note that STAR denes the junction start/end as intronic bases, while many other software dene them as exonic bases. The columns have the following meaning:
> - column 1: chromosome
> - column 2: first base of the intron (1-based)
> - column 3: last base of the intron (1-based)
> - column 4: strand (0: undened, 1: +, 2: -)
> - column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
> - column 6: 0: unannotated, 1: annotated in the splice junctions database. Note that in 2-pass mode, junctions detected in the 1st pass are reported as annotated, in addition to annotated junctions from GTF.
> - column 7: number of uniquely mapping reads crossing the junction
> - column 8: number of multi-mapping reads crossing the junction
> - column 9: maximum spliced alignment overhang

``` markdown
chr1	10131	10422	2	2	0	0	1	35
chr1	14441	17223	2	2	0	0	1	25
chr1	14830	14969	2	2	0	0	12	41
chr1	14830	185490	2	2	0	0	31	41
chr1	15039	15795	2	2	1	0	2	8
chr1	15948	16606	2	2	1	0	1	24
chr1	16766	16853	2	2	0	0	1	40
chr1	16766	16857	2	2	1	0	8	31
chr1	17056	17232	2	2	1	0	16	44
chr1	17056	187754	2	2	0	0	16	44
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/682a994a-c54e-434f-b1ac-c162bd9caf63)

**`cat results/star_alignment/Solo.out/Barcodes.stats`**

``` markdown
         noNoAdapter              0
             noNoUMI              0
              noNoCB              0
             noNinCB              0
            noNinUMI              0
    noUMIhomopolymer          27667
         noNoWLmatch        2002734
         noTooManyMM              0
  noTooManyWLmatches              0
     yesWLmatchExact       70976190
 yesOneWLmatchWithMM         547334
yesMultWLmatchWithMM        2116754
```

**`cat results/star_alignment/Solo.out/Gene/Summary.csv`**
``` markdown
Number of Reads,75670679
Reads With Valid Barcodes,0.971093
Sequencing Saturation,0.761387
Q30 Bases in CB+UMI,0.955356
Q30 Bases in RNA read,0.888375
Reads Mapped to Genome: Unique+Multiple,0.925104
Reads Mapped to Genome: Unique,0.675028
Reads Mapped to Gene: Unique+Multiple Gene,NoMulti
Reads Mapped to Gene: Unique Gene,0.385351
Estimated Number of Cells,564
Unique Reads in Cells Mapped to Gene,25755892
Fraction of Unique Reads in Cells,0.883268
Mean Reads per Cell,45666
Median Reads per Cell,35806
UMIs in Cells,5855011
Mean UMI per Cell,10381
Median UMI per Cell,7934
Mean Gene per Cell,2868
Median Gene per Cell,2413
Total Gene Detected,24212
```

**`cat results/star_alignment/Solo.out/Gene/Features.stats`**

``` markdown
                 noUnmapped        5102468
                noNoFeature       29214965
               MultiFeature       10006110
subMultiFeatureMultiGenomic        9238505
         noTooManyWLmatches         156976
       noMMtoWLwithoutExact              0
                 yesWLmatch       29159759
         yessubWLmatchExact       28466637
yessubWLmatch_UniqueFeature       29159759
            yesCellBarcodes          22909
                    yesUMIs        6957895
```

> - nNinBarcode: number of reads with more than 2 Ns in cell barcode (CB)
> - nUMIhomopolymer: number of reads with homopolymer in CB
> - nTooMany: not used at the moment
> - nNoMatch: number of reads with CBs that do not match whitelist even with one mismatch
> - All of the above reads are discarded from Solo output. Remaining reads are checked for overlap with features (e.g. genes):
> - nUnmapped: number of reads unmapped to the genome
> - nNoFeature: number of reads that map to the genome but do not belong to a feature
> - nAmbigFeature: number of reads that belong to more than one feature
> - nAmbigFeatureMultimap: number of reads that belong to more than one feature and are also multimapping to the genome (this is a subset of the nAmbigFeature)
> - nTooMany: number of reads with ambiguous CB (i.e. CB matches whitelist with one mismatch but with posterior probability <0.95)
> - nNoExactMatch: number of reads with CB that matches a whitelist barcode with 1 mis-match, but this whitelist barcode does not get any other reads with exact matches of CB
> - nCellBarcodes: number of distinct CBs detected
> - nUMIs: number of distinct UMIs detected


#### row file: features.tsv, barcodes.tsv, matrix.mtx - 60649 features with 6794880 cell barcodes
**`head results/star_alignment/Solo.out/Gene/raw/features.tsv`**

``` markdown
ENSG00000223972.5	DDX11L1	Gene Expression
ENSG00000227232.5	WASH7P	Gene Expression
ENSG00000278267.1	MIR6859-1	Gene Expression
ENSG00000243485.5	MIR1302-2HG	Gene Expression
ENSG00000284332.1	MIR1302-2	Gene Expression
ENSG00000237613.2	FAM138A	Gene Expression
ENSG00000268020.3	OR4G4P	Gene Expression
ENSG00000240361.2	OR4G11P	Gene Expression
ENSG00000186092.7	OR4F5	Gene Expression
ENSG00000238009.6	RP11-34P13.7	Gene Expression
```

**`head results/star_alignment/Solo.out/Gene/raw/barcodes.tsv`**

``` markdown
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCG
AAACCCAAGAAACCTG
AAACCCAAGAAACGAA
AAACCCAAGAAACGTC
AAACCCAAGAAACTAC
AAACCCAAGAAACTCA
AAACCCAAGAAACTGC
```

**`head results/star_alignment/Solo.out/Gene/raw/matrix.mtx`**

``` markdown
%%MatrixMarket matrix coordinate integer general
%
60649 6794880 2502872
60631 4704 1
60635 10432 1
16804 16886 1
51809 16933 1
25141 17937 1
10546 19742 1
52780 20161 1
```

#### filtered file: features.tsv, barcodes.tsv, matrix.mtx - 60649 features with 564 cell barcodes
**`head results/star_alignment/Solo.out/Gene/filtered/matrix.mtx`**

``` markdown
%%MatrixMarket matrix coordinate integer general
%
60649 564 1617942
63 1 1
69 1 3
106 1 1
123 1 2
131 1 1
133 1 1
152 1 2
```

#### log files: Log.final.out, Log.out, Log.progress.out
**`cat results/star_alignment/Log.final.out`**

summary mapping statistics after mapping job is complete, very useful for quality control. The statistics are calculated for each read (single- or paired-end) and then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, (unlike the samtools agstat/idxstats, which count each mate separately). Most of the information is collected about the UNIQUE mappers (unlike samtools agstat/idxstats which does not separate unique or multi-mappers). Each splicing is counted in the numbers of splices, which would correspond to summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/da92224f-8e5d-4972-9444-3fd4446069a4)

**`head results/star_alignment/Log.out`**

**`head results/star_alignment/Log.progress.out`**

reports job progress statistics, such as the number of processed reads, % of mapped reads etc. It is updated in 1 minute intervals.

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/6e94efe5-2bb0-433d-9912-7f21ac7b82d5)

## 3.2 Stellarscope Cellsort

The second step to prepare Stellarscope’s input is **sorting the alignment BAM file so that all the cell barcodes bundled together, all the alignments for each read are together**, and **the cell barcodes that are
not in the whitelist are filtered out**. This is achieved using the `stellarscope cellsort [samfile] [whitelist]` command. Additionally, the argument `--outfile` is used to indicate the output sorted BAM file (in this case it’ll be called `results/stellarscope/Aligned.sortedByCB.bam`), and a `--tempdir` (e.g. `/tmp`, most modern Unix-like systems, including Linux distributions and macOS, typically include a /tmp directory as part of
their file system hierarchy). Using 15 CPU cores, `stellarscope cellsort` completed within 13 minutes.

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

>#### error
> ![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/488494e2-9e66-49e6-9519-c4c0e11529eb)
> ![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/0c9da6b1-d31e-464b-a718-40b9d000a1dc)

``` bash
# reinstall the stellarscope package
pip uninstall stellarscope
pip install -e .
```

# 4 Basic Stellarscope analysis

The basic usage of Stellarscope is `stellarscope assign [samfile] [gtffile]`. Where the `samfile` is the alignment containing multimapping reads (e.g. the CB-sorted `results/stellarscope/Aligned.sortedByCB.bam` file), and `gtffile` is the TE annotation in GTF format. The file `resources/retro.hg38.v1.gtf` defines the transcriptional unit of each TE to be quantified, for more details see the Stellarscope paper.

>#### file check: retro.hg38.v1.gtf
> ![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/3da4f926-800b-4cdc-846f-4b59b4a7ce13)

``` bash
# you can use this command to explore the annotation
head resources/retro.hg38.v1.gtf 
```

The `stellarscope assign` command will do the following operations: filter alignments by the user-provided list of passing barcodes (`--whitelist`), create an initial weight matrix with fragments and candidate assigned features, perform our multimapper-aware UMI deduplication, fit a bayesian model to the deduplicated weight matrix using an EM, once the EM is complete, the posterior probability matrix is used to reassign each ambiguous fragment to their generating locus.

Outputs from `stellarscope assign` include a sparse counts matrix with values for the TEs in the annotation, a stats file with a report of the run, and (if requested) an updated BAM. Additionally, at specific checkpoints throughout Stellarscope’s execution, Python pickle files are produced. These files mark the operations described above ensuring that key stages of the process are preserved and enabling easy retrival and continuation in case of interruptions, or for future reference.

In this Stellarscope analysis, the argument `--pooling_mode` is set to `individual`, thus **fitting one model for each cell barcode (i.e resolving ambiguous alignments to TEs within each cell)**.

Using 15 CPU cores, `stellarscope assign` completes within ~2.5 hour.

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
  --max_iter 500 --seed 240626 \
  --updated_sam \
  results/stellarscope/Aligned.sortedByCB.bam \
  resources/retro.hg38.v1.gtf > results/stellarscope/log/stellarscope_assign_individual_240626.log 2>&1 &
```

>#### error
>![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/3ef8fb05-4b03-4a4e-9328-b3bed1eb873e)
>
>set a seed for stellarscope assign

>#### arguments check: stellarscope assign
>- The string provided to `--exp_tag` sets the basename for all the Stellarscope output files.
>- The argument `--stranded_mode` is used to consider feature strand when assigning reads, and here it is set to F as [10x libraries are stranded](https://kb.10xgenomics.com/hc/en-us/articles/360004396971-Are-10x-Single-Cell-gene-expression-libraries-strand-specific) (you can find more information on each protocol’s sheet).
>- In this case the `--reassign_mode` argument is set to “`best_exclude`”. This means that fragments with multiple best assignments (i.e. their values in the posterior probability matrix are equal) are excluded from the final counts.
>- The argument `--max_iter` controls the maximum number of iterations for the EM algorithm.
>- The argument `--updated_sam` is useful to obtain an updated alignment file (`results/stellarscope/pbmc500_individual-updated.bam`)

### output
These are all the outputs in the `results/stellarscope/individual`

``` bash
# inspect the analysis results
tree -t results/stellarscope
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/a7a27142-e6de-45de-9823-c0dfb812862e)

``` markdown
results/stellarscope/individual
├── pbmc500_individual-other.bam
├── pbmc500_individual-tmp_tele.bam
├── pbmc500_individual-checkpoint.load_alignment.pickle
├── pbmc500_individual-umi_tracking.txt
├── pbmc500_individual-checkpoint.dedup_umi.pickle
├── pbmc500_individual-barcodes.tsv
├── pbmc500_individual-features.tsv
├── pbmc500_individual-TE_counts.mtx
├── pbmc500_individual-stats.final.tsv
├── pbmc500_individual-updated.bam
└── pbmc500_individual-checkpoint.final.pickle

0 directories, 11 files
```

### check the output
> #### file check： results files-individual
> - **pbmc500_individual-other.bam​, ​pbmc500_individual-tmp_tele.bam**​: contain the alignments that do and do not overlap the TE annotation, respectively、
> - ​**pbmc500_individual-checkpoint.load_alignment.pickle​**: Stellarscope’s loading of alignments that intersect with the TE annotation
> - ​**pbmc500_individual-umi_tracking.txt, ​pbmc500_individual-checkpoint.dedup_umi.pickle​**​: a record and checkpoint of the PCR duplicates removal process
> - the ​**sparse matrix files​** ( ​*-barcodes.tsv​ ,  ​*-features.tsv​ , and  ​*-TE_counts.mtx​ )
> - the ​**stats report file​** ( ​*-stats.final.tsv​ ), the ​**updated BAM​** ( ​*-updated.bam​ ) and ​**final checkpoint​** ( ​*-checkpoint.final.pickle​ )

#### pickle file: pbmc500_individual-checkpoint.load_alignment.pickle, pbmc500_individual-checkpoint.dedup_umi.pickle, pbmc500_individual-checkpoint.final.pickle

``` markdown
>>> import pickle
>>> with open("/home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-checkpoint.load_alignment.pickle", "rb") as file:
...     data1 = pickle.load(file)
...     print(data1)
... 
<Stellarscope samfile=results/stellarscope/Aligned.sortedByCB.bam, gtffile=resources/retro.hg38.v1.gtf>
>>> with open("/home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-checkpoint.dedup_umi.pickle", "rb") as file:
...     data2 = pickle.load(file)
...     print(data2)
... 
<Stellarscope samfile=results/stellarscope/Aligned.sortedByCB.bam, gtffile=resources/retro.hg38.v1.gtf>
>>> with open("/home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-checkpoint.final.pickle", "rb") as file:
...     data3 = pickle.load(file)
...     print(data3)
... 
<Stellarscope samfile=results/stellarscope/Aligned.sortedByCB.bam, gtffile=resources/retro.hg38.v1.gtf>
```

#### Results of TE annotation for 564 pbmc cells with 27813 features: pbmc500_individual-barcodes.tsv, pbmc500_individual-features.tsv, pbmc500_individual-TE_counts.mtx - the same as results in celltype analysis 

`head /home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-barcodes.tsv -n 20`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/0c4d5ece-990c-463a-8a07-205b3eda7ae6)

`head /home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-features.tsv -n 20`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/54b17008-ed91-4bc6-8faa-18475f19f0b3)

`head -n 20 /home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-TE_counts.mtx`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/7c955335-6f5e-4919-857e-7a517314f57e)

#### other files: pbmc500_individual-umi_tracking.txt, pbmc500_individual-stats.final.tsv

`head -n 20 /home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-umi_tracking.txt`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/b93cde23-79a5-40f9-8c6e-03a02f134bb0)

  `less /home/liumy/software/stellarscope/results/stellarscope/individual/pbmc500_individual-stats.final.tsv`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/6e761892-4fa8-414c-ab73-24689862d506)

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/4227ffa9-1d32-423d-af62-145e08b91cdf)


# 5 Advanced Stellarscope analysis: multiple pooling modes

For certain analyses, the need to conduct Stellarscope’s reassigning using multiple pooling modes may arise. In this scenario it will be beneficial to start by loading the data and generating the checkpoint file corresponding to the UMI deduplication stage. This will provide the flexibility to employ different pooling modes and arguments, facilitating reassignment procedures.

You can achieve this by combining the commands `stellarscope assign` with the flag `--skip_em` and `stellarscope resume`.

## 5.1 Stellarscope load

**In this case (`--skip_em`), there is no pooling mode argument, no reassign mode argument, and no iterations argument.**

Stellarscope will load the alignments intersecting TE annotation and remove PCR duplicates in around 1.5 hours.

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

Upon inspection of the `results/stellarscope` directory, you can see that the two checkpoints `pbmc500_stload-checkpoint.load_alignment.pickle` and `pbmc500_stload-checkpoint.dedup_umi.pickle` have been created.

In this case, the argument `--logfile` was provided, and thus, the file `results/stellarscope/pbmc500_stload.log` is created. This file captures all the information that would have been displayed in the terminal during the Stellarscope’s execution. It includes messages, warnings, and errors (if any) generated as it ran.

### Output

``` bash
# the stellarscope files are output with the prefix `pbmc500_stload`
tree -t results/stellarscope
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/32ad6190-9e98-4faa-a1ef-a500ac073356)

``` markdown
├── pbmc500_stload-other.bam
├── pbmc500_stload-tmp_tele.bam
├── pbmc500_stload-checkpoint.load_alignment.pickle
├── pbmc500_stload-umi_tracking.txt
├── pbmc500_stload-checkpoint.dedup_umi.pickle
├── pbmc500_stload.log
├── pbmc500_stload-stats.final.tsv
```

### check the output
#### bam file: pbmc500_stload-other.bam, pbmc500_stload-tmp_tele.bam, contain the alignments that do and do not overlap the TE annotation, respectively
`samtools view pbmc500_stload-other.bam | head`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/0d701c21-1e98-4964-a89c-ee915dbf1530)

#### pickle file: pbmc500_stload-checkpoint.load_alignment.pickle, pbmc500_stload-checkpoint.dedup_umi.pickle

``` markdown
>>> import pickle
>>> with open("/home/liumy/software/stellarscope/results/stellarscope/pbmc500_stload-checkpoint.load_alignment.pickle", "rb") as file:
...     data = pickle.load(file)
...     print(data)
... 
<Stellarscope samfile=results/stellarscope/Aligned.sortedByCB.bam, gtffile=resources/retro.hg38.v1.gtf>

>>> with open("/home/liumy/software/stellarscope/results/stellarscope/pbmc500_stload-checkpoint.dedup_umi.pickle", "rb") as file:
...     data = pickle.load(file)
...     print(data)
... 
<Stellarscope samfile=results/stellarscope/Aligned.sortedByCB.bam, gtffile=resources/retro.hg38.v1.gtf>
```

#### other files: pbmc500_stload-umi_tracking.txt, pbmc500_stload.log, pbmc500_stload-stats.final.tsv

`head -n 20 /home/liumy/software/stellarscope/results/stellarscope/pbmc500_stload-umi_tracking.txt`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/cb8b5f43-4843-4400-b10e-ac2688d0df1d)

`less /home/liumy/software/stellarscope/results/stellarscope/pbmc500_stload.log`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/66438766-7473-48ff-b165-36b4f2e95d33)

`less /home/liumy/software/stellarscope/results/stellarscope/pbmc500_stload-stats.final.tsv`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/663c4010-3381-4c82-95a8-fa9c0edd9f03)


## 5.2 Stellarscope resume

Now, two Stellarscope analyses will demonstrate how to use the `stellarscope_resume` command.
- In the first one, `--pooling_mode` pseudobulk is used, and matrices employing every reassignment mode are generated.
- In the second one, cell type labels from annotating the PBMCs using the reference transcriptomes from [azimuth](https://azimuth.hubmapconsortium.org/) are used to perform Stellarscope’s reassignment with `--pooling_mode` celltype, and the final counts matrix is generated by setting a threshold and considering assignments that exceed it.

### 5.2.1 Stellarscope Pseudobulk

In the pooling mode `pseudobulk` a single model is fitted, meaning that **ambiguous alignments are resolved sharing information across all cells**.

For this analysis, the `stellarscope resume` command accepts reassignment arguments (pooling mode, iterations number, etc.) and a `*-checkpoint.dedup_umi.pickle` file (instead of a BAM file).

The reassignment mode is the method that Stellarscope will use to resolve which one is the “best” reassignment for fragments that have multiple possible reassignments. In this case, the argument `--use_every_reasign_mode` is added so that Stellarscope outputs seven `.mtx` files in total.
- `best_exclude` (fragments with multiple best assignments are excluded)
- `best_conf` (assignments with values above a threshold set by `--conf_prob` are counted)
- `best_random` (the best assignment is randomly chosen from the best assignments set)
- `best_average` (the fragment is divided evenly among the best assignments)
- `initial_unique` (only the originally uniquely aligned reads are included; EM is not considered)
- `initial_random` (an assignment is randomly selected from the initial alignments; EM is not considered）
- `total_hits` (all assignments are included; EM is not considered)
**Note**: The final three modes do not perform reassignment or model fitting but are included for comparison. Note that the matrix corresponding to `best_exclude` will be called `pbmc500_pseudobulk-TE_counts.mtx`, as opposed to the other matrices, which will have the reassignment mode indicated in their name (e.g. `pbmc500_pseudobulk-TE_counts.best_conf.mtx`). This is because `best_exclude` is the default option.

This Stellarscope resume analysis took less than 2 minutes in total.

``` bash
# create a directory for stellarscope pseudobulk results
mkdir -p results/stellarscope/pseudobulk

# use stellarscope resume to continue from the deduplication checkpoint
nohup stellarscope resume \
            --exp_tag pbmc500_pseudobulk \
            --outdir results/stellarscope/pseudobulk \
            --nproc 15 \
            --pooling_mode pseudobulk \
            --use_every_reassign_mode \
            --max_iter 500\
            --updated_sam --seed 240626 \
            results/stellarscope/pbmc500_stload-checkpoint.dedup_umi.pickle \
            --logfile results/stellarscope/pseudobulk/pbmc500_pseudobulk.log > results/stellarscope/log/stellarscope_resume_pseudobulk_240626.log 2>&1 &
```

### Output
``` bash
# now the resutls from pooling mode pseudobulk are in the `pseudobulk` directory
tree -t results/stellarscope/pseudobulk
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/00dbb5e4-34c7-4cec-9433-57c81c74c52a)

``` markdown
results/stellarscope/pseudobulk
├── pbmc500_pseudobulk-barcodes.tsv
├── pbmc500_pseudobulk-features.tsv
├── pbmc500_pseudobulk-TE_counts.mtx
├── pbmc500_pseudobulk-TE_counts.best_conf.mtx
├── pbmc500_pseudobulk-TE_counts.best_random.mtx
├── pbmc500_pseudobulk-TE_counts.best_average.mtx
├── pbmc500_pseudobulk-TE_counts.initial_unique.mtx
├── pbmc500_pseudobulk-TE_counts.initial_random.mtx
├── pbmc500_pseudobulk-stats.final.tsv
├── pbmc500_pseudobulk-TE_counts.total_hits.mtx
├── pbmc500_pseudobulk-checkpoint.final.pickle
└── pbmc500_pseudobulk.log

0 directories, 12 files
```

### check the output

#### Results of TE annotation for 564 pbmc cells with 27813 features: pbmc500_pseudobulk-barcodes.tsv, pbmc500_pseudobulk-features.tsv, pbmc500_pseudobulk-TE_counts.mtx - the same as results in celltype analysis 

`head /home/liumy/software/stellarscope/results/stellarscope/pseudobulk/pbmc500_pseudobulk-TE_counts.mtx -n 20`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/9c94a4c6-25cb-4533-9881-0e909267da5d)

##### other files: pbmc500_pseudobulk-stats.final.tsv, pbmc500_pseudobulk.log, pbmc500_pseudobulk-checkpoint.final.pickle, pbmc500_pseudobulk-TE_counts.total_hits.mtx

`less /home/liumy/software/stellarscope/results/stellarscope/pseudobulk/pbmc500_pseudobulk-stats.final.tsv`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/2d5f0deb-b5b8-4b1b-8c66-a51c8e833081)

`less /home/liumy/software/stellarscope/results/stellarscope/pseudobulk/pbmc500_pseudobulk.log`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/e9af1a8e-f13a-4bf9-a115-22a10e814061)

``` markdown
>>> with open("/home/liumy/software/stellarscope/results/stellarscope/pseudobulk/pbmc500_pseudobulk-checkpoint.final.pickle", "rb") as file:
...     data = pickle.load(file)
...     print(data)
... 
<Stellarscope checkpoint=results/stellarscope/pbmc500_stload-checkpoint.dedup_umi.pickle>
```

`head -n 20 /home/liumy/software/stellarscope/results/stellarscope/pseudobulk/pbmc500_pseudobulk-TE_counts.total_hits.mtx`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/5718217f-0b13-4159-a5dc-e5d6ce616a08)


### 5.2.2 Stellarscope Celltype

This is an example of fitting one model per cell type. This strategy depends on user-provided cell identities. Using this `celltype` pooling mode carries the assumption that the (retro)transcriptome is shared more closely between cells that are classified as the same cell type. Annotating/solving cell identity is an entire research line, so users should be aware of their goals and the strategies they use to determine cell types in their data.

To exemplify this Stellarscope analysis, we have provided within the resources a file `resources/celltypes_tsv/pbmc500_azimuth_l1.tsv` where each of the PBMC’s cell barcodes has an assigned cell types obtained from azimuth. This is provided to Stellarscope using the argument `--celltype_tsv`.

This Stellarscope analysis took around 10 minutes in total.

``` bash
# create a directory for stellarscope celltype results
mkdir -p results/stellarscope/celltype

# use stellarscope resume to continue from the deduplication checkpoint
nohup stellarscope resume \
            --exp_tag pbmc500_l1 \
            --outdir results/stellarscope/celltype \
            --nproc 15 \
            --pooling_mode celltype \
            --celltype_tsv resources/celltypes_tsv/pbmc500_azimuth_l1.tsv \
            --reassign_mode best_conf \
            --conf_prob 0.95 \
            --max_iter 500\
            --updated_sam --seed 240626 \
            results/stellarscope/pbmc500_stload-checkpoint.dedup_umi.pickle \
            --logfile results/stellarscope/celltype/pbmc500_l1.log > results/stellarscope/log/stellarscope_resume_celltype_240626.log 2>&1 &
```

### Output
``` bash
# stellarscope pooling mode celltype results
tree -t results/stellarscope/celltype
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/fc450a29-9a71-4734-b4d0-a1655e9872a7)

``` markdown
results/stellarscope/celltype
├── pbmc500_l1-barcodes.tsv
├── pbmc500_l1-features.tsv
├── pbmc500_l1-stats.final.tsv
├── pbmc500_l1-TE_counts.mtx
├── pbmc500_l1-checkpoint.final.pickle
└── pbmc500_l1.log

0 directories, 6 files
```

### check the output

##### Running process

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/c27cc36b-9493-4d3d-82e2-4823eb52c737)

##### cell type annotation file from [Azimuth](https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC)

`head resources/celltypes_tsv/pbmc500_azimuth_l1.tsv`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/99fa77d7-57c4-46d6-8826-d6391120dab3)

`head resources/celltypes_tsv/pbmc500_azimuth_l2.tsv`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/fc74f0bb-aaaf-4077-8353-bcc94089064f)

##### Results of TE annotation for 564 pbmc cells with 27813 features: pbmc500_l1-barcodes.tsv, pbmc500_l1-features.tsv, pbmc500_l1-TE_counts.mtx

`head /home/liumy/software/stellarscope/results/stellarscope/celltype/pbmc500_l1-barcodes.tsv -n 20`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/10529326-37d8-460d-88ec-9dd8ccd22c5f)

`head /home/liumy/software/stellarscope/results/stellarscope/celltype/pbmc500_l1-features.tsv -n 20`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/89ad2733-766f-4d60-8448-6f7bc05968c9)

`head /home/liumy/software/stellarscope/results/stellarscope/celltype/pbmc500_l1-TE_counts.mtx -n 20`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/d5f97eda-51f5-4ceb-b64d-18ebd1fea4a2)

##### other files: pbmc500_l1-stats.final.tsv, pbmc500_l1-checkpoint.final.pickle, pbmc500_l1.log

`less /home/liumy/software/stellarscope/results/stellarscope/celltype/pbmc500_l1-stats.final.tsv`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/0f4602e2-074a-4fc1-aaa0-d3e7c37b1307)

``` python
with open('/home/liumy/software/stellarscope/results/stellarscope/celltype/pbmc500_l1-checkpoint.final.pickle', 'rb') as file:
     data = pickle.load(file)
     print(data)
```

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/e2336acf-f063-4a88-8e58-0230bda6c840)


`less /home/liumy/software/stellarscope/results/stellarscope/celltype/pbmc500_l1.log`

![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/b7135686-2876-4c03-906d-8219d0f11705)
