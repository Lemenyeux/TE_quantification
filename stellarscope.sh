#1. make and change working directory
mkdir stellarscope
cd /home/liumy/software/stellarscope

#2. download datasets
## download scRNA-seq data
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_fastqs.tar
tar -xf 500_PBMC_3p_LT_Chromium_X_fastqs.tar
### sudo apt-get install tree
tree 500_PBMC_3p_LT_Chromium_X_fastqs

#500_PBMC_3p_LT_Chromium_X_fastqs
#├── 500_PBMC_3p_LT_Chromium_X_S4_L003_I1_001.fastq.gz
#├── 500_PBMC_3p_LT_Chromium_X_S4_L003_I2_001.fastq.gz
#├── 500_PBMC_3p_LT_Chromium_X_S4_L003_R1_001.fastq.gz
#├── 500_PBMC_3p_LT_Chromium_X_S4_L003_R2_001.fastq.gz
#├── 500_PBMC_3p_LT_Chromium_X_S4_L004_I1_001.fastq.gz
#├── 500_PBMC_3p_LT_Chromium_X_S4_L004_I2_001.fastq.gz
#├── 500_PBMC_3p_LT_Chromium_X_S4_L004_R1_001.fastq.gz
#└── 500_PBMC_3p_LT_Chromium_X_S4_L004_R2_001.fastq.gz
#0 directories, 8 files

## Obtain resources (30mins)
## curl -OJL https://zenodo.org/records/10671822/files/resources.tar.gz?download=1
## curl: (56) Recv failure: Connection timed out
wget -c -O resources.tar.gz https://zenodo.org/records/10671822/files/resources.tar.gz?download=1

tar -xvf resources.tar.gz
tree resources

#resources
#├── celltypes_tsv
#│   ├── pbmc500_azimuth_l1.tsv
#│   └── pbmc500_azimuth_l2.tsv
#├── retro.hg38.v1.gtf
#├── STAR_GRCh38.d1.vd1_gencode.v38
#│   ├── chrLength.txt
#│   ├── chrNameLength.txt
#│   ├── chrName.txt
#│   ├── chrStart.txt
#│   ├── exonGeTrInfo.tab
#│   ├── exonInfo.tab
#│   ├── geneInfo.tab
#│   ├── Genome
#│   ├── genomeParameters.txt
#│   ├── SA
#│   ├── SAindex
#│   ├── sjdbInfo.txt
#│   ├── sjdbList.fromGTF.out.tab
#│   ├── sjdbList.out.tab
#│   └── transcriptInfo.tab
#├── stellarscope_protocol_env.yaml
#└── whitelist_10x
#    └── 3M-february-2018.txt
#3 directories, 20 files


#3. Create a Stellarscope conda environment
#conda env create --yes -f resources/stellarscope_protocol_env.yaml
#conda activate stellarscope
##name: stellarscope_protocol_env
##channels:
##  - conda-forge
##  - bioconda
##  - defaults
##dependencies:
##  - python >=3.6
##  - future
##  - pip
##  - pyyaml
##  - cython
##  - numpy
##  - scipy
##  - pysam >=0.19
##  - htslib
##  - intervaltree
##  - pandas
##  - samtools >=1.16
##  - packaging
##  - anndata
##  - star ==2.7.10b
##  - pip:
##      - -e git+ssh://git@github.com/nixonlab/stellarscope.git@dev#egg=stellarscope

conda create -n stellarscope
conda activate stellarscope 

conda install python samtools 'star==2.7.10b'
pip install future pyyaml cython numpy scipy 'pysam>=0.19' htslib intervaltree pandas packaging anndata 
git clone https://github.com/nixonlab/stellarscope.git
cd stellarscope
pip install -e .


#4. Prepare Stellarscope’s input-BAM file
## STAR alignment: align the raw reads from PBMCs to the reference genome
## the alignment includes multimapping reads
## output a BAM file

### create directory for STAR results
mkdir -p results/star_alignment

### run STAR alignment
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

## --genoneDir: The STAR reference genome index
## --readFilesIn: STAR requires that the first FASTQ file(s) contain the cDNA read and the second FASTQ file(s) contain the cell barcode+UMI read, and that they’re separated by a blank space
## --soloCBwhitelist: A whitelist provides the list of known cell barcodes.
## --outSAMunmapped: Include unmapped reads in the output BAM file.
## --outSAMattributes: Specify which SAM attributes to include.
## --outSAMtype BAM SortedByCoordinate: Output sorted BAM file.
## In this alignment multimapping reads are allowed by setting the argument  ​--outFilterMultimapNmax​  to  ​500​ 
## --outFilterMultimapScoreRange​  is set to 5​ so that for each multimapping read, its alignments additional to the best one will actually be included in the output BAM file.

### output
tree -t results/star_alignment

#5 Stellarscope Cellsort
## sorting the alignment BAM file so that all the cell barcodes bundled together, all the
## alignments for each read are together, and the cell barcodes that are not in the whitelist are filtered out.

## create directory for Stellarscope results
mkdir -p results/stellarscope

## sort reads alignments by CB
## ​stellarscope cellsort [samfile] [whitelist]​ command
nohup stellarscope cellsort \
  --nproc 15 \
  --tempdir /tmp \
  --outfile results/stellarscope/Aligned.sortedByCB.bam \
  results/star_alignment/Aligned.sortedByCoord.out.bam \
  results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv > results/stellarscope/log/stellarscope_cellsort_240624.log 2>&1 &
#[804 OK]

#6. Basic Stellarscope analysis
## ​stellarscope assign [samfile] [gtffile]​ 
## ​samfile​: the alignment containing multimapping reads (e.g. the CB-sorted results/stellarscope/Aligned.sortedByCB.bam​ file)
## ​gtffile：the TE annotation in GTF format.
###1. filter alignments by the user-provided list of passing barcodes ( ​--whitelist​ )
###2. create an initial weight matrix with fragments and candidate assigned features
###3. perform our multimapper-aware UMI deduplication
###4. fit a bayesian model to the deduplicated weight matrix using an EM
###5. once the EM is complete, the posterior probability matrix is used to reassign each ambiguous fragment to their generating locus. 

## explore the annotation file
head resources/retro.hg38.v1.gtf

## create a directory for stellarscope individual results
mkdir -p results/stellarscope/individual
mkdir -p results/stellarscope/log

## stellarscope analysis (pooling mode individual)

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

### [989 240624 17:13 error] 
### [9808 240624 17:50]
### ​--exp_tag: the basename for all the Stellarscope output files
### ​--stranded_mode​: consider feature strand when assigning reads, and here it is set to F as 10x libraries are stranded.
### --pooling_mode​ ​individual:​ fitting one model for each cell barcode (i.e resolving ambiguous alignments to TEs within each cell). 
### ​--reassign_mode ​best_exclude: fragments with multiple best assignments (i.e. their values in the posterior probability matrix are equal) are excluded from the final counts.
### ​--max_iter: the maximum number of iterations for the EM algorithm.
### ​--updated_sam​: obtain an updated alignment file

## Output
###1. a sparse counts matrix with values for the TEs in the annotation
###2. a stats file with a report of the run 
###3. an updated BAM file (if requested)
###4. Python pickle files are produced at specific checkpoints throughout Stellarscope’s execution
## inspect the analysis results
tree -t results/stellarscope

### pbmc500_individual-other.bam​, ​pbmc500_individual-tmp_tele.bam​: contain the alignments that do and do not overlap the TE annotation, respectively
### ​pbmc500_individual-checkpoint.load_alignment.pickle: Stellarscope’s loading of alignments that intersect with the TE annotation
### ​pbmc500_individual-umi_tracking.txt, ​pbmc500_individual-checkpoint.dedup_umi.pickle​: a record and checkpoint of the PCR duplicates removal process
### the sparse matrix files ( ​*-barcodes.tsv​ ,  ​*-features.tsv​ , and  ​*-TE_counts.mtx​ )
### the stats report file ( ​*-stats.final.tsv​ ), the updated BAM ( ​*-updated.bam​ ) and final checkpoint ( ​*-checkpoint.final.pickle​ )

#7. Advanced Stellarscope analysis: multiple pooling modes
##7.1 Stellarscope load
### In this case ( ​--skip_em​ ), there is no pooling mode argument, no reassign mode argument, and no iterations argument.
### Stellarscope will load the alignments intersecting TE annotation and remove PCR duplicates in around 45 minutes.

## load the alignment data into stellarscope
nohup stellarscope assign \
  --exp_tag pbmc500_stload \
  --outdir results/stellarscope \
  --skip_em \
  --stranded_mode F \
  --whitelist results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
  --updated_sam \
  results/stellarscope/Aligned.sortedByCB.bam \
  resources/retro.hg38.v1.gtf \
  --logfile results/stellarscope/pbmc500_stload.log > results/stellarscope/log/stellarscope_assign_stload_240624.log 2>&1 &

##chmod +x /home/liumy/software/stellarscope/stellarscope/stellarscope/stellarscope_assign.py
##nohup /home/liumy/software/stellarscope/stellarscope/stellarscope/stellarscope_assign.py \
##  --exp_tag pbmc500_stload \
##  --outdir results/stellarscope \
##  --skip_em \
##  --stranded_mode F \
##  --whitelist results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
##  --updated_sam \
##  results/stellarscope/Aligned.sortedByCB.bam \
##  resources/retro.hg38.v1.gtf \
##  --logfile results/stellarscope/pbmc500_stload.log > results/stellarscope/log/stellarscope_assign_stload_240624.log 2>&1 &

###[1188]
### output
### two checkpoints: ​pbmc500_stload-checkpoint.load_alignment.pickle​; ​pbmc500_stload-checkpoint.dedup_umi.pickle
###​​ --logfile​: ​results/stellarscope/pbmc500_stload.log; capture all the information that would have been displayed in the terminal during the Stellarscope’s execution
### It includes messages, warnings, and errors (if any) generated as it ran.

## the stellarscope files are output with the prefix `pbmc500_stload`
tree -t results/stellarscope

##7.2 Stellarscope resume
### ​--pooling_mode​  pseudobulk: matrices employing every reassignment mode are generated
### --pooling_mode​  celltype: cell type labels from annotating the PBMCs using the reference transcriptomes from azimuth are used
### the final counts matrix is generated by setting a threshold and considering assignments that exceed it.

##7.2.1 Stellarscope Pseudobulk
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



