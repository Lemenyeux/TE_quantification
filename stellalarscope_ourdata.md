# 1. Setup

``` bash
#1. make and change working directory
mkdir /home/liumy/TEQTL/scRNA/stellarscope
cd /home/liumy/TEQTL/scRNA/stellarscope

#2. check the fastq datasets
tree -t /home/liumy/TEQTL/scRNA/raw_data/testing_sample
fastqfile="/home/liumy/TEQTL/scRNA/raw_data/testing_sample"

#3. check the resources file
tree -t /home/liumy/software/stellarscope/resources
resource_file="/home/liumy/software/stellarscope/resources"

log_file="/home/liumy/TEQTL/scRNA/stellarscope/log"
```

``` markdown
/home/liumy/TEQTL/scRNA/raw_data/testing_sample
├── D19-4295_S1_L001_R1_001.fastq.gz
├── D19-4295_S1_L001_R2_001.fastq.gz
├── D19-4295_S1_L002_R1_001.fastq.gz
├── D19-4295_S1_L003_R1_001.fastq.gz
├── D19-4295_S1_L002_R2_001.fastq.gz
├── D19-4295_S1_L004_R1_001.fastq.gz
├── D19-4295_S1_L003_R2_001.fastq.gz
└── D19-4295_S1_L004_R2_001.fastq.gz

0 directories, 8 files
```

# 2. STAR alignment
``` bash
# activate the envrionment
conda activate stellarscope

# create directory for STAR results
mkdir -p results/star_alignment

# run STAR alignment
nohup STAR \
  --runThreadN 15 \
  --genomeDir ${resource_file}/STAR_GRCh38.d1.vd1_gencode.v38 \
  --readFilesIn ${fastqfile}/D19-4295_S1_L001_R2_001.fastq.gz,${fastqfile}/D19-4295_S1_L002_R2_001.fastq.gz,${fastqfile}/D19-4295_S1_L003_R2_001.fastq.gz,${fastqfile}/D19-4295_S1_L004_R2_001.fastq.gz ${fastqfile}/D19-4295_S1_L001_R1_001.fastq.gz,${fastqfile}/D19-4295_S1_L002_R1_001.fastq.gz,${fastqfile}/D19-4295_S1_L003_R1_001.fastq.gz,${fastqfile}/D19-4295_S1_L004_R1_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --soloCBwhitelist ${resource_file}/whitelist_10x/3M-february-2018.txt \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS NM nM MD CR CY UR UY CB UB GX GN sS sQ sM \
  --outSAMtype BAM SortedByCoordinate \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --limitOutSJcollapsed 5000000 \
  --outFilterMultimapNmax 500 \
  --outFilterMultimapScoreRange 5 \
  --outFileNamePrefix results/star_alignment/ > ${log_file}/STAR_240628.log 2>&1 &

# [31019 16:37-17:25]

# check the output
tree -t results/star_alignment
```

# 3. Stellarscope Cellsort
``` bash
# create directory for Stellarscope results
mkdir -p results/stellarscope

## sort reads alignments by CB
nohup stellarscope cellsort \
  --nproc 15 \
  --tempdir /tmp \
  --outfile results/stellarscope/Aligned.sortedByCB.bam \
  results/star_alignment/Aligned.sortedByCoord.out.bam \
  results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv > ${log_file}/stellarscope_cellsort_240628.log 2>&1 &

# [31227 17:27 - 17:56]
```

# 4. Basic Stellarscope analysis
``` bash
## create a directory for stellarscope individual results
mkdir -p results/stellarscope/individual

## stellarscope analysis (pooling mode individual)
nohup stellarscope assign \
  --exp_tag D19_4295_individual \
  --outdir results/stellarscope/individual \
  --nproc 15 \
  --stranded_mode F \
  --whitelist results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
  --pooling_mode individual \
  --reassign_mode best_exclude \
  --max_iter 500 --debug --seed 240628 \
  --updated_sam \
  results/stellarscope/Aligned.sortedByCB.bam \
  ${resource_file}/retro.hg38.v1.gtf > ${log_file}/stellarscope_assign_individual_240628.log 2>&1 &

# [31475 18:19 2h30min]

## inspect the analysis results
tree -t results/stellarscope/individual
```

# 5. Advanced Stellarscope analysis: multiple pooling modes
## 5.1 Stellarscope load
``` bash
# load the alignment data into stellarscope
nohup stellarscope assign \
  --exp_tag D19_4295_stload \
  --outdir results/stellarscope \
  --skip_em \
  --stranded_mode F \
  --whitelist results/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
  --updated_sam --seed 240628 \
  results/stellarscope/Aligned.sortedByCB.bam \
  resources/retro.hg38.v1.gtf \
  --logfile results/stellarscope/D19_4295_stload.log > ${log_file}/stellarscope_assign_stload_240628.log 2>&1 &
[18600]

# the stellarscope files are output with the prefix `pbmc500_stload`
tree -t results/stellarscope
```

## 5.2 Stellarscope resume - Pseudobulk
``` bash
# create a directory for stellarscope pseudobulk results
mkdir -p results/stellarscope/pseudobulk

# use stellarscope resume to continue from the deduplication checkpoint
nohup stellarscope resume \
            --exp_tag D19_4295_pseudobulk \
            --outdir results/stellarscope/pseudobulk \
            --nproc 15 \
            --pooling_mode pseudobulk \
            --reassign_mode best_conf \
            --max_iter 500\
            --updated_sam --seed 240626 \
            results/stellarscope/D19_4295_stload-checkpoint.dedup_umi.pickle \
            --logfile results/stellarscope/pseudobulk/pbmc500_pseudobulk.log > ${log_file}/stellarscope_resume_pseudobulk_240626.log 2>&1 &

#[20707 less than 2 mins]

# stellarscope pooling mode pseudobulk results
tree -t results/stellarscope
```

## 5.3 Stellarscope resume - Celltype - no celltype annotation file - no run
``` bash
# create a directory for stellarscope celltype results
mkdir -p results/stellarscope/celltype

# use stellarscope resume to continue from the deduplication checkpoint
nohup stellarscope resume \
            --exp_tag D19_4295_celltype \
            --outdir results/stellarscope/celltype \
            --nproc 15 \
            --pooling_mode celltype \
            --celltype_tsv resources/celltypes_tsv/pbmc500_azimuth_l1.tsv \
            --reassign_mode best_conf \
            --conf_prob 0.95 \
            --max_iter 500\
            --updated_sam --seed 240626 \
            results/stellarscope/D19_4295_stload-checkpoint.dedup_umi.pickle \
            --logfile results/stellarscope/celltype/D19_4295_celltype.log > results/stellarscope/log/stellarscope_resume_celltype_240626.log 2>&1 &

#[21021]

# stellarscope pooling mode celltype results
tree -t results/stellarscope/celltype
```
