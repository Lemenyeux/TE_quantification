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

# run STAR alignment - named results/star_alignment_star
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
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR \
  --limitOutSJcollapsed 5000000 \
  --outFilterMultimapNmax 500 \
  --outFilterMultimapScoreRange 5 \
  --outFileNamePrefix results/star_alignment/ > ${log_file}/STAR_240628.log 2>&1 &

# [31019 16:37-17:25 Estimated Number of Cells,6447 without --soloCellFilter EmptyDrops_CR]
# [Estimated Number of Cells,6447 without --soloCellFilter EmptyDrops_CR]

# genomeGenerate: generate STAR genome index identical to CellRanger's
STAR  --runMode genomeGenerate \
    --runThreadN 15 \
    --genomeDir /home/liumy/software/stellarscope/resources/CR_refdata-gex-GRCh38-2024-A/ \
    --genomeFastaFiles /home/liumy/software/cellranger/refdata-gex-GRCh38-2024-A/fasta/genome.fa  \
    --sjdbGTFfile /home/liumy/software/stellarscope/resources/genes.gtf

# rerun STAR alignment using STAR genome index identical to CellRanger's
nohup STAR \
  --runThreadN 15 \
  --genomeDir ${resource_file}/CR_refdata-gex-GRCh38-2024-A \
  --readFilesIn ${fastqfile}/D19-4295_S1_L001_R2_001.fastq.gz,${fastqfile}/D19-4295_S1_L002_R2_001.fastq.gz,${fastqfile}/D19-4295_S1_L003_R2_001.fastq.gz,${fastqfile}/D19-4295_S1_L004_R2_001.fastq.gz ${fastqfile}/D19-4295_S1_L001_R1_001.fastq.gz,${fastqfile}/D19-4295_S1_L002_R1_001.fastq.gz,${fastqfile}/D19-4295_S1_L003_R1_001.fastq.gz,${fastqfile}/D19-4295_S1_L004_R1_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --soloCBwhitelist ${resource_file}/whitelist_10x/3M-february-2018.txt \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS NM nM MD CR CY UR UY CB UB GX GN sS sQ sM \
  --outSAMtype BAM SortedByCoordinate \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR \
  --limitOutSJcollapsed 5000000 \
  --outFilterMultimapNmax 500 \
  --outFilterMultimapScoreRange 5 \
  --outFileNamePrefix results/star_alignment/ > ${log_file}/STAR_cellranger240705.log 2>&1 &

# [111427 Estimated Number of Cells,8515]


# check the output
tree -t results/star_alignment
```

## check the codes
> `--soloUMIfiltering MultiGeneUMI_CR`: basic + remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0. Only works with --soloUMIdedup 1MM_CR
> `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`: same as 1MM_multi_pseudocounts, multimatching to WL is allowed for CBs with N-bases. This option matches best with CellRanger >= 3.0.0
> `--soloUMIdedup 1MM_CR`: CellRanger2-4 algorithm for 1MM UMI collapsing.
> `--soloCellFilter EmptyDrops_CR`: EmptyDrops filtering in CellRanger flavor. Please cite the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y.
>> Can be followed by 10 numeric parameters:
>> nExpectedCells   maxPercentile   maxMinRatio   indMin   indMax   umiMin   umiMinFracMedian   candMaxN   FDR   simN
>> The harcoded values are from CellRanger:
>> 3000             0.99            10            45000    90000    500      0.01                20000     0.01  10000
> `clipAdapterType CellRanger4`: 5p and 3p adapter clipping similar to CellRanger4. Utilizes Opal package by Martin Šošić: https://github.com/Martinsos/opal

## How to make STARsolo raw gene counts (almost) identical to CellRanger's
CellRanger uses its own "filtered" version of annotations (GTF file) which is a subset of ENSEMBL annotations, with several gene biotypes removed (mostly small non-coding RNA). Annotations affect the counts, and to match CellRanger counts CellRanger annotations have to be used. 10X provides several versions of the CellRanger annotations: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest. For the best match, the annotations in CellRanger run and STARsolo run should be exactly the same.

## check the results------------------
> **STAR using the index identical to Cellranger's**
> **/home/liumy/TEQTL/scRNA/stellarscope/results/star_alignment/Log.final.out  STAR using the index identical to Cellranger's**
``` markdown


                                 Started job on |       Jul 05 10:51:47
                             Started mapping on |       Jul 05 10:56:47
                                    Finished on |       Jul 05 12:02:29
       Mapping speed, Million of reads per hour |       176.11

                          Number of input reads |       192840779
                      Average input read length |       48
                                    UNIQUE READS:
                   Uniquely mapped reads number |       122935083
                        Uniquely mapped reads % |       63.75%
                          Average mapped length |       54.59
                       Number of splices: Total |       1221691
            Number of splices: Annotated (sjdb) |       1072500
                       Number of splices: GT/AG |       1117112
                       Number of splices: GC/AG |       20942
                       Number of splices: AT/AC |       1675
               Number of splices: Non-canonical |       81962
                      Mismatch rate per base, % |       0.38%
                         Deletion rate per base |       0.01%
                        Deletion average length |       1.30
                        Insertion rate per base |       0.01%
                       Insertion average length |       1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       22671591
             % of reads mapped to multiple loci |       11.76%
        Number of reads mapped to too many loci |       8
             % of reads mapped to too many loci |       0.00%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       40318965
                 % of reads unmapped: too short |       20.91%
                Number of reads unmapped: other |       6915132
                     % of reads unmapped: other |       3.59%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```
> **/home/liumy/TEQTL/scRNA/stellarscope/results/star_alignment/Solo.out/Gene/Summary.csv STAR using the index identical to Cellranger's**
``` marndown
Number of Reads,192840779
Reads With Valid Barcodes,0.981961
Sequencing Saturation,0.144525
Q30 Bases in CB+UMI,0.963877
Q30 Bases in RNA read,0.911485
Reads Mapped to Genome: Unique+Multiple,0.755062
Reads Mapped to Genome: Unique,0.637495
Reads Mapped to Gene: Unique+Multiple Gene,NoMulti
Reads Mapped to Gene: Unique Gene,0.193852
Estimated Number of Cells,8515
Unique Reads in Cells Mapped to Gene,26308021
Fraction of Unique Reads in Cells,0.703752
Mean Reads per Cell,3089
Median Reads per Cell,2379
UMIs in Cells,22289943
Mean UMI per Cell,2617
Median UMI per Cell,2025
Mean Gene per Cell,1568
Median Gene per Cell,1337
Total Gene Detected,29575
```
> **/home/liumy/TEQTL/scRNA/stellarscope/results/star_alignment/Solo.out/Gene/Summary.csv STAR using STAR genome index withou -soloCellFilter EmptyDrops_CR**
``` markdown
Number of Reads,192840779
Reads With Valid Barcodes,0.981986
Sequencing Saturation,0.141748
Q30 Bases in CB+UMI,0.963877
Q30 Bases in RNA read,0.911485
Reads Mapped to Genome: Unique+Multiple,0.755183
Reads Mapped to Genome: Unique,0.637804
Reads Mapped to Gene: Unique+Multiple Gene,NoMulti
Reads Mapped to Gene: Unique Gene,0.179157
Estimated Number of Cells,6447
Unique Reads in Cells Mapped to Gene,22846101
Fraction of Unique Reads in Cells,0.661271
Mean Reads per Cell,3543
Median Reads per Cell,2896
UMIs in Cells,19389259
Mean UMI per Cell,3007
Median UMI per Cell,2472
Mean Gene per Cell,1820
Median Gene per Cell,1624
Total Gene Detected,35381
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

# 6 other tests
``` bash
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
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR \
  --limitOutSJcollapsed 5000000 \
  --outFilterMultimapNmax 500 \
  --outFilterMultimapScoreRange 5 \
  --outFileNamePrefix results/star_alignment/cellfilter > ${log_file}/STAR_testcellfilter_240630.log 2>&1 &

# [40163 Estimated Number of Cells,8396 ]


# soloCellFiltering
STAR --runMode soloCellFiltering \
    results/star_alignment/Solo.out/Gene/raw/ \
    results/star_alignment/Solo.out/Gene/filtered_emptyDrop/ \
    --soloCellFilter EmptyDrops_CR \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts 

# [Estimated Number of Cells: 8396]
# [Estimated Number of Cells: 8396 with --soloStrand Reverse]
```
