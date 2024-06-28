# fastq to bam file using cell ranger
``` bash
nohup cellranger count --id=d19_4295 --fastqs=/home/liumy/TEQTL/scRNA/raw_data/testing_sample --sample=D19-4295 --output-dir=/home/liumy/TEQTL/scRNA/results/d19_4295 --transcriptome=/home/liumy/software/cellranger/refdata-gex-GRCh38-2024-A --create-bam=true > /home/liumy/TEQTL/scRNA/results/d19_4295/fastqtobam_d19_4295.log 2>&1 &
# 536650 [240609 17:00-19:22]
```

# exclude lines with blank CB
``` bash
samtools view possorted_genome_bam.bam -h | awk '/^@/ || /CB:/' | samtools view -h -b > possorted_genome_bam.clean.bam
```

``` bash
# creat conda env
conda create -n scte
conda activate scte

# installation
git clone https://github.com/JiekaiLab/scTE.git
cd scTE
conda install python
python setup.py install

# Building genome indices
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
# Main file, gene annotation on reference chromosomes in GTF and GFF3 file formats. These are the main GENCODE gene annotation files. They contain annotation (genes, transcripts, exons, start_codon, stop_codon, UTRs, CDS) on the reference chromosomes,   which are chr1-22, X, Y, M in human and chr1-19, X, Y, M in human.

# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
# description of rmsk.txt.gz file: http://genome.ucsc.edu/cgi-bin/hgTables


scTE_build -g hg38 # Human
scTE_build -g hg38 -m inclusive


vim ~/.bashrc
export PATH=/home/liumy/software/scTE/bin:$PATH
source ~/.bashrc

# Analysis of 10x style scRNA-seq data
cd /home/liumy/TEQTL/scRNA/results/d19_4295_te
nohup scTE -i possorted_genome_bam.clean.bam -o d19_4295 -x /home/liumy/software/scTE/hg38.inclusive.idx -p 10 -CB CB -UMI UB > /home/liumy/TEQTL/scRNA/log/scTE_d19_4295.log 2>&1 &

[544955 240610-11:27-13:16]
[10000 obs. Of 58390 variables]

nohup scTE -i possorted_genome_bam.clean.bam -o d19_4295 -x /home/liumy/software/scTE/hg38.inclusive.idx -p 10 --hdf5 True -CB CB -UMI UB --expect-cell  > /home/liumy/TEQTL/scRNA/log/scTE_d19_4295.log 2>&1 &
[545414]
```
