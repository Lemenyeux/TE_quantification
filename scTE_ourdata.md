# 1. Fastq to bam file using Cell Ranger
``` bash
nohup cellranger count --id=d19_4295 \
   --fastqs=/home/liumy/TEQTL/scRNA/raw_data/testing_sample \
   --sample=D19-4295 \
   --output-dir=/home/liumy/TEQTL/scRNA/results/d19_4295 \
   --transcriptome=/home/liumy/software/cellranger/refdata-gex-GRCh38-2024-A \
   --create-bam=true > /home/liumy/TEQTL/scRNA/results/d19_4295/fastqtobam_d19_4295.log 2>&1 &

# 536650 [240609 17:00-19:22]
```

# 2. Exclude lines with blank CB
``` bash
samtools view possorted_genome_bam.bam -h | awk '/^@/ || /CB:/' | samtools view -h -b > possorted_genome_bam.clean.bam
```

# 3. Setup

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
```

# 4. Analysis of 10x style scRNA-seq data
``` bash
cd /home/liumy/TEQTL/scRNA/results/d19_4295_te
nohup scTE -i possorted_genome_bam.clean.bam -o d19_4295 -x /home/liumy/software/scTE/hg38.inclusive.idx -p 10 -CB CB -UMI UB > /home/liumy/TEQTL/scRNA/log/scTE_d19_4295.log 2>&1 &

[544955 240610-11:27-13:16]
[10000 obs. Of 58390 variables]

nohup scTE -i possorted_genome_bam.clean.bam -o d19_4295 -x /home/liumy/software/scTE/hg38.inclusive.idx -p 10 --hdf5 True -CB CB -UMI UB --expect-cell  > /home/liumy/TEQTL/scRNA/log/scTE_d19_4295.log 2>&1 &
[545414]
```

# 5. Output
# 5.1 results from Cell Ranger - filtered_feature_bc_matrix
``` bash
cd /home/liumy/TEQTL/scRNA/scTE/results/d19_4295/outs/filtered_feature_bc_matrix
zcat barcodes.tsv.gz | head
zcat features.tsv.gz | head
zcat matrix.mtx.gz | head
```

**barcodes.tsv.gz: 280 cells**
``` markdown
AAACCCAAGAACGCGT-1
AAACCCAAGTACTCGT-1
AAACCCACAAGTGTCT-1
AAACCCACAATGAACA-1
AAACCCAGTATAGCTC-1
AAACCCAGTCACTCTC-1
AAACCCAGTGGTAACG-1
AAACCCATCATGAAAG-1
AAACCCATCCACCTGT-1
AAACCCATCCTTCACG-1
```

**features.tsv.gz: 1328 feature**
``` markdown
ENSG00000290825	DDX11L2	Gene Expression
ENSG00000243485	MIR1302-2HG	Gene Expression
ENSG00000237613	FAM138A	Gene Expression
ENSG00000290826	ENSG00000290826	Gene Expression
ENSG00000186092	OR4F5	Gene Expression
ENSG00000238009	ENSG00000238009	Gene Expression
ENSG00000239945	ENSG00000239945	Gene Expression
ENSG00000239906	ENSG00000239906	Gene Expression
ENSG00000241860	ENSG00000241860	Gene Expression
ENSG00000241599	ENSG00000241599	Gene Expression
```

**matrix.mtx.gz**
``` markdown
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-8.0.1", "format_version": 2}
38606 10538 29034127
23 1 1
52 1 1
81 1 1
84 1 1
95 1 1
99 1 1
205 1 1
```

## 5.2 input of scTE 
**rmsk.txt.gz: 5683690 feautres**
`zcat ~/TEQTL/scRNA/scTE/other_data/rmsk.txt.gz | head`
``` markdown
585  463   13   6    17  chr1  10000  10468  -248945954  +  (TAACCC)n  Simple_repeat  Simple_repeat  1      471   0     1
585  3612  114  215  13  chr1  10468  11447  -248944975  -  TAR1       Satellite      telo           -399   1712  483   2
585  484   251  132  0   chr1  11504  11675  -248944747  -  L1MC5a     LINE           L1             -2382  395   199   3
585  239   294  19   10  chr1  11677  11780  -248944642  -  MER5B      DNA            hAT-Charlie    -74    104   1     4
585  318   230  37   0   chr1  15264  15355  -248941067  -  MIR3       SINE           MIR            -119   143   49    5
585  18    232  0    19  chr1  15797  15849  -248940573  +  (TGCTCC)n  Simple_repeat  Simple_repeat  1      52    0     6
585  18    137  0    0   chr1  16712  16744  -248939678  +  (TGG)n     Simple_repeat  Simple_repeat  1      32    0     7
585  239   338  129  0   chr1  18906  19048  -248937374  +  L2a        LINE           L2             2942   3104  -322  8
585  994   312  60   25  chr1  19971  20405  -248936017  +  L3         LINE           CR1            2680   3129  -970  9
585  270   331  7    27  chr1  20530  20679  -248935743  +  Plat_L3    LINE           CR1            2802   2947  -639  1
```
![image](https://github.com/Lemenyeux/TE_quantification/assets/87812974/5635beef-46ce-406c-8e4e-7dc8c332b7a9)

# 5.3 results from scTE
**d19_4295.csv: 10000 cells**
``` python
import anndata
import numpy as np
import pandas as pd

# Configure logging to write to a file (optional)
import logging
logging.basicConfig(filename='anndata.log', level=logging.INFO)

# Set verbosity level to errors only for Scanpy
import scanpy as sc
sc.settings.verbosity = 0

# Read the .h5ad file
adata = anndata.read_h5ad('d19_4295.h5ad', backed='r')
print(adata)
# AnnData object with n_obs × n_vars = 10000 × 58389 backed at 'd19_4295.h5ad'

# Load TE names
te_feature = pd.read_csv('/home/liumy/TEQTL/scRNA/scTE/other_data/rmsk.txt.gz', compression='gzip', header=None, sep='\t')
# Extract the 11th column (note that pandas uses 0-based indexing, so the 11th column is at index 10)
te_names = te_feature.iloc[:, 10].tolist()

# Ensure that the genes_to_select are in the var_names of the AnnData object
te_names = [gene for gene in te_names if gene in adata.var_names]
print(len(te_names))
print(te_names[:5])
# ['TAR1', 'L1MC5a', 'MER5B', 'MIR3', 'L2a']

# Subset the AnnData object to include only selected TE names
adata_subset = adata[:, te_names]
adata_subset

# View of AnnData object with n_obs × n_vars = 10000 × 4819290 backed at 'd19_4295.h5ad'

# Save the subset AnnData object to an h5ad file
adata_subset.write('/home/liumy/TEQTL/scRNA/scTE/results/d19_4295_te/d19_4295_te.h5ad')

# Convert the expression matrix to a DataFrame and save to a CSV file
expression_matrix_df = pd.DataFrame(adata_subset.X, index=adata_subset.obs_names, columns=adata_subset.var_names)
expression_matrix_df.to_csv('/home/liumy/TEQTL/scRNA/scTE/results/d19_4295_te/d19_4295_te.csv')

```

