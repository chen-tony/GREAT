# GREAT

The Genomic Informed Care for Motivating High Risk Individuals Eligible for Evidence-based Prevention (GREAT) framework is a novel approach to incorporate polygenic risk scores (PRS) in the clinic. Here, we describe our analysis procedure for performing principal components analysis (PCA), computing multi-ancestry PRS, and calibrating PRS across the ancestry continuum. 

## Principal Components Analysis in 1000 Genomes 
Genotype and ancestry data from the 1000 Genomes Project (Phase 3) is publicly available. 

First, we pull the genotype data (reference: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html):
```{bash}
wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
mv 'all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen

wget https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1
mv 'all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1' all_hg38.pvar.zst

gwet https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1
mv 'hg38_corrected.psam?dl=1' all_hg38.psam

```

For population information, click "Download List" from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38. 

Next, we will convert the data to build37 using liftOver:
```{bash}

```

Then, we will perform PCA using a chosen set of SNPs:
```{bash}
plink2 \
--bfile \
--pca allele-wts 5 \
--out 1000G_pca
```

Finally, we will prepare the output for computing PC scores in external data
```{R}
library(dplyr)
library(data.table)
pca_loadings = 

```

## Computing PRS
We use data from recently-published multi-ancestry genome-wide association studies (GWAS) to compute PRS.
```{R}

```

## Calibrating PRS
```{R}

```

## Risk stratification
```{R}

```
