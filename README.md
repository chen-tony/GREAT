# GREAT

The Genomic Informed Care for Motivating High Risk Individuals Eligible for Evidence-based Prevention (GREAT) framework is a novel approach to incorporate polygenic risk scores (PRS) in the clinic. Here, we describe our analysis procedure for performing principal components analysis (PCA), computing multi-ancestry PRS, and calibrating PRS across the ancestry continuum. We will focus on lung cancer analysis, with application to data from the UK Biobank (coded generally as ukb_hg37). 

## Principal Components Analysis in 1000 Genomes 
Genotype and ancestry data from the 1000 Genomes Project (Phase 3) is publicly available. 

First, we pull the genotype data (reference: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html). For population information, click "Download the list" under "3202 samples" from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38. 
```{bash}
wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
mv 'all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen

wget https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1
mv 'all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1' all_hg38.pvar.zst

gwet https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1
mv 'hg38_corrected.psam?dl=1' all_hg38.psam

```

Next, we will perform PCA using a chosen set of SNPs:
```{bash}
plink2 \
--bfile all_hg38 \
--extract pca_snps_hg38.txt \
--make-bed \
--out all_hg38_pca

plink2 \
--bfile all_hg38_pca \
--pca allele-wts 5 \
--out 1000G_pca

```

Finally, we will prepare the PCA output for computing PC scores in external data
```{R}
library(dplyr)
library(data.table)

pca_snp_map = fread('pca_snps.txt')
pca_loadings = fread('1000G_pca.eigenvec.allele')

# plink gives duplicate rows for either reference or alternate allele as the effect allele
# we will just use the case where the alternate allele is the effect allele
pca_loadings_clean = pca_loadings %>%
filter(ALT == A1) %>%
left_join(pca_snp_map, by=c('#CHROM'='Chromosome', 'ID'='SNP (build 38)')) %>%
select(CHR=`#CHROM`, SNP37=`SNP (build 37)`, POS37=`Position (build 37)`, SNP38=`SNP (build 38)`, POS38=`Position (build 38)`, REF, ALT, A1, paste0('PC', 1:5))

write.table(pca_loadings_clean, '1000G_pca.eigenvec.allele.bim', row.names=F, quote=F)

```

## Computing PC Scores
We apply the PC loadings to the 1000 Genomes data (build 38), and our UK Biobank data (build 37)
```{bash}
plink2 \
--bfile all_hg38_pca \
--score 1000G_pca_eigenvec.allele.bim 4 8 header \
cols=+scoresums,-scoreavgs,-dosagesum,-nallele \
--score-col-nums 9-13 \
--out 1000G_pca

plink2 \
--bfile ukb_hg37 \
--score 1000G_pca_eigenvec.allele.bim 2 8 header \
cols=+scoresums,-scoreavgs,-dosagesum,-nallele \
--score-col-nums 9-13 \
--out ukb_pca

```

## Computing PRS
We first subset the UK Biobank and 1000 Genomes data to the 101 SNPs used for lung cancer
```{bash}
plink2 \
--bfile all_hg38 \
--extract prs_snps_hg38.txt \
--make-bed \
--out 1000G_101

plink2 \
--bfile ukb_hg37 \
--extract prs_snps_hg37.txt \
--make-bed \
--out UKB_101

```

We use data from recently-published multi-ancestry genome-wide association studies (GWAS) to compute PRS. Subsetted genotypes are directly read into R, where we can more carefully compute PRS. 
```{R}
library(genio)

```

## Calibrating PRS
```{R}


```

## Risk stratification
```{R}


```
