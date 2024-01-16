# GREAT

The Genomic Informed Care for Motivating High Risk Individuals Eligible for Evidence-based Prevention (GREAT) framework is a novel approach to incorporate polygenic risk scores (PRS) in the clinic. Here, we describe our analysis procedure for performing principal components analysis (PCA), computing multi-ancestry PRS, and calibrating PRS across the ancestry continuum. We will focus on lung cancer analysis, with application to data from the UK Biobank (coded generally as ukb_hg37). 

## Principal Components Analysis in 1000 Genomes 
Genotype and ancestry data from the 1000 Genomes Project (Phase 3) is publicly available. 

First, we pull the genotype data (reference: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html). Population information has been processed in `1000GP_hg38.sample`.
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

write.table(pca_loadings_clean, 'pca.eigenvec.allele.bim', row.names=F, quote=F)

```

## Computing PC Scores
We apply the PC loadings to the 1000 Genomes data (build 38), and our UK Biobank data (build 37)
```{bash}
plink2 \
--bfile all_hg38_pca \
--score pca_eigenvec.allele.bim 4 8 header \
cols=+scoresums,-scoreavgs,-dosagesum,-nallele \
--score-col-nums 9-13 \
--out 1000G_pca

plink2 \
--bfile ukb_hg37 \
--score pca_eigenvec.allele.bim 2 8 header \
cols=+scoresums,-scoreavgs,-dosagesum,-nallele \
--score-col-nums 9-13 \
--out ukb_pca

```

## Computing and calibrating PRS
We first subset the UK Biobank and 1000 Genomes data to the 101 SNPs used for lung cancer
```{bash}
plink2 \
--bfile all_hg38 \
--extract prs_snps_hg38.txt \
--make-bed \
--out 1000GP_101

plink2 \
--bfile ukb_hg37 \
--extract prs_snps_hg37.txt \
--make-bed \
--out UKB_101

```

We use data from recently-published multi-ancestry genome-wide association studies (GWAS) to compute PRS. Subsetted genotypes are directly read into R, where we can more carefully compute PRS. First, we apply the PRS weights to the 1000 Genomes data so that we can perform ancestry calibration. 
```{R}
library(genio)

prs_weights = read.table('prs_weights.txt', header=T)

# read in genotypes
full_map = read_bim('1000GP_101.bim', verbose=F) %>%
  mutate(chr=as.numeric(chr))
full_fam = read_fam('1000GP_101.fam', verbose=F)
full_G = read_bed('1000GP_101.bed', verbose=F,
                  m_loci = nrow(full_map), n_ind = nrow(full_fam))

# flip genotype to line up with PRS weights
ix_101 = match(prs_weights$rsid_hg38, full_map$id)

map = full_map[ix_101,]
G = full_G[ix_101,]

flip = which(prs_weights$REF != map$ref) 

G[flip,] = 2-G[flip,]

# compute PRS
prs_1000G = c(prs_weights$Beta %*% G)

# ancestry calibration
pc_1000G = read.table('1000G_pca.sscore') %>%
  select(ID=`#FID`, PC1=SCORE1_SUM, PC2=SCORE2_SUM, PC3=SCORE3_SUM, PC4=SCORE4_SUM, PC5=SCORE5_SUM)

data = data.frame(pc_1000G, prs=prs_1000G)

## remove mean-effects of ancestry
fit_mean =lm(prs ~ PC1+PC2+PC3+PC4+PC5, data=data)
pc5_mean_coefs = signif(coef(fit), 4)

pc5_comb = cbind(1, data.matrix(pc_1000G)[,paste0('PC', 1:5)]) %*% pc5_coefs
prs_1000G_adj = prs_1000G - pc5_comb

## remove variance-effects of ancestry
data$var_1000G_adj = prs_1000G_adj^2
fit_var = lm(var_1000G_adj ~ PC1+PC2+PC3+PC4+PC5, data=data)
pc5_coefs_var = signif(coef(fit_var), 4)

pc5_comb_var = cbind(1, data.matrix(pc_1000G)[,paste0('PC', 1:5)]) %*% pc5_coefs_var
prs_1000G_adj_var = (prs_1000G - pc5_comb) / sqrt(pc5_comb_var)

# check calibration by ancestry
sample_1000G = read.table('1000GP_hg38.sample')

data_1000G = cbind(data[,1:6], 
                   `Raw PRS`=prs_1000G,
                   `Calibrated PRS`=prs_1000G_adj_var) %>%
  left_join(sample_1000G)

data_1000G %>%
  group_by(GROUP) %>%
  summarize(raw_mean=mean(`Raw PRS`), raw_sd=sd(`Raw PRS`),
            cal_mean=mean(`Calibrated PRS`), cal_sd=sd(`Calibrated PRS`))

```

Now we apply results to the UK Biobank:
```{R}
full_map = read_bim('UKB_101.bim', verbose=F) %>%
  mutate(chr=as.numeric(chr))
full_fam = read_fam('UKB_101.fam', verbose=F)
full_G = read_bed('UKB_101.bed', verbose=F,
                  m_loci = nrow(full_map), n_ind = nrow(full_fam))

# flip genotype to line up with PRS weights
ix_101 = match(prs_weights$rsid_hg37, full_map$id)

map = full_map[ix_101,]
G = full_G[ix_101,]

flip = which(prs_weights$REF != map$ref) 

G[flip,] = 2-G[flip,]

# impute with SNP means after flipping
G_means = rowMeans(G, na.rm=T) # impute w/ EUR control mean

for (i in 1:101) {
  G[i, is.na(G[i,])] = G_means[i]
}

# compute PRS
prs_UKB = c(prs_weights$Beta %*% G)

# ancestry calibration
pc_UKB = read.table('ukb_pca.sscore') %>%
  select(ID=`#FID`, PC1=SCORE1_SUM, PC2=SCORE2_SUM, PC3=SCORE3_SUM, PC4=SCORE4_SUM, PC5=SCORE5_SUM)

pc5_comb = cbind(1, data.matrix(pc_UKB)[,paste0('PC', 1:5)]) %*% pc5_coefs
pc5_comb_var = cbind(1, data.matrix(pc_UKB)[,paste0('PC', 1:5)]) %*% pc5_coefs_var

prs_UKB_adj_var = (prs_UKB - pc5_comb) / sqrt(pc5_comb_var)

```

## Risk stratification
For lung cancer, we use the 20th and 80th percentiles to divide 1000 Genomes PRS distribution into three risk groups and compute odds ratios. Since we cannot directly share UK Biobank data, we will use variable placeholders. 
```{R}
percentiles_1000G = percentile(data_1000G$`Calibrated PRS`, c(0.2, 0.8))

data_UKB = data.frame(
  Lung_Cancer,
  Ancestry,
  PRS = prs_UKB_adj_var
)

data_UKB = data_UKB %>%
  mutate(Risk = case_when(
    PRS > percentiles_1000G['80%'] ~ 'VeryHighRisk',
    PRS > percentiles_1000G['20%'] ~ 'HighRisk',
    PRS <= percentiles_1000G['20%'] ~ 'AtRisk'
  ))

# odds ratio for very high risk group 
exp(coef(glm(Lung_Cancer ~ Risk, data_UKB %>% filter(adj_risk != 'HighRisk'), family='binomial')))['riskVeryHighRisk']

# odds ratio for high risk group
exp(coef(glm(Lung_Cancer ~ Risk, data_UKB %>% filter(adj_risk != 'VeryHighRisk'), family='binomial')))['riskHighRisk']

```
