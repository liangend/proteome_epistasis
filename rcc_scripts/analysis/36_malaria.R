library(data.table)
library(ggplot2)
library(plink2R)
setwd('/project/xuanyao/jinghui')

malaria = fread('/project2/xuanyao/jinghui/MalariaGEN/MalariaGEN_GWAS.sample')
malaria = malaria[-1, c(1:3, 7, 9:11, 13, 16:21)]
malaria$severe_malaria = as.numeric(malaria$severe_malaria)

malaria$rs334_dosage = 0
malaria$rs334_dosage[malaria$rs334_genotype == 'AT'] = 1
malaria$rs334_dosage[malaria$rs334_genotype == 'AA'] = 2
malaria$popn_PC1 = as.numeric(malaria$popn_PC1)
malaria$popn_PC2 = as.numeric(malaria$popn_PC2)
malaria$popn_PC3 = as.numeric(malaria$popn_PC3)
malaria$popn_PC4 = as.numeric(malaria$popn_PC4)
malaria$popn_PC5 = as.numeric(malaria$popn_PC5)
malaria = na.omit(malaria)

fit0 = glm(severe_malaria ~ country + clinical_sex + popn_PC1 +
             popn_PC2 + popn_PC3 + popn_PC4 + popn_PC5, 
           family = binomial(), data = malaria)
malaria$severe_res = residuals(fit0)

gene_meta = fread('mediation_pathway/00_ref/genecode.GRCh38.gene.meta.gtf')

## MBL2 * ABO effect on malaria


#### ABO
## ABO FUSION SNPs
load('cis_trans_inter/06_fusion_model/ukb_eur/ABO.wgt.RDat')
abo_wgt = wgt.matrix[wgt.matrix[,1] != 0, ]
abo_wgt = abo_wgt[order(abs(abo_wgt[,1]), decreasing = T), ]
abo_wgt = as.data.frame(abo_wgt)
abo_cis_geno = read_plink('cis_trans_inter/00_ref/cis_geno/ABO',
                          impute="avg")
abo_wgt$bp = abo_cis_geno$bim$V4[match(rownames(abo_wgt), 
                                       abo_cis_geno$bim$V2)]
## ABO SNPs in MalariaGen
abo_chr = gene_meta$chr[which(gene_meta$gene_name == 'ABO')]
abo_geno = read_plink(paste0('/project2/xuanyao/jinghui/MalariaGEN/genotype/geno_by_chr/',
                             abo_chr), impute="avg")
abo_snp = abo_geno$bim$V2[which(abo_geno$bim$V4 %in% abo_wgt$bp)]
abo_geno = abo_geno$bed[, abo_snp]
rownames(abo_geno) = sapply(strsplit(rownames(abo_geno), ':', fixed = T), 
                            '[', 1)
abo_geno = abo_geno[match(malaria$chip_id, rownames(abo_geno)), ]

#### MBL2
## MBL2 FUSION SNPs
load('cis_trans_inter/06_fusion_model/ukb_eur/MBL2.wgt.RDat')
mbl2_wgt = wgt.matrix[wgt.matrix[,1] != 0, ]
mbl2_wgt = mbl2_wgt[order(abs(mbl2_wgt[,1]), decreasing = T), ]
mbl2_wgt = as.data.frame(mbl2_wgt)
mbl2_cis_geno = read_plink('cis_trans_inter/00_ref/cis_geno/MBL2',
                           impute="avg")
mbl2_wgt$bp = mbl2_cis_geno$bim$V4[match(rownames(mbl2_wgt), 
                                         mbl2_cis_geno$bim$V2)]
## MBL2 SNPs in MalariaGen
mbl2_chr = gene_meta$chr[which(gene_meta$gene_name == 'MBL2')]
mbl2_geno = read_plink(paste0('/project2/xuanyao/jinghui/MalariaGEN/genotype/geno_by_chr/',
                             mbl2_chr), impute="avg")
mbl2_snp = mbl2_geno$bim$V2[which(mbl2_geno$bim$V4 %in% mbl2_wgt$bp)]
mbl2_geno = mbl2_geno$bed[, mbl2_snp]
rownames(mbl2_geno) = sapply(strsplit(rownames(mbl2_geno), ':', fixed = T), 
                            '[', 1)
mbl2_geno = mbl2_geno[match(malaria$chip_id, rownames(mbl2_geno)), ]

#### CD209
## CD209 FUSION SNPs
load('cis_trans_inter/06_fusion_model/ukb_eur/CD209.wgt.RDat')
cd209_wgt = wgt.matrix[wgt.matrix[,1] != 0, ]
cd209_wgt = cd209_wgt[order(abs(cd209_wgt[,1]), decreasing = T), ]
cd209_wgt = as.data.frame(cd209_wgt)
cd209_cis_geno = read_plink('cis_trans_inter/00_ref/cis_geno/CD209',
                            impute="avg")
cd209_wgt$bp = cd209_cis_geno$bim$V4[match(rownames(cd209_wgt), 
                                           cd209_cis_geno$bim$V2)]
## CD209 SNPs in MalariaGen
cd209_chr = gene_meta$chr[which(gene_meta$gene_name == 'CD209')]
cd209_geno = read_plink(paste0('/project2/xuanyao/jinghui/MalariaGEN/genotype/geno_by_chr/',
                              cd209_chr), impute="avg")
cd209_snp = cd209_geno$bim$V2[which(cd209_geno$bim$V4 %in% cd209_wgt$bp)]
cd209_geno = cd209_geno$bed[, cd209_snp]
rownames(cd209_geno) = sapply(strsplit(rownames(cd209_geno), ':', fixed = T), 
                             '[', 1)
cd209_geno = cd209_geno[match(malaria$chip_id, rownames(cd209_geno)), ]


#### ABO * MBL2 effect on malaria
abo_mbl2_inter = matrix(0, nrow = ncol(mbl2_geno), ncol = ncol(abo_geno))
for (i in 1:ncol(abo_geno)) {
  abo_i = abo_geno[, i]
  for (j in 1:ncol(mbl2_geno)) {
    mbl2_j = mbl2_geno[, j]
    fit_ij = lm(malaria$severe_res ~ abo_i * mbl2_j)
    coef_ij = summary(fit_ij)$coefficients
    abo_mbl2_inter[j, i] = coef_ij[4, 4]
  }
}

#### ABO * CD209 effect on malaria
abo_cd209_inter = matrix(0, nrow = ncol(cd209_geno), ncol = ncol(abo_geno))
for (i in 1:ncol(abo_geno)) {
  abo_i = abo_geno[, i]
  for (j in 1:ncol(cd209_geno)) {
    cd209_j = cd209_geno[, j]
    fit_ij = lm(malaria$severe_res ~ abo_i * mbl2_j)
    coef_ij = summary(fit_ij)$coefficients
    abo_cd209_inter[j, i] = coef_ij[4, 4]
  }
}


