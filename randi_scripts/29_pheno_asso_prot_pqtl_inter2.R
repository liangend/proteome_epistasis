library(data.table)
library(sandwich)
library(plink2R)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

genotype_dir = 'cis_trans_inter/09_inter_pheno_asso/all_white_geno/'
plink_dir = 'software/'

ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait_share.txt') # 56 traits
sig_inter = fread('cis_trans_inter/16_prot_trans_pqtl_inter/inter_all.txt')
sig_inter = sig_inter[sig_inter$qval < 0.05, ]
all_cov = fread('cis_trans_inter/00_ref/all_norm_cov.txt')
all_cov = na.omit(all_cov)

## phenotype
trait_i = ukb_trait$abbrev[i]
file_i = ukb_trait$file[i]
type_i = ukb_trait$Type[i]

pheno_i = fread(paste0('/gpfs/data/ukb-share/extracted_phenotypes/', trait_i, 
                       '/', file_i))
pheno_i = na.omit(pheno_i)
colnames(pheno_i)[3] = 'pheno'
## remove negative values representing missing data
pheno_i = pheno_i[pheno_i$pheno >= 0, ]

if (type_i == 'Continuous') {
  pheno_i$pheno = (pheno_i$pheno - mean(pheno_i$pheno)) / sd(pheno_i$pheno)
} else {
  pheno_i$pheno[pheno_i$pheno > 0] = 1
}

common_id = intersect(all_cov$FID, pheno_i$FID)
pheno_i = pheno_i[match(common_id, pheno_i$FID), ]
all_cov = all_cov[match(common_id, all_cov$FID), -(1:2)]

## residualize phenotypes by regressing on covariates
if (type_i == 'Continuous') {
  res_i = residuals(lm(unlist(pheno_i$pheno) ~ as.matrix(all_cov)))
} else {
  res_i = residuals(glm(unlist(pheno_i$pheno) ~ as.matrix(all_cov), family = binomial))
}

## regression between phenotype and interaction
reg_coef = list()
for (n in 1:nrow(sig_inter)) {
  target_n = sig_inter$prot[n]
  target_pred = fread(paste0('cis_trans_inter/16_prot_trans_pqtl_inter/pred_prot/',
                             target_n, '_pred.txt'))
  chr_n = sig_inter$chr[n]
  if (!file.exists(paste0("cis_trans_inter/06_fusion_inter/tmp/", target_n, '_', n, ".bed"))) {
    fwrite(data.frame(chr = chr_n, bp1 = sig_inter$bp_hg19[n],
                      bp2 = sig_inter$bp_hg19[n], bp3 = sig_inter$bp_hg19[n]), 
           paste0("cis_trans_inter/06_fusion_inter/tmp/", target_n, '_', n, ".txt"), 
           col.names = F, sep = '\t')
    system(paste0(plink_dir, "plink --bfile ", genotype_dir, "chr", chr_n,
                  "_QC --extract cis_trans_inter/06_fusion_inter/tmp/", target_n, '_', n, 
                  ".txt --make-bed --range --thread-num 16 --out cis_trans_inter/06_fusion_inter/tmp/", 
                  target_n, '_', n), intern=T)
    if (!file.exists(paste0("cis_trans_inter/06_fusion_inter/tmp/", target_n, '_', n, ".bed"))) {
      reg_coef[[n]] = NA
      next
    }
  }
  geno_n = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/', 
                             target_n, '_', n), impute="avg")
  
  ## remove ids used in FUSION training
  target_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                                 target_n, '.txt'))
  common_id2 = target_pred$id[!target_pred$id %in% target_train_id$FID]
  
  res_n = res_i[match(common_id2, common_id)]
  pred_n = target_pred$pred[match(common_id2, target_pred$id)]
  pqtl_n = geno_n$bed[match(common_id2, geno_n$fam$V1), 1]
  
  fit1 = lm(res_n ~ pred_n * pqtl_n)
  fit_coef = summary(fit1)$coefficients
  sandwich_se = diag(vcovHC(fit1, type = "HC"))^0.5
  sandwich_t = fit_coef[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  fit_coef = cbind(fit_coef, sandwich_t = sandwich_t, sandwich_p = sandwich_p)
  
  reg_coef[[n]] = fit_coef
}
names(reg_coef) = 1:nrow(sig_inter)
saveRDS(reg_coef, paste0('cis_trans_inter/16_prot_trans_pqtl_inter/asso_pheno/', 
                         trait_i, '.rds'))




