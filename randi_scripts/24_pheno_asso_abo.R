library(data.table)
library(sandwich)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait.csv')
all_cov = fread('cis_trans_inter/00_ref/all_norm_cov.txt')
all_cov = na.omit(all_cov)

## phenotype
trait_i = ukb_trait$abbrev[i]
file_i = list.files(paste0('cis_trans_inter/00_ref/ukb_traits/', trait_i))
pheno_i = fread(paste0('cis_trans_inter/00_ref/ukb_traits/', trait_i, 
                       '/', file_i))
type_i = ukb_trait$Type[i]
pheno_i = na.omit(pheno_i)
colnames(pheno_i)[3] = 'pheno'

common_id = intersect(all_cov$FID, pheno_i$FID)
pheno_i = pheno_i[match(common_id, pheno_i$FID), ]
all_cov = all_cov[match(common_id, all_cov$FID), -(1:2)]

## residualize phenotypes by regressing on covariates
if (type_i == 'Continuous') {
  pheno_i$pheno = (pheno_i$pheno - mean(pheno_i$pheno)) / sd(pheno_i$pheno)
  res_i = residuals(lm(unlist(pheno_i$pheno) ~ as.matrix(all_cov)))
} else {
  if (trait_i == 'Bipolar') {
    pheno_i$pheno[pheno_i$pheno > 0] = 1
  } else if (trait_i == 'Hypertension'){
    pheno_i$pheno[pheno_i$pheno == 1065] = 1
    pheno_i$pheno[pheno_i$pheno != 1065] = 0
  }
  res_i = residuals(glm(unlist(pheno_i$pheno) ~ as.matrix(all_cov), family = binomial))
}

## regression between phenotype and ABO
target_n = 'ABO'
target_pred = fread(paste0('cis_trans_inter/09_inter_pheno_asso/pred_prot/',
                           target_n, '_pred.txt'))

## remove ids used in FUSION training and coef estimating
target_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                               target_n, '.txt'))
target_coef_est = fread(paste0('cis_trans_inter/05_fusion/pred_expr/',
                               target_n, '_good_pred.txt'))

rm_id = unique(c(target_train_id$FID, target_coef_est$id))
common_id2 = target_pred$id[!target_pred$id %in% rm_id]

res_n = res_i[match(common_id2, common_id)]
target_n = target_pred$pred[match(common_id2, target_pred$id)]

## use cis-predicted target directly
fit0 = lm(res_n ~ target_n)
fit0_coef = summary(fit0)$coefficients

saveRDS(fit0_coef, paste0('cis_trans_inter/13_pheno_asso_abo/',
                          trait_i, '.rds'))

