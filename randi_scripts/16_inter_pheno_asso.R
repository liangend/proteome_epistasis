library(data.table)
library(sandwich)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait.csv')
sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter.txt')
ukb_cov = fread('cis_trans_inter/00_ref/gcta_cov_less.qcovar')

## phenotype
trait_i = ukb_trait$abbrev[i]
file_i = list.files(paste0('cis_trans_inter/00_ref/ukb_traits/', trait_i))
pheno_i = fread(paste0('cis_trans_inter/00_ref/ukb_traits/', trait_i, 
                       '/', file_i))
type_i = ukb_trait$Type[i]
pheno_i = na.omit(pheno_i)
colnames(pheno_i)[3] = 'pheno'

common_id = intersect(ukb_cov$V1, pheno_i$FID)
pheno_i = pheno_i[match(common_id, pheno_i$FID), ]
ukb_cov = ukb_cov[match(common_id, ukb_cov$V1), -(1:2)]


## residualize phenotypes by regressing on covariates
if (type_i == 'Continuous') {
  pheno_i$pheno = (pheno_i$pheno - mean(pheno_i$pheno)) / sd(pheno_i$pheno)
  res_i = residuals(lm(unlist(pheno_i$pheno) ~ as.matrix(ukb_cov)))
} else {
  if (trait_i == 'Bipolar') {
    pheno_i$pheno[pheno_i$pheno > 0] = 1
  } else if (trait_i == 'Hypertension'){
    pheno_i$pheno[pheno_i$pheno == 1065] = 1
    pheno_i$pheno[pheno_i$pheno != 1065] = 0
  }
  res_i = residuals(glm(unlist(pheno_i$pheno) ~ as.matrix(ukb_cov), family = binomial))
}

## regression between phenotype and interaction
reg_coef = list()
n_sample = c()
for (n in 1:nrow(sig_inter)) {
  target_n = sig_inter$target[n]
  reg_n = sig_inter$regulator[n]
  
  target_pred = fread(paste0('cis_trans_inter/05_fusion/pred_expr/',
                             target_n, '_good_pred.txt'))
  reg_pred = fread(paste0('cis_trans_inter/05_fusion/pred_expr/',
                          reg_n, '_good_pred.txt'))
  common_id2 = intersect(target_pred$id, reg_pred$id)
  res_n = res_i[match(common_id2, common_id)]
  target_n = target_pred$pred[match(common_id2, target_pred$id)]
  reg_n = reg_pred$pred[match(common_id2, reg_pred$id)]
  
  fit1 = lm(res_n ~ target_n * reg_n)
  fit_coef = summary(fit1)$coefficients
  sandwich_se = diag(vcovHC(fit1, type = "HC"))^0.5
  sandwich_t = fit_coef[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  fit_coef = cbind(fit_coef, sandwich_t = sandwich_t, sandwich_p = sandwich_p)

  reg_coef[[n]] = fit_coef
  n_sample[n] = length(summary(fit1)$residuals)
}
names(reg_coef) = paste0(sig_inter$target, '*', sig_inter$regulator)
reg_coef$n_sample = n_sample
saveRDS(reg_coef, paste0('cis_trans_inter/09_inter_pheno_asso/asso_ukb_ppp_ind/', trait_i, '.rds'))


