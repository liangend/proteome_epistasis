library(data.table)
library(sandwich)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait.csv')
sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter_w_cis.txt')
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

## regression between phenotype and interaction
add_coef = list()
for (n in 1:nrow(sig_inter)) {
  target_n = sig_inter$target[n]
  reg_n = sig_inter$regulator[n]
  
  target_pred = fread(paste0('cis_trans_inter/09_inter_pheno_asso/pred_prot/',
                             target_n, '_pred.txt'))
  reg_pred = fread(paste0('cis_trans_inter/09_inter_pheno_asso/pred_prot/',
                          reg_n, '_pred.txt'))
  
  ## remove ids used in FUSION training
  target_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                                 target_n, '.txt'))
  reg_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                              reg_n, '.txt'))
  train_id = unique(c(target_train_id$FID, reg_train_id$FID))
  common_id2 = target_pred$id[!target_pred$id %in% train_id]
  
  res_n = res_i[match(common_id2, common_id)]
  target_n = target_pred$pred[match(common_id2, target_pred$id)]
  reg_n = reg_pred$pred[match(common_id2, reg_pred$id)]
  
  fit2 = lm(res_n ~ target_n)
  fit2_coef = summary(fit2)$coefficients
  sandwich_se2 = diag(vcovHC(fit2, type = "HC"))^0.5
  sandwich_t2 = fit2_coef[, 1]/sandwich_se2
  sandwich_p2 = pchisq(sandwich_t2^2, 1, lower.tail = FALSE)
  fit2_coef = cbind(fit2_coef, sandwich_t = sandwich_t2, sandwich_p = sandwich_p2)
  
  add_coef[[n]] = fit2_coef
}
names(add_coef) = sig_inter$target

saveRDS(add_coef, paste0('cis_trans_inter/10_add_pheno_asso/asso_all_ind/', 
                         trait_i, '.rds'))


