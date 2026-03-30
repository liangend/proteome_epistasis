library(data.table)
library(sandwich)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait_share.txt')
sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter_w_cis.txt')

inter_coef = readRDS('cis_trans_inter/12_pheno_asso_pred_w_inter/inter_coef_after_center.rds')
add_coef = readRDS('cis_trans_inter/12_pheno_asso_pred_w_inter/add_coef_after_center.rds')

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
add1_pwas = list()
add2_pwas = list()
add1_inter_pwas = list()
add2_inter_pwas = list()

for (n in 1:nrow(sig_inter)) {
  inter_coef_n = inter_coef[[n]]
  add_coef_n = add_coef[[n]]
  
  target_n = sig_inter$target[n]
  reg_n = sig_inter$regulator[n]
  
  target_pred = fread(paste0('cis_trans_inter/09_inter_pheno_asso/pred_prot/',
                             target_n, '_pred.txt'))
  reg_pred = fread(paste0('cis_trans_inter/09_inter_pheno_asso/pred_prot/',
                          reg_n, '_pred.txt'))
  
  ## remove ids used in FUSION training and coef estimating
  target_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                                 target_n, '.txt'))
  reg_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                              reg_n, '.txt'))
  target_coef_est = fread(paste0('cis_trans_inter/05_fusion/pred_expr/',
                                 target_n, '_good_pred.txt'))
  reg_coef_est = fread(paste0('cis_trans_inter/05_fusion/pred_expr/',
                              reg_n, '_good_pred.txt'))
  rm_id = unique(c(target_train_id$FID, reg_train_id$FID,
                   target_coef_est$id, reg_coef_est$id))
  common_id2 = target_pred$id[!target_pred$id %in% rm_id]
  
  res_n = res_i[match(common_id2, common_id)]
  target_n = target_pred$pred[match(common_id2, target_pred$id)]
  target_n = target_n - mean(target_n)
  reg_n = reg_pred$pred[match(common_id2, reg_pred$id)]
  reg_n = reg_n - mean(reg_n)

  ## use cis-predicted target directly
  fit0 = lm(res_n ~ target_n)
  fit0_coef = summary(fit0)$coefficients
  
  ## predict target expression using two addivite coef
  target_est1 = add_coef_n[1,1] + add_coef_n[2,1] * target_n +
    add_coef_n[3,1] * reg_n
  fit1 = lm(res_n ~ target_est1)
  fit1_coef = summary(fit1)$coefficients
  
  ## predict target expression using one additive and one interaction coef
  target_est2 = inter_coef_n[1,1] + inter_coef_n[2,1] * target_n +
    inter_coef_n[4,1] * target_n * reg_n
  fit2 = lm(res_n ~ target_est2)
  fit2_coef = summary(fit2)$coefficients
  
  ## predict target expression using both additive and interaction coef
  target_est3 = inter_coef_n[1,1] + inter_coef_n[2,1] * target_n + 
    inter_coef_n[3,1] * reg_n + inter_coef_n[4,1] * target_n * reg_n
  
  fit3 = lm(res_n ~ target_est3)
  fit3_coef = summary(fit3)$coefficients
  
  add1_pwas[[n]] = fit0_coef
  add2_pwas[[n]] = fit1_coef
  add1_inter_pwas[[n]] = fit2_coef
  add2_inter_pwas[[n]] = fit3_coef
}
names(add1_pwas) = sig_inter$target
names(add2_pwas) = sig_inter$target
names(add1_inter_pwas) = sig_inter$target
names(add2_inter_pwas) = sig_inter$target

saveRDS(add1_pwas, paste0('cis_trans_inter/12_pheno_asso_pred_w_inter/center/asso_1_add/',
                          trait_i, '.rds'))
saveRDS(add2_pwas, paste0('cis_trans_inter/12_pheno_asso_pred_w_inter/center/asso_2_add/',
                          trait_i, '.rds'))
saveRDS(add1_inter_pwas, paste0('cis_trans_inter/12_pheno_asso_pred_w_inter/center/asso_1_add_inter/',
                                trait_i, '.rds'))
saveRDS(add2_inter_pwas, paste0('cis_trans_inter/12_pheno_asso_pred_w_inter/center/asso_2_add_inter/', 
                                trait_i, '.rds'))



