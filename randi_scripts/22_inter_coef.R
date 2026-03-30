library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')

sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter_w_cis.txt')
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')
st_cov = na.omit(st_cov)
coef_all1 = list()
coef_all2 = list()
for (i in 1:nrow(sig_inter)) {
  target_i = sig_inter$target[i]
  reg_i = sig_inter$regulator[i]
  target_pred = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
			     target_i, '_good_pred.txt'))
  reg_pred = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
			  reg_i, '_good_pred.txt'))
  target_pred$reg_pred = reg_pred$pred[match(target_pred$id, reg_pred$id)]
  target_pred = na.omit(target_pred)
  common_id = intersect(target_pred$id, st_cov$V1)
  target_pred = target_pred[match(common_id, target_pred$id), ]
  st_cov_sub = st_cov[match(common_id, st_cov$V1), -(1:2)]
  
  target_res = residuals(lm(target_pred$obs ~ as.matrix(st_cov_sub)))
  target_pred$pred = target_pred$pred - mean(target_pred$pred)
  target_pred$reg_pred = target_pred$reg_pred - mean(target_pred$reg_pred)
  model_fit1 = lm(target_res ~ target_pred$pred + target_pred$reg_pred)
  model_fit2 = lm(target_res ~ target_pred$pred * target_pred$reg_pred)
  coef_all1[[i]] = summary(model_fit1)$coefficients
  coef_all2[[i]] = summary(model_fit2)$coefficients
  names(coef_all1)[i] = paste0(target_i, '+', reg_i)
  names(coef_all2)[i] = paste0(target_i, '*', reg_i)
}
saveRDS(coef_all1, 'cis_trans_inter/12_pheno_asso_pred_w_inter/add_coef_after_center.rds')
saveRDS(coef_all2, 'cis_trans_inter/12_pheno_asso_pred_w_inter/inter_coef_after_center.rds')
