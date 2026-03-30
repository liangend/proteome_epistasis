setwd('/project/xuanyao/jinghui')
library(data.table)

args = commandArgs(trailingOnly = T)
target_index = as.numeric(args[1])

## pairwise interactions to be tested
inter_test = fread('cis_trans_inter/07_prot_inter/prot_perm_test.txt')
uniq_target = unique(inter_test$target)
target_prot = uniq_target[target_index]
inter_test = inter_test[inter_test$target == target_prot, ]

# standardized covariates
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

# target protein prediction and observation
target_pred = fread(paste0('cis_trans_inter/06_fusion_pred/', 
                           target_prot, '_good_pred.txt'))
target_pred$pred = (target_pred$pred - mean(target_pred$pred))/
  sd(target_pred$pred)
target_pred = na.omit(target_pred)

reg_add_p = c()
inter_p = c()
for(i in 1:nrow(inter_test)){
  reg_i = inter_test$regulator[i]

  # cis prediction of regulator i
  reg_pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/', 
                            reg_i, '_good_pred.txt'))
  reg_pred_i$pred = (reg_pred_i$pred - mean(reg_pred_i$pred))/
    sd(reg_pred_i$pred)
  common_id = intersect(reg_pred_i$id, target_pred$id)
  
  target_ex_i = target_pred$obs[match(common_id, target_pred$id)]
  target_pred_i = target_pred$pred[match(common_id, target_pred$id)]
  reg_pred_i = reg_pred_i$pred[match(common_id, reg_pred_i$id)]
  
  st_cov_i = st_cov[match(common_id, st_cov$V1), ]
  st_cov_i = as.matrix(st_cov_i[,-(1:2)])
  
  ## get the residualized prot expr
  lm_res = lm(target_ex_i ~ st_cov_i)
  res_i = residuals(lm_res)
  
  ## permute regulator for 1000 times
  reg_add_p_i = c()
  inter_p_i = c()
  for (rep in 1:1000) {
    reg_pred_shuf = sample(reg_pred_i)
    test = lm(res_i ~ target_pred_i + reg_pred_shuf + target_pred_i * reg_pred_shuf)
    test_coeff = summary(test)$coefficients
    reg_add_p_i[rep] = test_coeff['reg_pred_shuf', 4]
    inter_p_i[rep] = test_coeff[nrow(test_coeff), 4]
  }
  reg_add_p = cbind(reg_add_p, reg_add_p_i)
  inter_p = cbind(inter_p, inter_p_i)
}

reg_add_p = as.data.frame(reg_add_p)
inter_p = as.data.frame(inter_p)
colnames(reg_add_p) = inter_test$regulator
colnames(inter_p) = inter_test$regulator
p_save = list(inter_p = inter_p, reg_add_p = reg_add_p)
saveRDS(p_save, paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/sig_perm/', 
                       target_prot, '.rds'))


