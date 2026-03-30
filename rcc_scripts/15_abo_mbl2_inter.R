setwd('/project/xuanyao/jinghui')
library(data.table)
library(sandwich)
args = commandArgs(trailingOnly = T)

load('cis_trans_inter/data/ukb_ppp_prot_expr.rdata')
ukb_cov = fread('cis_trans_inter/00_ref/gcta_cov_full.qcovar')

ABO_pred = fread('cis_trans_inter/06_fusion_pred/eur/ABO_good_pred.txt')
MBL2_pred = fread('cis_trans_inter/06_fusion_pred/eur/MBL2_good_pred.txt')

ABO_id = as.character(ABO_pred$id)
MBL2_id = as.character(MBL2_pred$id)

## 29 batches, each with 100 proteins
batch_num = as.numeric(args[1])
gene_num = ((batch_num-1)*100+1):min(batch_num*100, length(gene_name))

inter_result = c()
for (i in gene_num) {
  prot_i = prot_ex[, i]
  names(prot_i) = sample_id
  prot_i = na.omit(prot_i)
  common_id = Reduce(intersect, list(ABO_id, MBL2_id, 
                                     as.character(ukb_cov$V1), 
                                     names(prot_i)))
  prot_i = prot_i[common_id]
  ABO_i = ABO_pred$pred[match(common_id, as.character(ABO_pred$id))]
  MBL2_i = MBL2_pred$pred[match(common_id, as.character(MBL2_pred$id))]
  cov_i = ukb_cov[match(common_id, as.character(ukb_cov$V1)), -(1:2)]
  res_i = residuals(lm(prot_i ~ as.matrix(cov_i)))
  
  ## regression to estimate interaction effect
  test_i = lm(res_i ~ ABO_i + MBL2_i + ABO_i * MBL2_i)
  coeff_i = summary(test_i)$coefficients
  
  ## using sandwich to correct for heteroscedasticity
  sandwich_se = diag(vcovHC(test_i, type = "HC"))^0.5
  sandwich_t = coef(summary(test_i))[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  
  result_i = c(coeff_i[-1, 1], coeff_i[-1, 4], sandwich_p[-1])
  names(result_i) = c('ABO_eff', 'MBL2_eff', 'ABO.MBL2_eff', 
                      'ABO_p', 'MBL2_p', 'ABO.MBL2_p',
                      'ABO_sand_p', 'MBL2_sand_p', 'ABO.MBL2_sand_p')
  
  inter_result = rbind(inter_result, result_i)
  print(i)
}

inter_result = data.frame(prot = gene_name[gene_num], inter_result)
fwrite(inter_result, 'cis_trans_inter/14_ABO_MBL2_inter/inter_p.txt', 
       sep = '\t', append = T)








