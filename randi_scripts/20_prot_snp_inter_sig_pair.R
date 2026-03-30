setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(plink2R)
library(sandwich)
args = commandArgs(trailingOnly = T)
inter_i = as.numeric(args[1])
genotype_dir = 'ukb_ppp/genotype/plink_geno/'
plink_dir = 'software/'

# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('ukb_ppp/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter.txt')
target_i = sig_inter$target[inter_i]
reg_i = sig_inter$regulator[inter_i]

# regulator prediction
reg_pred = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
                          reg_i, '_good_pred.txt'))
reg_pred$pred = (reg_pred$pred - mean(reg_pred$pred)) / 
  sd(reg_pred$pred)

# target observation
target_ex = prot_ex[, target_i]
names(target_ex) = sample_id
target_ex = na.omit(target_ex)

# target cis qtls used for prediction
target_qtl = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/', target_i),
                        impute="avg")
  
# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

# extract common id
common_id = Reduce(intersect, list(names(target_ex), reg_pred$id,
                                   target_qtl$fam$V1, st_cov$V1))
target_ex = target_ex[common_id]
reg_pred = reg_pred$pred[match(common_id, reg_pred$id)]
target_snp = target_qtl$bed[match(common_id, target_qtl$fam$V1), ]
reg_cov = st_cov[match(common_id, st_cov$V1), -(1:2)]

target_res_ex = residuals(lm(target_ex ~ as.matrix(reg_cov)))

prot_snp_result = c()
for (j in 1:ncol(target_snp)) {
  snp_j = target_snp[, j]
  snp_j = (snp_j - mean(snp_j, na.rm = T)) / sd(snp_j, na.rm = T)
  
  test_j = lm(target_res_ex ~ reg_pred + snp_j + reg_pred * snp_j)
  coeff_j = summary(test_j)$coefficients
  
  ## using sandwich to correct for heteroscedasticity
  sandwich_se = diag(vcovHC(test_j, type = "HC"))^0.5
  sandwich_t = coef(summary(test_j))[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  
  result_j = c(coeff_j[-1, 1], coeff_j[-1, 4], sandwich_p[-1])
  names(result_j) = c('reg_eff', 'snp_eff', 'reg.snp_eff', 
                      'reg_p', 'snp_p', 'reg.snp_p',
                      'reg_sand_p', 'snp_sand_p', 'reg.snp_sand_p')

  prot_snp_result = rbind(prot_snp_result, result_j)
}

prot_snp_result = data.frame(snp = colnames(target_snp), prot_snp_result)
fwrite(prot_snp_result, paste0('cis_trans_inter/11_prot_snp_inter/', 
                               target_i, '_', reg_i, '.txt'), sep = '\t')



