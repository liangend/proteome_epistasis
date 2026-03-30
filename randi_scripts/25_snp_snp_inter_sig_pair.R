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

sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter_w_cis.txt')
rm_pair = which((sig_inter$target_chr == 6 & sig_inter$reg_chr == 6 &
                   sig_inter$target_tss > 28477797 &  sig_inter$target_tss < 33448354 &
                   sig_inter$reg_tss > 28477797 & sig_inter$reg_tss < 33448354) | 
                  sig_inter$target %in% c('BTN3A2', 'HLA-A') |
                  sig_inter$regulator == 'BTN3A2')
sig_inter = sig_inter[-rm_pair, ]

target_i = sig_inter$target[inter_i]
reg_i = sig_inter$regulator[inter_i]

# target observation
target_ex = prot_ex[, target_i]
names(target_ex) = sample_id
target_ex = na.omit(target_ex)

# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

# target cis qtls used for prediction
target_qtl = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/', target_i),
                        impute="avg")

# regulator cis qtls used for prediction
reg_qtl = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/', reg_i),
                     impute="avg")

# extract common id
common_id = Reduce(intersect, list(names(target_ex), st_cov$V1,
                                   target_qtl$fam$V1))
target_ex = target_ex[common_id]
reg_cov = st_cov[match(common_id, st_cov$V1), -(1:2)]

target_snp = target_qtl$bed[match(common_id, target_qtl$fam$V1), ]
reg_snp = reg_qtl$bed[match(common_id, reg_qtl$fam$V1), ]

target_res_ex = residuals(lm(target_ex ~ as.matrix(reg_cov)))

snp_snp_inter = list()
snp_snp_inter_wo_add = list()
for (i in 1:ncol(reg_snp)) {
  reg_snp_i = reg_snp[, i]
  reg_snp_i = (reg_snp_i - mean(reg_snp_i, na.rm = T)) / 
    sd(reg_snp_i, na.rm = T)
  snp_i = colnames(reg_snp)[i]
  for (j in 1:ncol(target_snp)) {
    target_snp_j = target_snp[, j]
    target_snp_j = (target_snp_j - mean(target_snp_j, na.rm = T)) / 
      sd(target_snp_j, na.rm = T)
    snp_j = colnames(target_snp)[j]
    
    test_ij = lm(target_res_ex ~ reg_snp_i + target_snp_j + 
                   reg_snp_i * target_snp_j)
    coef_ij = summary(test_ij)$coefficients
    
    ## using sandwich to correct for heteroscedasticity
    sandwich_se = diag(vcovHC(test_ij, type = "HC"))^0.5
    sandwich_t = coef(summary(test_ij))[, 1]/sandwich_se
    sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
    
    coef_ij = cbind(coef_ij, sandwich_p)
    
    ## test using interaction only
    inter_snp = reg_snp_i * target_snp_j
    test_ij2 = lm(target_res_ex ~ inter_snp)
    coef_ij2 = summary(test_ij2)$coefficients
    
    snp_snp_inter[[paste0(snp_i,'*',snp_j)]] = coef_ij
    snp_snp_inter_wo_add[[paste0(snp_i,'*',snp_j)]] = coef_ij2
  }
}

saveRDS(snp_snp_inter, paste0('cis_trans_inter/14_snp_snp_inter/inter_w_add/', 
                              target_i, '_', reg_i, '.rds'))
saveRDS(snp_snp_inter_wo_add, paste0('cis_trans_inter/14_snp_snp_inter/inter_wo_add/', 
                                     target_i, '_', reg_i, '.rds'))



