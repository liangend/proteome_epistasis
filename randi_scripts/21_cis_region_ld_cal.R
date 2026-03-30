setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(plink2R)

sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter.txt')

for (i in 1:6) {
  target_i = sig_inter$target[i]
  target_qtl = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/', target_i),
                          impute="avg")
  snp_ld = cor(target_qtl$bed)
  saveRDS(snp_ld, paste0('cis_trans_inter/11_prot_snp_inter/', target_i, '_ld.rds'))
}

