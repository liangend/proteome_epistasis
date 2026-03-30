library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')

for (i in 1:22) {
  cov_i = fread(paste0('ukb_ppp/expr_cov_by_chr/cov_chr', i, '.txt'), header = T)
  cov_i = as.data.frame(cov_i)
  rownames(cov_i) = cov_i$var
  cov_i = t(cov_i[, -1])
  cov_i = cbind(FID = rownames(cov_i), IID = rownames(cov_i),
                cov_i)
  fwrite(cov_i, paste0('cis_trans_inter/05_fusion/cov_fusion/chr', i, '.txt'), sep = '\t',
         quote = F)
}

