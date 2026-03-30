args = commandArgs(trailingOnly = T)

setwd('/project/xuanyao/jinghui')
library(data.table)
library(phenix)
source('cis_trans_inter/gbat/scripts/make_smartsva.R')

load('cis_trans_inter/data/dgn_2qt.rdata')
colnames(ex) = gnames

chr_index = as.numeric(args[1])
sub_index = as.numeric(args[2])

g = fread(paste0('cis_trans_inter/data/ge_quantnorm/chr', chr_index, 
                 '_sub', sub_index, '.txt'))
g = as.data.frame(g)
pr = fread(paste0('cis_trans_inter/data/ge_quantnorm/pearsonR_chr', chr_index, 
                  '_sub', sub_index, '.txt'))
pr = as.data.frame(pr)

gene_meta = fread('gtex/00_ref/gene_meta_dgn.txt')
gene_chr1 = gene_meta$gene_id[which(gene_meta$chr == 'chr1')]
ex = ex[, gene_chr1]

## cv cis prediction of genes on chr1
cvgp_chr1 = c()
for (i in 1:10) {
  cvgp_chr1_i = fread(paste0('cis_trans_inter/data/ge_quantnorm/chr1', 
                             '_sub', i, '.txt'))
  cvgp_chr1 = cbind(cvgp_chr1, cvgp_chr1_i)
}
cvgp_chr1 = as.data.frame(cvgp_chr1)

trans_p = c()
inter_p = c()
gene_cal = c()
for(i in 1:ncol(g)){
  if(!is.na(g[1,i]) & (pr[i,2]>0.1)){
    tryCatch({
      g_sva=make_sva(ex,100,g[,i]) # sva conditional on trans gene i
      g_norm=(g[,i]-mean(g[,i]))/sd(g[,i])
      if(!is.na(g_sva[1,1])){
        rlm=lm(ex~g_sva[,1:20])
        rQN=quantnorm(as.matrix(resid(rlm)))
      } else{
        rQN=ex
      }
      colnames(rQN) = colnames(ex)
      trans_p_i = c()
      inter_p_i = c()
      for (j in 1:ncol(rQN)) {
        gene_j = rQN[, j]
        if (colnames(rQN)[j] %in% colnames(cvgp_chr1)) {
          ## cv cis prediction of gene j
          cvgp_j = cvgp_chr1[, colnames(rQN)[j]]
          test=lm(gene_j ~ g_norm + cvgp_j + g_norm * cvgp_j)
          test_coeff = summary(test)$coefficients
          trans_p_i = c(trans_p_i, test_coeff[2,4])
          inter_p_i = c(inter_p_i, test_coeff[4,4])
          names(trans_p_i)[length(trans_p_i)] = colnames(rQN)[j]
          names(inter_p_i)[length(inter_p_i)] = colnames(rQN)[j]
        }
      }
      trans_p = rbind(trans_p, trans_p_i)
      inter_p = rbind(inter_p, inter_p_i)
      gene_cal = c(gene_cal, colnames(g)[i])
    }, error=function(e){})
  }
}
trans_p = as.data.frame(trans_p)
inter_p = as.data.frame(inter_p)
rownames(trans_p) = gene_cal
rownames(inter_p) = gene_cal
fwrite(trans_p, paste0('cis_trans_inter/02_cis_cvgp_inter/trans_p/chr', chr_index, 
                       '_sub', sub_index, '.txt'), sep = '\t', row.names = T)
fwrite(inter_p, paste0('cis_trans_inter/02_cis_cvgp_inter/inter_p/chr', chr_index, 
                       '_sub', sub_index, '.txt'), sep = '\t', row.names = T)

