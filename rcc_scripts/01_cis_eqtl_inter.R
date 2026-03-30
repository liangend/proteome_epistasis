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

dgn_cis = fread('cis_trans_inter/data/dgn_cis_snp.traw')
# check if sample IDs are matched
# samples = fread('cis_trans_inter/data/samples.txt')
# cis_sample = colnames(dgn_cis)[-(1:6)]
# cis_sample = sapply(strsplit(cis_sample, '_', fixed = T), '[', 1)
# sum(cis_sample == samples$sample)
dgn_cis_qtl = fread('gtex/06_qtl_z/dgn_cis/all.cis_gene.txt.gz')
gene_meta = fread('gtex/00_ref/gene_meta_dgn.txt')
gene_chr1 = gene_meta$gene_id[which(gene_meta$chr == 'chr1')]
ex = ex[, gene_chr1]

trans_p = c()
inter_p = c()
gene_cal = c()
for(i in 1:ncol(g)){
  if(!is.na(g[1,i]) & (pr[i,2]>0.1)){
    tryCatch({
      g_sva=make_sva(ex,100,g[,i])
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
        eqtl_j = dgn_cis_qtl$variant_id[which(dgn_cis_qtl$phenotype_id == colnames(rQN)[j])]
        if (length(eqtl_j) > 0) {
          snp_j = unlist(dgn_cis[which(dgn_cis$SNP == eqtl_j), -(1:6)])
          test=lm(gene_j ~ g_norm + snp_j + g_norm * snp_j)
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
fwrite(trans_p, paste0('cis_trans_inter/output/trans_p/chr', chr_index, 
                       '_sub', sub_index, '.txt'), sep = '\t', row.names = T)
fwrite(inter_p, paste0('cis_trans_inter/output/inter_p/chr', chr_index, 
                       '_sub', sub_index, '.txt'), sep = '\t', row.names = T)

