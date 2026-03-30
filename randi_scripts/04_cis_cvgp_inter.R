args = commandArgs(trailingOnly = T)
setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(phenix)
source('cis_trans_inter/99_scripts/gbat/scripts/make_smartsva.R')

load('ukb_ppp/ukb_ppp_prot_expr.rdata') # expression matrix prot_ex; protein name gene_name; and sample id smaple_id
colnames(prot_ex) = gene_name

chr_index = as.numeric(args[1])

pos_gene = fread('ukb_ppp/ukb_ppp_gene_pos.txt')
cvgp_all = list.files('cis_trans_inter/02_cvgp', pattern = '.cvblp')
cvgp_all = sub('.indi.cvblp', '', cvgp_all)
# genes having cv cis prediction
pos_gene = pos_gene[pos_gene$gene %in% cvgp_all, ]

## genes on chr1
gene_chr1 = pos_gene$gene[which(pos_gene$chr == 1)]
prot_ex = prot_ex[, gene_chr1]

# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

## genes on chr_i
gene_chr_i = pos_gene$gene[which(pos_gene$chr == chr_index)]

trans_p = c()
inter_p = c()
gene_cal = c()
for(i in 1:length(gene_chr_i)){
  gene_i = gene_chr_i[i]
  cvpg_i = fread(paste0('cis_trans_inter/02_cvgp/', gene_i, '.indi.cvblp'))
  obs_i = fread(paste0('cis_trans_inter/01_obs_phe/', gene_i, '.txt'))
  cvpg_i = cvpg_i$V4[match(obs_i$V1, cvpg_i$V1)]
  names(cvpg_i) = obs_i$V1
  pred_cor = cor(cvpg_i, obs_i$V3)
  if(pred_cor > 0.1){
    cvpg_norm_i = (cvpg_i - mean(cvpg_i))/sd(cvpg_i)

    trans_p_i = c()
    inter_p_i = c()
    # loop over every protein on chr1
    for (j in 1:length(gene_chr1)) {
      gene_j = gene_chr1[j]
      prot_ex_j = prot_ex[, gene_j]
      names(prot_ex_j) = sample_id
      prot_ex_j = na.omit(prot_ex_j)
      # cv cis prediction of prot j
      cvpg_j = fread(paste0('cis_trans_inter/02_cvgp/', gene_j, '.indi.cvblp'))
      cvpg_j$V1 = as.character(cvpg_j$V1)
      common_id = Reduce(intersect, list(names(prot_ex_j), cvpg_j$V1, names(cvpg_norm_i)))
      prot_ex_j = prot_ex_j[common_id]
      cvpg_norm_i = cvpg_norm_i[common_id]
      cvpg_j = cvpg_j$V4[match(common_id, cvpg_j$V1)]
      cvgp_norm_j = (cvpg_j - mean(cvpg_j))/sd(cvpg_j)
      st_cov_ij = st_cov[match(common_id, st_cov$V1), ]
      st_cov_ij = as.matrix(st_cov_ij[,-(1:2)])
      
      test=lm(prot_ex_j ~ cvpg_norm_i + cvgp_norm_j + cvpg_norm_i * cvgp_norm_j + st_cov_ij)
      test_coeff = summary(test)$coefficients
      trans_p_i = c(trans_p_i, test_coeff[2,4])
      inter_p_i = c(inter_p_i, test_coeff[nrow(test_coeff),4])
      names(trans_p_i)[length(trans_p_i)] = gene_j
      names(inter_p_i)[length(inter_p_i)] = gene_j
    }
    trans_p = rbind(trans_p, trans_p_i)
    inter_p = rbind(inter_p, inter_p_i)
    gene_cal = c(gene_cal, gene_i)
  }
}
trans_p = as.data.frame(trans_p)
inter_p = as.data.frame(inter_p)
rownames(trans_p) = gene_cal
rownames(inter_p) = gene_cal
fwrite(trans_p, paste0('cis_trans_inter/04_cis_cvgp_inter/trans_p/chr', chr_index, 
                       '.txt'), sep = '\t', row.names = T)
fwrite(inter_p, paste0('cis_trans_inter/04_cis_cvgp_inter/inter_p/chr', chr_index, 
                       '.txt'), sep = '\t', row.names = T)
