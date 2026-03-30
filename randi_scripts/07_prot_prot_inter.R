args = commandArgs(trailingOnly = T)
setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)

# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('ukb_ppp/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

cis_chr = as.numeric(args[1])
trans_chr = as.numeric(args[2])
if (cis_chr == trans_chr) {
  stop("same chr")
} else {
  pos_gene = fread('ukb_ppp/ukb_ppp_gene_pos.txt')
  
  pred_all = list.files('cis_trans_inter/05_fusion/pred_expr', pattern = '_good_pred.txt')
  pred_all = sub('_good_pred.txt', '', pred_all)
  # genes having good cis prediction (cor > 0.1)
  pos_gene = pos_gene[pos_gene$gene %in% pred_all, ]
  
  ## genes on cis chr
  gene_cis = pos_gene$gene[which(pos_gene$chr == cis_chr)]
  prot_ex = prot_ex[, gene_cis]
  
  # standardized covarites
  st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')
  
  ## genes on trans chr
  gene_trans = pos_gene$gene[which(pos_gene$chr == trans_chr)]
  
  trans_p = c()
  inter_p = c()
  gene_cal = c()
  for(i in 1:length(gene_trans)){
    trans_gene_i = gene_trans[i]
    trans_pred_i = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
                                trans_gene_i, '_good_pred.txt'))
    trans_pred_i$pred = (trans_pred_i$pred - mean(trans_pred_i$pred))/
      sd(trans_pred_i$pred)
    
    trans_p_i = c()
    inter_p_i = c()
    # loop over every protein on cis chr
    for (j in 1:length(gene_cis)) {
      cis_gene_j = gene_cis[j]
      cis_ex_j = prot_ex[, cis_gene_j]
      names(cis_ex_j) = sample_id
      cis_ex_j = na.omit(cis_ex_j)
      # cis prediction of prot j
      cis_pred_j = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
                                cis_gene_j, '_good_pred.txt'))
      cis_pred_j$pred = (cis_pred_j$pred - mean(cis_pred_j$pred))/
        sd(cis_pred_j$pred)
      common_id = Reduce(intersect, list(names(cis_ex_j), cis_pred_j$id, trans_pred_i$id))
      cis_ex_j = cis_ex_j[common_id]
      trans_prot_i = trans_pred_i$pred[match(common_id, trans_pred_i$id)]
      cis_prot_j = cis_pred_j$pred[match(common_id, cis_pred_j$id)]
      
      st_cov_ij = st_cov[match(common_id, st_cov$V1), ]
      st_cov_ij = as.matrix(st_cov_ij[,-(1:2)])
      
      test=lm(cis_ex_j ~ cis_prot_j + trans_prot_i + cis_prot_j * trans_prot_i + st_cov_ij)
      test_coeff = summary(test)$coefficients
      trans_p_i = c(trans_p_i, test_coeff['trans_prot_i', 4])
      inter_p_i = c(inter_p_i, test_coeff[nrow(test_coeff), 4])
      names(trans_p_i)[length(trans_p_i)] = cis_gene_j
      names(inter_p_i)[length(inter_p_i)] = cis_gene_j
    }
    trans_p = rbind(trans_p, trans_p_i)
    inter_p = rbind(inter_p, inter_p_i)
    gene_cal = c(gene_cal, trans_gene_i)
    
  }
  trans_p = as.data.frame(trans_p)
  inter_p = as.data.frame(inter_p)
  rownames(trans_p) = gene_cal
  rownames(inter_p) = gene_cal
  p_save = list(inter_p = inter_p, trans_p = trans_p)
  saveRDS(p_save, paste0('cis_trans_inter/06_fusion_inter/prot_prot_inter/cis_chr', 
                         cis_chr, '_trans_chr', trans_chr, '.rds'))
}






