args = commandArgs(trailingOnly = T)
setwd('/project/xuanyao/jinghui')
library(data.table)
library(sandwich)

# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('cis_trans_inter/00_ref/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

cis_chr = as.numeric(args[1])
trans_chr = as.numeric(args[1])
# if (cis_chr == trans_chr) {
#   stop("same chr")
if (cis_chr != trans_chr) {
  stop("not same chr")
} else {
  pos_gene = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')
  
  pred_all = list.files('cis_trans_inter/06_fusion_pred/eur/', pattern = '_good_pred.txt')
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
  for(i in 1:length(gene_cis)){
    cis_gene_i = gene_cis[i]
    cis_ex_i = prot_ex[, cis_gene_i]
    names(cis_ex_i) = sample_id
    cis_ex_i = na.omit(cis_ex_i)
    # cis prediction of prot i
    cis_pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', 
                              cis_gene_i, '_good_pred.txt'))
    cis_pred_i$pred = (cis_pred_i$pred - mean(cis_pred_i$pred))/
      sd(cis_pred_i$pred)
    
    # use common individuals between cis and trans prot
    common_id = Reduce(intersect, list(names(cis_ex_i), cis_pred_i$id, st_cov$V1))
    cis_ex_i = cis_ex_i[as.character(common_id)]
    st_cov_ij = st_cov[match(common_id, st_cov$V1), -(1:2)]
    
    cis_ex_res = residuals(lm(cis_ex_i ~ as.matrix(st_cov_ij)))
    
    trans_p_i = c()
    inter_p_i = c()
    # loop over every protein on trans chr
    for (j in 1:length(gene_trans)) {
      trans_gene_j = gene_trans[j]
      trans_pred_j = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', 
                                  trans_gene_j, '_good_pred.txt'))
      trans_pred_j$pred = (trans_pred_j$pred - mean(trans_pred_j$pred))/
        sd(trans_pred_j$pred)
      
      common_id2 = intersect(trans_pred_j$id, common_id)
      
      cis_prot_j = cis_pred_i$pred[match(common_id2, cis_pred_i$id)]
      trans_prot_j = trans_pred_j$pred[match(common_id2, trans_pred_j$id)]
      cis_ex_res_j = cis_ex_res[match(common_id2, common_id)]
      ## OLS of linear model
      test=lm(cis_ex_res_j ~ cis_prot_j + trans_prot_j + cis_prot_j * trans_prot_j)
      ## using sandwich to correct for heteroscedasticity
      sandwich_se = diag(vcovHC(test, type = "HC"))^0.5
      sandwich_t = coef(summary(test))[, 1]/sandwich_se
      sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
      
      # test_coeff = summary(test)$coefficients
      # trans_p_i = c(trans_p_i, test_coeff['trans_prot_i', 4])
      # inter_p_i = c(inter_p_i, test_coeff[nrow(test_coeff), 4])
      trans_p_i = c(trans_p_i, sandwich_p['trans_prot_j'])
      inter_p_i = c(inter_p_i, sandwich_p['cis_prot_j:trans_prot_j'])
      
      names(trans_p_i)[length(trans_p_i)] = trans_gene_j
      names(inter_p_i)[length(inter_p_i)] = trans_gene_j
    }
    trans_p = cbind(trans_p, trans_p_i)
    inter_p = cbind(inter_p, inter_p_i)
    gene_cal = c(gene_cal, cis_gene_i)
    
  }
  trans_p = as.data.frame(trans_p)
  inter_p = as.data.frame(inter_p)
  colnames(trans_p) = gene_cal
  colnames(inter_p) = gene_cal
  p_save = list(inter_p = inter_p, trans_p = trans_p)
  saveRDS(p_save, paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/sandwich_p/cis_chr', 
                         cis_chr, '_trans_chr', trans_chr, '.rds'))
}


