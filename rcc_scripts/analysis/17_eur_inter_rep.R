setwd('/project/xuanyao/jinghui')
library(data.table)
library(sandwich)
library(qvalue)
# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('cis_trans_inter/00_ref/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

pos_gene = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')

pred_all = list.files('cis_trans_inter/06_fusion_pred/afr_pred_expr/')
prot_all = sapply(strsplit(pred_all, '_', fixed = T), '[', 1)
pos_gene = pos_gene[pos_gene$gene %in% prot_all, ]

prot_ex = prot_ex[, prot_all]

# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov_full.qcovar')

# significant interactions need to be tested
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter_w_cis.txt')
inter_p = c()
inter_p_sand = c()

for(i in 1:nrow(sig_inter)){
  target_i = sig_inter$target[i]
  reg_i = sig_inter$regulator[i]
  target_file = list.files('cis_trans_inter/06_fusion_pred/afr_pred_expr/', 
                           pattern = target_i)
  target_pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/afr_pred_expr/', 
                               target_file))
  reg_file = list.files('cis_trans_inter/06_fusion_pred/afr_pred_expr/', 
                        pattern = reg_i)
  reg_pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/afr_pred_expr/', 
                            reg_file))
  target_obs = prot_ex[, target_i]
  names(target_obs) = sample_id
      
  # use common individuals between cis and trans prot
  common_id = Reduce(intersect, list(names(target_obs), reg_pred_i$id, 
                                     target_pred_i$id))
  target_obs = target_obs[as.character(common_id)]
  target_pred = target_pred_i$pred[match(common_id, target_pred_i$id)]
  reg_pred = reg_pred_i$pred[match(common_id, reg_pred_i$id)]
  
  st_cov_i = st_cov[match(common_id, st_cov$V1), ]
  st_cov_i = as.matrix(st_cov_i[,c(3, 130:150)])
  
  ## OLS of linear model
  test=lm(target_obs ~ target_pred + reg_pred + target_pred * reg_pred + st_cov_i)
  ## using sandwich to correct for heteroscedasticity
  sandwich_se = diag(vcovHC(test, type = "HC"))^0.5
  sandwich_t = coef(summary(test))[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  
  test_coeff = summary(test)$coefficients
  if ("target_pred:reg_pred" %in% rownames(test_coeff)) {
    inter_p[i] = test_coeff["target_pred:reg_pred", 4]
  } else {
    inter_p[i] = NA
  }
  inter_p_sand[i] = sandwich_p['target_pred:reg_pred']
}

sig_inter$afr_p = inter_p
sig_inter$afr_sand_p = inter_p_sand
1 - qvalue(inter_p_sand)$pi0

## If significant interactions are in PPI
ppi_all = readRDS('pqtl/11_prot_complex/ppi_list.rds')
sig_reg = c('ABO','ABO','APOH','ABO',
            'ABO','ABO','OBP2B')
sig_target = c('MBL2','CD209','NPTXR',
               'ICAM5','NPTX1','FCGR2B','MBL2')
sig_pair = paste(sig_reg, sig_target, sep = ',')
sum(sig_pair %in% ppi_all)


## pval of pQTLs of ABO on MBL2
ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
abo_qtl = ukb_qtl[ukb_qtl$target_gene == 'ABO', ]

mbl2_file = list.files('/project2/xuanyao/jinghui/UKB_PPP_combined/MBL2_P11226_OID30759_v1_Inflammation_II/')
log10_p = c()
for (i in 1:nrow(abo_qtl)) {
  chr_i = abo_qtl$CHROM[i]
  bp_i = abo_qtl$bp_hg38[i]
  file_i = mbl2_file[grep(paste0('chr', chr_i, '_'), mbl2_file)]
  summ_stat_i = fread(paste0('/project2/xuanyao/jinghui/UKB_PPP_combined/MBL2_P11226_OID30759_v1_Inflammation_II/',
                             file_i))
  log10_p[i] = summ_stat_i$LOG10P[which(summ_stat_i$GENPOS == bp_i)]
  print(i)
}
abo_qtl = abo_qtl[,c(1:2, 14)]
abo_qtl$log10p_mbl2 = log10_p



