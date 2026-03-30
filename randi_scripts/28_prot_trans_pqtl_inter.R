setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(plink2R)
library(sandwich)
args = commandArgs(trailingOnly = T)
genotype_dir = 'ukb_ppp/genotype/plink_geno/'
plink_dir = 'software/'

# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('ukb_ppp/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

trans_pqtl = fread('cis_trans_inter/00_ref/ukb_trans_pqtl.txt')
good_pred = list.files('cis_trans_inter/05_fusion/pred_expr', pattern = 'good_pred.txt')
good_pred = sub('_good_pred.txt', '', good_pred)
trans_pqtl = trans_pqtl[trans_pqtl$prot %in% good_pred, ]

prot_index = as.numeric(args[1])
uniq_prot = unique(trans_pqtl$prot)
prot_i = uniq_prot[prot_index]
trans_pqtl_i = trans_pqtl[trans_pqtl$prot == prot_i, ]

# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

# FUSION predicted protein
prot_pred_i = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
                           prot_i, '_good_pred.txt'))
prot_pred_i$pred = (prot_pred_i$pred - mean(prot_pred_i$pred))/
  sd(prot_pred_i$pred)
prot_pred_i = na.omit(prot_pred_i)

# match up ids between covariates and protein
common_id = intersect(st_cov$V1, prot_pred_i$id)
st_cov = st_cov[match(common_id, st_cov$V1), ]
st_cov = as.matrix(st_cov[,-(1:2)])

prot_pred_i = prot_pred_i[match(common_id, prot_pred_i$id), ]
prot_pred_i$prot_res = residuals(lm(prot_pred_i$obs ~ st_cov))

# save the id file for genotype extraction
fwrite(data.frame(FID = common_id, IID = common_id), 
       paste0('cis_trans_inter/06_fusion_inter/tmp/', prot_i, '_id.txt'), 
       sep = '\t', col.names = F)

inter_p = c()
inter_p_wo_cis_prot = c()
# loop over every trans-pqtl
for (j in 1:nrow(trans_pqtl_i)) {
  trans_qtl_j = trans_pqtl_i[j, ] 
  ## extrat trans-pqtl genotype from plink file
  fwrite(data.frame(chr = trans_qtl_j$chr, bp1 = trans_qtl_j$bp_hg19,
                    bp2 = trans_qtl_j$bp_hg19, bp3 = trans_qtl_j$bp_hg19), 
         paste0("cis_trans_inter/06_fusion_inter/tmp/", prot_i, '_', j, ".txt"), 
         col.names = F, sep = '\t')
  system(paste0(plink_dir, "plink --bfile ", genotype_dir, "chr", trans_qtl_j$chr,
                "_QC --extract cis_trans_inter/06_fusion_inter/tmp/",
                prot_i, '_', j, ".txt --keep cis_trans_inter/06_fusion_inter/tmp/", prot_i, 
                "_id.txt --make-bed --range --thread-num 16 --out cis_trans_inter/06_fusion_inter/tmp/", 
                prot_i, '_', j), intern=T)
  if (!file.exists(paste0("cis_trans_inter/06_fusion_inter/tmp/", prot_i, '_', j, ".bed"))) {
    next
  }
  trans_qtl = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/', prot_i, '_', j),
                         impute="avg")

  trans_qtl = trans_qtl$bed[match(common_id, trans_qtl$fam$V1), 1]
  trans_qtl = (trans_qtl - mean(trans_qtl, na.rm = T))/sd(trans_qtl, na.rm = T)
  
  test = lm(prot_pred_i$prot_res ~ prot_pred_i$pred * trans_qtl)
  test_coeff = summary(test)$coefficients
  
  sandwich_se = diag(vcovHC(test, type = "HC"))^0.5
  sandwich_t = test_coeff[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  test_coeff = cbind(test_coeff, sandwich_t = sandwich_t, sandwich_p = sandwich_p)
  
  trans_qtl_j$inter_p = test_coeff[nrow(test_coeff), 4]
  trans_qtl_j$inter_p_sand = test_coeff[nrow(test_coeff), 6]
  fwrite(trans_qtl_j, 
         paste0('cis_trans_inter/16_prot_trans_pqtl_inter/', prot_i, '.txt'), 
         sep = '\t', append = T)
}





