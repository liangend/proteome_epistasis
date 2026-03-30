library(data.table)
library(plink2R)
setwd('/gpfs/data/xliu-lab/jinghui')

genotype_dir='cis_trans_inter/09_inter_pheno_asso/all_white_geno/'
plink_dir = 'software/plink'
output_dir = 'cis_trans_inter/16_prot_trans_pqtl_inter/pred_prot/'

fusion_pos = fread('cis_trans_inter/00_ref/ukb_fusion.pos')
sig_inter = fread('cis_trans_inter/16_prot_trans_pqtl_inter/inter_all.txt')
sig_inter = sig_inter[sig_inter$qval < 0.05, ]
pred_prot = fusion_pos[fusion_pos$ID %in% sig_inter$prot, ]
# pred_prot = fread('cis_trans_inter/09_inter_pheno_asso/prot_in_sig_inter.txt')

for (i in 1:nrow(pred_prot)) {
  prot_i = pred_prot$ID[i]
  chr_i = pred_prot$CHR[i]
  tss_i = pred_prot$P0[i]
  model_i = pred_prot$WGT[i]
  
  ## extract cis region genotype
  fwrite(data.frame(chr = chr_i, start = max(0, tss_i-500000),
                    end = tss_i+500000, prot = prot_i),
         paste0("cis_trans_inter/09_inter_pheno_asso/cis_geno/", prot_i, ".txt"), 
         col.names = F, sep = '\t')
  system(paste0(plink_dir, " --bfile ", genotype_dir, "chr", chr_i,
                "_QC --extract range cis_trans_inter/09_inter_pheno_asso/cis_geno/", prot_i,
                ".txt --make-bed --out cis_trans_inter/09_inter_pheno_asso/cis_geno/", prot_i))
  
  ## predict prot expr using fusion weights
  load(model_i)
  weight_i = wgt.matrix
  weight_i = weight_i[weight_i[,1] != 0, ]
  cis_geno = read_plink(paste0('cis_trans_inter/09_inter_pheno_asso/cis_geno/', prot_i),
                        impute="avg")
  shared_snp = intersect(rownames(weight_i), colnames(cis_geno$bed))
  pred_expr = cis_geno$bed[, shared_snp] %*% weight_i[shared_snp, 1]
  
  ## save the prediction
  pred_save = data.frame(id = cis_geno$fam$V1, pred = pred_expr)
  fwrite(pred_save, paste0(output_dir, prot_i, '_pred.txt'), sep = '\t')
  print(i)
}



