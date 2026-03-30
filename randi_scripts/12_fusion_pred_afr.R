library(data.table)
library(plink2R)
setwd('/gpfs/data/xliu-lab/jinghui')

gene_pos_dir='ukb_ppp/ukb_ppp_gene_pos.txt'
genotype_dir='ukb_ppp/genotype/afr_geno/'
expr_dir='cis_trans_inter/05_fusion/fusion_black_data/'
plink_dir = 'software/plink'
output_dir = 'cis_trans_inter/05_fusion/afr_pred_expr/'

pos_gene = fread(gene_pos_dir)
all_model = list.files('cis_trans_inter/05_fusion/afr_output/')

load('ukb_ppp/ukb_ppp_prot_expr.rdata')

for (i in all_model) {
  prot_i = sub('.wgt.RDat', '', i)
  chr_i = pos_gene$chr[which(pos_gene$gene == prot_i)]
  tss_i = pos_gene$start_hg19[which(pos_gene$gene == prot_i)]
  ## extract ind needed to be predicted
  pred_id = fread(paste0(expr_dir, 'pred_id.txt'))
  fwrite(data.frame(FID = pred_id$id,
                    IID = pred_id$id),
         paste0('cis_trans_inter/05_fusion/pred_id/', prot_i, '.txt'), sep = '\t')

  ## extract cis region genotype
  fwrite(data.frame(chr = chr_i, start = max(0, tss_i-500000),
                    end = tss_i+500000, prot = prot_i),
         paste0("cis_trans_inter/05_fusion/cis_geno/", prot_i, ".txt"), col.names = F, sep = '\t')
  system(paste0(plink_dir, " --bfile ", genotype_dir, "chr", chr_i,
                "_QC --keep cis_trans_inter/05_fusion/pred_id/", prot_i,
                ".txt --extract range cis_trans_inter/05_fusion/cis_geno/", prot_i,
                ".txt --make-bed --out cis_trans_inter/05_fusion/cis_geno/", prot_i))

  ## predict prot expr using fusion weights
  load(paste0('cis_trans_inter/05_fusion/afr_output/', i))
  weight_i = wgt.matrix
  cis_geno = read_plink(paste0('cis_trans_inter/05_fusion/cis_geno/', prot_i),
                        impute="avg")
  if (sum(rownames(weight_i) %in% colnames(cis_geno$bed)) != nrow(weight_i)) {
	next
  }
  cis_geno$bed = cis_geno$bed[, rownames(weight_i)]
  pred_expr = cis_geno$bed %*% weight_i[, 1]

  ## save the prediction if cor between obs and pred > 0.1
  obs_expr = prot_ex[match(cis_geno$fam$V1, sample_id), which(gene_name == prot_i)]
  pred_save = data.frame(id = cis_geno$fam$V1, pred = pred_expr, obs = obs_expr)
  cor_i = cor(obs_expr, pred_expr, use = 'na.or.complete')
  if (!is.na(cor_i) & cor_i > 0.1) {
    fwrite(pred_save, paste0(output_dir, prot_i, '_good_pred.txt'), sep = '\t')
  } else {
    fwrite(pred_save, paste0(output_dir, prot_i, '_bad_pred.txt'), sep = '\t')
  }
}

