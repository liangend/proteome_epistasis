setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(plink2R)
args = commandArgs(trailingOnly = T)
genotype_dir = 'ukb_ppp/genotype/plink_geno/'
plink_dir = 'software/'

# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('ukb_ppp/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

cis_pqtl = fread('cis_trans_inter/00_ref/ukb_cis_pqtl_in_geno.txt')

cis_chr = as.numeric(args[1])
trans_chr = as.numeric(args[2])
print(paste0('cis chr', cis_chr, ' trans chr', trans_chr))

if (cis_chr == trans_chr) {
  stop("same chr")
} else {
  pos_gene = fread('ukb_ppp/ukb_ppp_gene_pos.txt')
  
  pred_all = list.files('cis_trans_inter/05_fusion/pred_expr', pattern = '_good_pred.txt')
  pred_all = sub('_good_pred.txt', '', pred_all)
  # genes having good cis prediction (cor > 0.1)
  pos_gene = pos_gene[pos_gene$gene %in% pred_all, ]
  pos_gene = pos_gene[pos_gene$gene %in% cis_pqtl$prot, ]
  
  ## genes on cis chr
  gene_cis = pos_gene$gene[which(pos_gene$chr == cis_chr)]
  prot_ex = prot_ex[, gene_cis]
  
  # standardized covarites
  st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')
  
  ## genes on trans chr
  gene_trans = pos_gene$gene[which(pos_gene$chr == trans_chr)]
  
  inter_p_w_cis_prot = c()
  inter_p_wo_cis_prot = c()
  gene_cal = c()
  for(i in 1:length(gene_trans)){
    trans_gene_i = gene_trans[i]
    trans_pred_i = fread(paste0('cis_trans_inter/05_fusion/pred_expr/', 
                                trans_gene_i, '_good_pred.txt'))
    trans_pred_i$pred = (trans_pred_i$pred - mean(trans_pred_i$pred))/
      sd(trans_pred_i$pred)
    
    inter_p_w_cis_prot_i = c()
    inter_p_wo_cis_prot_i = c()
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
      
      cis_qtl_j = cis_pqtl[cis_pqtl$prot == cis_gene_j, ]
      cis_qtl_j = cis_qtl_j[1, ] ## one prot has two cis-pqtls, and I only use the first one here
      
      common_id = Reduce(intersect, list(names(cis_ex_j), cis_pred_j$id, trans_pred_i$id))
      
      ## extrat cis-pqtl genotype from plink file
      fwrite(data.frame(FID = common_id, IID = common_id), 
             paste0('cis_trans_inter/06_fusion_inter/tmp/cis_', cis_chr, '_trans_', trans_chr, '_id.txt'), 
             sep = '\t', col.names = F)
      fwrite(data.frame(chr = cis_qtl_j$chr, bp1 = cis_qtl_j$bp_hg19,
                        bp2 = cis_qtl_j$bp_hg19, bp3 = cis_qtl_j$bp_hg19), 
             paste0("cis_trans_inter/06_fusion_inter/tmp/cis_", cis_chr, '_trans_', trans_chr, "_qtl.txt"), 
             col.names=F, sep = '\t')
      system(paste0(plink_dir, "plink --bfile ", genotype_dir, "chr", cis_chr,
                    "_QC --extract cis_trans_inter/06_fusion_inter/tmp/cis_",
                    cis_chr, '_trans_', trans_chr, "_qtl.txt --keep cis_trans_inter/06_fusion_inter/tmp/cis_", 
		    cis_chr, '_trans_', trans_chr, 
                    "_id.txt --make-bed --range --thread-num 16 --out cis_trans_inter/06_fusion_inter/tmp/cis_", 
                    cis_chr, '_trans_', trans_chr), intern=T)
      cis_qtl = read_plink(paste0('cis_trans_inter/06_fusion_inter/tmp/cis_', cis_chr, '_trans_', trans_chr),
                           impute="avg")
      
      cis_ex_j = cis_ex_j[common_id]
      trans_prot_i = trans_pred_i$pred[match(common_id, trans_pred_i$id)]
      cis_prot_j = cis_pred_j$pred[match(common_id, cis_pred_j$id)]
      
      st_cov_ij = st_cov[match(common_id, st_cov$V1), ]
      st_cov_ij = as.matrix(st_cov_ij[,-(1:2)])
      
      cis_qtl = cis_qtl$bed[match(common_id, cis_qtl$fam$V1), 1]
      cis_qtl = (cis_qtl - mean(cis_qtl, na.rm = T))/sd(cis_qtl, na.rm = T)
      
      test1=lm(cis_ex_j ~ st_cov_ij + trans_prot_i + cis_qtl + cis_qtl * trans_prot_i)
      test2=lm(cis_ex_j ~ st_cov_ij + cis_prot_j + trans_prot_i + cis_qtl + cis_qtl * trans_prot_i)
      test_coeff1 = summary(test1)$coefficients
      test_coeff2 = summary(test2)$coefficients
      inter_p_wo_cis_prot_i = c(inter_p_wo_cis_prot_i, test_coeff1[nrow(test_coeff1), 4])
      inter_p_w_cis_prot_i = c(inter_p_w_cis_prot_i, test_coeff2[nrow(test_coeff2), 4])

      names(inter_p_wo_cis_prot_i)[length(inter_p_wo_cis_prot_i)] = cis_gene_j
      names(inter_p_w_cis_prot_i)[length(inter_p_w_cis_prot_i)] = cis_gene_j
    }
    inter_p_wo_cis_prot = rbind(inter_p_wo_cis_prot, inter_p_wo_cis_prot_i)
    inter_p_w_cis_prot = rbind(inter_p_w_cis_prot, inter_p_w_cis_prot_i)
    gene_cal = c(gene_cal, trans_gene_i)
    print(paste0('cis j ',j, ', trans i ', i))
  }
  inter_p_wo_cis_prot = as.data.frame(inter_p_wo_cis_prot)
  inter_p_w_cis_prot = as.data.frame(inter_p_w_cis_prot)
  rownames(inter_p_wo_cis_prot) = gene_cal
  rownames(inter_p_w_cis_prot) = gene_cal
  p_save = list(inter_p_wo_cis_prot = inter_p_wo_cis_prot, 
                inter_p_w_cis_prot = inter_p_w_cis_prot)
  saveRDS(p_save, paste0('cis_trans_inter/06_fusion_inter/prot_snp_inter/cis_chr', 
                         cis_chr, '_trans_chr', trans_chr, '.rds'))
}


