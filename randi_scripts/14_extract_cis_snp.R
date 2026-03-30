setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(plink2R)

genotype_dir = 'ukb_ppp/genotype/plink_geno/'
plink_dir = 'software/'
gene_pos_dir='ukb_ppp/ukb_ppp_gene_pos.txt'
load('ukb_ppp/ukb_ppp_prot_expr.rdata')
colnames(prot_ex) = gene_name
pos_gene = fread(gene_pos_dir)

sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter_w_cis.txt')
genes = unique(c(sig_inter$target, sig_inter$regulator))
for (i in 1:length(genes)) {
  gene_i = genes[i]
  chr_i = pos_gene$chr[which(pos_gene$gene == gene_i)]
  
  load(paste0('cis_trans_inter/05_fusion/output/', 
              gene_i, '.wgt.RDat'))
  weight_i = wgt.matrix
  weight_i = weight_i[weight_i[,1] != 0, ]
  
  fwrite(data.frame(FID = sample_id, IID = sample_id), 
         paste0('cis_trans_inter/06_fusion_inter/tmp/', gene_i, '.txt'), 
         sep = '\t', col.names = F)
  fwrite(data.frame(snp = rownames(weight_i)), 
         paste0("cis_trans_inter/06_fusion_inter/tmp/", gene_i, "_qtl.txt"), 
         sep = '\t', col.names = F)
  system(paste0(plink_dir, "plink --bfile ", genotype_dir, "chr", chr_i,
                "_QC --extract cis_trans_inter/06_fusion_inter/tmp/",
                gene_i, "_qtl.txt --keep cis_trans_inter/06_fusion_inter/tmp/", 
                gene_i, ".txt --make-bed --thread-num 16 --out cis_trans_inter/06_fusion_inter/tmp/", 
                gene_i), intern=T)
}
