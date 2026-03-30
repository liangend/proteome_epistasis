library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
batch = as.numeric(args[1])

gene_pos_dir='ukb_ppp/ukb_ppp_gene_pos.txt'
genotype_dir='ukb_ppp/genotype/plink_geno/'
expr_dir='ukb_ppp/expr_cov_by_chr/'
cov_dir='cis_trans_inter/05_fusion/cov_fusion/'
plink_dir = 'software/plink'
fusion_dir = 'software/fusion_twas-master/'
gcta_dir = 'software/gcta/gcta64'
tmp_dir = 'cis_trans_inter/05_fusion/tmp/'
output_dir = 'cis_trans_inter/05_fusion/output/'

pos_gene=fread(gene_pos_dir)
prot_done = list.files('cis_trans_inter/05_fusion/pred_expr/')
prot_done = sapply(strsplit(prot_done, '_', fixed = T), '[', 1)
pos_gene = pos_gene[!(pos_gene$gene %in% prot_done), ]

i = batch
prot_i = pos_gene$gene[i]
chr = pos_gene$chr[i]

expr_chr = fread(paste0(expr_dir,'prot_expr_chr',chr,'.txt'), header=T)

if(length(prot_i)>0){
    expr_chr = expr_chr[which(expr_chr$NAME == prot_i), ]
    ## make phenotype data
    fwrite(data.frame(FID = colnames(expr_chr)[-1], 
                      IID = colnames(expr_chr)[-1], 
                      prot = unlist(expr_chr[1, -1])),
           paste0('cis_trans_inter/05_fusion/obs_phe/', prot_i, '.txt'), sep = '\t')
    
    ## extract cis regions
    tss1=max(0, pos_gene$start_hg19[i]-500000)
    tss2=pos_gene$start_hg19[i]+500000
    ## make cis region plink file
    fwrite(data.frame(chr = chr, start = tss1, end = tss2, prot = prot_i),
           paste0("cis_trans_inter/05_fusion/cis_geno/", prot_i, ".txt"), col.names = F, sep = '\t')
    system(paste0(plink_dir, " --bfile ", genotype_dir, "chr",chr,
                  "_QC --pheno cis_trans_inter/05_fusion/obs_phe/", prot_i, 
                  ".txt --keep cis_trans_inter/05_fusion/obs_phe/", prot_i, 
                  ".txt --extract range cis_trans_inter/05_fusion/cis_geno/", prot_i,
                  ".txt --make-bed --out cis_trans_inter/05_fusion/cis_geno/", prot_i))

    ## calculate cis snp weight using fusion. Set h2 = 0.2 to skip h2 estimation since it takes too long
    system(paste0('mkdir ', tmp_dir, prot_i))
    system(paste0('Rscript ', fusion_dir, 'FUSION.compute_weights.R --bfile cis_trans_inter/05_fusion/cis_geno/',
                  prot_i, ' --tmp ', tmp_dir, prot_i, ' --out ', output_dir, prot_i, ' --covar ', cov_dir, 
                  'chr', chr, '.txt --PATH_plink ', plink_dir, ' --PATH_gcta ', gcta_dir, 
                  ' --hsq_set 0.2 --model enet'))
}
     
