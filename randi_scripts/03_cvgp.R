library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
chr = as.numeric(args[1])

gene_pos_dir='ukb_ppp/ukb_ppp_gene_pos.txt'
genotype_dir='ukb_ppp/genotype/'
rdata_dir='ukb_ppp/ukb_ppp_prot_expr.rdata'
cov_dir='ukb_ppp/covar.txt'
script_dir='cis_trans_inter/99_scripts/gbat/scripts/'
plink_dir = 'software/'
sample_dir='ukb_ppp/genotype/ukb_ppp_chr1_QC.fam'
gcta_dir = 'software/gcta/'

load(rdata_dir)
pos_gene=fread(gene_pos_dir)
allmap=fread(paste(genotype_dir,"ukb_ppp_chr",chr,"_QC.bim",sep=""))
covar=fread(cov_dir)
all_sample = fread(sample_dir)

## keep ind with both covariates and genotypes
common_id = intersect(covar$ID, all_sample$V1)
covar = covar[match(common_id, covar$ID), ]

chr_flag=which(pos_gene$chr==chr)
print(paste0('total number of gene: ' , length(chr_flag)))

for(flag_i in chr_flag){
  prot_i = pos_gene$gene[flag_i]
  exp_gflag=which(gene_name == prot_i)
  if(!( length(exp_gflag)==0)){
    ## extract cis regions
    tss1=max(0,pos_gene$start_hg19[flag_i]-100000)
    tss2=pos_gene$start_hg19[flag_i]+100000
    cis_flag=which(allmap[,4]>tss1 & allmap[,4]<tss2)
    
    if(length(cis_flag)>0){
      phe=prot_ex[,exp_gflag]
      names(phe) = sample_id
      phe = na.omit(phe)
      reml_id = intersect(common_id, names(phe))
      
      ## make phenotype data. Subset the first 10000 inds to reduce computational burden
      phe = phe[reml_id]
      if (length(phe) > 10000) {
        phe = phe[1:10000]
      }
      fwrite(data.frame(FID = names(phe), IID = names(phe), prot = phe), 
             paste0('cis_trans_inter/01_obs_phe/', prot_i, '.txt'), sep = '\t', col.names = F)
      ## make cis region plink file
      fwrite(allmap[cis_flag,c(1,4,4,4)], 
             paste0("cis_trans_inter/03_geno_grm/chr", chr, ".txt"), col.names=F, sep = '\t')
      system(paste0(plink_dir, "plink --bfile ",genotype_dir,"ukb_ppp_chr",chr,
                   "_QC --extract cis_trans_inter/03_geno_grm/chr",
                   chr,".txt --keep cis_trans_inter/01_obs_phe/", prot_i, 
                   ".txt --make-bed --range --thread-num 16 --out cis_trans_inter/03_geno_grm/chr", chr))
      ## calculate grm using gcta
      system(paste(gcta_dir, "gcta64 --bfile cis_trans_inter/03_geno_grm/chr", chr,
                   " --make-grm --out cis_trans_inter/03_geno_grm/chr", chr, sep = ""))
      ## calculate cvgp and h2
      system(paste0(gcta_dir, "gcta64 --reml --cvblup --thread-num 16 --grm cis_trans_inter/03_geno_grm/chr", chr,
                    " --pheno cis_trans_inter/01_obs_phe/", prot_i, ".txt", 
                    " --qcovar cis_trans_inter/00_ref/gcta_cov.qcovar --out cis_trans_inter/02_cvgp/", 
                    prot_i))
    }
  }
  print(which(chr_flag == flag_i))
}





