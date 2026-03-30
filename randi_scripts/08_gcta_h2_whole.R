library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
batch = as.numeric(args[1])

gene_pos_dir='ukb_ppp/ukb_ppp_gene_pos.txt'
genotype_dir='ukb_ppp/genotype/plink_geno/'
rdata_dir='ukb_ppp/ukb_ppp_prot_expr.rdata'
cov_dir='ukb_ppp/covar.txt'
script_dir='cis_trans_inter/99_scripts/gbat/scripts/'
plink_dir = 'software/'
sample_dir='ukb_ppp/genotype/plink_geno/chr1_QC.fam'
gcta_dir = 'software/gcta/'

load(rdata_dir)
pos_gene=fread(gene_pos_dir)
#gene_done = list.files('cis_trans_inter/07_cis_h2', pattern = '.hsq')
#gene_done = sub('.hsq', '', gene_done)
#pos_gene = pos_gene[!(pos_gene$gene %in% gene_done), ]

covar=fread(cov_dir)
all_sample = fread(sample_dir)

## keep ind with both covariates and genotypes
common_id = intersect(covar$ID, all_sample$V1)
covar = covar[match(common_id, covar$ID), ]

pos_gene = pos_gene[((batch-1)*100+1):min(batch*100, nrow(pos_gene))]

for(flag_i in 1:nrow(pos_gene)){
  chr = pos_gene$chr[flag_i]
  prot_i = pos_gene$gene[flag_i]
  
  allmap=fread(paste(genotype_dir,"chr",chr,"_QC.bim",sep=""))
  exp_gflag=which(gene_name == prot_i)
  if(!(length(exp_gflag)==0)){
      tryCatch({
        phe=prot_ex[,exp_gflag]
        names(phe) = sample_id
        phe = na.omit(phe)
        reml_id = intersect(common_id, names(phe))
        
        ## make phenotype data. Subset the first 1000 inds to reduce computational burden
        phe = phe[reml_id]
        if (length(phe) > 1000) {
          phe = phe[sample(1:length(phe), 1000)]
        }
        # fwrite(data.frame(FID = names(phe), IID = names(phe), prot = phe), 
        #        paste0('cis_trans_inter/01_obs_phe/', prot_i, '.txt'), sep = '\t', col.names = F)
        ## make plink file
	system(paste0(plink_dir, "plink --bfile ",genotype_dir,
		      "chr_all --keep cis_trans_inter/01_obs_phe/", prot_i,
                      ".txt --make-bed --thread-num 16 --out cis_trans_inter/03_geno_grm/", prot_i))
        
        ## calculate grm using gcta
        system(paste(gcta_dir, "gcta64 --bfile cis_trans_inter/03_geno_grm/", prot_i,
                     " --make-grm --out cis_trans_inter/03_geno_grm/", prot_i, sep = ""))
        ## calculate cvgp and h2
        system(paste0(gcta_dir, 
                      "gcta64 --reml --thread-num 16 --reml-no-constrain --grm cis_trans_inter/03_geno_grm/", 
                      prot_i, " --pheno cis_trans_inter/01_obs_phe/", prot_i, ".txt", 
                      " --qcovar cis_trans_inter/00_ref/gcta_cov_less.qcovar --out cis_trans_inter/07_h2/whole_h2/", 
                      prot_i))
      }, error=function(e){})
  }
  print(flag_i)
}




