library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')
ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait.csv')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

id_i = ukb_trait$UKB_ID[i]
trait_i = ukb_trait$abbrev[i]
system(paste0("python3 /gpfs/data/ukb-share/extraction_scripts/extract_pheno.py ",
              id_i, " -n ", trait_i, " -c -t cis_trans_inter/00_ref/ukb_traits/"))
