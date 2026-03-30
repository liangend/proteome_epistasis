library(data.table)
setwd('/gpfs/data/xliu-lab/jinghui')

cov_dir='ukb_ppp/covar.txt'
sample_dir='ukb_ppp/genotype/ukb_ppp_chr1_QC.fam'

## keep cov that have genotypes
all_sample = fread(sample_dir)
covar=fread(cov_dir)
common_id = intersect(covar$ID, all_sample$V1)
covar = covar[match(common_id, covar$ID), ]

# calculate age
covar$age = 2008 - covar$birth_year
covar = covar[, -3]

covar$ukb_centre = as.factor(covar$ukb_centre)
covar$geno_batch = as.factor(covar$geno_batch)

stsv = model.matrix(~., data = covar[,-1])
stsv = stsv[, -1]
# standardize continuous covariates
stsv[, (ncol(stsv)-20):ncol(stsv)] = apply(stsv[, (ncol(stsv)-20):ncol(stsv)], 2,
                                           function(x){(x-mean(x))/sd(x)})
# attach family ID and ind ID (they are the same)
stsv = as.data.frame(stsv)
stsv = cbind(covar$ID, covar$ID, stsv)

fwrite(stsv, 'cis_trans_inter/00_ref/st_gcta_cov.qcovar', sep = '\t', col.names = F)





