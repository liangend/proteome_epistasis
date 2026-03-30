library(data.table)
library(sandwich)
setwd('/gpfs/data/xliu-lab/jinghui')
args=commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

ukb_trait = fread('cis_trans_inter/00_ref/ukb_trait_share.txt') # 56 traits
sig_inter = fread('cis_trans_inter/06_fusion_inter/sig_prot_prot_inter.txt')
all_cov = fread('cis_trans_inter/00_ref/all_norm_cov.txt')
all_cov = na.omit(all_cov)

## phenotype
trait_i = ukb_trait$abbrev[i]
file_i = ukb_trait$file[i]
type_i = ukb_trait$Type[i]

pheno_i = fread(paste0('/gpfs/data/ukb-share/extracted_phenotypes/', trait_i, 
                       '/', file_i))
pheno_i = na.omit(pheno_i)
colnames(pheno_i)[3] = 'pheno'

## remove negative values representing missing data
pheno_i = pheno_i[pheno_i$pheno >= 0, ]

if (type_i == 'Continuous') {
  pheno_i$pheno = (pheno_i$pheno - mean(pheno_i$pheno)) / sd(pheno_i$pheno)
} else {
  pheno_i$pheno[pheno_i$pheno > 0] = 1
}

abo_blood = fread('cis_trans_inter/00_ref/ukb_traits/blood_type/blood_type674178.pheno')
colnames(abo_blood)[3] = 'blood_type'
abo_blood = na.omit(abo_blood)

abo_blood$o_type = 'OO'
abo_blood$o_type[which(abo_blood$blood_type %in% c('BO', 'AO'))] = 'OX'
abo_blood$o_type[which(abo_blood$blood_type %in% c('AA', 'BB', 'AB'))] = 'XX'

common_id = Reduce(intersect, list(all_cov$FID, pheno_i$FID, abo_blood$FID))
pheno_i = pheno_i[match(common_id, pheno_i$FID), ]
all_cov = all_cov[match(common_id, all_cov$FID), -(1:2)]


## residualize phenotypes by regressing on covariates
if (type_i == 'Continuous') {
  res_i = residuals(lm(unlist(pheno_i$pheno) ~ as.matrix(all_cov)))
} else {
  res_i = residuals(glm(unlist(pheno_i$pheno) ~ as.matrix(all_cov), family = binomial))
}
names(res_i) = common_id

## regression between phenotype and interaction
inter_coef = list()
for (n in 1:nrow(sig_inter)) {
  target_n = sig_inter$target[n]
  target_pred = fread(paste0('cis_trans_inter/09_inter_pheno_asso/pred_prot/',
                             target_n, '_pred.txt'))
  
  ## remove ids used in FUSION training
  target_train_id = fread(paste0('cis_trans_inter/05_fusion/obs_phe/eur/',
                                 target_n, '.txt'))
  train_id = target_train_id$FID
  common_id2 = target_pred$id[!target_pred$id %in% train_id]
  
  res_n = res_i[match(common_id2, names(res_i))]
  target_n = target_pred$pred[match(common_id2, target_pred$id)]
  o_type_n = abo_blood$o_type[match(common_id2, abo_blood$FID)]
  
  fit1 = lm(res_n ~ target_n * o_type_n)
  fit_coef = summary(fit1)$coefficients
  sandwich_se = diag(vcovHC(fit1, type = "HC"))^0.5
  sandwich_t = fit_coef[, 1]/sandwich_se
  sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  fit_coef = cbind(fit_coef, sandwich_t = sandwich_t, sandwich_p = sandwich_p)
  
  inter_coef[[n]] = fit_coef
}
names(inter_coef) = paste0(sig_inter$target, '*', sig_inter$regulator)
saveRDS(inter_coef, paste0('cis_trans_inter/15_pheno_asso_inter_otype/', 
                         trait_i, '.rds'))
