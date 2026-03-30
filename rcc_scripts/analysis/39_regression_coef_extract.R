setwd('/project/xuanyao/jinghui')
library(data.table)
library(sandwich)
library(plink2R)
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')
pos_gene = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')

# expression matrix prot_ex; protein name gene_name; and sample id smaple_id
load('cis_trans_inter/00_ref/ukb_ppp_prot_expr.rdata') 
colnames(prot_ex) = gene_name

# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

### check the detailed regression of a protein-pair
i = 6
cis_gene_i = sig_inter$target[i]
trans_gene_j = sig_inter$regulator[i]

cis_ex_i = prot_ex[, cis_gene_i]
names(cis_ex_i) = sample_id
cis_ex_i = na.omit(cis_ex_i)
# cis prediction of target prot
cis_pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', 
                          cis_gene_i, '_good_pred.txt'))
cis_pred_i$pred = (cis_pred_i$pred - mean(cis_pred_i$pred))/
  sd(cis_pred_i$pred)
  
# use common individuals between cis and trans prot
common_id = Reduce(intersect, list(names(cis_ex_i), cis_pred_i$id, st_cov$V1))
cis_ex_i = cis_ex_i[as.character(common_id)]
st_cov_ij = st_cov[match(common_id, st_cov$V1), -(1:2)]

cis_ex_res = residuals(lm(cis_ex_i ~ as.matrix(st_cov_ij)))

# cis prediction of regulator prot
trans_pred_j = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', 
                            trans_gene_j, '_good_pred.txt'))
trans_pred_j$pred = (trans_pred_j$pred - mean(trans_pred_j$pred))/
  sd(trans_pred_j$pred)
common_id2 = intersect(trans_pred_j$id, common_id)

cis_prot_j = cis_pred_i$pred[match(common_id2, cis_pred_i$id)]
trans_prot_j = trans_pred_j$pred[match(common_id2, trans_pred_j$id)]
cis_ex_res_j = cis_ex_res[match(common_id2, common_id)]
    
## OLS of linear model
test=lm(cis_ex_res_j ~ cis_prot_j + trans_prot_j + cis_prot_j * trans_prot_j)
## using sandwich to correct for heteroscedasticity
sandwich_se = diag(vcovHC(test, type = "HC"))^0.5
sandwich_t = coef(summary(test))[, 1]/sandwich_se
sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)

coef_test = summary(test)$coefficients
coef_test = cbind(coef_test, sandwich_se, sandwich_p)
rownames(coef_test) = c('intercept', cis_gene_i, trans_gene_j, 'interaction')
coef_test

### FUT2 secretor status
fut2_cis = read_plink('cis_trans_inter/00_ref/cis_geno/FUT2/FUT2', impute="avg")
fut2_cis = fut2_cis$bed
fut2_cis = fut2_cis[,1]
names(fut2_cis) = sapply(strsplit(names(fut2_cis), ':', fixed = T), '[', 1)
fut2_cis = fut2_cis[names(cis_ex_res_j)]

cix_ex_tab = data.frame(cis_ex = cis_ex_res_j, 
                        cis_prot = cis_prot_j,
                        trans_prot = trans_prot_j,
                        fut2 = fut2_cis)

test2=lm(cis_ex ~ cis_prot + trans_prot + cis_prot * trans_prot, 
         data = cix_ex_tab[cix_ex_tab$fut2 == 2, ])
## using sandwich to correct for heteroscedasticity
sandwich_se = diag(vcovHC(test2, type = "HC"))^0.5
sandwich_t = coef(summary(test2))[, 1]/sandwich_se
sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)

coef_test2 = summary(test2)$coefficients
coef_test2 = cbind(coef_test2, sandwich_se, sandwich_p)
rownames(coef_test2) = c('intercept', cis_gene_i, trans_gene_j, 'interaction')
coef_test2





