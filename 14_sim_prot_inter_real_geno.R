setwd('/project/xuanyao/jinghui')
library(simtrait)
library(data.table)
library(sandwich)
library(plink2R)

args = commandArgs(trailingOnly = T)
rep_index = as.numeric(args[1])

gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta = gene_meta[gene_meta$gene_type == 'protein_coding', ]
gene_meta = gene_meta[gene_meta$chr %in% paste0('chr', 1:22), ]

### simulate genotypes

### Parameters of simulation

# correlation between SNP1 and SNP2
# sample size
N = 3000
# training testing split
train_set = 1:1500
test_set = 1501:N
# number of simulations
nsim=1000
# heritability
h2_all = c(0.05, 0.1, 0.2, 0.5)


p_inter1 = matrix(NA, nrow = nsim, ncol = length(h2_all))
p_inter1_sand = matrix(NA, nrow = nsim, ncol = length(h2_all))
ld_mat = matrix(NA, nrow = nsim, ncol = length(h2_all))

set.seed(101)
gene_meta_sample = gene_meta[sample(nrow(gene_meta), nsim), ]
gene_meta_sample = gene_meta_sample[order(gene_meta_sample$chr), ]

for (j in 1:length(h2_all)) {
  h2 = h2_all[j]
  chr_ind = 'chrx'
  for(i in 1:nsim){
    set.seed(rep_index + i + j)
    ##### simulate cis gene expr
    ## randomly pick a gene
    chr1 = gene_meta_sample$chr[i]
    tss1 = gene_meta_sample$start[i]
    if (chr1 != chr_ind) {
      ## 1000G genotype 
      geno1 = read_plink(paste0('cis_trans_inter/00_ref/ukb_eur_geno/', chr1))
      geno1_bim = geno1$bim
      geno1_bed = geno1$bed
      chr_ind = chr1
      gc()
    }
    
    ## cis region (20 kb) of the gene
    geno1_bim_cis = geno1_bim[abs(geno1_bim$V4 - tss1) < 100000, ]
    if (nrow(geno1_bim_cis) < 10) {
      next
    }
    cis_geno = geno1_bed[, geno1_bim_cis$V2]
    cis_geno = apply(cis_geno, 2, function(x) {
      x[is.na(x)] = mean(x, na.rm = T)
      return(x)
    })
    cor_cis = abs(cor(cis_geno))
    ave_cor = apply(cor_cis, 1, mean)
    
    snp_pool = colnames(cis_geno)
    
    ## 5 causal SNPs
    causal_snp = sample(names(snp_pool), 5)
    causal_geno = cis_geno[, causal_snp]
    causal_geno = apply(causal_geno, 2, scale)
    
    ## 5 tagged SNPs
    cor_causal = cor_cis[-which(rownames(cor_cis) %in% causal_snp), 
                         causal_snp]
    max_cor = apply(cor_causal, 1, max)
    max_cor = max_cor[max_cor > 0.2]
    
    ld_snp_cand = intersect(names(max_cor), names(snp_pool))
    
    if (length(ld_snp_cand) < 5) {
      next
    }
    ld_snp = sample(ld_snp_cand, 5)
    ld_geno = cis_geno[, ld_snp]
    ld_geno = apply(ld_geno, 2, scale)
    
    cor_geno = cis_geno[, c(causal_snp, ld_snp)]
    
    mean_ld = mean(abs(cor(cor_geno)[6:10, 1:5])) # average cor between tagged and causal snps
    ld_mat[i, j] = mean_ld
    
    ## cis gene expr
    beta = rnorm(5)
    bv = causal_geno %*% beta
    sigma_e = sqrt((1 - h2)/h2 * var(bv))
    y = bv + rnorm(N, 0, sigma_e)
    y = as.vector(y)
    
    ##### simulate trans gene expr
    trans_x1 = rbinom(length(test_set), 2, 0.5)
    trans_x2 = rbinom(length(test_set), 2, 0.5)
    trans_x3 = rbinom(length(test_set), 2, 0.5)
    trans_x1 = scale(trans_x1)
    trans_x2 = scale(trans_x2)
    trans_x3 = scale(trans_x3)
    beta_trans = rnorm(3)
    trans_pred = beta_trans[1] * trans_x1 + beta_trans[2] * trans_x2 + 
      beta_trans[3] * trans_x3
    
    ###### train the model using the training set
    train_fit = lm(y[train_set] ~ ld_geno[train_set, ])
    coef_est = summary(train_fit)$coefficients[-1, 1]
    ## make singular coefficient 0
    names(coef_est) = sub("ld_geno\\[train_set, \\]", '', names(coef_est))
    coef_est = coef_est[match(ld_snp, names(coef_est))]
    coef_est[is.na(coef_est)] = 0
    
    # predict the cis gene expr using the trained model
    cis_pred = ld_geno[test_set, ] %*% coef_est
    cis_pred = as.vector(cis_pred)
    
    # interaction between cis gene expr and trans gene expr
    fit1 = lm(y[test_set] ~ cis_pred * trans_pred)
    fit_coef = summary(fit1)$coefficients
    if (!all(dim(fit_coef) == c(4,4))) {
      next
    }  
    p_inter1[i, j] = fit_coef[4, 4]
    
    sandwich_se1 = diag(vcovHC(fit1, type = "HC"))^0.5
    sandwich_t1 = coef(summary(fit1))[, 1]/sandwich_se1
    sandwich_p1 = pchisq(sandwich_t1^2, 1, lower.tail = FALSE)
    p_inter1_sand[i, j] = sandwich_p1[4]
  }
  gc()
  print(j)
}
all_p = list(p_inter1, p_inter1_sand)
p_infl = sapply(all_p, function(x){apply(x, 2, function(x){pval_infl(x, df = 1)})})


p_infl=as.data.frame(p_infl)
colnames(p_infl) = c('inter', 'inter_sve')
rownames(p_infl) = c('h2_1', 'h2_2', 'h2_3', 'h2_4')

saveRDS(all_p, paste0('cis_trans_inter/04_false_positive_test/expr_inter_ukb/p_', 
                      rep_index, '.rds'))
saveRDS(ld_mat, paste0('cis_trans_inter/04_false_positive_test/expr_inter_ukb/ld_', 
                       rep_index, '.rds'))
fwrite(p_infl, paste0('cis_trans_inter/04_false_positive_test/expr_inter_ukb/p_infl', 
                      rep_index, '.txt'), sep = '\t', row.names = T)





