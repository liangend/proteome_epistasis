library(mvtnorm)
library(gaston)
library(bindata)
library(sandwich)
library(simtrait)
library(data.table)
setwd('/project/xuanyao/jinghui')

args = commandArgs(trailingOnly = T)
rep_index = as.numeric(args[1])

# Function to normalize genotype
stand_norm = function(x){return((x - mean(x)) / sd(x))}

### Parameters of simulation
# two alleles
size = 2
# MAF of 4 SNP
maf1 = 0.2
maf2 = 0.2

maf3 = 0.3
maf4 = 0.3

# sample size
N = 1000

# training testing split
train_set = 1:(N/2)
test_set = (N/2 + 1):1000

# number of simulations
nsim = 1000

# heritability
h2_all = c(0.05, 0.1, 0.15, 0.2)
inter_h2_all = c(0.18, 0.13, 0.08, 0.03)
trans_h2 = 0.2

for (n in 1:5) {
  rep_n = (rep_index-1)*5 + n
  set.seed(rep_n)
  inter_p1 = matrix(0, nrow = nsim, ncol = length(h2_all))
  inter_p2 = matrix(0, nrow = nsim, ncol = length(h2_all))
  inter_p3 = matrix(0, nrow = nsim, ncol = length(h2_all))
  
  inter_p_sand1 = matrix(0, nrow = nsim, ncol = length(h2_all))
  inter_p_sand2 = matrix(0, nrow = nsim, ncol = length(h2_all))
  inter_p_sand3 = matrix(0, nrow = nsim, ncol = length(h2_all))
  for (j in 1:4) {
    h2 = h2_all[j]
    inter_h2 = inter_h2_all[j]
    
    for(i in 1:nsim){
      # simulate SNPs x1 and x2 with cor = rho
      x1 = rbinom(N, 2, maf1)
      x2 = rbinom(N, 2, maf2)
      x3 = rbinom(N, 2, maf3)
      x4 = rbinom(N, 2, maf4)
      
      x1 = stand_norm(x1)
      x2 = stand_norm(x2)
      x3 = stand_norm(x3)
      x4 = stand_norm(x4)
      
      # 4 causal snps of cis gene
      beta = rnorm(4)
      bv = beta[1] * x1 + beta[2] * x2 + beta[3] * x3 + beta[4] * x4
      
      # trans gene expr
      trans_x1 = rbinom(N, 2, 0.5)
      trans_x2 = rbinom(N, 2, 0.5)
      trans_x3 = rbinom(N, 2, 0.5)
      trans_x1 = stand_norm(trans_x1)
      trans_x2 = stand_norm(trans_x2)
      trans_x3 = stand_norm(trans_x3)
      beta_trans = rnorm(3)
      
      trans_bv = beta_trans[1] * trans_x1 + beta_trans[2] * trans_x2 + 
        beta_trans[3] * trans_x3
      
      sigma_e_trans = sqrt((1 - trans_h2)/trans_h2 * var(trans_bv))
      trans_y = trans_bv + rnorm(N, 0, sigma_e_trans)
      
      inter_bv = trans_y * bv
      var_inter = var(bv) * (inter_h2 / h2)
      inter_bv = inter_bv * sqrt(var_inter / var(inter_bv))
      sigma_e = sqrt((1 - h2 - inter_h2)/(h2 + inter_h2) * var(bv))
      # cis gene expr
      y = bv + inter_bv + rnorm(N, 0, sigma_e) 
      
      # train the model using the training set
      cis_fit = lm(y[train_set] ~ x1[train_set] + x2[train_set] + x3[train_set])
      coef_cis = summary(cis_fit)$coefficients
      
      trans_fit = lm(trans_y[train_set] ~ trans_x1[train_set] + trans_x2[train_set] + trans_x3[train_set])
      coef_trans = summary(trans_fit)$coefficients
      
      # predict the cis gene expr using the trained model
      cis_pred = coef_cis[2,1]*x1[test_set] + coef_cis[3,1]*x2[test_set] + 
        coef_cis[4,1]*x3[test_set]
      trans_pred = coef_trans[2,1]*trans_x1[test_set] + coef_trans[3,1]*trans_x2[test_set] + 
        coef_trans[4,1]*trans_x3[test_set]
      
      # interaction between cis gene expr and trans gene expr
      reg1 = lm(y[test_set] ~ cis_pred * trans_pred)
      # interaction between cis snp and trans gene expr
      reg2 = lm(y[test_set] ~ trans_pred * x1[test_set])
      # interaction between cis snp and trans gene expr with cis_pred as an additive term
      reg3 = lm(y[test_set] ~ cis_pred + trans_pred + trans_pred * x1[test_set])
      
      ## regular interaction p value
      inter_p1[i,j] = summary(reg1)$coefficients[4, 4]
      inter_p2[i,j] = summary(reg2)$coefficients[4, 4]
      inter_p3[i,j] = summary(reg3)$coefficients[5, 4]
      
      ## using sandwich to correct for heteroscedasticity
      sandwich_se1 = diag(vcovHC(reg1, type = "HC"))^0.5
      sandwich_t1 = coef(summary(reg1))[, 1]/sandwich_se1
      sandwich_p1 = pchisq(sandwich_t1^2, 1, lower.tail = FALSE)
      
      sandwich_se2 = diag(vcovHC(reg2, type = "HC"))^0.5
      sandwich_t2 = coef(summary(reg2))[, 1]/sandwich_se2
      sandwich_p2 = pchisq(sandwich_t2^2, 1, lower.tail = FALSE)
      
      sandwich_se3 = diag(vcovHC(reg3, type = "HC"))^0.5
      sandwich_t3 = coef(summary(reg3))[, 1]/sandwich_se3
      sandwich_p3 = pchisq(sandwich_t3^2, 1, lower.tail = FALSE)
      
      inter_p_sand1[i,j] = sandwich_p1[4]
      inter_p_sand2[i,j] = sandwich_p2[4]
      inter_p_sand3[i,j] = sandwich_p3[5]
    }
  }
  all_p = list(inter_p1, inter_p2, inter_p3, 
               inter_p_sand1, inter_p_sand2, inter_p_sand3)
  p_infl = sapply(all_p, function(x){apply(x, 2, function(x){pval_infl(x, df = 1)})})
  p_infl=as.data.frame(p_infl)
  colnames(p_infl) = c('inter1', 'inter2', 'inter3', 
                       'inter_sand1', 'inter_sand2', 'inter_sand3')
  rownames(p_infl) = c('h2_1', 'h2_2', 'h2_3', 'h2_4')
  
  saveRDS(all_p, paste0('cis_trans_inter/05_power_test/gene_gene_inter/obs_sand/p_', rep_n, '.rds'))
  fwrite(p_infl, paste0('cis_trans_inter/05_power_test/gene_gene_inter/obs_sand/p_infl', rep_n, '.txt'), 
         sep = '\t', row.names = T)
}



