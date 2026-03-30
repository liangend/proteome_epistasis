library(mvtnorm)
library(gaston)
library(bindata)
library(simtrait)
library(data.table)
setwd('/project/xuanyao/jinghui')

args = commandArgs(trailingOnly = T)
rep_index = as.numeric(args[1])

### simulate genotypes
# Function to create n correlated pairs
rmvBinomial <- function(n, size, p1, p2, rho) {
  X <- replicate(n, {
    colSums(rmvbin(size, c(p1,p2), bincorr=(1-rho)*diag(2)+rho))
  })
  t(X)
}

# Function to normalize genotype
stand_norm = function(x){return((x - mean(x)) / sd(x))}

### Parameters of simulation
# two alleles
size <- 2
# MAF two SNPs
maf1 <- 0.1
maf2 <- 0.1

maf3 <- 0.2
maf4 <- 0.2

maf5 <- 0.3
maf6 <- 0.3

# correlation between SNP1 and SNP2
rho_all = c(0, 0.2, 0.4, 0.6, 0.8)
# sample size
N = 1000
# training testing split
train_set = 1:(N/2)
test_set = (N/2 + 1):1000
# number of simulations
nsim=1000
# heritability
h2_all = c(0.05, 0.1, 0.2, 0.5)

for (n in 1:2) {
  rep_n = (rep_index-1)*2 + n
  set.seed(rep_n)
  for (rho_index in 1:length(rho_all)) {
    rho = rho_all[rho_index]
    p_inter1 = matrix(0, nrow = nsim, ncol = length(h2_all))
    p_inter2 = matrix(0, nrow = nsim, ncol = length(h2_all))
    p_inter3 = matrix(0, nrow = nsim, ncol = length(h2_all))
    
    p_trans_add1 = matrix(0, nrow = nsim, ncol = length(h2_all))
    p_trans_add2 = matrix(0, nrow = nsim, ncol = length(h2_all))
    p_trans_add3 = matrix(0, nrow = nsim, ncol = length(h2_all))
    
    for (j in 1:length(h2_all)) {
      h2 = h2_all[j]
      for(i in 1:nsim){
        # simulate SNPs x1 and x2 with cor = rho
        X1 = rmvBinomial(N, size=size, p1=maf1, p2=maf2, rho=rho)
        X2 = rmvBinomial(N, size=size, p1=maf3, p2=maf4, rho=rho)
        X3 = rmvBinomial(N, size=size, p1=maf5, p2=maf6, rho=rho)
        x4 = rbinom(N, 2, 0.5)
        
        X1 = apply(X1, 2, stand_norm)
        X2 = apply(X2, 2, stand_norm)
        X3 = apply(X3, 2, stand_norm)
        x4 = stand_norm(x4)
        # 4 causal snps of cis gene
        beta = rnorm(4)
        bv = beta[1] * X1[,1] + beta[2] * X2[,1] + beta[3] * X3[,1] + beta[4] * x4
        sigma_e = sqrt((1 - h2)/h2 * var(bv))
        
        # cis gene expr
        y = bv + rnorm(N, 0, sigma_e)
        
        # trans gene expr
        trans_x1 = rbinom(N/2, 2, 0.5)
        trans_x2 = rbinom(N/2, 2, 0.5)
        trans_x3 = rbinom(N/2, 2, 0.5)
        trans_x1 = stand_norm(trans_x1)
        trans_x2 = stand_norm(trans_x2)
        trans_x3 = stand_norm(trans_x3)
        beta_trans = rnorm(3)
        trans_pred = beta_trans[1] * trans_x1 + beta_trans[2] * trans_x2 + 
          beta_trans[3] * trans_x3
        
        # train the model using the training set
        train_fit = lm(y[train_set] ~ X1[train_set, 2] + X2[train_set, 2] + X3[train_set, 2])
        coef_est = summary(train_fit)$coefficients
        
        # predict the cis gene expr using the trained model
        cis_pred = coef_est[2,1]*X1[test_set, 2] + coef_est[3,1]*X2[test_set, 2] + 
          coef_est[4,1]*X3[test_set, 2]
        
        # interaction between cis gene expr and trans gene expr
        reg1 = lm(y[test_set] ~ cis_pred * trans_pred)
        # interaction between cis snp and trans gene expr
        reg2 = lm(y[test_set] ~ trans_pred * X1[test_set, 2])
        # interaction between cis snp and trans gene expr with cis_pred as an additive term
        reg3 = lm(y[test_set] ~ cis_pred + trans_pred + trans_pred * X1[test_set, 2])
        
        p_inter1[i, j] = summary(reg1)$coefficients[4, 4]
        p_inter2[i, j] = summary(reg2)$coefficients[4, 4]
        p_inter3[i, j] = summary(reg3)$coefficients[5, 4]
        
        p_trans_add1[i, j] = summary(reg1)$coefficients['trans_pred', 'Pr(>|t|)']
        p_trans_add2[i, j] = summary(reg2)$coefficients['trans_pred', 'Pr(>|t|)']
        p_trans_add3[i, j] = summary(reg3)$coefficients['trans_pred', 'Pr(>|t|)']
      }
      print(j)
    }
    all_p = list(p_inter1, p_inter2, p_inter3, 
                 p_trans_add1, p_trans_add2, p_trans_add3)
    p_infl = sapply(all_p, function(x){apply(x, 2, function(x){pval_infl(x, df = 1)})})
    p_infl=as.data.frame(p_infl)
    colnames(p_infl) = c('inter1', 'inter2', 'inter3', 
                         'trans_add1', 'trans_add2', 'trans_add3')
    rownames(p_infl) = c('h2_1', 'h2_2', 'h2_3', 'h2_4')
    
    saveRDS(all_p, paste0('cis_trans_inter/04_false_positive_test/expr_inter/ld_',
                          rho_index, '/p_', rep_n, '.rds'))
    fwrite(p_infl, paste0('cis_trans_inter/04_false_positive_test/expr_inter/ld_',
                          rho_index, '/p_infl', rep_n, '.txt'), sep = '\t', row.names = T)
  }
}





