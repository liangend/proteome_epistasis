library(mvtnorm)
library(gaston)
library(bindata)
setwd('/project/xuanyao/jinghui')

args = commandArgs(trailingOnly = T)
rho_index = as.numeric(args[1])
h2_index = as.numeric(args[2])

### simulate genotypes
# Function to create n correlated SNP sets
rmvBinomial <- function(n, n_snp, size, maf, rho) {
  X <- replicate(n, {
    cov_mat = matrix(0, nrow = n_snp, ncol = n_snp)
    cov_mat[1,] = rho
    cov_mat[,1] = rho
    cov_mat[-1,-1] = rho^2
    diag(cov_mat) = 1
    colSums(rmvbin(size, rep(maf, n_snp), bincorr = cov_mat))
  })
  t(X)
}

# Function to plot a QQ plot
qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

# Function to normalize genotype
stand_norm = function(x){return((x - mean(x)) / sd(x))}

### Parameters of simulation
# two alleles
size = 2
# number of SNPs
n_snp = 4
# MAF of each SNP
maf = 0.3
# correlation between SNP1 and SNP2
rho_all = c(0, 0.2, 0.4, 0.6, 0.8)
rho = rho_all[rho_index]
# sample size
N = 1000
# training testing split
train_set = 1:(N/2)
test_set = (N/2 + 1):1000
# number of simulations
nsim=1000
# heritability
h2_all = c(0.2, 0.4, 0.6, 0.8)
h2 = h2_all[h2_index]

pval1 = c()
pval2 = c()
pval3 = c()
for(i in 1:nsim){
  # simulate SNPs x1 and x2 with cor = rho
  X = rmvBinomial(N, n_snp, size=size, maf, rho)
  
  X = apply(X, 2, stand_norm)

  # 1 causal snps of cis gene
  beta = sqrt(h2)
  y = beta * X[,1] + rnorm(N, 0, sqrt(1 - h2))
  
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
  train_fit = lm(y[train_set] ~ X[train_set, 2] + X[train_set, 3] + X[train_set, 4])
  coef_est = summary(train_fit)$coefficients
  
  # predict the cis gene expr using the trained model
  cis_pred = coef_est[2,1]*X[test_set, 2] + coef_est[3,1]*X[test_set, 3] + 
    coef_est[4,1]*X[test_set, 4]
  
  # interaction between cis gene expr and trans gene expr
  reg1 = lm(y[test_set] ~ cis_pred * trans_pred)
  # interaction between cis snp and trans gene expr
  reg2 = lm(y[test_set] ~ X[test_set, 2] * trans_pred)
  # interaction between cis snp and trans gene expr with cis_pred as an additive term
  reg3 = lm(y[test_set] ~ cis_pred + trans_pred + trans_pred * X[test_set, 2])
  
  pval1[i] = summary(reg1)$coefficients[4, 4]
  pval2[i] = summary(reg2)$coefficients[4, 4]
  pval3[i] = summary(reg3)$coefficients[5, 4]
}

all_p = list(pval1, pval2, pval3)
saveRDS(all_p, paste0('cis_trans_inter/04_false_positive_test/snp_in_ld/p_r2_', rho_index, 
                      '_h2_', h2_index, '.rds'))



