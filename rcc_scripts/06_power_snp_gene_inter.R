library(mvtnorm)
library(gaston)
library(bindata)
setwd('/project/xuanyao/jinghui')

args = commandArgs(trailingOnly = T)
rho_index = as.numeric(args[1])
h2_index = as.numeric(args[2])

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
# MAF of SNP 1 and 2
maf1 <- 0.2
maf2 <- 0.2

maf3 <- 0.2
maf4 <- 0.2

maf5 <- 0.3
maf6 <- 0.3

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
h2_all = c(0.2, 0.3, 0.4, 0.5)
h2 = h2_all[h2_index]

inter_h2_all = c(0.5, 0.4, 0.3, 0.2)
inter_h2 = inter_h2_all[h2_index]

trans_h2 = 0.3

pval1 = c()
pval2 = c()
pval3 = c()
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
  
  inter_bv = trans_y * X1[,1]
  var_inter = var(bv) * (inter_h2 / h2)
  inter_bv = inter_bv * sqrt(var_inter / var(inter_bv))
  sigma_e = sqrt((1 - h2 - inter_h2)/(h2 + inter_h2) * var(bv))
  # cis gene expr
  y = bv + inter_bv + rnorm(N, 0, sigma_e) 
  
  # train the model using the training set
  cis_fit = lm(y[train_set] ~ X1[train_set, 2] + X2[train_set, 2] + X3[train_set, 2])
  coef_cis = summary(cis_fit)$coefficients
  
  trans_fit = lm(trans_y[train_set] ~ trans_x1[train_set] + trans_x2[train_set] + trans_x3[train_set])
  coef_trans = summary(trans_fit)$coefficients
  
  # predict the cis gene expr using the trained model
  cis_pred = coef_cis[2,1]*X1[test_set, 2] + coef_cis[3,1]*X2[test_set, 2] + 
    coef_cis[4,1]*X3[test_set, 2]
  trans_pred = coef_trans[2,1]*trans_x1[test_set] + coef_trans[3,1]*trans_x2[test_set] + 
    coef_trans[4,1]*trans_x3[test_set]
  
  # interaction between cis gene expr and trans gene expr
  reg1 = lm(y[test_set] ~ cis_pred * trans_pred)
  # interaction between cis snp and trans gene expr
  reg2 = lm(y[test_set] ~ trans_pred * X1[test_set, 2])
  # interaction between cis snp and trans gene expr with cis_pred as an additive term
  reg3 = lm(y[test_set] ~ cis_pred + trans_pred + trans_pred * X1[test_set, 2])
  
  pval1[i] = summary(reg1)$coefficients[4, 4]
  pval2[i] = summary(reg2)$coefficients[4, 4]
  pval3[i] = summary(reg3)$coefficients[5, 4]
  if (i %% 100 == 0) {
    print(i)
  }
}

all_p = list(pval1, pval2, pval3)
saveRDS(all_p, paste0('cis_trans_inter/05_power_test/snp_gene_inter/p_r2_', rho_index, 
                      '_h2_', h2_index, '.rds'))



