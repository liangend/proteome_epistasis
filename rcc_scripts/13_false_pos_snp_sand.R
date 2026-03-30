library(mvtnorm)
library(bindata)
library(simtrait)
library(data.table)
library(sandwich)
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
size = 2
# MAF of SNP 1 and 2
maf1 = 0.2
maf2 = 0.2
# correlation between SNP1 and SNP2
rho_all = c(0, 0.3, 0.6, 0.9)

# sample size
N = 1000
# number of simulations
nsim = 1000
# heritability
h2 = c(0.05, 0.1, 0.2, 0.5)

for (n in 1:5) {
  rep_n = (rep_index-1)*5 + n
  #set.seed(rep_n)
  for (rho_index in 1:length(rho_all)) {
    rho = rho_all[rho_index]
    p_inter1 = matrix(0, nrow = nsim, ncol = length(h2))
    p_inter2 = matrix(0, nrow = nsim, ncol = length(h2))
    p_inter3 = matrix(0, nrow = nsim, ncol = length(h2))
    
    p_x3_add1 = matrix(0, nrow = nsim, ncol = length(h2))
    p_x3_add2 = matrix(0, nrow = nsim, ncol = length(h2))
    p_x3_add3 = matrix(0, nrow = nsim, ncol = length(h2))
    
    p_x2_add1 = matrix(0, nrow = nsim, ncol = length(h2))
    p_x2_add3 = matrix(0, nrow = nsim, ncol = length(h2))
    
    p_inter_sand1 = matrix(0, nrow = nsim, ncol = length(h2))
    p_inter_sand2 = matrix(0, nrow = nsim, ncol = length(h2))
    p_inter_sand3 = matrix(0, nrow = nsim, ncol = length(h2))

    for (j in 1:length(h2)) {
      h2_j = h2[j]
      for(i in 1:nsim){
        # simulate SNPs x1 and x2 with cor = rho
        X <- rmvBinomial(N, size=size, p1=maf1, p2=maf2, rho=rho)
        # independent SNP x3 (SNP in trans)
        x3 = rbinom(N, 2, 0.4)
        
        x1_norm = stand_norm(X[,1])
        x2_norm = stand_norm(X[,2])
        x3_norm = stand_norm(x3)
        
        # effect size of x1
        b1 = sqrt(h2_j / var(x1_norm))
        # phenotype
        y = b1 * x1_norm + rnorm(N, 0, 1 - h2_j)
        
        reg1 = lm(y ~ x2_norm * x3_norm)
        reg2 = lm(y ~ x1_norm * x3_norm)
        reg3 = lm(y ~ x1_norm + x2_norm * x3_norm)
        
        p_inter1[i, j] = summary(reg1)$coefficients["x2_norm:x3_norm", "Pr(>|t|)"]
        p_inter2[i, j] = summary(reg2)$coefficients["x1_norm:x3_norm", "Pr(>|t|)"]
        p_inter3[i, j] = summary(reg3)$coefficients["x2_norm:x3_norm", "Pr(>|t|)"]
        
        p_x3_add1[i, j] = summary(reg1)$coefficients["x3_norm", "Pr(>|t|)"]
        p_x3_add2[i, j] = summary(reg2)$coefficients["x3_norm", "Pr(>|t|)"]
        p_x3_add3[i, j] = summary(reg3)$coefficients["x3_norm", "Pr(>|t|)"]
        
        p_x2_add1[i, j] = summary(reg1)$coefficients["x2_norm", "Pr(>|t|)"]
        p_x2_add3[i, j] = summary(reg3)$coefficients["x2_norm", "Pr(>|t|)"]
        
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
        
        p_inter_sand1[i, j] = sandwich_p1[4]
        p_inter_sand2[i, j] = sandwich_p2[4]
        p_inter_sand3[i, j] = sandwich_p3[5]
        
      }
      print(j)
    }
    all_p = list(p_inter1, p_inter2, p_inter3, p_x3_add1, p_x3_add2, p_x3_add3, 
                 p_x2_add1, p_x2_add3, 
                 p_inter_sand1, p_inter_sand2, p_inter_sand3)
    p_infl = sapply(all_p, function(x){apply(x, 2, function(x){pval_infl(x, df = 1)})})
    p_infl=as.data.frame(p_infl)
    colnames(p_infl) = c('inter1', 'inter2', 'inter3', 'x3_add1', 'x3_add2',
                         'x3_add3', 'x2_add1', 'x2_add3',
                         'inter_sand1', 'inter_sand2', 'inter_sand3')
    rownames(p_infl) = c('h2_1', 'h2_2', 'h2_3', 'h2_4')
    
    saveRDS(all_p, paste0('cis_trans_inter/04_false_positive_test/snp_sandwich_inter/ld_',
                          rho_index, '/p_', rep_n, '.rds'))
    fwrite(p_infl, paste0('cis_trans_inter/04_false_positive_test/snp_sandwich_inter/ld_',
                          rho_index, '/p_infl', rep_n, '.txt'), sep = '\t', row.names = T)
    
  }
}








