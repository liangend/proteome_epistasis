library(mvtnorm)
library(bindata)
library(simtrait)
library(data.table)
library(sandwich)
library(ggplot2)

setwd('/project/xuanyao/jinghui')

### simulate genotypes
# Function to create n correlated pairs
rmvBinomial <- function(n, size, p1, p2, rho) {
  X <- replicate(n, {
    colSums(rmvbin(size, c(p1,p2), bincorr=(1-rho)*diag(2)+rho))
  })
  t(X)
}

### Parameters of simulation
# two alleles
size <- 2
# MAF of SNP 1 and 2
maf1 <- 0.2
maf2 <- 0.2

# sample size
N = 1000

rho = 0.8
h2_j = 0.5

inter_p = c()
sand_p = c()
var_x1 = c()
var_x2 = c()
for (i in 1:30) {
  # simulate SNPs x1 and x2 with cor = rho
  X <- rmvBinomial(N, size=size, p1=maf1, p2=maf2, rho=rho)
  x1 = X[,1]
  x2 = X[,2]
  
  # independent SNP x3 (SNP in trans)
  x3 = rbinom(N, 2, 0.4)
  
  # effect size of x1
  b1 = sqrt(h2_j / var(x1))
  # phenotype
  y = b1 * x1 + rnorm(N, 0, 1 - h2_j)
  
  sim_tab = data.frame(x1 = x1, x2 = x2, 
                       x3 = x3, 
                       y = y)
  var_x1 = rbind(var_x1, aggregate(y ~ x1, data = sim_tab, var)[, 2])
  var_x2 = rbind(var_x2, aggregate(y ~ x2, data = sim_tab, var)[, 2])
  reg1 = lm(y ~ x2 * x3)
  inter_p[i] = summary(reg1)$coefficients[4,4]
  
  sandwich_se = diag(vcovHC(reg1, type = "HC"))^0.5
  sandwich_t = coef(summary(reg1))[, 1]/sandwich_se
  sand_p[i] = pchisq(sandwich_t^2, 1, lower.tail = FALSE)[4]
  print(i)
}
plot(inter_p, sand_p, 
     xlab = "p value without SVE", ylab = 'p value with SVE',
     main = 'MAF = 0.2')
abline(0, 1, lty = 2)

var_x1_plt = data.frame(var = unlist(as.data.frame(var_x1)), 
                        snp = rep(c(0,1,2), each = 30))
var_x2_plt = data.frame(var = unlist(as.data.frame(var_x2)), 
                        snp = rep(c(0,1,2), each = 30))

ggplot(var_x1_plt, aes(x = as.factor(snp), y = var)) + 
  geom_boxplot() +
  labs(title = 'MAF = 0.5', x = "SNP1", y = "Var(y)", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggplot(var_x2_plt, aes(x = as.factor(snp), y = var)) + 
  geom_boxplot() +
  labs(title = 'MAF = 0.5', x = "SNP2", y = "Var(y)", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

table(x1, x2)
aggregate(y ~ x2, data = sim_tab, var)
ggplot(sim_tab, aes(x = as.factor(x1), y = y)) + 
  geom_violin() + geom_point() + 
  labs(title = '', x = "SNP1", y = "Var(y)", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



