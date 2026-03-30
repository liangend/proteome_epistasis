setwd('/project/xuanyao/jinghui')
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gaston)
library(simtrait)

qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))

  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

### Inflation of PWAS y ~ protA + protA*protB when y is only affected by protB
nround = 100
nsim = 1000
N = 20000 # number of individuals
h2_reg = 0.2 # regulator h2
h2_target_cis = 0.1 # target cis var is 0.1 of the total genetic var
h2_target_trans = 0.8 # target trans var is 0.8 of the total genetic var
h2_target_inter = 0.1 # target interaction var is 0.1 of the total genetic var
res_coef = 4 # residual var is 4 times genetic var
train_set = 1:5000
test_set = 5001:10000
pwas_set = 10001:20000

lambda_raw = c()
lambda_wo_reg = c()
lambda_w_reg = c()
lambda_wo_reg_center = c()
lambda_w_reg_center = c()

fpr_raw = c()
fpr_wo_reg = c()
fpr_w_reg = c()
fpr_wo_reg_center = c()
fpr_w_reg_center = c()

for (n in 1:nround) {
  pwas_coef_raw = c()
  pwas_coef_wo_reg = c()
  pwas_p_wo_reg = c()
  pwas_coef_w_reg = c()
  pwas_p_w_reg = c()
 
  pwas_p_raw = c()
  pwas_coef_wo_reg_center = c()
  pwas_p_wo_reg_center = c()
  pwas_coef_w_reg_center = c()
  pwas_p_w_reg_center = c()
  
  for (i in 1:nsim) {
    reg_snp1 = rbinom(N, 2, 0.2) # regulator cis SNP1
    reg_snp2 = rbinom(N, 2, 0.2) # regulator cis SNP2
    reg_snp3 = rbinom(N, 2, 0.2) # regulator cis SNP3
    
    target_snp1 = rbinom(N, 2, 0.2) # target cis SNP1
    target_snp2 = rbinom(N, 2, 0.2) # target cis SNP2
    target_snp3 = rbinom(N, 2, 0.2) # target cis SNP3
    
    beta_reg = rnorm(3)
    reg_cis = reg_snp1*beta_reg[1] + reg_snp2*beta_reg[2] + reg_snp3*beta_reg[3]
    var_reg_e = (1 - h2_reg)/h2_reg * var(reg_cis)
    reg = reg_cis + rnorm(N, 0, sqrt(var_reg_e))
    
    beta_target = rnorm(3)
    target_cis = target_snp1*beta_target[1] + target_snp2*beta_target[2] +
      target_snp3*beta_target[3]
    
    target_cis_var = var(target_cis)
    add_eff = sqrt(target_cis_var * (h2_target_trans/h2_target_cis) / var(reg))  # regulator effect on target
    inter_eff = sqrt(target_cis_var * (h2_target_inter/h2_target_cis) / var(reg * target_cis)) # interaction effect
    var_target_e = res_coef * var(target_cis + add_eff * reg + inter_eff * reg * target_cis)
    target = target_cis + add_eff * reg +
      inter_eff * reg * target_cis + rnorm(N, 0, sqrt(var_target_e)) # target expression
    ## normalize expression
    reg = (reg - mean(reg))/sd(reg)
    target = (target - mean(target))/sd(target)
    
    trait_eff = 5 # regulator effect on the phenotype. In real data, slope of ALP on ABO expression is -5
    Y = trait_eff * reg + rnorm(N,0,10) # Y is the phenotype. In real data, the regression residual var is 26
    Y_center = Y - mean(Y)
    
    ## training testing split
    reg_train = reg[train_set]
    reg_snp1_train = reg_snp1[train_set]
    reg_snp2_train = reg_snp2[train_set]
    reg_snp3_train = reg_snp3[train_set]
    
    reg_snp1_test = reg_snp1[test_set]
    reg_snp2_test = reg_snp2[test_set]
    reg_snp3_test = reg_snp3[test_set]
    
    reg_snp1_pwas = reg_snp1[pwas_set]
    reg_snp2_pwas = reg_snp2[pwas_set]
    reg_snp3_pwas = reg_snp3[pwas_set]
    
    target_train = target[train_set]
    target_snp1_train = target_snp1[train_set]
    target_snp2_train = target_snp2[train_set]
    target_snp3_train = target_snp3[train_set]
    
    target_test = target[test_set]
    target_snp1_test = target_snp1[test_set]
    target_snp2_test = target_snp2[test_set]
    target_snp3_test = target_snp3[test_set]
    
    target_snp1_pwas = target_snp1[pwas_set]
    target_snp2_pwas = target_snp2[pwas_set]
    target_snp3_pwas = target_snp3[pwas_set]
    
    ## predicting cis-regulatory component
    fit_reg = lm(reg_train ~ reg_snp1_train + reg_snp2_train + reg_snp3_train)
    reg_coef = summary(fit_reg)$coefficients
    
    reg_hat = reg_coef[1,1] + reg_coef[2,1] * reg_snp1_test + reg_coef[3,1] * reg_snp2_test +
      reg_coef[4,1] * reg_snp3_test
    reg_hat_center = reg_hat - mean(reg_hat)
    
    reg_hat_pwas = reg_coef[1,1] + reg_coef[2,1] * reg_snp1_pwas + reg_coef[3,1] * reg_snp2_pwas +
      reg_coef[4,1] * reg_snp3_pwas
    reg_hat_pwas_center = reg_hat_pwas - mean(reg_hat_pwas)
    
    fit_target = lm(target_train ~ target_snp1_train + target_snp2_train + target_snp3_train)
    target_coef = summary(fit_target)$coefficients
    
    target_hat = target_coef[1,1] + target_coef[2,1] * target_snp1_test + target_coef[3,1] * target_snp2_test +
      target_coef[4,1] * target_snp3_test
    target_hat_center = target_hat - mean(target_hat)
    
    target_hat_pwas = target_coef[1,1] + target_coef[2,1] * target_snp1_pwas + target_coef[3,1] * target_snp2_pwas +
      target_coef[4,1] * target_snp3_pwas
    target_hat_pwas_center = target_hat_pwas - mean(target_hat_pwas)
    
    ## predicting target by incorporating epistasis
    fit_target_epi = lm(target_test ~ target_hat * reg_hat)
    fit_target_epi_center = lm(target_test ~ target_hat_center * reg_hat_center)
    target_epi_coef = summary(fit_target_epi)$coefficients
    target_epi_center_coef = summary(fit_target_epi_center)$coefficients
    # without regulator additive part
    target_pwas1 = target_epi_coef[1,1] + target_epi_coef[2,1] * target_hat_pwas +
      target_epi_coef[4,1] * target_hat_pwas * reg_hat_pwas
    target_pwas1_center = target_epi_center_coef[1,1] + target_epi_center_coef[2,1] * target_hat_pwas_center +
      target_epi_center_coef[4,1] * target_hat_pwas_center * reg_hat_pwas_center
    
    # with regulator additive part
    target_pwas2 = target_epi_coef[1,1] + target_epi_coef[2,1] * target_hat_pwas +
      target_epi_coef[3,1] * reg_hat_pwas + target_epi_coef[4,1] * target_hat_pwas * reg_hat_pwas
    target_pwas2_center = target_epi_center_coef[1,1] + target_epi_center_coef[2,1] * target_hat_pwas_center +
      target_epi_center_coef[3,1] * reg_hat_pwas_center +
      target_epi_center_coef[4,1] * target_hat_pwas_center * reg_hat_pwas_center
    
    ## PWAS of two models
    Y_pwas = Y[pwas_set]
    
    pwas0 = lm(Y_pwas ~ target_hat_pwas)
    pwas1 = lm(Y_pwas ~ target_pwas1)
    pwas2 = lm(Y_pwas ~ target_pwas2)
    
    pwas1_center = lm(Y_pwas ~ target_pwas1_center)
    pwas2_center = lm(Y_pwas ~ target_pwas2_center)
    
    pwas_coef_raw[i] = summary(pwas0)$coefficients[2,1]
    pwas_p_raw[i] = summary(pwas0)$coefficients[2,4]
    
    pwas_coef_wo_reg[i] = summary(pwas1)$coefficients[2,1]
    pwas_p_wo_reg[i] = summary(pwas1)$coefficients[2,4]
    
    pwas_coef_w_reg[i] = summary(pwas2)$coefficients[2,1]
    pwas_p_w_reg[i] = summary(pwas2)$coefficients[2,4]
    
    pwas_coef_wo_reg_center[i] = summary(pwas1_center)$coefficients[2,1]
    pwas_p_wo_reg_center[i] = summary(pwas1_center)$coefficients[2,4]
    
    pwas_coef_w_reg_center[i] = summary(pwas2_center)$coefficients[2,1]
    pwas_p_w_reg_center[i] = summary(pwas2_center)$coefficients[2,4]
  }
  lambda_raw[n] = pval_infl(pwas_p_raw, df = 1)
  lambda_wo_reg[n] = pval_infl(pwas_p_wo_reg, df = 1)
  lambda_w_reg[n] = pval_infl(pwas_p_w_reg, df = 1)
  lambda_wo_reg_center[n] = pval_infl(pwas_p_wo_reg_center, df = 1)
  lambda_w_reg_center[n] = pval_infl(pwas_p_w_reg_center, df = 1)
  
  fpr_raw[n] = mean(pwas_p_raw < 0.05)
  fpr_wo_reg[n] = mean(pwas_p_wo_reg < 0.05)
  fpr_w_reg[n] = mean(pwas_p_w_reg < 0.05)
  fpr_wo_reg_center[n] = mean(pwas_p_wo_reg_center < 0.05)
  fpr_w_reg_center[n] = mean(pwas_p_w_reg_center < 0.05)
  print(n)
}

lambda_all = list(lambda_wo_reg, 
                  lambda_w_reg,
                  lambda_wo_reg_center,
                  lambda_w_reg_center,
                  lambda_raw,
                  fpr_wo_reg,
                  fpr_w_reg,
                  fpr_wo_reg_center,
                  fpr_w_reg_center,
                  fpr_raw)
saveRDS(lambda_all, 'cis_trans_inter/11_inter_pheno_asso/lambda_sim.rds')

lambda_all = readRDS('cis_trans_inter/11_inter_pheno_asso/lambda_sim.rds')
lambda_plt = data.frame(lambda = c(lambda_all[[2]], lambda_all[[1]], lambda_all[[3]]),
                        group = rep(c('Adding trans additive effects',
                                      'Excluding trans additive effects without centering',
                                      'Excluding trans additive effects after centering'), each = 100))
lambda_plt$group = factor(lambda_plt$group,
                          levels = c('Adding trans additive effects',
                                     'Excluding trans additive effects without centering',
                                     'Excluding trans additive effects after centering'))
library(ggbreak)
ggplot(lambda_plt, aes(x = group, y = lambda, fill = group)) +
  geom_boxplot(width = 0.5) +
  labs(x = "", y = expression(lambda),
       title = "Centering removes inflation") +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_break(breaks = c(1.5, 151)) +
  scale_y_continuous(limits = c(0.8, 151.2)) +
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(1,3,2)]) +
  theme(text = element_text(size=14, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13,
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

fdr_plt = data.frame(fdr = c(lambda_all[[6]], lambda_all[[5]], lambda_all[[7]]),
                     group = rep(c('Adding trans additive effects',
                                   'Excluding trans additive effects without centering',
                                   'Excluding trans additive effects after centering'), each = 100))
fdr_plt$group = factor(fdr_plt$group,
                       levels = c('Adding trans additive effects',
                                  'Excluding trans additive effects without centering',
                                  'Excluding trans additive effects after centering'))
ggplot(fdr_plt, aes(x = group, y = fdr, fill = group)) +
  geom_boxplot(width = 0.5) +
  labs(x = "", y = 'False dicovery rate',
       title = "Centering controls false dicovery rate") +
  geom_hline(yintercept = 0.05, linetype = 2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(1,3,2)]) +
  theme(text = element_text(size=14, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13,
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

###### when the true effect size of target is not zero
nsim = 1000
N = 20000 # number of individuals
h2_reg = 0.2 # regulator h2
h2_target_cis = 0.1 # target cis var is 0.1 of the total genetic var
h2_target_trans = 0.8 # target trans var is 0.8 of the total genetic var
h2_target_inter = 0.1 # target interaction var is 0.1 of the total genetic var
res_coef = 4 # residual var is 4 times genetic var
train_set = 1:5000
test_set = 5001:10000
pwas_set = 10001:20000

pwas_coef_raw = c()
pwas_coef_wo_reg = c()
pwas_p_wo_reg = c()
pwas_coef_w_reg = c()
pwas_p_w_reg = c()

pwas_p_raw = c()
pwas_coef_wo_reg_center = c()
pwas_p_wo_reg_center = c()
pwas_coef_w_reg_center = c()
pwas_p_w_reg_center = c()

for (i in 1:nsim) {
  reg_snp1 = rbinom(N, 2, 0.1) # regulator cis SNP1
  reg_snp2 = rbinom(N, 2, 0.3) # regulator cis SNP2
  reg_snp3 = rbinom(N, 2, 0.5) # regulator cis SNP3
  
  target_snp1 = rbinom(N, 2, 0.1) # target cis SNP1
  target_snp2 = rbinom(N, 2, 0.3) # target cis SNP2
  target_snp3 = rbinom(N, 2, 0.5) # target cis SNP3
  
  beta_reg = rnorm(3)
  reg_cis = reg_snp1*beta_reg[1] + reg_snp2*beta_reg[2] + reg_snp3*beta_reg[3]
  var_reg_e = (1 - h2_reg)/h2_reg * var(reg_cis)
  reg = reg_cis + rnorm(N, 0, sqrt(var_reg_e))
  
  beta_target = rnorm(3)
  target_cis = target_snp1*beta_target[1] + target_snp2*beta_target[2] +
    target_snp3*beta_target[3]
  
  target_cis_var = var(target_cis)
  add_eff = sqrt(target_cis_var * (h2_target_trans/h2_target_cis) / var(reg))  # regulator effect on target
  inter_eff = sqrt(target_cis_var * (h2_target_inter/h2_target_cis) / var(reg * target_cis)) # interaction effect
  var_target_e = res_coef * var(target_cis + add_eff * reg + inter_eff * reg * target_cis)
  target = target_cis + add_eff * reg +
    inter_eff * reg * target_cis + rnorm(N, 0, sqrt(var_target_e)) # target expression
  ## normalize expression
  reg = (reg - mean(reg))/sd(reg)
  target = (target - mean(target))/sd(target)
  
  trait_eff = 5 # regulator effect on the phenotype. In real data, slope of ALP on ABO expression is -5
  target_eff = 2 # target effect on the phenotype.
  Y = trait_eff * reg + target_eff * target + rnorm(N,0,20) # Y is the phenotype. In real data, the regression residual var is 26
  Y_center = Y - mean(Y)
  
  ## training testing split
  reg_train = reg[train_set]
  reg_snp1_train = reg_snp1[train_set]
  reg_snp2_train = reg_snp2[train_set]
  reg_snp3_train = reg_snp3[train_set]
  
  reg_snp1_test = reg_snp1[test_set]
  reg_snp2_test = reg_snp2[test_set]
  reg_snp3_test = reg_snp3[test_set]
  
  reg_snp1_pwas = reg_snp1[pwas_set]
  reg_snp2_pwas = reg_snp2[pwas_set]
  reg_snp3_pwas = reg_snp3[pwas_set]
  
  target_train = target[train_set]
  target_snp1_train = target_snp1[train_set]
  target_snp2_train = target_snp2[train_set]
  target_snp3_train = target_snp3[train_set]
  
  target_test = target[test_set]
  target_snp1_test = target_snp1[test_set]
  target_snp2_test = target_snp2[test_set]
  target_snp3_test = target_snp3[test_set]
  
  target_snp1_pwas = target_snp1[pwas_set]
  target_snp2_pwas = target_snp2[pwas_set]
  target_snp3_pwas = target_snp3[pwas_set]
  
  ## predicting cis-regulatory component
  fit_reg = lm(reg_train ~ reg_snp1_train + reg_snp2_train + reg_snp3_train)
  reg_coef = summary(fit_reg)$coefficients
  
  reg_hat = reg_coef[1,1] + reg_coef[2,1] * reg_snp1_test + reg_coef[3,1] * reg_snp2_test +
    reg_coef[4,1] * reg_snp3_test
  reg_hat_center = reg_hat - mean(reg_hat)
  
  reg_hat_pwas = reg_coef[1,1] + reg_coef[2,1] * reg_snp1_pwas + reg_coef[3,1] * reg_snp2_pwas +
    reg_coef[4,1] * reg_snp3_pwas
  reg_hat_pwas_center = reg_hat_pwas - mean(reg_hat_pwas)
  
  fit_target = lm(target_train ~ target_snp1_train + target_snp2_train + target_snp3_train)
  target_coef = summary(fit_target)$coefficients
  
  target_hat = target_coef[1,1] + target_coef[2,1] * target_snp1_test + target_coef[3,1] * target_snp2_test +
    target_coef[4,1] * target_snp3_test
  target_hat_center = target_hat - mean(target_hat)
  
  target_hat_pwas = target_coef[1,1] + target_coef[2,1] * target_snp1_pwas + target_coef[3,1] * target_snp2_pwas +
    target_coef[4,1] * target_snp3_pwas
  target_hat_pwas_center = target_hat_pwas - mean(target_hat_pwas)
  
  ## predicting target by incorporating epistasis
  fit_target_epi = lm(target_test ~ target_hat * reg_hat)
  fit_target_epi_center = lm(target_test ~ target_hat_center * reg_hat_center)
  target_epi_coef = summary(fit_target_epi)$coefficients
  target_epi_center_coef = summary(fit_target_epi_center)$coefficients
  # without regulator additive part
  target_pwas1 = target_epi_coef[1,1] + target_epi_coef[2,1] * target_hat_pwas +
    target_epi_coef[4,1] * target_hat_pwas * reg_hat_pwas
  target_pwas1_center = target_epi_center_coef[1,1] + target_epi_center_coef[2,1] * target_hat_pwas_center +
    target_epi_center_coef[4,1] * target_hat_pwas_center * reg_hat_pwas_center
  
  # with regulator additive part
  target_pwas2 = target_epi_coef[1,1] + target_epi_coef[2,1] * target_hat_pwas +
    target_epi_coef[3,1] * reg_hat_pwas + target_epi_coef[4,1] * target_hat_pwas * reg_hat_pwas
  target_pwas2_center = target_epi_center_coef[1,1] + target_epi_center_coef[2,1] * target_hat_pwas_center +
    target_epi_center_coef[3,1] * reg_hat_pwas_center +
    target_epi_center_coef[4,1] * target_hat_pwas_center * reg_hat_pwas_center
  
  ## PWAS of two models
  Y_pwas = Y[pwas_set]
  
  pwas0 = lm(Y_pwas ~ target_hat_pwas)
  pwas1 = lm(Y_pwas ~ target_pwas1)
  pwas2 = lm(Y_pwas ~ target_pwas2)
  
  pwas1_center = lm(Y_pwas ~ target_pwas1_center)
  pwas2_center = lm(Y_pwas ~ target_pwas2_center)
  
  pwas_coef_raw[i] = summary(pwas0)$coefficients[2,1]
  pwas_p_raw[i] = summary(pwas0)$coefficients[2,4]
  
  pwas_coef_wo_reg[i] = summary(pwas1)$coefficients[2,1]
  pwas_p_wo_reg[i] = summary(pwas1)$coefficients[2,4]
  
  pwas_coef_w_reg[i] = summary(pwas2)$coefficients[2,1]
  pwas_p_w_reg[i] = summary(pwas2)$coefficients[2,4]
  
  pwas_coef_wo_reg_center[i] = summary(pwas1_center)$coefficients[2,1]
  pwas_p_wo_reg_center[i] = summary(pwas1_center)$coefficients[2,4]
  
  pwas_coef_w_reg_center[i] = summary(pwas2_center)$coefficients[2,1]
  pwas_p_w_reg_center[i] = summary(pwas2_center)$coefficients[2,4]
}

coef_plt2 = data.frame(coef = c(pwas_coef_raw, pwas_coef_wo_reg, pwas_coef_wo_reg_center,
                                pwas_coef_w_reg_center),
                       group = rep(c('raw', 'Excluding trans additive effects without centering',
                                     'Excluding trans additive effects after centering',
                                     'Adding trans additive effects'),
                                   each = nsim))
coef_plt2$group = factor(coef_plt2$group,
                         levels = c('raw', 'Adding trans additive effects',
                                    'Excluding trans additive effects without centering',
                                    'Excluding trans additive effects after centering'))
ggplot(coef_plt2, aes(x = group, y = coef, fill = group)) +
  geom_boxplot(width = 0.5) +
  labs(x = "", y = "PWAS coefficient estimate",
       title = "Centering is unbiased") +
  geom_hline(yintercept = 2, linetype = 2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(1,3,2,4)]) +
  theme(text = element_text(size=14, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13,
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


coef_plt2 = data.frame(coef = c(pwas_coef_raw, pwas_coef_wo_reg, pwas_coef_wo_reg_center,
                                pwas_coef_w_reg_center),
                       group = rep(c('raw', 'Excluding trans additive effects without centering',
                                     'Excluding trans additive effects after centering',
                                     'Adding trans additive effects'),
                                   each = nsim))
coef_plt2$group = factor(coef_plt2$group,
                         levels = c('raw', 'Adding trans additive effects',
                                    'Excluding trans additive effects without centering',
                                    'Excluding trans additive effects after centering'))
ggplot(coef_plt2, aes(x = group, y = coef, fill = group)) +
  geom_boxplot(width = 0.5) +
  labs(x = "", y = "PWAS coefficient estimate",
       title = "Adding trans additive effects causes severe PWAS bias") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(1,3,2,4)]) +
  theme(text = element_text(size=14, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13,
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
# 
# p_infl1 = pval_infl(pwas_p_wo_reg_center, df = 1)
# p_infl2 = pval_infl(pwas_p_wo_reg, df = 1)
# p_infl3 = pval_infl(pwas_p_w_reg, df = 1)
# 
# qqplot.pvalues(pwas_p_wo_reg_center, col.abline = "black", pch = 19, 
#                col.CB = "gray80", col = brewer.pal(8, "Set2")[2], 
#                main = expression('Q-Q plot of simulated PWAS '*italic(P)),
#                xlab = expression('Expected -log'[10]*'('*italic(P)*')'),
#                ylab = expression('Observed -log'[10]*'('*italic(P)*')'))
# qq_plot(pwas_p_wo_reg, col = brewer.pal(8, "Set2")[3], pch = 19)
# qq_plot(pwas_p_w_reg, col = brewer.pal(8, "Set2")[1], pch = 19)





