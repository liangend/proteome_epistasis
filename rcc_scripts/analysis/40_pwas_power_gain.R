setwd('/project/xuanyao/jinghui')
library(ggplot2)
library(RColorBrewer)

###### when the true effect size of target is not zero
nsim = 1000
N = 20000 # number of individuals
h2_reg = 0.2 # regulator h2

res_coef = 4 # residual var is 4 times genetic var
train_set = 1:5000
test_set = 5001:10000
pwas_set = 10001:20000

h2_target_cis = 0.1 # target cis var is 0.1 of the total genetic var
h2_target_trans_all = c(0.8, 0.6, 0.4, 0.2) # target trans var is 0.8 of the total genetic var
h2_target_inter_all = c(0.1, 0.3, 0.5, 0.7) # target interaction var is 0.1 of the total genetic var

power_raw = c()
power_epi_adj = c()

for (n in 1:4) {
  h2_target_trans = h2_target_trans_all[n]
  h2_target_inter = h2_target_inter_all[n]
  
  pwas_coef_raw = c()
  pwas_p_raw = c()
  pwas_coef_wo_reg_center = c()
  pwas_p_wo_reg_center = c()
  
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
    target_pwas1_center = target_epi_center_coef[1,1] + target_epi_center_coef[2,1] * target_hat_pwas_center +
      target_epi_center_coef[4,1] * target_hat_pwas_center * reg_hat_pwas_center
    
    ## PWAS of two models
    Y_pwas = Y[pwas_set]
    
    pwas0 = lm(Y_pwas ~ target_hat_pwas)
    pwas1_center = lm(Y_pwas ~ target_pwas1_center)
    
    pwas_coef_raw[i] = summary(pwas0)$coefficients[2,1]
    pwas_p_raw[i] = summary(pwas0)$coefficients[2,4]
    
    pwas_coef_wo_reg_center[i] = summary(pwas1_center)$coefficients[2,1]
    pwas_p_wo_reg_center[i] = summary(pwas1_center)$coefficients[2,4]
  }
  power_raw[n] = sum(pwas_p_raw < 0.05)
  power_epi_adj[n] = sum(pwas_p_wo_reg_center < 0.05)
}


power_plt = data.frame(power = c(power_raw, power_epi_adj)/nsim,
                       group = rep(c('Raw', 'Epistasis adjusted'),
                                   each = 4),
                       epi_h2 = c(0.1, 0.3, 0.5, 0.7))
ggplot(power_plt, aes(x = as.factor(epi_h2), y = power, fill = group)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  labs(x = expression("Proportion of h"[g]^2 *" explained by epistasis"), y = "Power",
       title = "PWAS power gain by including epistasis effects", fill = '') +
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(2,1)]) +
  theme(text = element_text(size=14, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


