setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(boot)
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')

alk_pho = fread('cis_trans_inter/data/Alkaline_phosphatase.pheno')
alk_pho = na.omit(alk_pho)
colnames(alk_pho)[3] = 'pheno'
alk_pho$pheno = (alk_pho$pheno - mean(alk_pho$pheno)) / sd(alk_pho$pheno)

### ALP association with interaction terms alone
ld_indep_inter = readRDS('cis_trans_inter/17_snp_snp_inter/inter_w_add_indep.rds')
ukb_cov = fread('cis_trans_inter/06_fusion_pred/all_eur/all_norm_cov.txt')
ukb_cov = na.omit(ukb_cov)

names(ld_indep_inter)

alp_p = c()
for (inter_i in 1:length(ld_indep_inter)) {
  inter_sig_i = ld_indep_inter[[inter_i]]
  target_i = unlist(strsplit(names(ld_indep_inter)[inter_i], '_'))[1]
  reg_i = unlist(strsplit(names(ld_indep_inter)[inter_i], '_'))[2]
  
  target_prot = fread(paste0('cis_trans_inter/06_fusion_pred/all_eur/',
                             target_i, '_pred.txt'))
  reg_prot = fread(paste0('cis_trans_inter/06_fusion_pred/all_eur/',
                          reg_i, '_pred.txt'))
  
  common_id = Reduce(intersect, list(alk_pho$FID, target_prot$id, reg_prot$id,
                                     ukb_cov$FID))
  
  alk_pho_sub = alk_pho$pheno[match(common_id, alk_pho$FID)]
  ukb_cov_sub = ukb_cov[match(common_id, ukb_cov$FID), -(1:2)]
  alk_pho_res = residuals(lm(alk_pho_sub ~ as.matrix(ukb_cov_sub)))
  
  target_prot_sub = target_prot$pred[match(common_id, target_prot$id)]
  reg_prot_sub = reg_prot$pred[match(common_id, reg_prot$id)]
  
  inter_sub = target_prot_sub * reg_prot_sub
  fit_i = lm(alk_pho_res ~ inter_sub)
  coef_i = summary(fit_i)
  alp_p[inter_i] = coef_i$coefficients[2,4]
}
names(alp_p) = names(ld_indep_inter)

##### ALP and blood types
# alk_pho$pheno = (alk_pho$pheno - mean(alk_pho$pheno)) / 
#   sd(alk_pho$pheno)
abo_blood = fread('cis_trans_inter/data/ukb_blood_type.pheno')
abo_blood = na.omit(abo_blood)

target_prot = 'CD209'
target = fread(paste0('cis_trans_inter/06_fusion_pred/all_eur/',
                      target_prot, '_pred.txt'))
abo = fread('cis_trans_inter/06_fusion_pred/all_eur/ABO_pred.txt')
ukb_cov = fread('cis_trans_inter/06_fusion_pred/all_eur/all_norm_cov.txt')
ukb_cov = na.omit(ukb_cov)

all_coef = readRDS('cis_trans_inter/15_pheno_asso_pred_w_inter/inter_coef.rds')
target_coef = all_coef$`CD209*ABO`

common_id = Reduce(intersect, list(alk_pho$FID, abo_blood$FID,
                                   target$id, ukb_cov$FID))
alk_pho_sub = alk_pho$pheno[match(common_id, alk_pho$FID)]
abo_blood_sub = abo_blood$`23165`[match(common_id, abo_blood$FID)]
target_sub = target$pred[match(common_id, target$id)]
ukb_cov_sub = ukb_cov[match(common_id, ukb_cov$FID), -(1:2)]

## effect of target on ALP by blood type
alk_pho_res = residuals(lm(alk_pho_sub ~ as.matrix(ukb_cov_sub)))
# AA_ind = which(abo_blood_sub == 'AA')
# AB_ind = which(abo_blood_sub == 'AB')
# AO_ind = which(abo_blood_sub == 'AO')
# BB_ind = which(abo_blood_sub == 'BB')
# BO_ind = which(abo_blood_sub == 'BO')
# OO_ind = which(abo_blood_sub == 'OO')

reg_tab = data.frame(alk_pho_res = alk_pho_res, target = target_sub, 
                     abo_type = abo_blood_sub)
reg_tab$abo_pred = abo$pred[match(common_id, abo$id)]
reg_tab$blood_type = 'O'
reg_tab$blood_type[which(reg_tab$abo_type %in% c('AA', 'AO'))] = 'A'
reg_tab$blood_type[which(reg_tab$abo_type %in% c('BB', 'BO'))] = 'B'
reg_tab$blood_type[which(reg_tab$abo_type %in% c('AB'))] = 'AB'

reg_tab$o_type = 'OO'
reg_tab$o_type[which(reg_tab$abo_type %in% c('BO', 'AO'))] = 'OX'
reg_tab$o_type[which(reg_tab$abo_type %in% c('AA', 'BB', 'AB'))] = 'XX'

reg_tab$is_o = 'O'
reg_tab$is_o[which(reg_tab$blood_type != 'O')] = 'not O'

reg_tab$pred_target = target_coef[1,1] + target_coef[2,1] * (reg_tab$target - mean(reg_tab$target))+
  target_coef[3,1] * (reg_tab$abo_pred - mean(reg_tab$abo_pred)) * (reg_tab$target - mean(reg_tab$target))

ggplot(reg_tab, aes(x=target, y=alk_pho_res, color = is_o)) + 
  geom_smooth(method='lm') + 
  #geom_point(alpha = 0.1) + 
  # xlim(c(-0.8,0.1)) +
  # coord_cartesian(ylim = c(-7, 6)) + 
  labs(x = paste0("Predicted ", target_prot, ' with interaction'), y = 'ALP', 
       title = paste0('Effect of ', target_prot, ' on ALP'), color = '') + 
  scale_color_manual(values = brewer.pal(8, "Set1")) +
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

fit1 = lm(alk_pho_res ~ target * is_o, data = reg_tab)
summary(fit1)


## number of individuals for each blood type
abo_n_tab = as.data.frame(table(abo_blood_sub))
ggplot(abo_n_tab, aes(x = reorder(abo_blood_sub, -Freq), y = Freq, 
                      label = Freq)) + 
  geom_bar(stat="identity", width = 0.8, fill = 'steelblue') +
  geom_text(vjust=-0.2, size = 4) + 
  labs(title = '', x = "", y = "# ind", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# A, B, AB, O blood types
boot_lm = function(data, indices){
  df <- data[indices, ]
  aggre_coef = sapply(split(df, df$blood_type), 
                      function(sub){summary(lm(alk_pho_res ~ target, data = sub))$coefficients[2,1]})
  return(aggre_coef)
}

boot_lm_by_o = function(data, indices){
  df <- data[indices, ]
  aggre_coef = sapply(split(df, df$o_type), 
                      function(sub){summary(lm(alk_pho_res ~ target, data = sub))$coefficients[2,1]})
  return(aggre_coef)
}

# coef_boot = boot(reg_tab, boot_lm, 1000)
# coef_tab = data.frame(beta = coef_boot$t0,
#                       ci_low = apply(coef_boot$t, 2, function(x){quantile(x, 0.025)}),
#                       ci_up = apply(coef_boot$t, 2, function(x){quantile(x, 0.975)}),
#                       type = c('A', 'AB', 'B', 'O'))

## FCGR2B coef
coef_tab = data.frame(beta = c(0.037748127, 0.018416359, 0.005926087, 0.043973324),
                      ci_low = c(0.01672009, -0.04642254, -0.03473172, 0.02418331),
                      ci_up = c(0.05910919, 0.07866973, 0.05076496, 0.06394378),
                      type = c('A', 'AB', 'B', 'O'))

## CD209 coef
coef_tab = data.frame(beta = c(0.007052894, 0.065317592, -0.043172327, -0.035355251),
                      ci_low = c(-0.01896265, -0.02881639, -0.10713492, -0.06286292),
                      ci_up = c(0.032669180, 0.150386767, 0.014193287, -0.005833564),
                      type = c('A', 'AB', 'B', 'O'))

ggplot(coef_tab, aes(x = type, y = beta, color = type)) + 
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width=.1, 
                position=position_dodge()) + 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = paste0('Effect of ', target_prot, ' on alkaline phosphatase by blood types'), 
       x = "", y = "beta", color = '') +
  scale_color_manual(values = brewer.pal(8, "Dark2")) +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

# O and non-O blood types
# coef_boot_o = boot(reg_tab, boot_lm_by_o, 100)
# coef_tab_o = data.frame(beta = coef_boot_o$t0,
#                       ci_low = apply(coef_boot_o$t, 2, function(x){quantile(x, 0.025)}),
#                       ci_up = apply(coef_boot_o$t, 2, function(x){quantile(x, 0.975)}),
#                       type = c('not O', 'O'))

coef_tab_o = data.frame(beta = c(0.004050432, -0.035355251),
                        ci_low = c(-0.02087068, -0.06221709),
                        ci_up = c(0.02538119, -0.00116622),
                        type = c('not O', 'O'))

ggplot(coef_tab_o, aes(x = type, y = beta, color = type)) + 
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width=.1, 
                position=position_dodge()) + 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = paste0('Effect of ', target_prot, ' on alkaline phosphatase by blood types'), 
       x = "", y = "beta", color = '') +
  scale_color_manual(values = brewer.pal(8, "Dark2")) +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')



