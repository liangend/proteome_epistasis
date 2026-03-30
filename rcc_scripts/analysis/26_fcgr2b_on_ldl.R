setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

ldl_coef = readRDS('cis_trans_inter/11_inter_pheno_asso/asso_all_ind/LDL.rds')
ldl_add_coef = readRDS('cis_trans_inter/12_add_pheno_asso/asso_all_ind/LDL.rds')

ldl = fread('cis_trans_inter/data/LDL674178.pheno')
ldl = na.omit(ldl)
ldl$`30780` = (ldl$`30780` - mean(ldl$`30780`))/sd(ldl$`30780`)
abo_blood = fread('cis_trans_inter/data/ukb_blood_type.pheno')
abo_blood = na.omit(abo_blood)
fcgr2b = fread('cis_trans_inter/06_fusion_pred/all_eur/FCGR2B_pred.txt')
abo = fread('cis_trans_inter/06_fusion_pred/all_eur/ABO_pred.txt')
ukb_cov = fread('cis_trans_inter/06_fusion_pred/all_eur/all_norm_cov.txt')

common_id = Reduce(intersect, list(ldl$FID, abo$id,
                                   fcgr2b$id, ukb_cov$FID))
ldl_sub = ldl$`30780`[match(common_id, ldl$FID)]
abo_sub = abo$pred[match(common_id, abo$id)]
fcgr2b_sub = fcgr2b$pred[match(common_id, fcgr2b$id)]
ukb_cov_sub = ukb_cov[match(common_id, ukb_cov$FID), -(1:2)]
ldl_res = residuals(lm(ldl_sub ~ as.matrix(ukb_cov_sub)))
fit1 = lm(ldl_res ~ fcgr2b_sub)
summary(fit1)

fit2 = lm(ldl_res ~ fcgr2b_sub * abo_sub)
summary(fit2)

common_id = Reduce(intersect, list(ldl$FID, abo_blood$FID,
                                   fcgr2b$id, ukb_cov$FID))
ldl_sub = ldl$`30780`[match(common_id, ldl$FID)]
abo_blood_sub = abo_blood$`23165`[match(common_id, abo_blood$FID)]
fcgr2b_sub = fcgr2b$pred[match(common_id, fcgr2b$id)]
ukb_cov_sub = ukb_cov[match(common_id, ukb_cov$FID), -(1:2)]

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

## effect of FCGR2B on LDL by blood type
ldl_res = residuals(lm(ldl_sub ~ as.matrix(ukb_cov_sub)))
AA_ind = which(abo_blood_sub == 'AA')
AB_ind = which(abo_blood_sub == 'AB')
AO_ind = which(abo_blood_sub == 'AO')
BB_ind = which(abo_blood_sub == 'BB')
BO_ind = which(abo_blood_sub == 'BO')
OO_ind = which(abo_blood_sub == 'OO')

# all blood types
fit_AA = lm(ldl_res[AA_ind] ~ fcgr2b_sub[AA_ind])
fit_AB = lm(ldl_res[AB_ind] ~ fcgr2b_sub[AB_ind])
fit_AO = lm(ldl_res[AO_ind] ~ fcgr2b_sub[AO_ind])
fit_BB = lm(ldl_res[BB_ind] ~ fcgr2b_sub[BB_ind])
fit_BO = lm(ldl_res[BO_ind] ~ fcgr2b_sub[BO_ind])
fit_OO = lm(ldl_res[OO_ind] ~ fcgr2b_sub[OO_ind])

coef_tab = data.frame(beta = c(summary(fit_AA)$coefficients[2,1],
                               summary(fit_AB)$coefficients[2,1],
                               summary(fit_AO)$coefficients[2,1],
                               summary(fit_BB)$coefficients[2,1],
                               summary(fit_BO)$coefficients[2,1],
                               summary(fit_OO)$coefficients[2,1]),
                      p = c(summary(fit_AA)$coefficients[2,4],
                            summary(fit_AB)$coefficients[2,4],
                            summary(fit_AO)$coefficients[2,4],
                            summary(fit_BB)$coefficients[2,4],
                            summary(fit_BO)$coefficients[2,4],
                            summary(fit_OO)$coefficients[2,4]),
                      type = c('AA', 'AB', 'AO', 'BB', 'BO', 'OO'))

ggplot(coef_tab, aes(x = type, y = beta, label = round(p, 3))) + 
  geom_bar(stat="identity", width = 0.8, fill = 'steelblue') +
  geom_text(vjust=-0.2, size = 4) + 
  labs(title = '', x = "", y = "beta", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# A, B, AB, O blood types
fit_A = lm(ldl_res[c(AA_ind, AO_ind)] ~ fcgr2b_sub[c(AA_ind, AO_ind)])
fit_B = lm(ldl_res[c(BB_ind, BO_ind)] ~ fcgr2b_sub[c(BB_ind, BO_ind)])
fit_AB = lm(ldl_res[AB_ind] ~ fcgr2b_sub[AB_ind])
fit_O = lm(ldl_res[OO_ind] ~ fcgr2b_sub[OO_ind])

coef_tab2 = data.frame(beta = c(summary(fit_A)$coefficients[2,1],
                               summary(fit_B)$coefficients[2,1],
                               summary(fit_AB)$coefficients[2,1],
                               summary(fit_O)$coefficients[2,1]),
                      p = c(summary(fit_A)$coefficients[2,4],
                            summary(fit_B)$coefficients[2,4],
                            summary(fit_AB)$coefficients[2,4],
                            summary(fit_O)$coefficients[2,4]),
                      type = c('A', 'B', 'AB', 'O'))

ggplot(coef_tab2, aes(x = type, y = beta, label = round(p, 3))) + 
  geom_bar(stat="identity", width = 0.8, fill = 'steelblue') +
  geom_text(vjust=-0.2, size = 4) + 
  labs(title = '', x = "", y = "beta", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# O and non-O blood types
fit_not_O = lm(ldl_res[-OO_ind] ~ fcgr2b_sub[-OO_ind])
fit_O = lm(ldl_res[OO_ind] ~ fcgr2b_sub[OO_ind])

coef_tab3 = data.frame(beta = c(summary(fit_not_O)$coefficients[2,1],
                                summary(fit_O)$coefficients[2,1]),
                       p = c(summary(fit_not_O)$coefficients[2,4],
                             summary(fit_O)$coefficients[2,4]),
                       type = c('not O', 'O'))

ggplot(coef_tab3, aes(x = type, y = beta, label = round(p, 3))) + 
  geom_bar(stat="identity", width = 0.8, fill = 'steelblue') +
  geom_text(vjust=-0.2, size = 4) + 
  labs(title = '', x = "", y = "beta", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
