setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(RColorBrewer)

alp = fread('cis_trans_inter/data/Alkaline_phosphatase.pheno')
alp = na.omit(alp)
colnames(alp)[3] = 'pheno'
alp$pheno = (alp$pheno - mean(alp$pheno)) / sd(alp$pheno)


target_prot = 'FCGR2B'
target = fread(paste0('cis_trans_inter/06_fusion_pred/all_eur/',
                      target_prot, '_pred.txt'))

abo = fread('cis_trans_inter/06_fusion_pred/all_eur/ABO_pred.txt')

ukb_cov = fread('cis_trans_inter/06_fusion_pred/all_eur/all_norm_cov.txt')
ukb_cov = na.omit(ukb_cov)

common_id = Reduce(intersect, list(alp$FID, abo$id,
                                   target$id, ukb_cov$FID))
alp_sub = alp$pheno[match(common_id, alp$FID)]
abo_sub = abo$pred[match(common_id, abo$id)]
target_sub = target$pred[match(common_id, target$id)]
ukb_cov_sub = ukb_cov[match(common_id, ukb_cov$FID), -(1:2)]

## effect of target on ALP by blood type
alp_res = residuals(lm(alp_sub ~ as.matrix(ukb_cov_sub)))

reg_tab = data.frame(alp = alp_res, target = (target_sub-mean(target_sub)), 
                     abo = (abo_sub - mean(abo_sub)))
fit1 = lm(alp ~ target * abo, data = reg_tab)
summary(fit1)

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

