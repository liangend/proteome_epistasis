setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(RColorBrewer)
inter_sig = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')

# cor_cis = c()
# cor_add = c()
# cor_inter = c()
# cor_inter_wo_add = c()
# n_sample = c()
# for (n in 1:nrow(inter_sig)) {
#   target_prot = inter_sig$target[n]
#   reg_prot = inter_sig$regulator[n]
#   target = fread(paste0('cis_trans_inter/06_fusion_pred/eur/',
#                         target_prot, '_good_pred.txt'))
#   reg = fread(paste0('cis_trans_inter/06_fusion_pred/eur/',
#                      reg_prot, '_good_pred.txt'))
#   target$reg = reg$pred[match(target$id, reg$id)]
#   target = na.omit(target)
#   n_sample[n] = nrow(target)
#   print(paste0('sample size = ', nrow(target)))
#   
#   target$pred = target$pred - mean(target$pred)
#   target$reg = target$reg - mean(target$reg)
#   
#   ### residualize target expression
#   ukb_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')
#   ukb_cov = na.omit(ukb_cov)
#   ukb_cov_sub = ukb_cov[match(target$id, ukb_cov$V1), -(1:2)]
#   target$target_res = residuals(lm(target$obs ~ as.matrix(ukb_cov_sub)))
#   
#   ### k-fold cross-validation
#   k = 10
#   fold_size = ceiling(nrow(target) / k)
#   pred_add = c()
#   pred_inter = c()
#   pred_inter_wo_add = c()
#   for (i in 1:k) {
#     ## training testing split
#     if (i < k) {
#       test_set = (1 + (i-1)*fold_size):(i*fold_size)
#     } else {
#       test_set = (1 + (i-1)*fold_size):nrow(target)
#     }
#     train_set = setdiff(1:nrow(target), test_set)
#     train_dat = target[train_set, ]
#     test_dat = target[test_set, ]
#     
#     ## model fitting
#     model_add = lm(target_res ~ pred + reg, data = train_dat)
#     model_inter = lm(target_res ~ pred * reg, data = train_dat)
#     coef_add = summary(model_add)$coefficients
#     coef_inter = summary(model_inter)$coefficients
#     
#     ## model prediction
#     pred_add_i = coef_add[1,1] + coef_add[2,1]*test_dat$pred + 
#       coef_add[3,1]*test_dat$reg
#     pred_inter_i = coef_inter[1,1] + coef_inter[2,1]*test_dat$pred + 
#       coef_inter[3,1]*test_dat$reg + coef_inter[4,1]*test_dat$pred*test_dat$reg
#     pred_inter_wo_add_i = coef_inter[1,1] + coef_inter[2,1]*test_dat$pred + 
#       coef_inter[4,1]*test_dat$pred*test_dat$reg
#     
#     pred_add = c(pred_add, pred_add_i)
#     pred_inter = c(pred_inter, pred_inter_i)
#     pred_inter_wo_add = c(pred_inter_wo_add, pred_inter_wo_add_i)
#   }
#   
#   cor_cis[n] = cor(target$pred, target$target_res)
#   cor_add[n] = cor(pred_add, target$target_res)
#   cor_inter[n] = cor(pred_inter, target$target_res)
#   cor_inter_wo_add[n] = cor(pred_inter_wo_add, target$target_res)
# }
# 
# pred_tab = data.frame(pred_cor = c(cor_cis, cor_add, cor_inter, cor_inter_wo_add),
#                       target = rep(inter_sig$target, 4), 
#                       regulator = rep(inter_sig$regulator, 4),
#                       group = rep(c('cis', 'target+reg', 'target*reg', 'target+interaction'),
#                                   each = nrow(inter_sig)))
# pred_tab$inter_pair = paste0(pred_tab$target, '*', 
#                              pred_tab$regulator)
# pred_tab$r2 = pred_tab$pred_cor^2
# fwrite(pred_tab, 'cis_trans_inter/06_fusion_pred/acc_cv.txt', sep = '\t')

pred_tab = fread('cis_trans_inter/06_fusion_pred/acc_cv.txt')
pred_tab = pred_tab[pred_tab$group %in% c('cis', 'target+interaction'), ]
pred_tab$group[pred_tab$group == 'cis'] = 'cis only'
pred_tab$group[pred_tab$group == 'target+interaction'] = 'cis + interaction'
ggplot(pred_tab, aes(x=reorder(inter_pair, -r2), y=r2, fill=group)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  labs(x = '', y = 'pred R2', title = '') + 
  #coord_cartesian(ylim = c(0.2, 0.3)) + 
  scale_fill_manual(values = brewer.pal(8, "Set1")) +
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

