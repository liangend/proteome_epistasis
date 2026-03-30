setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
h2_tab = fread('cis_trans_inter/08_prot_h2/cis_h2_summ.txt')
h2_tab_sub = h2_tab[h2_tab$R2 > 0.01, ]
h2_tab_sub$r2_cat = cut(h2_tab_sub$R2, 
                        breaks = quantile(h2_tab_sub$R2, seq(0,1,0.2)),
                        include.lowest = T)
h2_tab_sub$inter_hot = h2_tab_sub$n_inter >= 10
h2_tab_sub$target_hot = h2_tab_sub$n_target >= 10

ukb_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')
ukb_cov = ukb_cov[, 2:3]
colnames(ukb_cov) = c('id', 'sex')

## proteins that are target hotspot
h2_hot = h2_tab_sub[h2_tab_sub$target_hot, ]
h2_hot = h2_hot[order(h2_hot$n_target), ]
par(mfrow = c(3,3))
prot_hot_sample = c(seq(1, nrow(h2_hot), 3)[1:7], 20, 22)
for (i in prot_hot_sample) {
  pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/', h2_hot$prot[i], '_good_pred.txt'))
  pred_i$sex = ukb_cov$sex[match(pred_i$id, ukb_cov$id)]
  plot(pred_i$pred, pred_i$obs, xlab = 'pred', ylab = 'obs',
       main = paste0(h2_hot$prot[i], ', R2 = ', round(h2_hot$R2[i],2), ', n target = ', h2_hot$n_target[i]))
}
ggplot(pred_i, aes(x=pred, y=obs, color=as.factor(sex))) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  labs(x = "pred", y = 'obs', title = '', color = 'sex') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## proteins that are not target hotspot
non_hot_sample = c(sample(h2_tab_sub$prot[!h2_tab_sub$target_hot & h2_tab_sub$R2 < 0.1], 3),
                   sample(h2_tab_sub$prot[!h2_tab_sub$target_hot & h2_tab_sub$R2 > 0.2 & h2_tab_sub$R2 < 0.5], 3),
                   sample(h2_tab_sub$prot[!h2_tab_sub$target_hot & h2_tab_sub$R2 > 0.7], 3))
for (i in non_hot_sample) {
  pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/', i, '_good_pred.txt'))
  plot(pred_i$pred, pred_i$obs, xlab = 'pred', ylab = 'obs', 
       main = paste0(i, ', R2 = ', round(h2_tab_sub$R2[which(h2_tab_sub$prot == i)], 2)))
}








