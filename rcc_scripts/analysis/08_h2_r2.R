setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
### comparison of cis h2 and R2
# prediction accuracy
pred_all = list.files('cis_trans_inter/06_fusion_pred/eur/')
pred_prot = sapply(strsplit(pred_all, '_', fixed = T), '[', 1)
na_prop = c()
cor_all = c()
for (i in 1:length(pred_all)) {
  pred_i = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', pred_all[i]))
  na_prop[i] = sum(is.na(pred_i))/nrow(pred_i)
  cor_all[i] = cor(pred_i$pred, pred_i$obs, use = 'na.or.complete')
}
pred_tab = data.frame(prot = pred_prot, na_prop = na_prop,
                      cor = cor_all)
median(pred_tab$cor, na.rm = T)
median((pred_tab$cor)^2, na.rm = T)

# prot h2
all_h2_file = list.files('cis_trans_inter/08_prot_cis_h2/', pattern = '.hsq')
h2_all = c()
se_all = c()
for (i in all_h2_file) {
  h2_i = fread(paste0('cis_trans_inter/08_prot_cis_h2/', i), nrows = 4)
  h2_all[i] = h2_i$Variance[4]
  se_all[i] = h2_i$SE[4]
}

h2_tab = data.frame(prot = sub('.hsq', '', all_h2_file), 
                    h2 = h2_all, se = se_all)
h2_tab$pred_cor = cor_all[match(h2_tab$prot, pred_prot)]
h2_tab$na_prop = na_prop[match(h2_tab$prot, pred_prot)]
h2_tab$R2 = h2_tab$pred_cor ^ 2
h2_tab = h2_tab[order(h2_tab$h2), ]
h2_tab = h2_tab[h2_tab$se > 0, ]
median(h2_tab$h2)
# n_inter = as.data.frame(table(sig_p_tab$target))
# h2_tab$n_inter = n_inter$Freq[match(h2_tab$prot, n_inter$Var1)]
# fwrite(h2_tab, 'cis_trans_inter/08_prot_cis_h2/h2_summ.txt', sep = '\t')

## plot R2 and h2
h2_tab = fread('cis_trans_inter/08_prot_h2/cis_h2_500kb_summ.txt')

h2_tab = h2_tab[!is.na(h2_tab$h2) & !is.na(h2_tab$R2), ]
h2_tab$rank = rank(h2_tab$h2)
ggplot(h2_tab) + 
  geom_point(aes(x = rank, y = R2), size = 0.8, color = 'darkorange') +
  geom_point(aes(x = rank, y = h2), size = 0.8, color = 'black') +
  labs(title = expression(paste("Protein prediction accuracy close to ", 
                                italic("cis"), "-h"^2)), 
       x = expression(paste("Proteins ordered by ", italic("cis"), "-h"^2)), 
       y = expression(paste(italic("cis"), "-h"^2, " or R"^2))) +
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# r2_ave = tapply(h2_tab$R2, (seq_along(h2_tab$R2) - 1) %/% 5, function(x){mean(x, na.rm = T)})
# plot(1:nrow(h2_tab), h2_tab$h2, pch = 20, xlab = 'Protein ordered by cis-h2', ylab = 'cis-h2 or R2')
# points(1:nrow(h2_tab), h2_tab$R2, col = 'grey', cex = 0.5)
# points(1:length(r2_ave) * 5 - 2, r2_ave, col = 'red', cex = 1, pch = 20)
# legend('topleft', legend = c("h2", "R2", "Mean R2 for every 5 prot"),
#        col=c("black", "grey", "red"), pch=c(20,1,20), cex=0.8)
# 
# mean(h2_tab$R2 / h2_tab$h2, na.rm = T)










