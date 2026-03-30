library(data.table)
library(boot)
library(ggplot2)
setwd('/project/xuanyao/jinghui')

cis_h2_prot = list.files('cis_trans_inter/08_prot_h2/cis_h2_1mb')
trans_h2_prot = list.files('cis_trans_inter/08_prot_h2/trans_h2_1mb/')
whole_h2_prot = list.files('cis_trans_inter/08_prot_h2/whole_h2/')

shared_prot = Reduce(intersect, list(cis_h2_prot, trans_h2_prot, whole_h2_prot))

cis_h2 = c()
cis_h2_se = c()
trans_h2 = c()
trans_h2_se = c()
whole_h2 = c()
whole_h2_se = c()
for (i in shared_prot) {
  cis_i = fread(paste0('cis_trans_inter/08_prot_h2/cis_h2_1mb/', i), nrows = 4)
  cis_h2[i] = cis_i$Variance[4]
  cis_h2_se[i] = cis_i$SE[4]
  
  trans_i = fread(paste0('cis_trans_inter/08_prot_h2/trans_h2_1mb/', i), nrows = 4)
  trans_h2[i] = trans_i$Variance[4]
  trans_h2_se[i] = trans_i$SE[4]
  
  whole_i = fread(paste0('cis_trans_inter/08_prot_h2/whole_h2/', i), nrows = 4)
  whole_h2[i] = whole_i$Variance[4]
  whole_h2_se[i] = whole_i$SE[4]
}

h2_tab = data.frame(cis_h2 = cis_h2,
                    cis_h2_se = cis_h2_se,
                    trans_h2 = trans_h2,
                    trans_h2_se = trans_h2_se,
                    whole_h2 = whole_h2,
                    whole_h2_se = whole_h2_se)

median(h2_tab$trans_h2) / 
  (median(h2_tab$trans_h2) + median(h2_tab$cis_h2))
mean(h2_tab$trans_h2) / 
  (mean(h2_tab$trans_h2) + mean(h2_tab$cis_h2))

ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
h2_ldsc = fread('pqtl/05_h2/04_h2_summ/prot_h2_ukb_1mb_5mb.txt')
h2_ldsc$gene_name = ukb_prot$gene_name[match(h2_ldsc$file, ukb_prot$file)]

shared_prot_name = sub('.hsq', '', shared_prot)
h2_ldsc = h2_ldsc[h2_ldsc$gene_name %in% shared_prot_name, ]

par(mfrow = c(2,2))
hist(h2_ldsc$cis_1mb_h2, breaks = 100, main = 'cis h2 LDSC', xlab = 'h2')
abline(v = mean(h2_ldsc$cis_1mb_h2), col = 'red', lwd = 2)
abline(v = median(h2_ldsc$cis_1mb_h2), col = 'blue', lwd = 2)

hist(h2_tab$cis_h2, breaks = 100, main = 'cis h2 REML', xlab = 'h2')
abline(v = mean(h2_tab$cis_h2), col = 'red', lwd = 2)
abline(v = median(h2_tab$cis_h2), col = 'blue', lwd = 2)

hist(h2_ldsc$trans_5mb_h2_ind, breaks = 100, main = 'trans h2 LDSC', xlab = 'h2')
abline(v = mean(h2_ldsc$trans_5mb_h2_ind), col = 'red', lwd = 2)
abline(v = median(h2_ldsc$trans_5mb_h2_ind), col = 'blue', lwd = 2)

hist(h2_tab$trans_h2, breaks = 100, main = 'trans h2 REML', xlab = 'h2')
abline(v = mean(h2_tab$trans_h2), col = 'red', lwd = 2)
abline(v = median(h2_tab$trans_h2), col = 'blue', lwd = 2)

ldsc_mean = mean(h2_ldsc$trans_5mb_h2_ind) / 
  (mean(h2_ldsc$trans_5mb_h2_ind) + mean(h2_ldsc$cis_1mb_h2))
ldsc_median1 = median(h2_ldsc$trans_5mb_h2_ind) / 
  (median(h2_ldsc$trans_5mb_h2_ind) + median(h2_ldsc$cis_1mb_h2))
ldsc_median2 = median(h2_ldsc$trans_5mb_h2_ind) / 
  (median(h2_ldsc$trans_5mb_h2_ind + h2_ldsc$cis_1mb_h2))

reml_mean = mean(h2_tab$trans_h2) / 
  (mean(h2_tab$trans_h2) + mean(h2_tab$cis_h2))
reml_median1 = median(h2_tab$trans_h2) / 
  (median(h2_tab$trans_h2) + median(h2_tab$cis_h2))
reml_median2 = median(h2_tab$trans_h2) / 
  (median(h2_tab$trans_h2 + h2_tab$cis_h2))

boot_prop_ldsc = function(data, i) {
  df = data[i, ]
  mean_prop = mean(df$trans_5mb_h2_ind)/
    mean(df$trans_5mb_h2_ind + df$cis_1mb_h2)
  median_prop1 = median(df$trans_5mb_h2_ind) / 
    (median(df$trans_5mb_h2_ind) + median(df$cis_1mb_h2))
  median_prop2 = median(df$trans_5mb_h2_ind) / 
    (median(df$trans_5mb_h2_ind + df$cis_1mb_h2))
  return(c(mean_prop, median_prop1, median_prop2))
}

boot_prop_reml = function(data, i) {
  df = data[i, ]
  mean_prop = mean(df$trans_h2)/
    mean(df$trans_h2 + df$cis_h2)
  median_prop1 = median(df$trans_h2) / 
    (median(df$trans_h2) + median(df$cis_h2))
  median_prop2 = median(df$trans_h2) / 
    (median(df$trans_h2 + df$cis_h2))
  return(c(mean_prop, median_prop1, median_prop2))
}
ldsc_boot = boot(h2_ldsc, boot_prop_ldsc, 1000)
reml_boot = boot(h2_tab, boot_prop_reml, 1000)

trans_prop_all = data.frame(trans_prop = c(ldsc_mean, ldsc_median1, 
                                           reml_mean, reml_median1),
                            lower = c(quantile(ldsc_boot$t[,1], 0.025), 
                                      quantile(ldsc_boot$t[,2], 0.025),
                                      quantile(reml_boot$t[,1], 0.025),
                                      quantile(reml_boot$t[,2], 0.025)), 
                            upper = c(quantile(ldsc_boot$t[,1], 0.975), 
                                      quantile(ldsc_boot$t[,2], 0.975),
                                      quantile(reml_boot$t[,1], 0.975),
                                      quantile(reml_boot$t[,2], 0.975)),
                            method = rep(c('LDSC', 'REML'), each = 2),
                            value = c(rep(c('mean', 'median') , 2)))
ggplot(trans_prop_all, aes(x= method,  y = trans_prop, 
                           color = value)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  geom_errorbar(position=position_dodge(0.3), aes(ymin=lower, ymax=upper), width=.1) +
  labs(x = "", y = bquote('trans/(cis + trans) h'^2), title = '', color = '') + 
  ylim(c(0.5,1.1)) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

