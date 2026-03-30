setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gaston)
library(simtrait)
library(RColorBrewer)

qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

inter_test = fread('cis_trans_inter/07_prot_inter/prot_perm_test.txt')
uniq_target = unique(inter_test$target)

infl_add = c()
infl_inter = c()
is_hot = c()
for (i in 1:length(uniq_target)) {
  target_i = uniq_target[i]
  pval_i = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/sig_perm/', 
                   target_i, '.rds'))
  add_p_i = pval_i$reg_add_p
  inter_p_i = pval_i$inter_p
  
  ## calculate p value inflation
  add_infl_i = apply(add_p_i, 2, pval_infl)
  inter_infl_i = apply(inter_p_i, 2, pval_infl)
  is_hot_i = inter_test$is_hot_target[which(inter_test$target == target_i)]
  
  infl_add = c(infl_add, add_infl_i)
  infl_inter = c(infl_inter, inter_infl_i)
  is_hot = c(is_hot, is_hot_i)
}

infl_tab = data.frame(infl_add = infl_add, infl_inter = infl_inter, is_hot = is_hot, 
                      target = inter_test$target)
ggplot(infl_tab, aes(x= is_hot,  y = infl_inter, color = is_hot)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "is target hotspot", y = 'pval inflation', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

## number of target vs. p value inflation
infl_tab_sub = infl_tab[infl_tab$is_hot, ]
infl_tab_sub$n_target = rep(unname(table(infl_tab_sub$target)), times = unname(table(infl_tab_sub$target)))

ggplot(infl_tab_sub, aes(x= as.factor(n_target),  y = infl_inter)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "# targets", y = 'interaction pval inflation', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

ggplot(infl_tab_sub, aes(x= as.factor(n_target),  y = infl_add)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "# targets", y = 'additive pval inflation', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

# target_i = readRDS('cis_trans_inter/07_prot_inter/prot_prot_inter/sig_perm/ASAH2.rds')
# qqplot.pvalues(target_i$inter_p[,1], col.abline = "black", pch = 19)






