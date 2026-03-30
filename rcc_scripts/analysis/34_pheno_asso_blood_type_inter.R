setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(openxlsx)
library(RColorBrewer)
library(gridExtra)
all_inter_file = list.files('cis_trans_inter/19_pheno_asso_inter_otype/')
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter_w_cis.txt')

ukb_trait = read.xlsx('cis_trans_inter/00_ref/ukb_trait.xlsx', sheet = 3)

inter_p = c() # interaction p in the interaction model
for (i in all_inter_file) {
  coef_i = readRDS(paste0('cis_trans_inter/19_pheno_asso_inter_otype/', i))
  inter_p = rbind(inter_p, sapply(coef_i, function(x){min(x[5:6,6])}))
}
rownames(inter_p) = sub('.rds', '', all_inter_file)

### significant interaction associations
inter_p = inter_p[!rownames(inter_p) %in% c('Hypertension', 'Snoring'), ]
inter_q = qvalue(inter_p)
inter_q = inter_q$qvalues
sum(inter_q < 0.05)
inter_q = as.data.frame(inter_q)
inter_q$trait = rownames(inter_q)
inter_q$trait_group = ukb_trait$Group[match(inter_q$trait,
                                            ukb_trait$abbrev)]
inter_q = na.omit(inter_q)
inter_q_long = reshape(inter_q, varying = names(inter_q)[1:31],
                       v.names = "q_val", 
                       timevar = "interaction", 
                       times = names(inter_q)[1:31],
                       direction = "long")
inter_q_long = inter_q_long[, 1:4]
inter_q_long = inter_q_long[inter_q_long$q_val < 0.05, ]
sig_inter_summ = as.data.frame(table(inter_q_long$trait_group, 
                                     inter_q_long$interaction))
ggplot(sig_inter_summ, aes(x = reorder(Var2, -Freq), y = Freq, 
                           fill = Var1)) + 
  geom_bar(stat="identity") +
  labs(x = "", 
       y = "# significant associations", fill = 'trait') +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## target add p comparison: model with interactions vs. without interactions
add_p_in_inter = as.data.frame(add_p_in_inter)
add_p_in_inter$trait = rownames(add_p_in_inter)
add_p_in_inter$trait_group = ukb_trait$Group[match(add_p_in_inter$trait, 
                                                   ukb_trait$abbrev)]

add_p_wo_inter = as.data.frame(add_p_wo_inter)
add_p_plt = cbind(add_p_in_inter, add_p_wo_inter)
add_p_plt = na.omit(add_p_plt)

comp_plt = list()
for (i in 1:31) {
  plt_tab_i = add_p_plt[, c(i, 32, 33, i+33)]
  inter_i = colnames(plt_tab_i)[1]
  colnames(plt_tab_i)[c(1,4)] = c('w_inter', 'wo_inter')
  comp_plt[[i]] = ggplot(plt_tab_i, aes(x = -log10(w_inter), y = -log10(wo_inter), 
                                        color = trait_group)) + 
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, lty = 2) + 
    labs(x = "-log10(p) with inter", 
         y = "-log10(p) without inter", color = '', title = inter_i) +
    theme(text = element_text(size=13, colour = "black"), 
          axis.text.x = element_text(colour = "black", size = 13),
          axis.text.y = element_text(colour = "black", size = 13),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none')
}

do.call(grid.arrange, c(comp_plt, nrow = 7))


## target add p comparison: target expression prediction w/wo interactions
add_p_w_inter = as.data.frame(add_p_w_inter)
add_p_w_inter$trait = rownames(add_p_w_inter)
add_p_w_inter$trait_group = ukb_trait$Group[match(add_p_w_inter$trait, 
                                                  ukb_trait$abbrev)]
colnames(add_p_w_inter) = colnames(add_p_in_inter)

add_p_plt2 = cbind(add_p_w_inter, add_p_wo_inter)
add_p_plt2 = na.omit(add_p_plt2)


comp_plt2 = list()
for (i in 1:31) {
  plt_tab_i = add_p_plt2[, c(i, 32, 33, i+33)]
  inter_i = colnames(plt_tab_i)[1]
  colnames(plt_tab_i)[c(1,4)] = c('w_inter', 'wo_inter')
  comp_plt2[[i]] = ggplot(plt_tab_i, aes(x = -log10(w_inter), y = -log10(wo_inter), 
                                         color = trait_group)) + 
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, lty = 2) + 
    labs(x = "-log10(p) with inter", 
         y = "-log10(p) wo inter", color = '', title = inter_i) +
    theme(text = element_text(size=13, colour = "black"), 
          axis.text.x = element_text(colour = "black", size = 13),
          axis.text.y = element_text(colour = "black", size = 13),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none')
}
do.call(grid.arrange, c(comp_plt2, nrow = 7))

