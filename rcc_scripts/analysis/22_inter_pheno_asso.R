setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(openxlsx)
library(RColorBrewer)
library(gridExtra)
all_reg_file = list.files('cis_trans_inter/12_add_pheno_asso/asso_all_ind/')
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter_w_cis.txt')
sig_inter0 = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')

ukb_trait = read.xlsx('cis_trans_inter/00_ref/ukb_trait.xlsx', sheet = 3)

add_p_in_inter = c()  # target additive p in the interaction model
reg_add_p = c() # regulator additive p in the interaction model
inter_p = c() # interaction p in the interaction model
n_sample = c()
add_p_w_inter = c() # target additive p when using interaction to predict target expression
add_p_wo_inter = c() # target additive p without interaction
for (i in all_reg_file) {
  ## interaction model
  coef_i = readRDS(paste0('cis_trans_inter/11_inter_pheno_asso/asso_all_ind/', i))
  add_p_in_inter = rbind(add_p_in_inter, sapply(coef_i[1:(length(coef_i)-1)], 
                                                function(x){x[2,4]}))
  reg_add_p = rbind(reg_add_p, sapply(coef_i[1:(length(coef_i)-1)], 
                                      function(x){x[3,4]}))
  inter_p = rbind(inter_p, sapply(coef_i[1:(length(coef_i)-1)], 
                                  function(x){x[4,6]}))
  n_sample = rbind(n_sample, coef_i$n_sample)
  
  ## additive model using interaction to predict target expression 
  coef_add_w_inter_i = readRDS(paste0('cis_trans_inter/15_pheno_asso_pred_w_inter/asso_2_add_inter/', i))
  add_p_w_inter = rbind(add_p_w_inter, sapply(coef_add_w_inter_i, function(x){x[2,4]}))
  
  ## additive model without interaction
  coef_add_i = readRDS(paste0('cis_trans_inter/12_add_pheno_asso/asso_all_ind/', i))
  add_p_wo_inter = rbind(add_p_wo_inter, sapply(coef_add_i, function(x){x[2,4]}))
}

rownames(add_p_in_inter) = sub('.rds', '', all_reg_file)
rownames(reg_add_p) = sub('.rds', '', all_reg_file)
rownames(inter_p) = sub('.rds', '', all_reg_file)

rownames(n_sample) = sub('.rds', '', all_reg_file)

rownames(add_p_w_inter) = sub('.rds', '', all_reg_file)

rownames(add_p_wo_inter) = sub('.rds', '', all_reg_file)

add_q = qvalue(add_p_wo_inter)
add_q = add_q$qvalues

add_q_w_inter = qvalue(add_p_w_inter)
add_q_w_inter = add_q_w_inter$qvalues

### significant interaction associations
inter_p = inter_p[!rownames(inter_p) %in% c('Hypertension', 'Snoring'), ]
inter_p_sub = inter_p[, colnames(inter_p) %in% paste0(sig_inter0$target, "*", sig_inter0$regulator)]
inter_q = qvalue(inter_p_sub)
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

