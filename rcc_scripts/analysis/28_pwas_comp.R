setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(openxlsx)
library(RColorBrewer)
library(gridExtra)
all_reg_file = list.files('cis_trans_inter/15_pheno_asso_pred_w_inter/asso_1_add_inter/')
inter_sig = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter_w_cis.txt')

ukb_trait = read.xlsx('cis_trans_inter/00_ref/ukb_trait.xlsx', sheet = 3)

### PWAS p value between target predicted expression and phenotype
add1_p_wo_inter = c() # target ~ target_cis_pred
add1_p_w_inter = c()  # target ~ target_cis_pred + target_cis_pred * regulator_cis_pred 
#add2_p_wo_inter = c() # target ~ target_cis_pred + regulator_cis_pred
#add2_p_w_inter = c()  # target ~ target_cis_pred + regulator_cis_pred + target_cis_pred * regulator_cis_pred 
abo_p = c()

for (i in 1:length(all_reg_file)) {
  trait_i = all_reg_file[i]
  # target ~ target_cis_pred
  coef_i = readRDS(paste0('cis_trans_inter/15_pheno_asso_pred_w_inter/center/asso_1_add/', trait_i))
  add1_p_wo_inter = rbind(add1_p_wo_inter, sapply(coef_i, function(x){x[2,4]}))
  
  # target ~ target_cis_pred + target_cis_pred * regulator_cis_pred 
  coef_i = readRDS(paste0('cis_trans_inter/15_pheno_asso_pred_w_inter/center/asso_1_add_inter/', trait_i))
  add1_p_w_inter = rbind(add1_p_w_inter, sapply(coef_i, function(x){x[2,4]}))
  
  # # target ~ target_cis_pred + regulator_cis_pred
  # coef_i = readRDS(paste0('cis_trans_inter/15_pheno_asso_pred_w_inter/asso_2_add/', trait_i))
  # add2_p_wo_inter = rbind(add2_p_wo_inter, sapply(coef_i, function(x){x[2,4]}))
  # 
  # # target ~ target_cis_pred + regulator_cis_pred + target_cis_pred * regulator_cis_pred 
  # coef_i = readRDS(paste0('cis_trans_inter/15_pheno_asso_pred_w_inter/asso_2_add_inter/', trait_i))
  # add2_p_w_inter = rbind(add2_p_w_inter, sapply(coef_i, function(x){x[2,4]}))
  
  # abo pwas value
  coef_i = readRDS(paste0('cis_trans_inter/16_pheno_asso_abo/', trait_i))
  abo_p[i] = coef_i[2,4]
}

rownames(add1_p_wo_inter) = sub('.rds', '', all_reg_file)
rownames(add1_p_w_inter) = sub('.rds', '', all_reg_file)
names(abo_p) = sub('.rds', '', all_reg_file)

## remove protein pairs that are located in HLA region
rm_pair = which((inter_sig$target_chr == 6 & inter_sig$reg_chr == 6 &
                   inter_sig$target_tss > 28477797 &  inter_sig$target_tss < 33448354 &
                   inter_sig$reg_tss > 28477797 & inter_sig$reg_tss < 33448354) | 
                  inter_sig$target %in% c('BTN3A2', 'HLA-A') |
                  inter_sig$regulator == 'BTN3A2')
keep_pair = setdiff(1:31, rm_pair)

## target ~ target_cis_pred vs. target ~ target_cis_pred + target_cis_pred * regulator_cis_pred 
add1_p_wo_inter = as.data.frame(add1_p_wo_inter)
add1_p_wo_inter$trait = rownames(add1_p_wo_inter)
add1_p_wo_inter$trait_group = ukb_trait$Group[match(add1_p_wo_inter$trait, 
                                                    ukb_trait$abbrev)]

add1_p_w_inter = as.data.frame(add1_p_w_inter)
add_p_plt = cbind(add1_p_wo_inter, add1_p_w_inter)
add_p_plt = na.omit(add_p_plt)

plt_tab_all = c()
for (i in 1:length(keep_pair)) {
  pair_i = keep_pair[i]
  plt_tab_i = add_p_plt[, c(pair_i, 32, 33, pair_i+33)]
  inter_i = paste0(inter_sig$target[pair_i], ' * ', inter_sig$regulator[pair_i])
  colnames(plt_tab_i)[c(1,4)] = c('wo_inter', 'w_inter')
  plt_tab_i$inter = inter_i
  plt_tab_all = rbind(plt_tab_all, plt_tab_i)
}
plt_tab_sub = plt_tab_all[plt_tab_all$trait_group %in% 
                            c('Anthropometric', 'Blood biochemistry', 'Blood cell',
                              'Body fat', 'Cardiometabolic', unique(plt_tab_all$trait_group)[1], 
                              'Mental health', 'Urine metabolite'), ]

sig_inter_diff_chr = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')
sig_pair = paste0(sig_inter_diff_chr$target, ' * ', sig_inter_diff_chr$regulator)
plt_tab_sub = plt_tab_sub[plt_tab_sub$inter %in% sig_pair, ]
ggplot(plt_tab_sub, aes(x = -log10(w_inter), y = -log10(wo_inter), 
                      color = trait_group)) + 
  geom_point(size = 3) +
  facet_wrap( ~ inter, nrow = 2, scales = "free") + 
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  labs(x = expression("-log"[10] * "(" *italic(P)[adj]* ") with interaction") , 
       y = expression("-log"[10] * "(" *italic(P)[raw]* ") without interaction"), color = '') +
  scale_color_manual(values = brewer.pal(8, "Dark2")) + 
  theme(text = element_text(size=13, colour = "black"), 
        strip.background = element_rect(color="black", fill="bisque"),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

sig_inter_pheno = plt_tab_all[plt_tab_all$w_inter < 0.05/136/7, ]
colnames(sig_inter_pheno)[c(1,4)] = c('pwas_p_wo_inter', 'pwas_p_w_inter')
# fwrite(sig_inter_pheno, 'cis_trans_inter/15_pheno_asso_pred_w_inter/sig_inter_pheno.txt', sep = '\t')

## target ~ target_cis_pred vs. target ~ target_cis_pred + regulator_cis_pred 
add2_p_wo_inter = as.data.frame(add2_p_wo_inter)
add_p_plt2 = cbind(add1_p_wo_inter, add2_p_wo_inter)
add_p_plt2 = na.omit(add_p_plt2)

comp_plt2 = list()
for (i in 1:length(keep_pair)) {
  pair_i = keep_pair[i]
  plt_tab_i = add_p_plt2[, c(pair_i, 32, 33, pair_i+33)]
  inter_i = paste0(inter_sig$target[pair_i], '*', inter_sig$regulator[pair_i])
  colnames(plt_tab_i)[c(1,4)] = c('wo_inter', 'w_inter')
  comp_plt2[[i]] = ggplot(plt_tab_i, aes(x = -log10(w_inter), y = -log10(wo_inter), 
                                        color = trait_group)) + 
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, lty = 2) + 
    labs(x = "-log10(p) with reg", 
         y = "-log10(p) without reg", color = '', title = inter_i) +
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
do.call(grid.arrange, c(comp_plt2, nrow = 5))

## target ~ target_cis_pred vs. target ~ target_cis_pred + regulator_cis_pred + target_cis_pred*regulator_cis_pred
add2_p_w_inter = as.data.frame(add2_p_w_inter)
add_p_plt3 = cbind(add1_p_wo_inter, add2_p_w_inter)
add_p_plt3 = na.omit(add_p_plt3)

comp_plt3 = list()
for (i in 1:length(keep_pair)) {
  pair_i = keep_pair[i]
  plt_tab_i = add_p_plt3[, c(pair_i, 32, 33, pair_i+33)]
  inter_i = paste0(inter_sig$target[pair_i], '*', inter_sig$regulator[pair_i])
  colnames(plt_tab_i)[c(1,4)] = c('wo_inter', 'w_inter')
  comp_plt3[[i]] = ggplot(plt_tab_i, aes(x = -log10(w_inter), y = -log10(wo_inter), 
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
do.call(grid.arrange, c(comp_plt3, nrow = 5))


## target ~ target_cis_pred + regulator_cis_pred vs. 
## target ~ target_cis_pred + regulator_cis_pred + target_cis_pred*regulator_cis_pred
add2_p_wo_inter = as.data.frame(add2_p_wo_inter)
add2_p_wo_inter$trait = rownames(add2_p_wo_inter)
add2_p_wo_inter$trait_group = ukb_trait$Group[match(add2_p_wo_inter$trait, 
                                                    ukb_trait$abbrev)]

add2_p_w_inter = as.data.frame(add2_p_w_inter)
add_p_plt4 = cbind(add2_p_wo_inter, add2_p_w_inter)
add_p_plt4 = na.omit(add_p_plt4)

comp_plt4 = list()
for (i in 1:length(keep_pair)) {
  pair_i = keep_pair[i]
  plt_tab_i = add_p_plt4[, c(pair_i, 32, 33, pair_i+33)]
  inter_i = paste0(inter_sig$target[pair_i], '*', inter_sig$regulator[pair_i])
  colnames(plt_tab_i)[c(1,4)] = c('wo_inter', 'w_inter')
  comp_plt4[[i]] = ggplot(plt_tab_i, aes(x = -log10(w_inter), y = -log10(wo_inter), 
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
do.call(grid.arrange, c(comp_plt4, nrow = 5))

## ABO PWAS
abo_p_plt = data.frame(p = abo_p, trait = sub('.rds', '', all_reg_file))
abo_p_plt$trait_group = ukb_trait$Group[match(abo_p_plt$trait, 
                                              ukb_trait$abbrev)]
abo_p_plt = abo_p_plt[order(abo_p_plt$p), ]
abo_p_sub = abo_p_plt[1:20, ]
ggplot(abo_p_sub, aes(x = reorder(trait, p), y = -log10(p))) + 
  geom_bar(stat="identity") +
  labs(x = "", 
       y = "PWAS -log10(p)", fill = '', title = 'ABO') +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 60, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')




