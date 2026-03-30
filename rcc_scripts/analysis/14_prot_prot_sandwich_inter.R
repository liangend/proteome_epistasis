setwd('/project/xuanyao/jinghui')
library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
library(gaston)
library(RColorBrewer)
library(simtrait)

qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

## read all interaction p values after sandwich variance correction
inter_p = list()
trans_p = list()
for (i in 1:22) {
  cis_chr = i
  inter_i = c()
  trans_i = c()
  for (j in 1:22) {
    trans_chr = j
    p_ij = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/sandwich_p/cis_chr', 
                          cis_chr, '_trans_chr', trans_chr, '.rds'))
    inter_i = rbind(inter_i, p_ij$inter_p)
    trans_i = rbind(trans_i, p_ij$trans_p)
  }
  for (k in 1:ncol(inter_i)) {
    prot_p = inter_i[, k]
    names(prot_p) = rownames(inter_i)
    inter_p[[length(inter_p) + 1]] = prot_p
    trans_p[[length(inter_p) + 1]] = trans_i[, k]
  }
  names(inter_p)[(length(inter_p) - ncol(inter_i) + 1):length(inter_p)] = colnames(inter_i)
}
all_inter_p = unlist(inter_p)
all_trans_p = unlist(trans_p)

inter_p_tab = data.frame(pval = unname(all_inter_p), 
                         target = sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 1),
                         regulator = sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 2))
trans_p_tab = data.frame(pval = unname(all_trans_p), 
                         target = inter_p_tab$target,
                         regulator = inter_p_tab$regulator)

## remove interactions between the same protein
inter_p_tab = inter_p_tab[inter_p_tab$target != inter_p_tab$regulator, ]

## remove interactions between proteins that are close to each other
cis_window = 1000000
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
inter_p_tab$target_chr = ukb_prot$chr[match(inter_p_tab$target, ukb_prot$gene_name)]
inter_p_tab$reg_chr = ukb_prot$chr[match(inter_p_tab$regulator, ukb_prot$gene_name)]
# inter_p_tab$target_tss = ukb_prot$gene_start[match(inter_p_tab$target, ukb_prot$gene_name)]
# inter_p_tab$reg_tss = ukb_prot$gene_start[match(inter_p_tab$regulator, ukb_prot$gene_name)]
# inter_p_tab$reg_target_dist = abs(inter_p_tab$target_tss - inter_p_tab$reg_tss)
inter_p_tab_sub = inter_p_tab[!(inter_p_tab$reg_chr == inter_p_tab$target_chr), ]

p_infl = pval_infl(inter_p_tab_sub$pval)
qqplot.pvalues(inter_p_tab_sub$pval, col.abline = "black", pch = 19, 
               main = expression('Q-Q plot of SVE corrected interaction p, '*lambda*' = 1.0'))

## define significant interactions as q < 0.05
inter_q = qvalue(inter_p_tab_sub$pval)
inter_p_tab_sub$qval = inter_q$qvalues
inter_p_tab_sub$sig = inter_p_tab_sub$qval < 0.05

# pval_infl(inter_p_tab_sub$pval, df = 1)
# qqplot.pvalues(inter_p_tab_sub$pval, col.abline = "black", pch = 19, 
#                main = 'QQ-plot of variance-corrected interaction pval, inflation factor = 1.01')

inter_sig = inter_p_tab_sub[inter_p_tab_sub$sig, ]
# fwrite(inter_sig, 'cis_trans_inter/07_prot_inter/sig_prot_prot_inter_w_cis.txt',
#        sep = '\t')
inter_sig = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')

## remove protein pairs that are both located in HLA region
# rm_pair = which(inter_sig$target_chr == 6 & inter_sig$reg_chr == 6 &
#                   inter_sig$target_tss > 28477797 &  inter_sig$target_tss < 33448354 &
#                   inter_sig$reg_tss > 28477797 & inter_sig$reg_tss < 33448354)
# inter_sig = inter_sig[-rm_pair, ]
# inter_sig = inter_sig[!inter_sig$target %in% c('BTN3A2', 'HLA-A') & 
#                         inter_sig$regulator != 'BTN3A2', ]

inter_sig_trans_p = c()
for (i in 1:nrow(inter_sig)) {
  inter_sig_trans_p[i] = trans_p_tab$pval[which(trans_p_tab$target == inter_sig$target[i] &
                                                  trans_p_tab$regulator == inter_sig$regulator[i])] 
}


## raindrop of significant interactions
# gene_pos = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')
# chr_pos = readRDS('pqtl/00_ref/chromosome_location.rds')
# 
# inter_sig$target_pos = inter_sig$target_tss + chr_pos$tot[match(inter_sig$target_chr, chr_pos$CHR)]
# inter_sig$reg_pos = inter_sig$reg_tss + chr_pos$tot[match(inter_sig$reg_chr, chr_pos$CHR)]
# ggplot(inter_sig) +
#   geom_point(aes(x = reg_pos, y = target_pos,size = -log10(pval)),
#              alpha = 1, shape = 1) +
#   labs(y = "Target chromosome", x = "Regulator chromosome", size = quote(-log[10](p))) +
#   scale_x_continuous(limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
#                      label = chr_pos$CHR, breaks = chr_pos$center, expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
#                      label = chr_pos$CHR, breaks = chr_pos$center, expand = c(0, 0)) + 
#   scale_size(guide = guide_legend(override.aes = list(alpha = 1)), range = c(1, 5)) +
#   geom_hline(data = chr_pos, 
#              aes(yintercept = xmax),
#              linetype = "dotted", color = "grey") +
#   geom_vline(data = chr_pos, 
#              aes(xintercept = xmax),
#              linetype = "dotted", color = "grey") +
#   geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey') + 
#   theme(text = element_text(size=15, colour = "black"), 
#         axis.text.x = element_text(colour = "black", size = 10),
#         axis.text.y = element_text(colour = "black", size = 10),
#         axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())


### plot interaction network
library(igraph)
inter_net = graph_from_data_frame(inter_sig[, c(3,2)], directed = TRUE)

# use different colors for targets and regulators
V(inter_net)$color = "grey"
edges = ends(inter_net, E(inter_net))

reg_nodes = unique(edges[,1])
target_nodes = unique(edges[,2])

# Assign colors
V(inter_net)$color[V(inter_net)$name %in% reg_nodes] = "lightsalmon"   # regulator nodes (from)
V(inter_net)$color[V(inter_net)$name %in% target_nodes] = "lightskyblue1"  # target nodes (to)
degrees = rep(30, 9)
degrees[1] = 40
degrees[4] = 35
plot(inter_net,
     vertex.label.cex = 0.8,
     vertex.label.color = 'black',
     vertex.size = degrees,
     node.size = 1,
     edge.arrow.size = 0.02,
     edge.color = "black",
     layout = layout_with_kk)

## interaction example
# 1. MBL2 * ABO
abo_tab = fread('cis_trans_inter/06_fusion_pred/eur/ABO_good_pred.txt')
mbl2_tab = fread('cis_trans_inter/06_fusion_pred/eur/MBL2_good_pred.txt')
blood_type = fread('cis_trans_inter/data/ukb_blood_type.pheno')

abo_mbl2 = data.frame(id = mbl2_tab$id, 
                      mbl2_obs = mbl2_tab$obs, mbl2_pred = mbl2_tab$pred)
abo_mbl2$abo_pred = abo_tab$pred[match(abo_mbl2$id, abo_tab$id)]
ggplot(abo_mbl2, 
       aes(x = mbl2_pred, y = (mbl2_pred - mbl2_obs)^2)) + 
  geom_point(alpha = 0.3) + 
  labs(x = "Predicted MBL2", y = expression("Residual"^2), 
       title = "Heteroscedasticity of predicted MBL2") + 
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

abo_mbl2$blood_type = blood_type$`23165`[match(abo_mbl2$id, blood_type$FID)]
abo_mbl2 = abo_mbl2[!is.na(abo_mbl2$blood_type), ]
abo_mbl2$blood_O = 'not O'
abo_mbl2$blood_O[which(abo_mbl2$blood_type == 'OO')] = 'O'

ggplot(abo_mbl2, 
       aes(x=mbl2_pred, y=mbl2_obs, color=blood_O)) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(method='lm') +
  labs(x = "Predicted MBL2", y = 'Observed MBL2', 
       title = expression("Different " * italic("cis") * 
                            "-effects of MBL2 across ABO blood types"), color = '') + 
  scale_color_manual(values = brewer.pal(8, "Paired")[c(2,8)]) +
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

## ABO interacting proteins with different effects based on different number of O alleles
abo_inter_prot = c('MBL2', 'CD209', 'ICAM5', 
                   'NPTX1', 'FCGR2B')
o_eff = c()
not_o_eff = c()
o_eff_sd = c()
not_o_eff_sd = c()
inter_eff_p = c()
for (i in abo_inter_prot) {
  prot_i = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', 
                        i, '_good_pred.txt'))
  prot_i$blood_type = blood_type$`23165`[match(prot_i$id, blood_type$FID)]
  prot_i = prot_i[!is.na(prot_i$blood_type), ]
  prot_i$blood_O = 'not O'
  prot_i$blood_O[which(prot_i$blood_type == 'OO')] = 'O'
  prot_i_o = prot_i[prot_i$blood_O == 'O', ]
  prot_i_not_o = prot_i[prot_i$blood_O == 'not O', ]
  fit_o = summary(lm(obs ~ pred, data = prot_i_o))
  fit_not_o = summary(lm(obs ~ pred, data = prot_i_not_o))
  o_eff = c(o_eff, fit_o$coefficients[2,1])
  not_o_eff = c(not_o_eff, fit_not_o$coefficients[2,1])
  o_eff_sd = c(o_eff_sd, fit_o$coefficients[2,2])
  not_o_eff_sd = c(not_o_eff_sd, fit_not_o$coefficients[2,2])
  fit_inter = summary(lm(obs ~ pred * blood_O, data = prot_i))
  inter_eff_p = c(inter_eff_p, fit_inter$coefficients[4,4])
  print(i)
}

eff_tab = data.frame(beta = c(o_eff, not_o_eff), 
                     sd = c(o_eff_sd, not_o_eff_sd),
                     prot = rep(abo_inter_prot, 2),
                     p = rep(inter_eff_p, 2), 
                     group = rep(c('O blood type', 'not O blood type'), each = 5))
ggplot(eff_tab, aes(x=prot,  y=beta, color = group)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  geom_errorbar(position=position_dodge(0.3), aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=0.5) +
  labs(x = "", y = 'cis-regultory effect', title = '', color = '') + 
  scale_color_manual(values = brewer.pal(8, "Paired")[c(2,8)]) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



