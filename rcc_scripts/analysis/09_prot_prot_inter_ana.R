setwd('/project/xuanyao/jinghui')
library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
library(gaston)
library(RColorBrewer)
library(tidyverse)
library(simtrait)
library(ggVennDiagram)

qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

## read all interaction p values
inter_p = list()
for (i in 1:22) {
  cis_chr = i
  inter_i = c()
  for (j in 1:22) {
    trans_chr = j
    if (i == j) {
      next
    }
    p_ij = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/obs_p/cis_chr', 
                          cis_chr, '_trans_chr', trans_chr, '.rds'))
    inter_ij = p_ij$inter_p
    inter_i = rbind(inter_i, inter_ij)
  }
  for (k in 1:ncol(inter_i)) {
    prot_p = inter_i[, k]
    names(prot_p) = rownames(inter_i)
    inter_p[[length(inter_p) + 1]] = prot_p
  }
  names(inter_p)[(length(inter_p) - ncol(inter_i) + 1):length(inter_p)] = colnames(inter_i)
}
all_inter_p = unlist(inter_p)
inter_p_tab = data.frame(pval = unname(all_inter_p), 
                         target = sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 1),
                         regulator = sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 2))
inter_p_tab$qval = qvalue(inter_p_tab$pval)$qvalues
inter_p_tab$sig = inter_p_tab$qval < 0.05
p_infl = pval_infl(inter_p_tab$pval)
qqplot.pvalues(inter_p_tab$pval, col.abline = "black", pch = 19, 
               main = expression('Q-Q plot of OLS interaction p, '*lambda*' = 1.1'))

n_sig_target = aggregate(sig ~ target, data = inter_p_tab, sum)
hot_target_prot = n_sig_target$target[n_sig_target$sig >= 10]
inter_p_tab$is_hot_target = inter_p_tab$target %in% hot_target_prot
inter_sig_tab = inter_p_tab[inter_p_tab$qval < 0.05, ]
inter_sig_tab_nonhot = inter_sig_tab[!inter_sig_tab$is_hot_target, ]
inter_sig_tab_nonhot = inter_sig_tab_nonhot[order(inter_sig_tab_nonhot$pval), ]
inter_sig_tab_test = inter_sig_tab[inter_sig_tab$is_hot_target, ]
inter_sig_tab_test = rbind(inter_sig_tab_test, inter_sig_tab_nonhot[1:10, ])
inter_sig_tab_test = inter_sig_tab_test[order(inter_sig_tab_test$target), ]
fwrite(inter_sig_tab_test, 'cis_trans_inter/07_prot_inter/prot_perm_test.txt', sep = '\t')

qqplot.pvalues(inter_p_tab$pval[!inter_p_tab$target %in% hot_target_prot], col.abline = "black", pch = 19, 
               main = 'observed prot * prot interaction p of hotspot target vs. nonhotspot')
qq_plot(inter_p_tab$pval[inter_p_tab$target %in% hot_target_prot], 'red', pch = 19)
legend('topleft', legend = c("non-hotspot", "hotspot"),
       col=c("black", "red"), pch=c(20,20), cex=0.8)

## read all interaction random p values
# inter_perm_p = list()
# for (i in 1:22) {
#   cis_chr = i
#   inter_i = c()
#   for (j in 1:22) {
#     trans_chr = j
#     if (i == j) {
#       next
#     }
#     p_ij = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/perm_p/cis_chr', 
#                           cis_chr, '_trans_chr', trans_chr, '.rds'))
#     inter_ij = p_ij$inter_p
#     inter_i = rbind(inter_i, inter_ij)
#   }
#   for (k in 1:ncol(inter_i)) {
#     prot_p = inter_i[, k]
#     names(prot_p) = rownames(inter_i)
#     inter_perm_p[[length(inter_perm_p) + 1]] = prot_p
#   }
#   names(inter_perm_p)[(length(inter_perm_p) - 
#                          ncol(inter_i) + 1):length(inter_perm_p)] = colnames(inter_i)
# }
# all_inter_perm_p = unlist(inter_perm_p)
# inter_perm_p_tab = data.frame(pval = unname(all_inter_perm_p), 
#                          target = sapply(strsplit(names(all_inter_perm_p), '.', fixed = T), '[', 1),
#                          regulator = sapply(strsplit(names(all_inter_perm_p), '.', fixed = T), '[', 2))
# inter_perm_p_tab$qval = qvalue(inter_perm_p_tab$pval)$qvalues
# inter_perm_p_tab$sig = inter_perm_p_tab$qval < 0.05
# n_sig_target_perm = aggregate(sig ~ target, data = inter_perm_p_tab, sum)
# hot_target_prot_perm = n_sig_target_perm$target[n_sig_target_perm$sig >= 10]

qqplot.pvalues(inter_perm_p_tab$pval[!inter_perm_p_tab$target %in% hot_target_prot_perm], 
               col.abline = "black", pch = 19, ylim = c(0,110),
               main = 'permute prot * prot interaction p of hotspot target vs. nonhotspot')
qq_plot(inter_perm_p_tab$pval[inter_perm_p_tab$target %in% hot_target_prot_perm], 'red', pch = 19)
legend('topleft', legend = c("non-hotspot", "hotspot"),
       col=c("black", "red"), pch=c(20,20), cex=0.8)

pval_infl(all_inter_perm_p, df = 1)
qqplot.pvalues(all_inter_p, col.abline = "black", pch = 19, 
               main = '')
qq_plot(all_inter_perm_p, 'red', pch = 19)
legend('topleft', legend = c("obs", "rand"),
       col=c("black", "red"), pch=c(20,20), cex=0.8)

## using q value (fdr) < 0.05 as the cutoff
q_val = qvalue(all_inter_p)
sig_p = all_inter_p[q_val$qvalues < 0.05]

# q_val_perm = qvalue(all_inter_perm_p)
# sig_perm_p = all_inter_perm_p[q_val_perm$qvalues < 0.05]

## using permuted p to correct obs p
# p_obs_rank = frank(all_inter_p)
# names(p_obs_rank) = names(all_inter_p)
# name.obs.smallp = names(all_inter_p)[all_inter_p < 10^(-4)]
# 
# names(all_inter_perm_p) = paste0(names(all_inter_perm_p), '.perm')
# all_rank = frank(c(all_inter_p, all_inter_perm_p))
# names(all_rank) = c(names(all_inter_p), names(all_inter_perm_p))
# 
# q = pmin(all_rank[name.obs.smallp]/p_obs_rank[name.obs.smallp]-1, 1)
# res = data.frame("snp" = name.obs.smallp, "p" = all_inter_p[all_inter_p < 10^(-4)], "q" = q)
# res_sig = res[res$q < 0.1, ]

## plot significant p values
gene_pos = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')
chr_pos = readRDS('pqtl/00_ref/chromosome_location.rds')
sig_p_tab = data.frame(p = sig_p, q = q_val$qvalues[q_val$qvalues < 0.05])
sig_p_tab$target = sapply(strsplit(names(sig_p), '.', fixed = T), '[', 1)
sig_p_tab$regulator = sapply(strsplit(names(sig_p), '.', fixed = T), '[', 2)
sig_p_tab$target_chr = gene_pos$chr[match(sig_p_tab$target, gene_pos$gene)]
sig_p_tab$target_tss = gene_pos$start_hg19[match(sig_p_tab$target, gene_pos$gene)]
sig_p_tab$regulator_chr = gene_pos$chr[match(sig_p_tab$regulator, gene_pos$gene)]
sig_p_tab$regulator_tss = gene_pos$start_hg19[match(sig_p_tab$regulator, gene_pos$gene)]
sig_p_tab$target_pos = sig_p_tab$target_tss + chr_pos$tot[match(sig_p_tab$target_chr, chr_pos$CHR)]
sig_p_tab$regulator_pos = sig_p_tab$regulator_tss + chr_pos$tot[match(sig_p_tab$regulator_chr, chr_pos$CHR)]

# sig_perm_p_tab = data.frame(p = sig_perm_p, q = q_val_perm$qvalues[q_val_perm$qvalues < 0.05])
# sig_perm_p_tab$target = sapply(strsplit(names(sig_perm_p), '.', fixed = T), '[', 1)
# sig_perm_p_tab$regulator = sapply(strsplit(names(sig_perm_p), '.', fixed = T), '[', 2)
# sig_perm_p_tab$target_chr = gene_pos$chr[match(sig_perm_p_tab$target, gene_pos$gene)]
# sig_perm_p_tab$target_tss = gene_pos$start_hg19[match(sig_perm_p_tab$target, gene_pos$gene)]
# sig_perm_p_tab$regulator_chr = gene_pos$chr[match(sig_perm_p_tab$regulator, gene_pos$gene)]
# sig_perm_p_tab$regulator_tss = gene_pos$start_hg19[match(sig_perm_p_tab$regulator, gene_pos$gene)]
# sig_perm_p_tab$target_pos = sig_perm_p_tab$target_tss + chr_pos$tot[match(sig_perm_p_tab$target_chr, chr_pos$CHR)]
# sig_perm_p_tab$regulator_pos = sig_perm_p_tab$regulator_tss + chr_pos$tot[match(sig_perm_p_tab$regulator_chr, chr_pos$CHR)]

ggplot(sig_p_tab) +
  geom_point(aes(x = regulator_pos, y = target_pos,size = -log10(p)),
    alpha = 1, shape = 1) +
  labs(y = "Target chromosome", x = "Regulator chromosome", size = quote(-log[10](p))) +
  scale_x_continuous(limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
    label = chr_pos$CHR, breaks = chr_pos$center, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
    label = chr_pos$CHR, breaks = chr_pos$center, expand = c(0, 0)) + 
  scale_size(guide = guide_legend(override.aes = list(alpha = 1)), range = c(1, 5)) +
  geom_hline(data = chr_pos, 
    aes(yintercept = xmax),
    linetype = "dotted", color = "grey") +
  geom_vline(data = chr_pos, 
             aes(xintercept = xmax),
             linetype = "dotted", color = "grey") +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey') + 
  theme(text = element_text(size=15, colour = "black"), 
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())

### investigate the regulator hotspots
## 1. if hotspots have higher h2
h2_tab = fread('cis_trans_inter/08_prot_h2/cis_h2_summ.txt')
n_inter_perm = as.data.frame(table(c(sig_perm_p_tab$target, sig_perm_p_tab$regulator)))
n_target_perm = as.data.frame(table(c(sig_perm_p_tab$target)))
n_reg_perm = as.data.frame(table(c(sig_perm_p_tab$regulator)))
h2_tab$n_inter_perm = n_inter_perm$Freq[match(h2_tab$prot, n_inter_perm$Var1)]
h2_tab$n_target_perm = n_target_perm$Freq[match(h2_tab$prot, n_target_perm$Var1)]
h2_tab$n_reg_perm = n_reg_perm$Freq[match(h2_tab$prot, n_reg_perm$Var1)]

h2_tab$n_inter[is.na(h2_tab$n_inter)] = 0
h2_tab$n_target[is.na(h2_tab$n_target)] = 0
h2_tab$n_reg[is.na(h2_tab$n_reg)] = 0

h2_tab$n_inter_perm[is.na(h2_tab$n_inter_perm)] = 0
h2_tab$n_target_perm[is.na(h2_tab$n_target_perm)] = 0
h2_tab$n_reg_perm[is.na(h2_tab$n_reg_perm)] = 0

h2_tab_sub = h2_tab[h2_tab$R2 > 0.01, ]
h2_tab_sub$r2_cat = cut(h2_tab_sub$R2, 
                        breaks = quantile(h2_tab_sub$R2, seq(0,1,0.2)),
                        include.lowest = T)
h2_tab_sub$inter_hot = h2_tab_sub$n_inter >= 10
h2_tab_sub$target_hot = h2_tab_sub$n_target >= 10
h2_tab_sub$inter_perm_hot = h2_tab_sub$n_inter_perm >= 10
h2_tab_sub$target_perm_hot = h2_tab_sub$n_target_perm >= 10
## target hotspot
ggplot(h2_tab_sub, aes(x= r2_cat,  y = n_target, color = r2_cat)) +
  geom_boxplot() +
  labs(x = "R2", y = '# targets', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

hotspot_target = aggregate(target_hot ~ r2_cat, data = h2_tab_sub, sum)
ggplot(hotspot_target, aes(x= r2_cat,  y = target_hot, fill = r2_cat, label = target_hot)) +
  geom_bar(stat="identity") +
  geom_text(position=position_dodge(0.9), vjust=-0.2) + 
  labs(x = "R2", y = '# hotspots', title = 'Obs target', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

ggplot(h2_tab_sub, aes(x= r2_cat,  y = n_target_perm, color = r2_cat)) +
  geom_boxplot() +
  labs(x = "R2", y = '# targets', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

hotspot_perm_target = aggregate(target_perm_hot ~ r2_cat, data = h2_tab_sub, sum)
ggplot(hotspot_perm_target, aes(x= r2_cat,  y = target_perm_hot, fill = r2_cat, 
                                label = target_perm_hot)) +
  geom_bar(stat="identity") +
  geom_text(position=position_dodge(0.9), vjust=-0.2) + 
  labs(x = "R2", y = '# hotspots', title = 'Perm target', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

## interaction hotspot
ggplot(h2_tab_sub, aes(x= r2_cat,  y = n_inter, color = r2_cat)) +
  geom_boxplot() +
  labs(x = "R2", y = '# interactions', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

hotspot_inter = aggregate(inter_hot ~ r2_cat, data = h2_tab_sub, sum)
ggplot(hotspot_inter, aes(x= r2_cat,  y = inter_hot, fill = r2_cat, label = inter_hot)) +
  geom_bar(stat="identity") +
  geom_text(position=position_dodge(0.9), vjust=-0.2) + 
  labs(x = "R2", y = '# hotspots', title = 'Obs interaction', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

ggplot(h2_tab_sub, aes(x= r2_cat,  y = n_inter_perm, color = r2_cat)) +
  geom_boxplot() +
  labs(x = "R2", y = '# interactions', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

hotspot_perm_inter = aggregate(inter_perm_hot ~ r2_cat, data = h2_tab_sub, sum)
ggplot(hotspot_perm_inter, aes(x= r2_cat,  y = inter_perm_hot, fill = r2_cat, 
                               label = inter_perm_hot)) +
  geom_bar(stat="identity") +
  geom_text(position=position_dodge(0.9), vjust=-0.2) + 
  labs(x = "R2", y = '# hotspots', title = 'Perm interaction', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

target_overlap = list(obs = h2_tab_sub$prot[h2_tab_sub$target_hot], 
                       rand = h2_tab_sub$prot[h2_tab_sub$target_perm_hot])
ggVennDiagram(target_overlap) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = 'Target hotspot overlap', color = '') + 
  theme(legend.position = "none")

inter_overlap = list(obs = h2_tab_sub$prot[h2_tab_sub$inter_hot], 
                     rand = h2_tab_sub$prot[h2_tab_sub$inter_perm_hot])
ggVennDiagram(inter_overlap) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = 'Interaction hotspot overlap', color = '') + 
  theme(legend.position = "none")



## regulator hotspot
ggplot(h2_tab_sub, aes(x= r2_cat,  y = n_reg, color = r2_cat)) +
  geom_boxplot() +
  labs(x = "R2", y = '# regulators', title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

prot_hotspot = data.frame(prot = h2_tab_sub$prot[h2_tab_sub$target_hot])
prot_non_hotspot = data.frame(prot = h2_tab_sub$prot[!h2_tab_sub$target_hot])
fwrite(prot_non_hotspot, 'cis_trans_inter/prot_non_hot.txt', sep = '\t')

## 2. if hotspots have stronger cis-pQTLs
ukb_cis = fread('cis_trans_inter/00_ref/ukb_cis_pqtl.txt')
h2_tab$cis_logP = ukb_cis$logP[match(h2_tab$prot, ukb_cis$prot)]
cor.test(h2_tab$n_inter, h2_tab$cis_logP, use = 'na.or.complete')
plot(h2_tab$n_inter, h2_tab$cis_logP, xlab = '# of sig regulator', ylab = 'cis-pQTL -log10(p)', 
     main = 'cor = 0.11 (p = 0.15)')


#### enrichment in TF and PPI
## tf enrichment
tf_gene_id = fread('pqtl/15_TF/humantf_ccbr/TFs_Ensembl_v_1.01.txt', header = F)
hg38_gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.unadj.gtf')
hg38_gene_meta$gene_id = sapply(strsplit(hg38_gene_meta$gene_id, '.', fixed = T), '[', 1)
hg38_gene_meta = hg38_gene_meta[hg38_gene_meta$gene_type == 'protein_coding']
tf_gene_id$gene_name = hg38_gene_meta$gene_name[match(tf_gene_id$V1, hg38_gene_meta$gene_id)]

gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta_sub = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
gene_meta_sub = gene_meta_sub[!grepl('gene_status', gene_meta_sub$gene_name), ]
gene_meta_sub = gene_meta_sub[!grepl('ENSG', gene_meta_sub$gene_name), ]

gene_meta_sub$gene_length = gene_meta_sub$end - gene_meta_sub$start
sig_p_tab_sub = sig_p_tab[sig_p_tab$regulator %in% gene_meta_sub$gene_name, ]
background_tf = matrix(rep(F, nrow(sig_p_tab_sub) * 1000), ncol = 1000)
for (i in 1:nrow(sig_p_tab_sub)) {
  gene_i = sig_p_tab_sub$regulator[i]
  gene_length_i = gene_meta_sub$gene_length[match(gene_i, gene_meta_sub$gene_name)]
  gene_meta_i = gene_meta_sub[gene_meta_sub$gene_length >= 0.9*gene_length_i & 
                                gene_meta_sub$gene_length <= 1.1*gene_length_i, ]
  background_gene_i = sample(gene_meta_i$gene_name, 1000, replace = T)
  background_tf[i, ] = background_gene_i %in% tf_gene_id$gene_name
  if (i %% 1000 == 0) {
    print(i)
  }
}
background_tf = apply(background_tf, 2, sum)
obs_tf = sum(sig_p_tab_sub$regulator %in% tf_gene_id$gene_name)
tf_enrich = obs_tf / background_tf


## PPI enrichment
ppi_list = readRDS('pqtl/11_prot_complex/ppi_list.rds')
all_prot = unique(sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 1))
background_ppi = c()
for (i in 1:1000) {
  gene_pair_i = paste0(sample(all_prot), ',', all_prot)
  background_ppi[i] = sum(gene_pair_i %in% ppi_list)
  if (i %% 100 == 0) {
    print(i)
  }
}

enrich_all = data.frame(enrich = c(obs_tf/mean(background_tf), obs_ppi/mean(background_ppi)),
                        ci_low = c(obs_tf/quantile(background_tf, 0.025), obs_ppi/quantile(background_ppi, 0.025)),
                        ci_high = c(obs_tf/quantile(background_tf, 0.975), obs_ppi/quantile(background_ppi, 0.975)),
                        type = c('TF', 'PPI'))
# p value of tf enrichment
tf_enrich_p = obs_tf/background_tf - 1
empirical_p_tf = mean(tf_enrich_p > 0) 
2*min(c(empirical_p_tf, 1-empirical_p_tf))

# p value of ppi enrichment
ppi_enrich_p = obs_ppi/background_ppi - 1
empirical_p_ppi = mean(ppi_enrich_p > 0) 
2*min(c(empirical_p_ppi, 1-empirical_p_ppi))

ggplot(enrich_all, aes(x = type, y = enrich, fill = type)) + 
  geom_bar(stat="identity", position=position_dodge(0.1), width = 0.5) +
  geom_errorbar(position=position_dodge(0.1), aes(ymin=ci_low, ymax=ci_high), width=.2) + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "", y = "Enrichment", fill = '',  title = '') + 
  scale_fill_manual(values = brewer.pal(3, "Set1")[1:2]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')




