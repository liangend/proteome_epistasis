setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(gaston)

qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

## read all interaction p values
inter_p_w_cis = list()
inter_p_wo_cis = list()
for (i in 1:22) {
  cis_chr = i
  inter_i_w_cis = c()
  inter_i_wo_cis = c()
  for (j in 1:22) {
    trans_chr = j
    if (i == j) {
      next
    }
    p_ij = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_snp_inter/cis_chr', 
                          cis_chr, '_trans_chr', trans_chr, '.rds'))
    inter_i_w_cis = rbind(inter_i_w_cis, p_ij$inter_p_w_cis_prot)
    inter_i_wo_cis = rbind(inter_i_wo_cis, p_ij$inter_p_wo_cis_prot)
  }
  for (k in 1:ncol(inter_i_w_cis)) {
    prot_p_w_cis = inter_i_w_cis[, k]
    names(prot_p_w_cis) = rownames(inter_i_w_cis)
    inter_p_w_cis[[length(inter_p_w_cis) + 1]] = prot_p_w_cis
    
    prot_p_wo_cis = inter_i_wo_cis[, k]
    names(prot_p_wo_cis) = rownames(inter_i_wo_cis)
    inter_p_wo_cis[[length(inter_p_wo_cis) + 1]] = prot_p_wo_cis
  }
  names(inter_p_w_cis)[(length(inter_p_w_cis) - 
                          ncol(inter_i_w_cis) + 1):length(inter_p_w_cis)] = colnames(inter_i_w_cis)
  names(inter_p_wo_cis)[(length(inter_p_wo_cis) - 
                          ncol(inter_i_wo_cis) + 1):length(inter_p_wo_cis)] = colnames(inter_i_wo_cis)
}
all_inter_p_w_cis = unlist(inter_p_w_cis)
all_inter_p_wo_cis = unlist(inter_p_wo_cis)

## using q value (fdr) < 0.05 as the cutoff
q_val_w_cis = qvalue(all_inter_p_w_cis)
sig_p_w_cis = all_inter_p_w_cis[q_val_w_cis$qvalues < 0.05]

q_val_wo_cis = qvalue(all_inter_p_wo_cis)
sig_p_wo_cis = all_inter_p_wo_cis[q_val_wo_cis$qvalues < 0.05]

sum(names(sig_p_w_cis) %in% names(sig_p_wo_cis))

## qq-plot of two interaction p values
qqplot.pvalues(all_inter_p_wo_cis, col.abline = "black", pch = 1, 
               main = 'Protein SNP interaction')
qq_plot(all_inter_p_w_cis, col = "red", pch = 1)
legend('topleft', legend = c("without target pred", "with target pred"),
       col=c("black", "red"), pch=c(1,1), cex=0.8)

## plot significant p values in the model with target pred
gene_pos = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')
chr_pos = readRDS('pqtl/00_ref/chromosome_location.rds')
sig_p_tab = data.frame(p = sig_p_w_cis, q = q_val_w_cis$qvalues[q_val_w_cis$qvalues < 0.05])
sig_p_tab$target = sapply(strsplit(names(sig_p_w_cis), '.', fixed = T), '[', 1)
sig_p_tab$regulator = sapply(strsplit(names(sig_p_w_cis), '.', fixed = T), '[', 2)
sig_p_tab$target_chr = gene_pos$chr[match(sig_p_tab$target, gene_pos$gene)]
sig_p_tab$target_tss = gene_pos$start_hg19[match(sig_p_tab$target, gene_pos$gene)]
sig_p_tab$regulator_chr = gene_pos$chr[match(sig_p_tab$regulator, gene_pos$gene)]
sig_p_tab$regulator_tss = gene_pos$start_hg19[match(sig_p_tab$regulator, gene_pos$gene)]
sig_p_tab$target_pos = sig_p_tab$target_tss + chr_pos$tot[match(sig_p_tab$target_chr, chr_pos$CHR)]
sig_p_tab$regulator_pos = sig_p_tab$regulator_tss + chr_pos$tot[match(sig_p_tab$regulator_chr, chr_pos$CHR)]

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
sort(table(sig_p_tab$target),decreasing = T)[1:10]
sort(table(sig_p_tab$regulator),decreasing = T)[1:10]


## Enrichment of cis-pQTL in 3'UTR
cis_pqtl_annot = fread('cis_trans_inter/09_snpEff/ukb_cis_annot.vcf')
ukb_cis_qtl = fread('cis_trans_inter/00_ref/ukb_cis_pqtl.txt')
cis_pqtl_annot$annot = sapply(strsplit(cis_pqtl_annot$INFO, '|', fixed = T),
                              '[', 2)
(table(cis_pqtl_annot$annot)/nrow(cis_pqtl_annot))[1]

sig_p_tab$annot = cis_pqtl_annot$annot[match(sig_p_tab$target, ukb_cis_qtl$prot)]
(table(sig_p_tab$annot)/nrow(sig_p_tab))[1]












