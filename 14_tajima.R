setwd('/project/xuanyao/jinghui')
library(data.table)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(plink2R)

### Tajima's D for each gene
tajima = fread('cis_trans_inter/21_tajima_D/tajdEd.txt.gz')
# gene_meta = fread('mediation_pathway/00_ref/genecode.GRCh38.gene.meta.gtf')
# gene_meta_sub = gene_meta[gene_meta$gene_type == 'protein_coding', ]
# gene_meta_sub = na.omit(gene_meta_sub)
# gene_meta_bed = data.frame(chr = gene_meta_sub$chr, 
#                            start = gene_meta_sub$start_hg37,
#                            end = gene_meta_sub$end_hg37,
#                            gene = gene_meta_sub$gene_name)
# fwrite(gene_meta_bed, 'cis_trans_inter/21_tajima_D/gene_meta.bed', 
#        sep = '\t', col.names = F)
gene_meta = fread('cis_trans_inter/21_tajima_D/gene_meta_hg18.bed')
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')
all_gene = list.files('cis_trans_inter/06_fusion_pred/eur/', pattern = 'good_pred.txt')
all_gene = sub('_good_pred.txt', '', all_gene)
uniq_gene = unique(c(sig_inter$target, sig_inter$regulator))
tajima_all = c()
for (i in 1:length(all_gene)) {
  gene_i = all_gene[i]
  chr_i = gene_meta$V1[which(gene_meta$V4 == gene_i)]
  start_i = gene_meta$V2[which(gene_meta$V4 == gene_i)]
  end_i = gene_meta$V3[which(gene_meta$V4 == gene_i)]
  if (length(start_i) == 0) {
    tajima_all[i] = NA
    next
  }
  if (is.na(start_i)) {
    tajima_all[i] = NA
    next
  }
 
  tajima_i = tajima[tajima$V2 == chr_i, ]
  index_start_i = which.min(abs(tajima_i$V3 - start_i))
  index_end_i = which.min(abs(tajima_i$V3 - end_i))
  tajima_all[i] = mean(tajima_i$V5[index_start_i:index_end_i])
}
names(tajima_all) = all_gene
tajima_all = na.omit(tajima_all)
tajima_sig = tajima_all[uniq_gene]
boxplot(tajima_all, tajima_sig[1:6], tajima_sig,
        names = c('whole genome', 'target protein', 'target + regulator'),
        ylab = 'Tajima D')

## ABO (chr9) interactor p 
inter_p = c()
for (i in 1:22) {
  if (i == 9) {
    next
  }
  inter_p_i = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/sandwich_p/cis_chr', 
                             i, '_trans_chr9.rds'))
  abo_p_i = inter_p_i$inter_p[which(rownames(inter_p_i$inter_p) == 'ABO'), ]
  inter_p = c(inter_p, unlist(abo_p_i))
}

abo_target = data.frame(gene = names(inter_p), inter_p = inter_p)
abo_target$tajima = tajima_all[match(abo_target$gene, names(tajima_all))]
abo_target = na.omit(abo_target)

### the enrichment of significant ABO target in high tajima'D (> 1.5) under different interaction p cutoffs
tajima_cut = 1.5

boot_enrich = function(data, i){
  df = data[i, ]
  return(sum(df$tajima > tajima_cut) / nrow(df))
}

inter_cut_all = c(0.1, 0.05, 0.01, 0.001, 1.7e-7)
tajima_enrich = c()
enrich_low = c()
enrich_upper = c()
for (i in 1:length(inter_cut_all)) {
  sig_i = abo_target[abo_target$inter_p < inter_cut_all[i], ]
  not_sig_i = abo_target[abo_target$inter_p > inter_cut_all[i], ]
  sig_boot = boot(sig_i, boot_enrich, R = 1000)
  not_sig_boot = boot(not_sig_i, boot_enrich, R = 1000)
  enrich_i = sig_boot$t / not_sig_boot$t
  tajima_enrich[i] = mean(enrich_i)
  enrich_low[i] = quantile(enrich_i, 0.025)
  enrich_upper[i] = quantile(enrich_i, 0.975)
}

enrich_tab = data.frame(inter_cut = c('< 0.1', '< 0.05', '< 0.01', '< 0.001', 
                                      '< 1.7e-7'),
                        tajima_enrich = tajima_enrich,
                        enrich_low = enrich_low, 
                        enrich_upper = enrich_upper)
ggplot(enrich_tab, aes(x = reorder(inter_cut, tajima_enrich), 
                       y = tajima_enrich, color = inter_cut)) + 
  geom_point(stat="identity", size = 5) +
  geom_errorbar(position=position_dodge(0.5), 
                aes(ymin=enrich_low, ymax=enrich_upper), width=.1) + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  scale_color_manual(values = c(rep('steelblue', 4), 'red')) + 
  labs(x = "Interaction P cutoff", y = "Fold enrichment", 
       title = "ABO target with high Tajima's D") + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
