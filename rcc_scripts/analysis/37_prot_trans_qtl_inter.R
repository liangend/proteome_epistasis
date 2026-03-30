library(data.table)
library(ggplot2)
library(boot)
setwd('/project/xuanyao/jinghui')
all_inter = fread('cis_trans_inter/20_prot_trans_pqtl_inter/inter_all.txt')
### enrichment of significant epistasis in PPI
gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.gtf')
gene_meta = gene_meta[gene_meta$gene_type == 'protein_coding', ]
near_gene = c()
for (i in 1:nrow(all_inter)) {
  chr_i = all_inter$chr[i]
  bp_i = all_inter$bp_hg38[i]
  gene_meta_i = gene_meta[gene_meta$chr == paste0('chr', chr_i), ]
  near_gene[i] = gene_meta_i$gene_name[which.min(abs(gene_meta_i$start - bp_i))]
}
all_inter$near_gene = near_gene

ppi_list = readRDS('pqtl/11_prot_complex/ppi_list.rds')
all_inter$is_ppi = paste0(all_inter$near_gene, ',', 
                          all_inter$prot) %in% ppi_list
all_inter$is_sig = all_inter$qval < 0.05

boot_mean = function(data, i){
  df = data[i, ]
  sig_ppi = sum(df$is_ppi & df$is_sig) / sum(df$is_sig)
  all_ppi = sum(df$is_ppi) / nrow(df)
  return(sig_ppi / all_ppi)
}

ppi_boot = boot(all_inter, boot_mean, R = 1000)
ppi_plot = data.frame(type = 'PPI',
                      enrich = mean(ppi_boot$t),
                      quant_lower = quantile(ppi_boot$t, 0.025),
                      quant_higher = quantile(ppi_boot$t, 0.975))
ggplot(ppi_plot, aes(x = type, y = enrich)) + 
  geom_bar(stat="identity", position=position_dodge(0.5), 
           width = 0.5, fill = 'steelblue') +
  geom_errorbar(position=position_dodge(0.5), 
                aes(ymin=quant_lower, ymax=quant_higher), width=.3) + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "", y = "Enrichment", fill = '',  title = 'Epistatic trans-pQTL in PPI') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## significant trans-pQTL * prot interactions
sig_inter = all_inter[all_inter$qval < 0.05, ]
sig_inter[50, ]
gene_meta = fread('mediation_pathway/00_ref/genecode.GRCh38.gene.meta.gtf')
gene_meta_sub = gene_meta[gene_meta$chr == 'chr20', ]
gene_meta_sub[which.min(abs(gene_meta_sub$start - 31223514)),]

## interaction with blood type
asso_file = list.files('cis_trans_inter/20_prot_trans_pqtl_inter/asso_pheno/')
inter_p = c()
for (i in 1:length(asso_file)) {
  asso_i = readRDS(paste0('cis_trans_inter/20_prot_trans_pqtl_inter/asso_pheno/', 
                          asso_file[i]))
  inter_p_i = sapply(asso_i, function(x){ifelse(is.null(nrow(x)), NA, x[4,6])})
  inter_p = rbind(inter_p, inter_p_i)
}
rownames(inter_p) = sub('.rds', '', asso_file)
inter_p = inter_p[!rownames(inter_p) %in% c('Hypertension', 'Snoring'), ]

min_p = apply(inter_p, 2, min)
