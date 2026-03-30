setwd('/project/xuanyao/jinghui')
library(data.table)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(plink2R)

### Tajima's D for each gene
tajima = fread('cis_trans_inter/21_tajima_D/tajdEd.txt.gz')
gene_meta = fread('mediation_pathway/00_ref/genecode.GRCh38.gene.meta.gtf')
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')
all_gene = list.files('cis_trans_inter/06_fusion_pred/eur/', pattern = 'good_pred.txt')
all_gene = sub('_good_pred.txt', '', all_gene)
uniq_gene = unique(c(sig_inter$target, sig_inter$regulator))
tajima_all = c()
for (i in 1:length(all_gene)) {
  gene_i = all_gene[i]
  chr_i = gene_meta$chr[which(gene_meta$gene_name == gene_i)]
  start_i = gene_meta$start_hg37[which(gene_meta$gene_name == gene_i)]
  end_i = gene_meta$end_hg37[which(gene_meta$gene_name == gene_i)]
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

# boot_enrich = function(data, i, inter_cut){
#   df = data[i, ]
#   sig_inter = sum(df$inter_p < inter_cut & df$tajima > tajima_cut) /
#     sum(df$inter_p < inter_cut)
#   not_inter = sum(df$inter_p > inter_cut & df$tajima > tajima_cut) / 
#     sum(df$inter_p > inter_cut)
#   return(sig_inter / not_inter)
# }

# for (i in 1:length(inter_cut_all)) {
#   tajima_boot = boot(abo_target, boot_enrich, R = 1000, 
#                      inter_cut = inter_cut_all[i])
#   enrich_i = tajima_boot$t
#   tajima_enrich[i] = mean(tajima_boot$t)
#   enrich_low[i] = quantile(tajima_boot$t, 0.025)
#   enrich_upper[i] = quantile(tajima_boot$t, 0.975)
# }

boot_enrich = function(data, i){
  df = data[i, ]
  return(sum(df$tajima > tajima_cut) / nrow(df))
}

inter_cut_all = c(0.1, 0.05, 0.01, 0.001, 1.7e-9)
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
                                      '< 1.7e-9'),
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

### examples of ABO target showing different variance under different O alleles
prot = 'NPTX1'
prot_expr = fread(paste0('cis_trans_inter/06_fusion_pred/eur/', 
                         prot, '_good_pred.txt'))

abo_type = fread('cis_trans_inter/data/ukb_blood_type.pheno')
colnames(abo_type)[3] = 'blood_type'

prot_expr$blood_type = abo_type$blood_type[match(prot_expr$id, 
                                                 abo_type$FID)]
prot_expr = na.omit(prot_expr)

var.test(prot_expr$obs[prot_expr$blood_type != 'OO'], 
         prot_expr$obs[prot_expr$blood_type == 'OO'])

## variance of protein expression level changing across ABO alleles
abo_cis = read_plink('cis_trans_inter/00_ref/cis_geno/ABO',
                     impute="avg")
prot_expr$o_allele = abo_cis$bed[match(prot_expr$id, abo_cis$fam$V1),
                                 'rs8176719']
prot_expr$o_allele = as.factor(prot_expr$o_allele)

boot_var_fun = function(formula, data, i){
  df = data[i, ]
  aggre_var = aggregate(formula, data = df, var)
  return(c(var(df[,3]), aggre_var[,2]))
}
var_boot2 = boot(prot_expr, boot_var_fun, 1000, formula = obs ~ o_allele)
var_dat2 = data.frame(var = apply(var_boot2$t, 2, mean),
                      ci_low = apply(var_boot2$t, 2, function(x){quantile(x, 0.025)}),
                      ci_high = apply(var_boot2$t, 2, function(x){quantile(x, 0.975)}),
                      group = c('all', '0', '1', '2'))
ggplot(var_dat2[var_dat2$group != 'all', ], aes(x = group, y = var)) + 
  geom_point(position = position_dodge(0.3), size = 5,  color = 'steelblue') +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .1, 
                position=position_dodge(0.3), color = 'steelblue') +
  labs(x = "rs8176719 (ABO) dosage", y = 'Variance', 
       title = paste0(prot, ' variation across O alleles')) + 
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


## variance of protein expression level changing across different blood type 
# var_boot = boot(prot_expr, boot_var_fun, 1000, formula = obs ~ blood_type)
# var_dat = data.frame(var = apply(var_boot$t, 2, mean),
#                      ci_low = apply(var_boot$t, 2, function(x){quantile(x, 0.025)}),
#                      ci_high = apply(var_boot$t, 2, function(x){quantile(x, 0.975)}),
#                      group = c('all', 'AA', 'AB', 'AO', 'BB', 'BO', 'OO'))
# ggplot(var_dat, aes(x = factor(group, levels = c('all', 'AA', 'AB', 'AO', 
#                                                  'BB', 'BO', 'OO')), 
#                     y = var, 
#                     color = factor(group, levels = c('all', 'AA', 'AB', 'AO', 
#                                                      'BB', 'BO', 'OO')))) + 
#   geom_point(position = position_dodge(0.3), size = 5) +
#   geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .1, 
#                 position=position_dodge(0.3)) +
#   labs(x = "blood type", y = 'Var', 
#        title = paste0('Var of ', prot, ' across blood types'), 
#        color = '') + 
#   scale_color_manual(values = c('steelblue', rep("black", 6))) +
#   theme(text = element_text(size=13, colour = "black"), 
#         axis.text.x = element_text(colour = "black", size = 13),
#         axis.text.y = element_text(colour = "black", size = 13),
#         axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         legend.position = 'none')

