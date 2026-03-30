setwd('/project/xuanyao/jinghui')
library(data.table)
library(qvalue)
library(boot)
library(ggplot2)
library(RColorBrewer)
### significant regulator pred by target snp interaction effect
sig_inter = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')
snp_prot_inter = c()
for (i in 1:nrow(sig_inter)) {
  target_i = sig_inter$target[i]
  reg_i = sig_inter$regulator[i]
  inter_i = fread(paste0('cis_trans_inter/13_prot_snp_inter/', target_i,
                         '_', reg_i, '.txt'))
  ## find the fusion effect for target SNPs
  load(paste0('cis_trans_inter/06_fusion_model/ukb_eur/', target_i, '.wgt.RDat'))
  weight_i = as.data.frame(wgt.matrix)
  weight_i = weight_i[weight_i$enet !=0, ]
  inter_i$fusion_eff = weight_i$enet[match(inter_i$snp, rownames(weight_i))]
  
  inter_i$pair = paste0(target_i, '_', reg_i)
  snp_prot_inter = rbind(snp_prot_inter, inter_i)
}

snp_prot_inter_q = qvalue(snp_prot_inter$reg.snp_sand_p)
snp_prot_inter_qval = snp_prot_inter_q$qvalues
snp_prot_inter$snp_prot_inter_qval = snp_prot_inter_qval

snp_prot_inter_sig = snp_prot_inter[snp_prot_inter$snp_prot_inter_qval < 0.05, ]

eff_dat = data.frame(beta = c(abs(snp_prot_inter$fusion_eff), 
                              abs(snp_prot_inter_sig$fusion_eff)),
                     group = rep(c('all SNPs', 'significant SNPs'), 
                                 times = c(nrow(snp_prot_inter), nrow(snp_prot_inter_sig))))
ggplot(eff_dat, aes(x = group, y = beta, fill = group)) + 
  geom_boxplot(width = 0.5) +
  labs(title = 'FUSION SNP effect', x = "", y = "|β|", fill = '') +
  scale_fill_manual(values = brewer.pal(8, "Dark2")) + 
  ylim(c(0,0.2)) + 
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

### LD pruning
# inter_sig_ld_prun = c()
# for (i in 1:nrow(sig_inter)) {
#   target_i = sig_inter$target[i]
#   reg_i = sig_inter$regulator[i]
#   inter_sub = snp_prot_inter_sig[snp_prot_inter_sig$pair == paste0( target_i, '_', reg_i), ]
#   ld_i = readRDS(paste0('cis_trans_inter/13_prot_snp_inter/', target_i, '_ld.rds'))
#   while(nrow(inter_sub) > 0) {
#     top_inter = which.min(inter_sub$snp_prot_inter_qval)[1]
#     top_snp = inter_sub$snp[top_inter]
#     top_ld = (ld_i[top_snp, inter_sub$snp])^2
#     rm_snp = names(which(top_ld > 0.1))
#     inter_sig_ld_prun = rbind(inter_sig_ld_prun, inter_sub[top_inter, ])
#     inter_sub = inter_sub[!(inter_sub$snp %in% c(rm_snp, top_snp)), ]
#   }
#   print(i)
# }

## make vcf input for snpEff
# cis_snp = fread('cis_trans_inter/00_ref/cis_geno/cis_snp_all.bim')
# inter_snp = data.frame(chr = cis_snp$V1[match(snp_prot_inter_sig$snp, cis_snp$V2)],
#                        POS = cis_snp$V4[match(snp_prot_inter_sig$snp, cis_snp$V2)],
#                        ID = '.',
#                        REF = cis_snp$V5[match(snp_prot_inter_sig$snp, cis_snp$V2)],
#                        ALT = cis_snp$V6[match(snp_prot_inter_sig$snp, cis_snp$V2)],
#                        QUAL = '.', FILTER = '.', INFO ='.')
# colnames(inter_snp)[1] = '#CHROM'
# fwrite(inter_snp, 'cis_trans_inter/13_prot_snp_inter/inter_sig_snp.vcf',
#        sep = '\t')
# 
# fusion_snp = data.frame(chr = cis_snp$V1[match(snp_prot_inter$snp, cis_snp$V2)],
#                        POS = cis_snp$V4[match(snp_prot_inter$snp, cis_snp$V2)],
#                        ID = '.',
#                        REF = cis_snp$V5[match(snp_prot_inter$snp, cis_snp$V2)],
#                        ALT = cis_snp$V6[match(snp_prot_inter$snp, cis_snp$V2)],
#                        QUAL = '.', FILTER = '.', INFO ='.')
# colnames(fusion_snp)[1] = '#CHROM'
# fwrite(fusion_snp, 'cis_trans_inter/13_prot_snp_inter/fusion_snp.vcf',
#        sep = '\t')

### make snp bed files
# cis_snp = fread('cis_trans_inter/00_ref/cis_geno/cis_snp_all.bim')
# inter_snp = fread('cis_trans_inter/13_prot_snp_inter/inter_sig_snp.vcf')
# inter_snp_bed = data.frame(chr = paste0('chr', inter_snp$`#CHROM`),
#                            start = inter_snp$POS, end = inter_snp$POS)
# fusion_snp = fread('cis_trans_inter/13_prot_snp_inter/fusion_snp.vcf')
# fusion_snp_bed = data.frame(chr = paste0('chr', fusion_snp$`#CHROM`),
#                            start = fusion_snp$POS, end = fusion_snp$POS)
# fwrite(inter_snp_bed, 'cis_trans_inter/13_prot_snp_inter/inter_sig_snp.bed',
#        col.names = F, sep = '\t')
# fwrite(fusion_snp_bed, 'cis_trans_inter/13_prot_snp_inter/fusion_snp.bed',
#        col.names = F, sep = '\t')


## snp annotation
inter_sig_annot = fread('cis_trans_inter/13_prot_snp_inter/inter_sig_snp_annot.vcf')
inter_sig_annot$annot = sapply(strsplit(inter_sig_annot$INFO, '|', fixed = T), '[', 2)
inter_sig_annot_tab = as.data.frame(table(inter_sig_annot$annot) / nrow(inter_sig_annot))

fusion_annot = fread('cis_trans_inter/13_prot_snp_inter/fusion_snp_annot.vcf')
fusion_annot$annot = sapply(strsplit(fusion_annot$INFO, '|', fixed = T), '[', 2)
fusion_annot_tab = as.data.frame(table(fusion_annot$annot) / nrow(fusion_annot))

inter_sig_annot_tab$fusion_freq = fusion_annot_tab$Freq[match(inter_sig_annot_tab$Var1,
                                                              fusion_annot_tab$Var1)]
annot_cat = as.character(inter_sig_annot_tab$Var1)
boot_annot = function(data, i){
  df = data[i, ]
  annot_freq = table(df$annot)/nrow(df)
  annot_freq = annot_freq[match(annot_cat, names(annot_freq))]
}
inter_sig_boot = boot(inter_sig_annot, boot_annot, R = 1000)
inter_sig_boot$t[is.na(inter_sig_boot$t)] = 0
trans_annot_enrich = t(apply(inter_sig_boot$t, 1, 
                             function(x){x/inter_sig_annot_tab$fusion_freq}))
colnames(trans_annot_enrich) = annot_cat

inter_sig_annot_tab$enrich = inter_sig_annot_tab$Freq / inter_sig_annot_tab$fusion_freq
inter_sig_annot_tab$ci_low = apply(trans_annot_enrich, 2, function(x){quantile(x, 0.05)})
inter_sig_annot_tab$ci_high = apply(trans_annot_enrich, 2, function(x){quantile(x, 0.95)})
ggplot(inter_sig_annot_tab, aes(x = reorder(Var1, enrich), y = enrich)) + 
  geom_point(stat="identity", size = 2, color = 'black') +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.2) +
  geom_hline(yintercept = 1, lty = 2) + 
  coord_flip() + 
  labs(title = '', x = "", y = "Fold enrichment of interacting SNP vs cis-pQTL") +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### SNP chromatin state
inter_snp_state = fread('cis_trans_inter/13_prot_snp_inter/inter_sig_snp_w_state.bed')
fusion_snp_state = fread('cis_trans_inter/13_prot_snp_inter/fusion_snp_w_state.bed')

inter_snp_state$fusion_eff = snp_prot_inter_sig$fusion_eff
fusion_snp_state$fusion_eff = snp_prot_inter$fusion_eff
inter_snp_state = na.omit(inter_snp_state)
background_prop = matrix(0, nrow = 1000, ncol = nrow(inter_snp_state))
for (i in 1:nrow(inter_snp_state)) {
  eff_i = abs(inter_snp_state$fusion_eff[i])
  fusion_snp_sub = fusion_snp_state[which(abs(fusion_snp_state$fusion_eff) > 0.9 * eff_i &
                                            abs(fusion_snp_state$fusion_eff) < 1.1 * eff_i), ]
  background_prop[, i] = sample(fusion_snp_sub$V7, 1000, replace = T)
}
inter_snp_prop = as.data.frame(table(inter_snp_state$V7)/nrow(inter_snp_state))
obs_count = as.data.frame(table(inter_snp_state$V7))
inter_snp_enrich = t(apply(background_prop, 1, function(x){
  background_i = as.data.frame(table(x))
  count_i = background_i$Freq[match(obs_count$Var1, background_i$x)]
  return(log2(obs_count$Freq/count_i))
}))


# boot_background = function(data, i){
#   df = data[i, ]
#   prop_all = as.data.frame(table(df$V7)/nrow(df))
#   return(prop_all$Freq[match(inter_snp_prop$Var1, prop_all$Var1)])
# }
# fusion_prop = boot(fusion_snp_state, boot_background, R = 1000)
# 
# fusion_prop = fusion_prop$t
# inter_snp_enrich = t(apply(fusion_prop, 1, function(x){log2(inter_snp_prop$Freq/x)}))

enrich_plt = data.frame(enrich = apply(inter_snp_enrich, 2, function(x){mean(x, na.rm = T)}),
                        quant_lower = apply(inter_snp_enrich, 2, function(x){quantile(x, 0.025, na.rm = T)}),
                        quant_higher = apply(inter_snp_enrich, 2, function(x){quantile(x, 0.975, na.rm = T)}),
                        state = obs_count$Var1)

state_annot = fread('pqtl/07_chr_state/roadmap_25_annot.txt')
enrich_plt$annot = state_annot$DESCRIPTION[match(enrich_plt$state, state_annot$`STATE NO.`)]
ggplot(enrich_plt, aes(x = reorder(annot, -enrich), y = enrich)) + 
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), 
                width=.2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "", y = expression(paste("log"[2], "fold enrichment")), 
       title = expression(paste('Functional enrichment of ', italic('cis'), 
                                '-SNPs with significant interactions'))) +
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())





