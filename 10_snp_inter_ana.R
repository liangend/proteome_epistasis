setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(plink2R)
library(boot)
library(RColorBrewer)

## SNP*SNP interactions on significant protein interactions
sig_inter_list = fread('cis_trans_inter/07_prot_inter/sig_prot_prot_inter.txt')
inter_snp_all = c()
for (i in 1:length(sig_inter_list)) {
  target_i = sig_inter_list$target[i]
  reg_i = sig_inter_list$regulator[i]
  snp_inter_i = readRDS(paste0('cis_trans_inter/17_snp_snp_inter/inter_w_add/',
                               target_i, '_', reg_i, '.rds'))
  inter_p_i = sapply(snp_inter_i, function(x){x[4,5]})
  inter_beta_i = sapply(snp_inter_i, function(x){x[4,1]})
  inter_se_i = inter_beta_i / qchisq(inter_p_i, 1, lower.tail = F)^0.5
  
  inter_i = data.frame(reg_prot = reg_i, target_prot = target_i,
                       reg_snp = sapply(strsplit(names(inter_p_i), '*', fixed = T), '[', 1),
                       target_snp = sapply(strsplit(names(inter_p_i), '*', fixed = T), '[', 2),
                       beta = unname(inter_beta_i), se = unname(inter_se_i),
                       p = unname(inter_p_i))
  
  ## find the fusion effect for target SNPs
  load(paste0('cis_trans_inter/06_fusion_model/ukb_eur/', target_i, '.wgt.RDat'))
  weight_i = as.data.frame(wgt.matrix)
  weight_i = weight_i[weight_i$enet !=0, ]
  inter_i$target_eff = weight_i$enet[match(inter_i$target_snp, rownames(weight_i))]
  inter_snp_all = rbind(inter_snp_all, inter_i)
}
snp_inter_q = qvalue(inter_snp_all$p)
inter_snp_all$qval = snp_inter_q$qvalues

inter_sig = inter_snp_all[inter_snp_all$qval < 0.05, ]
inter_not_sig = inter_snp_all[inter_snp_all$qval > 0.05, ]
# fwrite(inter_sig, 'cis_trans_inter/17_snp_snp_inter/snp_inter.txt', sep = '\t')

eff_dat = data.frame(beta = c(abs(inter_not_sig$target_eff), 
                              abs(inter_sig$target_eff)),
                     group = rep(c('non-significant SNPs', 'significant SNPs'), 
                                 times = c(nrow(inter_not_sig), nrow(inter_sig))))

ggplot(eff_dat, aes(x = group, y = beta, fill = group)) + 
  geom_boxplot(width = 0.5) +
  labs(title = 'SNPs with epistasis have a larger effect', x = "", y = "|β|", fill = '') +
  scale_fill_manual(values = brewer.pal(8, "Dark2")) + 
  ylim(c(0,0.3)) + 
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
wilcox.test(abs(inter_not_sig$target_eff), 
            abs(inter_sig$target_eff))

# make SNP vcf bed file
inter_sig_full = c()
for (i in 1:length(sig_inter_list)) {
  target_i = sig_inter_list$target[i]
  reg_i = sig_inter_list$regulator[i]
  
  inter_sig_i = inter_sig[inter_sig$reg_prot == reg_i &
                            inter_sig$target_prot == target_i, ]
  
  reg_geno = read_plink(paste0('cis_trans_inter/00_ref/cis_geno/', reg_i))
  reg_bim = reg_geno$bim
  inter_sig_i = cbind(inter_sig_i, 
                      reg_bim[match(inter_sig_i$reg_snp, 
                                    reg_bim$V2), c(1, 4:6)])

  target_geno = read_plink(paste0('cis_trans_inter/00_ref/cis_geno/', target_i))
  target_bim = target_geno$bim
  inter_sig_i = cbind(inter_sig_i, 
                      target_bim[match(inter_sig_i$target_snp, 
                                       target_bim$V2), c(1, 4:6)])
  
  inter_sig_full = rbind(inter_sig_full, inter_sig_i)
}
colnames(inter_sig_full)[8:15] = c('reg_chr', 'reg_bp', 'reg_ref', 'reg_alt',
                                   'target_chr', 'target_bp', 
                                   'target_ref', 'target_alt')

target_snp_vcf = data.frame(CHROM = inter_sig_full$target_chr, 
                            POS = inter_sig_full$target_bp,
                            ID = '.',
                            REF = inter_sig_full$target_ref, 
                            ALT = inter_sig_full$target_alt,
                            QUAL = '.', FILTER = '.', INFO = '.')
reg_snp_vcf = data.frame(CHROM = inter_sig_full$reg_chr, 
                         POS = inter_sig_full$reg_bp,
                         ID = '.',
                         REF = inter_sig_full$reg_ref, 
                         ALT = inter_sig_full$reg_alt,
                         QUAL = '.', FILTER = '.', INFO = '.')
colnames(target_snp_vcf)[1] = '#CHROM'
colnames(reg_snp_vcf)[1] = '#CHROM'
target_snp_vcf = unique(target_snp_vcf)
reg_snp_vcf = unique(reg_snp_vcf)
fwrite(target_snp_vcf,
       'cis_trans_inter/18_sig_inter_snp_annot/target_snp_all.vcf', sep = '\t')
fwrite(reg_snp_vcf,
       'cis_trans_inter/18_sig_inter_snp_annot/reg_snp_all.vcf', sep = '\t')

target_snp_bed = inter_sig_full[, c('target_chr','target_bp','target_bp','target_eff')]
reg_snp_bed = inter_sig_full[, c('reg_chr','reg_bp','reg_bp')]
target_snp_bed$target_chr = paste0('chr', target_snp_bed$target_chr)
reg_snp_bed$reg_chr = paste0('chr', reg_snp_bed$reg_chr)
fwrite(target_snp_bed, 'cis_trans_inter/18_sig_inter_snp_annot/target_snp_all.bed',
       col.names = F, sep = '\t')
fwrite(reg_snp_bed, 'cis_trans_inter/18_sig_inter_snp_annot/reg_snp_all.bed',
       col.names = F, sep = '\t')

### SNP chromatin state after matching the effect size
inter_snp_state = fread('cis_trans_inter/18_sig_inter_snp_annot/target_snp_all_w_state.bed')
inter_snp_state$fusion_eff = inter_sig$target_eff
inter_snp_state = unique(inter_snp_state)
inter_snp_state = na.omit(inter_snp_state)

fusion_snp_state = fread('cis_trans_inter/13_prot_snp_inter/fusion_snp_w_state.bed')
fusion_snp_state = unique(fusion_snp_state)

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
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(l = 40, r = 10, t = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




