setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(plink2R)
library(RColorBrewer)

ld_indep_inter = readRDS('cis_trans_inter/17_snp_snp_inter/inter_w_add_indep.rds')

## UKB protein level
load('cis_trans_inter/data/ukb_ppp_prot_expr.rdata')
colnames(prot_ex) = gene_name

ukb_cov = fread('cis_trans_inter/06_fusion_pred/all_eur/all_norm_cov.txt')
ukb_cov = na.omit(ukb_cov)

names(ld_indep_inter)
inter_i = 1 # CD209*ABO
inter_i = 8 # MBL2*ABO
inter_sig_i = ld_indep_inter[[inter_i]]
target_i = unlist(strsplit(names(ld_indep_inter)[inter_i], '_'))[1]
reg_i = unlist(strsplit(names(ld_indep_inter)[inter_i], '_'))[2]
inter_snp_pair = names(ld_indep_inter[[inter_i]])
target_snp = sapply(strsplit(inter_snp_pair, '*', fixed = T), '[', 2)
reg_snp = sapply(strsplit(inter_snp_pair, '*', fixed = T), '[', 1)

## target geno
target_geno = read_plink(paste0('cis_trans_inter/00_ref/cis_geno/', target_i))
geno_id = target_geno$fam$V1
target_bed = target_geno$bed

## regulator geno
reg_geno = read_plink(paste0('cis_trans_inter/00_ref/cis_geno/', reg_i))
reg_bed = reg_geno$bed

## asso with target protein
prot_target = prot_ex[, target_i]
names(prot_target) = sample_id
prot_target = na.omit(prot_target)

common_id = Reduce(intersect, list(names(prot_target), geno_id, ukb_cov$FID))

prot_sub = prot_target[match(common_id, names(prot_target))]
ukb_cov_sub = ukb_cov[match(common_id, ukb_cov$FID), -(1:2)]
prot_res = residuals(lm(prot_sub ~ as.matrix(ukb_cov_sub)))

## specific inter
extr_pair = names(which.min(inter_sig_i)) # MBL2*ABO
extr_pair = names((inter_sig_i)[2]) # CD209*ABO
extr_target_snp = unlist(strsplit(extr_pair, '*', fixed = T))[2]
extr_reg_snp = unlist(strsplit(extr_pair, '*', fixed = T))[1]
extr_reg_snp = 'rs8176719' # snp determining O blood type

target_snp1 = target_bed[match(common_id, geno_id), extr_target_snp]
reg_snp1 = reg_bed[match(common_id, geno_id), extr_reg_snp]

fit1 = lm(prot_res~target_snp1*reg_snp1)
summary(fit1)

prot_dat = data.frame(prot = prot_sub, target_snp = target_snp1,
                      reg_snp = reg_snp1)
prot_dat = na.omit(prot_dat)

### CD209*ABO
ggplot(prot_dat, aes(x = as.factor(target_snp), y = prot, 
                     fill = as.factor(reg_snp))) + 
  geom_boxplot() + 
  geom_smooth(aes(group = 1), method='lm', color = 'red', 
              linetype = 2, se = F) +
  facet_wrap(~ reg_snp, nrow = 1) + ylim(c(-1.7, 2.2)) + 
  labs(title = paste0(extr_reg_snp, ' (', reg_i, ')'), 
       x = paste0(extr_target_snp, ' (', target_i, ')'), 
       y = paste0('Residualized ', target_i)) +
  scale_fill_manual(values = brewer.pal(8, "Blues")[4:6]) +
  theme(text = element_text(size=13, colour = "black"), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(color="black", fill="bisque"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

fit1 = lm(prot~target_snp, data = prot_dat[prot_dat$reg_snp == 2, ])
summary(fit1)$coefficients

### MBL2*ABO
ggplot(prot_dat, aes(x = as.factor(reg_snp), y = prot, 
                     fill = as.factor(target_snp))) + 
  geom_boxplot() +
  geom_smooth(aes(group = 1), method='lm', color = 'red', 
              linetype = 2, se = F) +
  facet_wrap(~ target_snp, nrow = 1) + ylim(c(-5.2, 3.9)) + 
  labs(title = paste0(extr_target_snp, ' (', target_i, ')'), 
       x = paste0(extr_reg_snp, ' (', reg_i, ')'), 
       y = paste0('Residualized ', target_i)) +
  scale_fill_manual(values = brewer.pal(8, "Blues")[4:6]) +
  theme(text = element_text(size=13, colour = "black"), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(color="black", fill="bisque"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
fit1 = lm(prot~reg_snp, data = prot_dat[prot_dat$target_snp == 2, ])
summary(fit1)$coefficients


