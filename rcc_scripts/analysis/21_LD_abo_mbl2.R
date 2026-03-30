setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(plink2R)
library(RColorBrewer)
library(pheatmap)
library(gaston)
## LD between MBL2 and ABO cis var
mbl2_cis = read_plink('cis_trans_inter/00_ref/cis_geno/MBL2',
                      impute="avg")
abo_cis = read_plink('cis_trans_inter/00_ref/cis_geno/ABO',
                     impute="avg")

load('cis_trans_inter/06_fusion_model/ABO.wgt.RDat')
abo_weight = wgt.matrix
abo_weight = abo_weight[abo_weight[,1]!=0, ]
abo_weight[,2] = abs(abo_weight[,1])
abo_weight = abo_weight[order(abo_weight[,2], decreasing = T), ]

load('cis_trans_inter/06_fusion_model/MBL2.wgt.RDat')
mbl2_weight = wgt.matrix
mbl2_weight = mbl2_weight[mbl2_weight[,1]!=0, ]
mbl2_weight[,2] = abs(mbl2_weight[,1])
mbl2_weight = mbl2_weight[order(mbl2_weight[,2], decreasing = T), ]

# extract top 20 SNPs of ABO and MBL2
abo_top_cis = abo_cis$bed[, colnames(abo_cis$bed) %in% rownames(abo_weight)[1:20]]
colnames(abo_top_cis) = paste0('ABO_', colnames(abo_top_cis))

mbl2_top_cis = mbl2_cis$bed[, colnames(mbl2_cis$bed) %in% rownames(mbl2_weight)[1:20]]
colnames(mbl2_top_cis) = paste0('MBL2_', colnames(mbl2_top_cis))

# make bed matrix based on two cis regions
abo_mbl2_bed = cbind(abo_top_cis, mbl2_top_cis)
abo_bim = abo_cis$bim[abo_cis$bim$V2 %in% rownames(abo_weight)[1:20], ]
abo_bim$V2 = paste0('ABO_', abo_bim$V2)
mbl2_bim = mbl2_cis$bim[mbl2_cis$bim$V2 %in% rownames(mbl2_weight)[1:20], ]
mbl2_bim$V2 = paste0('MBL2_', mbl2_bim$V2)
abo_mbl2_bim = rbind(abo_bim, mbl2_bim)
colnames(abo_cis$fam) = c('famid', 'id', 'father', 'mother', 'sex', 'pheno')
colnames(abo_mbl2_bim) = c('chr', 'id', 'dist', 'pos', 'A1', 'A2')


ld_abo_mbl2 = LD(as.bed.matrix(abo_mbl2_bed, abo_cis$fam, abo_mbl2_bim), 
                lim = c(1,ncol(abo_mbl2_bed)))
LD.plot(ld_abo_mbl2,
        graphical.par = list(mar = c(0,0,1,0)), 
        cex.ld = 0.4, cex.snp = 0.4)

### combination of ABO O allele and MBL2 allele
abo_mbl2 = data.frame(o_allele = abo_top_cis[, 'ABO_rs8176719'],
                      mbl2 = mbl2_top_cis[, 'MBL2_rs11003131'])
allele_freq = as.data.frame(table(abo_mbl2$o_allele, abo_mbl2$mbl2))
allele_freq = allele_freq[allele_freq$Var2 %in% 0:2, ]
allele_freq$Freq = allele_freq$Freq / sum(allele_freq$Freq)
allele_freq$ABO_MBL2 = paste0(allele_freq$Var1, ':', allele_freq$Var2)
ggplot(allele_freq, aes(x = reorder(ABO_MBL2, -Freq), y = Freq)) + 
  geom_bar(stat="identity", width=0.5, fill = 'steelblue') + 
  labs(x = "ABO:MBL2", y = 'Frequency', title = 'O allele and MBL2 allele comb', 
       color = '') + 
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())