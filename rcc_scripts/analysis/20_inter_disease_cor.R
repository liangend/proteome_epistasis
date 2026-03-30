setwd('/project/xuanyao/jinghui')
library(data.table)
library(plink2R)
library(ggplot2)
library(ivreg)
library(gaston)
library(MendelianRandomization)

mbl2_expr = fread('cis_trans_inter/06_fusion_pred/eur/MBL2_good_pred.txt')
abo_expr = fread('cis_trans_inter/06_fusion_pred/eur/ABO_good_pred.txt')
pred_tab = data.frame(ID = intersect(mbl2_expr$id, abo_expr$id))
pred_tab$abo_pred = abo_expr$pred[match(pred_tab$ID, abo_expr$id)]
pred_tab$mbl2_pred = mbl2_expr$pred[match(pred_tab$ID, mbl2_expr$id)]
pred_tab = na.omit(pred_tab)

## trait: oral_trait/oral_issue674206, Rheumatoid_arthritis_pain674178, HDL674178, LDL674178
trait = ''
## oral ulcer
ulcer = fread('cis_trans_inter/data/oral_trait/oral_issue674206.pheno')
ulcer = ulcer[ulcer$`6149` %in% c(1, -7)]
ulcer$`6149`[ulcer$`6149` == -7] = 0

## RA
ra = fread('cis_trans_inter/data/Rheumatoid_arthritis_pain674178.pheno')
ra = ra[ra$`120001` %in% c(1, 0)]

## LDL
ldl = fread('cis_trans_inter/data/LDL674178.pheno')

## HDL
hdl = fread('cis_trans_inter/data/HDL674178.pheno')

pred_tab$ulcer = ulcer$`6149`[match(pred_tab$ID, ulcer$FID)]
pred_tab$ra = ra$`120001`[match(pred_tab$ID, ra$FID)]
pred_tab$ldl = ldl$`30780`[match(pred_tab$ID, ldl$FID)]
pred_tab$hdl = hdl$`30760`[match(pred_tab$ID, hdl$FID)]
pred_tab$inter = pred_tab$abo_pred * pred_tab$mbl2_pred

summary(glm(pred_tab$ulcer ~ pred_tab$inter, family = binomial))
summary(glm(pred_tab$ra ~ pred_tab$inter, family = binomial))

summary(lm(pred_tab$ldl ~ pred_tab$abo_pred))
summary(lm(pred_tab$ldl ~ pred_tab$mbl2_pred))
summary(lm(pred_tab$ldl ~ pred_tab$inter))
cor(pred_tab$ldl, pred_tab$inter, use = 'na.or.complete')
plot(pred_tab$ldl, pred_tab$inter)
pred_tab$inter_bin = cut(pred_tab$inter, breaks = 5)
boxplot(ldl ~ inter_bin, data = pred_tab, xlab = 'MBL2 * ABO', y = 'LDL')





