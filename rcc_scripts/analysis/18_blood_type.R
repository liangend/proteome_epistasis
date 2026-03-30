setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(plink2R)
library(sandwich)
library(RColorBrewer)
library(gaston)
library(gridExtra)

## h2 and R2 of ABO and MBL2
pos_gene = fread('cis_trans_inter/00_ref/ukb_ppp_gene_pos.txt')
prot_h2 = fread('cis_trans_inter/08_prot_h2/cis_h2_summ.txt')

## Elastic net of ABO and MBL2
load('cis_trans_inter/06_fusion_model/ABO.wgt.RDat')
enet_abo = as.data.frame(wgt.matrix)
enet_abo = enet_abo[enet_abo$enet!=0, ]
enet_abo$top1 = abs(enet_abo$enet)
colnames(enet_abo)[2] = 'abs_coef'

load('cis_trans_inter/06_fusion_model/MBL2.wgt.RDat')
enet_mbl2 = as.data.frame(wgt.matrix)
enet_mbl2 = enet_mbl2[enet_mbl2$enet!=0, ]
enet_mbl2$top1 = abs(enet_mbl2$enet)
colnames(enet_mbl2)[2] = 'abs_coef'


## Plot the LD matrix of ABO cis non-zero SNPs
abo_cis = read_plink('cis_trans_inter/00_ref/cis_geno/ABO', impute="avg")
colnames(abo_cis$fam) = c('famid', 'id', 'father', 'mother', 'sex', 'pheno')
colnames(abo_cis$bim) = c('chr', 'id', 'dist', 'pos', 'A1', 'A2')
ld_abo_cis = LD(as.bed.matrix(abo_cis$bed, abo_cis$fam, abo_cis$bim), 
                c(1,ncol(abo_cis$bed)))
LD.plot(ld_abo_cis,
        graphical.par = list(mar = c(0,0,1,0)), 
        cex.ld = 0.2, cex.snp = 0.4)
barplot(enet_abo$enet, ylab = 'Elastic net weight')

cor(abo_cis$bed[,'rs505922'], abo_cis$bed[,'rs8176719'])

## ABO raw expr by blood type
load('cis_trans_inter/data/ukb_ppp_prot_expr.rdata')
blood_type = fread('cis_trans_inter/data/ukb_blood_type.pheno')
white_id = fread('cis_trans_inter/00_ref/white_id.txt')
colnames(prot_ex) = gene_name
abo_ex = data.frame(ID = sample_id, ABO = prot_ex[, 'ABO'])
abo_ex$blood_type = blood_type$`23165`[match(abo_ex$ID, blood_type$FID)]
abo_ex = abo_ex[abo_ex$ID %in% white_id$FID, ]
abo_ex$MBL2 = prot_ex[match(abo_ex$ID, sample_id), 'MBL2']

p1 = ggplot(abo_ex[!is.na(abo_ex$blood_type), ], 
       aes(x=reorder(blood_type, -ABO), y=ABO, fill=blood_type)) + 
  geom_boxplot() + 
  labs(x = "Blood type haplotype", y = 'ABO', 
       title = 'Normalized expression', color = '') + 
  scale_fill_manual(values = brewer.pal(8, "Dark2")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

## ABO raw expr by cis-pQTLs
abo_ex$rs8176719 = abo_cis$bed[match(abo_ex$ID, abo_cis$fam$id), 'rs8176719']
abo_ex$rs687289 = abo_cis$bed[match(abo_ex$ID, abo_cis$fam$id), 'rs687289']
abo_ex$rs527210 = abo_cis$bed[match(abo_ex$ID, abo_cis$fam$id), 'rs527210']
abo_ex$rs505922 = abo_cis$bed[match(abo_ex$ID, abo_cis$fam$id), 'rs505922']
cor(abo_cis$bed[, c('rs8176719', 'rs687289', 'rs527210', 'rs505922')])

ggplot(abo_ex[abo_ex$rs505922 %in% 0:2, ], aes(x=as.factor(rs505922), y=ABO)) + 
  geom_boxplot() + 
  labs(x = "SNP dosage", y = 'ABO prot expr', 
       title = 'rs505922', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

## ABO cis pred expr by blood type
abo_pred_tab = fread('cis_trans_inter/06_fusion_pred/eur/ABO_good_pred.txt')
abo_pred_tab$pred = (abo_pred_tab$pred - mean(abo_pred_tab$pred))/
  sd(abo_pred_tab$pred)
abo_pred_tab$blood_type = blood_type$`23165`[match(abo_pred_tab$id, blood_type$FID)]
abo_pred_tab$abo_cis = abo_cis$bed[match(abo_pred_tab$id, abo_cis$fam$id), 'rs505922']

p2= ggplot(abo_pred_tab[!is.na(abo_pred_tab$blood_type), ], 
       aes(x=reorder(blood_type, -pred), y=pred, fill=blood_type)) + 
  geom_boxplot() + 
  labs(x = "Blood type haplotype", y = 'ABO', 
       title = 'Predicted expression', color = '') + 
  scale_fill_manual(values = brewer.pal(8, "Dark2")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
grid.arrange(p1, p2, nrow = 1)

## raw and predicted expression together
abo_ex_all = na.omit(abo_pred_tab)
abo_ex_all = data.frame(abo = c(abo_ex_all$pred, abo_ex_all$obs),
                        group = rep(c('Predicted', 'Raw'), each = nrow(abo_ex_all)),
                        blood = rep(abo_ex_all$blood_type, 2))
ggplot(abo_ex_all, 
       aes(x=reorder(blood, -abo), y=abo, fill=group)) + 
  geom_boxplot() + 
  labs(x = "Blood type haplotype", y = 'ABO expression', 
       title = '', fill = '') + 
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(2,8)]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## MBL2 expression vs. ABO expression
abo_ex$abo_pred = abo_pred_tab$pred[match(abo_ex$ID, abo_pred_tab$id)]
ggplot(abo_ex[!is.na(abo_ex$blood_type), ], 
       aes(x=abo_pred, y=MBL2, color=blood_type)) + 
  geom_point() + 
  labs(x = "ABO cis pred", y = 'MBL2 expr', 
       title = '', color = 'Blood type') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
cor(abo_ex$MBL2, abo_ex$abo_pred, use = 'na.or.complete')

ggplot(abo_ex[abo_ex$blood_type == 'OO', ], 
       aes(x=abo_pred, y=MBL2, color=blood_type)) + 
  geom_point() + 
  labs(x = "ABO cis pred", y = 'MBL2 expr', 
       title = 'OO type only', color = 'Blood type') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = 'none')
cor.test(abo_ex[abo_ex$blood_type == 'OO', ]$MBL2, 
    abo_ex[abo_ex$blood_type == 'OO', ]$abo_pred, use = 'na.or.complete')

summary(lm(MBL2 ~ abo_pred * blood_type, data = abo_ex))


## MBL2 vs MBL2_pred by blood type
mbl2_tab = fread('cis_trans_inter/06_fusion_pred/eur/MBL2_good_pred.txt')
mbl2_tab$blood_type = blood_type$`23165`[match(mbl2_tab$id, blood_type$FID)]
mbl2_tab = na.omit(mbl2_tab)
mbl2_tab$abo_pred = abo_pred_tab$pred[match(mbl2_tab$id, abo_pred_tab$id)]
mbl2_tab$pred_interval = cut(mbl2_tab$pred, 7)
ggplot(mbl2_tab, 
       aes(x=pred, y=obs, color=blood_type)) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(method='lm') +
  labs(x = "Predicted MBL2", y = 'Observed MBL2', 
       title = '', color = '') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(mbl2_tab, 
       aes(x=pred, y=obs, color=pred_interval)) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(method='lm') +
  labs(x = "Predicted MBL2", y = 'Observed MBL2', 
       title = '', color = '') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
summary(lm(obs ~ pred + pred_interval, data = mbl2_tab))

ggplot(mbl2_tab, 
       aes(x=pred, y=obs, shape=pred_interval, color = blood_type)) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(method='lm') +
  labs(x = "Predicted MBL2", y = 'Observed MBL2', 
       title = '') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
summary(lm(obs ~ pred + pred_interval + blood_type, data = mbl2_tab))

ggplot(mbl2_tab, 
       aes(x=abo_pred, y=obs, color=pred_interval)) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(method='lm') +
  labs(x = "Predicted MBL2", y = 'Observed MBL2', 
       title = '', color = '') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(mbl2_tab, 
       aes(x=blood_type, y=obs, color=blood_type)) + 
  geom_boxplot() + 
  labs(x = "Blood type", y = 'Observed MBL2', 
       title = '', color = '') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
summary(lm(obs ~ blood_type, data = mbl2_tab))

## MBL2 vs MBL2_pred by top cis-pQTL of ABO
mbl2_tab$rs505922 = abo_cis$bed[match(mbl2_tab$id, abo_cis$fam$id), 'rs505922']
ggplot(mbl2_tab[!is.na(mbl2_tab$rs505922), ], 
       aes(x=pred, y=obs, color=as.factor(rs505922))) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(method='lm') +
  labs(x = "Predicted MBL2", y = 'Observed MBL2', 
       title = '', color = 'rs505922') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggplot(mbl2_tab[!is.na(mbl2_tab$rs505922), ], 
       aes(x=as.factor(rs505922), y=obs, color=as.factor(rs505922))) + 
  geom_boxplot() + 
  labs(x = "rs505922", y = 'Observed MBL2', 
       title = '', color = '') + 
  scale_color_brewer(palette="Dark2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
summary(lm(obs ~ rs505922, data = mbl2_tab))

### coefficient of ABO*MBL2
# standardized covarites
st_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

mbl2_ex = prot_ex[, 'MBL2']
names(mbl2_ex) = sample_id
mbl2_ex = na.omit(mbl2_ex)

mbl2_pred = fread('cis_trans_inter/06_fusion_pred/eur/MBL2_good_pred.txt')
mbl2_pred$pred = (mbl2_pred$pred - mean(mbl2_pred$pred))/
  sd(mbl2_pred$pred)

abo_pred = fread('cis_trans_inter/06_fusion_pred/eur/ABO_good_pred.txt')
abo_pred$pred = (abo_pred$pred - mean(abo_pred$pred))/
  sd(abo_pred$pred)

# use common individuals between cis and trans prot
common_id = Reduce(intersect, list(names(mbl2_ex), abo_pred$id, mbl2_pred$id))
mbl2_ex = mbl2_ex[as.character(common_id)]
abo_pred = abo_pred$pred[match(common_id, abo_pred$id)]
mbl2_pred = mbl2_pred$pred[match(common_id, mbl2_pred$id)]

st_cov_ij = st_cov[match(common_id, st_cov$V1), ]
st_cov_ij = as.matrix(st_cov_ij[,-(1:2)])

#blood_type_fit = blood_type$`23165`[match(common_id, blood_type$FID)]
## OLS of linear model
fit1 = lm(mbl2_ex ~ mbl2_pred + abo_pred + mbl2_pred * abo_pred + st_cov_ij)
## using sandwich to correct for heteroscedasticity
sandwich_se = diag(vcovHC(fit1, type = "HC"))^0.5
sandwich_t = coef(summary(fit1))[, 1]/sandwich_se
sandwich_p = pchisq(sandwich_t^2, 1, lower.tail = FALSE)
coef(summary(fit1))[c(1:3, 152), ]
sandwich_p[c(1:3, 152)]







