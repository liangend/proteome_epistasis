library(data.table)
library(ggplot2)
setwd('/project/xuanyao/jinghui')

## interaction with blood type
asso_file = list.files('cis_trans_inter/19_pheno_asso_inter_otype/')
inter_p_otype = c()
for (i in 1:length(asso_file)) {
  asso_i = readRDS(paste0('cis_trans_inter/19_pheno_asso_inter_otype/', 
                          asso_file[i]))
  inter_p_i = sapply(asso_i, function(x){min(x[5:6,6])})
  inter_p_otype = rbind(inter_p_otype, inter_p_i)
}
rownames(inter_p_otype) = sub('.rds', '', asso_file)
inter_p_otype = inter_p_otype[,c(1,2,4:6)]
inter_p_otype = inter_p_otype[!rownames(inter_p_otype) %in% c('Hypertension', 'Snoring'), ]

## CD209 example
abo_type = fread('cis_trans_inter/data/ukb_blood_type.pheno')
colnames(abo_type)[3] = 'blood_type'
cd209 = fread('cis_trans_inter/06_fusion_pred/all_eur/CD209_pred.txt')
abo = fread('cis_trans_inter/06_fusion_pred/all_eur/ABO_pred.txt')
baso = fread('cis_trans_inter/data/BASO.pheno')
colnames(baso)[3] = 'baso'
cov_all = fread('cis_trans_inter/data/all_norm_cov.txt')
cov_all = na.omit(cov_all)
baso$blood_type = abo_type$blood_type[match(baso$FID, abo_type$FID)]
baso$cd209 = cd209$pred[match(baso$FID, cd209$id)]
baso$o_type = 'OO'
baso$o_type[which(baso$blood_type %in% c('AO', 'BO'))] = 'OX'
baso$o_type[which(baso$blood_type %in% c('AA', 'AB', 'BB'))] = 'XX'

baso_sub = na.omit(baso)
baso_sub$baso = (baso_sub$baso - mean(baso_sub$baso))/sd(baso_sub$baso)

common_id = intersect(cov_all$FID, baso_sub$FID)
all_cov_sub = cov_all[match(common_id, cov_all$FID), -(1:2)]
baso_sub = baso_sub[match(common_id, baso_sub$FID), ]
baso_sub$baso_res = residuals(lm(baso_sub$baso ~ as.matrix(all_cov_sub)))
baso_sub$abo = abo$pred[match(baso_sub$FID, abo$id)]
baso_sub$inter = baso_sub$cd209 * baso_sub$abo

ggplot(baso_sub, aes(x = cd209, y = baso_res)) + 
  #geom_point() + 
  facet_wrap(~ o_type, scale = "free_x") + 
  geom_smooth(aes(color = o_type), method = "lm") + 
  labs(x = "CD209", y = 'BASO', title = '', 
       color = '') + 
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(baso_sub, aes(x = inter, y = baso_res)) + 
  #geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "CD209*ABO", y = 'BASO', title = 'p = 0.007', 
       color = '') + 
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

fit0 = lm(baso_res ~ cd209, data = baso_sub[baso_sub$o_type == 'XX', ])
summary(fit0)








