setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
infl_file = list.files('cis_trans_inter/04_false_positive_test/snp_norm_inter/ld_5/', 
                       pattern = '.txt')
n_sim = length(infl_file)
inter1_infl = list()
inter2_infl = list()
inter3_infl = list()

x3_add1_infl = list()
x3_add2_infl = list()
x3_add3_infl = list()

x2_add1_infl = list()
x2_add3_infl = list()

for (ld in 1:5) {
  inter1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  x3_add1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  x3_add2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  x3_add3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  x2_add1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  x2_add3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  for (j in 1:length(infl_file)) {
    file_j = infl_file[j]
    infl_j = fread(paste0('cis_trans_inter/04_false_positive_test/snp_norm_inter/ld_', 
                          ld, '/', file_j))
    inter1_infl_ld[j,] = unlist(infl_j[,2])
    inter2_infl_ld[j,] = unlist(infl_j[,3])
    inter3_infl_ld[j,] = unlist(infl_j[,4])
    
    x3_add1_infl_ld[j,] = unlist(infl_j[,5])
    x3_add2_infl_ld[j,] = unlist(infl_j[,6])
    x3_add3_infl_ld[j,] = unlist(infl_j[,7])
    
    x2_add1_infl_ld[j,] = unlist(infl_j[,8])
    x2_add3_infl_ld[j,] = unlist(infl_j[,9])
  }
  inter1_infl[[ld]] = inter1_infl_ld
  inter2_infl[[ld]] = inter2_infl_ld
  inter3_infl[[ld]] = inter3_infl_ld
  
  x3_add1_infl[[ld]] = x3_add1_infl_ld
  x3_add2_infl[[ld]] = x3_add2_infl_ld
  x3_add3_infl[[ld]] = x3_add3_infl_ld
  
  x2_add1_infl[[ld]] = x2_add1_infl_ld
  x2_add3_infl[[ld]] = x2_add3_infl_ld
}

## inflation factor of interaction p between x2 and x3
inter1_all = data.frame(infl = c(unlist(inter1_infl), unlist(inter2_infl[[1]])),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 6)),
                        ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8, 1), each = n_sim * 4)))
ggplot(inter1_all[inter1_all$ld != 0.4, ], aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = brewer.pal(4, "Dark2")) + 
  labs(x = "LD", y = 'inflation of SNP2 * SNP3', 
       title = ' Truth: y = x1 + e \n Model: y ~ x2 + x3 + x2*x3', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## inflation factor of interaction p between x2 and x3 when x1 is in the model
inter2_all = data.frame(infl = c(unlist(inter3_infl), unlist(inter2_infl[[1]])),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 6)),
                        # ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8, 1), each = n_sim * 4)))
                        ld = as.factor(rep(c(0, 0.3, 0.4, 0.6, 0.9, 1), each = n_sim * 4)))
p3 = ggplot(inter2_all[inter2_all$ld != 0.4, ], aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = brewer.pal(4, "Dark2")) + 
  ylim(c(0.7,1.6)) +
  labs(x = "LD", y = 'inflation of SNP2 * SNP3', 
       title = 'Without SVE, SNP1 included', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## inflation factor of additive p for x3
add_all = data.frame(infl = c(unlist(x3_add1_infl), unlist(x3_add2_infl[[1]])),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 6)),
                        ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8, 1), each = n_sim * 4)))
ggplot(add_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD (cor between x1 and x2)", y = 'pval inflation of x3', 
       title = ' Truth: y = x1 + e \n Model: y ~ x2 + x3 + x2*x3', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


