setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
infl_file = list.files('cis_trans_inter/04_false_positive_test/expr_inter/ld_5/', 
                       pattern = '.txt')
n_sim = length(infl_file)
inter1_infl = list()
inter2_infl = list()
inter3_infl = list()

x3_add1_infl = list()
x3_add2_infl = list()
x3_add3_infl = list()

for (ld in 1:5) {
  inter1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  x3_add1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  x3_add2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  x3_add3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  for (j in 1:length(infl_file)) {
    file_j = infl_file[j]
    infl_j = fread(paste0('cis_trans_inter/04_false_positive_test/expr_inter/ld_', 
                          ld, '/', file_j))
    inter1_infl_ld[j,] = unlist(infl_j[,2])
    inter2_infl_ld[j,] = unlist(infl_j[,3])
    inter3_infl_ld[j,] = unlist(infl_j[,4])
    
    x3_add1_infl_ld[j,] = unlist(infl_j[,5])
    x3_add2_infl_ld[j,] = unlist(infl_j[,6])
    x3_add3_infl_ld[j,] = unlist(infl_j[,7])
  }
  inter1_infl[[ld]] = inter1_infl_ld
  inter2_infl[[ld]] = inter2_infl_ld
  inter3_infl[[ld]] = inter3_infl_ld
  
  x3_add1_infl[[ld]] = x3_add1_infl_ld
  x3_add2_infl[[ld]] = x3_add2_infl_ld
  x3_add3_infl[[ld]] = x3_add3_infl_ld
}

## inflation factor of interaction p of predA * predB
inter1_all = data.frame(infl = unlist(inter1_infl),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 5)),
                        ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8), each = n_sim * 4)))
ggplot(inter1_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD", y = 'pval inflation of predA*predB', 
       title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## inflation factor of interaction p of predB * x5
inter2_all = data.frame(infl = unlist(inter2_infl),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 5)),
                        ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8), each = n_sim * 4)))
ggplot(inter2_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD", y = 'pval inflation of x5*predB', 
       title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## inflation factor of interaction p of predB * x5 when predA is included in the model
inter3_all = data.frame(infl = unlist(inter3_infl),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 5)),
                        ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8), each = n_sim * 4)))
ggplot(inter3_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD", y = 'pval inflation of x5*predB', 
       title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## inflation of additive effect for pred B
add_all1 = data.frame(infl = unlist(x3_add1_infl),
                        h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 5)),
                        ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8), each = n_sim * 4)))
ggplot(add_all1, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD", y = 'pval inflation of predB', 
       title = 'model: y ~ pred A + pred B + pred A * pred B', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

add_all2 = data.frame(infl = unlist(x3_add2_infl),
                      h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 5)),
                      ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8), each = n_sim * 4)))
ggplot(add_all2, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD", y = 'pval inflation of predB', 
       title = 'model: y ~ x5 + pred B + x5 * pred B', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

add_all3 = data.frame(infl = unlist(x3_add3_infl),
                      h2 = as.factor(rep(rep(c(0.2, 0.4, 0.6, 0.8), each = n_sim), 5)),
                      ld = as.factor(rep(c(0, 0.2, 0.4, 0.6, 0.8), each = n_sim * 4)))
ggplot(add_all3, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "LD", y = 'pval inflation of predB', 
       title = 'model: y ~ pred A + pred B + x5 + x5 * pred B', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



