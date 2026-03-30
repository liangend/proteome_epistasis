setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
infl_file = list.files('cis_trans_inter/04_false_positive_test/snp_sandwich_inter/ld_4/', 
                       pattern = '.txt')
n_sim = length(infl_file)

inter1_infl = list()
inter2_infl = list()
inter3_infl = list()
sand1_infl = list()
sand2_infl = list()
sand3_infl = list()

for (ld in 1:4) {
  inter1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  sand1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  sand2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  sand3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)

  for (j in 1:length(infl_file)) {
    file_j = infl_file[j]
    infl_j = fread(paste0('cis_trans_inter/04_false_positive_test/snp_sandwich_inter/ld_', 
                          ld, '/', file_j))
    inter1_infl_ld[j,] = unlist(infl_j[,2])
    inter2_infl_ld[j,] = unlist(infl_j[,3])
    inter3_infl_ld[j,] = unlist(infl_j[,4])
    
    sand1_infl_ld[j,] = unlist(infl_j[,10])
    sand2_infl_ld[j,] = unlist(infl_j[,11])
    sand3_infl_ld[j,] = unlist(infl_j[,12])

  }
  inter1_infl[[ld]] = inter1_infl_ld
  inter2_infl[[ld]] = inter2_infl_ld
  inter3_infl[[ld]] = inter3_infl_ld
  
  sand1_infl[[ld]] = sand1_infl_ld
  sand2_infl[[ld]] = sand2_infl_ld
  sand3_infl[[ld]] = sand3_infl_ld
  
}

## inflation factor of interaction p between x2 and x3
inter1_all = data.frame(infl = c(unlist(inter1_infl), unlist(inter2_infl[[1]])),
                        h2 = as.factor(rep(rep(c(0.05, 0.1, 0.2, 0.5), each = n_sim), 5)),
                        ld = as.factor(rep(c(0, 0.3, 0.6, 0.9, 1), each = n_sim * 4)))
p1 = ggplot(inter1_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = brewer.pal(4, "Dark2")) + 
  ylim(c(0.7,1.6)) +
  labs(x = "LD", y = 'inflation of SNP2 * SNP3', 
       title = 'Without SVE', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# sandwich corrected p value
sand1_all = data.frame(infl = c(unlist(sand1_infl), unlist(sand2_infl[[1]])),
                        h2 = as.factor(rep(rep(c(0.05, 0.1, 0.2, 0.5), each = n_sim), 5)),
                        ld = as.factor(rep(c(0, 0.3, 0.6, 0.9, 1), each = n_sim * 4)))
p2 = ggplot(sand1_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = brewer.pal(4, "Dark2")) + 
  ylim(c(0.7,1.6)) +
  labs(x = "LD", y = 'inflation of SNP2 * SNP3', 
       title = 'With SVE', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### put p3 figure from 04_infl_snp_inter.R together
grid.arrange(p1, p2, p3, nrow = 2)



