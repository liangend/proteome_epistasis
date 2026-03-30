setwd('/project/xuanyao/jinghui')
library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
library(gaston)
library(RColorBrewer)
library(simtrait)

qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

infl_file = list.files('cis_trans_inter/05_power_test/gene_gene_inter/obs_sand', pattern = 'txt')
n_sim = length(infl_file)
inter1_infl = matrix(0, nrow = n_sim, ncol = 4)
inter2_infl = matrix(0, nrow = n_sim, ncol = 4)
inter3_infl = matrix(0, nrow = n_sim, ncol = 4)

inter1_sand_infl = matrix(0, nrow = n_sim, ncol = 4)
inter2_sand_infl = matrix(0, nrow = n_sim, ncol = 4)
inter3_sand_infl = matrix(0, nrow = n_sim, ncol = 4)
  
for (j in 1:length(infl_file)) {
  file_j = infl_file[j]
  infl_j = fread(paste0('cis_trans_inter/05_power_test/gene_gene_inter/obs_sand/', 
                        file_j))
  inter1_infl[j,] = unlist(infl_j[,2])
  inter2_infl[j,] = unlist(infl_j[,3])
  inter3_infl[j,] = unlist(infl_j[,4])
  
  inter1_sand_infl[j,] = unlist(infl_j[,5])
  inter2_sand_infl[j,] = unlist(infl_j[,6])
  inter3_sand_infl[j,] = unlist(infl_j[,7])
}
inter1_infl = as.data.frame(inter1_infl)
inter1_sand_infl = as.data.frame(inter1_sand_infl)
inter1_all = data.frame(infl = c(unlist(inter1_infl), unlist(inter1_sand_infl)),
                        h2 = as.factor(rep(rep(c(0.18, 0.13, 0.08, 0.03), each = n_sim), 2)),
                        method = rep(c('OLS', 'sandwich'), each = 100 * 4))

ggplot(inter1_all, aes(x=h2, y=infl, fill = method)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = "inter h2", y = 'pval inflation factor of predA*predB', 
       title = '', color = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

pval_i = readRDS('cis_trans_inter/05_power_test/gene_gene_inter/obs_sand/p_2.rds')
par(mfrow = c(2,2))
qqplot.pvalues(pval_i[[1]][,1], col.abline = "black", pch = 19, 
               main = 'add h2 = 0.05, inter h2 = 0.18')
qq_plot(pval_i[[4]][,1], col = "red", pch = 19)
legend('topleft', legend = c('OLS', 'sandwich'), col = c('black', 'red'),
       pch = 19, cex = 0.7)

qqplot.pvalues(pval_i[[1]][,2], col.abline = "black", pch = 19, 
               main = 'add h2 = 0.1, inter h2 = 0.13')
qq_plot(pval_i[[4]][,2], col = "red", pch = 19)

qqplot.pvalues(pval_i[[1]][,3], col.abline = "black", pch = 19, 
               main = 'add h2 = 0.15, inter h2 = 0.08')
qq_plot(pval_i[[4]][,3], col = "red", pch = 19)

qqplot.pvalues(pval_i[[1]][,4], col.abline = "black", pch = 19, 
               main = 'add h2 = 0.2, inter h2 = 0.03')
qq_plot(pval_i[[4]][,4], col = "red", pch = 19)



