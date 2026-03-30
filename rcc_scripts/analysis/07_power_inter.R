library(data.table)
library(ggplot2)
library(gridExtra)
library(gaston)
setwd('/project/xuanyao/jinghui')


qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}
p_all = list()
for (i in 1:5) {
  p_all[[i]] = list()
  for (j in 1:4) {
    p_all[[i]][[j]] = readRDS(paste0('cis_trans_inter/05_power_test/snp_in_ld_gene_inter/p_r2_', i, '_h2_', j, '.rds'))
  }
}

par(mfrow = c(5,4))
cor_all = c(0,0.2,0.4,0.6,0.8)
h2_all = c(0.2, 0.3, 0.4, 0.5)
trans_h2_all = c(0.5, 0.4, 0.3, 0.2)
for (i in 1:5) {
  for (j in 1:4) {
    qqplot.pvalues(p_all[[i]][[j]][[1]], col.abline = "black", pch = 19, 
                   main = paste0('h2_a=', h2_all[j], ', h2_i=', trans_h2_all[j], 
                                 ', cor of SNPs=', cor_all[i]))
    qq_plot(p_all[[i]][[j]][[2]], col = "red", pch = 19)
    qq_plot(p_all[[i]][[j]][[3]], col = "green", pch = 19)
  }
}

legend("right", legend = c("A ~ predA+predB+predA*predB", 
                           "A ~ snp5+predB+snp5*predB", 
                           "A ~ predA+snp5+predB+snp5*predB"),
       col = c('black', "red", "green"), pch = c(19, 19, 19), cex = 0.8)


