setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gridExtra)
library(gaston)

## gene expr interaction and gene expr * cis snp interaction
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
    p_all[[i]][[j]] = readRDS(paste0('cis_trans_inter/04_false_positive_test/snp_in_ld/p_r2_', i, '_h2_', j, '.rds'))
  }
}

par(mfrow = c(3,2))
cor_all = c(0,0.2,0.4,0.6,0.8)
for (i in 1:5) {
  qqplot.pvalues(p_all[[i]][[1]][[1]], col.abline = "black", pch = 19, 
                 main = paste0('A ~ predA + predB +  SNPA*predB, cor of SNPs=', 
                               cor_all[i]))
  qq_plot(p_all[[i]][[2]][[1]], col = "red", pch = 19)
  qq_plot(p_all[[i]][[3]][[1]], col = "green", pch = 19)
  qq_plot(p_all[[i]][[4]][[1]], col = "purple", pch = 19)
}
dev.off()
par(mfrow = c(3,2))
for (i in 1:5) {
  qqplot.pvalues(p_all[[i]][[1]][[2]], col.abline = "black", pch = 19, 
                 main = paste0('A ~ predB + SNP2 +  SNP2*predB, cor of SNPs=', 
                               cor_all[i]))
  qq_plot(p_all[[i]][[2]][[2]], col = "red", pch = 19)
  qq_plot(p_all[[i]][[3]][[2]], col = "green", pch = 19)
  qq_plot(p_all[[i]][[4]][[2]], col = "purple", pch = 19)
}

dev.off()
par(mfrow = c(3,2))
for (i in 1:5) {
  qqplot.pvalues(p_all[[i]][[1]][[3]], col.abline = "black", pch = 19, 
                 main = paste0('A ~ predA + predB + SNP2 + SNP2*predB, cor of SNPs=', 
                               cor_all[i]))
  qq_plot(p_all[[i]][[2]][[3]], col = "red", pch = 19)
  qq_plot(p_all[[i]][[3]][[3]], col = "green", pch = 19)
  qq_plot(p_all[[i]][[4]][[3]], col = "purple", pch = 19)
}

legend("bottomright", legend = c("h2 = 0.2", "h2 = 0.4", "h2 = 0.6", "h2 = 0.8"),
       col = c('black', "red", "green", "purple"), pch = c(19, 19, 19, 19), cex = 0.8)


