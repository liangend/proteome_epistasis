setwd('/project/xuanyao/jinghui')
library(gaston)
## p1 to p5 with cor between x1 and x2 ranging between 0, 0.2, 0.4, 0.6, 0.8
all_p1 = readRDS('cis_trans_inter/04_false_positive_test/snp_inter/all_p1.rds')
all_p2 = readRDS('cis_trans_inter/04_false_positive_test/snp_inter/all_p2.rds')
all_p3 = readRDS('cis_trans_inter/04_false_positive_test/snp_inter/all_p3.rds')
all_p4 = readRDS('cis_trans_inter/04_false_positive_test/snp_inter/all_p4.rds')
all_p5 = readRDS('cis_trans_inter/04_false_positive_test/snp_inter/all_p5.rds')
# Function to plot a QQ plot
qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

## each list include 3 interaction p values, each having 4 levels of h2:
# 1. y ~ x2 + x3 + x2*x3
# 2. y ~ x1 + x3 + x2*x3
# 3. y ~ x1 + x2 + x3 + x2*x3
qqplot.pvalues(all_p5[[1]][,1], col.abline = "black", pch = 19, 
               main = 'y ~ x1 + x2 + x3 + x2*x3, cor(x1, x2) = 0.8')
qq_plot(all_p5[[1]][,2], col = "red", pch = 19)
qq_plot(all_p5[[1]][,3], col = "green", pch = 19)
qq_plot(all_p5[[1]][,4], col = "purple", pch = 19)
legend("bottomright", legend = c("h2 = 0.2", "h2 = 0.4", "h2 = 0.6", "h2 = 0.8"),
       col = c('black', "red", "green", "purple"), pch = c(19, 19, 19, 19), cex = 0.8)







