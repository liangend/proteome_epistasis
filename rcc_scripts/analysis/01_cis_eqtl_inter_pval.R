library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
p_trans = fread('cis_trans_inter/output/trans_p_all.txt', header = T)
p_inter = fread('cis_trans_inter/output/inter_p_all.txt', header = T)

p_trans = p_trans[-which(p_trans$WASH7P == 'WASH7P'), ]
p_inter = p_inter[-which(p_inter$WASH7P == 'WASH7P'), ]
p_trans = as.data.frame(p_trans)
p_inter = as.data.frame(p_inter)

rownames(p_trans) = p_trans$V1
rownames(p_inter) = p_inter$V1
p_trans = p_trans[, -1]
p_inter = p_inter[, -1]

p_trans = apply(p_trans, 2, as.numeric)
p_inter = apply(p_inter, 2, as.numeric)


gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(mapping = aes(x = expected, ymin = clower, ymax = cupper),
                alpha = 0.1) +
    geom_point(aes(expected, observed), size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) + ylab(log10Po)
}

them_plot = theme(text = element_text(size=15, colour = "black"), 
                  axis.text.x = element_text(colour = "black", size = 15),
                  axis.text.y = element_text(colour = "black", size = 15),
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank())

gene_chr1 = colnames(p_inter)
minp_inter = apply(p_inter, 2, min)
minp_trans = apply(p_trans, 2, min)

minp_inter = sort(minp_inter)
minp_trans = sort(minp_trans)

i = 1
p1 = gg_qqplot(p_inter[, names(minp_inter)[1]]) + them_plot + labs(title = names(minp_inter)[1])
p2 = gg_qqplot(p_inter[, names(minp_inter)[2]]) + them_plot + labs(title = names(minp_inter)[2])
p3 = gg_qqplot(p_inter[, names(minp_inter)[3]]) + them_plot + labs(title = names(minp_inter)[3])
p4 = gg_qqplot(p_inter[, names(minp_inter)[4]]) + them_plot + labs(title = names(minp_inter)[4])
p5 = gg_qqplot(p_inter[, names(minp_inter)[5]]) + them_plot + labs(title = names(minp_inter)[5])
p6 = gg_qqplot(p_inter[, names(minp_inter)[6]]) + them_plot + labs(title = names(minp_inter)[6])
p7 = gg_qqplot(p_inter[, names(minp_inter)[7]]) + them_plot + labs(title = names(minp_inter)[7])
p8 = gg_qqplot(p_inter[, names(minp_inter)[8]]) + them_plot + labs(title = names(minp_inter)[8])
p9 = gg_qqplot(p_inter[, names(minp_inter)[9]]) + them_plot + labs(title = names(minp_inter)[9])
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)

p1 = gg_qqplot(p_trans[, names(minp_trans)[1]]) + them_plot + labs(title = names(minp_trans)[1])
p2 = gg_qqplot(p_trans[, names(minp_trans)[2]]) + them_plot + labs(title = names(minp_trans)[2])
p3 = gg_qqplot(p_trans[, names(minp_trans)[3]]) + them_plot + labs(title = names(minp_trans)[3])
p4 = gg_qqplot(p_trans[, names(minp_trans)[4]]) + them_plot + labs(title = names(minp_trans)[4])
p5 = gg_qqplot(p_trans[, names(minp_trans)[5]]) + them_plot + labs(title = names(minp_trans)[5])
p6 = gg_qqplot(p_trans[, names(minp_trans)[6]]) + them_plot + labs(title = names(minp_trans)[6])
p7 = gg_qqplot(p_trans[, names(minp_trans)[7]]) + them_plot + labs(title = names(minp_trans)[7])
p8 = gg_qqplot(p_trans[, names(minp_trans)[8]]) + them_plot + labs(title = names(minp_trans)[8])
p9 = gg_qqplot(p_trans[, names(minp_trans)[9]]) + them_plot + labs(title = names(minp_trans)[9])
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)




