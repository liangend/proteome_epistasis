library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
p_trans = fread('cis_trans_inter/03_prot_cis_cvgp_inter/trans_p_all.txt', header = T)
p_inter = fread('cis_trans_inter/03_prot_cis_cvgp_inter/inter_p_all.txt', header = T)

p_trans = as.data.frame(p_trans)
p_inter = as.data.frame(p_inter)

rownames(p_trans) = p_trans$V1
rownames(p_inter) = p_inter$V1
p_trans = p_trans[, -1]
p_inter = p_inter[, -1]


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

minp_inter = apply(p_inter, 2, min)
minp_trans = apply(p_trans, 2, min)

minp_inter = sort(minp_inter)
minp_trans = sort(minp_trans)

## most significant examples
inter_qq = list()
inter_hist = list()
for (i in 1:9) {
  inter_qq[[i]] = gg_qqplot(p_inter[, names(minp_inter)[i]]) + them_plot + labs(title = names(minp_inter)[i])
  inter_hist[[i]] = ggplot(data.frame(p = p_inter[, names(minp_inter)[i]]), aes(x=p)) +
    geom_histogram(binwidth=0.01) + them_plot +
    labs(title = names(minp_inter)[i]) 
}
grid.arrange(inter_qq[[1]], inter_qq[[2]], inter_qq[[3]], 
             inter_qq[[4]], inter_qq[[5]], inter_qq[[6]], 
             inter_qq[[7]], inter_qq[[8]], inter_qq[[9]], nrow = 3)
grid.arrange(inter_hist[[1]], inter_hist[[2]], inter_hist[[3]], 
             inter_hist[[4]], inter_hist[[5]], inter_hist[[6]], 
             inter_hist[[7]], inter_hist[[8]], inter_hist[[9]], nrow = 3)


trans_qq = list()
trans_hist = list()
for (i in 1:9) {
  trans_qq[[i]] = gg_qqplot(p_trans[, names(minp_trans)[i]]) + them_plot + labs(title = names(minp_trans)[i])
  trans_hist[[i]] = ggplot(data.frame(p = p_trans[, names(minp_trans)[i]]), aes(x=p)) +
    geom_histogram(binwidth=0.01) + them_plot +
    labs(title = names(minp_trans)[i]) 
}
grid.arrange(trans_qq[[1]], trans_qq[[2]], trans_qq[[3]], 
             trans_qq[[4]], trans_qq[[5]], trans_qq[[6]], 
             trans_qq[[7]], trans_qq[[8]], trans_qq[[9]], nrow = 3)
grid.arrange(trans_hist[[1]], trans_hist[[2]], trans_hist[[3]], 
             trans_hist[[4]], trans_hist[[5]], trans_hist[[6]], 
             trans_hist[[7]], trans_hist[[8]], trans_hist[[9]], nrow = 3)


## random examples
inter_qq_rand = list()
inter_hist_rand = list()
set.seed(101)
rand_gene = sample(colnames(p_inter), 9)
for (i in 1:9) {
  inter_qq_rand[[i]] = gg_qqplot(p_inter[, rand_gene[i]]) + them_plot + labs(title = rand_gene[i])
  inter_hist_rand[[i]] = ggplot(data.frame(p = p_inter[, rand_gene[i]]), aes(x=p)) +
    geom_histogram(binwidth=0.01) + them_plot +
    labs(title = rand_gene[i]) 
}
grid.arrange(inter_qq_rand[[1]], inter_qq_rand[[2]], inter_qq_rand[[3]], 
             inter_qq_rand[[4]], inter_qq_rand[[5]], inter_qq_rand[[6]], 
             inter_qq_rand[[7]], inter_qq_rand[[8]], inter_qq_rand[[9]], nrow = 3)
grid.arrange(inter_hist_rand[[1]], inter_hist_rand[[2]], inter_hist_rand[[3]], 
             inter_hist_rand[[4]], inter_hist_rand[[5]], inter_hist_rand[[6]], 
             inter_hist_rand[[7]], inter_hist_rand[[8]], inter_hist_rand[[9]], nrow = 3)


trans_qq_rand = list()
trans_hist_rand = list()
for (i in 1:9) {
  trans_qq_rand[[i]] = gg_qqplot(p_trans[, rand_gene[i]]) + them_plot + labs(title = rand_gene[i])
  trans_hist_rand[[i]] = ggplot(data.frame(p = p_trans[, rand_gene[i]]), aes(x=p)) +
    geom_histogram(binwidth=0.01) + them_plot +
    labs(title = rand_gene[i]) 
}
grid.arrange(trans_qq_rand[[1]], trans_qq_rand[[2]], trans_qq_rand[[3]], 
             trans_qq_rand[[4]], trans_qq_rand[[5]], trans_qq_rand[[6]], 
             trans_qq_rand[[7]], trans_qq_rand[[8]], trans_qq_rand[[9]], nrow = 3)
grid.arrange(trans_hist_rand[[1]], trans_hist_rand[[2]], trans_hist_rand[[3]], 
             trans_hist_rand[[4]], trans_hist_rand[[5]], trans_hist_rand[[6]], 
             trans_hist_rand[[7]], trans_hist_rand[[8]], trans_hist_rand[[9]], nrow = 3)


