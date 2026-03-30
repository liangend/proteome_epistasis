setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(boot)
infl_file = list.files('cis_trans_inter/04_false_positive_test/expr_sand_inter/ld_4/', 
                       pattern = '.txt')
n_sim = length(infl_file)
inter1_infl = list()
inter2_infl = list()
inter3_infl = list()

inter_sand1_infl = list()
inter_sand2_infl = list()
inter_sand3_infl = list()

for (ld in 1:4) {
  inter1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  inter_sand1_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter_sand2_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  inter_sand3_infl_ld = matrix(0, nrow = n_sim, ncol = 4)
  
  for (j in 1:length(infl_file)) {
    file_j = infl_file[j]
    infl_j = fread(paste0('cis_trans_inter/04_false_positive_test/expr_sand_inter/ld_', 
                          ld, '/', file_j))
    inter1_infl_ld[j,] = unlist(infl_j[,2])
    inter2_infl_ld[j,] = unlist(infl_j[,3])
    inter3_infl_ld[j,] = unlist(infl_j[,4])
    
    inter_sand1_infl_ld[j,] = unlist(infl_j[,5])
    inter_sand2_infl_ld[j,] = unlist(infl_j[,6])
    inter_sand3_infl_ld[j,] = unlist(infl_j[,7])
  }
  inter1_infl[[ld]] = inter1_infl_ld
  inter2_infl[[ld]] = inter2_infl_ld
  inter3_infl[[ld]] = inter3_infl_ld
  
  inter_sand1_infl[[ld]] = inter_sand1_infl_ld
  inter_sand2_infl[[ld]] = inter_sand2_infl_ld
  inter_sand3_infl[[ld]] = inter_sand3_infl_ld
}

## inflation factor of interaction p of predA * predB
inter1_all = data.frame(infl = unlist(inter1_infl),
                        h2 = as.factor(rep(rep(c(0.05, 0.1, 0.2, 0.5), each = n_sim), 4)),
                        ld = as.factor(rep(c(0, 0.3, 0.6, 0.9), each = n_sim * 4)),
                        group = 'wo SVE')
p1 = ggplot(inter1_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = brewer.pal(4, "Dark2")) + 
  labs(x = "LD", y = 'pval inflation of pred A * pred B', title = 'Without SVE') + 
  ylim(c(0.8,1.4)) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

inter_sand1_all = data.frame(infl = unlist(inter_sand1_infl),
                        h2 = as.factor(rep(rep(c(0.05, 0.1, 0.2, 0.5), each = n_sim), 4)),
                        ld = as.factor(rep(c(0, 0.3, 0.6, 0.9), each = n_sim * 4)),
                        group = 'w SVE')
p2 = ggplot(inter_sand1_all, aes(x=ld, y=infl, fill=h2)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = brewer.pal(4, "Dark2")) + 
  labs(x = "LD", y = 'pval inflation of pred A * pred B', title = 'With SVE') + 
  ylim(c(0.8,1.4)) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)


### plot two groups together
inter_plt = rbind(inter1_all, inter_sand1_all)
boot_mean = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean[,3])
}
infl_boot = boot(inter_plt, boot_mean, 1000, formula = infl~ld*group)

infl_boot_plt = data.frame(infl = apply(infl_boot$t, 2, mean), 
                           quant_lower = apply(infl_boot$t, 2, function(x){quantile(x, 0.025)}), 
                           quant_higher = apply(infl_boot$t, 2, function(x){quantile(x, 0.975)}),
                           group = rep(c('With SVE', 'Without SVE'), each = 4), 
                           ld = as.factor(rep(c(0, 0.3, 0.6, 0.9), 2)))
ggplot(infl_boot_plt, 
       aes(x=ld, y=infl, 
           color=factor(group, levels = c('Without SVE', 'With SVE')))) + 
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) + 
  geom_point(position=position_dodge(0.3), size = 5) +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, 
                position=position_dodge(0.3)) + 
  scale_color_manual(values = c('darkgrey', brewer.pal(8, "Set1")[5])) + 
  labs(x = expression('LD (r'^2*')'), y = expression(lambda), 
       title = 'SVE eliminates LD-induced inflation in \n interaction tests', color = '') + 
  ylim(c(0.9,1.2)) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

