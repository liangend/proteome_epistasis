setwd('/project/xuanyao/jinghui')
library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

## read all additive and interaction p values 
add_p = list()
inter_p = list()
for (i in 1:22) {
  cis_chr = i
  add_i = c()
  inter_i = c()
  for (j in 1:22) {
    trans_chr = j
    if (i == j) {
      next
    }
    p_ij = readRDS(paste0('cis_trans_inter/07_prot_inter/prot_prot_inter/sandwich_p/cis_chr', 
                          cis_chr, '_trans_chr', trans_chr, '.rds'))
    add_ij = p_ij$trans_p
    add_i = rbind(add_i, add_ij)
    inter_ij = p_ij$inter_p
    inter_i = rbind(inter_i, inter_ij)
  }
  for (k in 1:ncol(inter_i)) {
    prot_p = add_i[, k]
    names(prot_p) = rownames(add_i)
    add_p[[length(add_p) + 1]] = prot_p
    
    prot_p = inter_i[, k]
    names(prot_p) = rownames(inter_i)
    inter_p[[length(inter_p) + 1]] = prot_p
  }
  names(add_p)[(length(add_p) - ncol(add_i) + 1):length(add_p)] = colnames(add_i)
  names(inter_p)[(length(inter_p) - ncol(inter_i) + 1):length(inter_p)] = colnames(inter_i)
}
all_add_p = unlist(add_p)
add_p_tab = data.frame(pval = unname(all_add_p), 
                         target = sapply(strsplit(names(all_add_p), '.', fixed = T), '[', 1),
                         regulator = sapply(strsplit(names(all_add_p), '.', fixed = T), '[', 2))

all_inter_p = unlist(inter_p)
inter_p_tab = data.frame(pval = unname(all_inter_p), 
                         target = sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 1),
                         regulator = sapply(strsplit(names(all_inter_p), '.', fixed = T), '[', 2))

## pi1 values by regulator
reg_add_pi1 = aggregate(pval ~ regulator, data = add_p_tab, 
                  function(x){1 - qvalue(x)$pi0})
reg_inter_pi1 = aggregate(pval ~ regulator, data = inter_p_tab, 
                          function(x){1 - qvalue(x)$pi0})
reg_add_pi1 = reg_add_pi1[order(reg_add_pi1$pval, decreasing = T), ]
reg_inter_pi1 = reg_inter_pi1[order(reg_inter_pi1$pval, decreasing = T), ]
p1 = ggplot(reg_add_pi1[1:10, ], aes(x = reorder(regulator, -pval), y = pval)) + 
  geom_bar(stat="identity", width = 0.8, fill = 'steelblue') +
  labs(title = 'pi1 of additive p values for each protein with other proteins', x = "", y = "pi1", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(reg_inter_pi1[1:10, ], aes(x = reorder(regulator, -pval), y = pval)) + 
  geom_bar(stat="identity", width = 0.8, fill = 'steelblue') +
  labs(title = 'pi1 of interactive p values for each protein with other proteins', x = "", y = "pi1", color = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 2)

