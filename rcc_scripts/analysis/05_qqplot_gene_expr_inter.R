setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gridExtra)
library(gaston)

#### cis regions of top interactions and random interactions
dgn_gene_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
dgn_gene_meta_chr1 = dgn_gene_meta[dgn_gene_meta$chr == 'chr1', ]

top_inter_gene = c('ZNF326', 'KLHL20', 'TAF13', 'DNTTIP2', 'FAM46C', 
                   'EVI5', 'SH3GLB1', 'CLIC4', 'CTNNBIP1')
inter_p = list()
for (i in 1:9) {
  gene_i = top_inter_gene[i]
  set_i = dgn_gene_meta$gene_set[which(dgn_gene_meta$gene_id == gene_i)]
  tss_i = dgn_gene_meta$start[which(dgn_gene_meta$gene_id == gene_i)]
  
  dgn_z = fread(paste0('gtex/06_qtl_z/dgn_z/geneSet', set_i, 
                       '.chr1.trans_qtl_pairs.txt.gz'))
  dgn_z_sub = dgn_z[dgn_z$phenotype_id == gene_i, ]
  dgn_z_sub$bp = as.numeric(sapply(strsplit(dgn_z_sub$variant_id, ':', fixed = T), '[', 2))
  dgn_z_cis = dgn_z_sub[which(abs(dgn_z_sub$bp - tss_i) < 1000000), ]
  
  dgn_manh = data.frame(SNP = dgn_z_cis$variant_id, P = dgn_z_cis$pval,
                        BP = dgn_z_cis$bp)
  inter_p[[i]] = ggplot() +
    geom_point(data = dgn_manh, aes(x = BP, y = -log10(P))) +
    labs(title = gene_i, y = "-log10(p)", x = "position") + 
    geom_vline(xintercept = tss_i, linetype="dashed", linewidth = 0.3) +
    theme(text = element_text(size=15, colour = "black"), 
          plot.margin = unit(c(1.5,1,0.5,0.5), "cm"),
          axis.text.x = element_text(colour = "black", size = 11),
          axis.text.y = element_text(colour = 'black', size = 11),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  print(i)
}

grid.arrange(inter_p[[1]], inter_p[[2]], inter_p[[3]], 
             inter_p[[4]], inter_p[[5]], inter_p[[6]], 
             inter_p[[7]], inter_p[[8]], inter_p[[9]], nrow = 3)

rand_inter_gene = c('WLS', 'TUFT1', 'SH2D1B', 'SOAT1', 'RASSF5', 
                   'MLLT11', 'DENND2C', 'ACADM', 'TRIM33')
rand_p = list()
for (i in 1:9) {
  gene_i = rand_inter_gene[i]
  set_i = dgn_gene_meta$gene_set[which(dgn_gene_meta$gene_id == gene_i)]
  tss_i = dgn_gene_meta$start[which(dgn_gene_meta$gene_id == gene_i)]
  
  dgn_z = fread(paste0('gtex/06_qtl_z/dgn_z/geneSet', set_i, 
                       '.chr1.trans_qtl_pairs.txt.gz'))
  dgn_z_sub = dgn_z[dgn_z$phenotype_id == gene_i, ]
  dgn_z_sub$bp = as.numeric(sapply(strsplit(dgn_z_sub$variant_id, ':', fixed = T), '[', 2))
  dgn_z_cis = dgn_z_sub[which(abs(dgn_z_sub$bp - tss_i) < 1000000), ]
  
  dgn_manh = data.frame(SNP = dgn_z_cis$variant_id, P = dgn_z_cis$pval,
                        BP = dgn_z_cis$bp)
  rand_p[[i]] = ggplot() +
    geom_point(data = dgn_manh, aes(x = BP, y = -log10(P))) +
    labs(title = gene_i, y = "-log10(p)", x = "position") + 
    geom_vline(xintercept = tss_i, linetype="dashed", linewidth = 0.3) +
    theme(text = element_text(size=15, colour = "black"), 
          plot.margin = unit(c(1.5,1,0.5,0.5), "cm"),
          axis.text.x = element_text(colour = "black", size = 11),
          axis.text.y = element_text(colour = 'black', size = 11),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  print(i)
}
grid.arrange(rand_p[[1]], rand_p[[2]], rand_p[[3]], 
             rand_p[[4]], rand_p[[5]], rand_p[[6]], 
             rand_p[[7]], rand_p[[8]], rand_p[[9]], nrow = 3)


## gene expr interaction and gene expr * cis snp interaction
qq_plot <- function(pvals, col, pch) {
  # Generate expected p-values
  expected <- -log10(ppoints(length(pvals)))
  observed <- -log10(sort(pvals))
  
  # Add to plot
  points(expected, observed, col = col, pch = pch)
}

p1_1 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_1_h2_1.rds')
p1_2 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_1_h2_2.rds')
p1_3 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_1_h2_3.rds')
p1_4 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_1_h2_4.rds')

p2_1 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_2_h2_1.rds')
p2_2 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_2_h2_2.rds')
p2_3 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_2_h2_3.rds')
p2_4 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_2_h2_4.rds')

p3_1 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_3_h2_1.rds')
p3_2 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_3_h2_2.rds')
p3_3 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_3_h2_3.rds')
p3_4 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_3_h2_4.rds')

p4_1 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_4_h2_1.rds')
p4_2 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_4_h2_2.rds')
p4_3 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_4_h2_3.rds')
p4_4 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_4_h2_4.rds')

p5_1 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_5_h2_1.rds')
p5_2 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_5_h2_2.rds')
p5_3 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_5_h2_3.rds')
p5_4 = readRDS('cis_trans_inter/04_false_positive_test/expr_inter/p_r2_5_h2_4.rds')

qqplot.pvalues(p4_1[[1]], col.abline = "black", pch = 19, 
               main = 'A ~ predB + snp1 + snp1 * predB, cor of SNPs = 0.2')
qq_plot(p4_2[[1]], col = "red", pch = 19)
qq_plot(p4_3[[1]], col = "green", pch = 19)
qq_plot(p4_4[[1]], col = "purple", pch = 19)
legend("bottomright", legend = c("h2 = 0.2", "h2 = 0.4", "h2 = 0.6", "h2 = 0.8"),
       col = c('black', "red", "green", "purple"), pch = c(19, 19, 19, 19), cex = 0.8)


