setwd('/gpfs/data/xliu-lab/jinghui')
input='ukb_ppp/ukb_ppp_prot_expr.rdata'
tool_dir='cis_trans_inter/99_scripts/gbat/scripts/'
library(sva)
library(data.table)

source(paste(tool_dir,'make_sva.R',sep=""))
load(input)
num_sv=20
n_na = apply(prot_ex, 2, function(x){sum(is.na(x))})
rm_col = which(n_na/nrow(prot_ex) > 0.2)
prot_ex = prot_ex[, -rm_col]
for (i in 1:ncol(prot_ex)) {
  prot_ex[,i][which(is.na(prot_ex[,i]))] = mean(prot_ex[,i], na.rm = T)
}
sv=make_sva(prot_ex,num_sv)
saveRDS(sv, 'ukb_ppp/sva.rds')
