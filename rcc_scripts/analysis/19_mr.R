setwd('/project/xuanyao/jinghui')
library(data.table)
library(plink2R)
library(ggplot2)
library(ivreg)
library(gaston)
library(MendelianRandomization)

load('cis_trans_inter/data/ukb_ppp_prot_expr.rdata')
colnames(prot_ex) = gene_name
ukb_cov = fread('cis_trans_inter/00_ref/gcta_cov.qcovar')

mbl2 = read_plink('cis_trans_inter/00_ref/cis_geno/MBL2',
                  impute="avg")

## trait: oral_trait/oral_issue674206, Rheumatoid_arthritis_pain674178, HDL674178, LDL674178
trait = 'LDL674178'
pheno = fread(paste0('cis_trans_inter/data/', trait, '.pheno'))

## oral ulcer
pheno = pheno[pheno$`6149` %in% c(1, -7)]
pheno = data.frame(ID = mbl2$fam$V1, 
                   pheno = pheno$`6149`[match(mbl2$fam$V1, pheno$FID)])
pheno$pheno[pheno$pheno == -7] = 0
## RA
pheno = pheno[pheno$`120001` %in% c(1, 0)]
pheno = data.frame(ID = mbl2$fam$V1, 
                   pheno = pheno$`120001`[match(mbl2$fam$V1, pheno$FID)])
## LDL
pheno = data.frame(ID = mbl2$fam$V1, 
                   pheno = pheno$`30780`[match(mbl2$fam$V1, pheno$FID)])


pheno = na.omit(pheno)
pheno$mbl2 = prot_ex[match(pheno$ID, sample_id), 'MBL2']
pheno$abo = prot_ex[match(pheno$ID, sample_id), 'ABO']
pheno = pheno[!is.na(pheno$mbl2), ]

mbl2_geno = mbl2$bed[match(pheno$ID, mbl2$fam$V1), ]

ukb_cov_sub = ukb_cov[match(pheno$ID, ukb_cov$V1), -(1:2)]
mbl2_res = residuals(lm(pheno$mbl2 ~ as.matrix(ukb_cov_sub)))

gwas_mbl2 = apply(mbl2_geno, 2, function(x){summary(lm(mbl2_res ~ x))$coef[2,4]})
mbl2_geno_sub = mbl2_geno[, which(gwas_mbl2 < 5e-8)]

mbl2_geno_cor = cor(mbl2_geno_sub)
gwas_mbl20 = gwas_mbl2[gwas_mbl2 < 5e-8]
mbl2_indep_sig_snp = c()
while (length(gwas_mbl20) > 0) {
  most_sig = which.min(gwas_mbl20)
  mbl2_indep_sig_snp = c(mbl2_indep_sig_snp, names(most_sig))
  cor_i = mbl2_geno_cor[, names(most_sig)]
  rm_snp = names(which(abs(cor_i) > 0.1))
  gwas_mbl20 = gwas_mbl20[-which(names(gwas_mbl20) %in% rm_snp)]
}

mbl2_geno_final = mbl2_geno_sub[, mbl2_indep_sig_snp]
mbl2_beta = apply(mbl2_geno_final, 2, function(x){summary(lm(mbl2_res ~ x))$coef[2,1]})
mbl2_se = apply(mbl2_geno_final, 2, function(x){summary(lm(mbl2_res ~ x))$coef[2,2]})

pheno_res = residuals(glm(pheno$pheno ~ as.matrix(ukb_cov_sub), family = gaussian))
pheno_beta = apply(mbl2_geno_final, 2, function(x){summary(lm(pheno_res ~ x))$coef[2,1]})
pheno_se = apply(mbl2_geno_final, 2, function(x){summary(lm(pheno_res ~ x))$coef[2,2]})

MRInputObject = mr_input(bx = mbl2_beta,
                         bxse = mbl2_se,
                         by = pheno_beta,
                         byse = pheno_se)
IVWObject = mr_ivw(MRInputObject)
IVWObject

abo = read_plink('cis_trans_inter/00_ref/cis_geno/ABO',
                 impute="avg")
mbl2_abo = mbl2_geno_final
mbl2_abo = cbind(mbl2_abo, rs505922 = abo$bed[match(rownames(mbl2_abo), rownames(abo$bed)), 
                                              'rs505922'])

inter_g = c()
abo_snp_center = mbl2_abo[,ncol(mbl2_abo)] - mean(mbl2_abo[,ncol(mbl2_abo)])
for (i in 1:ncol(mbl2_abo)) {
  snp_i = mbl2_abo[, i]
  snp_i = snp_i - mean(snp_i)
  inter_g = cbind(inter_g, snp_i * abo_snp_center)
}
tsls.out <- summary(ivreg(pheno_res ~ mbl2_res | inter_g), 
                    vcov = sandwich::sandwich, diagnostics = F)
tsls.out


MRSq = function(Y,A,G.mat,k=5){
  G.hat<- apply(G.mat,2,mean) ## MAF
  ## Generate Z, the new IV matrix
  K = dim(G.mat)[2]
  n = dim(G.mat)[1]
  Z<- numeric();
  for (L in  (K-k+1):K) {
    print("Creating the Z matrix")
    if (L==(K-k+1)) {
      S<- t(combn(1:K,L))
      Z.L<- matrix(nrow=n,ncol=dim(S)[1])
      for (j in 1:dim(S)[1]) {
        if (length(S[j,])==1) {
          Z.L[,j]<- G.mat[,S[j,]]-G.hat[S[j,]];
        } else {
          Z.L[,j]<- apply(sweep(G.mat[,S[j,]],2,G.hat[S[j,]]),1,prod);
        }
      }
    } else {
      S<- t(combn(1:K,L))
      Z.L<- matrix(nrow=n,ncol=dim(S)[1])
      for (j in 1:dim(S)[1]) {
        k.L  <- S[j,];
        if (K==k) {
          P.L <- 1;
        } else if ((K-k)==1) {
          P.L  <- G.mat[,k.L[1:(K-k)]]-G.hat[k.L[1:(K-k)]];
        } else {
          P.L  <- apply(sweep(G.mat[,k.L[1:(K-k)]],2,G.hat[k.L[1:(K-k)]]),1,prod)
        }
        Z.L[,j]<- P.L*(apply(G.mat[,k.L[(K-k+1):L]],1,prod)-prod(G.hat[k.L[(K-k+1):L]]));
      }
    }
    Z<- cbind(Z,Z.L);
  }
  print("The Z matrix created!, running regressions")
  print(dim(Z))
  tsls.out <- summary(ivreg(Y ~ A | Z), vcov = sandwich::sandwich, diagnostics = F)
  print(tsls.out)
  est <- c(tsls.out$coef[2,1], tsls.out$coef[2,2])
  return(est)
}







