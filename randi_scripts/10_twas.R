setwd('/gpfs/data/xliu-lab/jinghui')
library(data.table)
library(plink2R)
library(optparse)
args=commandArgs(trailingOnly=TRUE)
chr = as.numeric(args[1])

sumstats_file = 'cis_trans_inter/00_ref/gwas_sum_stats/LDL_30780.tsv.gz'
weights_pos = 'cis_trans_inter/00_ref/ukb_fusion.pos'
weights_dir = 'cis_trans_inter/05_fusion/output/'
ref_geno = paste0('ukb_ppp/genotype/plink_geno/chr', chr)
out = paste0('cis_trans_inter/chr', chr, '.dat')
force_model = 'enet'
min_r2pred = 0.7

allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)
  
  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip
  
  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  
  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  
  return(snp)
}

# Load in summary stats
sumstat = fread(sumstats_file)
colnames(sumstat)[3:4] = c('A1', 'A2')
sumstat$Z = sumstat$BETA / sumstat$SE

# Load in list of weights
# TODO : TEST FOR NO HEADER HERE
wgtlist = fread(weights_pos)
wgtlist = wgtlist[wgtlist$CHR == chr, ]

N = nrow(wgtlist)
out.tbl = data.frame( "PANEL" = rep(NA,N) , "FILE" = character(N) , "ID" = character(N) , 
                      "CHR" = numeric(N) , "P0" = character(N) , "P1" = character(N) ,
                      "HSQ" = numeric(N) , "BEST.GWAS.ID" = character(N) , 
                      "BEST.GWAS.Z" = numeric(N) , "EQTL.ID" = character(N) , 
                      "EQTL.R2" = numeric(N) , "EQTL.Z" = numeric(N) , 
                      "EQTL.GWAS.Z" = numeric(N) , "NSNP" = numeric(N) , 
                      "NWGT" = numeric(N) , "MODEL" = character(N) , 
                      "MODELCV.R2" = character(N) , "MODELCV.PV" = character(N) , 
                      "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , 
                      stringsAsFactors=FALSE )

# Load in reference data
genos = read_plink(paste(ref_geno,'_QC',sep=''),impute="avg")
#genos = read_plink(ref_geno, impute="avg")

# Match summary data to input, record NA where summary data is missing
m = match(genos$bim[,2] , sumstat$SNP)
sum.missing = is.na(m)
sumstat = sumstat[m,]
sumstat$SNP = genos$bim[,2]
sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]

# QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )

# Flip Z-scores for mismatching alleles
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
  genos$bim = genos$bim[qc$keep,]
  genos$bed = genos$bed[,qc$keep]
  sumstat = sumstat[qc$keep,]
}

# TODO: WARNING if too many NAs in summary stats

FAIL.ctr = 0

## For each wgt file:
for ( w in 1:nrow(wgtlist) ) {
  #cat( unlist(wgtlist[w,]) , '\n' )
  # Load weights
  wgt.file = paste(wgtlist$WGT[w],sep='')
  load(wgt.file)
  # Remove NAs (these should not be here)
  wgt.matrix[is.na(wgt.matrix)] = 0
  
  # Match up the SNPs and weights
  m = match( snps[,2] , genos$bim[,2] )
  m.keep = !is.na(m)
  snps = snps[m.keep,]
  wgt.matrix = wgt.matrix[m.keep,,drop=F]
  cur.genos = scale(genos$bed[,m[m.keep]])
  cur.bim = genos$bim[m[m.keep],]
  # Flip WEIGHTS for mismatching alleles
  qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
  wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
  
  cur.FAIL = FALSE
  
  # Match up the SNPs and the summary stats
  m = match(cur.bim[,2] , sumstat$SNP)
  cur.Z = sumstat$Z[m]
  
  # which rows have rsq
  row.rsq = grep( "rsq" , rownames(cv.performance) )
  # which rows have p-values
  row.pval = grep( "pval" , rownames(cv.performance) )	
  
  # Identify the best model
  if ( !is.na(force_model) ) {
    mod.best = which( colnames(wgt.matrix) == force_model )
    if ( length(mod.best) == 0 ) {
      cat( "WARNING : --force_model" , mod.best ,"does not exist for", unlist(wgtlist[w,]) , "\n")
      cur.FAIL = TRUE
    }	
  } else {
    # get the most significant model
    mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
  }
  if ( length(mod.best) == 0 ) {
    cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
    FAIL.ctr = FAIL.ctr + 1
    next
  }
  
  if ( sum(wgt.matrix[, mod.best] != 0) == 0 ) {
    cat( "WARNING : " , unlist(wgtlist[w,]) , names(cv.performance)[ mod.best ] , "had", length(cur.Z) , "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
    cur.FAIL = TRUE
  }
  
  # if this is a top1 model, clear out all the other weights
  if ( substr( (colnames(cv.performance))[ mod.best ],1,4) == "top1" ) wgt.matrix[ -which.max(wgt.matrix[,mod.best]^2)  , mod.best] = 0
  
  # Compute LD matrix
  if ( length(cur.Z) == 0 ) {
    cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping SNPs\n")
    cur.FAIL = TRUE
    out.tbl$NSNP[w] = NA
  } else if ( !cur.FAIL ) {
    cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)	
    out.tbl$NSNP[w] = nrow(cur.LD)
    cur.miss = is.na(cur.Z)
    # Impute missing Z-scores
    if ( sum(cur.miss) != 0 ) {
      if ( sum(!cur.miss) == 0 ) {
        cat( "WARNING : " , unlist(wgtlist[w,]) , "had no overlapping GWAS Z-scores, skipping this gene\n")
        cur.FAIL = TRUE
      } else if ( mean(cur.miss) > opt$max_impute ) {
        cat( "WARNING : " , unlist(wgtlist[w,]) , "had" , sum(cur.miss) , "/" , length(cur.miss) , "non-overlapping GWAS Z-scores, skipping this gene.\n")
        cur.FAIL = TRUE
      } else {
        cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
        cur.impz = cur.wgt %*% cur.Z[!cur.miss]
        cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
        cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
        
        all.r2pred = rep(1,length(cur.Z))
        all.r2pred[ cur.miss ] = cur.r2pred
        if ( sum(is.na(all.r2pred)) != 0 ) {
          cat( "WARNING : " , unlist(wgtlist[w,]) , "had missing GWAS Z-scores that could not be imputed, skipping this gene.\n" )
          cur.FAIL = TRUE
        } else if ( mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) < min_r2pred ) {
          cat( "WARNING : " , unlist(wgtlist[w,]) , "had mean GWAS Z-score imputation r2 of" , mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) , "at expression weight SNPs, skipping this gene.\n")
          cur.FAIL = TRUE
        }
      }
    }
    
    if ( !cur.FAIL ) {
      # Compute TWAS Z-score
      cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
      cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]
      
      if ( cur.twasr2pred > 0 ) {
        cur.twas = cur.twasz / sqrt(cur.twasr2pred)
      } else {
        cur.FAIL=T
        cat( "WARNING : " , unlist(wgtlist[w,]) , " had zero predictive accuracy, try a different model.\n")
      }
    }
  }
  
  # populate the output
  out.tbl$FILE[w] = wgt.file
  out.tbl$CHR[w] = wgtlist$CHR[w]
  out.tbl$P0[w] = wgtlist$P0[w]
  out.tbl$P1[w] = wgtlist$P1[w]
  out.tbl$ID[w] = wgtlist$ID[w]
  if ( exists("hsq") ) {
    out.tbl$HSQ[w] = hsq[1]
  }
  out.tbl$MODEL[w] = colnames( cv.performance )[ mod.best ]
  out.tbl$MODELCV.R2[w] = paste(format(cv.performance[row.rsq,mod.best],digits=2,trim=T),collapse=',')
  out.tbl$MODELCV.PV[w] = paste(format(cv.performance[row.pval,mod.best],digits=2,trim=T),collapse=',')
  
  eqtlmod = colnames(wgt.matrix) == "top1"
  topeqtl = which.max( wgt.matrix[,eqtlmod]^2 )
  
  if ( cur.FAIL || sum(eqtlmod) == 0 || length(topeqtl) == 0 || is.na(topeqtl) ) {
    out.tbl$EQTL.ID[w] = NA
    out.tbl$EQTL.R2[w] = NA
    out.tbl$EQTL.Z[w] =  NA
    out.tbl$EQTL.GWAS.Z[w] = NA
  } else {
    out.tbl$EQTL.ID[w] = rownames(wgt.matrix)[topeqtl]
    out.tbl$EQTL.R2[w] = cv.performance[1,eqtlmod]
    out.tbl$EQTL.Z[w] = wgt.matrix[ topeqtl , eqtlmod ]
    out.tbl$EQTL.GWAS.Z[w] = cur.Z[ topeqtl ]
    
  }
  
  topgwas = which.max( cur.Z^2 )
  if ( !cur.FAIL && length(topgwas) != 0 && !is.na(topgwas) ) {
    out.tbl$BEST.GWAS.ID[w] = snps[ topgwas , 2 ]
    out.tbl$BEST.GWAS.Z[w] = cur.Z[ topgwas ]
  } else {
    out.tbl$BEST.GWAS.ID[w] = NA
    out.tbl$BEST.GWAS.Z[w] = NA
  }
  
  if ( !cur.FAIL ) {
    out.tbl$NWGT[w] = sum( wgt.matrix[,mod.best] != 0 )
    out.tbl$TWAS.Z[w] = cur.twas
    out.tbl$TWAS.P[w] = 2*(pnorm( abs(out.tbl$TWAS.Z[w]) , lower.tail=F))
  } else {
    out.tbl$TWAS.Z[w] = NA
    out.tbl$TWAS.P[w] = NA
  }
  
  if ( cur.FAIL ) FAIL.ctr = FAIL.ctr + 1
}

cat("Analysis completed.\n")
cat("NOTE:",FAIL.ctr,"/",nrow(wgtlist),"genes were skipped\n")
if ( FAIL.ctr / nrow(wgtlist) > 0.1 ) {
  cat("If a large number of genes were skipped, verify that your GWAS Z-scores, expression weights, and LDREF data use the same SNPs (or nearly)\n")
  cat("Or consider pre-imputing your summary statistics to the LDREF markers using summary-imputation software such as [https://github.com/bogdanlab/fizi]\n")
}

# WRITE MHC TO SEPARATE FILE
mhc = as.numeric(out.tbl$CHR) == 6 & as.numeric(out.tbl$P0) > 26e6 & as.numeric(out.tbl$P1) < 34e6

out.tbl$P0 = apply( as.matrix(out.tbl$P0) , 1 , toString )
out.tbl$P1 = apply( as.matrix(out.tbl$P1) , 1 , toString )

if ( sum( mhc ) > 0 ) {
  cat("Results in the MHC are written to",paste(out,".MHC",sep=''),", evaluate with caution due to complex LD structure\n")
  write.table( format( out.tbl[mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=paste(out,".MHC",sep='') )
}
write.table( format( out.tbl[!mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=out )


