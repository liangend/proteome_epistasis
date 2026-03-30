#!/usr/bin/env bash 
#SBATCH --job-name=predDb
#SBATCH --time=24:00:00
#SBATCH --mem=90GB
#SBATCH --ntasks=1
#SBATCH --partition=tier2q
#SBATCH --cpus-per-task=16
#SBATCH --output=../98_logs/job%j.out
#SBATCH --error=../98_logs/job%j.err

source ~/.bashrc
conda activate pred_nextflow
module load gcc/12.1.0
module load R/4.3.1
cd /gpfs/data/xliu-lab/jinghui/software/PredictDb-nextflow
nextflow run main.nf \
	--gene_annotation '/gpfs/data/xliu-lab/jinghui/ukb_ppp/prot_annot.gtf' \
	--snp_annotation '/gpfs/data/xliu-lab/jinghui/ukb_ppp/genotype/snp_annot/chr1.vcf' \
	--genotype '/gpfs/data/xliu-lab/jinghui/ukb_ppp/genotype/geno_012/chr1.test' \
	--gene_exp '/gpfs/data/xliu-lab/jinghui/ukb_ppp/expr_cov_by_chr/prot_expr_chr1.txt' \
	--covariates '/gpfs/data/xliu-lab/jinghui/ukb_ppp/expr_cov_by_chr/cov_chr1.txt' \
	--outdir '/gpfs/data/xliu-lab/jinghui/cis_trans_inter/05_predDb' \
	-resume \
	--prefix 'chr1'
