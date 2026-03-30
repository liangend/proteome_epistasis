#!/bin/bash
#SBATCH --job-name=trans_eqtl
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10GB
#SBATCH --time=30:00:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=../98_logs/ana_%j.out
#SBATCH --error=../98_logs/ana_%j.err

module load java
java -Xmx8g -jar /project/xuanyao/jinghui/software/snpEff/snpEff.jar \
	GRCh37.75 \
	../09_snpEff/ukb_cis_qtl.vcf \
	> ../09_snpEff/ukb_cis_annot.vcf
