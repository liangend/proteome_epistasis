#!/usr/bin/env bash 
#SBATCH --job-name=ukb_job
#SBATCH --time=24:00:00
#SBATCH --mem=5GB
#SBATCH --partition=tier2q
#SBATCH --output=../98_logs/job%j.out
#SBATCH --error=../98_logs/job%j.err

python3 /gpfs/data/ukb-share/extraction_scripts/extract_pheno.py 41270 \
	-n diag1 -c \
	-t ../00_ref/ukb_traits_more/
