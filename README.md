# Robust cis-by-trans epistasis in the human plasma proteome highlights an ABO-centered interaction network
Code used in the study.
### Code details
01_fusion.R: FUSION model training

02_fusion_pred.R: Protein prediction based on FUSION model

03_gcta_h2_cis.R: Calculate the protein cis-h2 using GCTA

04_sim_snp_inter_ols.R: Simulation of interaction P value at the SNP level calculated by OLS

04_sim_snp_inter_sve.R: Simulation of interaction P value at the SNP level calculated by SVE

05_sim_prot_inter.R: Simulation of interaction P value at the protein level using simulated genotypes

06_sim_prot_inter_real_geno.R: Simulation of interaction P value at the protein level using real genotypes

07_prot_prot_inter_dif_chr.R: Pairwise epistasis analysis at the protein level for inter-chr pairs

08_prot_prot_inter_same_chr.R: Pairwise epistasis analysis at the protein level for intra-chr pairs

09_snp_snp_inter.R: Pairwise epistasis analysis at the SNP level

10_snp_inter_ana.R: Chromatin states of SNPs with significant epistatic effects

11_pwas_sim.R: PWAS simulation

12_ukb_epi_pwas.R: Regular and epistasis-adjusted PWAS on UKB traits

13_pred_acc.R: Protein prediction accuracy with and without epistatic effects

14_tajima.R: Enrichment of ABO-interacting genes in high Tajima'D

15_pwas_power_gain.R: Simulation of PWAS power gain with different epistatic effect sizes
