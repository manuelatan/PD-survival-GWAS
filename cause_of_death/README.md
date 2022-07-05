# Code for cause of death analysis

This repository contains scripts for performing cause of death analysis for the top 3 Parkinson's disease (PD) mortality genome-wide association study (GWAS) loci. 

In the subset of cases who had died and had cause of death available, we classified cause of death as either:
1. Related to PD and end of life, or
2. 'Interrupted' and likey unrelated to PD

In all these cases, we looked at the interaction between death cause and each of the top SNPs. We also conducted stratified analyses, where we separately analysed patients with PD-related cause of death, and patients with non-PD cause of death.

# Code contents

1. [extract_mortality_snps_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/cause_of_death/extract_mortality_snps_script.sh): In plink, this script extracts the top 10 loci from the PD mortality GWAS from the datasets of interest and recodes them into .raw format which is then analysed in R.

2. [cause_of_death_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/cause_of_death/cause_of_death_script.R): This script performs the competing risk analysis for cause of death, using the Fine-Gray method, for each of the top mortality SNPs. It also plots the cumulative incidence function plots for the top 4 mortality SNPs.
