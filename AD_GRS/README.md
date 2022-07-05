# Code for Alzheimer's Disease Genetic Risk Score analysis

This repository contains scripts used to analyse the Alzheimer's Disease (AD) Genetic Risk Score (GRS) in relation to Parkinson's Disease progression. We created AD GRSs based on the most recent AD case-control GWAS by Wightman et al. [https://pubmed.ncbi.nlm.nih.gov/34493870/](https://pubmed.ncbi.nlm.nih.gov/34493870/), excluding the APOE region. 

Initially, I downloaded the full summary statistics for the AD GWAS from [https://ctg.cncr.nl/software/summary_statistics](https://ctg.cncr.nl/software/summary_statistics), however this data contained z-scores rather than betas, and the z-scores have been calculated dependent on minor allele frequency, which was also not available in the summary statistics. Therefore I used the betas for the genome-wide significant loci published in the preprint [https://www.medrxiv.org/content/10.1101/2020.11.20.20235275v1](https://www.medrxiv.org/content/10.1101/2020.11.20.20235275v1).

# Code contents

1. [AD_PRS_script_new.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/AD_GRS/AD_PRS_script_new.R): This script formats the AD summary statistics, and analyses the AD GRS vs. mortality, Hoehn and Yahr stage 3+, and cognitive impairment in each cohort. I also conduct random-effects meta-analysis to combine results across cohorts.

2. [AD_PRS_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/AD_GRS/AD_PRS_script.sh): This script creates the AD GRS in each cohort, both with and without the APOE region.
