# Code for candiate gene analysis

This repository contains scripts used to extract the results for candidate genes/variants from the Parkinson's Disease (PD) progression genome-wide association studies (GWASs). 

We looked at the 90 PD risk SNPs from the most recent PD case-control GWAS meta-analysis by Nalls et al. (2019) [https://pubmed.ncbi.nlm.nih.gov/31701892/](https://pubmed.ncbi.nlm.nih.gov/31701892/).

We also looked at the PD Genetic Risk Score (GRS), a cumulative score created from the 90 risk SNPs.

Finally, we looked at other candidate variants that have been reported by other PD progression GWASs and candidate gene studies:
* SLC44A1 [(Iwaki et al., 2019)](https://pubmed.ncbi.nlm.nih.gov/31505070/)
* RIMS2 [(Liu et al., 2021)](https://pubmed.ncbi.nlm.nih.gov/33958783/)
* WWOX [(Liu et al., 2021)](https://pubmed.ncbi.nlm.nih.gov/33958783/)
* TMEM108 [(Liu et al., 2021)](https://pubmed.ncbi.nlm.nih.gov/33958783/)
* APOE e2
* MAPT H1 haplotype [(Williams-Gray et al., 2013)](https://pubmed.ncbi.nlm.nih.gov/23781007/)
* rs2242367 adjacent to LRRK2 [(Jabbari et al, 2021)](https://pubmed.ncbi.nlm.nih.gov/33341150/)


# Code contents

1. [candidate_genes_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/800e8bacb2bf44163c2ed9b8608b2aed5a23d737/candidate_genes/candidate_genes_script.R): This script looks at the results from the PD progression GWASs for the 90 PD risk loci (except for 2 which were not present due to minor allele filtering), and the other candidate variants. Creates heatmap and results tables.

2. [GRS_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/800e8bacb2bf44163c2ed9b8608b2aed5a23d737/candidate_genes/GRS_script.R): This script analyses the PD GRS vs mortality, Hoehn and Yahr stage 3+, and cognitive impairment in each cohort. Then I perform meta-analysis to combine results across cohorts.
