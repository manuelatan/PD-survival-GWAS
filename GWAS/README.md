# Code for running GWAS

This repository contains scripts used to perform genome-wide association studies (GWAS) of Parkinson's Disease (PD) progression to specified clinical milestones: mortality, Hoehn and Yahr stage 3 or greater, and cognitive impairment.

# Code contents

These scripts are run within each cohort separately.


1. [make_rscripts.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/GWAS/make_rscripts.sh): This script will make a number of R scripts, one for each subset of the genome (50k SNPs). Each R script runs Cox proportional hazards model for the outcome of interest, e.g. mortality.

Run using
```
sh make_rscripts.sh
```

The R scripts were then submitted to the UCL kronos High Performance Computing job scheduler and run in parallel, e.g.
```
qsub -pe make 2 -cwd GWASscript_0.r
qsub -pe make 2 -cwd GWASscript_1.r
qsub -pe make 2 -cwd GWASscript_2.r
```

2. [interpret_results_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main//GWAS/interpret_results_script.sh): This script will combine GWAS results from all the subsets of the genome, calculate genomic inflation/lambda, and format for METAL. 

* For meta-analysis of mortality, I excluded PPMI so merged all other cohorts by chr:bp in hg19
* For meta-analysis of cognitive impairment and HY3, this included PPMI (in hg38) so merged by rsID


3. [metaanalysis_mortality_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/GWAS/metaanalysis_mortality_script.sh): Run meta-analysis in METAL. This is for the mortality GWAS, the other GWASs follow the same format just different cohorts are included.

4. [postprocessing_metaanalysis_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/GWAS/postprocessing_metaanalysis_script.sh): This script will tidy meta-analysis results, filter meta-analysis results (e.g. for heterogeneity criteria), calculate genomic inflation/lambda, format for FUMA, and format for LocusZoom plots.
