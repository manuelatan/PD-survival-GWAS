# Code for running GWAS

This repository contains scripts used to perform genome-wide association studies (GWAS) of Parkinson's Disease (PD) progression to specified clinical milestones: mortality, Hoehn and Yahr stage 3 or greater, and cognitive impairment.

# Code contents

These scripts are run within each cohort separately.


1. [make_rscripts.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/2707bd3510eeaae40c31360ad6ba7c33eca7852a/GWAS/make_rscripts.sh): This script will make a number of R scripts, one for each subset of the genome (50k SNPs). Each R script runs Cox proportional hazards model for the outcome of interest, e.g. mortality.

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

2. [interpret_results_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/25a84193491803ba0a422de1d208d867056a5d0c/GWAS/interpret_results_script.sh): This script will combine GWAS results from all the subsets of the genome, calculate genomic inflation/lambda, format for FUMA, and format for METAL. 

* For meta-analysis of mortality, I excluded PPMI so merged all other cohorts by chr:bp in hg19
* For meta-analysis of cognitive impairment and HY3, this included PPMI (in hg38) so merged by rsID


3. [metaanalysis_mortality_script.sh](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/8ee01580ea695e8f1b0f4c44f364c008535e48da/GWAS/metaanalysis_mortality_script.sh): Run meta-analysis in METAL

4.
