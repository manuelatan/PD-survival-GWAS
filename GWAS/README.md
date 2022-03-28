# Code for running GWAS

This repository contains scripts used to perform genome-wide association studies (GWAS) of Parkinson's Disease (PD) progression to specified clinical milestones: mortality, Hoehn and Yahr stage 3 or greater, and cognitive impairment.

# Code contents

These scripts are run within each cohort separately.


1. make_rscripts.sh: This script will make a number of R scripts, one for each subset of the genome (50k SNPs). Each R script runs Cox proportional hazards model for the outcome of interest, e.g. mortality.

Run using
```
sh make_rscripts.sh
```

The R scripts were then submitted to the kronos High Performance Computing job scheduler and run in parallel.

2. interpret_results_script.sh: 
