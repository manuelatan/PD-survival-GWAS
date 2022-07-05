# Code for clinical data summaries and power calculations

This repository contains scripts to summarise clinical data across all cohorts (e.g. the final number of individuals included in the mortality GWAS, how many individuals met the outcome or did not meet the outcome, and median time to event). 
This also includes scripts for performing power calculations using the R package survSNP.

# Code contents

1. [clinical_data_summary_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/summaries_and_power/clinical_data_summary_script.R): This script reads in clinical data across all cohorts, filters for just the individuals included in the final GWAS, and combines this data. In this script I then summarise clinical data for each outcome.

2. [power_calculations_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/summaries_and_power/power_calculations_script.R): This script performs power calculations for the survival GWAS and creates plots.
