#------------CAUSE OF DEATH SUB-ANALYSIS------------#
#Sub-analysis in PDdeath causes vs. interrupted
#In QSBB and UKB data only

#---Load packages---####
library(ggplot2)
library(survival)
library(survminer)
library(lme4)
library(data.table)
library(tidyverse)

#---Read in genotype data for top 10 SNPs from each cohort---####

QSBB_geno <- fread("QSBB.snpqc.mortality_snps.recodeA.raw")
UKB_geno <- fread("UKB.snpqc.mortality_snps.recodeA.raw")

#---Read in clinical data for UKB and QSBB---####

#Read in clinical data for UKB and QSBB - this has the cause of death data
UKBclinical <- read.table("clinical/UKB_survival_all_2020-09-17.txt", header = T)

#Select relevant columns and just keep UKB cases that have died, only incident and prevalent cases
UKBclinical <- UKBclinical %>% 
  filter(event_death == 1) %>% 
  filter(cohort == "UKB_incident" | cohort == "UKB_prevalent") %>% 
  select(FID, IID, cohort, event_death, timeToEvent_death, age_onset_imput, gender, PDdeath) %>% 
  mutate(PDdeath = ifelse(is.na(PDdeath), "interrupted", "yes"))



QSBBclinical <- read.table("clinical/QSBB_survival_all_2020-07-30.txt", header = T)

QSBBclinical <- QSBBclinical %>% 
  mutate(cohort == "QSBB") %>% 
  mutate(gender = tolower(gender)) %>% 
  select(FID, IID, cohort, event_death, timeToEvent_death, age_onset_imput, gender, PDdeath)

#---Read in Principal Components---####

QSBB_PCs <- fread("../QSBB/PCA.eigenvec")
colnames(QSBB_PCs) <- c("FID", "IID", paste("PC", 1:20, sep = ""))

UKB_PCs <- fread("../UKB/PCA.eigenvec")
colnames(UKB_PCs) <- c("FID", "IID", paste("PC", 1:20, sep = ""))

#---Merge each clinical dataset with the genotype dataset---####

QSBB_merged <- QSBBclinical %>% 
  inner_join(QSBB_geno, by = c("FID", "IID")) %>% 
  inner_join(QSBB_PCs, by = c("FID", "IID"))

UKB_merged <- UKBclinical %>% 
  inner_join(UKB_geno, by = c("FID", "IID")) %>% 
  inner_join(UKB_PCs, by = c("FID", "IID"))

#---Merge datasets---####

match <- names(UKB_merged) %in% names(QSBB_merged)
UKB_merged_selected <- UKB_merged[, match]

#Merge clinical and genotype data from each cohort
merged <- rbind(QSBB_merged, UKB_merged_selected)

#Count of each cohort
merged %>% 
  group_by(cohort) %>% 
  summarise(count = n())

merged %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

#Filter just for cases that have cause of death (i.e. remove the QSBB with NA cause of death)
merged_deathcause <- merged %>% 
  filter(!is.na(PDdeath))

merged_deathcause %>% 
  group_by(cohort, PDdeath) %>% 
  summarise(count = n())

merged_deathcause %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

#Calculate mean and median time to event
merged_deathcause %>% 
  summarise(median = median(timeToEvent_death, na.rm = TRUE))

#---Rename SNP column names with rsIDs---####

#Read in SNP positions and rsIDs (from FUMA independent SNPs)
snps <- read.table("mortality_snps_rsids")
colnames(snps) <- c("chr", "bp", "rsid", "effect_allele", "other_allele")

#Arrange SNPs by chr and bp
snps_sorted <- snps %>% 
  arrange(chr, bp) %>% 
  mutate(snp = paste(chr, ":", bp, ":", other_allele, ":", effect_allele, "_", effect_allele, sep = ""))

merged_colnames <- colnames(merged_deathcause)[13:19]

#Keep only SNPs that are in the main dataframe
rsids <- snps_sorted %>% 
  filter(snp %in% merged_colnames)

rsids <- as.vector(rsids[['rsid']])

colnames(merged_deathcause)[13:19] <- rsids

#---Make functions to plot Kaplan-Meier surves---####


### Function to plot Kaplan-Meier curve stratified by PD death cause ###
# This will select/filter for only PD deaths OR interrupted deaths
# PDdeath_val must be either "yes" (PD death) or "interrupted"
# SNP value must match the column names in the main dataset
# Will also return results from Cox model with covariates (age onset, sex, PC1-PC5)
plot_KM_stratified <- function(df, SNP, PDdeath_val){
  
  #Filter by PDdeath group
  filtered <- df %>% 
    filter(PDdeath == PDdeath_val)
  
  #Fit mortality vs. SNP
  #Note we are using the survminer::surv_fit() function rather than survfit()
  #Otherwise ggsurvplot has scoping issues
  fit <- surv_fit(Surv(time = filtered$timeToEvent_death,
                      event = filtered$event_death) ~ filtered[[SNP]],
                 data = filtered)
  
  if (PDdeath_val == "yes") {
    filename <- paste("./plots/", SNP, "_PDdeath.png", sep = "")
  } else {
    filename <- paste("./plots/", SNP, "_interrupted.png", sep = "")
  }
 
  #Make legend labels - this depends on whether there are homozygotes for the alt allele
  if (length(unique(filtered[[SNP]])) == 3) {
    legend.labs = c("0", "1", "2")
  } else {
    legend.labs = c("0", "1")
  }

  legend_title <- paste(SNP, " allele count", sep = "")
  
  #Plot Kaplan Meier curve for SNP vs. mortality in stratified data
  ggsurvplot(fit, data = filtered, pval = FALSE,
             legend.title = legend_title,
             legend.labs = legend.labs) +
    ggsave(filename, width = 6, height = 7)
  
  #Fit Cox Proportional Hazards model for mortality vs.SNP with covariates
  cox <- coxph(Surv(time = filtered$timeToEvent_death, 
                    event = filtered$event_death) ~ filtered[[SNP]]
               + age_onset_imput + gender
               + PC1 + PC2 + PC3 + PC4 + PC5, 
               data = filtered)
  
  return(summary(cox))
  
}


### Function to plot Kaplan-Meier curve faceted by PD death cause ###
# This will include the whole dataset and save a single plot faceted by PD death group
# The Kaplan-Meier curve will show SNP allele counts in a dominant model for simplicity
# I.e. SNP 0 vs. SNP 1 / 2
# PDdeath_val must be either "yes" (PD death) or "interrupted"
# SNP value must match the column names in the main dataset
# Will also return results from Cox model with interaction between SNP x death
# Covariates (age onset, sex, PC1-PC5)

plot_KM_faceted <- function(input, SNP){
  
  #Recode rs429358 SNP in a dominant model
  input <- input %>% 
    mutate(SNP_dominant = ifelse(input[[SNP]] == 0, "0",
                                 ifelse(input[[SNP]] == 1 | input[[SNP]] == 2, "1/2", NA)))

  #Fit SNP vs. mortality
  fit <- surv_fit(Surv(time = timeToEvent_death,
                       event = event_death) ~ SNP_dominant,
                  data = input)
  
  #Plot Kaplan Meier curve for SNP vs. mortality
  ggsurvplot_facet(fit = fit, data = input, 
                   pval = FALSE,
                   facet.by = "PDdeath",
                   panel.labs = list(PDdeath = c("interrupted death", "PD death")),
                   short.panel.labs = TRUE,
                   legend.title = paste(SNP, " allele count", sep = "")) +
    ggsave(paste("./plots/", SNP, "_faceted.png", sep = ""))
  
  #Fit Cox Proportional Hazards model for mortality vs. SNP with covariates
  cox <- coxph(Surv(time = timeToEvent_death,
                    event = event_death) ~ input[[SNP]]*PDdeath
                 + age_onset_imput + gender
                 + PC1 + PC2 + PC3 + PC4 + PC5,
                 data = input)

  return(summary(cox))
  
}

#---Plot Kaplan-Meier curves for APOE rs429358 SNP---####

plot_KM_stratified(merged_deathcause, "rs429358", "yes")
plot_KM_stratified(merged_deathcause, "rs429358", "interrupted")

plot_KM_faceted(input = merged_deathcause, "rs429358")

#---Plot Kaplan-Meier curves for TBXAS1 rs4726467 SNP---####

merged_deathcause %>% 
  group_by(rs4726467) %>% 
  summarise(count = n())

plot_KM_stratified(merged_deathcause, "rs4726467", "yes")
plot_KM_stratified(merged_deathcause, "rs4726467", "interrupted")

plot_KM_faceted(merged_deathcause, "rs4726467")


#---Plot Kaplan-Meier curves for TBXAS1 rs144889025 SNP---####

merged_deathcause %>% 
  group_by(rs144889025) %>% 
  summarise(count = n())

plot_KM_stratified(merged_deathcause, "rs144889025", "yes")
plot_KM_stratified(merged_deathcause, "rs144889025", "interrupted")

plot_KM_faceted(merged_deathcause, "rs144889025")

#---Plot Kaplan-Meier curves for SYT10 rs10437796 SNP---####

merged_deathcause %>%
  group_by(rs10437796) %>% 
  summarise(count = n())

plot_KM_stratified(merged_deathcause, "rs10437796", "yes")
plot_KM_stratified(merged_deathcause, "rs10437796", "interrupted")

plot_KM_faceted(merged_deathcause, "rs10437796")

