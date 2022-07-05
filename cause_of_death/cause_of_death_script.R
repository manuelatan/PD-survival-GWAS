#------------CAUSE OF DEATH SUB-ANALYSIS------------#
#Sub-analysis in PDdeath causes vs. interrupted
#In QSBB and UKB data only

#---Load packages---####
library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(tidyverse)
library(cmprsk)

#---Load cmprsk function---####

#Downloaded from http://www.stat.unipg.it/luca/R
#From https://www.nature.com/articles/1705727
#I edited the CumIncidence function to save the plots in the ./plots/ directory
source("./cmprsk/CumIncidence_MT.r")


#---Read in genotype data for top 10 SNPs from each cohort---####

QSBB_geno <- fread("QSBB.snpqc.mortality_snps.recodeA.raw")
UKB_geno <- fread("UKB.snpqc.mortality_snps.recodeA.raw")

#---Read in clinical data for UKB and QSBB---####

#Read in clinical data for UKB and QSBB - this has the cause of death data
UKBclinical <- read.table("clinical/UKB_survival_all_2020-09-17.txt", header = T)

#PDdeath is coded as PDdeath == yes, or NA for non-PD
UKBclinical %>% 
  group_by(event_death, PDdeath) %>% 
  summarise(count = n())

#Select relevant columns and only incident and prevalent cases
UKBclinical <- UKBclinical %>% 
  filter(cohort == "UKB_incident" | cohort == "UKB_prevalent") %>% 
  select(FID, IID, cohort, event_death, timeToEvent_death, age_onset_imput, gender, PDdeath) %>% 
  mutate(PDdeath = ifelse(event_death == 0 & is.na(PDdeath), 0, #0 indicates censored, still alive
                          ifelse(event_death == 1 & is.na(PDdeath), 2, #Competing cause of death (non-PD) coded as 2
                                 ifelse(event_death == 1 & PDdeath == "yes", 1, NA)))) #Main cause of death (PD) coded as 1

#Check coding
UKBclinical %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

#Read in QSBB data
QSBBclinical <- read.table("clinical/QSBB_survival_all_2020-07-30.txt", header = T)

QSBBclinical %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

QSBBclinical <- QSBBclinical %>% 
  mutate(cohort == "QSBB") %>% 
  mutate(gender = tolower(gender),
         PDdeath = ifelse(PDdeath == "yes", 1,
                          ifelse(PDdeath == "interrupted", 2, NA))) %>% 
  select(FID, IID, cohort, event_death, timeToEvent_death, age_onset_imput, gender, PDdeath)

QSBBclinical %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

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

#Match column names
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

#---Rename SNP column names with rsIDs---####

#Read in SNP positions and rsIDs (from FUMA independent SNPs)
snps <- read.table("mortality_snps_rsids")
colnames(snps) <- c("chr", "bp", "rsid", "effect_allele", "other_allele")

#Arrange SNPs by chr and bp
snps_sorted <- snps %>% 
  arrange(chr, bp) %>% 
  mutate(snp = paste(chr, ":", bp, ":", other_allele, ":", effect_allele, "_", effect_allele, sep = ""))

merged_colnames <- colnames(merged)[13:19]

#Keep only SNPs that are in the main dataframe
rsids <- snps_sorted %>% 
  filter(snp %in% merged_colnames)

rsids <- as.vector(rsids[['rsid']])

colnames(merged)[13:19] <- rsids

#---Remove cases missing time to event or PD death data---####

merged %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

#Remove cases missing time to event data
#Otherwise the competing risk analysis cannot produce time point estimates and SEs
merged_nomiss <- merged %>% 
  filter(!is.na(timeToEvent_death), !is.na(PDdeath))

#---Summaries---####

merged_nomiss %>% 
  group_by(PDdeath) %>% 
  summarise(count = n())

merged_nomiss %>% 
  group_by(PDdeath) %>% 
  summarise(mean_time = mean(timeToEvent_death),
            median_time = median(timeToEvent_death))

#Check mean and median survival time in died vs. surviving cases (not separating PD death vs. non-PD death)
merged_nomiss %>% 
  mutate(died = ifelse(PDdeath == 1 | PDdeath == 2, 1,
                       ifelse(PDdeath == 0, 0, NA))) %>% 
  group_by(died) %>% 
  summarise(mean_time = mean(timeToEvent_death),
            median_time = median(timeToEvent_death))


#---Recode gender as numeric---####

#Recode gender as numeric - needed for crr function
merged_nomiss <- merged_nomiss %>% 
  mutate(SEX = ifelse(gender == "male", 1,
                      ifelse(gender == "female", 2, NA)))


#---Cumulative incidence curves for APOE rs429358---####


#Summary of SNP allele carriers
merged_nomiss %>% 
  group_by(rs429358) %>% 
  summarise(count = n())

#Code SNP in dominant model
merged_nomiss <- merged_nomiss %>% 
  mutate(rs429358_dominant = ifelse(rs429358 == 0, 0,
                                    ifelse(rs429358 == 1 | rs429358 == 2, 1, NA)))

#Check if any missing APOE SNP data
merged_nomiss %>% 
  filter(is.na(rs429358_dominant)) %>% 
  summarise(count = n())


#Make SNP into a factor
rs429358_dominant_fact <- factor(merged_nomiss$rs429358_dominant, levels=c(0,1), labels= c("rs429358 non-carriers", "rs429358 carriers"))

#Run analysis
fit <- CumIncidence_MT(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, rs429358_dominant_fact, cencode = 0, xlab="Years from onset")

#---Subdistribution hazard models for APOE rs429358---####

#Make covariate matrix
cov.mat <- cbind(merged_nomiss$rs429358_dominant, merged_nomiss$age_onset_imput, merged_nomiss$SEX,
                 merged_nomiss$PC1, merged_nomiss$PC2, merged_nomiss$PC3, merged_nomiss$PC4, merged_nomiss$PC5)

# Subdistribution hazard model for PD death (failcode = 1) 
crr.1 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 1, cencode = 0)
summary(crr.1)
#The coefficients from the summary(crr) are in the order that we specified in the cov.mat

# Subdistribution hazard model for nonPD death (failcode = 2)
crr.2 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 2, cencode = 0)
summary(crr.2)



#---Cumulative incidence curves for TBXAS1 rs4726467---####

#Summary of SNP allele carriers
merged_nomiss %>% 
  group_by(rs4726467) %>% 
  summarise(count = n())

#Code SNP in dominant model
merged_nomiss <- merged_nomiss %>% 
  mutate(rs4726467_dominant = ifelse(rs4726467 == 0, 0,
                                     ifelse(rs4726467 == 1 | rs4726467 == 2, 1, NA)))

#Check if any missing SNP data
merged_nomiss %>% 
  filter(is.na(rs4726467_dominant)) %>% 
  summarise(count = n())


#Make SNP into a factor
rs4726467_dominant_fact <- factor(merged_nomiss$rs4726467_dominant, levels=c(0,1), labels= c("rs4726467 non-carriers", "rs4726467 carriers"))

#Run analysis
fit <- CumIncidence_MT(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, rs4726467_dominant_fact, cencode = 0, xlab="Years from onset")

#---Subdistribution hazard models for TBXAS1 SNP rs4726467---####

#Make covariate matrix
cov.mat <- cbind(merged_nomiss$rs4726467_dominant, merged_nomiss$age_onset_imput, merged_nomiss$SEX,
                 merged_nomiss$PC1, merged_nomiss$PC2, merged_nomiss$PC3, merged_nomiss$PC4, merged_nomiss$PC5)

# Subdistribution hazard model for PD death (failcode = 1) 
crr.1 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 1, cencode = 0)
summary(crr.1)
#The coefficients from the summary(crr) are in the order that we specified in the cov.mat

# Subdistribution hazard model for nonPD death (failcode = 2)
crr.2 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 2, cencode = 0)
summary(crr.2)





#---Cumulative incidence curves for TBXAS1 rs144889025---####

#Summary of SNP allele carriers
merged_nomiss %>% 
  group_by(rs144889025) %>% 
  summarise(count = n())

#Code SNP in dominant model
merged_nomiss <- merged_nomiss %>% 
  mutate(rs144889025_dominant = ifelse(rs144889025 == 0, 0,
                                       ifelse(rs144889025 == 1 | rs144889025 == 2, 1, NA)))

merged_nomiss %>% 
  group_by(rs144889025_dominant) %>% 
  summarise(count = n())

#Check if any missing SNP data
merged_nomiss %>% 
  filter(is.na(rs144889025_dominant)) %>% 
  summarise(count = n())


#Make SNP into a factor
rs144889025_dominant_fact <- factor(merged_nomiss$rs144889025_dominant, levels=c(0,1), labels= c("rs144889025 non-carriers", "rs144889025 carriers"))

#Run analysis
fit <- CumIncidence_MT(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, rs144889025_dominant_fact, cencode = 0, xlab="Years from onset")

#---Subdistribution hazard models for TBXAS1 SNP rs144889025---####

#Make covariate matrix
cov.mat <- cbind(merged_nomiss$rs144889025_dominant, merged_nomiss$age_onset_imput, merged_nomiss$SEX,
                 merged_nomiss$PC1, merged_nomiss$PC2, merged_nomiss$PC3, merged_nomiss$PC4, merged_nomiss$PC5)

# Subdistribution hazard model for PD death (failcode = 1) 
crr.1 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 1, cencode = 0)
summary(crr.1)
#The coefficients from the summary(crr) are in the order that we specified in the cov.mat

# Subdistribution hazard model for nonPD death (failcode = 2)
crr.2 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 2, cencode = 0)
summary(crr.2)





#---Cumulative incidence curves for SYT10 rs10437796---####

#Summary of SNP allele carriers
merged_nomiss %>% 
  group_by(rs10437796) %>% 
  summarise(count = n())

#Code SNP in dominant model
merged_nomiss <- merged_nomiss %>% 
  mutate(rs10437796_dominant = ifelse(rs10437796 == 0, 0,
                                      ifelse(rs10437796 == 1 | rs10437796 == 2, 1, NA)))

merged_nomiss %>% 
  group_by(rs10437796_dominant) %>% 
  summarise(count = n())

#Check if any missing SNP data
merged_nomiss %>% 
  filter(is.na(rs144889025_dominant)) %>% 
  summarise(count = n())


#Make SNP into a factor
rs10437796_dominant_fact <- factor(merged_nomiss$rs10437796_dominant, levels=c(0,1), labels= c("rs10437796 non-carriers", "rs10437796 carriers"))

#Run analysis
fit <- CumIncidence_MT(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, rs10437796_dominant_fact, cencode = 0, xlab="Years from onset")

#---Subdistribution hazard models for SYT10 SNP rs10437796---####

#Make covariate matrix
cov.mat <- cbind(merged_nomiss$rs10437796_dominant, merged_nomiss$age_onset_imput, merged_nomiss$SEX,
                 merged_nomiss$PC1, merged_nomiss$PC2, merged_nomiss$PC3, merged_nomiss$PC4, merged_nomiss$PC5)

# Subdistribution hazard model for PD death (failcode = 1) 
crr.1 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 1, cencode = 0)
summary(crr.1)
#The coefficients from the summary(crr) are in the order that we specified in the cov.mat

# Subdistribution hazard model for nonPD death (failcode = 2)
crr.2 <- crr(merged_nomiss$timeToEvent_death, merged_nomiss$PDdeath, cov.mat, failcode = 2, cencode = 0)
summary(crr.2)




