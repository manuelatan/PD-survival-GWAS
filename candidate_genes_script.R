### PD RISK SNPS AND CANDIDATE VARIANT ANALYSIS ###

#---Load packages---####

library(readxl)
library(tidyverse)
library(data.table)

#---Read in PD GWAS data---####

#Read in PD GWAS summary stats in GRCh37/hg19
PD_SNPs <- read_excel("./PD_SNPs/Table S2. Detailed summary statistics.xlsx")

#Remove SNPs that failed final filtering
PD_SNPs_pass <- PD_SNPs %>% 
  filter(`Failed final filtering and QC` ==  0)

#---Export PD GWAS risk SNPs for GWAS analysis---####


#Export by chr:bp in hg19/GRCh37
export_GRS_chrbp <- PD_SNPs_pass %>% 
  mutate(CHRBP = paste(CHR, ":", BP, sep = ""),
         A1 = toupper(`Effect allele`)) %>% 
  select(CHRBP, A1, `Beta, all studies`) %>% 
  rename(beta = `Beta, all studies`)

#Export by rsID
export_GRS_rsid <- PD_SNPs_pass %>% 
  mutate(A1 = toupper(`Effect allele`)) %>% 
  select(SNP, A1, `Beta, all studies`) %>% 
  rename(beta = `Beta, all studies`)

write.table(export_GRS_chrbp, "./PD_GRS/Nalls_score_chrbp.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

write.table(export_GRS_rsid, "./PD_GRS/Nalls_score_rsid.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

#---Read in results from survival GWAS meta-analyses---####

mortality <- fread("../metaanalysis_mortality/covars_aao_gender/SURVIVAL_MORTALITY_META_20210907_hg191.tbl")

mortality_filtered <- mortality %>% 
  filter(TotalSampleSize >= 1000) %>% 
  select(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, `P-value`, ) %>% 
  separate(MarkerName, into = c("CHR", "BP"), convert = TRUE) %>% 
  rename(mortality_beta = Effect,
         mortality_pval = `P-value`)


HY3 <- fread("../metaanalysis_HY3/covars_aao_gender/SURVIVAL_HY3_META1.tbl")

HY3_filtered <- HY3 %>% 
  filter(TotalSampleSize >= 1000) %>% 
  select(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, `P-value`) %>%
  rename(SNP = MarkerName,
         HY3_beta = Effect,
         HY3_pval = `P-value`)

dementia <- fread("../metaanalysis_dementia/covars_aao_gender/SURVIVAL_DEMENTIA_META1.tbl")

dementia_filtered <- dementia %>%
  filter(TotalSampleSize >= 1000) %>% 
  select(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, `P-value`) %>% 
  rename(SNP = MarkerName,
         dementia_beta = Effect,
         dementia_pval = `P-value`)

#---Merge PD GWAS SNPs with survival GWAS SNPs---####

#Make CHR and BP into correct format in mortality data
mortality_filtered$CHR <- as.double(mortality_filtered$CHR)
mortality_filtered$BP <- as.double(mortality_filtered$BP)

#Merge PD risk SNPS with survival GWAS results
merged <- PD_SNPs_pass %>% 
  left_join(mortality_filtered, by = c("CHR", "BP")) %>%
  left_join(HY3_filtered, by = "SNP") %>% 
  left_join(dementia_filtered, by = "SNP")


#Check which SNPs are missing from the survival GWAS results

merged_nonmiss <- merged %>%
  filter(!is.na(mortality_beta), !is.na(HY3_beta), !is.na(dementia_beta))

#Looks like the only missing SNPs are those that are the rare GBA and LRRK2 variants

#---Write table of results---####

#Check that the effect alleles match in PD GWAS vs progression GWAS
merged_nonmiss %>%
  filter(`Effect allele`!=Allele1.x)

PDSNPs_export <- merged_nonmiss %>%
  mutate(effect_allele = toupper(Allele1.x),
         other_allele = toupper(Allele2.x)) %>% 
  select(`Nearest Gene`, SNP, CHR, BP, effect_allele, other_allele, Freq1.x, 
         mortality_beta, mortality_pval,
         HY3_beta, HY3_pval,
         dementia_beta, dementia_pval)

write.table(PDSNPs_export, "./PD_SNPs/PDSNPs_results.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

#---Correction for multiple testing---####

#This is within each domain, divided by the number of genes tested
analysis_wide_significance <- 0.05/88

#This is correcting across the 3 domains
testwide_significance <- analysis_wide_significance/3



#---Create variables for significance level---####
#p > 0.05
#p < 0.05
#p < 0.01
#p < analysis wide significance
#p < test wide significance

#Set scipen
options(scipen = 5)

merged_nonmiss <- merged_nonmiss %>%
  mutate(Pvalue_mortality = ifelse(mortality_pval > 0.05, "p>0.05",
                                            ifelse(mortality_pval < 0.05 & mortality_pval > 0.01, "p<0.05",
                                                   ifelse(mortality_pval < 0.01 & mortality_pval > analysis_wide_significance, "p<0.01",
                                                          ifelse(mortality_pval < analysis_wide_significance & mortality_pval > testwide_significance, "p<0.00057",
                                                                 ifelse(mortality_pval < testwide_significance, "p<0.00019", NA)))))) %>% 
  mutate(Pvalue_HY3 = ifelse(HY3_pval > 0.05, "p>0.05",
                                        ifelse(HY3_pval < 0.05 & HY3_pval > 0.01, "p<0.05",
                                               ifelse(HY3_pval < 0.01 & HY3_pval > analysis_wide_significance, "p<0.01",
                                                      ifelse(HY3_pval < analysis_wide_significance & HY3_pval > testwide_significance, "p<0.00057",
                                                             ifelse(HY3_pval < testwide_significance, "p<0.00019", NA)))))) %>% 
  mutate(Pvalue_dementia = ifelse(dementia_pval > 0.05, "p>0.05",
                                      ifelse(dementia_pval < 0.05 & dementia_pval > 0.01, "p<0.05",
                                             ifelse(dementia_pval < 0.01 & dementia_pval > analysis_wide_significance, "p<0.01",
                                                    ifelse(dementia_pval < analysis_wide_significance & dementia_pval > testwide_significance, "p<0.00057",
                                                           ifelse(dementia_pval < testwide_significance, "p<0.00019", NA)))))) 


#Filter out variants that have no significant associations
merged_nonmiss_filtered <- merged_nonmiss %>% 
  filter(mortality_pval < 0.05 | HY3_pval < 0.05 | dementia_pval < 0.05)

#Save intermediate merged dataset
write.table(merged_nonmiss, "./PD_SNPs/merged_file.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#---Convert to long format---####

merged_filtered_long <- gather(data = merged_nonmiss_filtered, key = phenotype, value = significance,
                               Pvalue_mortality:Pvalue_dementia)

#Check values for significance
merged_filtered_long %>% 
  group_by(significance) %>% 
  summarise(count = n())

#---Plot heat map for PD GWAS risk loci---####

#Colours for heatmap (purple)
colours <- c('#6a51a3', '#9e9ac8', '#cbc9e2', '#f2f0f7')

#Plot heatmap
ggplot(data = merged_filtered_long, mapping = aes(x = phenotype, y = SNP, fill = significance)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = colours) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_x_discrete(limits = c("Pvalue_dementia",  "Pvalue_HY3", "Pvalue_mortality"),
                   labels=c("Pvalue_mortality" = "mortality", 
                            "Pvalue_HY3" = "HY3", 
                            "Pvalue_dementia" = "dementia")) +
  coord_flip() + #rotate plot so that SNPs are on the x-axis and GWAS phenotype is on the y-object
  ggsave("./PD_SNPs/heatmap_PD_SNPs.png", width = 7.5, height = 2.5)



### OTHER CANDIDATE VARIANTS ####
#---Make list of candidate variants from other studies---####

#Read in candidate variant list and positions
#SLC44A1 from Iwaki, RIMS2 from Liu, TMEM108 from Liu, WWOX from Liu, APOE e2, MAPT H1
candidate_vars <- read.table("./other_candidate_vars/candidate_vars.txt", header = T)

#Separate chr and bp into separate columns
candidate_vars <- candidate_vars %>% 
  separate(`chr.bp`, into = c("CHR", "BP"), sep = ":")

#Make into double format (numeric)
candidate_vars$CHR <- as.double(candidate_vars$CHR)
candidate_vars$BP <- as.double(candidate_vars$BP)

#---Merge candidate variant list with survival GWAS results---####

#Merge PD risk SNPS with survival GWAS results
candidate_vars_merged <- candidate_vars %>% 
  left_join(mortality_filtered, by = c("CHR", "BP")) %>%
  left_join(HY3_filtered, by = c("rsid" = "SNP")) %>% 
  left_join(dementia_filtered, by = c("rsid" = "SNP"))


#---Correction for multiple testing---####

#This is within each domain, divided by the number of variants tested
CV_analysis_wide_significance <- 0.05/7

#This is correcting across the 3 domains
CV_testwide_significance <- analysis_wide_significance/3

#---Create variables for significance level---####
#p > 0.05
#p < 0.05
#p < 0.01
#p < analysis wide significance
#p < test wide significance

candidate_vars_merged <- candidate_vars_merged %>%
  mutate(Pvalue_mortality = ifelse(mortality_pval > 0.05, "p>0.05",
                                   ifelse(mortality_pval < 0.05 & mortality_pval > 0.01, "p<0.05",
                                          ifelse(mortality_pval < 0.01 & mortality_pval > CV_analysis_wide_significance, "p<0.01",
                                                 ifelse(mortality_pval < CV_analysis_wide_significance & mortality_pval > CV_testwide_significance, "p<0.0071",
                                                        ifelse(mortality_pval < CV_testwide_significance, "p<0.0028", NA)))))) %>% 
  mutate(Pvalue_HY3 = ifelse(HY3_pval > 0.05, "p>0.05",
                             ifelse(HY3_pval < 0.05 & HY3_pval > 0.01, "p<0.05",
                                    ifelse(HY3_pval < 0.01 & HY3_pval > CV_analysis_wide_significance, "p<0.01",
                                           ifelse(HY3_pval < CV_analysis_wide_significance & HY3_pval > CV_testwide_significance, "p<0.0071",
                                                  ifelse(HY3_pval < CV_testwide_significance, "p<0.0028", NA)))))) %>% 
  mutate(Pvalue_dementia = ifelse(dementia_pval > 0.05, "p>0.05",
                                  ifelse(dementia_pval < 0.05 & dementia_pval > 0.01, "p<0.05",
                                         ifelse(dementia_pval < 0.01 & dementia_pval > CV_analysis_wide_significance, "p<0.01",
                                                ifelse(dementia_pval < CV_analysis_wide_significance & dementia_pval > CV_testwide_significance, "p<0.0071",
                                                       ifelse(dementia_pval < CV_testwide_significance, "p<0.0028", NA)))))) 


#---Convert to long format for heatmap---####

candidate_vars_merged_long <- gather(data = candidate_vars_merged, key = phenotype, value = significance,
                               Pvalue_mortality:Pvalue_dementia)

#Check values for significance
candidate_vars_merged_long %>% 
  group_by(significance) %>% 
  summarise(count = n())


#---Export results---####

#If frequency is > 0.5, flip alleles and effect direction
candidate_vars_export <- candidate_vars_merged %>% 
  mutate(effect_allele = ifelse(Freq1.x > 0.5, toupper(Allele2.x), toupper(Allele1.x)),
         other_allele = ifelse(Freq1.x > 0.5, toupper(Allele1.x), toupper(Allele2.x)),
         freq = ifelse(Freq1.x > 0.5, 1-Freq1.x, Freq1.x),
         beta_mortality = ifelse(Freq1.x > 0.5, -mortality_beta, mortality_beta),
         beta_HY3 = ifelse(Freq1.x > 0.5, -HY3_beta, HY3_beta),
         beta_dementia = ifelse(Freq1.x > 0.5, -dementia_beta, dementia_beta)) %>% 
  select(rsid, CHR, BP, effect_allele, other_allele, freq, 
         beta_mortality, mortality_pval, 
         beta_HY3, HY3_pval,
         beta_dementia, dementia_pval)

write.table(candidate_vars_export, "./other_candidate_vars/candidate_vars_results.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")
