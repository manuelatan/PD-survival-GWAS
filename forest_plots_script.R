###### MAKE FOREST PLOTS OF TOP SNP EFFECTS IN EACH COHORT #####

#---Load packages---####
library(dplyr)
library(stringr)
library(meta)
library(grid)
library(meta)

#---Extract APOE top SNP effects from each cohort and plot forest plot---####

# Top SNP is rs429358, 19:45411941
# #Excluding PPMI because not enough people met the outcome

#rs429358, 19:45411941
grep 19:45411941 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/PROBAND_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/Oxford_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/QSBB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/UKB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/UKB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/Calypso_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Cambridge/mortality/covars_aao_gender/Cambridge_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/DIGPD_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/Aasly_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/Oslo_survival_mortality_METAL_hg19.tab > rs429358_cohorts.txt

#Get raw data for top APOE SNP
#This is just to check the Cox proportional hazards assumption
#Note that the SNP names are slightly different e.g. 19:45411941, 19.45411941
#So we are grepping with just the bp position
grep "45411941" \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Cambridge/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/allGWAS_results.txt > rs429358_raw_cohorts.txt

#Read in results into R
rs429358 <- read.table("./rs429358_cohorts.txt")

#Rename columns
rs429358 <- rs429358 %>% 
  rename(file = "V1",
         effect_allele = "V2",
         noneffect_allele = "V3",
         beta = "V4",
         se = "V5",
         Pvalue = "V6",
         N = "V7",
         MAF = "V8")

#Read in raw data
rs429358_raw <- read.table("./rs429358_raw_cohorts.txt")

rs429358_raw <- rs429358_raw %>% 
  rename(file = "V1",
         Coeff = "V2",
         se = "V3",
         Pvalue	= "V4",
         Cox.zphPVal = "V5",
         N = "V6",
         ov.lik.ratio	= "V7",
         logrank = "V8",
         r2 = "V9")

#Filter for exact SNP matching, grep does not seem to be matching exact string
rs429358_raw <- rs429358_raw %>% 
  filter(grepl("19.45411941|19:45411941", file))

#The noneffect allele T is reading in as TRUE
rs429358 <- rs429358 %>% 
  mutate(noneffect_allele = ifelse(noneffect_allele == TRUE, "T", NA))

#Add cohort column and calculate SD from SE
rs429358 <- rs429358 %>% 
  mutate(cohort = ifelse(str_detect(file, "PROBAND"), "Tracking Parkinson's", 
                         ifelse(str_detect(file, "Oxford"), "Oxford Discovery",
                                ifelse(str_detect(file, "QSBB"), "QSBB",
                                       ifelse(str_detect(file, "incident"), "UKB_incident",
                                              ifelse(str_detect(file, "prevalent"), "UKB_prevalent",
                                                     ifelse(str_detect(file, "Calypso"), "Calypso",
                                                            ifelse(str_detect(file, "Cambridge"), "CamPaiGN",
                                                                   ifelse(str_detect(file, "DIGPD"), "DIGPD",
                                                                          ifelse(str_detect(file, "Aasly"), "Trondheim",
                                                                                 ifelse(str_detect(file, "Oslo"), "Oslo", NA))))))))))) %>% 
  mutate(sd = se * sqrt(N),
         HR = exp(beta)) 

#Remove Oxford as the confidence intervals are too wide
#rs61871952_noOxford <- rs61871952 %>% 
#  filter(cohort!="Oxford Discovery")

#Do meta-analysis using random effects model
#Using Knapp-Hartung(-Sidik-Jonkman) adjustment
meta_rs429358 <- metagen(beta,
                           se,
                           data = rs429358,
                           studlab = cohort,
                           comb.fixed = TRUE,
                           comb.random = TRUE,
                           prediction = TRUE,
                           sm = "HR")

meta_rs429358

#Make forest plot
png(file = './plots/forestplot_rs429358.png', width = 700, height = 300) 
forest(meta_rs429358,
       layout = "meta",
       leftcols = c("studlab"),
       rightcols = c("effect.ci"),
       rightlabs = c("HR", "[95% CI]"),
       # print.tau2 = FALSE,
       # smlab = "",
       colgap.forest.left = unit(10,"mm"),
       print.tau2 = FALSE,
       prediction = FALSE
       # colgap.forest.right = unit(5,"mm"),
       #comb.fixed = TRUE
)
grid.text("Effect size for rs429358 on mortality", 0.5, .98, 
          gp=gpar(fontsize = 6, cex=2, fontface = 2))
dev.off()



#---Extract TBXAS1 top SNP effects from each cohort and plot forest plot---####

# Top SNP is rs4726467, 7:139637422
# #Excluding PPMI because not enough people met the outcome

grep 7:139637422 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/PROBAND_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/Oxford_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/QSBB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/UKB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/UKB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/Calypso_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Cambridge/mortality/covars_aao_gender/Cambridge_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/DIGPD_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/Aasly_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/Oslo_survival_mortality_METAL_hg19.tab > rs4726467_cohorts.txt

#Get raw data for top TBXAS1 SNP
grep 7:139637422 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Cambridge/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/allGWAS_results.txt > rs4726467_raw_cohorts.txt

#Read in results into R
rs4726467 <- read.table("./rs4726467_cohorts.txt")

#Rename columns
rs4726467 <- rs4726467 %>% 
  rename(file = "V1",
         effect_allele = "V2",
         noneffect_allele = "V3",
         beta = "V4",
         se = "V5",
         Pvalue = "V6",
         N = "V7",
         MAF = "V8")

#Read in raw data
rs4726467_raw <- read.table("./rs4726467_raw_cohorts.txt")

rs4726467_raw <- rs4726467_raw %>% 
  rename(file = "V1",
         Coeff = "V2",
         se = "V3",
         Pvalue	= "V4",
         Cox.zphPVal = "V5",
         N = "V6",
         ov.lik.ratio	= "V7",
         logrank = "V8",
         r2 = "V9")

#The effect allele T is reading in as TRUE
rs4726467 <- rs4726467 %>% 
  mutate(effect_allele = ifelse(effect_allele == TRUE, "T", NA))

#Add cohort column and calculate SD from SE
rs4726467 <- rs4726467 %>% 
  mutate(cohort = ifelse(str_detect(file, "PROBAND"), "Tracking Parkinson's", 
                         ifelse(str_detect(file, "Oxford"), "Oxford Discovery",
                                ifelse(str_detect(file, "QSBB"), "QSBB",
                                       ifelse(str_detect(file, "incident"), "UKB_incident",
                                              ifelse(str_detect(file, "prevalent"), "UKB_prevalent",
                                                     ifelse(str_detect(file, "Calypso"), "Calypso",
                                                            ifelse(str_detect(file, "Cambridge"), "CamPaiGN",
                                                                   ifelse(str_detect(file, "DIGPD"), "DIGPD",
                                                                          ifelse(str_detect(file, "Aasly"), "Trondheim",
                                                                                 ifelse(str_detect(file, "Oslo"), "Oslo", NA))))))))))) %>% 
  mutate(sd = se * sqrt(N),
         HR = exp(beta)) 


#Do meta-analysis using random effects model
meta_rs4726467 <- metagen(beta,
                         se,
                         data = rs4726467,
                         studlab = cohort,
                         comb.fixed = TRUE,
                         comb.random = TRUE,
                         prediction = TRUE,
                         sm = "HR")

meta_rs4726467

#Make forest plot
png(file = './plots/forestplot_rs4726467.png', width = 700, height = 300) 
forest(meta_rs4726467,
       layout = "meta",
       leftcols = c("studlab"),
       rightcols = c("effect.ci"),
       rightlabs = c("HR", "[95% CI]"),
       # print.tau2 = FALSE,
       # smlab = "",
       colgap.forest.left = unit(10,"mm"),
       print.tau2 = FALSE,
       prediction = FALSE
       # colgap.forest.right = unit(5,"mm"),
       #comb.fixed = TRUE
)
grid.text("Effect size for rs4726467 (TBXAS1) on mortality", 0.5, .98, 
          gp=gpar(fontsize = 6, cex=2, fontface = 2))
dev.off()



#---Extract SYT10 top SNP effects from each cohort and plot forest plot---####

# Top SNP is rs10437796, 12:33635494
# #Excluding PPMI because not enough people met the outcome

grep 12:33635494 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/PROBAND_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/Oxford_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/QSBB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/UKB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/UKB_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/Calypso_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Cambridge/mortality/covars_aao_gender/Cambridge_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/DIGPD_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/Aasly_survival_mortality_METAL_hg19.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/Oslo_survival_mortality_METAL_hg19.tab > rs10437796_cohorts.txt

#Get raw data for top SYT10 SNP
grep 33635494 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Cambridge/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/allGWAS_results.txt \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/allGWAS_results.txt > rs10437796_raw_cohorts.txt

#Read in results into R
rs10437796 <- read.table("./rs10437796_cohorts.txt")

#Rename columns
rs10437796 <- rs10437796 %>% 
  rename(file = "V1",
         effect_allele = "V2",
         noneffect_allele = "V3",
         beta = "V4",
         se = "V5",
         Pvalue = "V6",
         N = "V7",
         MAF = "V8")

#Read in raw data
rs10437796_raw <- read.table("./rs10437796_raw_cohorts.txt")

rs10437796_raw <- rs10437796_raw %>% 
  rename(file = "V1",
         Coeff = "V2",
         se = "V3",
         Pvalue	= "V4",
         Cox.zphPVal = "V5",
         N = "V6",
         ov.lik.ratio	= "V7",
         logrank = "V8",
         r2 = "V9")

#Filter for exact SNP matching, grep does not seem to be matching exact string
rs10437796_raw <- rs10437796_raw %>% 
  filter(grepl("12:33635494|12.33635494", file))

#Add cohort column and calculate SD from SE
rs10437796 <- rs10437796 %>% 
  mutate(cohort = ifelse(str_detect(file, "PROBAND"), "Tracking Parkinson's", 
                         ifelse(str_detect(file, "Oxford"), "Oxford Discovery",
                                ifelse(str_detect(file, "QSBB"), "QSBB",
                                       ifelse(str_detect(file, "incident"), "UKB_incident",
                                              ifelse(str_detect(file, "prevalent"), "UKB_prevalent",
                                                     ifelse(str_detect(file, "Calypso"), "Calypso",
                                                            ifelse(str_detect(file, "Cambridge"), "CamPaiGN",
                                                                   ifelse(str_detect(file, "DIGPD"), "DIGPD",
                                                                          ifelse(str_detect(file, "Aasly"), "Trondheim",
                                                                                 ifelse(str_detect(file, "Oslo"), "Oslo", NA))))))))))) %>% 
  mutate(sd = se * sqrt(N),
         HR = exp(beta)) 


#Do meta-analysis using random effects model
meta_rs10437796 <- metagen(beta,
                          se,
                          data = rs10437796,
                          studlab = cohort,
                          comb.fixed = TRUE,
                          comb.random = TRUE,
                          prediction = TRUE,
                          sm = "HR")

meta_rs10437796

#Make forest plot
png(file = './plots/forestplot_rs10437796.png', width = 700, height = 300) 
forest(meta_rs10437796,
       layout = "meta",
       leftcols = c("studlab"),
       rightcols = c("effect.ci"),
       rightlabs = c("HR", "[95% CI]"),
       # print.tau2 = FALSE,
       # smlab = "",
       colgap.forest.left = unit(10,"mm"),
       print.tau2 = FALSE,
       prediction = FALSE
       # colgap.forest.right = unit(5,"mm"),
       #comb.fixed = TRUE
)
grid.text("Effect size for rs10437796 (SYT10) on mortality", 0.5, .98, 
          gp=gpar(fontsize = 6, cex=2, fontface = 2))
dev.off()


