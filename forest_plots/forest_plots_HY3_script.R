###### MAKE FOREST PLOTS OF TOP SNP EFFECTS IN EACH COHORT #####

#---Load packages---####
library(dplyr)
library(stringr)
library(meta)
library(grid)

#---Extract GGT top SNP effects from each cohort and plot forest plot---####

# Top SNP is rs112809886, 22:24641838

#rs112809886, 22:24641838
grep rs112809886 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/HY3/covars_aao_gender/PROBAND_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/HY3/covars_aao_gender/Oxford_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/PPMI/HY3/covars_aao_gender/PPMI_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/HY3/covars_aao_gender/DIGPD_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oslo/HY3/covars_aao_gender/Oslo_survival_HY3_METAL.tab > rs112809886_cohorts.txt

#Read in results into R
rs112809886 <- read.table("./rs112809886_cohorts.txt")

#Rename columns
rs112809886 <- rs112809886 %>% 
  rename(file = "V1",
         effect_allele = "V2",
         noneffect_allele = "V3",
         beta = "V4",
         se = "V5",
         Pvalue = "V6",
         N = "V7",
         MAF = "V8")


#Add cohort column and calculate SD from SE
rs112809886 <- rs112809886 %>% 
  mutate(cohort = ifelse(str_detect(file, "PROBAND"), "Tracking Parkinson's", 
                         ifelse(str_detect(file, "Oxford"), "Oxford Discovery",
                                ifelse(str_detect(file, "QSBB"), "QSBB",
                                       ifelse(str_detect(file, "incident"), "UKB_incident",
                                              ifelse(str_detect(file, "prevalent"), "UKB_prevalent",
                                                     ifelse(str_detect(file, "Calypso"), "Calypso",
                                                            ifelse(str_detect(file, "Cambridge"), "CamPaiGN",
                                                                   ifelse(str_detect(file, "DIGPD"), "DIGPD",
                                                                          ifelse(str_detect(file, "Aasly"), "Trondheim",
                                                                                 ifelse(str_detect(file, "Oslo"), "Oslo", 
                                                                                        ifelse(str_detect(file, "PPMI"), "PPMI",NA)))))))))))) %>% 
  mutate(sd = se * sqrt(N),
         HR = exp(beta)) 

#Remove PPMI because of the wide confidence intervals
rs112809886_noPPMI <- rs112809886 %>% 
  filter(cohort!="PPMI")


#Remove Oslo because lambda > 1.2
rs112809886_noOslo <- rs112809886 %>% 
  filter(cohort!="Oslo")

#Do meta-analysis using random effects model
#Using Knapp-Hartung(-Sidik-Jonkman) adjustment
meta_rs112809886 <- metagen(beta,
                           se,
                           data = rs112809886_noOslo,
                           studlab = cohort,
                           comb.fixed = TRUE,
                           comb.random = TRUE,
                           prediction = TRUE,
                           sm = "HR")

meta_rs112809886

#Try to adjust the values for PPMI as the conf intervals are too wide
meta_rs112809886$lower <- c(0.349, -3, 1.441)
meta_rs112809886$upper <- c(1.75, 3, 2.98)

#Make forest plot
png(file = './plots/forestplot_rs112809886.png', width = 700, height = 300)
forest(meta_rs112809886,
       layout = "meta",
       leftcols = c("studlab"),
       rightcols = c("effect.ci"),
       rightlabs = c("HR", "[95% CI]"),
       # print.tau2 = FALSE,
       # smlab = "",
       colgap.forest.left = unit(10,"mm"),
       print.tau2 = FALSE,
       prediction = FALSE,
       # colgap.forest.right = unit(5,"mm"),
       #comb.fixed = TRUE
       )
grid.text("Effect size for rs112809886 on HY3", 0.5, .98, 
          gp=gpar(fontsize = 6, cex=2, fontface = 2))
dev.off()



#---Including proxy in Tracking Parkinson's (PROBAND)---####

#Proxy for rs112809886 in PROBAND = rs148492066 (D' = 1 and R2 = 1)
#grep rs148492066 /data/kronos/kronos/mtan/survival_GWAS/PROBAND/HY3/covars_aao_gender/PROBAND_survival_HY3_METAL.tab > rs148492066_PROBAND_proxy.txt

#Read in PROBAND proxy results
rs112809886_PROBAND_proxy <- read.table("./rs148492066_PROBAND_proxy.txt")

#Rename columns
rs112809886_PROBAND_proxy <- rs112809886_PROBAND_proxy %>% 
  rename(file = "V1",
         effect_allele = "V2",
         noneffect_allele = "V3",
         beta = "V4",
         se = "V5",
         Pvalue = "V6",
         N = "V7",
         MAF = "V8")

#Merge with main data
rs112809886_with_PROBAND_proxy <- rbind(rs112809886, rs112809886_PROBAND_proxy)


#Add cohort column and calculate SD from SE
rs112809886_with_PROBAND_proxy <- rs112809886_with_PROBAND_proxy %>% 
  mutate(cohort = ifelse(str_detect(file, "PROBAND"), "Tracking Parkinson's", 
                         ifelse(str_detect(file, "Oxford"), "Oxford Discovery",
                                ifelse(str_detect(file, "QSBB"), "QSBB",
                                       ifelse(str_detect(file, "incident"), "UKB_incident",
                                              ifelse(str_detect(file, "prevalent"), "UKB_prevalent",
                                                     ifelse(str_detect(file, "Calypso"), "Calypso",
                                                            ifelse(str_detect(file, "Cambridge"), "CamPaiGN",
                                                                   ifelse(str_detect(file, "DIGPD"), "DIGPD",
                                                                          ifelse(str_detect(file, "Aasly"), "Trondheim",
                                                                                 ifelse(str_detect(file, "Oslo"), "Oslo", 
                                                                                        ifelse(str_detect(file, "PPMI"), "PPMI",
                                                                                               ifelse(str_detect(file, "rs148492066"), "Tracking Parkinson's proxy", NA))))))))))))) %>% 
  mutate(sd = se * sqrt(N),
         HR = exp(beta))

#Remove Oslo because lambda > 1.2
rs112809886_with_PROBAND_proxy_noOslo <- rs112809886_with_PROBAND_proxy %>% 
  filter(cohort!="Oslo")


#Do meta-analysis using random effects model
#Using Knapp-Hartung(-Sidik-Jonkman) adjustment
meta_rs112809886_wPROBAND_proxy <- metagen(beta,
                            se,
                            data = rs112809886_with_PROBAND_proxy_noOslo,
                            studlab = cohort,
                            comb.fixed = TRUE,
                            comb.random = TRUE,
                            prediction = TRUE,
                            sm = "HR")

meta_rs112809886_wPROBAND_proxy

#Try to adjust the values for PPMI as the conf intervals are too wide
meta_rs112809886_wPROBAND_proxy$lower <- c(0.349, -3, 1.441, -0.587)
meta_rs112809886_wPROBAND_proxy$upper <- c(1.745, 3, 2.987, 0.527)

#Make forest plot
png(file = './plots/forestplot_rs112809886_wPROBAND_proxy.png', width = 700, height = 300)
forest(meta_rs112809886_wPROBAND_proxy,
       layout = "meta",
       leftcols = c("studlab"),
       rightcols = c("effect.ci"),
       rightlabs = c("HR", "[95% CI]"),
       # print.tau2 = FALSE,
       # smlab = "",
       colgap.forest.left = unit(10,"mm"),
       print.tau2 = FALSE,
       prediction = FALSE,
       # colgap.forest.right = unit(5,"mm"),
       #comb.fixed = TRUE
)
grid.text("Effect size for rs112809886 on HY3", 0.5, .98, 
          gp=gpar(fontsize = 6, cex=2, fontface = 2))
dev.off()
