###### MAKE FOREST PLOTS OF TOP SNP EFFECTS IN EACH COHORT #####

#---Load packages---####
library(dplyr)
library(stringr)
library(grid)
library(forestplot)
library(meta)

#---Make functions---####

#Make forest plot function
#Using forestplot package and function
#SNP is the SNP name
#cliplower and clipupper are the limits of the x-axis (the hazard ratio) - as this varies between SNPs
make_forest_plot <- function(SNP, cliplower, clipupper) {
  
  #Filename
  filename <- paste("./", SNP, "_cohorts.txt", sep = "")
  
  #Read in results into R
  data <- read.table(filename)
  
  #Rename columns
  data <- data %>% 
    rename(file = "V1",
           effect_allele = "V2",
           noneffect_allele = "V3",
           beta = "V4",
           se = "V5",
           Pvalue = "V6",
           N = "V7",
           MAF = "V8")
  
  #Correct allele if noneffect allele T is reading in as TRUE
  data <- data %>% 
    mutate(noneffect_allele = ifelse(noneffect_allele == TRUE, "T", NA))
  
  #Add cohort column and calculate SD from SE
  data <- data %>% 
    mutate(cohort = ifelse(str_detect(file, "PROBAND"), "Tracking Parkinson's", 
                           ifelse(str_detect(file, "Oxford"), "Oxford Discovery",
                                  ifelse(str_detect(file, "QSBB"), "QSBB",
                                         ifelse(str_detect(file, "incident"), "UKB_incident",
                                                ifelse(str_detect(file, "prevalent"), "UKB_prevalent",
                                                       ifelse(str_detect(file, "Calypso"), "Calypso",
                                                              ifelse(str_detect(file, "CamPaIGN"), "CamPaIGN",
                                                                     ifelse(str_detect(file, "DIGPD"), "DIGPD",
                                                                            ifelse(str_detect(file, "Aasly"), "Trondheim",
                                                                                   ifelse(str_detect(file, "Oslo"), "Oslo",
                                                                                          ifelse(str_detect(file, "PPMI"), "PPMI", NA))))))))))),
           hazard_ratio = exp(beta), #Calculate hazard ratio
           CI_95_lower = exp(beta - 1.96*se), #Calculate 95% CI lower bound
           CI_95_upper =  exp(beta + 1.96*se), #Calculate 95% CI upper bound
           HR = sprintf("%.2f", round(hazard_ratio, 2)), #Create new text column with HR, keeping trailing 0s
           CI_95 = paste("[", sprintf("%.2f", round(CI_95_lower, 2)), "-", sprintf("%.2f", round(CI_95_upper, 2)), "]", sep = ""))
  #Make column for 95% CI combining lower and upper values
  
  
  #Select relevant columns for forest plot
  base <- data %>% 
    select(study = cohort, 
           mean = hazard_ratio, 
           lower = CI_95_lower,
           upper = CI_95_upper,
           HR = HR,
           CI_95 = CI_95)
  
  #Read in METAL meta-analysis results
  metal_filename <- paste(SNP, "_metal.txt", sep = "")
  
  metal <- read.table(metal_filename)
  colnames(metal) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "Pvalue", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize")
  
  #Swap alleles and effect direction if the freq > 0.5
  metal <- metal %>% 
    mutate(A1 = ifelse(Freq1 > 0.5, Allele2, Allele1),
           A2 = ifelse(Freq1 > 0.5, Allele1, Allele2),
           freq = ifelse(Freq1 > 0.5, 1-Freq1, Freq1),
           beta = ifelse(Freq1 > 0.5, -Effect, Effect),
           hazard_ratio = exp(beta),
           CI_95_lower = exp(beta - 1.96*StdErr),
           CI_95_upper =  exp(beta + 1.96*StdErr),
           HR = sprintf("%.2f", round(hazard_ratio, 2)), #Create column with HR, keeping trailing 0s
           CI_95 = paste("[", sprintf("%.2f", round(CI_95_lower, 2)), "-", sprintf("%.2f", round(CI_95_upper, 2)), "]", sep = ""),
           study = "Meta-analysis")
  
  #Format summary meta analysis results for plot
  summary <- metal %>%
    select(mean = hazard_ratio,
           lower = CI_95_lower,
           upper = CI_95_upper,
           study = study,
           HR = HR,
           CI_95 = CI_95) %>% 
    tibble(summary = TRUE)
  
  #Make header dataframe
  header <- tibble(study = c("Study"),
                   HR = c("HR"),
                   CI_95 = c("95% CI"),
                   summary = TRUE)
  
  #Combine
  forestplot_df <- bind_rows(header, base, summary)
  
  #Make forest plot file name
  plotname <- paste("./plots/", SNP, "_forest.png", sep = "")
  
  #Make forest plot filename
  png(plotname, width = 7, height = 5, units = "in", res = 300)
  
  #Calculate the number of rows in the forestplot_df - this changes depending on the number of cohorts included
  #This will determine where the horizontal line is placed. We want it to be just before the meta-analysis results
  list_lines <- list("2" = gpar(lwd = 1), "A" = gpar(lwd = 1)) 
  names(list_lines) <- c("2", paste(nrow(forestplot_df)))
  
  #Make list of xticks
  xticks_list <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 5, 7, 10)
  
  #Filter according to specified limits
  xticks_list_filtered <- xticks_list[xticks_list <= clipupper]
  xticks_list_filtered <- xticks_list_filtered[xticks_list_filtered >= cliplower]
  
  #Make forest plot
  plot <- forestplot_df %>% 
    forestplot(labeltext = c(study, HR, CI_95), 
               is.summary = summary,
               hrzl_lines = list_lines,
               boxsize = 0.2,
               xlog = TRUE, 
               clip = c(cliplower, clipupper),
               xticks = xticks_list_filtered,
               line.margin = .1,
               xlab = "Hazard Ratio",
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "red",
                              hrz_lines = "#444444"),
               txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 0.8),
                                xlab  = gpar(fontfamily = "", cex = 1.2)))
  print(plot)
  dev.off()
}


#---Extract GGT top SNP effects from each cohort and plot forest plot---####

# Top SNP is rs112809886, 22:24641838

#rs112809886, 22:24641838
#Excluding Oslo as this had genomic inflation factor > 1.2
grep rs112809886 \
/data/kronos/kronos/mtan/survival_GWAS/PROBAND/HY3/covars_aao_gender/PROBAND_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/Oxford/HY3/covars_aao_gender/Oxford_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/PPMI/HY3/covars_aao_gender/PPMI_survival_HY3_METAL.tab \
/data/kronos/kronos/mtan/survival_GWAS/DIGPD/HY3/covars_aao_gender/DIGPD_survival_HY3_METAL.tab > rs112809886_cohorts.txt


#Get meta-analysis results
grep rs112809886 /data/kronos/kronos/mtan/survival_GWAS/metaanalysis_HY3/covars_aao_gender/SURVIVAL_HY3_META1.tbl > rs112809886_metal.txt

make_forest_plot("rs112809886", 0.5, 10)


#---Including proxy in Tracking Parkinson's (PROBAND)---####

#Proxy for rs112809886 in PROBAND = rs148492066 (D' = 1 and R2 = 1)
#grep rs148492066 /data/kronos/kronos/mtan/survival_GWAS/PROBAND/HY3/covars_aao_gender/PROBAND_survival_HY3_METAL.tab > rs148492066_PROBAND_proxy.txt

#Read in PROBAND proxy results
rs112809886_PROBAND_proxy <- read.table("./rs148492066_PROBAND_proxy.txt")

#Read in results into R
data <- read.table("./rs112809886_cohorts.txt")

#Combine
data_w_proxy <- rbind(data, rs112809886_PROBAND_proxy)

#Rename columns
data_w_proxy <- data_w_proxy %>% 
  rename(file = "V1",
         effect_allele = "V2",
         noneffect_allele = "V3",
         beta = "V4",
         se = "V5",
         Pvalue = "V6",
         N = "V7",
         MAF = "V8")

#Add cohort column and calculate SD from SE
data_w_proxy <- data_w_proxy %>% 
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
  mutate(hazard_ratio = exp(beta), #Calculate hazard ratio
         CI_95_lower = exp(beta - 1.96*se), #Calculate 95% CI lower bound
         CI_95_upper =  exp(beta + 1.96*se), #Calculate 95% CI upper bound
         HR = sprintf("%.2f", round(hazard_ratio, 2)), #Create new text column with HR, keeping trailing 0s
         CI_95 = paste("[", sprintf("%.2f", round(CI_95_lower, 2)), "-", sprintf("%.2f", round(CI_95_upper, 2)), "]", sep = ""))
#Make column for 95% CI combining lower and upper values


#Do meta-analysis using random effects model
#Using Knapp-Hartung(-Sidik-Jonkman) adjustment
meta_rs112809886_wPROBAND_proxy <- metagen(beta,
                            se,
                            data = data_w_proxy,
                            studlab = cohort,
                            comb.fixed = TRUE,
                            comb.random = TRUE,
                            prediction = TRUE,
                            sm = "HR")

meta_rs112809886_wPROBAND_proxy

#Format summary meta analysis results for plot
summary <- data.frame(matrix(ncol = 6))
colnames(summary) <- c("hazard_ratio", "CI_95_lower", "CI_95_upper", "study", "HR", "CI_95")

summary[1, "hazard_ratio"] <- exp(meta_rs112809886_wPROBAND_proxy$TE.random) # Hazard Ratio which is exp of beta
summary[1, "CI_95_lower"] <- exp(meta_rs112809886_wPROBAND_proxy$lower.random) # Hazard Ratio lower 95% CI bound
summary[1, "CI_95_upper"] <- exp(meta_rs112809886_wPROBAND_proxy$upper.random) # Hazard Ratio upper 95% CI bound

summary <- summary %>% 
  mutate(HR = sprintf("%.2f", round(hazard_ratio, 2)), #Create column with HR, keeping trailing 0s
         CI_95 = paste("[", sprintf("%.2f", round(CI_95_lower, 2)), "-", sprintf("%.2f", round(CI_95_upper, 2)), "]", sep = ""),
         study = "Meta-analysis")

#Prepare for plot
summary <- summary %>%
  select(mean = hazard_ratio,
         lower = CI_95_lower,
         upper = CI_95_upper,
         study = study,
         HR = HR,
         CI_95 = CI_95) %>% 
  tibble(summary = TRUE)

#Make header dataframe
header <- tibble(study = c("Study"),
                 HR = c("HR"),
                 CI_95 = c("95% CI"),
                 summary = TRUE)

#Make base dataframe - individual cohort results.Select relevant columns for forest plot
base <- data_w_proxy %>% 
  select(study = cohort, 
         mean = hazard_ratio, 
         lower = CI_95_lower,
         upper = CI_95_upper,
         HR = HR,
         CI_95 = CI_95)

#Combine
forestplot_df <- bind_rows(header, base, summary)

#Make forest plot filename
png("./plots/rs112809886_forest_w_proxy.png", width = 7, height = 5, units = "in", res = 300)

#Calculate the number of rows in the forestplot_df - this changes depending on the number of cohorts included
#This will determine where the horizontal line is placed. We want it to be just before the meta-analysis results
list_lines <- list("2" = gpar(lwd = 1), "A" = gpar(lwd = 1)) 
names(list_lines) <- c("2", paste(nrow(forestplot_df)))

#Make list of xticks
xticks_list <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 5, 7, 10, 12)

#Make forest plot
plot <- forestplot_df %>% 
  forestplot(labeltext = c(study, HR, CI_95), 
             is.summary = summary,
             hrzl_lines = list_lines,
             boxsize = 0.2,
             xlog = TRUE, 
             clip = c(0.1, 10.5),
             xticks = xticks_list,
             line.margin = .1,
             xlab = "Hazard Ratio",
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "red",
                            hrz_lines = "#444444"),
             txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 0.8),
                              xlab  = gpar(fontfamily = "", cex = 1.2)))
print(plot)
dev.off()

