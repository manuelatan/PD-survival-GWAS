### PD GRS ANALYSIS ####


#---Load packages---####
library(tidyverse)
library(data.table)
library(survival)
library(meta)

#---Analyse GRS vs mortality in each cohort---####

cohorts_list <- c("Aasly", "Calypso", "Cambridge", "DIGPD", "Oslo", "Oxford", "PROBAND", "QSBB", "PPMI", "UKB")

#Make results dataframe
coefficients<-as.data.frame(matrix(ncol= 9))
names(coefficients) <- c("cohort","Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio","logrank", "r2" )


for (i in 1:length(cohorts_list)) {
  
  name <- cohorts_list[i]
  
  #Read in GRS
  GRS <- fread(paste("./PD_GRS/", name, ".GRS.profile", sep = ""))
  
  #Read in principal components
  PCs_location <- paste("../", name, "/PCA.eigenvec", sep = "")
  PCs <- fread(PCs_location)

  #Need to change the FID and IID in the PCs as these have been calculated in unimputed data
  
  if (name == "Calypso" | name == "Cambridge" | name == "QSBB" | name == "PPMI" | name == "UKB") {
    
    colnames(PCs) <- c("V1", "V2", paste("PC", 1:20, sep = ""))
    
    PCs <- PCs %>% 
      mutate(FID = V1,
             IID = FID) %>% 
      select(-V1, -V2)
  } else if (name == "Oslo") {
    #If reading in the Oslo data, the PC file is a bit different
    #Colnames are already in the file
    PCs <- PCs %>% 
      mutate(FID = IID)
  
  } else {
    
    colnames(PCs) <- c("V1", "V2", paste("PC", 1:20, sep = ""))
    
    PCs <- PCs %>% 
      mutate(FID = paste(V1, "_", V2, sep = ""),
             IID = FID) %>% 
      select(-V1, -V2)
    
  }
  
  #Read in clinical data for mortality
  #First get the filename of the mortality clinical data file
  mortality_filename <- list.files(path =  paste("../", name, "/mortality/covars_aao_gender/", sep = ""),
                      pattern = "2021-|2020-")
  
  #Read in mortality data
  mortality <- fread(paste("../", name, "/mortality/covars_aao_gender/", mortality_filename, sep = ""))
  
  #Merge all data
  merged <- mortality %>% 
    inner_join(GRS, by = c("FID", "IID")) %>% 
    inner_join(PCs, by = c("FID", "IID"))
  
  #The age at onset variables are sometimes named age_onset or age_onset_imput
  #Rename to age_onset
  merged <- merged %>% 
    rename_all(recode, age_onset_imput = "age_onset", age_diagnosis = "age_onset")
  
  #If UKB cohort, need to split into incident and prevalent before analysing GRS
  if (name == "UKB") { 
    
    UKB_subcohorts <- c("prevalent", "incident")
    
    for (j in 1:length(UKB_subcohorts)) {
    
      UK_subcohort <- merged %>% 
        filter(PD_status == UKB_subcohorts[j])
    
      #Standardise GRS
      mean_GRS <- mean(UK_subcohort$SCORE)
      sd_GRS <- sd(UK_subcohort$SCORE)
      
      UK_subcohort <- UK_subcohort %>% 
        mutate(zGRS = (SCORE - mean_GRS)/sd_GRS)
      
      #Analyse mortality vs. GRS
      model.cox <- coxph(Surv(UK_subcohort$timeToEvent_death, UK_subcohort$event_death) ~ zGRS + UK_subcohort$age_onset + UK_subcohort$gender 
                         + UK_subcohort$PC1+ UK_subcohort$PC2 + UK_subcohort$PC3 + UK_subcohort$PC4 + UK_subcohort$PC5, data=UK_subcohort)
      
      summary(model.cox)
      
      kmz <- cox.zph(model.cox, transform = "km")
      
      k = j+9
      
      coefficients[k,1]<- paste("UKB", UKB_subcohorts[j], sep = "")
      coefficients[k,2]<- summary(model.cox)$coefficients[1,1]
      coefficients[k,3]<- summary(model.cox)$coefficients[1,3]
      coefficients[k,4]<- summary(model.cox)$coefficients[1,5]
      coefficients[k,5]<- kmz$table[1,3]
      coefficients[k,6]<- model.cox$n
      coefficients[k,7]<- summary(model.cox)$logtest[[1]]
      coefficients[k,8]<- summary(model.cox)$sctest[[1]]
      coefficients[k,9]<- summary(model.cox)$rsq[[1]] 
      } 
    
    }
  else { 
      
      #FOR NON UKB COHORTS
      
      #Standardise GRS
      mean_GRS <- mean(merged$SCORE)
      sd_GRS <- sd(merged$SCORE)
      
      merged <- merged %>% 
        mutate(zGRS = (SCORE - mean_GRS)/sd_GRS)
      
      #Analyse mortality vs. GRS
      model.cox <- coxph(Surv(merged$timeToEvent_death, merged$event_death) ~ zGRS + merged$age_onset + merged$gender 
                         + merged$PC1+ merged$PC2 + merged$PC3 + merged$PC4 + merged$PC5, data=merged)
      
      summary(model.cox)
      
      kmz <- cox.zph(model.cox, transform = "km")
      
      coefficients[i,1]<- name
      coefficients[i,2]<- summary(model.cox)$coefficients[1,1]
      coefficients[i,3]<- summary(model.cox)$coefficients[1,3]
      coefficients[i,4]<- summary(model.cox)$coefficients[1,5]
      coefficients[i,5]<- kmz$table[1,3]
      coefficients[i,6]<- model.cox$n
      coefficients[i,7]<- summary(model.cox)$logtest[[1]]
      coefficients[i,8]<- summary(model.cox)$sctest[[1]]
      coefficients[i,9]<- summary(model.cox)$rsq[[1]] # nagelkerke r square
    }
    
  }




#---Meta analyse GRS vs mortality---####

#Random effects meta-analysis
meta_mortality <- metagen(TE = Coeff,
                 seTE = se,
                 studlab = cohort,
                 data = coefficients,
                 sm = "HR",
                 comb.fixed = TRUE,
                 comb.random = TRUE)

meta_mortality

#---Analyse GRS vs HY3 in each cohort---####

cohorts_HY3 <- c("DIGPD", "Oslo", "Oxford", "PROBAND", "PPMI")

#Make results dataframe
coefficients_HY3 <-as.data.frame(matrix(ncol= 9))
names(coefficients_HY3) <- c("cohort","Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio","logrank", "r2" )


for (i in 1:length(cohorts_HY3)) {
  
  name <- cohorts_HY3[i]
  
  #Read in GRS
  GRS <- fread(paste("./PD_GRS/", name, ".GRS.profile", sep = ""))
  
  #Read in principal components
  PCs_location <- paste("../", name, "/PCA.eigenvec", sep = "")
  PCs <- fread(PCs_location)
  
  #Need to change the FID and IID in the PCs as these have been calculated in unimputed data
  
  if (name == "PPMI") {
    
    colnames(PCs) <- c("V1", "V2", paste("PC", 1:20, sep = ""))
    
    PCs <- PCs %>% 
      mutate(FID = V1,
             IID = FID) %>% 
      select(-V1, -V2)
  } else if (name == "Oslo") {
    #If reading in the Oslo data, the PC file is a bit different
    #Colnames are already in the file
    PCs <- PCs %>% 
      mutate(FID = IID)
    
  } else {
    
    colnames(PCs) <- c("V1", "V2", paste("PC", 1:20, sep = ""))
    
    PCs <- PCs %>% 
      mutate(FID = paste(V1, "_", V2, sep = ""),
             IID = FID) %>% 
      select(-V1, -V2)
    
  }
  
  #Read in clinical data for mortality
  #First get the filename of the mortality clinical data file
  HY3_filename <- list.files(path =  paste("../", name, "/HY3/covars_aao_gender/", sep = ""),
                             pattern = "2021-|2020-")
  
  #Read in HY3 data
  HY3 <- fread(paste("../", name, "/HY3/covars_aao_gender/", HY3_filename, sep = ""))
  
  #Merge all data
  merged <- HY3 %>% 
    inner_join(GRS, by = c("FID", "IID")) %>% 
    inner_join(PCs, by = c("FID", "IID"))
  
  #The age at onset variables are sometimes named age_onset or age_onset_imput
  #Rename to age_onset
  merged <- merged %>% 
    rename_all(recode, age_onset_imput = "age_onset", age_diagnosis = "age_onset")
  
  #Standardise GRS
  mean_GRS <- mean(merged$SCORE)
  sd_GRS <- sd(merged$SCORE)
  
  merged <- merged %>% 
    mutate(zGRS = (SCORE - mean_GRS)/sd_GRS)
  
  #Analyse mortality vs. GRS
  model.cox <- coxph(Surv(merged$timeToEvent_HY3, merged$event_HY3) ~ zGRS + merged$age_onset + merged$gender 
                     + merged$PC1+ merged$PC2 + merged$PC3 + merged$PC4 + merged$PC5, data=merged)
  
  summary(model.cox)
  
  kmz <- cox.zph(model.cox, transform = "km")
  
  coefficients_HY3[i,1]<- name
  coefficients_HY3[i,2]<- summary(model.cox)$coefficients[1,1]
  coefficients_HY3[i,3]<- summary(model.cox)$coefficients[1,3]
  coefficients_HY3[i,4]<- summary(model.cox)$coefficients[1,5]
  coefficients_HY3[i,5]<- kmz$table[1,3]
  coefficients_HY3[i,6]<- model.cox$n
  coefficients_HY3[i,7]<- summary(model.cox)$logtest[[1]]
  coefficients_HY3[i,8]<- summary(model.cox)$sctest[[1]]
  coefficients_HY3[i,9]<- summary(model.cox)$rsq[[1]] # nagelkerke r square
  
}

#---Meta-analyse GRS vs HY3---####


#Random effects meta-analysis
meta_HY3 <- metagen(TE = Coeff,
                          seTE = se,
                          studlab = cohort,
                          data = coefficients_HY3,
                          sm = "HR",
                          comb.fixed = TRUE,
                          comb.random = TRUE)

meta_HY3



#---Analyse GRS vs dementia in each cohort---####

cohorts_dementia <- c("Aasly", "DIGPD", "Oslo", "Oxford", "PROBAND", "PPMI")

#Make results dataframe
coefficients_dementia <-as.data.frame(matrix(ncol= 9))
names(coefficients_dementia) <- c("cohort","Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio","logrank", "r2" )


for (i in 1:length(cohorts_dementia)) {
  
  name <- cohorts_dementia[i]
  
  #Read in GRS
  GRS <- fread(paste("./PD_GRS/", name, ".GRS.profile", sep = ""))
  
  #Read in principal components
  PCs_location <- paste("../", name, "/PCA.eigenvec", sep = "")
  PCs <- fread(PCs_location)
  
  #Need to change the FID and IID in the PCs as these have been calculated in unimputed data
  
  if (name == "PPMI") {
    
    colnames(PCs) <- c("V1", "V2", paste("PC", 1:20, sep = ""))
    
    PCs <- PCs %>% 
      mutate(FID = V1,
             IID = FID) %>% 
      select(-V1, -V2)
  } else if (name == "Oslo") {
    #If reading in the Oslo data, the PC file is a bit different
    #Colnames are already in the file
    PCs <- PCs %>% 
      mutate(FID = IID)
    
  } else {
    
    colnames(PCs) <- c("V1", "V2", paste("PC", 1:20, sep = ""))
    
    PCs <- PCs %>% 
      mutate(FID = paste(V1, "_", V2, sep = ""),
             IID = FID) %>% 
      select(-V1, -V2)
    
  }
  
  #Read in clinical data for dementia
  #First get the filename of the dementia clinical data file
  dementia_filename <- list.files(path =  paste("../", name, "/dementia/covars_aao_gender/", sep = ""),
                                  pattern = "2021-|2020-")
  
  #Read in dementia data
  dementia <- fread(paste("../", name, "/dementia/covars_aao_gender/", dementia_filename, sep = ""))
  
  #Merge all data
  merged <- dementia %>% 
    inner_join(GRS, by = c("FID", "IID")) %>% 
    inner_join(PCs, by = c("FID", "IID"))
  
  #The age at onset variables are sometimes named age_onset or age_onset_imput
  #Rename to age_onset
  merged <- merged %>% 
    rename_all(recode, age_onset_imput = "age_onset", age_diagnosis = "age_onset")
  
  #Standardise GRS
  mean_GRS <- mean(merged$SCORE)
  sd_GRS <- sd(merged$SCORE)
  
  merged <- merged %>% 
    mutate(zGRS = (SCORE - mean_GRS)/sd_GRS)
  
  #Analyse mortality vs. GRS
  model.cox <- coxph(Surv(merged$timeToEvent_dementia, merged$event_dementia) ~ zGRS + merged$age_onset + merged$gender 
                     + merged$PC1+ merged$PC2 + merged$PC3 + merged$PC4 + merged$PC5, data=merged)
  
  summary(model.cox)
  
  kmz <- cox.zph(model.cox, transform = "km")
  
  coefficients_dementia[i,1]<- name
  coefficients_dementia[i,2]<- summary(model.cox)$coefficients[1,1]
  coefficients_dementia[i,3]<- summary(model.cox)$coefficients[1,3]
  coefficients_dementia[i,4]<- summary(model.cox)$coefficients[1,5]
  coefficients_dementia[i,5]<- kmz$table[1,3]
  coefficients_dementia[i,6]<- model.cox$n
  coefficients_dementia[i,7]<- summary(model.cox)$logtest[[1]]
  coefficients_dementia[i,8]<- summary(model.cox)$sctest[[1]]
  coefficients_dementia[i,9]<- summary(model.cox)$rsq[[1]] # nagelkerke r square
  
}




#---Meta-analyse GRS vs dementia---####


#Random effects meta-analysis
meta_dementia <- metagen(TE = Coeff,
                    seTE = se,
                    studlab = cohort,
                    data = coefficients_dementia,
                    sm = "HR",
                    comb.fixed = TRUE,
                    comb.random = TRUE)

meta_dementia



