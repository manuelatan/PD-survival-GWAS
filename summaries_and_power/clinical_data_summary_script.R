### CLINICAL OUTCOME SUMMARIES ###
#Summarise the number of individuals who met the outcome across all cohorts

#---Load packages---####
library(tidyverse)
library(data.table)


#---Read in mortality clinical data from each cohort and merge---####

#Make list of cohorts - all for mortality analysis. Excluding PPMI
#Cambridge is included here even though it was exclued from meta-analysis - will remove later
cohorts_all <- c("Aasly", "Calypso", "Cambridge", "DIGPD",
                 "Oslo", "Oxford", "PROBAND", "QSBB", "UKB")

#Make empty data list
datalist <- list()

#Loop to read in all cohorts, get list of final individuals, and add to data list
for (i in 1:length(cohorts_all)){
  
  #Cohort name
  name <- cohorts_all[i]
  
  #Directory where clinical data is stored
  directory <- paste("../", name, "/mortality/covars_aao_gender/", sep = "")
  
  #Filename which contains 2020 or 2021
  mortality_filename <- list.files(path = directory,
                                   pattern = "2021-|2020-")
  
  #Read in mortality data
  data <- fread(paste("../", name, "/mortality/covars_aao_gender/",
                      mortality_filename,
                      sep = ""))
  
  if (name == "UKB") {
    prevalent <- fread(paste("../", name, "/mortality/", name, "prevalent_final_keep.txt", sep = ""))
    incident <-  fread(paste("../", name, "/mortality/", name, "incident_final_keep.txt", sep = ""))
    
    prevalent <- prevalent %>% 
      mutate(cohort = "UKB_prevalent")
    
    incident <- incident %>% 
      mutate(cohort = "UKB_incident")
    
    final_individuals <- rbind(prevalent, incident)
    
    #Read in genetic data fam file
    genetic <- fread("../UKB/UKB_PD.snpqc.fam")
    colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")

    data <- data %>%
      inner_join(final_individuals, by = "IID") %>% 
      inner_join(genetic, by = c("FID", "IID"))
    
  } else {
    #Read in final list of individuals - if it exists
    individuals_filename <- paste("../", name, "/mortality/", name, "_final_keep.txt", sep = "")
    
    #Merge to get final individuals in analysis
    if (file.exists(individuals_filename)) {
      
      #Read in final individuals list
      final_individuals <- fread(individuals_filename)
      
      #Read in genetic fam file
      genetic <- fread(paste("../", name, "/", name, "_PD.snpqc.fam", sep = ""))
      colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")
      
      data <- data %>%
        inner_join(final_individuals, by = "IID") %>% 
        inner_join(genetic, by = c("FID", "IID")) %>% 
        mutate(cohort = name)
      
    } else {
      
      #Read in genetic fam file
      genetic <- fread(paste("../", name, "/", name, "_PD.snpqc.fam", sep = ""))
      colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")
      
      data <- data %>% 
        inner_join(genetic, by = c("FID", "IID")) %>% 
        mutate(cohort = name)
      
    }
  
  }
  
  #The age at onset variables are sometimes named age_onset or age_onset_imput
  #Rename to age_onset
  data <- data %>% 
    rename_all(recode, age_onset_imput = "age_onset", age_diagnosis = "age_onset")
  
  #Select just relevant columns
  data_selected <- data %>% 
    select(FID, IID, timeToEvent_death, event_death,
           age_onset, gender, cohort) %>% 
    mutate(gender = tolower(gender)) #make gender values lowercase
  
  #Add to datalist
  datalist[[i]] <- data_selected
 
}


#Merge data in datalist
merged <- bind_rows(datalist)

#Remove PPMI and Campaign as these were excluded from the final GWAS results
merged_mortality <- merged %>% 
  filter(cohort != "Cambridge", cohort !="PPMI") %>% 
  filter(IID != "PACCFCX") #Remove one individual from Oslo cohort which was overlapping with Aasly

merged_mortality %>% 
  group_by(cohort) %>% 
  summarise(count = n())

#Remove rows missing any data
merged_mortality_final <- na.omit(merged_mortality)

#---Mortality summaries---####

#Total number of individuals who died or survived, mean time to event
merged_mortality_final %>% 
  group_by(event_death) %>% 
  summarise(count = n(),
            mean_timeToEvent = mean(timeToEvent_death),
            median_timeToEvent = median(timeToEvent_death))



#---Read in H&Y3 clinical data from each cohort and merge---####

#Make list of cohorts in H&Y3 analysis
cohorts_HY3 <- c("DIGPD", "Oslo", "Oxford", "PROBAND", "PPMI")

#Make empty data list
datalist_HY3 <- list()

#Loop to read in all cohorts, get list of final individuals, and add to data list
for (i in 1:length(cohorts_HY3)){
  
  #Cohort name
  name <- cohorts_HY3[i]
  
  #Directory where clinical data is stored
  directory <- paste("../", name, "/HY3/covars_aao_gender/", sep = "")
  
  #Filename which contains 2020 or 2021
  HY3_filename <- list.files(path = directory,
                                   pattern = "2021-|2020-")
  
  #Read in mortality data
  data <- fread(paste("../", name, "/HY3/covars_aao_gender/",
                      HY3_filename,
                      sep = ""))
  
  #Read in final list of individuals - if it exists
  individuals_filename <- paste("../", name, "/HY3/", name, "_final.txt", sep = "")
    
  #Merge to get final individuals in analysis
  if (file.exists(individuals_filename)) {
      
      #Read in final individuals list
      final_individuals <- fread(individuals_filename)
      
      #Read in genetic fam file
      
      if (name == "PPMI") {
        genetic <- fread(paste("../", name, "/", name, "_july2018.snpqc.PD.fam", sep = ""))
      } else {
        genetic <- fread(paste("../", name, "/", name, "_PD.snpqc.fam", sep = ""))
      }
      
      colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")
      
      data <- data %>%
        inner_join(final_individuals, by = "IID") %>% 
        inner_join(genetic, by = c("FID", "IID")) %>% 
        mutate(cohort = name)
      
    } else {
      
      #Read in genetic fam file
      genetic <- fread(paste("../", name, "/", name, "_PD.snpqc.fam", sep = ""))
      colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")
      
      data <- data %>% 
        inner_join(genetic, by = c("FID", "IID")) %>% 
        mutate(cohort = name)
      
    }
  
  #The age at onset variables are sometimes named age_onset or age_onset_imput
  #Rename to age_onset
  data <- data %>% 
    rename_all(recode, age_onset_imput = "age_onset", age_diagnosis = "age_onset")
  
  #Select just relevant columns
  data_selected <- data %>% 
    select(FID, IID, timeToEvent_HY3, event_HY3,
           age_onset, gender, cohort) %>% 
    mutate(gender = tolower(gender)) #make gender values lowercase
  
  #Add to datalist
  datalist_HY3[[i]] <- data_selected
  
}


#Merge data in datalist
merged_HY3 <- bind_rows(datalist_HY3)

#Remove Oslo as these were excluded from the final GWAS results
merged_HY3_final <- merged_HY3 %>% 
  filter(cohort != "Oslo")

merged_HY3_final %>% 
  group_by(cohort) %>% 
  summarise(count = n())

#---H&Y3 summaries---####

#Total number of individuals who died or survived, mean time to event
merged_HY3_final %>% 
  group_by(event_HY3) %>% 
  summarise(count = n(),
            mean_timeToEvent = mean(timeToEvent_HY3),
            median_timeToEvent = median(timeToEvent_HY3))





#---Read in cognitive impairment clinical data from each cohort and merge---####

#Make list of cohorts in H&Y3 analysis
cohorts_dementia <- c("Aasly", "DIGPD", "Oslo", "Oxford", "PROBAND", "PPMI")

#Make empty data list
datalist_dementia <- list()

#Loop to read in all cohorts, get list of final individuals, and add to data list
for (i in 1:length(cohorts_dementia)){
  
  #Cohort name
  name <- cohorts_dementia[i]
  
  #Directory where clinical data is stored
  directory <- paste("../", name, "/dementia/covars_aao_gender/", sep = "")
  
  #Filename which contains 2020 or 2021
  dementia_filename <- list.files(path = directory,
                             pattern = "2021-|2020-")
  
  #Read in mortality data
  data <- fread(paste("../", name, "/dementia/covars_aao_gender/",
                      dementia_filename,
                      sep = ""))
  
  #Read in final list of individuals - if it exists
  individuals_filename <- paste("../", name, "/dementia/", name, "_final.txt", sep = "")
  
  #Merge to get final individuals in analysis
  if (file.exists(individuals_filename)) {
    
    #Read in final individuals list
    final_individuals <- fread(individuals_filename)
    
    #Read in genetic fam file
    
    if (name == "PPMI") {
      genetic <- fread(paste("../", name, "/", name, "_july2018.snpqc.PD.fam", sep = ""))
    } else {
      genetic <- fread(paste("../", name, "/", name, "_PD.snpqc.fam", sep = ""))
    }
    
    colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")
    
    data <- data %>%
      inner_join(final_individuals, by = "IID") %>% 
      inner_join(genetic, by = c("FID", "IID")) %>% 
      mutate(cohort = name)
    
  } else {
    
    #Read in genetic fam file
    genetic <- fread(paste("../", name, "/", name, "_PD.snpqc.fam", sep = ""))
    colnames(genetic) <- c("FID", "IID", "f", "m", "sex", "pheno")
    
    data <- data %>% 
      inner_join(genetic, by = c("FID", "IID")) %>% 
      mutate(cohort = name)
    
  }
  
  #The age at onset variables are sometimes named age_onset or age_onset_imput
  #Rename to age_onset
  data <- data %>% 
    rename_all(recode, age_onset_imput = "age_onset", age_diagnosis = "age_onset")
  
  #Select just relevant columns
  data_selected <- data %>% 
    select(FID, IID, timeToEvent_dementia, event_dementia,
           age_onset, gender, cohort) %>% 
    mutate(gender = tolower(gender)) #make gender values lowercase
  
  #Add to datalist
  datalist_dementia[[i]] <- data_selected
  
}


#Merge data in datalist
merged_dementia <- bind_rows(datalist_dementia)

#Remove individuals missing any datapoints
merged_dementia_final <- na.omit(merged_dementia)

merged_dementia_final %>% 
  group_by(cohort) %>% 
  summarise(count = n())

#---Dementia summaries---####

#Total number of individuals who met outcome, mean time to event
merged_dementia_final %>% 
  group_by(event_dementia) %>% 
  summarise(count = n(),
            mean_timeToEvent = mean(timeToEvent_dementia),
            median_timeToEvent = median(timeToEvent_dementia))




