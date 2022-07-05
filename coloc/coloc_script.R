##### COLOCALISATION ANALYSIS #####

#Following Regina Reynolds script for coloc analysis of RBD GWAS
#https://rhreynolds.github.io/RBD-GWAS-analysis/
#https://github.com/RHReynolds/RBD-GWAS-analysis

#---Load packages---####

library(remotes)
library(tidyverse)
library(coloc)
library(data.table)
library(colochelpR) #coloc helper functions package 
library(tidyr)
library(MafDb.1Kgenomes.phase3.hs37d5)
library(here)

#coloc package installed from
#install_github("chr1swallace/coloc",build_vignettes=TRUE)
#https://chr1swallace.github.io/coloc/index.html

#colochelpR
#install_github("RHReynolds/colochelpR")
#https://github.com/RHReynolds/colochelpR

#if there are problems installing colochelpR because of other packages
#Install them from BiocManager
#BiocManager::install("GenomicRanges")


#---Read in meta-analysis results---####

#Read in mortality results
mortality <- fread("../metaanalysis_mortality/covars_aao_gender/SURVIVAL_MORTALITY_META_20220513_hg191.tbl")

#Filter for SNPs that pass QC
mortality_qc <- mortality %>% 
  filter(TotalSampleSize >= 1000) %>% #Filter for SNPs that are in at least 1000 individuals
  filter(HetPVal > 0.05) %>% #Filter for SNPs that pass heterozygosity criteria
  filter(HetISq < 80) %>% 
  mutate(MAF_variability = MaxFreq - MinFreq) %>%
  filter(MAF_variability <= 0.15) %>% #Remove SNPs that have >0.15 MAF variability
  mutate(Al1 = toupper(Allele1),
         Al2 = toupper(Allele2))

#Tidy GWAS data and rename columns
mortality_tidy <- mortality_qc %>%
  dplyr::rename(SNP = MarkerName,
         beta = Effect,
         se = StdErr,
         p.value = `P-value`,
         maf = Freq1) %>% 
  dplyr::mutate(GWAS = "PD_mortality") %>% 
  dplyr::select(GWAS, SNP, beta, se, p.value, Al1, Al2, maf)

#---Format GWAS data for coloc---####

#As some MAFs are > 0.5, need to switch alleles and inverse beta
mortality_tidy_new <- mortality_tidy %>% 
  dplyr::mutate(beta_new = ifelse(maf > 0.5, -beta, beta),
         Al1_new = ifelse(maf > 0.5, Al2, Al1),
         Al2_new = ifelse(maf > 0.5, Al1, Al2),
         maf_new = ifelse(maf > 0.5, 1 - maf, maf)) %>% 
  dplyr::select(GWAS, SNP, beta_new, se, p.value, Al1_new, Al2_new, maf_new) %>% 
  dplyr::rename(beta = beta_new,
         Al1 = Al1_new,
         Al2 = Al2_new,
         maf = maf_new)

#Run helper functions to check data and calculate varbeta
mortality_tidy_w_varbeta <- mortality_tidy_new %>%
  get_varbeta() %>% 
  check_coloc_data_format(beta_or_pval = "beta", check_maf = T)

#Save as output
fwrite(mortality_tidy_w_varbeta, "./GWAS/PDmortality_tidy_varbeta.txt",
       quote = F, row.names = F, col.names = T, sep = "\t")


#Sort out biallelic SNPs (SNPs that have the same chr:bp but multiple alleles)
#We will use the SNP with the highest MAF in the data

#Find duplicates and choose by highest MAF
duplicates_top_maf <- mortality_tidy_w_varbeta %>% 
  dplyr::filter(!(duplicated(SNP) | duplicated(SNP, fromLast = TRUE))) %>% 
  dplyr::group_by(SNP) %>% 
  dplyr::top_n(1, maf)

# Remove any remaining duplicates which do not differentiate on beta/p-value and maf
duplicates_top_maf <- duplicates_top_maf %>% 
  dplyr::anti_join(duplicates_top_maf %>% 
                     dplyr::filter(duplicated(SNP) | duplicated(SNP, fromLast = TRUE)))

# Above steps results in removal of all duplicate SNPs 
# That is all SNP duplicates are duplicates (and not multiallelic SNPs)
nrow(duplicates_top_maf)

# Thus simply use dplyr::distinct() to find all unique entries
mortality_tidy_w_varbeta <- mortality_tidy_w_varbeta %>% 
  dplyr::distinct(across(everything()))

#Remove any SNPs with beta = 0, as this produces NaNs in coloc output otherwise
mortality_tidy_w_varbeta <- mortality_tidy_w_varbeta %>% 
  filter(beta!=0)

#Save as output
fwrite(mortality_tidy_w_varbeta, "./GWAS/PDmortality_tidy_varbeta_noduplicates.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")

#---Extract all genes within +/- 1 Mb of all significant hits in the GWAS---####
#Using GRCh37 gene definitions from ensembl
#GWAS data is in GRCh37 build

#First just have a look at the GWAS significant hits
mortality_tidy_w_varbeta %>% 
  filter(p.value < 5e-8)

# Extract all genes within +/- 1Mb of significant GWAS hits
ensembl_gene_ids_overlapping_1Mb_window_hit <- mortality_tidy_w_varbeta %>% 
  tidyr::separate(col = SNP, into = c("CHR", "BP"), sep = ":", remove = FALSE) %>% 
  dplyr::mutate(CHR = as.integer(CHR),
                BP = as.integer(BP)) %>% 
  get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                      CHR_column = "CHR",
                                      BP_column = "BP",
                                      mart = 37)


# Number of genes within +/- 1Mb of all significant hits in PD progression GWAS
ensembl_gene_ids_overlapping_1Mb_window_hit <-
  tibble(ensembl_id = ensembl_gene_ids_overlapping_1Mb_window_hit) %>% 
  colochelpR::biomart_df(.,
                         columnToFilter = "ensembl_id", 
                         attributes = c("ensembl_gene_id", "hgnc_symbol",
                                        "chromosome_name", "start_position", "end_position"),
                         filter = "ensembl_gene_id",
                         mart = 37)


#Save as RDS
saveRDS(ensembl_gene_ids_overlapping_1Mb_window_hit, "ensembl_gene_ids_overlapping_1Mb_window_hit.RDS")

#---eQTL datasets used---####

##eQTLGen
#31,684 blood samples from 37 cohorts (assayed using 3 gene expression platforms)
#Data downloaded from: https://www.eqtlgen.org/cis-eqtls.html. 
#Full cis-eQTL summary statistics
#Downloaded on 22/03/2022

##Psychencode
#eQTL data derived from 1387 individuals (number is retrieved from FAQ section of psychencode resource: https://faq.gersteinlab.org/category/capstone4/page/1/; see post entitled: "Inquiry regarding PsychENCODE eQTL resource" from May 3, 2019).
#Data downloaded from: http://resource.psychencode.org/. 
#Full set of cis-eQTLs with no p-value or FDR filtering: Full_hg19_cis-eQTL.txt.gz
#Downloaded on 22/03/2022.


#---Setting up coloc: Look up MAFs---####

#Look up MAFS
#MAFs were not available in both eQTL datasets so need to use p-value and MAF instead
#MAFs not available in files so used reference database to look up MAFs
#This can be performed using the package MafDb.1Kgenomes.phase3.hs37d5, which contains MAFs for a number of populations.

mafdb <- MafDb.1Kgenomes.phase3.hs37d5
populations(mafdb)


#Used EUR_AF for all lookups

example_df <- data.frame(rs_id = c("rs12921634", "rs1476958", "rs56189750"))

# Example code
mafs <- GenomicScores::gscores(x = mafdb, ranges = example_df$rs_id %>% as.character(), pop = "EUR_AF")

mafs <- mafs %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "SNP") %>%
  dplyr::rename(maf = EUR_AF) %>% 
  dplyr::select(SNP, maf) %>% 
  dplyr::mutate(maf = as.numeric(maf)) %>% 
  dplyr::filter(maf > 0 & maf < 1) 

example_df <- example_df %>%
  inner_join(mafs, by = c("rs_id" = "SNP"))

#---Run coloc analysis with eQTLGen---####

#Using script: coloc_PDmortality_eQTLGen.R
#WD: /data/kronos/kronos/mtan/survival_GWAS/coloc/results
#Run using command:
#nohup Rscript /data/kronos/kronos/mtan/survival_GWAS/coloc/coloc_PDmortality_eQTLGen.R &>/data/kronos/kronos/mtan/survival_GWAS/coloc/logs/coloc_PDmortality_eQTLGen.log&


#---Run coloc analysis with Psychencode---####

#Using script: coloc_PDmortality_psychencode.R
#WD: /data/kronos/kronos/mtan/survival_GWAS/coloc/results
#Run using command:
#nohup Rscript /data/kronos/kronos/mtan/survival_GWAS/coloc/coloc_PDmortality_psychencode.R &>/data/kronos/kronos/mtan/survival_GWAS/coloc/logs/coloc_PDmortality_psychencode.log&




#---Load results---####

results_dir <- tibble(dir = list.files(here::here("results"), recursive = T, full.names = T, pattern = ".rda") %>% 
           str_replace(., "/[^/]*$", "") %>%
           str_replace(., "/[^/]*$", ""),
file_path = list.files(here::here("results"), recursive = T, full.names = T, pattern = ".rda"),
         file_name = list.files(here::here("results"), recursive = T, full.names = F, pattern = ".rda") %>% 
           str_replace(., ".rda", "")) %>%
  tidyr::separate(file_name, into = c("dataset", "prior_type", "gene"), sep = "/", remove = F)

# Split into liberal/robust results
dataset_dir <- setNames(results_dir %>% dplyr::group_split(prior_type),
                        c(results_dir %>% .[["prior_type"]] %>% unique() %>% sort()))

# Priors df
priors_df <- tibble(prior_type = c("liberal", "robust"),
                    p12 = c(1e-05, 5e-06))

all_results <- setNames(vector(mode = "list", length = length(dataset_dir)),
                        names(dataset_dir))

for(i in 1:length(dataset_dir)){
  
  priors_df_filtered <- 
    priors_df %>% 
    dplyr::filter(prior_type == names(dataset_dir[i]))
  
  dataset_file_paths <- dataset_dir[[i]]
  
  dataset_names <- dataset_file_paths$dataset %>% unique()
  
  results_list <- vector(mode = "list", length = length(dataset_names))
  
  for(j in 1:length(dataset_names)){
    
    dataset <- dataset_names[j]
    
    print(dataset)
    
    dir_to_load <- 
      dataset_file_paths %>% 
      dplyr::filter(dataset == dataset_names[j]) %>% 
      .[["dir"]] %>% 
      unique()
    
    dir_to_load <- str_c(dir_to_load, "/", priors_df_filtered$prior_type)
    
    for(k in 1:length(dir_to_load)){
      
      print(dir_to_load[k])
      
      if(dataset %in% c("eQTLGen", "psychencode")){
        
        results <- 
          colochelpR::merge_coloc_summaries(dir_to_load[k], add_signif_SNP = F, recursive = T, pattern = ".rda") %>% 
          dplyr::select(GWAS_1, gene_2, everything(), -eQTL_dataset_2)
        
      }
      
      results <- 
        results %>% 
        dplyr::mutate(dataset = dataset,
                      p12 = priors_df_filtered$p12) %>% 
        dplyr::select(GWAS_1, dataset, gene_2, everything())
      
      if(k == 1){
        
        results_list[[j]] <- results 
        
      } else{
        
        results_list[[j]] <- 
          results_list[[j]] %>% 
          dplyr::bind_rows(results)
        
      }
      
    }
    
  }
  
  all_results[[i]] <- results_list
  
}

results <- 
  all_results %>% 
  lapply(., function(x){
    
    x %>% 
      qdapTools::list_df2df() %>% 
      biomart_df(columnToFilter = "gene_2",
                 mart = 37,
                 attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = c("ensembl_gene_id")) %>% 
      dplyr::select(GWAS_1, dataset, gene_2, hgnc_symbol, everything(), -X1)
    
  }) 

saveRDS(results, file = here::here("results", "coloc_summary.Rds"))


#---Filtering by PP.H4 > 0.75 (i.e. colocalisation)---####
#Typically, a PP.H4 > 0.75 is considered to be strong evidence of a shared SNP.

#Read in results file if not loaded
results <- readRDS(here::here("results", "coloc_summary.Rds"))

coloc <- 
  results %>% 
  qdapTools::list_df2df(., col1 = "prior_type") %>% 
  dplyr::mutate(hgnc_symbol = case_when(gene_2 == "ENSG00000247775" ~ "SNCA-AS1",
                                        TRUE ~ hgnc_symbol)) %>%
  dplyr::inner_join(results %>% 
                      qdapTools::list_df2df(., col1 = "prior_type") %>% 
                      dplyr::filter(prior_type == "liberal", PP.H4.abf >= 0.75) %>%
                      dplyr::distinct(dataset, gene_2)) 

coloc %>% 
  dplyr::arrange(dataset, prior_type, hgnc_symbol) %>% 
  dplyr::select(prior_type, p12, GWAS_1, dataset, gene_2, hgnc_symbol, everything()) %>% 
  DT::datatable(rownames = FALSE,
                options = list(scrollX = TRUE),
                class = 'white-space: nowrap') 


#Make results dataframe with no filtering of PP.H4.abf - just to check max value
full_results <- results %>% 
  qdapTools::list_df2df(., col1 = "prior_type") %>% 
  arrange(-PP.H4.abf) #Arrange by PP.H4

#Look at genes with sum PP.H3 + PP.H4 > 0.8
#Cumulative probability of 80% that SNps associate with both trait 1 and trait 2 within the region in question
full_results %>% 
  dplyr::filter(sum_PPH3_PPH4 > 0.8) %>% 
  dplyr::arrange(-PP.H4.abf) %>% #Arrange by PP.H4.abf largest to smallest
  DT::datatable(rownames = FALSE,
                options = list(scrollX = TRUE),
                class = 'white-space: nowrap')
#The highest PP.H4 is 0.028
#No evidence that SNPs are shared

#Save full results
write.table(full_results, "./outputs/coloc_full_results.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")


#----Conclusion---####

#Conclude: There are no genes with PP.H4 > 0.75
#No evidence of colocalisation
#PP.H3 appears to be high, indicating they appear to be independent SNPs
