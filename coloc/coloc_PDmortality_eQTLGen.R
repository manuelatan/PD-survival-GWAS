# Description: script to run coloc analyses using eQTLGen

#---Import data and libraries-----------------------------------------

#Path for where additional packages are installed (those not already installed on kronos)
pack <- "/data/kronos/kronos/mtan/R_packages"

library(BSgenome)
library(coloc)
library(colochelpR, lib.loc = pack)
library(data.table)
library(MafDb.1Kgenomes.phase3.hs37d5, lib.loc = pack)
library(tidyverse)


mortality_tidy_w_varbeta <- 
  fread("/data/kronos/kronos/mtan/survival_GWAS/coloc/GWAS/PDmortality_tidy_varbeta_noduplicates.txt") %>% 
  as_tibble()

mafdb <- MafDb.1Kgenomes.phase3.hs37d5

GWAS_path <- "/data/kronos/kronos/mtan/survival_GWAS/metaanalysis_mortality/covars_aao_gender/SURVIVAL_MORTALITY_META_20220513_hg191.tbl"
eQTL_path <- "/data/kronos/kronos/mtan/survival_GWAS/coloc/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
results_path <- "/data/kronos/kronos/mtan/survival_GWAS/coloc/results"

#---GWAS--------------------------------------------------------------

# Extract all genes within +/- 1Mb of significant PD hits
ensembl_gene_ids_overlapping_1Mb_window_hit <- 
  mortality_tidy_w_varbeta %>% 
  tidyr::separate(col = SNP, into = c("CHR", "BP"), sep = ":", remove = FALSE) %>% 
  dplyr::mutate(CHR = as.integer(CHR),
                BP = as.integer(BP)) %>% 
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP",
                                                  mart = 37)

# PD sample size 
df1_N <- 5744

# PD proportion cases
df_1_propor_cases <- 1

#---eQTLGen--------------------------------------------------------------

print(str_c(Sys.time(), " - eQTLGen"))

# Load in eQTL file
eQTL <- fread(eQTL_path)

# Create results directory
results_path_GWAS_eQTLGen <- make_results_dir(results_path = results_path, folder_name = "eQTLGen")

# Within results directory, create a folder for "liberal" and "robust" coloc p12 prior
results_path_priors <- setNames(vector(mode = "list", length = 2),
                                c("liberal", "robust"))

results_path_priors$liberal <- make_results_dir(results_path = results_path_GWAS_eQTLGen, folder_name = "liberal")
results_path_priors$robust <- make_results_dir(results_path = results_path_GWAS_eQTLGen, folder_name = "robust")

# Set up loop
for(j in seq_along(ensembl_gene_ids_overlapping_1Mb_window_hit)){
  
  ensembl_gene_id_to_search <- ensembl_gene_ids_overlapping_1Mb_window_hit[j]
  
  print(str_c(Sys.time(), " - PD - eQTLGen - ", ensembl_gene_id_to_search))
  
  # Filter eQTLs for gene, tidy eQTLs
  eQTL_tidy_gene_filtered <- 
    eQTL %>% 
    dplyr::filter(Gene == ensembl_gene_id_to_search) %>%
    colochelpR::tidy_eQTL_eQTLGen() %>% 
    colochelpR::check_coloc_data_format(beta_or_pval = "pval", check_maf = F) %>% 
    dplyr::distinct(SNP, .keep_all = T)
  
  
  if (nrow(eQTL_tidy_gene_filtered) == 0) {
    print(str_c("No QTLs overlapping: ", ensembl_gene_id_to_search))
    next
  }
  
  # Finding mafs and joining
  mafs <- 
    GenomicScores::gscores(x = mafdb, 
                           ranges = eQTL_tidy_gene_filtered$SNP_rs %>% 
                             as.character(), 
                           pop = "EUR_AF")
  mafs <- 
    mafs %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "SNP_rs") %>%
    dplyr::rename(maf = EUR_AF) %>% 
    dplyr::select(SNP_rs, maf) %>% 
    dplyr::mutate(maf = as.numeric(maf)) %>% 
    dplyr::filter(maf > 0 & maf < 1) #Added this line as was getting coloc error
    # dataset 2: MAF should be a numeric, strictly >0 & <1
  
  eQTL_tidy_gene_filtered_maf <- 
    eQTL_tidy_gene_filtered %>%
    dplyr::inner_join(mafs) %>%
    dplyr::select(-SNP_rs) %>% 
    na.omit() %>%
    as_tibble() 
  
  # Run coloc
  p12 <- setNames(c(1e-05, 5e-06),
                  c("liberal", "robust"))
  
  for(k in 1:length(p12)){
    
    print(str_c("Results for '", names(p12[k]), "' p12 prior;  p12 = ", p12[k]))
    
    coloc_results_annotated <-
      colochelpR::get_coloc_results(df1 = mortality_tidy_w_varbeta, df2 = eQTL_tidy_gene_filtered_maf, 
                                    # Harmonise set to false as all it does is flip the beta value 
                                    # As we are using p-values for second dataset this is not necessary
                                    harmonise = F, 
                                    df1_type = "quant", df2_type = "quant", 
                                    df1_beta_or_pval = "beta", df2_beta_or_pval = "pval",
                                    df_1_propor_cases = df_1_propor_cases,
                                    df1_N = df1_N, df2_N = max(eQTL_tidy_gene_filtered_maf$N),
                                    annotate_signif_SNP_df1_df2 = T, 
                                    key_cols = c("GWAS_1", "eQTL_dataset_2", "gene_2"), 
                                    df_1_name = "GWAS", df_2_name = "eQTL", 
                                    df1_path = GWAS_path, df2_path = eQTL_path,
                                    p1 = 1e-04, p2 = 1e-04, p12 = as.numeric(p12[k]))
    
    colochelpR::save_coloc_results(coloc_results_annotated, results_dir_path = results_path_priors[[names(p12[k])]])
    
  }
  
}

print("Done!")