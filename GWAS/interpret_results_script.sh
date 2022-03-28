######SCRIPT TO COMBINE GWAS RESULTS AND FORMAT FOR FUMA AND META-ANALYSIS######
#Created: 24/02/2020
#Created by Manuela Tan
#Last updated: 08/06/2021
#Including an export for METAL using chr:bp not rsID
#This is because we have excluded PPMI from the meta-analysis of mortality and all the other datasets are in GRCh37/hg19 so we can merge by chr:bp, no need to merge with rsIDs, may be losing some SNPs

#Make combined file
cat COHORT_survival_mortality_GWASresults_* | grep -v 'SNP' > allGWAS_results.txt

#Read into R
R
library(dplyr)
library(data.table)
library(tidyr)

#Read in data
data <- fread("allGWAS_results.txt")

#Label columns
colnames(data) <- c("SNP", "Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio", "logrank", "r2")


#Sort by p value
data <- data %>% arrange(Pvalue)

#Check Cox.zphPVal
#This is really, really important to get right. If this p-value is very small, then it means that the variable is time dependent and therefore needs to be treated separately. 
#Check this column if you get any significant hits.
 
#Calculate lambda
p <- data$Pvalue
n <- length(data$Pvalue)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda

#Remove X from the SNP names
data <- data %>%
	mutate(SNP = substring(SNP, 2))

#Split SNP name into chr and position and alleles
data_split <- data %>%
	separate(SNP, into = c("chr", "bp", "REF", "ALT", "A1"))

#Check some instances where the A1 coded by plink --recode is not the same as the ALT column (from VCF)
#We need to use the plink A1 column because this is the allele count that is used for the survival model
#Everything else is the non-effect allele
data_split  <- data_split %>%
  mutate(effect = ifelse(A1 == ALT, ALT,
                          ifelse(A1 == REF, REF, NA)),
         noneffect = ifelse(A1 == ALT, REF,
                          ifelse(A1 == REF, ALT, NA)))

#Remove indels
data_split_noindels <- data_split %>%
	filter(!is.na(effect)) %>%
	filter(!is.na(noneffect))


#Export for FUMA
fwrite(data_split_noindels, "COHORT_survival_mortality_FUMA.txt", quote = F, row.names = F, col.names = T, sep = "\t")



#Read in frequency file - needed for METAL
freq <- fread("COHORT_PD.snpqc.freq.frq")

freq <- freq %>%
	select(-CHR, -NCHROBS) %>%
	rename(A1_freq = A1,
		A2_freq = A2)

#Merge with results file
data_split_noindels_freq <- data_split_noindels %>%
	mutate(SNP = paste(chr,bp,REF,ALT, sep = ":")) %>%
	inner_join(freq, by = "SNP")

#Check that alleles match (there should be no mismatches)
data_split_noindels_freq %>%
	filter(A1!=A1_freq) %>%
	summarise(count = n())

#Export in GRCh37/hg19 with chr:pos for METAL
export_METAL_hg19 <- data_split_noindels_freq %>%
	mutate(SNP_new = paste(chr, bp, sep = ":")) %>%
	select(SNP_new, effect, noneffect, Coeff, se, Pvalue, N, MAF) %>%
	rename(SNP = SNP_new,
	effect_allele = effect,
		noneffect_allele = noneffect,
		beta = Coeff)

#Export for METAL
fwrite(export_METAL_hg19, "COHORT_survival_mortality_METAL_hg19.tab", quote = F, sep = "\t", row.names = F)

#Read in rsIDs file
rsids <- fread("HRC_rs_ids_GRCh37.txt")

#Merge results with rsIDs
data_split_noindels_freq_rsids <- data_split_noindels_freq %>%
	mutate(chrbp = paste(chr,bp, sep =":")) %>%
	inner_join(rsids, by = "chrbp")

#Check allele matches against HRC rsIDs
data_split_noindels_freq_rsids <- data_split_noindels_freq_rsids %>%
	mutate(allele_match = ifelse(effect == ALT.y & noneffect == REF.y, "match1",
				ifelse(effect == REF.y & noneffect == ALT.y, "match2", "mismatch")))

#Count allele matches and mismatches
data_split_noindels_freq_rsids %>%
	group_by(allele_match) %>%
	summarise(count = n())

#Filter out allele mismatches
data_split_noindels_freq_rsids_alleleMatch <- data_split_noindels_freq_rsids %>%
	filter(allele_match == "match1" | allele_match == "match2")

#Check that rsIDs are unique
data_split_noindels_freq_rsids_alleleMatch[duplicated(data_split_noindels_freq_rsids_alleleMatch$ID), ]

#There are duplicated elements where the rsIDs are missing
#Remove these
data_split_noindels_freq_rsids_alleleMatch <- data_split_noindels_freq_rsids_alleleMatch %>%
	filter(!is.na(ID)) %>%
	filter(ID!=".")

#Also still duplicated rsIDs where there are multiple alleles e.g. 1:96656801:C:A and 1:96656801:C:T
#Remove these 
#Have ordered data according to p value so should keep one of the pair with the smaller p value
data_split_noindels_freq_rsids_alleleMatch_unique <- data_split_noindels_freq_rsids_alleleMatch %>% 
	distinct(ID, .keep_all = TRUE)

#Format for METAL
export_METAL <- data_split_noindels_freq_rsids_alleleMatch_unique %>%
	select(ID, effect, noneffect, Coeff, se, Pvalue, N, MAF) %>%
	rename(rsid = ID,
		effect_allele = effect,
		noneffect_allele = noneffect,
		beta = Coeff)

#Export for METAL with rsIDs
fwrite(export_METAL, "COHORT_survival_mortality_METAL.tab", quote = F, sep = "\t", row.names = F)

q()
n