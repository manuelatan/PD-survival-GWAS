######POST META-ANALYSIS FILTERING######
#Created 17/03/2020
#Last updated: 07/09/2021
#Created by: Manuela Tan
#NOTE THAT THE EFFECT DIRECTION REFERS TO ALLELE1 - which is sometimes the REF/MAJOR ALLELE


R
library(data.table)
library(dplyr)
library(tidyr)

#Read in meta analysis results
data <- fread("SURVIVAL_MORTALITY_META_20210907_hg191.tbl")
#8573053 SNPs

#Arrange by p-value
data_sorted <- data %>%
	arrange(`P-value`)

#Filter for SNPs that are in at least 1000 individuals
data_filtered_sorted <- data_sorted %>%
	filter(TotalSampleSize >= 1000)

#Get the maximum sample size and total SNPs
data_filtered_sorted %>%
	summarise(max_samplesize = max(TotalSampleSize),
	count = n())


#Calculate lambda before applying any other filters
p <- data_filtered_sorted$`P-value`
n <- length(data_filtered_sorted$`P-value`)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda

#Filter out SNPs with HetPVal < 0.05 (Cochran's Q-test for heterogeneity)
#Also filter out SNPs with HetISq > 80
data_filtered_sorted_het <-data_filtered_sorted %>%
	filter(HetPVal > 0.05) %>%
	filter(HetISq < 80)

#Check MAF variability - remove variants with MAF variability > 15%
data_filtered_sorted_het_MAF <- data_filtered_sorted_het %>%
	mutate(MAF_variability = MaxFreq - MinFreq) %>%
	filter(MAF_variability <= 0.15)

#Calculate lambda after applying  filters
p <- data_filtered_sorted_het_MAF$`P-value`
n <- length(data_filtered_sorted_het_MAF$`P-value`)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda

#Separate chr and bp
data_filtered_sorted_het_MAF <- data_filtered_sorted_het_MAF %>%
	separate(MarkerName, into = c("chr", "bp"))

#Make chr and bp  numeric
data_filtered_sorted_het_MAF <- data_filtered_sorted_het_MAF %>%
	mutate(chr = as.numeric(chr),
			bp = as.numeric(bp))

#Disable scientific notation (this doesnt work with the version of fwrite on kronos)
options(scipen=999)

#Export for FUMA
export_FUMA <- data_filtered_sorted_het_MAF %>%
	select(chr, bp, `P-value`, Allele1, Allele2, Effect, StdErr, TotalSampleSize) %>%
	rename(pval = `P-value`)
	
write.table(export_FUMA, "metaanalysis_mortality_FUMA_20210907.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Filter for regional association plots for 3 GWAS significant loci
	
	#Fwrite is making mistakes with scientific notation when writing bp so try write.table

	#Filter for just chr19
	#Alleles must be in capital letters for LocusZoom
	export_LZ_chr19 <- data_filtered_sorted_het_MAF %>%
		arrange(chr, bp) %>%
		select(chr, bp, `P-value`, Allele1, Allele2, Effect, StdErr, TotalSampleSize) %>%
		mutate(Allele1 = toupper(Allele1),
				Allele2 = toupper(Allele2)) %>%
		rename(pval = `P-value`) %>%
		filter(chr == 19)

	write.table(export_LZ_chr19, "./locuszoom/metaanalysis_mortality_LZ_chr19.txt", quote = F, row.names = F, col.names = T, sep = "\t")

	#Filter for just chr7
	export_LZ_chr7 <- data_filtered_sorted_het_MAF %>%
		arrange(chr, bp) %>%
		select(chr, bp, `P-value`, Allele1, Allele2, Effect, StdErr, TotalSampleSize) %>%
		mutate(Allele1 = toupper(Allele1),
				Allele2 = toupper(Allele2)) %>%
		rename(pval = `P-value`) %>%
		filter(chr == 7)

	write.table(export_LZ_chr7, "./locuszoom/metaanalysis_mortality_LZ_chr7.txt", quote = F, row.names = F, col.names = T, sep = "\t")

	#Filter for just chr12
	export_LZ_chr12 <- data_filtered_sorted_het_MAF %>%
		arrange(chr, bp) %>%
		select(chr, bp, `P-value`, Allele1, Allele2, Effect, StdErr, TotalSampleSize) %>%
		mutate(Allele1 = toupper(Allele1),
				Allele2 = toupper(Allele2)) %>%
		rename(pval = `P-value`) %>%
		filter(chr == 12)

	write.table(export_LZ_chr12, "./locuszoom/metaanalysis_mortality_LZ_chr12.txt", quote = F, row.names = F, col.names = T, sep = "\t")

	q()
	n

	#Zip files for LocusZoom
	bgzip metaanalysis_mortality_LZ_chr19.txt
	bgzip metaanalysis_mortality_LZ_chr7.txt
	bgzip metaanalysis_mortality_LZ_chr12.txt


	#Tabix
	tabix -b2 -e2 -S1 -f metaanalysis_mortality_LZ_chr19.txt.gz
	tabix -b2 -e2 -S1 -f metaanalysis_mortality_LZ_chr7.txt.gz
	tabix -b2 -e2 -S1 -f metaanalysis_mortality_LZ_chr12.txt.gz


#Zip for FUMA
gzip metaanalysis_mortality_FUMA_20210907.txt 




