######POST META-ANALYSIS FILTERING######
#Created 17/03/2020
#Last updated: 13/05/2022
#Created by: Manuela Tan
#NOTE THAT THE EFFECT DIRECTION REFERS TO ALLELE1 - which is sometimes the REF/MAJOR ALLELE
#I don't know why METAL does this when I have specified effect vs. non-effect in the input files

R
library(data.table)
library(dplyr)
library(tidyr)

#Read in meta analysis results
data <- fread("SURVIVAL_MORTALITY_META_20220513_hg191.tbl")
#8613171 SNPs

#If using stats from meta-analysis without Genomic Control correction (for LDSC)
#data <- fread("SURVIVAL_MORTALITY_UKBSEP_NOPPMI_META_hg19_GCoff1.tbl")

#Arrange by p-value
data_sorted <- data %>%
	arrange(`P-value`)

#Filter for SNPs that are in at least 1000 individuals
data_filtered_sorted <- data_sorted %>%
	filter(TotalSampleSize >= 1000)
#7696389 SNPs

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
#7313918 SNPs remaining for mortality

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
	
write.table(export_FUMA, "metaanalysis_mortality_FUMA_20220513.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Filter for regional association plots
	
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

#Zip files for LocusZoom
bgzip ./locuszoom/metaanalysis_mortality_LZ_chr19.txt
bgzip ./locuszoom/metaanalysis_mortality_LZ_chr7.txt
bgzip ./locuszoom/metaanalysis_mortality_LZ_chr12.txt



#Tabix
tabix -b2 -e2 -S1 -f ./locuszoom/metaanalysis_mortality_LZ_chr19.txt.gz
tabix -b2 -e2 -S1 -f ./locuszoom/metaanalysis_mortality_LZ_chr7.txt.gz
tabix -b2 -e2 -S1 -f ./locuszoom/metaanalysis_mortality_LZ_chr12.txt.gz

#Format for GCTA-COJO
	
	#SNP names must match those in the reference sample (merged_hg19 all cohorts)

	#Read in merged SNP names
	reference <- fread("/data/kronos/kronos/mtan/survival_GWAS/merged_hg19_all/merged_hg19.snpqc.bim")
	colnames(reference) <- c("chr", "SNP_ref", "x", "bp", "A1_ref", "A2_ref")

	#If the alleles match (including flipped)
	COJO <- data_filtered_sorted_het_MAF %>%
		inner_join(reference, by = c("chr", "bp"))

	#Make new variables for COJO to match the reference sample
	#If the alleles are flipped, then flip the beta effect and frequency
	COJO <- COJO %>%
	 	mutate(A1 = toupper(Allele1),
	 			A2 = toupper(Allele2)) %>%
		mutate(SNP_final = ifelse(A1_ref == A1 & A2_ref == A2, SNP_ref,
								ifelse(A1_ref == A2 & A2_ref == A1, SNP_ref, NA))) %>%
		mutate(A1_final = ifelse(A1_ref == A1 & A2_ref == A2, A1,
								ifelse(A1_ref == A2 & A2_ref == A1, A2, NA))) %>%
		mutate(A2_final = ifelse(A1_ref == A1 & A2_ref == A2, A2,
								ifelse(A1_ref == A2 & A2_ref == A1, A1, NA))) %>%
		mutate(b = ifelse(A1_ref == A1 & A2_ref == A2, Effect,
								ifelse(A1_ref == A2 & A2_ref == A1, -Effect, NA))) %>%
		mutate(freq = ifelse(A1_ref == A1 & A2_ref == A2, Freq1,
								ifelse(A1_ref == A2 & A2_ref == A1, 1 - Freq1, NA)))

	export_COJO <- COJO %>%
		select(SNP_final, A1_final, A2_final, freq, b, StdErr, `P-value`, TotalSampleSize) %>%
		rename(SNP = SNP_final,
			A1 = A1_final,
			A2 = A2_final,
			se = StdErr,
			p = `P-value`,
			N = TotalSampleSize)

	fwrite(export_COJO, "metaanalysis_mortality_COJO.ma", quote = F, row.names = F, col.names = T, sep = "\t")

	#To run COJO - using the merged data as reference sample
	/data/kronos/kronos/mtan/software/gcta_1.92.3beta2/gcta64 \
	--bfile /data/kronos/kronos/mtan/survival_GWAS/merged_hg19_all/merged_hg19.snpqc \
	--maf 0.01 \
	--cojo-file metaanalysis_mortality_COJO.ma \
	--cojo-slct \
	--out COJO_results_slct

	#Because we have only kept SNPs that are in the merged_hg19.snpqc file
	#there are some SNPs in the meta-analysis results that are left out
	#Not sure if should left_join rather than inner_join or whether this would make any difference
	#e.g. if GCTA-COJO only takes the SNPs overlapping between the reference and the meta-analysis results

	#To run COJO - using the merged data as reference sample with larger p-value threshold
	/data/kronos/kronos/mtan/software/gcta_1.92.3beta2/gcta64 \
	--bfile /data/kronos/kronos/mtan/survival_GWAS/merged_hg19_all/merged_hg19.snpqc \
	--maf 0.01 \
	--cojo-file metaanalysis_mortality_COJO.ma \
	--cojo-slct \
	--cojo-p 5e-5 \
	--out COJO_results_slct_5e-5


	
#Annotate with rsIDs for LDSC - should be done with Genomic Control OFF in meta-analysis
		
	#Read in rsIDs file
	rsids <- fread("/data/kronos/kronos/mtan/HRC_rs_ids_GRCh37.txt")

	#Merge results with rsIDs
	data_filtered_sorted_het_MAF_rsids <- data_filtered_sorted_het_MAF %>%
		mutate(chrbp = paste(chr,":",bp, sep = "")) %>%
		inner_join(rsids, by = "chrbp")

	#Convert REF and ALT to lowercase
	data_filtered_sorted_het_MAF_rsids <- data_filtered_sorted_het_MAF_rsids %>%
		mutate(REF = tolower(REF),
				ALT = tolower(ALT))

	#Check allele matches against HRC rsIDs
	data_filtered_sorted_het_MAF_rsids <- data_filtered_sorted_het_MAF_rsids %>%
		mutate(allele_match = ifelse(Allele1 == ALT & Allele2 == REF, "match1",
					ifelse(Allele1 == REF & Allele2 == ALT, "match2", "mismatch")))

	#Count allele matches and mismatches
	data_filtered_sorted_het_MAF_rsids %>%
		group_by(allele_match) %>%
		summarise(count = n())
	#3996934 match1
	#3316984 match2
	#25869 mismatch

	#Filter out allele mismatches
	data_filtered_sorted_het_MAF_rsids_alleleMatch <- data_filtered_sorted_het_MAF_rsids %>%
		filter(allele_match == "match1" | allele_match == "match2")
	#7313918 SNPs

	#Check that rsIDs are unique
	data_filtered_sorted_het_MAF_rsids_alleleMatch[duplicated(data_filtered_sorted_het_MAF_rsids_alleleMatch$ID), ]

	#There are duplicated elements where the rsIDs are missing
	#Remove these
	data_filtered_sorted_het_MAF_rsids_alleleMatch <- data_filtered_sorted_het_MAF_rsids_alleleMatch %>%
		filter(!is.na(ID)) %>%
		filter(ID!=".")
	#7313559 SNPs

	#Also still duplicated rsIDs where there are multiple alleles e.g. 1:96656801:C:A and 1:96656801:C:T
	#Remove these 
	#Have ordered data according to p value so should keep one of the pair with the smaller p value
	data_filtered_sorted_het_MAF_rsids_alleleMatch_unique <- data_filtered_sorted_het_MAF_rsids_alleleMatch %>% distinct(ID, .keep_all = TRUE)
	#7313559 SNPs remaining

	#Export for LDSC
	#Note that munge_sumstats.py interprets A1 as the reference allele and that the A1 column in the .sumstats file format refers to the reference allele.
	#Allele2 is the effect allele
	export_LDSC <- data_filtered_sorted_het_MAF_rsids_alleleMatch_unique %>%
		rename(a1 = Allele1,
				a2 = Allele2,
				snpid = ID,
				N = TotalSampleSize,
				pval = `P-value`,
				beta = Effect,
				se = StdErr) %>%
		select(snpid, a1, a2, N, pval, beta, se)

	fwrite(export_LDSC, "metaanalysis_mortality_LDSC_GCoff.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Format for final results table in paper, annotate nearest gene

	#Annotate with closest genes (just for hits with p < 0.00001)

	#Convert first to df - as otherwise it is data.table and df
	data_filtered_sorted_het_MAF_rsids_alleleMatch_unique_df <- as.data.frame(data_filtered_sorted_het_MAF_rsids_alleleMatch_unique)

	topsnps <- data_filtered_sorted_het_MAF_rsids_alleleMatch_unique_df %>%
		filter(`P-value` < 0.00001)

	#Read in gene coordinates table 
	gene.coords <- read.table("NCBI37.3.gene.loc")

	### Loop over each line in results file

	res.annotated=as.data.frame(do.call(rbind,lapply(1:nrow(topsnps),function(x){
  
 	### Get SNP chr and bp and corresponding genes
  	snp.chr=topsnps[x,1]
	snp.bp=topsnps[x,2]
	gene.coords.chr=gene.coords[which(gene.coords$V2==snp.chr),]
	  
	  ### Calculate distance between snp bp and all genes, select gene with minimum value, combine columns and print distance between SNP and gene
	  ### If distance is 0 then SNP is within gene coordinates. Distance is in BP.
	  b=cbind(topsnps[x,],gene.coords.chr[which.min(abs(snp.bp-((gene.coords.chr$V3+gene.coords.chr$V4)/2))),])
	  if(b$bp<b$V3){
	    d=cbind(b,as.character(b$V3-b$bp))
	  } else if(b$V4<b$bp){
	    d=cbind(b,as.character(b$bp-b$V4))
	  } else if(b$bp>b$V3 & b$bp<b$V4){
	    d=cbind(b,as.character("0"))
	  }
	  names(d)[30]=c("Distance.To.Gene(BP)")
	  d
	})))

	res.annotated=res.annotated[order(res.annotated$`P-value`),]

	res.annotated <- res.annotated %>%
		rename(gene.chr = V2,
				gene.start = V3,
				gene.end = V4,
				gene.name = V6) %>%
		select(-V1)
	
	res.annotated.final <- res.annotated %>%
		select(chr, bp, ID, Allele2, Allele1, gene.name, `Distance.To.Gene(BP)`, Effect, StdErr, `P-value`)

	write.table(res.annotated.final, "final_top10snps_table.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
q()
n


#Zip for FUMA
gzip metaanalysis_mortality_FUMA_20220513.txt




