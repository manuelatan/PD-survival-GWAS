### MAKE FINAL TABLE FOR PAPER ###

#Created 10/09/2021
#Last updated 18/05/2022
#Created by Manuela Tan
#WD: /data/kronos/kronos/mtan/survival_GWAS/metaanalysis_mortality/covars_aao_gender




R
library(dplyr)
library(data.table)
library(tidyr)

#Read in FUMA independent significant SNPs
FUMA_ind_SNPs <- fread("./FUMA/FUMA_job184549/IndSigSNPs.txt")
#106 SNPs

#Arrange by p-value
FUMA_ind_SNPs_sorted <- FUMA_ind_SNPs %>%
arrange(p)

#Read in meta-analysis results to get alleles and effect sizes
metal <- fread("SURVIVAL_MORTALITY_META_20220513_hg191.tbl")

metal <- metal %>%
separate(MarkerName, into = c("chr", "pos")) %>%
mutate(chr = as.integer(chr),
pos = as.integer(pos))

FUMA_METAL <- FUMA_ind_SNPs_sorted %>%
inner_join(metal, by = c("chr", "pos"))

#Read in gene coordinates table 
gene.coords <- read.table("NCBI37.3.gene.loc")


#Annotate with closest genes (just for hits with p < 0.00001)

	#Convert first to df - as otherwise it is data.table and df
	FUMA_METAL_df <- as.data.frame(FUMA_METAL)

	### Loop over each line in results file

	res.annotated=as.data.frame(do.call(rbind,lapply(1:nrow(FUMA_METAL_df),function(x){
  
 	### Get SNP chr and bp and corresponding genes
  	snp.chr=FUMA_METAL_df[x,'chr']
	snp.bp=FUMA_METAL_df[x,'pos']
	gene.coords.chr=gene.coords[which(gene.coords$V2==snp.chr),]
	  
	  ### Calculate distance between snp bp and all genes, select gene with minimum value, combine columns and print distance between SNP and gene
	  ### If distance is 0 then SNP is within gene coordinates. Distance is in BP.
	  b=cbind(FUMA_METAL_df[x,],gene.coords.chr[which.min(abs(snp.bp-((gene.coords.chr$V3+gene.coords.chr$V4)/2))),])
	  if(b$pos<b$V3){
	    d=cbind(b,as.character(b$V3-b$pos))
	  } else if(b$V4<b$pos){
	    d=cbind(b,as.character(b$pos-b$V4))
	  } else if(b$pos>b$V3 & b$pos<b$V4){
	    d=cbind(b,as.character("0"))
	  }
	  names(d)[31]=c("distance_to_gene")
	  d
	})))

	#Make final table
	#Sometimes METAL takes the major (more common allele) as the effect
	#If this is the case and Freq1 > 0.5, switch alleles, effect direction and frequency
	#For the 95% CI, round to 2 decimal places and use the sprintf function to keep trailing zeros
	final_table <- res.annotated %>%
		select(chr, pos, rsID, Allele1, Allele2, Freq1, V6, distance_to_gene, Effect, StdErr, p) %>%
		mutate(Allele1 = toupper(Allele1)) %>%
		mutate(Allele2 = toupper(Allele2)) %>%
		mutate(A1 = ifelse(Freq1 > 0.5, Allele2, Allele1),
				A2 = ifelse(Freq1 > 0.5, Allele1, Allele2),
				freq = ifelse(Freq1 > 0.5, 1-Freq1, Freq1),
				beta = ifelse(Freq1 > 0.5, -Effect, Effect)) %>%
		mutate(HR = exp(beta),
				CI_95_lower = sprintf("%.2f", round(exp(beta - 1.96*StdErr),2)),
				CI_95_upper = sprintf("%.2f", round(exp(beta + 1.96*StdErr),2)),
				CI_95 = paste(CI_95_lower, "-", CI_95_upper, sep = "")) %>%
		select(chr, pos, rsID, 
			effect = A1, 
			noneffect = A2, 
			freq, 
			nearest_gene = V6, distance_to_gene, beta, StdErr, 
			HR, CI_95, p)

	#Make final table with top 10 SNPs
	final_table_top10 <- final_table %>%
	head(10)


	write.table(final_table_top10, "final_top10snps_table.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
