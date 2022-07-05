### AD PRS FOR PD PROGRESSION MORTALITY AND DEMENTIA ###
#Created by: Manuela Tan
#Created 19/10/2021
#Last updated 30/06/2022
#Last udpated: Separated Cambridge into Cambridgeclinic and CamPaIGN
#WD: /data/kronos/kronos/mtan/survival_GWAS/AD_PRS
#Goal is to create AD GRS and PRSs
#Excluding APOE region
#Using new AD GWAS Wightman 2021


### Format AD sumstats ###
	
	#See R script
	#12659833 SNPs remaining

### Update SNP names in each PD dataset to the format chr:pos ###

	#Done in PD GRS
	#wd ../candidate_genes/PD_GRS

### PPMI liftover hg38 to hg19 ###
	### PPMI - liftover from hg38 to hg19

	#Make BED file
	R
	library(data.table)
	library(dplyr)
	bim <- fread("PPMI_july2018.snpqc.PD.extra_qc.bim")
	#SNPs should be in chrN:start-end formats
	bim <- bim %>%
		mutate(chr = paste("chr", V1, sep = ""),
				start = V4-1,
				end = V4, 
				name = V2)	

	liftover <- bim %>%
		select(chr, start, end, name)
	#Disable scientific notation
	options(scipen=999)

	#Save file
	write.table(liftover, "PPMI.liftover_hg38_to_hg19.txt", quote = F, row.names = F, col.names = F, sep = "\t")
	q()
	n

	#Run liftOver on mac desktop - because Kronos does not have updated software to run liftover
	./liftOver PPMI.liftover_hg38_to_hg19.txt hg38ToHg19.over.chain.gz output.bed unlifted.bed

	R
	library(data.table)
	library(dplyr)
	output <- fread("output.bed")
	old_map <- fread("PPMI.liftover_hg38_to_hg19.txt")

	#Check whether chromosomes match
	#Merge by rsID
	merged <- output %>%
		inner_join(old_map, by = "V4")

	#Select matching chromosome SNPs only and export
	export_hg19 <- merged %>%
		filter(V1.x == V1.y) %>%
		select(V4, V3.x)

	write.table(export_hg19, "export_hg19.txt", quote = F, col.names = F, row.names = F)


	q()
	n



	#Update positions with hg19
	/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile PPMI_july2018.snpqc.PD.extra_qc \
	--extract export_hg19.txt \
	--update-map export_hg19.txt \
	--make-bed \
	--out PPMI_july2018.snpqc.PD.extra_qc.hg19


	#Recode cohort SNP names to chr:bp
	/data/kronos/kronos/mtan/software/plink2 \
	--bfile /data/kronos/kronos/mtan/survival_GWAS/PPMI/PPMI_july2018.snpqc.PD.extra_qc.hg19 \
	--set-all-var-ids @:# \
	--make-bed \
	--out ../candidate_genes/PD_GRS/PPMI_PD.snpqc.chrbp_ids

	#Remove duplicate SNP IDs otherwise this causes problems when making GRS
	#Need to use plink2 function for this because plink1 only removes duplicates by ID and alleles
	#We want to remove duplicates by ID only
	/data/kronos/kronos/mtan/software/plink2 --bfile ../candidate_genes/PD_GRS/PPMI_PD.snpqc.chrbp_ids \
	--rm-dup force-first \
	--make-bed \
	--out ../candidate_genes/PD_GRS/PPMI_PD.snpqc.chrbp_ids.nodupli
	

### In plink make AD GRS just using GWAS-significant loci ###

	for COHORT in Aasly Calypso CamPaIGN DIGPD Oslo Oxford PROBAND QSBB UKB PPMI
	do
		/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile ../candidate_genes/PD_GRS/"$COHORT"_PD.snpqc.chrbp_ids.nodupli \
		--score Wightman_score.txt 3 5 7 header \
		--out "$COHORT".GRS
	done

### In plink make AD GRS just using GWAS-significant loci excluding APOE ###

	for COHORT in Aasly Calypso CamPaIGN DIGPD Oslo Oxford PROBAND QSBB UKB PPMI
	do
		/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile ../candidate_genes/PD_GRS/"$COHORT"_PD.snpqc.chrbp_ids.nodupli \
		--exclude range APOE_exclude.txt \
		--score Wightman_score.txt 3 5 7 header \
		--out "$COHORT".GRS_noAPOE
	done


#Now see Rproj and script to process and analyse PRSs





