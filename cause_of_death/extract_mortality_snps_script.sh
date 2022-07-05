### EXTRACT TOP MORTALITY SNPS ###
#Created by: Manuela Tan
#Last updated: 23/05/2022
#WD: /data/kronos/kronos/mtan/survival_GWAS/cause_of_death
#Goal: extract top 10 mortality SNPs from each cohort so they can be analysed in R by cause of death

#Top SNPs are from FUMA independent SNPs 
#Data is in hg19/GRCh37

### Extract top 10 SNPs from all cohorts###
	
	#Excluding PPMI and Cambridge PD Research Clinic as they were not included in the final mortality GWAS meta-analysis
	sh

	for COHORT in Aasly Calypso CamPaIGN DIGPD Oslo Oxford PROBAND QSBB UKB

	do
		/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile ../"$COHORT"/"$COHORT"_PD.snpqc \
		--extract range top10_mortality_snps.txt \
		--recodeA \
		--out "$COHORT".snpqc.mortality_snps.recodeA
	done


