##### PD GRS VS PROGRESSION #####
#Goal: created PD GRS in each cohort and analyse in relation to mortality, H&Y3, dementia
#Created by: Manuela Tan
#Created: 15/10/2021
#Last updated: 15/10/2021
#WD: /data/kronos/kronos/mtan/survival_GWAS/candidate_genes/PD_GRS


### Get allele weights from Nalls 2019 GWAS ###
	
	#See candidate_genes_script.R


### Make GRS ###
	
	#For all cohorts except PPMI and UKB

	sh
	for COHORT in Aasly Calypso Cambridge DIGPD Oslo Oxford PROBAND QSBB

	do

		#Recode cohort SNP names to chr:bp
		/data/kronos/kronos/mtan/software/plink2 \
		--bfile /data/kronos/kronos/mtan/survival_GWAS/$COHORT/"$COHORT"_PD.snpqc \
		--set-all-var-ids @:# \
		--make-bed \
		--out "$COHORT"_PD.snpqc.chrbp_ids

		#Remove duplicate SNP IDs otherwise this causes problems when making GRS
		#Need to use plink2 function for this because plink1 only removes duplicates by ID and alleles
		#We want to remove duplicates by ID only
		/data/kronos/kronos/mtan/software/plink2 --bfile "$COHORT"_PD.snpqc.chrbp_ids \
		--rm-dup force-first \
		--make-bed \
		--out "$COHORT"_PD.snpqc.chrbp_ids.nodupli

		#Make GRS
		#Can also use the header option as the score file has headers
		/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile "$COHORT"_PD.snpqc.chrbp_ids.nodupli \
		--score Nalls_score_chrbp.txt \
		--out "$COHORT".GRS

	done


	#For UKB - generate GRSs on whole cohort as it should be independent. Will split later 

		#Recode cohort SNP names to chr:bp
		/data/kronos/kronos/mtan/software/plink2 \
		--bfile /data/kronos/kronos/mtan/survival_GWAS/UKB/UKB_PD.snpqc \
		--set-all-var-ids @:# \
		--make-bed \
		--out UKB_PD.snpqc.chrbp_ids

		#Remove duplicate SNP IDs otherwise this causes problems when making GRS
		#Need to use plink2 function for this because plink1 only removes duplicates by ID and alleles
		#We want to remove duplicates by ID only
		/data/kronos/kronos/mtan/software/plink2 --bfile UKB_PD.snpqc.chrbp_ids \
		--rm-dup force-first \
		--make-bed \
		--out UKB_PD.snpqc.chrbp_ids.nodupli

		#Make GRS
		/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile UKB_PD.snpqc.chrbp_ids.nodupli \
		--score Nalls_score_chrbp.txt header\
		--out UKB.GRS

	done

	#For PPMI - using rsIDs

		#Make GRS
		/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile /data/kronos/kronos/mtan/survival_GWAS/PPMI/PPMI_july2018.snpqc.PD.extra_qc \
		--score Nalls_score_rsid.txt \
		--out PPMI.GRS


