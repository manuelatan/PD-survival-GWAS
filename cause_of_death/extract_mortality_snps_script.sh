### EXTRACT TOP MORTALITY SNPS ###

#Top SNPs are from FUMA independent SNPs 
#Data is in hg19/GRCh37

### Extract top 10 SNPs from QSBB and UKB ###

	#QSBB
	/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile ../QSBB/QSBB_PD.snpqc \
	--extract range top10_mortality_snps.txt \
	--recodeA \
	--out QSBB.snpqc.mortality_snps.recodeA


	#UKB
	/data/kronos/kronos/mtan/software/plink_1-9/plink --bfile ../UKB/UKB_PD.snpqc \
	--extract range top10_mortality_snps.txt \
	--recodeA \
	--out UKB.snpqc.mortality_snps.recodeA

