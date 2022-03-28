#Shell script to write a separate r script for each subset of the genotype dataset
#There are 134 subsets of the data (up to 50k SNPs in each)
#So this will create 134 R scripts
for n in $(seq 0 133)
do
	echo "#!/usr/bin/Rscript
	library(data.table)
	library(survival)
	library(dplyr)

	#Load genetic data
	load(\"./subset_$n.RData\")

	#Load clinical data
	clinical <- fread(\"./mortality/covars_aao_gender/COHORT_survival_mortality.txt\")

	#Load Principal Components
	PCs <- fread(\"./PCA.eigenvec\")

	#Select just the first 5 principal components
	PCs <- PCs %>%
		select(V1:V7) %>%
		mutate(IID = paste(V1, \"_\", V2, sep = \"\"),
          FID = IID) %>%
		rename(PC1 = V3,
			PC2 = V4,
			PC3 = V5,
			PC4 = V6,
			PC5 = V7) %>%
   		select(-V1, -V2)
	
	
	#Read in list of final individuals to keep
	#After removing related individuals from merged dataset
	COHORT_final <- fread(\"./mortality/COHORT_final_keep.txt\")

	#Inner join all datasets (only keep individuals with values in all)
	#Individuals missing the clinical data have been removed from the clinical dataset so they should be removed after the join
	TABLE <- clinical %>%
		inner_join(PCs, by = c(\"IID\", \"FID\")) %>%
		inner_join(subset, by = c(\"IID\" = \"FID\")) %>%
		inner_join(COHORT_final, by = c(\"IID\"))


	coefficients<-as.data.frame(matrix(ncol= 9))
 
	names(coefficients) <- c(\"SNP\",\"Coeff\", \"se\", \"Pvalue\", \"Cox.zphPVal\", \"N\", \"ov.lik.ratio\",\"logrank\", \"r2\" )
 
	for (i in 12:ncol(TABLE)) {
	print(colnames(TABLE)[i])
  	cols <- c(i, 1:11)
  	snp <- TABLE[, ..cols]
  	model.cox<- coxph(Surv(snp\$timeToEvent_death, snp\$event_death) ~ snp[[1]] + snp\$age_onset_imput + snp\$gender + snp\$PC1+ snp\$PC2 + snp\$PC3 + snp\$PC4 + snp\$PC5, data=snp)
  	kmz<- cox.zph(model.cox, transform = \"km\")
  	j= i-11
  	coefficients[j,1]<- paste(colnames(TABLE)[i])
  	coefficients[j,2]<- summary(model.cox)\$coefficients[1,1]
  	coefficients[j,3]<- summary(model.cox)\$coefficients[1,3]
  	coefficients[j,4]<- summary(model.cox)\$coefficients[1,5]
  	coefficients[j,5]<- kmz\$table[1,3]
  	coefficients[j,6]<- model.cox\$n
  	coefficients[j,7]<- summary(model.cox)\$logtest[[1]]
  	coefficients[j,8]<- summary(model.cox)\$sctest[[1]]
  	coefficients[j,9]<- summary(model.cox)\$rsq[[1]] # nagelkerke r square
  	}
 
	fwrite(coefficients, \"./mortality/covars_aao_gender/COHORT_survival_mortality_GWASresults_$n.txt\", row.names=FALSE, sep=\"\t\", quote= FALSE)" > GWASscript_$n.r
done