########## ANNOTATED MANHATTAN PLOT ###########
#Based on code in https://github.com/GP2-TNC-WG/GP2-Bioinformatics-course/blob/master/Module_III.md#6
#I have modified the script slightly based on the data and also to include all SNPs not just log10P > 3.114074) 



#---Load packages---####
library(tidyverse)
library(data.table)
library(reshape2)
library(ggrepel)
library(svglite)

#---Read in GWAS summary statistics---####

#Read in METAL meta-analysis results
gwas <- fread("../SURVIVAL_dementia_META1.tbl")

#Remove SNPs not passing QC
gwas_qc <- gwas %>%
  filter(TotalSampleSize >= 1000) %>% #Remove SNPs not present in at least 1000 people
  filter(HetPVal > 0.05) %>% #Remove SNPs not passing heterogeneity criteria
  filter(HetISq < 80) %>% #Remove SNPs not passing heterogeneity criteria
  mutate(MAF_variability = MaxFreq - MinFreq) %>%
  filter(MAF_variability <= 0.15) #Remove SNPs with MAF variability > 15%


#Sumstats are in rsID format
#Read in GRCh37/hg19 positions to convert to chr:bp
hg19 <- fread("../../../../HRC_rs_ids_GRCh37.txt")

#Match GWAS sumstats to positions
gwas_qc_pos <- gwas_qc %>% 
  inner_join(hg19, by = c("MarkerName" = "ID")) %>% 
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2),
         allele_match = ifelse(Allele1 == REF & Allele2 == ALT, "match", #check alleles match
                               ifelse(Allele1 == ALT & Allele2 == REF, "match", "mismatch")))
#Note that the Allele1 from METAL results is not always the minor allele (freq > 0.5)

#Check how many allele matches/mismatches
gwas_qc_pos %>% 
  group_by(allele_match) %>% 
  summarise(count = n())

#Remove allele mismatches
#Format GWAS sumstats
gwas_qc_formatted <- gwas_qc_pos %>% 
  filter(allele_match!="mismatch") %>% 
  mutate(SNP = paste("chr", chrbp, sep = "")) %>% 
  separate(chrbp, into = c("CHR", "BP")) %>% 
  select(SNP, CHR, BP, P = `P-value`) %>% 
  mutate(CHR = as.numeric(CHR), #Format CHR and BP as numeric
         BP = as.numeric(BP))



#---Read in gene coordinates and annotate hits in GWAS sumstats---####

#Gene coordinates for hg37 downloaded from MAGMA website
#https://ctg.cncr.nl/software/magma

#Read in gene coordinates
gene.coords <- read.table("../NCBI37.3.gene.loc")


#Select SNPs to annotate
hits <- gwas_qc_formatted %>% 
  filter(P < 5e-08) #Only GWAS significant SNPs will be annotated

#Annotate
### Loop over each line in hits file

res.annotated=as.data.frame(do.call(rbind,lapply(1:nrow(hits),function(x){
  
  ### Get SNP chr and bp and corresponding genes
  snp.chr=hits[[x,2]] #Need to index with double brackets [[]] so that it is stored as a value not dataframe
  snp.bp=hits[[x,3]] #Need to index with double brackets [[]] so that it is stored as a value not dataframe
  gene.coords.chr=gene.coords[which(gene.coords$V2==snp.chr),] #Filter gene locations file to just the chromosome matching the SNP to annotate
  
  ### Calculate distance between snp bp and all genes, select gene with minimum value, combine columns and print distance between SNP and gene
  ### If distance is 0 then SNP is within gene coordinates. Distance is in BP.
  b=cbind(hits[x,],gene.coords.chr[which.min(abs(snp.bp-((gene.coords.chr$V3+gene.coords.chr$V4)/2))),])
  
  #Calculate distance between SNP and gene
  if(b$BP<b$V3){
    d=cbind(b,as.character(b$V3-b$BP))
  } else if(b$V4<b$BP){
    d=cbind(b,as.character(b$BP-b$V4))
  } else if(b$BP>b$V3 & b$BP<b$V4){
    d=cbind(b,as.character("0"))
  }
  names(d)[11]=c("Distance.To.Gene(BP)")
  d
})))

#Select relevant columns from annotated hits dataframe
# STATUS is used for color coding, where 0 and 1 will code in red or orange. This is a chance to indicate confidence of the SNPs if you desire - but can be left as all 0 or all 1

hits_annotated <- res.annotated %>% 
  select(SNP, GENE = V6, P) %>%
  mutate(STATUS = 0) %>% #Going to code all as STATUS = 0, red
  arrange(P) %>% #Arrange by p value
  distinct(GENE, .keep_all = TRUE) %>% #Keep only unique gene names - don't need to annotate twice
  select(SNP, STATUS, GENE) 


#---Munge/tidy GWAS summary statistics for plot---####

# Mung the GWAS summary statistics
# Factor code the p-values as per Pulit et al., 2016
gwas_qc_formatted$log10Praw <- -1*log(gwas_qc_formatted$P, base = 10)
gwas_qc_formatted$log10P <- ifelse(gwas_qc_formatted$log10Praw > 40, 40, gwas_qc_formatted$log10Praw)
gwas_qc_formatted$Plevel <- NA
gwas_qc_formatted$Plevel[gwas_qc_formatted$P < 5E-08] <- "possible"
gwas_qc_formatted$Plevel[gwas_qc_formatted$P < 5E-09] <- "likely"

# Reduction of the GWAS object
# This is for more efficient plotting
# This drops everything not useful in future AUC calcs estiamte in Nalls et al., 2019

#gwasFiltered <- subset(gwas_qc_formatted, log10P > 3.114074)
#I have not filtered the data - want to plot all SNPs


# Highlight the hits of interest to annotate
snpsOfInterest <- hits_annotated$SNP

# Merge the filtered GWAS with the annotated hits of interest
gwasToPlotUnsorted <- gwas_qc_formatted %>% 
  left_join(hits_annotated, by = "SNP") #modified original code slightly as it wsas not keeping GENE column

gwasToPlot <- gwasToPlotUnsorted[order(gwasToPlotUnsorted$CHR,gwasToPlotUnsorted$BP),]


#---Prepare the daraset to plot---####

# Prepare the dataset to plot
plotting <- gwasToPlot %>%
  group_by(CHR) %>% # Space out the chromosomes accordingly
  summarize(chr_len=max(BP)) %>%
  
  mutate(tot=cumsum(chr_len)-chr_len) %>% # Calculate the cumulative position of each chromosome
  select(-chr_len) %>%
  
  left_join(gwasToPlot, ., by=c("CHR"="CHR")) %>% # Have this information added to the original dataset so you can subset it for plotting
  
  arrange(ordered(CHR), BP) %>%
  mutate(BPcum=BP+tot) %>% # Space out the SNPs accordingly
  
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>% # Highlight the hits
  mutate(is_annotate=ifelse(log10P>7.30103, "yes", "no"))

#---Plot---####

# Have the x-axis accomodate all chromosome sizes
axisdf <- plotting %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot panel
thisManhattan <- ggplot(plotting, aes(x=BPcum, y=log10P)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) + # Show all the points and color depending on chromosome
  scale_color_manual(values = rep(c("light grey", "dark grey"), 22 )) +
  
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) + # Custom X axis that removes the spaces between X axis and SNPs
  scale_y_continuous(expand = c(0, 1) ) + #increased expand, this increases the Y axis scale around the highest point otherwise it was getting cut off
  geom_point(data=subset(plotting, is_highlight=="yes" & Plevel == "likely"), color = "red", size = 2) + # add highlighted points, these highlighted points are adding an estimate of confidence for "genome-wide significant hits"
  geom_point(data=subset(plotting, is_highlight=="yes" & Plevel == "possible"), color = "orange", size = 2) + # red = more likely to replicate than orange -- related to Pulit et al. 2016
  
  geom_label_repel(data=subset(plotting, is_annotate=="yes"), aes(label=GENE, fill=factor(STATUS)), alpha = 0.5,  size=2) + # # add label using ggrepel to avoid overlapping, here the label color coding is STATUS from the hits file and the text is the GENE from that file
  scale_fill_manual(values = c("aquamarine","cyan")) +
  
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  xlab("BP") +
  ylab("-log10P") +
  geom_hline(yintercept = -log10(5e-08), #Add horizontal line at p = 5e-08 to indicate GWAS significance
             linetype = "dashed",
             color = "blue")




# Export it
ggsave("./plots/dementia_manhattan.pdf", thisManhattan, width = 8, height = 4, dpi=300, units = "in")
ggsave("./plots/dementia_manhattan.png", thisManhattan, width = 7, height = 3.5, dpi=300, units = "in")
ggsave("./plots/dementia_manhattan.svg", thisManhattan, width = 7, height = 3.5, dpi=300, units = "in")

# Warning message for rows missing values for geom_label_repel
#These are the SNPs passing GWAS significance but we are not annotating with gene names
#Because another SNP in the locus with the same gene has already been annotated
#This is fine - the SNPs are still plotted just not annotated with gene names