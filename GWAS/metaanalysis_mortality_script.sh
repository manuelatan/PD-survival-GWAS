######META-ANALYSIS SCRIPT######
#Created: 16/03/2020
#Last updated: 13/05/2022
#Created by: Manuela Tan
#WD /data/kronos/kronos/mtan/survival_GWAS/metaanalysis_mortality/covars_aao_gender
#Meta-analysing UKB incident vs. prevalent separately


#Excluding PPMI as not enough people met the outcome
#Including CamPaIGN, separately from Cambridge PD research clinic 
#Excluding Cambridge PD research clinic as genomic inflation factor > 1.2

#Using chr:bp in hg19 (not rsIDs)
#After excluding related individuals from merged dataset

#COHORTS INCLUDED: PROBAND, Oxford, QSBB, UKB incident, UKB prevalent, Calypso, CamPaIGN, DIGPD, Aasly, Oslo

#Genomic control ON

##Download metal software
#http://csg.sph.umich.edu/abecasis/Metal/download/
#Downloaded to /data/kronos/mtan/software
#tar xvzf Linux-metal.tar.gz 

#Run script using:
#/data/kronos/kronos/mtan/software/metal/generic-metal/metal metaanalysis_mortality_script.sh

SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N 

# LOAD INPUT FILES

# Enable Genomic control correction (comment out line as needed)
GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/PROBAND/mortality/covars_aao_gender/PROBAND_survival_mortality_METAL_hg19.tab

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/Oxford/mortality/covars_aao_gender/Oxford_survival_mortality_METAL_hg19.tab


# === DESCRIBE AND PROCESS THE THIRD INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/QSBB/mortality/covars_aao_gender/QSBB_survival_mortality_METAL_hg19.tab


# === DESCRIBE AND PROCESS THE FOURTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/prevalent/UKB_survival_mortality_METAL_hg19.tab


# === DESCRIBE AND PROCESS THE FIFTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/UKB/mortality/covars_aao_gender/incident/UKB_survival_mortality_METAL_hg19.tab

# === DESCRIBE AND PROCESS THE SIXTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/Calypso/mortality/covars_aao_gender/Calypso_survival_mortality_METAL_hg19.tab

# === DESCRIBE AND PROCESS THE SEVENTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/CamPaIGN/mortality/covars_aao_gender/CamPaIGN_survival_mortality_METAL_hg19.tab


# === DESCRIBE AND PROCESS THE EIGTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/DIGPD/mortality/covars_aao_gender/DIGPD_survival_mortality_METAL_hg19.tab


# === DESCRIBE AND PROCESS THE NINTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/Aasly/mortality/covars_aao_gender/Aasly_survival_mortality_METAL_hg19.tab


# === DESCRIBE AND PROCESS THE TENTH INPUT FILE ===
MARKER SNP
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS /data/kronos/kronos/mtan/survival_GWAS/Oslo/mortality/covars_aao_gender/Oslo_survival_mortality_METAL_hg19.tab


OUTFILE SURVIVAL_MORTALITY_META_20220513_hg19 .tbl
ANALYZE HETEROGENEITY

QUIT
