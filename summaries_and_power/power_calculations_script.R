##### Power calculations for survival GWAS #####
#Created by: Manuela Tan

#---Load packages---####
library(survSNP)
library(dplyr)
library(ggplot2)

#Single power calculation - using stats for APOE SNP
#beta for APOE SNP = 0.2945
#lm is the landmark time used for powering the study - this is the median time to the outcome (not censoring)
res_single <- sim.snp.expsurv.power(1.34245, n=5744, raf = 0.16, erate = 0.321,
                                    pilm = 0.5, lm = 10.6, B=0,
                                    model="additive",test="additive",alpha=5e-8)
res_single


#Table with different values for hazard ratios, allele frequency

GHRs<-seq(1.05,1.5,by=0.05)
rafs<-c(0.1, 0.2, 0.3,0.4)
erates<-c(0.3, 0.5,0.7,0.9)

res1<-survSNP.power.table(GHRs, n=5744, rafs, erates,
                          pilm=0.5, lm=10.6, B=0,
                          model="additive",test="additive",alpha=5e-8)

#Note that we are assuming that the median for the survival function in the population is 1 unit of time 
#(that is why we have set pilm=0.5 and lm=1). 
#If we had desired to set power the study based on a population
#whose 0.7 quantile is say 2 units of time, we would have set pilm=0.7 and lm=2.

#Set lm as 10.6 years = median time to death

str(res1)

#Plot
plot <- ggplot(data = res1, mapping  = aes(x = GHR, y = pow0, group = as.factor(raf), color = as.factor(raf))) + 
  geom_line() +
  facet_wrap(~erate) +
  labs(color = "Allele frequency") +
  labs(x= "Hazard Ratio",
       y= "power") +
  theme_bw()

ggsave("plots/power_survivalGWAS.png", width = 9, height = 7)
