# Code for colocalization analysis

This repository contains scripts used to perform a colocalisation analysis using the Parkinson's Disease (PD) mortality GWAS and two eQTL datasets:
* [eQTLGen](https://pubmed.ncbi.nlm.nih.gov/34475573/). Downloaded from [https://www.eqtlgen.org/cis-eqtls.html](https://www.eqtlgen.org/cis-eqtls.html)

* [PsychENCODE](https://pubmed.ncbi.nlm.nih.gov/30545857/). Downloaded from [http://resource.psychencode.org/](http://resource.psychencode.org/)


# Code contents

1. [coloc_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/5df101b6816f0fde138779845f61c951d1cb0ac4/coloc/coloc_script.R): This script involves setting up coloc and gathering results
2. [coloc_PDmortality_eQTLGen.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/2f93271100267345859577c96367988b54b68b75/coloc/coloc_PDmortality_eQTLGen.R): This script runs coloc for PD mortality and eQTLGen
3. [coloc_PDmortality_psychencode.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/2f93271100267345859577c96367988b54b68b75/coloc/coloc_PDmortality_psychencode.R): This script runs coloc for PD mortality and PsychENCODE

