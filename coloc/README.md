# Code for colocalization analysis

This repository contains scripts used to perform a colocalisation analysis using the Parkinson's Disease (PD) mortality GWAS and two eQTL datasets:
* [eQTLGen](https://pubmed.ncbi.nlm.nih.gov/34475573/). Downloaded from [https://www.eqtlgen.org/cis-eqtls.html](https://www.eqtlgen.org/cis-eqtls.html)

* [PsychENCODE](https://pubmed.ncbi.nlm.nih.gov/30545857/). Downloaded from [http://resource.psychencode.org/](http://resource.psychencode.org/)


This code has been adapted from Regina Reynolds' very nice scripts for coloc analysis in the RBD GWAS (https://www.medrxiv.org/content/10.1101/2021.09.08.21254232v2) [https://www.medrxiv.org/content/10.1101/2021.09.08.21254232v2]:
* [https://rhreynolds.github.io/RBD-GWAS-analysis/](https://rhreynolds.github.io/RBD-GWAS-analysis/)
* [https://github.com/RHReynolds/RBD-GWAS-analysis](https://github.com/RHReynolds/RBD-GWAS-analysis)

# Code contents

1. [coloc_script.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/coloc/coloc_script.R): This script involves setting up coloc and gathering results
2. [coloc_PDmortality_eQTLGen.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/coloc/coloc_PDmortality_eQTLGen.R): This script runs coloc for PD mortality and eQTLGen
3. [coloc_PDmortality_psychencode.R](https://github.com/huw-morris-lab/PD-survival-GWAS/blob/main/coloc/coloc_PDmortality_psychencode.R): This script runs coloc for PD mortality and PsychENCODE

