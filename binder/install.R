## General dependencies:
install.packages("seqinr")
install.packages("stringr")
install.packages("devtools")

## BioC dependencies

#install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "sangerseqR"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "sangerseqR"))

devtools::install_github("ipb-halle/GoldenMutagenesis", ref="binder")
