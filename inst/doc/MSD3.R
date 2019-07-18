## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

## ------------------------------------------------------------------------
library("GoldenMutagenesis")

## ------------------------------------------------------------------------
input_sequence<-"ATGTCTCAGGTTCAGAGTGGCATTTTGCCAGAACATTGCCGCGCGGCGATTTGGATCGAAGCCAACGTGAAAGGGGAAGTTGACGCCCTGCGTGCGGCCAGTAAAACATTTGCCGACAAACTGGCAACTTTTGAAGCGAAATTCCCGGACGCGCATCTTGGTGCGGTGGTTGCCTTTGGTAACAACACCTGGCGCGCTCTGAGCGGCGGCGTTGGGGCAGAAGAGCTGAAAGATTTTCCGGGCTACGGTAAAGGCCTTGCGCCGACGACCCAGTTCGATGTGTTGATCCACATTCTTTCTCTGCGTCACGACGTAAACTTCTCTGTCGCCCAGGCGGCGATGGAAGCCTTTGGTGACTGCATTGAAGTGAAAGAAGAGATCCACGGCTTCCGTTGGGTTGAAGAGCGTGACCTGAGCGGCTTTGTTGACGGTACGGAAAACCCGGCGGGTGAAGAGACGCGTCGCGAAGTGGCGGTTATCAAAGACGGCGTGGATGCGGGCGGCAGCTATGTGTTTGTCCAGCGTTGGGAACACAACCTGAAGCAGCTCAACCGGATGAGCGTTCACGATCAGGAGATGGTGATCGGGCGCACCAAAGAGGCCAACGAAGAGATCGACGGCGACGAACGTCCGGAAACCTCTCACCTCACCCGCGTTGATCTGAAAGAAGATGGCAAAGGGCTGAAGATTGTTCGCCAGAGCCTGCCGTACGGCACTGCCAGTGGCACTCACGGTCTGTACTTCTGCGCCTACTGCGCGCGTCTGCATAACATTGAGCAGCAACTGCTGAGCATGTTTGGCGATACCGATGGTAAGCGTGATGCGATGTTGCGTTTCACCAAACCGGTAACCGGCGGCTATTATTTCGCACCGTCGCTGGACAAGTTGATGGCGCTGTAA"
recognition_site_bbsi<-"GAAGAC"
recognition_site_bsai<-"GGTCTC"
cuf<-"e_coli_316407.csv"



## ------------------------------------------------------------------------
mutations_bbsi<-domesticate(input_sequence, recognition_site_bbsi, cuf=cuf)
mutations_bbsi
mutations_bsai<-domesticate(input_sequence, recognition_site_bsai, cuf=cuf)
mutations_bsai


## ------------------------------------------------------------------------
#If domestication is necessary follow the workflow of the Point Mutagenesis vignette
#We will use a structure similar to the Point Mutagenesis workflow to insert different kinds of Multiple Site Saturations here.
mutations<-list(c(137,"NDT"), c(143,"NNN"), c( 147, "DBK"),c(232, "NDT"), c(234, "NDT"))
primers<-msd_mutate(input_sequence, prefix="TT" ,restriction_enzyme=recognition_site_bsai, suffix="A", vector=c("AATG", "AAGC"), replacements=mutations, replacement_range=5, binding_min_length=4 ,primer_length=9, target_temp=60, fragment_min_size=60 )
primers

## ------------------------------------------------------------------------
primers_lvl0<-primer_add_level(primers,  prefix="TT", restriction_enzyme=recognition_site_bbsi, suffix="AA", vector=c("CTCA", "CTCG"))
primers_lvl0

## ------------------------------------------------------------------------
print_primer(primers_lvl0)

## ----eval=FALSE----------------------------------------------------------
#  sink("primers.txt", append=FALSE, split=FALSE)
#  print_primer(primers_lvl0)
#  sink()

