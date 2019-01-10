## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library("GoldenMutagenesis")

## ------------------------------------------------------------------------
input_sequence<-"ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA"
recognition_site_bbsi<-"GAAGAC"
recognition_site_bsai<-"GGTCTC"
cuf<-"e_coli_316407.csv"


## ------------------------------------------------------------------------
mutations_bbsi<-domesticate(input_sequence, recognition_site_bbsi, cuf)
mutations_bbsi
mutations_bsai<-domesticate(input_sequence, recognition_site_bsai, cuf)
mutations_bsai


## ------------------------------------------------------------------------
mutations<-c(list(c(66, "V")), mutations_bbsi)
primers<-mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = mutations, binding_min_length=4 ,primer_length=9, target_temp=60, cuf=cuf)
primers

## ------------------------------------------------------------------------
print_primer(primers)

## ----eval=FALSE----------------------------------------------------------
#  sink("primers.txt", append=FALSE, split=FALSE)
#  print_primer(primers)
#  sink()

