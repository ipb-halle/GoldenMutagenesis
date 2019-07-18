## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library("GoldenMutagenesis")

## ------------------------------------------------------------------------
input_sequence<-"ATGTCTCAGGTTCAGAGTGGCATTTTGCCAGAACATTGCCGCGCGGCGATTTGGATCGAAGCCAACGTGAAAGGGGAAGTTGACGCCCTGCGTGCGGCCAGTAAAACATTTGCCGACAAACTGGCAACTTTTGAAGCGAAATTCCCGGACGCGCATCTTGGTGCGGTGGTTGCCTTTGGTAACAACACCTGGCGCGCTCTGAGCGGCGGCGTTGGGGCAGAAGAGCTGAAAGATTTTCCGGGCTACGGTAAAGGCCTTGCGCCGACGACCCAGTTCGATGTGTTGATCCACATTCTTTCTCTGCGTCACGACGTAAACTTCTCTGTCGCCCAGGCGGCGATGGAAGCCTTTGGTGACTGCATTGAAGTGAAAGAAGAGATCCACGGCTTCCGTTGGGTTGAAGAGCGTGACCTGAGCGGCTTTGTTGACGGTACGGAAAACCCGGCGGGTGAAGAGACGCGTCGCGAAGTGGCGGTTATCAAAGACGGCGTGGATGCGGGCGGCAGCTATGTGTTTGTCCAGCGTTGGGAACACAACCTGAAGCAGCTCAACCGGATGAGCGTTCACGATCAGGAGATGGTGATCGGGCGCACCAAAGAGGCCAACGAAGAGATCGACGGCGACGAACGTCCGGAAACCTCTCACCTCACCCGCGTTGATCTGAAAGAAGATGGCAAAGGGCTGAAGATTGTTCGCCAGAGCCTGCCGTACGGCACTGCCAGTGGCACTCACGGTCTGTACTTCTGCGCCTACTGCGCGCGTCTGCATAACATTGAGCAGCAACTGCTGAGCATGTTTGGCGATACCGATGGTAAGCGTGATGCGATGTTGCGTTTCACCAAACCGGTAACCGGCGGCTATTATTTCGCACCGTCGCTGGACAAGTTGATGGCGCTGTAA"

mutations<-c(137,143,147,232,234)


## ----eval=FALSE----------------------------------------------------------
#  abfile<-"sequences/Yfex_0activesite_for_EF01147142.ab1")
#  base_distribution(input_sequence=input_sequence, ab1file=abfile, replacements=mutations)

## ---- echo=F-------------------------------------------------------------
abfile<-system.file("sequences", "Yfex_0activesite_for_EF01147142.ab1", package="GoldenMutagenesis")
base_distribution(input_sequence=input_sequence, ab1file=abfile, replacements=mutations)

## ----eval=FALSE----------------------------------------------------------
#  abfile<-"sequences/Yfex_activesite_rev_EF01147143.ab1")
#  base_distribution(input_sequence=input_sequence, ab1file=abfile, replacements=mutations)

## ---- echo=F-------------------------------------------------------------
abfile<-system.file("sequences", "Yfex_activesite_rev_EF01147143.ab1", package="GoldenMutagenesis")
base_distribution(input_sequence=input_sequence, ab1file=abfile, replacements=mutations)

