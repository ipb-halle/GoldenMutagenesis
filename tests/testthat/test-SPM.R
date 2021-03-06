context("Single Point Mutagenesis")


#load("../../data/Point_Mutagenesis_BbsI_setup.RData")
load(file.path(system.file("data", package="GoldenMutagenesis"), "Point_Mutagenesis_BbsI_setup.RData"))

primer_test<-mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = mutations, binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf)
#load("../../data/Point_Mutagenesis_BbsI_result.RData")
load(file.path(system.file("data", package="GoldenMutagenesis"), "Point_Mutagenesis_BbsI_result.RData"))

test_that("Known primers are calculated correctly (SPM)", {
  expect_that(primer_test, is_identical_to(primers))
})

primer_test_lc<-mutate_spm(input_sequence, prefix=str_to_lower("TT"), restriction_enzyme = str_to_lower(recognition_site_bbsi), suffix = str_to_lower("AA"), vector=str_to_lower(c("CTCA", "CTCG")), replacements = mutations, binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf)

test_that("Known primers are calculated correctly (SPM) for lowercase sequences", {
  expect_that(primer_test_lc, is_identical_to(primer_test))
})

load("SPM_BbsI_result_lv2.RData")
primer_test_lvl2<-primer_prepare_level(primers)

test_that("Preparing for level2 works", {
  expect_that(primer_test_lvl2, is_identical_to(primers_test))
})

rm(primers)
load("SPM_5V.RData")
load("SPM_5V_240E.RData")
load("SPM_5V_240E_212E.RData")

test_that("Mutations on start and stop positions are calculated correctly", {
  expect_that(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf), is_identical_to(primers_5v))
  expect_that(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V"), c("240", "E")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf), is_identical_to(primers_5v_240e))
  expect_that(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V"), c("240", "E"), c("212", "E")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf), is_identical_to(primers_5v_240e_212e))
  expect_that(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf)@primers[[1]][[1]], is_identical_to(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V"), c("240", "E"), c("212", "E")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf)@primers[[1]][[1]]))
  expect_that(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V"), c("240", "E")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf)@primers[[1]][[2]]@binding_sequence, is_identical_to(mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, suffix = "AA", vector=c("CTCA", "CTCG"), replacements = list(c("5", "V"), c("240", "E"), c("212", "E")), binding_min_length=4 ,binding_max_length=9, target_temp=60, cuf=cuf)@primers[[2]][[2]]@binding_sequence))
  
})

test_that("Domestication is working for known sequences", {
  expect_that(domesticate(input_sequence, recognition_site_bbsi, cuf) ,is_identical_to(list(c(143, "K"))))
}) 

test_that("Domestication is working for lowercase sequences", {
  expect_that(domesticate(str_to_lower(input_sequence), str_to_lower(recognition_site_bbsi), cuf) ,is_identical_to(list(c(143, "K"))))
})

load("SPM_complex_input.RData")
test_that("Domestication is working for complex sequences", {
  expect_that(domesticate(sequence, restriction_enzyme = "GAAGAC"), is_identical_to(list(c(1247, "K"), c(116, "K"), c(1010, "S"))))
  expect_that(domesticate(sequence, restriction_enzyme = "GGTCTC"), is_identical_to(list(c(131, "L"), c(358, "L"), c(988, "R"), c(716, "D"), c(877, "L"))))
  expect_that(domesticate(sequence, restriction_enzyme = "CGTCTC"), is_identical_to(list(c(248, "E"))))
})

#Test for issue #29

test_that("Avoiding the introduction of unwanted restriction sites", {
  expect_that(length(domesticate(input_sequence = mutate_spm(input_sequence = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCGAGCTGCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA", restriction_enzyme = "GAAGAC", prefix = "TT", suffix="AA", vector = c("CTCA", "CTCG"), replacements = list(c("60", "T")))@newsequence, restriction_enzyme = "GGTCTC")), equals(0))
})
