context("Multiple Site Saturation Mutagenesis")

#load("../../data/MSD_BsaI_setup_lv2.RData")
load(file.path(system.file("data", package="GoldenMutagenesis"), "MSD_BsaI_setup_lv2.RData"))
primer_test<-msd_mutate(input_sequence = input_sequence, replacements = mutations, restriction_enzyme = recognition_site_bsai)
#load("../../data/MSD_BsaI_result_lv2.RData")
load(file.path(system.file("data", package="GoldenMutagenesis"), "MSD_BsaI_result_lv2.RData"))

test_that("Known primers are calculated correctly (MSD)", {
  expect_that(primer_test, is_identical_to(primers))
})

primer_test_lc<-msd_mutate(input_sequence = str_to_lower(input_sequence), replacements = mutations, restriction_enzyme = str_to_lower(recognition_site_bsai))

test_that("Known primers are calculated correctly (MSD) for lowercase sequences", {
  expect_that(primer_test_lc, is_identical_to(primer_test))
})

load("MSD_BsaI_result_lv0.RData")
primer_test_lvl0<-primer_add_level(primers,  prefix="TT", restriction_enzyme="GAAGAC", suffix="AA", vector=c("CTCA", "CTCG"))

test_that("Adding level0 to a set of primers works", {
  expect_that(primer_test_lvl0, is_identical_to(primers_lvl0))
})

primer_test_lvl0_lc<-primer_add_level(primers,  prefix=str_to_lower("TT"), restriction_enzyme=str_to_lower("GAAGAC"), suffix=str_to_lower("AA"), vector=str_to_lower(c("CTCA", "CTCG")))

test_that("Adding level0 to a set of primers works for lowercase sequences", {
  expect_that(primer_test_lvl0, is_identical_to(primer_test_lvl0_lc))
})

rm(primers)