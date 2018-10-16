context("Multiple Site Saturation Mutagenesis")

load("../data/MSD_BsaI_setup_lv2.RData")
primer_test<-msd_mutate(input_sequence = input_sequence, replacements = mutations, restriction_enzyme = recognition_site_bsai)
load("../data/MSD_BsaI_result_lv2.RData")

test_that("Known primers are calculated correctly (MSD)", {
  expect_that(primer_test, is_identical_to(primers))
})

primer_test_lc<-msd_mutate(input_sequence = str_to_lower(input_sequence), replacements = mutations, restriction_enzyme = str_to_lower(recognition_site_bsai))

test_that("Known primers are calculated correctly (MSD) for lowercase sequences", {
  expect_that(primer_test_lc, is_identical_to(primers_test))
})

rm(primers)