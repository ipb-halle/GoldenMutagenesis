context("Class Tests")

load("../../data/MSD_BsaI_result_lv2.RData")
output<-capture.output(print_primer(primers), file=NULL)
test_that("Printing of primers (MSD) works", {
expect_that(output, is_identical_to(readLines('primers_MSD.txt')))
})
rm(primers)

load("../../data/Point_Mutagenesis_BbsI_result.RData")
output<-capture.output(print_primer(primers), file=NULL)
test_that("Printing of primers (SPM) works", {
expect_that(output, is_identical_to(readLines('primers_SPM.txt', warn = F)))
})
rm(primers)