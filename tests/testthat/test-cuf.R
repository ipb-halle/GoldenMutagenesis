context("Codon Usage Tables")
test_that("Codon Usage Tables can be found", {
expect_match(list_cu_table(),"Salinibacter_ruber_DSM_13855.csv", fixed=T, all=F)
expect_match(list_cu_table(),"arabidopsis.csv", fixed=T, all=F)
expect_match(list_cu_table(),"e_coli_316407.csv", fixed=T, all=F)
expect_match(list_cu_table(),"s_cerevisiae_4932.csv", fixed=T, all=F)
})
