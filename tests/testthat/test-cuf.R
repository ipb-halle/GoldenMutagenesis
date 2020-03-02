context("Codon Usage Tables")
skip('skip')
test_that("Codon Usage Tables can be found", {
expect_that(list_cu_table(),is_equivalent_to(c( "Salinibacter_ruber_DSM_13855.csv","arabidopsis.csv","e_coli_316407.csv","s_cerevisiae_4932.csv")))
})
