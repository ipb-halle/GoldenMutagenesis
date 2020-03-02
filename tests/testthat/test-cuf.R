context("Codon Usage Tables")
test_that("Codon Usage Tables can be found", {
expect_that(list_cu_table(),is_equivalent_to(c( "arabidopsis.csv","e_coli_316407.csv","s_cerevisiae_4932.csv", "Salinibacter_ruber_DSM_13855.csv")))
})
