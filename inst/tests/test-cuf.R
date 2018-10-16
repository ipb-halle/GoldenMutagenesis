context("Codon Usage Tables")

expect_that(list_cu_table(),is_identical_to(c( "arabidopsis.csv","e_coli_316407.csv","s_cerevisiae_4932.csv")))