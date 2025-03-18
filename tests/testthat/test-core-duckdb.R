library(testthat)
library(phyloseq_duckdb)

test_that("Core DuckDB object creation works", {
  # Test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  tax_mat <- matrix(sample(letters, 60, replace = TRUE), nrow = 10)
  sam_df <- data.frame(Sample = paste0("S", 1:10), 
                      Treatment = rep(c("A", "B"), each = 5))
  
  # Test object creation
  otu <- otu_table_duckdb(otu_mat, taxa_are_rows = TRUE)
  expect_s4_class(otu, "otu_table_duckdb")
  
  tax <- tax_table_duckdb(tax_mat)
  expect_s4_class(tax, "tax_table_duckdb")
  
  sam <- sample_data_duckdb(sam_df)
  expect_s4_class(sam, "sample_data_duckdb")
  
  # Test phyloseq object creation
  ps <- phyloseq_duckdb(otu, tax, sam)
  expect_s4_class(ps, "phyloseq_duckdb")
  
  # Clean up
  close_phyloseq(ps)
})

test_that("Function aliases work correctly", {
  # Test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  
  # Test alias functionality
  otu_duck <- otu_table_duckdb(otu_mat, taxa_are_rows = TRUE)
  otu_alias <- otu_table(otu_mat, taxa_are_rows = TRUE)
  
  expect_identical(class(otu_duck), class(otu_alias))
  
  # Clean up
  close_otu_table(otu_duck)
  close_otu_table(otu_alias)
})

test_that("Connection management works", {
  # Test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  otu <- otu_table_duckdb(otu_mat, taxa_are_rows = TRUE)
  
  # Test connection closing
  expect_true(is(otu, "otu_table_duckdb"))
  close_otu_table(otu)
  expect_error(access_connection(otu), "Connection closed")
})
