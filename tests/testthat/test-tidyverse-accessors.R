library(testthat)
library(phyloseq_duckdb)
library(tibble)
library(dplyr)

test_that("OTU table conversion to tibble works", {
  # Create test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  rownames(otu_mat) <- paste0("OTU", 1:10)
  colnames(otu_mat) <- paste0("Sample", 1:10)
  
  # Create OTU table
  otu <- otu_table_duckdb(otu_mat, taxa_are_rows = TRUE)
  
  # Convert to tibble
  tbl <- as_tibble(otu)
  
  # Tests
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), ncol(otu_mat) + 1)  # +1 for OTU_ID column
  expect_equal(nrow(tbl), nrow(otu_mat))
  expect_equal(tbl$OTU_ID, rownames(otu_mat))
  
  # Test custom ID name
  tbl2 <- as_tibble(otu, name = "MyOTUs")
  expect_true("MyOTUs" %in% colnames(tbl2))
  
  # Clean up
  close_otu_table(otu)
})

test_that("Taxonomy table conversion to tibble works", {
  # Create test data
  tax_mat <- matrix(sample(letters, 50, replace = TRUE), nrow = 10)
  rownames(tax_mat) <- paste0("OTU", 1:10)
  colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family")
  
  # Create tax table
  tax <- tax_table_duckdb(tax_mat)
  
  # Convert to tibble
  tbl <- as_tibble(tax)
  
  # Tests
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), ncol(tax_mat) + 1)  # +1 for Taxa_ID column
  expect_equal(nrow(tbl), nrow(tax_mat))
  expect_equal(tbl$Taxa_ID, rownames(tax_mat))
  
  # Test custom ID name
  tbl2 <- as_tibble(tax, name = "MyTaxa")
  expect_true("MyTaxa" %in% colnames(tbl2))
  
  # Clean up
  close_tax_table(tax)
})

test_that("Sample data conversion to tibble works", {
  # Create test data
  sam_df <- data.frame(
    Treatment = rep(c("Control", "Treatment"), each = 5),
    Time = rep(c("0h", "24h", "48h", "72h", "96h"), 2),
    stringsAsFactors = FALSE
  )
  rownames(sam_df) <- paste0("Sample", 1:10)
  
  # Create sample data
  sam <- sample_data_duckdb(sam_df)
  
  # Convert to tibble
  tbl <- as_tibble(sam)
  
  # Tests
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), ncol(sam_df) + 1)  # +1 for Sample_ID column
  expect_equal(nrow(tbl), nrow(sam_df))
  expect_equal(tbl$Sample_ID, rownames(sam_df))
  
  # Test custom ID name
  tbl2 <- as_tibble(sam, name = "MySamples")
  expect_true("MySamples" %in% colnames(tbl2))
  
  # Clean up
  close_sample_data(sam)
})

test_that("Full phyloseq object component extraction works", {
  # Create test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  rownames(otu_mat) <- paste0("OTU", 1:10)
  colnames(otu_mat) <- paste0("Sample", 1:10)
  
  tax_mat <- matrix(sample(letters, 50, replace = TRUE), nrow = 10)
  rownames(tax_mat) <- rownames(otu_mat)
  colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family")
  
  sam_df <- data.frame(
    Treatment = rep(c("Control", "Treatment"), each = 5),
    Time = rep(c("0h", "24h", "48h", "72h", "96h"), 2),
    stringsAsFactors = FALSE
  )
  rownames(sam_df) <- colnames(otu_mat)
  
  # Create phyloseq object
  ps <- phyloseq_duckdb(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat),
    sample_data(sam_df)
  )
  
  # Test component extraction
  otu_tbl <- ps %>% as_tibble("otu")
  expect_s3_class(otu_tbl, "tbl_df")
  expect_true("OTU_ID" %in% colnames(otu_tbl))
  
  tax_tbl <- ps %>% as_tibble("tax", name = "MyTaxa")
  expect_s3_class(tax_tbl, "tbl_df")
  expect_true("MyTaxa" %in% colnames(tax_tbl))
  
  sam_tbl <- ps %>% as_tibble("sample")
  expect_s3_class(sam_tbl, "tbl_df")
  expect_true("Sample_ID" %in% colnames(sam_tbl))
  
  # Clean up
  close_phyloseq(ps)
})
