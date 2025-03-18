library(testthat)
library(phyloseq_duckdb)
library(ggplot2)

test_that("psmelt_duckdb works correctly", {
  # Test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  colnames(otu_mat) <- paste0("Sample", 1:10)
  rownames(otu_mat) <- paste0("OTU", 1:10)
  
  tax_mat <- matrix(c(rep("Firmicutes", 5), rep("Bacteroidetes", 5),
                     rep(c("A", "B"), each = 5)), 
                   nrow = 10)
  colnames(tax_mat) <- c("Phylum", "Class")
  rownames(tax_mat) <- rownames(otu_mat)
  
  sam_df <- data.frame(
    Sample = paste0("Sample", 1:10),
    Treatment = rep(c("Control", "Treatment"), each = 5)
  )
  rownames(sam_df) <- sam_df$Sample
  
  # Create phyloseq object
  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat),
    sample_data(sam_df)
  )
  
  # Test melting
  melted <- psmelt_duckdb(ps)
  expect_true(is.data.frame(melted))
  expect_true(all(c("Abundance", "OTU", "Sample", "Phylum", "Treatment") %in% colnames(melted)))
  
  # Clean up
  close_phyloseq(ps)
})

test_that("plot_bar_duckdb works", {
  # Create test phyloseq object
  ps <- phyloseq(
    otu_table(matrix(sample(0:100, 100), nrow = 10), taxa_are_rows = TRUE),
    tax_table(matrix(c(rep("Firmicutes", 5), rep("Bacteroidetes", 5)), nrow = 10, dimnames = list(paste0("OTU", 1:10), "Phylum")))
  )
  
  # Test basic bar plot
  p1 <- plot_bar_duckdb(ps, fill = "Phylum")
  expect_s3_class(p1, "ggplot")
  
  # Test stacked bar plot
  p2 <- plot_bar_duckdb(ps, fill = "Phylum", normalize = TRUE)
  expect_s3_class(p2, "ggplot")
  
  # Clean up
  close_phyloseq(ps)
})

test_that("plot_ordination_duckdb works", {
  skip_if_not_installed("vegan")
  
  # Create test phyloseq object
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  colnames(otu_mat) <- paste0("Sample", 1:10)
  rownames(otu_mat) <- paste0("OTU", 1:10)
  
  sam_df <- data.frame(
    Sample = paste0("Sample", 1:10),
    Treatment = rep(c("Control", "Treatment"), each = 5)
  )
  rownames(sam_df) <- sam_df$Sample
  
  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    sample_data(sam_df)
  )
  
  # Test PCoA plot
  p1 <- plot_ordination_duckdb(ps, method = "PCoA", distance = "bray", color = "Treatment")
  expect_s3_class(p1, "ggplot")
  
  # Test NMDS plot
  p2 <- plot_ordination_duckdb(ps, method = "NMDS", distance = "bray", color = "Treatment")
  expect_s3_class(p2, "ggplot")
  
  # Clean up
  close_phyloseq(ps)
})
