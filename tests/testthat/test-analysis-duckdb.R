library(testthat)
library(phyloseq_duckdb)

test_that("Taxonomic operations work", {
  # Test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  colnames(otu_mat) <- paste0("Sample", 1:10)
  rownames(otu_mat) <- paste0("OTU", 1:10)
  
  tax_mat <- matrix(c(rep("Firmicutes", 5), rep("Bacteroidetes", 5),
                     rep(c("A", "B"), each = 5),
                     rep(letters[1:5], 2)), 
                   nrow = 10)
  colnames(tax_mat) <- c("Phylum", "Class", "Genus")
  rownames(tax_mat) <- rownames(otu_mat)
  
  # Create phyloseq object
  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat)
  )
  
  # Test tax_glom
  ps_glom <- tax_glom(ps, taxrank = "Phylum")
  expect_equal(ntaxa(ps_glom), 2)  # Should have 2 phyla
  
  # Test prune_taxa
  ps_pruned <- prune_taxa(taxa_names(ps)[1:5], ps)
  expect_equal(ntaxa(ps_pruned), 5)
  
  # Clean up
  close_phyloseq(ps)
  close_phyloseq(ps_glom)
  close_phyloseq(ps_pruned)
})

test_that("Sample filtering works", {
  # Test data
  otu_mat <- matrix(sample(0:100, 100), nrow = 10)
  sam_df <- data.frame(
    Sample = paste0("S", 1:10),
    Treatment = rep(c("A", "B"), each = 5),
    Reads = sample(1000:5000, 10)
  )
  
  # Create phyloseq object
  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    sample_data(sam_df)
  )
  
  # Test filtering by metadata
  ps_treat <- filter_samples(ps, Treatment == "A")
  expect_equal(nsamples(ps_treat), 5)
  
  # Test filtering by abundance
  ps_high <- filter_samples(ps, Reads > 2000)
  expect_true(nsamples(ps_high) < nsamples(ps))
  
  # Clean up
  close_phyloseq(ps)
  close_phyloseq(ps_treat)
  close_phyloseq(ps_high)
})
