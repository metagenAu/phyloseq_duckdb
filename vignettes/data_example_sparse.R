#' Generate Example Sparse Microbiome Dataset
#'
#' Creates a large sparse microbiome dataset (10000 samples x 50000 OTUs) with 90% sparsity,
#' designed to demonstrate the benefits of DuckDB optimization. The dataset includes taxonomic
#' assignments and rich sample metadata.
#'
#' @keywords internal
generate_example_sparse <- function() {
  # Create sparse OTU matrix (90% sparsity)
  set.seed(42)
  n_samples <- 10000
  n_otus <- 50000
  sparsity <- 0.90

  otu_mat <- matrix(0, nrow = n_otus, ncol = n_samples)
  non_zero_entries <- sample(1:(n_otus * n_samples), 
                           size = floor((n_otus * n_samples) * (1 - sparsity)),
                           replace = FALSE)
  otu_mat[non_zero_entries] <- stats::rpois(length(non_zero_entries), lambda = 10)

  rownames(otu_mat) <- paste0("OTU", 1:n_otus)
  colnames(otu_mat) <- paste0("Sample", 1:n_samples)

  # Create taxonomy table
  phyla <- c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria")
  classes <- paste0("Class", 1:8)
  orders <- paste0("Order", 1:16)
  families <- paste0("Family", 1:32)
  genera <- paste0("Genus", 1:64)

  tax_mat <- matrix(
    c(sample(phyla, n_otus, replace = TRUE),
      sample(classes, n_otus, replace = TRUE),
      sample(orders, n_otus, replace = TRUE),
      sample(families, n_otus, replace = TRUE),
      sample(genera, n_otus, replace = TRUE)),
    nrow = n_otus,
    ncol = 5
  )

  rownames(tax_mat) <- rownames(otu_mat)
  colnames(tax_mat) <- c("Phylum", "Class", "Order", "Family", "Genus")

  # Create sample metadata
  sam_df <- data.frame(
    Sample = paste0("Sample", 1:n_samples),
    Treatment = rep(c("Control", "DrugA", "DrugB", "DrugC"), each = n_samples/4),
    TimePoint = rep(rep(c("0h", "24h", "48h", "72h"), each = n_samples/16), 4),
    Subject = rep(1:(n_samples/8), each = 8),
    Reads = sample(5000:50000, n_samples, replace = TRUE),
    pH = stats::runif(n_samples, 6.0, 8.0),
    Temperature = stats::rnorm(n_samples, 37, 0.5)
  )
  rownames(sam_df) <- sam_df$Sample

  list(
    otu_table = otu_mat,
    tax_table = tax_mat,
    sample_data = sam_df
  )
}

#' Example Sparse Microbiome Dataset
#'
#' A large sparse microbiome dataset (10000 samples x 50000 OTUs) with 90% sparsity,
#' designed to demonstrate the benefits of DuckDB optimization. The dataset includes
#' taxonomic assignments and rich sample metadata.
#'
#' @format A list containing three components:
#' \describe{
#'   \item{otu_table}{A sparse OTU abundance matrix (50000 OTUs x 10000 samples)}
#'   \item{tax_table}{A taxonomy assignment matrix (50000 OTUs x 5 ranks)}
#'   \item{sample_data}{A data frame of sample metadata (10000 samples x 7 variables)}
#' }
#'
#' @examples
#' data(example_sparse)
#' ps <- phyloseq(
#'   otu_table(example_sparse$otu_table, taxa_are_rows = TRUE),
#'   tax_table(example_sparse$tax_table),
#'   sample_data(example_sparse$sample_data)
#' )
"example_sparse"
