# Function aliases for DuckDB-based implementations
#' @title Function Aliases for DuckDB Implementation
#' @description
#' This file provides drop-in replacements for all core phyloseq functions,
#' allowing seamless use of the DuckDB-based implementation. Each function
#' is aliased to its DuckDB counterpart without the "_duckdb" suffix,
#' maintaining compatibility with existing phyloseq workflows while
#' providing the benefits of DuckDB's optimized storage and queries.
#'
#' @section Function Mapping:
#' Core Object Creation:
#' \itemize{
#'   \item \code{phyloseq()} -> \code{phyloseq_duckdb()}
#'   \item \code{otu_table()} -> \code{otu_table_duckdb()}
#'   \item \code{tax_table()} -> \code{tax_table_duckdb()}
#'   \item \code{sample_data()} -> \code{sample_data_duckdb()}
#' }
#'
#' Analysis Functions:
#' \itemize{
#'   \item \code{alpha_diversity()} -> \code{alpha_diversity_duckdb()}
#'   \item \code{beta_diversity()} -> \code{beta_diversity_duckdb()}
#'   \item \code{filter_taxa()} -> \code{filter_taxa_phyloseq_duckdb()}
#'   \item \code{filter_samples()} -> \code{filter_samples_duckdb()}
#'   \item \code{tax_glom()} -> \code{tax_glom_duckdb()}
#'   \item \code{prune_taxa()} -> \code{prune_taxa_duckdb()}
#'   \item \code{transform_sample_counts()} -> \code{transform_sample_counts_duckdb()}
#' }
#'
#' Visualization Functions:
#' \itemize{
#'   \item \code{psmelt()} -> \code{psmelt_duckdb()}
#'   \item \code{plot_bar()} -> \code{plot_bar_duckdb()}
#'   \item \code{plot_ordination()} -> \code{plot_ordination_duckdb()}
#' }
#'
#' Resource Management:
#' \itemize{
#'   \item \code{close_phyloseq()} -> \code{close_phyloseq_duckdb()}
#'   \item \code{close_otu_table()} -> \code{close_otu_table_duckdb()}
#'   \item \code{close_tax_table()} -> \code{close_tax_table_duckdb()}
#'   \item \code{close_sample_data()} -> \code{close_sample_data_duckdb()}
#' }
#'
#' @note
#' It's important to use the close_* functions to properly manage DuckDB
#' connections and clean up resources when you're done with an object.
#' @name duckdb-aliases

#' @rdname duckdb-aliases
#' @export
phyloseq <- function(otu_table = NULL, 
                    tax_table = NULL, 
                    sample_data = NULL,
                    taxa_are_rows = TRUE,
                    db_path = tempfile()) {
    phyloseq_duckdb(otu_table, tax_table, sample_data, taxa_are_rows, db_path)
}

#' @rdname duckdb-aliases
#' @export
otu_table <- function(object, taxa_are_rows = NULL, db_path = tempfile()) {
    otu_table_duckdb(object, taxa_are_rows, db_path)
}

#' @rdname duckdb-aliases
#' @export
tax_table <- function(object, db_path = tempfile()) {
    tax_table_duckdb(object, db_path)
}

#' @rdname duckdb-aliases
#' @export
sample_data <- function(object, db_path = tempfile()) {
    sample_data_duckdb(object, db_path)
}

#' @rdname duckdb-aliases
#' @export
alpha_diversity <- function(physeq, measures = c("observed", "shannon", "simpson")) {
    alpha_diversity_duckdb(physeq, measures)
}

#' @rdname duckdb-aliases
#' @export
beta_diversity <- function(physeq, method = "bray") {
    beta_diversity_duckdb(physeq, method)
}

#' @rdname duckdb-aliases
#' @export
filter_taxa <- function(physeq, min_abundance = 1, min_prevalence = 1) {
    filter_taxa_phyloseq_duckdb(physeq, min_abundance, min_prevalence)
}

#' @rdname duckdb-aliases
#' @export
transform_sample_counts <- function(x, transform = "log10") {
    transform_sample_counts_duckdb(x, transform)
}

#' @rdname duckdb-aliases
#' @export
close_phyloseq <- function(physeq) {
    close_phyloseq_duckdb(physeq)
}

#' @rdname duckdb-aliases
#' @export
close_otu_table <- function(x) {
    close_otu_table_duckdb(x)
}

#' @rdname duckdb-aliases
#' @export
close_tax_table <- function(x) {
    close_tax_table_duckdb(x)
}

#' @rdname duckdb-aliases
#' @export
close_sample_data <- function(x) {
    close_sample_data_duckdb(x)
}

#' @rdname duckdb-aliases
#' @export
filter_sample_data <- function(x, variable, value, condition = "=") {
    filter_sample_data_duckdb(x, variable, value, condition)
}

#' @rdname duckdb-aliases
#' @export
unique_taxa <- function(x, rank) {
    unique_taxa_duckdb(x, rank)
}

#' @rdname duckdb-aliases
#' @export
unique_sample_values <- function(x, variable) {
    unique_sample_values_duckdb(x, variable)
}

#' @rdname duckdb-aliases
#' @export
merge_sample_otu <- function(sample_data, otu_table) {
    merge_sample_otu_duckdb(sample_data, otu_table)
}

#' @rdname duckdb-aliases
#' @export
taxa_sums <- function(x) {
    taxa_sums_duckdb(x)
}

#' @rdname duckdb-aliases
#' @export
sample_sums <- function(x) {
    sample_sums_duckdb(x)
}

#' @rdname duckdb-aliases
#' @export
tax_glom <- function(physeq, taxrank) {
    tax_glom_duckdb(physeq, taxrank)
}

#' @rdname duckdb-aliases
#' @export
prune_taxa <- function(taxa, physeq, keep = TRUE) {
    prune_taxa_duckdb(taxa, physeq, keep)
}

#' @rdname duckdb-aliases
#' @export
filter_samples <- function(physeq, variable = NULL, value = NULL, 
                         condition = "=", min_abundance = NULL, 
                         min_taxa = NULL) {
    filter_samples_duckdb(physeq, variable, value, condition, 
                         min_abundance, min_taxa)
}

#' @rdname duckdb-aliases
#' @export
psmelt <- function(physeq) {
    psmelt_duckdb(physeq)
}

#' @rdname duckdb-aliases
#' @export
plot_bar <- function(physeq, x = "Sample", y = "Abundance", 
                    fill = NULL, facet_grid = NULL,
                    normalize = FALSE) {
    plot_bar_duckdb(physeq, x, y, fill, facet_grid, normalize)
}

#' @rdname duckdb-aliases
#' @export
plot_ordination <- function(physeq, method = "PCoA", 
                          distance = "bray",
                          color = NULL, shape = NULL, 
                          size = NULL, title = NULL) {
    plot_ordination_duckdb(physeq, method, distance,
                          color, shape, size, title)
}
