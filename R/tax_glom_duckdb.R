# DuckDB-based implementation for taxonomic agglomeration
#' Agglomerate taxa by taxonomic rank using DuckDB
#'
#' This function combines taxa at a specified taxonomic rank, summing their abundances.
#' Uses DuckDB's efficient aggregation capabilities for memory-efficient processing of
#' large datasets.
#'
#' @param physeq A phyloseq_duckdb object
#' @param taxrank Character string specifying the taxonomic rank to glom by
#' @return A new phyloseq_duckdb object with agglomerated taxa
#' @examples
#' # Agglomerate at Phylum level
#' ps_phylum <- tax_glom_duckdb(ps, "Phylum")
#'
#' # Agglomerate at Class level
#' ps_class <- tax_glom_duckdb(ps, "Class")
#' @seealso \code{\link{filter_taxa_duckdb}}
#' @export
tax_glom_duckdb <- function(physeq, taxrank) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Check if tax_table exists
    if (!DBI::dbExistsTable(physeq, "tax_table")) {
        stop("No taxonomy table found in the phyloseq object")
    }
    
    # Create temporary tables for the agglomeration
    DBI::dbExecute(physeq, sprintf("
        CREATE TEMP TABLE tax_glom AS
        SELECT t.%s as taxa_id, o.sample_id, SUM(o.abundance) as abundance
        FROM otu_table o
        JOIN tax_table t ON o.taxa_id = t.taxa_id
        GROUP BY t.%s, o.sample_id
    ", taxrank, taxrank))
    
    # Create new taxonomy table with only the specified rank
    DBI::dbExecute(physeq, sprintf("
        CREATE TEMP TABLE tax_table_new AS
        SELECT DISTINCT %s as taxa_id, %s as %s
        FROM tax_table
    ", taxrank, taxrank, taxrank))
    
    # Replace original tables
    DBI::dbExecute(physeq, "DROP TABLE otu_table")
    DBI::dbExecute(physeq, "DROP TABLE tax_table")
    DBI::dbExecute(physeq, "ALTER TABLE tax_glom RENAME TO otu_table")
    DBI::dbExecute(physeq, "ALTER TABLE tax_table_new RENAME TO tax_table")
    
    # Recreate indices
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_sample ON otu_table(sample_id)")
    DBI::dbExecute(physeq, "CREATE INDEX idx_tax_taxa ON tax_table(taxa_id)")
    
    return(physeq)
}
