# DuckDB-based implementation for pruning taxa
#' Prune specific taxa from a phyloseq object using DuckDB
#'
#' This function removes specified taxa from a phyloseq object, efficiently
#' using DuckDB's indexing and filtering capabilities. This is particularly
#' useful for removing unwanted taxa or focusing on a specific subset.
#'
#' @param taxa Character vector of taxa IDs to keep (if keep = TRUE) or remove (if keep = FALSE)
#' @param physeq A phyloseq_duckdb object
#' @param keep Logical, if TRUE keep only the specified taxa, if FALSE remove them
#' @return A new phyloseq_duckdb object with pruned taxa
#' @examples
#' # Keep only specific taxa
#' taxa_to_keep <- c("Firmicutes", "Bacteroidetes")
#' ps_pruned <- prune_taxa_duckdb(taxa_to_keep, ps, keep = TRUE)
#'
#' # Remove specific taxa
#' taxa_to_remove <- c("Unclassified", "Unknown")
#' ps_pruned <- prune_taxa_duckdb(taxa_to_remove, ps, keep = FALSE)
#' @seealso \code{\link{filter_taxa_duckdb}}, \code{\link{tax_glom_duckdb}}
#' @export
prune_taxa_duckdb <- function(taxa, physeq, keep = TRUE) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Convert taxa to string for SQL
    taxa_str <- paste0("'", paste(taxa, collapse = "','"), "'")
    
    # Create temporary OTU table with pruned taxa
    if (keep) {
        DBI::dbExecute(physeq, sprintf("
            CREATE TEMP TABLE otu_table_pruned AS
            SELECT *
            FROM otu_table
            WHERE taxa_id IN (%s)
        ", taxa_str))
    } else {
        DBI::dbExecute(physeq, sprintf("
            CREATE TEMP TABLE otu_table_pruned AS
            SELECT *
            FROM otu_table
            WHERE taxa_id NOT IN (%s)
        ", taxa_str))
    }
    
    # Update taxonomy table if it exists
    if (DBI::dbExistsTable(physeq, "tax_table")) {
        DBI::dbExecute(physeq, "
            CREATE TEMP TABLE tax_table_pruned AS
            SELECT t.*
            FROM tax_table t
            INNER JOIN (
                SELECT DISTINCT taxa_id
                FROM otu_table_pruned
            ) o ON t.taxa_id = o.taxa_id
        ")
        DBI::dbExecute(physeq, "DROP TABLE tax_table")
        DBI::dbExecute(physeq, "ALTER TABLE tax_table_pruned RENAME TO tax_table")
    }
    
    # Replace OTU table
    DBI::dbExecute(physeq, "DROP TABLE otu_table")
    DBI::dbExecute(physeq, "ALTER TABLE otu_table_pruned RENAME TO otu_table")
    
    # Recreate indices
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_sample ON otu_table(sample_id)")
    if (DBI::dbExistsTable(physeq, "tax_table")) {
        DBI::dbExecute(physeq, "CREATE INDEX idx_tax_taxa ON tax_table(taxa_id)")
    }
    
    return(physeq)
}
