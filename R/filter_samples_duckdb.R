# DuckDB-based implementation for filtering samples
#' Filter samples based on criteria using DuckDB
#'
#' This function filters samples based on metadata variables or abundance criteria,
#' efficiently using DuckDB's optimized query engine. Particularly useful for
#' subsetting large datasets based on experimental conditions or quality metrics.
#'
#' @param physeq A phyloseq_duckdb object
#' @param variable Name of the sample metadata variable to filter on
#' @param value Value to filter for
#' @param condition Comparison condition ('=', '>', '<', '>=', '<=', '!=')
#' @param min_abundance Minimum total abundance threshold for samples
#' @param min_taxa Minimum number of taxa present in samples
#' @return A new phyloseq_duckdb object with filtered samples
#' @examples
#' # Filter by metadata variable
#' ps_treated <- filter_samples_duckdb(ps, "Treatment", "Control")
#'
#' # Filter by abundance threshold
#' ps_high_counts <- filter_samples_duckdb(ps, 
#'                                        min_abundance = 1000,
#'                                        min_taxa = 10)
#'
#' # Filter by time point
#' ps_late <- filter_samples_duckdb(ps, "Time", 24, condition = ">=")
#' @seealso \code{\link{filter_taxa_duckdb}}, \code{\link{prune_taxa_duckdb}}
#' @export
filter_samples_duckdb <- function(physeq, variable = NULL, value = NULL, 
                                condition = "=", min_abundance = NULL, 
                                min_taxa = NULL) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Start building the sample filter conditions
    conditions <- character()
    
    # Add metadata condition if specified
    if (!is.null(variable) && !is.null(value)) {
        if (!DBI::dbExistsTable(physeq, "sample_data")) {
            stop("No sample data found in the phyloseq object")
        }
        
        # Handle different value types
        if (is.character(value) || is.factor(value)) {
            value_str <- sprintf("'%s'", value)
        } else {
            value_str <- as.character(value)
        }
        
        # Create temporary table with filtered samples
        DBI::dbExecute(physeq, sprintf("
            CREATE TEMP TABLE sample_filter AS
            SELECT sample_id
            FROM sample_data
            WHERE %s %s %s
        ", variable, condition, value_str))
        
        conditions <- c(conditions, "sample_id IN (SELECT sample_id FROM sample_filter)")
    }
    
    # Add abundance condition if specified
    if (!is.null(min_abundance)) {
        abundance_cond <- sprintf("
            sample_id IN (
                SELECT sample_id
                FROM otu_table
                GROUP BY sample_id
                HAVING SUM(abundance) >= %f
            )
        ", min_abundance)
        conditions <- c(conditions, abundance_cond)
    }
    
    # Add taxa count condition if specified
    if (!is.null(min_taxa)) {
        taxa_cond <- sprintf("
            sample_id IN (
                SELECT sample_id
                FROM otu_table
                WHERE abundance > 0
                GROUP BY sample_id
                HAVING COUNT(DISTINCT taxa_id) >= %d
            )
        ", min_taxa)
        conditions <- c(conditions, taxa_cond)
    }
    
    # Combine all conditions
    where_clause <- if (length(conditions) > 0) {
        paste("WHERE", paste(conditions, collapse = " AND "))
    } else {
        ""
    }
    
    # Create filtered OTU table
    DBI::dbExecute(physeq, sprintf("
        CREATE TEMP TABLE otu_table_filtered AS
        SELECT *
        FROM otu_table
        %s
    ", where_clause))
    
    # Update sample data if it exists
    if (DBI::dbExistsTable(physeq, "sample_data")) {
        DBI::dbExecute(physeq, "
            CREATE TEMP TABLE sample_data_filtered AS
            SELECT s.*
            FROM sample_data s
            INNER JOIN (
                SELECT DISTINCT sample_id
                FROM otu_table_filtered
            ) o ON s.sample_id = o.sample_id
        ")
        DBI::dbExecute(physeq, "DROP TABLE sample_data")
        DBI::dbExecute(physeq, "ALTER TABLE sample_data_filtered RENAME TO sample_data")
    }
    
    # Replace OTU table
    DBI::dbExecute(physeq, "DROP TABLE otu_table")
    DBI::dbExecute(physeq, "ALTER TABLE otu_table_filtered RENAME TO otu_table")
    
    # Clean up temporary tables
    if (DBI::dbExistsTable(physeq, "sample_filter")) {
        DBI::dbExecute(physeq, "DROP TABLE sample_filter")
    }
    
    # Recreate indices
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_sample ON otu_table(sample_id)")
    if (DBI::dbExistsTable(physeq, "sample_data")) {
        DBI::dbExecute(physeq, "CREATE INDEX idx_sample_id ON sample_data(sample_id)")
    }
    
    return(physeq)
}
