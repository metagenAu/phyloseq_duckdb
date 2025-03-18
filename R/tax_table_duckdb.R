# DuckDB-based implementation for handling taxonomy tables
#' Create a DuckDB-based taxonomy table
#'
#' This function creates and manages taxonomy tables using DuckDB for efficient storage
#' and querying of taxonomic classifications. The data is stored in a columnar format
#' optimized for fast filtering and joining operations with OTU tables.
#'
#' @param object A matrix, data.frame, or existing tax_table object
#' @param db_path Path to store the DuckDB database. Default is temporary
#' @return A DuckDB connection with taxonomy table
#' @examples
#' # Create from data frame
#' tax_df <- data.frame(
#'   taxa_id = paste0("OTU", 1:5),
#'   Phylum = c("Firmicutes", "Bacteroidetes", "Firmicutes", 
#'              "Proteobacteria", "Actinobacteria"),
#'   Class = c("Bacilli", "Bacteroidia", "Clostridia",
#'             "Gammaproteobacteria", "Actinobacteria")
#' )
#' tax <- tax_table_duckdb(tax_df)
#'
#' # Get unique phyla
#' phyla <- unique_taxa_duckdb(tax, "Phylum")
#'
#' # Don't forget to close the connection
#' close_tax_table_duckdb(tax)
#' @seealso \code{\link{filter_tax_table_duckdb}}, \code{\link{unique_taxa_duckdb}}
#' @export
tax_table_duckdb <- function(object, db_path = tempfile()) {
    # Initialize DuckDB connection
    con <- DBI::dbConnect(duckdb::duckdb(), db_path)
    
    # Convert input to data frame
    if (inherits(object, "taxonomyTable")) {
        tax_df <- as.data.frame(object)
    } else if (is.matrix(object)) {
        tax_df <- as.data.frame(object)
    } else if (is.data.frame(object)) {
        tax_df <- object
    } else {
        stop("Invalid input type for tax_table_duckdb")
    }
    
    # Add taxa_id column if not present
    if (!"taxa_id" %in% colnames(tax_df)) {
        tax_df$taxa_id <- rownames(tax_df)
    }
    
    # Create DuckDB table
    DBI::dbWriteTable(con, "tax_table", tax_df, overwrite = TRUE)
    
    # Create index on taxa_id
    DBI::dbExecute(con, "CREATE INDEX idx_tax_taxa ON tax_table(taxa_id)")
    
    # Return connection
    class(con) <- c("tax_table_duckdb", class(con))
    return(con)
}

#' Filter taxonomy table based on taxonomic ranks
#'
#' Efficiently filter the taxonomy table based on specific taxonomic ranks and values.
#' Uses DuckDB's indexing for fast filtering operations.
#'
#' @param x A tax_table_duckdb connection
#' @param rank Taxonomic rank to filter on (e.g., "Phylum", "Class", "Order")
#' @param value Value to filter for
#' @return Filtered tax_table_duckdb connection
#' @examples
#' # Filter for Firmicutes
#' firmicutes <- filter_tax_table_duckdb(tax, "Phylum", "Firmicutes")
#' @seealso \code{\link{unique_taxa_duckdb}}
#' @export
filter_tax_table_duckdb <- function(x, rank, value) {
    if (!inherits(x, "tax_table_duckdb")) {
        stop("Input must be a tax_table_duckdb object")
    }
    
    # Create filtered table
    DBI::dbExecute(x, sprintf("
        CREATE TABLE tax_table_filtered AS
        SELECT *
        FROM tax_table
        WHERE %s = '%s'
    ", rank, value))
    
    # Replace original table
    DBI::dbExecute(x, "DROP TABLE tax_table")
    DBI::dbExecute(x, "ALTER TABLE tax_table_filtered RENAME TO tax_table")
    
    # Recreate index
    DBI::dbExecute(x, "CREATE INDEX idx_tax_taxa ON tax_table(taxa_id)")
    
    return(x)
}

#' Get unique values for a taxonomic rank
#'
#' Efficiently retrieve all unique values for a specific taxonomic rank.
#' Takes advantage of DuckDB's columnar storage for fast distinct value queries.
#'
#' @param x A tax_table_duckdb connection
#' @param rank Taxonomic rank to get unique values for (e.g., "Phylum", "Class")
#' @return Character vector of unique values
#' @examples
#' # Get all unique phyla
#' phyla <- unique_taxa_duckdb(tax, "Phylum")
#' @seealso \code{\link{filter_tax_table_duckdb}}
#' @export
unique_taxa_duckdb <- function(x, rank) {
    if (!inherits(x, "tax_table_duckdb")) {
        stop("Input must be a tax_table_duckdb object")
    }
    
    result <- DBI::dbGetQuery(x, sprintf("
        SELECT DISTINCT %s
        FROM tax_table
        WHERE %s IS NOT NULL
        ORDER BY %s
    ", rank, rank, rank))
    
    return(result[[1]])
}

#' Close DuckDB connection and cleanup
#'
#' @param x A tax_table_duckdb connection
#' @export
close_tax_table_duckdb <- function(x) {
    if (!inherits(x, "tax_table_duckdb")) {
        stop("Input must be a tax_table_duckdb object")
    }
    DBI::dbDisconnect(x, shutdown = TRUE)
}
