#' DuckDB-based OTU table class
#'
#' @slot con DuckDB connection
#' @slot taxa_are_rows Logical indicating if taxa are rows
#' @exportClass otu_table_duckdb
setClass("otu_table_duckdb",
         slots = c(
           con = "DBIConnection",
           taxa_are_rows = "logical"
         ))

#' Create a DuckDB-based OTU table
#'
#' This function creates and manages OTU (Operational Taxonomic Unit) tables using DuckDB,
#' which is optimized for handling large sparse datasets. The data is stored in a columnar format
#' that is particularly efficient for analytical queries on large datasets.
#'
#' @param object A matrix, data.frame, or existing otu_table object to convert
#' @param taxa_are_rows Logical indicating if taxa are represented as rows (TRUE) or columns (FALSE)
#' @param db_path Path to store the DuckDB database. Default is temporary
#' @return An otu_table_duckdb object
#' @examples
#' # Create from matrix
#' mat <- matrix(sample(0:100, 100), 10, 10)
#' otu <- otu_table_duckdb(mat, taxa_are_rows = TRUE)
#'
#' # Calculate taxa sums
#' sums <- taxa_sums_duckdb(otu)
#'
#' # Don't forget to close the connection
#' close_otu_table_duckdb(otu)
#' @seealso \code{\link{taxa_sums_duckdb}}, \code{\link{sample_sums_duckdb}}
#' @export
otu_table_duckdb <- function(object, taxa_are_rows, db_path = tempfile()) {
    # Initialize DuckDB connection
    con <- DBI::dbConnect(duckdb::duckdb(), db_path)
    
    # Convert input to appropriate format
    if (inherits(object, "otu_table")) {
        mat <- as(object, "matrix")
    } else if (inherits(object, "data.frame")) {
        mat <- as.matrix(object)
    } else if (is.matrix(object)) {
        mat <- object
    } else {
        DBI::dbDisconnect(con, shutdown = TRUE)
        stop("Invalid input type for otu_table_duckdb")
    }
    
    # Create long-format data frame for sparse representation
    if (taxa_are_rows) {
        taxa_ids <- rownames(mat)
        sample_ids <- colnames(mat)
    } else {
        taxa_ids <- colnames(mat)
        sample_ids <- rownames(mat)
        mat <- t(mat)  # Transpose to make taxa rows for consistent DB structure
    }
    
    # Convert to long format (sparse representation)
    sparse_data <- data.frame(
        taxa_id = rep(taxa_ids, ncol(mat)),
        sample_id = rep(sample_ids, each = nrow(mat)),
        abundance = as.vector(mat)
    )
    
    # Remove zero abundances to maintain sparsity
    sparse_data <- sparse_data[sparse_data$abundance != 0, ]
    
    # Create DuckDB table
    DBI::dbWriteTable(con, "otu_table", sparse_data, overwrite = TRUE)
    
    # Create indices for faster querying
    DBI::dbExecute(con, "CREATE INDEX idx_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(con, "CREATE INDEX idx_sample ON otu_table(sample_id)")
    
    # Store metadata
    DBI::dbWriteTable(con, "metadata", data.frame(
        key = c("taxa_are_rows"),
        value = c(taxa_are_rows)
    ), overwrite = TRUE)
    
    # Create and return S4 object
    new("otu_table_duckdb", con = con, taxa_are_rows = taxa_are_rows)
}

#' Close DuckDB connection for OTU table
#'
#' @param x An otu_table_duckdb object
#' @export
close_otu_table_duckdb <- function(x) {
    if (!is(x, "otu_table_duckdb")) {
        stop("Input must be an otu_table_duckdb object")
    }
    DBI::dbDisconnect(x@con, shutdown = TRUE)
}

#' Calculate taxa sums for DuckDB OTU table
#'
#' Efficiently calculates the sum of abundances for each taxa across all samples.
#' Takes advantage of DuckDB's columnar storage for fast aggregation.
#'
#' @param x An otu_table_duckdb object
#' @return Named numeric vector of taxa sums
#' @examples
#' # Calculate taxa sums
#' sums <- taxa_sums_duckdb(otu)
#' @seealso \code{\link{sample_sums_duckdb}}
#' @export
taxa_sums_duckdb <- function(x) {
    if (!is(x, "otu_table_duckdb")) {
        stop("Input must be an otu_table_duckdb object")
    }
    
    result <- DBI::dbGetQuery(x@con, "
        SELECT taxa_id, SUM(abundance) as sum
        FROM otu_table
        GROUP BY taxa_id
        ORDER BY taxa_id
    ")
    
    sums <- setNames(result$sum, result$taxa_id)
    return(sums)
}

#' Calculate sample sums for DuckDB OTU table
#'
#' @param x An otu_table_duckdb object
#' @return Named numeric vector of sample sums
#' @export
sample_sums_duckdb <- function(x) {
    if (!is(x, "otu_table_duckdb")) {
        stop("Input must be an otu_table_duckdb object")
    }
    
    result <- DBI::dbGetQuery(x@con, "
        SELECT sample_id, SUM(abundance) as sum
        FROM otu_table
        GROUP BY sample_id
        ORDER BY sample_id
    ")
    
    sums <- setNames(result$sum, result$sample_id)
    return(sums)
}

#' Transform OTU abundances using DuckDB
#'
#' Apply common transformations to abundance data, such as log or square root.
#' Operations are performed efficiently using DuckDB's SQL engine.
#'
#' @param x An otu_table_duckdb object
#' @param transform The transformation to apply ('log10', 'log2', 'sqrt')
#' @return A new otu_table_duckdb object with transformed values
#' @examples
#' # Log transform abundances
#' otu_log <- transform_sample_counts_duckdb(otu, "log10")
#' @seealso \code{\link{filter_taxa_duckdb}}
#' @export
transform_sample_counts_duckdb <- function(x, transform = "log10") {
    if (!is(x, "otu_table_duckdb")) {
        stop("Input must be an otu_table_duckdb object")
    }
    
    # Create new temporary table with transformation
    transform_sql <- switch(transform,
        "log10" = "log10(abundance + 1)",
        "log2" = "log2(abundance + 1)",
        "sqrt" = "sqrt(abundance)",
        stop("Unsupported transformation")
    )
    
    DBI::dbExecute(x@con, sprintf("
        CREATE TABLE otu_table_transformed AS
        SELECT taxa_id, sample_id, %s as abundance
        FROM otu_table
    ", transform_sql))
    
    # Replace original table
    DBI::dbExecute(x@con, "DROP TABLE otu_table")
    DBI::dbExecute(x@con, "ALTER TABLE otu_table_transformed RENAME TO otu_table")
    
    # Recreate indices
    DBI::dbExecute(x@con, "CREATE INDEX idx_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(x@con, "CREATE INDEX idx_sample ON otu_table(sample_id)")
    
    return(x)
}

#' Filter OTU table based on abundance thresholds using DuckDB
#'
#' @param x An otu_table_duckdb object
#' @param min_abundance Minimum abundance threshold
#' @param min_samples Minimum number of samples where taxa must be present
#' @return Filtered otu_table_duckdb object
#' @export
filter_taxa_duckdb <- function(x, min_abundance = 1, min_samples = 1) {
    if (!is(x, "otu_table_duckdb")) {
        stop("Input must be an otu_table_duckdb object")
    }
    
    # Create filtered table
    DBI::dbExecute(x@con, sprintf("
        CREATE TABLE otu_table_filtered AS
        SELECT o.*
        FROM otu_table o
        INNER JOIN (
            SELECT taxa_id
            FROM otu_table
            WHERE abundance >= %f
            GROUP BY taxa_id
            HAVING COUNT(DISTINCT sample_id) >= %d
        ) f ON o.taxa_id = f.taxa_id
    ", min_abundance, min_samples))
    
    # Replace original table
    DBI::dbExecute(x@con, "DROP TABLE otu_table")
    DBI::dbExecute(x@con, "ALTER TABLE otu_table_filtered RENAME TO otu_table")
    
    # Recreate indices
    DBI::dbExecute(x@con, "CREATE INDEX idx_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(x@con, "CREATE INDEX idx_sample ON otu_table(sample_id)")
    
    return(x)
}
