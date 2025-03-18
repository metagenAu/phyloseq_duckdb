# DuckDB-based implementation for handling sample metadata
#' Create a DuckDB-based sample data table
#'
#' This function creates and manages sample metadata using DuckDB for efficient storage
#' and querying of sample information. The data is stored in a columnar format
#' optimized for fast filtering and joining operations with OTU tables.
#'
#' @param object A data.frame or existing sample_data object
#' @param db_path Path to store the DuckDB database. Default is temporary
#' @return A DuckDB connection with sample data table
#' @examples
#' # Create from data frame
#' sample_df <- data.frame(
#'   sample_id = paste0("Sample", 1:5),
#'   Treatment = c("Control", "Treatment", "Control", 
#'                 "Treatment", "Control"),
#'   Time = c(0, 0, 24, 24, 48)
#' )
#' sdata <- sample_data_duckdb(sample_df)
#'
#' # Get unique treatments
#' treatments <- unique_sample_values_duckdb(sdata, "Treatment")
#'
#' # Don't forget to close the connection
#' close_sample_data_duckdb(sdata)
#' @seealso \code{\link{filter_sample_data_duckdb}}, \code{\link{unique_sample_values_duckdb}}
#' @export
sample_data_duckdb <- function(object, db_path = tempfile()) {
    # Initialize DuckDB connection
    con <- DBI::dbConnect(duckdb::duckdb(), db_path)
    
    # Convert input to data frame
    if (inherits(object, "sample_data")) {
        sample_df <- as.data.frame(object)
    } else if (is.data.frame(object)) {
        sample_df <- object
    } else {
        stop("Invalid input type for sample_data_duckdb")
    }
    
    # Add sample_id column if not present
    if (!"sample_id" %in% colnames(sample_df)) {
        sample_df$sample_id <- rownames(sample_df)
    }
    
    # Create DuckDB table
    DBI::dbWriteTable(con, "sample_data", sample_df, overwrite = TRUE)
    
    # Create index on sample_id
    DBI::dbExecute(con, "CREATE INDEX idx_sample_id ON sample_data(sample_id)")
    
    # Return connection
    class(con) <- c("sample_data_duckdb", class(con))
    return(con)
}

#' Filter sample data based on metadata variables
#'
#' Efficiently filter the sample data based on specific metadata variables and values.
#' Uses DuckDB's indexing for fast filtering operations.
#'
#' @param x A sample_data_duckdb connection
#' @param variable Variable name to filter on
#' @param value Value to filter for
#' @param condition Comparison condition ('=', '>', '<', etc.)
#' @return Filtered sample_data_duckdb connection
#' @examples
#' # Filter for control samples
#' controls <- filter_sample_data_duckdb(sdata, "Treatment", "Control")
#'
#' # Filter for samples after 24 hours
#' late_samples <- filter_sample_data_duckdb(sdata, "Time", 24, ">")
#' @seealso \code{\link{unique_sample_values_duckdb}}
#' @export
filter_sample_data_duckdb <- function(x, variable, value, condition = "=") {
    if (!inherits(x, "sample_data_duckdb")) {
        stop("Input must be a sample_data_duckdb object")
    }
    
    # Handle different value types
    if (is.character(value) || is.factor(value)) {
        value_sql <- sprintf("'%s'", value)
    } else {
        value_sql <- as.character(value)
    }
    
    # Create filtered table
    DBI::dbExecute(x, sprintf("
        CREATE TABLE sample_data_filtered AS
        SELECT *
        FROM sample_data
        WHERE %s %s %s
    ", variable, condition, value_sql))
    
    # Replace original table
    DBI::dbExecute(x, "DROP TABLE sample_data")
    DBI::dbExecute(x, "ALTER TABLE sample_data_filtered RENAME TO sample_data")
    
    # Recreate index
    DBI::dbExecute(x, "CREATE INDEX idx_sample_id ON sample_data(sample_id)")
    
    return(x)
}

#' Get unique values for a sample metadata variable
#'
#' Efficiently retrieve all unique values for a specific sample metadata variable.
#' Takes advantage of DuckDB's columnar storage for fast distinct value queries.
#'
#' @param x A sample_data_duckdb connection
#' @param variable Variable to get unique values for
#' @return Vector of unique values
#' @examples
#' # Get all unique treatments
#' treatments <- unique_sample_values_duckdb(sdata, "Treatment")
#'
#' # Get all timepoints
#' times <- unique_sample_values_duckdb(sdata, "Time")
#' @seealso \code{\link{filter_sample_data_duckdb}}
#' @export
unique_sample_values_duckdb <- function(x, variable) {
    if (!inherits(x, "sample_data_duckdb")) {
        stop("Input must be a sample_data_duckdb object")
    }
    
    result <- DBI::dbGetQuery(x, sprintf("
        SELECT DISTINCT %s
        FROM sample_data
        WHERE %s IS NOT NULL
        ORDER BY %s
    ", variable, variable, variable))
    
    return(result[[1]])
}

#' Merge sample data with OTU table
#'
#' Efficiently merge sample metadata with OTU abundance data using DuckDB's
#' optimized join operations. This is particularly useful for analyzing
#' abundance patterns across different experimental conditions.
#'
#' @param sample_data A sample_data_duckdb connection
#' @param otu_table An otu_table_duckdb connection
#' @return A DuckDB connection with merged data
#' @examples
#' # Merge sample metadata with OTU abundances
#' merged <- merge_sample_otu_duckdb(sdata, otu)
#' @seealso \code{\link{sample_data_duckdb}}, \code{\link{otu_table_duckdb}}
#' @export
merge_sample_otu_duckdb <- function(sample_data, otu_table) {
    if (!inherits(sample_data, "sample_data_duckdb") || !inherits(otu_table, "otu_table_duckdb")) {
        stop("Inputs must be sample_data_duckdb and otu_table_duckdb objects")
    }
    
    # Create merged view
    DBI::dbExecute(otu_table, "
        CREATE VIEW merged_data AS
        SELECT o.*, s.*
        FROM otu_table o
        JOIN sample_data s USING (sample_id)
    ")
    
    return(otu_table)
}

#' Close DuckDB connection and cleanup
#'
#' @param x A sample_data_duckdb connection
#' @export
close_sample_data_duckdb <- function(x) {
    if (!inherits(x, "sample_data_duckdb")) {
        stop("Input must be a sample_data_duckdb object")
    }
    DBI::dbDisconnect(x, shutdown = TRUE)
}
