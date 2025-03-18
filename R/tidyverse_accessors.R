#' Convert phyloseq_duckdb components to tibbles
#'
#' These functions convert phyloseq_duckdb components (OTU table, taxonomy table,
#' or sample data) into tibbles, preserving rownames as a proper column.
#'
#' @param obj An object of class otu_table_duckdb, tax_table_duckdb, or sample_data_duckdb
#' @param name Optional name for the ID column (default varies by component type)
#' @return A tibble containing the data with rownames as a proper column
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom DBI dbGetQuery
#' @export
as_tibble.otu_table_duckdb <- function(obj, name = "OTU_ID") {
  # Get data from DuckDB connection
  data <- dbGetQuery(obj@conn, paste("SELECT * FROM", obj@table_name))
  
  # Convert to tibble, preserving rownames
  tbl <- tibble::as_tibble(data, rownames = name)
  
  # Return tibble
  tbl
}

#' @rdname as_tibble.otu_table_duckdb
#' @export
as_tibble.tax_table_duckdb <- function(obj, name = "Taxa_ID") {
  # Get data from DuckDB connection
  data <- dbGetQuery(obj@conn, paste("SELECT * FROM", obj@table_name))
  
  # Convert to tibble, preserving rownames
  tbl <- tibble::as_tibble(data, rownames = name)
  
  # Return tibble
  tbl
}

#' @rdname as_tibble.otu_table_duckdb
#' @export
as_tibble.sample_data_duckdb <- function(obj, name = "Sample_ID") {
  # Get data from DuckDB connection
  data <- dbGetQuery(obj@conn, paste("SELECT * FROM", obj@table_name))
  
  # Convert to tibble, preserving rownames
  tbl <- tibble::as_tibble(data, rownames = name)
  
  # Return tibble
  tbl
}

#' Convert phyloseq_duckdb object components to tibbles
#'
#' This function provides a convenient way to extract components from a phyloseq_duckdb
#' object as tibbles. It works with the tidyverse pipe operator.
#'
#' @param obj A phyloseq_duckdb object
#' @param component Which component to extract: "otu", "tax", or "sample"
#' @param name Optional name for the ID column
#' @return A tibble containing the requested component data
#' @examples
#' \dontrun{
#' # Create example phyloseq object
#' data(example_sparse)
#' ps <- phyloseq(
#'   otu_table(example_sparse$otu_table, taxa_are_rows = TRUE),
#'   tax_table(example_sparse$tax_table),
#'   sample_data(example_sparse$sample_data)
#' )
#'
#' # Extract components as tibbles
#' ps %>% as_tibble("otu")
#' ps %>% as_tibble("tax", name = "MyTaxaID")
#' ps %>% as_tibble("sample")
#' }
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @export
as_tibble.phyloseq_duckdb <- function(obj, component = c("otu", "tax", "sample"), name = NULL) {
  component <- match.arg(component)
  
  # Set default names if not provided
  if (is.null(name)) {
    name <- switch(component,
                  otu = "OTU_ID",
                  tax = "Taxa_ID",
                  sample = "Sample_ID")
  }
  
  # Extract and convert component
  switch(component,
         otu = as_tibble(otu_table(obj), name = name),
         tax = as_tibble(tax_table(obj), name = name),
         sample = as_tibble(sample_data(obj), name = name))
}

#' Register S3 methods
#' @importFrom tibble as_tibble
#' @importFrom methods setOldClass
#' @export
.onLoad <- function(libname, pkgname) {
  # Register S3 methods
  methods::setOldClass(c("otu_table_duckdb", "tax_table_duckdb", "sample_data_duckdb", "phyloseq_duckdb"))
  registerS3method("as_tibble", "otu_table_duckdb", as_tibble.otu_table_duckdb)
  registerS3method("as_tibble", "tax_table_duckdb", as_tibble.tax_table_duckdb)
  registerS3method("as_tibble", "sample_data_duckdb", as_tibble.sample_data_duckdb)
  registerS3method("as_tibble", "phyloseq_duckdb", as_tibble.phyloseq_duckdb)
}
