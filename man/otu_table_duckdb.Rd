% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otu_table_duckdb.R
\name{otu_table_duckdb}
\alias{otu_table_duckdb}
\title{Create a DuckDB-based OTU table}
\usage{
otu_table_duckdb(object, taxa_are_rows, db_path = tempfile())
}
\arguments{
\item{object}{A matrix, data.frame, or existing otu_table object to convert}

\item{taxa_are_rows}{Logical indicating if taxa are represented as rows (TRUE) or columns (FALSE)}

\item{db_path}{Path to store the DuckDB database. Default is temporary}
}
\value{
A DuckDB connection with OTU table
}
\description{
This function creates and manages OTU (Operational Taxonomic Unit) tables using DuckDB,
which is optimized for handling large sparse datasets. The data is stored in a columnar format
that is particularly efficient for analytical queries on large datasets.
}
\examples{
# Create from matrix
mat <- matrix(sample(0:100, 100), 10, 10)
otu <- otu_table_duckdb(mat, taxa_are_rows = TRUE)

# Calculate taxa sums
sums <- taxa_sums_duckdb(otu)

# Don't forget to close the connection
close_otu_table_duckdb(otu)
}
\seealso{
\code{\link{taxa_sums_duckdb}}, \code{\link{sample_sums_duckdb}}
}
