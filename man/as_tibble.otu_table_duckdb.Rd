% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyverse_accessors.R
\name{as_tibble.otu_table_duckdb}
\alias{as_tibble.otu_table_duckdb}
\alias{as_tibble.tax_table_duckdb}
\alias{as_tibble.sample_data_duckdb}
\title{Convert phyloseq_duckdb components to tibbles}
\usage{
\method{as_tibble}{otu_table_duckdb}(obj, name = "OTU_ID")

\method{as_tibble}{tax_table_duckdb}(obj, name = "Taxa_ID")

\method{as_tibble}{sample_data_duckdb}(obj, name = "Sample_ID")
}
\arguments{
\item{obj}{An object of class otu_table_duckdb, tax_table_duckdb, or sample_data_duckdb}

\item{name}{Optional name for the ID column (default varies by component type)}
}
\value{
A tibble containing the data with rownames as a proper column
}
\description{
These functions convert phyloseq_duckdb components (OTU table, taxonomy table,
or sample data) into tibbles, preserving rownames as a proper column.
}
