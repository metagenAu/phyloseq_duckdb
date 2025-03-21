% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_duckdb.R
\name{filter_taxa_phyloseq_duckdb}
\alias{filter_taxa_phyloseq_duckdb}
\title{Filter phyloseq object based on abundance and prevalence}
\usage{
filter_taxa_phyloseq_duckdb(physeq, min_abundance = 1, min_prevalence = 1)
}
\arguments{
\item{physeq}{A phyloseq_duckdb object}

\item{min_abundance}{Minimum abundance threshold}

\item{min_prevalence}{Minimum number of samples where taxa must be present}
}
\value{
Filtered phyloseq_duckdb object
}
\description{
Efficiently filter taxa based on abundance and prevalence criteria using DuckDB's
optimized query engine. This is particularly useful for removing rare taxa or
focusing on core microbiome members in large datasets.
}
\examples{
# Remove rare taxa (present in less than 5 samples)
ps_filtered <- filter_taxa_phyloseq_duckdb(ps, 
                                          min_abundance = 1,
                                          min_prevalence = 5)

# Focus on abundant taxa (at least 100 reads in 10 samples)
ps_abundant <- filter_taxa_phyloseq_duckdb(ps, 
                                          min_abundance = 100,
                                          min_prevalence = 10)
}
\seealso{
\code{\link{transform_sample_counts_duckdb}}
}
