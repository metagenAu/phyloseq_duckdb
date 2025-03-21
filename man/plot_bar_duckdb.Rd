% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bar_duckdb.R
\name{plot_bar_duckdb}
\alias{plot_bar_duckdb}
\title{Create bar plots from phyloseq objects using DuckDB}
\usage{
plot_bar_duckdb(
  physeq,
  x = "Sample",
  y = "Abundance",
  fill = NULL,
  facet_grid = NULL,
  normalize = FALSE
)
}
\arguments{
\item{physeq}{A phyloseq_duckdb object}

\item{x}{Character, variable to map to x-axis}

\item{y}{Character, variable to map to y-axis (default "Abundance")}

\item{fill}{Character, variable to map to fill aesthetic (default NULL)}

\item{facet_grid}{Character vector of length 2 for faceting (default NULL)}

\item{normalize}{Logical, whether to normalize abundances (default FALSE)}
}
\value{
A ggplot2 object
}
\description{
This function efficiently creates bar plots from phyloseq objects using DuckDB's
optimized data handling. It supports various types of bar plots including
stacked, grouped, and faceted visualizations.
}
\examples{
# Basic abundance plot by sample
p1 <- plot_bar(ps, x = "Sample")

# Stacked bar plot by phylum
p2 <- plot_bar(ps, x = "Sample", fill = "Phylum")

# Normalized abundances with faceting
p3 <- plot_bar(ps, x = "Sample", fill = "Phylum",
               facet_grid = c("Time", "Treatment"),
               normalize = TRUE)
}
\seealso{
\code{\link{psmelt_duckdb}}, \code{\link{plot_ordination_duckdb}}
}
