% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ordination_duckdb.R
\name{plot_ordination_duckdb}
\alias{plot_ordination_duckdb}
\title{Create ordination plots from phyloseq objects using DuckDB}
\usage{
plot_ordination_duckdb(
  physeq,
  method = "PCoA",
  distance = "bray",
  color = NULL,
  shape = NULL,
  size = NULL,
  title = NULL
)
}
\arguments{
\item{physeq}{A phyloseq_duckdb object}

\item{method}{Ordination method ('PCoA', 'NMDS', 'CCA', 'RDA')}

\item{distance}{Distance method for PCoA/NMDS ('bray', 'jaccard')}

\item{color}{Character, variable to map to color aesthetic}

\item{shape}{Character, variable to map to shape aesthetic}

\item{size}{Character, variable to map to size aesthetic}

\item{title}{Plot title}
}
\value{
A ggplot2 object
}
\description{
This function efficiently performs ordination analysis and creates plots using
DuckDB's optimized data handling. It supports various ordination methods and
visualization options, particularly suited for large sparse datasets.
}
\examples{
# PCoA with Bray-Curtis distances
p1 <- plot_ordination(ps, method = "PCoA", distance = "bray",
                     color = "Treatment")

# NMDS with sample metadata
p2 <- plot_ordination(ps, method = "NMDS", distance = "jaccard",
                     color = "Treatment", shape = "Time")

# RDA with multiple aesthetic mappings
p3 <- plot_ordination(ps, method = "RDA",
                     color = "Treatment",
                     shape = "Time",
                     size = "TotalAbundance")
}
\seealso{
\code{\link{psmelt_duckdb}}, \code{\link{plot_bar_duckdb}}
}
