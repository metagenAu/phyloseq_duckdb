# DuckDB-based implementation for bar plots
#' Create bar plots from phyloseq objects using DuckDB
#'
#' This function efficiently creates bar plots from phyloseq objects using DuckDB's
#' optimized data handling. It supports various types of bar plots including
#' stacked, grouped, and faceted visualizations.
#'
#' @param physeq A phyloseq_duckdb object
#' @param x Character, variable to map to x-axis
#' @param y Character, variable to map to y-axis (default "Abundance")
#' @param fill Character, variable to map to fill aesthetic (default NULL)
#' @param facet_grid Character vector of length 2 for faceting (default NULL)
#' @param normalize Logical, whether to normalize abundances (default FALSE)
#' @return A ggplot2 object
#' @examples
#' # Basic abundance plot by sample
#' p1 <- plot_bar(ps, x = "Sample")
#'
#' # Stacked bar plot by phylum
#' p2 <- plot_bar(ps, x = "Sample", fill = "Phylum")
#'
#' # Normalized abundances with faceting
#' p3 <- plot_bar(ps, x = "Sample", fill = "Phylum",
#'                facet_grid = c("Time", "Treatment"),
#'                normalize = TRUE)
#' @seealso \code{\link{psmelt_duckdb}}, \code{\link{plot_ordination_duckdb}}
#' @export
plot_bar_duckdb <- function(physeq, x = "Sample", y = "Abundance", 
                           fill = NULL, facet_grid = NULL,
                           normalize = FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needed for plotting. Please install it.",
             call. = FALSE)
    }
    
    # Melt the phyloseq object
    melted_data <- psmelt_duckdb(physeq)
    
    # Normalize if requested
    if (normalize) {
        # Calculate total abundance per sample
        sample_totals <- tapply(melted_data$Abundance, 
                              melted_data$Sample, 
                              sum)
        # Normalize abundances
        melted_data$Abundance <- melted_data$Abundance / 
            sample_totals[match(melted_data$Sample, names(sample_totals))]
    }
    
    # Create base plot
    p <- ggplot2::ggplot(melted_data, 
                        ggplot2::aes_string(x = x, y = y))
    
    # Add bars with fill if specified
    if (!is.null(fill)) {
        p <- p + ggplot2::geom_bar(
            ggplot2::aes_string(fill = fill),
            stat = "identity",
            position = "stack"
        )
    } else {
        p <- p + ggplot2::geom_bar(
            stat = "identity",
            position = "stack"
        )
    }
    
    # Add faceting if specified
    if (!is.null(facet_grid) && length(facet_grid) == 2) {
        facet_formula <- stats::as.formula(
            paste(facet_grid[1], "~", facet_grid[2])
        )
        p <- p + ggplot2::facet_grid(facet_formula)
    }
    
    # Add theme and labels
    p <- p + ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    
    if (normalize) {
        p <- p + ggplot2::ylab("Relative Abundance")
    }
    
    return(p)
}
