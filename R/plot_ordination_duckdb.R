# DuckDB-based implementation for ordination plots
#' Create ordination plots from phyloseq objects using DuckDB
#'
#' This function efficiently performs ordination analysis and creates plots using
#' DuckDB's optimized data handling. It supports various ordination methods and
#' visualization options, particularly suited for large sparse datasets.
#'
#' @param physeq A phyloseq_duckdb object
#' @param method Ordination method ('PCoA', 'NMDS', 'CCA', 'RDA')
#' @param distance Distance method for PCoA/NMDS ('bray', 'jaccard')
#' @param color Character, variable to map to color aesthetic
#' @param shape Character, variable to map to shape aesthetic
#' @param size Character, variable to map to size aesthetic
#' @param title Plot title
#' @return A ggplot2 object
#' @examples
#' # PCoA with Bray-Curtis distances
#' p1 <- plot_ordination(ps, method = "PCoA", distance = "bray",
#'                      color = "Treatment")
#'
#' # NMDS with sample metadata
#' p2 <- plot_ordination(ps, method = "NMDS", distance = "jaccard",
#'                      color = "Treatment", shape = "Time")
#'
#' # RDA with multiple aesthetic mappings
#' p3 <- plot_ordination(ps, method = "RDA",
#'                      color = "Treatment",
#'                      shape = "Time",
#'                      size = "TotalAbundance")
#' @seealso \code{\link{psmelt_duckdb}}, \code{\link{plot_bar_duckdb}}
#' @export
plot_ordination_duckdb <- function(physeq, method = "PCoA", 
                                 distance = "bray",
                                 color = NULL, shape = NULL, 
                                 size = NULL, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needed for plotting. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("Package 'vegan' needed for ordination. Please install it.",
             call. = FALSE)
    }
    
    # Get abundance matrix
    otu_matrix <- DBI::dbGetQuery(physeq, "
        SELECT sample_id, taxa_id, abundance
        FROM otu_table
    ")
    otu_matrix <- reshape2::dcast(otu_matrix, 
                                 sample_id ~ taxa_id, 
                                 value.var = "abundance", 
                                 fill = 0)
    row.names(otu_matrix) <- otu_matrix$sample_id
    otu_matrix$sample_id <- NULL
    
    # Get sample data if available
    sample_data <- NULL
    if (DBI::dbExistsTable(physeq, "sample_data")) {
        sample_data <- DBI::dbGetQuery(physeq, "SELECT * FROM sample_data")
        row.names(sample_data) <- sample_data$sample_id
        sample_data$sample_id <- NULL
    }
    
    # Perform ordination
    ord <- switch(method,
        "PCoA" = {
            # Calculate distance matrix
            dist_matrix <- vegan::vegdist(otu_matrix, method = distance)
            # Perform PCoA
            ord <- stats::cmdscale(dist_matrix, k = 2, eig = TRUE)
            # Extract scores
            data.frame(
                PC1 = ord$points[,1],
                PC2 = ord$points[,2],
                sample_id = row.names(ord$points)
            )
        },
        "NMDS" = {
            # Perform NMDS
            ord <- vegan::metaMDS(otu_matrix, distance = distance)
            # Extract scores
            data.frame(
                NMDS1 = ord$points[,1],
                NMDS2 = ord$points[,2],
                sample_id = row.names(ord$points)
            )
        },
        "RDA" = {
            if (is.null(sample_data)) {
                stop("Sample data required for RDA")
            }
            # Perform RDA
            ord <- vegan::rda(otu_matrix ~ ., data = sample_data)
            # Extract scores
            data.frame(
                RDA1 = vegan::scores(ord, display = "sites")[,1],
                RDA2 = vegan::scores(ord, display = "sites")[,2],
                sample_id = row.names(vegan::scores(ord, display = "sites"))
            )
        },
        "CCA" = {
            if (is.null(sample_data)) {
                stop("Sample data required for CCA")
            }
            # Perform CCA
            ord <- vegan::cca(otu_matrix ~ ., data = sample_data)
            # Extract scores
            data.frame(
                CCA1 = vegan::scores(ord, display = "sites")[,1],
                CCA2 = vegan::scores(ord, display = "sites")[,2],
                sample_id = row.names(vegan::scores(ord, display = "sites"))
            )
        },
        stop("Unsupported ordination method")
    )
    
    # Merge with sample data if available
    if (!is.null(sample_data)) {
        ord <- merge(ord, 
                    data.frame(sample_data, 
                              sample_id = row.names(sample_data)),
                    by = "sample_id")
    }
    
    # Create plot
    x_axis <- colnames(ord)[2]  # First ordination axis
    y_axis <- colnames(ord)[3]  # Second ordination axis
    
    p <- ggplot2::ggplot(ord, ggplot2::aes_string(x = x_axis, y = y_axis))
    
    # Add aesthetic mappings
    aes_params <- list()
    if (!is.null(color)) aes_params$color <- color
    if (!is.null(shape)) aes_params$shape <- shape
    if (!is.null(size)) aes_params$size <- size
    
    p <- p + do.call(ggplot2::geom_point, aes_params)
    
    # Add theme and labels
    p <- p + ggplot2::theme_bw()
    
    if (!is.null(title)) {
        p <- p + ggplot2::ggtitle(title)
    }
    
    # Add axis labels based on ordination method
    x_lab <- switch(method,
        "PCoA" = "PCoA1",
        "NMDS" = "NMDS1",
        "RDA" = "RDA1",
        "CCA" = "CCA1"
    )
    y_lab <- switch(method,
        "PCoA" = "PCoA2",
        "NMDS" = "NMDS2",
        "RDA" = "RDA2",
        "CCA" = "CCA2"
    )
    
    p <- p + ggplot2::xlab(x_lab) + ggplot2::ylab(y_lab)
    
    return(p)
}
