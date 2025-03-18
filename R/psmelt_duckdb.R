# DuckDB-based implementation for melting phyloseq objects
#' Melt phyloseq object to data.frame using DuckDB
#'
#' This function efficiently converts a phyloseq object into a long-format data frame
#' using DuckDB's optimized join operations. The resulting data frame is suitable
#' for plotting with ggplot2 and other visualization tools.
#'
#' @param physeq A phyloseq_duckdb object
#' @return A data.frame with all phyloseq data in long format
#' @examples
#' # Melt phyloseq object for plotting
#' ps_melted <- psmelt(ps)
#'
#' # Create bar plot
#' library(ggplot2)
#' ggplot(ps_melted, aes(x = Sample, y = Abundance, fill = Phylum)) +
#'   geom_bar(stat = "identity")
#' @seealso \code{\link{plot_bar_duckdb}}, \code{\link{plot_ordination_duckdb}}
#' @export
psmelt_duckdb <- function(physeq) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Start with OTU table
    query <- "SELECT taxa_id, sample_id, abundance AS Abundance FROM otu_table"
    
    # Add taxonomy if available
    if (DBI::dbExistsTable(physeq, "tax_table")) {
        # Get taxonomy column names
        tax_cols <- DBI::dbGetQuery(physeq, "
            SELECT column_name 
            FROM information_schema.columns 
            WHERE table_name = 'tax_table' 
            AND column_name != 'taxa_id'
        ")$column_name
        
        tax_cols_str <- paste(tax_cols, collapse = ", ")
        query <- sprintf("
            SELECT o.sample_id as Sample, 
                   o.taxa_id as OTU,
                   o.abundance as Abundance,
                   %s
            FROM otu_table o
            LEFT JOIN tax_table t ON o.taxa_id = t.taxa_id
        ", tax_cols_str)
    }
    
    # Add sample data if available
    if (DBI::dbExistsTable(physeq, "sample_data")) {
        # Get sample data column names
        sample_cols <- DBI::dbGetQuery(physeq, "
            SELECT column_name 
            FROM information_schema.columns 
            WHERE table_name = 'sample_data' 
            AND column_name != 'sample_id'
        ")$column_name
        
        sample_cols_str <- paste(
            sprintf("s.%s", sample_cols),
            collapse = ", "
        )
        
        if (DBI::dbExistsTable(physeq, "tax_table")) {
            query <- sprintf("
                SELECT m.*, %s
                FROM (%s) m
                LEFT JOIN sample_data s ON m.Sample = s.sample_id
            ", sample_cols_str, query)
        } else {
            query <- sprintf("
                SELECT o.sample_id as Sample,
                       o.taxa_id as OTU,
                       o.abundance as Abundance,
                       %s
                FROM otu_table o
                LEFT JOIN sample_data s ON o.sample_id = s.sample_id
            ", sample_cols_str)
        }
    }
    
    # Execute query and return as data frame
    melted <- DBI::dbGetQuery(physeq, query)
    
    # Convert to data.frame and ensure proper column types
    melted <- as.data.frame(melted, stringsAsFactors = FALSE)
    
    return(melted)
}
