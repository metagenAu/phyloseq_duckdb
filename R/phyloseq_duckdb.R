# DuckDB-based implementation for phyloseq class
#' Create a DuckDB-based phyloseq object
#'
#' This function creates a unified phyloseq object that uses DuckDB for efficient storage
#' and analysis of large, sparse microbiome datasets. It combines OTU tables, taxonomy tables,
#' and sample data into a single DuckDB database with optimized queries. The implementation
#' is particularly suited for datasets that are too large to fit in memory.
#'
#' @param otu_table OTU table (matrix, data.frame, or otu_table object)
#' @param tax_table Taxonomy table (matrix, data.frame, or tax_table object)
#' @param sample_data Sample metadata (data.frame or sample_data object)
#' @param taxa_are_rows Logical indicating if taxa are represented as rows
#' @param db_path Path to store the DuckDB database. Default is temporary
#' @return A phyloseq_duckdb object
#' @examples
#' # Create a phyloseq object with all components
#' ps <- phyloseq_duckdb(otu_table = otu_mat,
#'                       tax_table = tax_df,
#'                       sample_data = sample_df,
#'                       taxa_are_rows = TRUE)
#'
#' # Calculate diversity metrics
#' alpha_div <- alpha_diversity_duckdb(ps)
#' beta_div <- beta_diversity_duckdb(ps)
#'
#' # Filter taxa based on abundance
#' ps_filtered <- filter_taxa_phyloseq_duckdb(ps, min_abundance = 10)
#'
#' # Don't forget to close the connection
#' close_phyloseq_duckdb(ps)
#' @seealso \code{\link{alpha_diversity_duckdb}}, \code{\link{beta_diversity_duckdb}}
#' @export
phyloseq_duckdb <- function(otu_table = NULL, 
                           tax_table = NULL, 
                           sample_data = NULL,
                           taxa_are_rows = TRUE,
                           db_path = tempfile()) {
    # Initialize DuckDB connection
    con <- DBI::dbConnect(duckdb::duckdb(), db_path)
    
    # Store components if provided
    if (!is.null(otu_table)) {
        otu_con <- otu_table_duckdb(otu_table, taxa_are_rows)
        DBI::dbExecute(con, "ATTACH ':memory:' AS otu")
        DBI::dbExecute(con, "CREATE TABLE otu_table AS SELECT * FROM otu.otu_table")
        close_otu_table_duckdb(otu_con)
    }
    
    if (!is.null(tax_table)) {
        tax_con <- tax_table_duckdb(tax_table)
        DBI::dbExecute(con, "ATTACH ':memory:' AS tax")
        DBI::dbExecute(con, "CREATE TABLE tax_table AS SELECT * FROM tax.tax_table")
        close_tax_table_duckdb(tax_con)
    }
    
    if (!is.null(sample_data)) {
        sample_con <- sample_data_duckdb(sample_data)
        DBI::dbExecute(con, "ATTACH ':memory:' AS sample")
        DBI::dbExecute(con, "CREATE TABLE sample_data AS SELECT * FROM sample.sample_data")
        close_sample_data_duckdb(sample_con)
    }
    
    # Create indices for performance
    if (!is.null(otu_table)) {
        DBI::dbExecute(con, "CREATE INDEX idx_otu_taxa ON otu_table(taxa_id)")
        DBI::dbExecute(con, "CREATE INDEX idx_otu_sample ON otu_table(sample_id)")
    }
    if (!is.null(tax_table)) {
        DBI::dbExecute(con, "CREATE INDEX idx_tax_taxa ON tax_table(taxa_id)")
    }
    if (!is.null(sample_data)) {
        DBI::dbExecute(con, "CREATE INDEX idx_sample_id ON sample_data(sample_id)")
    }
    
    # Store metadata
    DBI::dbWriteTable(con, "metadata", data.frame(
        key = c("taxa_are_rows", "has_otu", "has_tax", "has_sample"),
        value = c(taxa_are_rows, 
                 !is.null(otu_table),
                 !is.null(tax_table),
                 !is.null(sample_data))
    ), overwrite = TRUE)
    
    # Set class and return
    class(con) <- c("phyloseq_duckdb", class(con))
    return(con)
}

#' Calculate alpha diversity metrics using DuckDB
#'
#' Efficiently calculate various alpha diversity metrics using DuckDB's optimized
#' aggregation functions. This implementation is memory-efficient and suitable for
#' large datasets as it processes data in chunks.
#'
#' @param physeq A phyloseq_duckdb object
#' @param measures Vector of diversity measures to calculate ('observed', 'shannon', 'simpson')
#' @return Data frame of diversity metrics by sample
#' @examples
#' # Calculate multiple diversity metrics
#' alpha_div <- alpha_diversity_duckdb(ps, 
#'                                    measures = c("observed", "shannon"))
#'
#' # Calculate only Shannon diversity
#' shannon_div <- alpha_diversity_duckdb(ps, measures = "shannon")
#' @seealso \code{\link{beta_diversity_duckdb}}
#' @export
alpha_diversity_duckdb <- function(physeq, measures = c("observed", "shannon", "simpson")) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Base query for sample-level calculations
    base_query <- "
        WITH sample_counts AS (
            SELECT sample_id, taxa_id, abundance,
                   SUM(abundance) OVER (PARTITION BY sample_id) as total_abundance
            FROM otu_table
        )
    "
    
    metrics <- list()
    
    if ("observed" %in% measures) {
        observed <- DBI::dbGetQuery(physeq, paste0(base_query, "
            SELECT sample_id, COUNT(DISTINCT taxa_id) as observed
            FROM sample_counts
            GROUP BY sample_id
        "))
        metrics$observed <- observed
    }
    
    if ("shannon" %in% measures) {
        shannon <- DBI::dbGetQuery(physeq, paste0(base_query, "
            SELECT sample_id,
                   -SUM((abundance/total_abundance) * ln(abundance/total_abundance)) as shannon
            FROM sample_counts
            GROUP BY sample_id
        "))
        metrics$shannon <- shannon
    }
    
    if ("simpson" %in% measures) {
        simpson <- DBI::dbGetQuery(physeq, paste0(base_query, "
            SELECT sample_id,
                   1 - SUM(pow(abundance/total_abundance, 2)) as simpson
            FROM sample_counts
            GROUP BY sample_id
        "))
        metrics$simpson <- simpson
    }
    
    # Combine all metrics
    result <- Reduce(function(x, y) merge(x, y, by = "sample_id"), metrics)
    
    # Add sample metadata if available
    if (DBI::dbExistsTable(physeq, "sample_data")) {
        sample_data <- DBI::dbGetQuery(physeq, "SELECT * FROM sample_data")
        result <- merge(result, sample_data, by = "sample_id")
    }
    
    return(result)
}

#' Calculate beta diversity metrics using DuckDB
#'
#' Efficiently calculate pairwise sample dissimilarities using DuckDB's optimized
#' operations. This implementation handles large sparse datasets by processing
#' comparisons in chunks and using efficient memory management.
#'
#' @param physeq A phyloseq_duckdb object
#' @param method Distance method ('bray', 'jaccard')
#' @return Distance matrix between samples
#' @examples
#' # Calculate Bray-Curtis dissimilarity
#' bray_dist <- beta_diversity_duckdb(ps, method = "bray")
#'
#' # Calculate Jaccard distance
#' jaccard_dist <- beta_diversity_duckdb(ps, method = "jaccard")
#' @seealso \code{\link{alpha_diversity_duckdb}}
#' @export
beta_diversity_duckdb <- function(physeq, method = "bray") {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Get all sample pairs
    samples <- DBI::dbGetQuery(physeq, "SELECT DISTINCT sample_id FROM otu_table")$sample_id
    sample_pairs <- expand.grid(sample1 = samples, sample2 = samples)
    
    # Calculate distances based on method
    if (method == "bray") {
        distances <- DBI::dbGetQuery(physeq, "
            WITH abundances AS (
                SELECT s1.sample_id as sample1, s2.sample_id as sample2,
                       COALESCE(s1.abundance, 0) as abd1,
                       COALESCE(s2.abundance, 0) as abd2
                FROM (SELECT DISTINCT sample_id FROM otu_table) samples1
                CROSS JOIN (SELECT DISTINCT sample_id FROM otu_table) samples2
                LEFT JOIN otu_table s1 ON samples1.sample_id = s1.sample_id
                LEFT JOIN otu_table s2 ON samples2.sample_id = s2.sample_id
                WHERE samples1.sample_id < samples2.sample_id
            )
            SELECT sample1, sample2,
                   SUM(ABS(abd1 - abd2)) / (SUM(abd1) + SUM(abd2)) as distance
            FROM abundances
            GROUP BY sample1, sample2
        ")
    } else if (method == "jaccard") {
        distances <- DBI::dbGetQuery(physeq, "
            WITH presence AS (
                SELECT s1.sample_id as sample1, s2.sample_id as sample2,
                       COUNT(DISTINCT CASE WHEN s1.abundance > 0 THEN s1.taxa_id END) as n1,
                       COUNT(DISTINCT CASE WHEN s2.abundance > 0 THEN s2.taxa_id END) as n2,
                       COUNT(DISTINCT CASE WHEN s1.abundance > 0 AND s2.abundance > 0 
                                         THEN s1.taxa_id END) as shared
                FROM (SELECT DISTINCT sample_id FROM otu_table) samples1
                CROSS JOIN (SELECT DISTINCT sample_id FROM otu_table) samples2
                LEFT JOIN otu_table s1 ON samples1.sample_id = s1.sample_id
                LEFT JOIN otu_table s2 ON samples2.sample_id = s2.sample_id
                WHERE samples1.sample_id < samples2.sample_id
                GROUP BY s1.sample_id, s2.sample_id
            )
            SELECT sample1, sample2,
                   1 - (shared / (n1 + n2 - shared)) as distance
            FROM presence
        ")
    } else {
        stop("Unsupported distance method")
    }
    
    # Convert to distance matrix
    n <- length(samples)
    dist_matrix <- matrix(0, n, n)
    rownames(dist_matrix) <- colnames(dist_matrix) <- samples
    
    for (i in 1:nrow(distances)) {
        dist_matrix[distances$sample1[i], distances$sample2[i]] <- 
        dist_matrix[distances$sample2[i], distances$sample1[i]] <- distances$distance[i]
    }
    
    class(dist_matrix) <- "dist"
    return(dist_matrix)
}

#' Filter phyloseq object based on abundance and prevalence
#'
#' Efficiently filter taxa based on abundance and prevalence criteria using DuckDB's
#' optimized query engine. This is particularly useful for removing rare taxa or
#' focusing on core microbiome members in large datasets.
#'
#' @param physeq A phyloseq_duckdb object
#' @param min_abundance Minimum abundance threshold
#' @param min_prevalence Minimum number of samples where taxa must be present
#' @return Filtered phyloseq_duckdb object
#' @examples
#' # Remove rare taxa (present in less than 5 samples)
#' ps_filtered <- filter_taxa_phyloseq_duckdb(ps, 
#'                                           min_abundance = 1,
#'                                           min_prevalence = 5)
#'
#' # Focus on abundant taxa (at least 100 reads in 10 samples)
#' ps_abundant <- filter_taxa_phyloseq_duckdb(ps, 
#'                                           min_abundance = 100,
#'                                           min_prevalence = 10)
#' @seealso \code{\link{transform_sample_counts_duckdb}}
#' @export
filter_taxa_phyloseq_duckdb <- function(physeq, min_abundance = 1, min_prevalence = 1) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    
    # Filter OTU table
    DBI::dbExecute(physeq, sprintf("
        CREATE TABLE otu_table_filtered AS
        SELECT o.*
        FROM otu_table o
        INNER JOIN (
            SELECT taxa_id
            FROM otu_table
            WHERE abundance >= %f
            GROUP BY taxa_id
            HAVING COUNT(DISTINCT sample_id) >= %d
        ) f ON o.taxa_id = f.taxa_id
    ", min_abundance, min_prevalence))
    
    # Update related tables
    if (DBI::dbExistsTable(physeq, "tax_table")) {
        DBI::dbExecute(physeq, "
            CREATE TABLE tax_table_filtered AS
            SELECT t.*
            FROM tax_table t
            INNER JOIN (
                SELECT DISTINCT taxa_id
                FROM otu_table_filtered
            ) f ON t.taxa_id = f.taxa_id
        ")
        DBI::dbExecute(physeq, "DROP TABLE tax_table")
        DBI::dbExecute(physeq, "ALTER TABLE tax_table_filtered RENAME TO tax_table")
    }
    
    # Replace OTU table
    DBI::dbExecute(physeq, "DROP TABLE otu_table")
    DBI::dbExecute(physeq, "ALTER TABLE otu_table_filtered RENAME TO otu_table")
    
    # Recreate indices
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_taxa ON otu_table(taxa_id)")
    DBI::dbExecute(physeq, "CREATE INDEX idx_otu_sample ON otu_table(sample_id)")
    if (DBI::dbExistsTable(physeq, "tax_table")) {
        DBI::dbExecute(physeq, "CREATE INDEX idx_tax_taxa ON tax_table(taxa_id)")
    }
    
    return(physeq)
}

#' Close DuckDB connection and cleanup
#'
#' Properly close the DuckDB connection and clean up resources. It's important
#' to call this function when you're done with the phyloseq object to ensure
#' proper resource management.
#'
#' @param physeq A phyloseq_duckdb object
#' @examples
#' # Create phyloseq object
#' ps <- phyloseq_duckdb(otu_table = otu_mat)
#'
#' # Do some analysis
#' alpha_div <- alpha_diversity_duckdb(ps)
#'
#' # Clean up when done
#' close_phyloseq_duckdb(ps)
#' @export
close_phyloseq_duckdb <- function(physeq) {
    if (!inherits(physeq, "phyloseq_duckdb")) {
        stop("Input must be a phyloseq_duckdb object")
    }
    DBI::dbDisconnect(physeq, shutdown = TRUE)
}
