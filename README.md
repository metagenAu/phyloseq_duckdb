<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# phyloseq_duckdb

A DuckDB-optimized implementation for handling and analyzing large-scale microbiome census data. This package maintains compatibility with the original phyloseq interface while providing significant performance improvements through columnar storage and sparse data optimization.

## Key Features

- Memory-efficient data transformation using DuckDB
- DuckDB-optimized joins for data melting operations
- Support for all standard phyloseq plotting options
- Full ggplot2 compatibility for customization
- Multiple ordination methods with vegan integration

## Core Functions

### Data Handling
- `otu_table_duckdb()`: DuckDB-based OTU table with sparse storage
- `tax_table_duckdb()`: Taxonomy table optimizations
- `sample_data_duckdb()`: Sample metadata handling
- `phyloseq_duckdb()`: Unified interface with optimized diversity calculations

### Analysis
- `tax_glom_duckdb()`: Agglomerate taxa by taxonomic rank
- `prune_taxa_duckdb()`: Keep or remove specific taxa
- `filter_samples_duckdb()`: Filter samples based on metadata or abundance

### Visualization
- `psmelt_duckdb()`: Efficiently melt phyloseq objects using DuckDB joins
- `plot_bar_duckdb()`: Create bar plots with support for stacking and faceting
- `plot_ordination_duckdb()`: Implement ordination methods (PCoA, NMDS, RDA, CCA)

## Function Aliases

All DuckDB functions have aliases without the "_duckdb" suffix for seamless transition:
- `otu_table()` → `otu_table_duckdb()`
- `tax_table()` → `tax_table_duckdb()`
- `sample_data()` → `sample_data_duckdb()`
- And more...

## Tidyverse Integration

`phyloseq_duckdb` provides seamless integration with the tidyverse ecosystem through `as_tibble` methods. You can easily convert OTU tables, taxonomy tables, and sample data into tibbles while preserving rownames:

```r
library(phyloseq_duckdb)
library(tidyverse)

# Load example data
data(example_sparse)
ps <- phyloseq(
  otu_table(example_sparse$otu_table, taxa_are_rows = TRUE),
  tax_table(example_sparse$tax_table),
  sample_data(example_sparse$sample_data)
)

# Convert components to tibbles
otu_tbl <- ps %>% as_tibble("otu")           # OTU table with OTU_ID column
tax_tbl <- ps %>% as_tibble("tax")           # Taxonomy table with Taxa_ID column
sam_tbl <- ps %>% as_tibble("sample")        # Sample data with Sample_ID column

# Customize ID column names
tax_tbl <- ps %>% as_tibble("tax", name = "MyTaxaID")

# Use with tidyverse functions
tax_tbl %>%
  filter(Phylum == "Firmicutes") %>%
  select(MyTaxaID, Family, Genus)
```

## Installation

```r
# Install from GitHub
devtools::install_github("metagenAu/phyloseq_duckdb")
```

## Usage

```r
library(phyloseq_duckdb)

# Create a phyloseq object with DuckDB backend
ps <- phyloseq(otu_table(counts_matrix),
               tax_table(taxonomy_matrix),
               sample_data(sample_info))

# Perform taxonomic agglomeration
ps_glom <- tax_glom(ps, taxrank = "Genus")

# Create a bar plot
plot_bar(ps_glom, fill = "Phylum")
```

## Performance Benefits

- Efficient handling of large, sparse datasets
- Reduced memory footprint through columnar storage
- Optimized joins and aggregations
- Stream processing for large datasets

## Dependencies

Core dependencies:
- DuckDB (>= 0.9.0)
- DBI (>= 1.1.0)
- ggplot2
- vegan
- reshape2
- patchwork
- Matrix

## License

AGPL-3

## Implementation Details

This package implements DuckDB-optimized versions of core phyloseq functionality:

1. Core Data Structures:
- Efficient columnar storage using DuckDB
- Sparse matrix optimization for OTU tables
- Automatic connection and resource management

2. Function Aliases:
- All DuckDB functions have drop-in replacements without the "_duckdb" suffix
- Maintains full compatibility with original phyloseq interface
- Proper cleanup with close_* functions

3. Performance Optimizations:
- Memory-efficient operations through DuckDB
- Optimized joins and aggregations
- Stream processing for large datasets
- Efficient indexing strategies

For detailed implementation information and benchmarks, please see the package vignettes.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Citation

If you use this package, please cite:

```
@software{phyloseq_duckdb2025,
  author = {{MetagenAu Team}},
  title = {phyloseq\_duckdb: DuckDB-optimized handling of high-throughput microbiome census data},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/metagenAu/phyloseq_duckdb}
}
