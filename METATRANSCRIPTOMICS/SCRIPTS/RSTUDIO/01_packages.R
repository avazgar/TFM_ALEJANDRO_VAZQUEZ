###########################################################################################################################
# Script: 01_packages.R
# Description:
#       - Installation of all the packages required for the transcriptomic analysis
#       - Define methods of installation (CRAN, Bioconductor, Github)
#       - Load all required libraries with brief description of their purpose
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metatranscriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###########################################################################################################################

# CRAN packages

cran_packages <- c(
# Data manipulation and management
  "tidyverse",      # Collection of packages for data manipulation and visualization
  "dplyr",          # Data manipulation (select, filter, mutate, summarise)
  "tidyr",          # Data cleaning and structuration
  "tibble",         # Creation and manipulation of Dataframes
  "stringr",        # Strings manipulation
  "data.table",     # Fast reading and manipulation of big tables
  "reshape2",       # Wide-long format conversion
  "here",           # File paths management
  "magrittr",       # Pipe operator (%>%) and fluid coding
  "forcats",        # Factors manipulation and ordination
  "purrr",          # Functional programming and iteration
  "rlang",          # Core tidyverse programming infrastructure
  "conflicted",     # Resolve conflicts between packages
# Import and export data
  "readr",	    # Fast import of text and tables
  "readxl",	    # Read Excel files
  "openxlsx",	    # Read/Write Excel files without Java
  "writexl",	    # Export Dataframes to Excel
# Visualization
  "ggplot2",	    # Flexible data visualization
  "ggpubr",	    # Utilities and improvements of ggplot graphs
"RColorBrewer",   # Color palettes
  "cowplot",        # Graph composition and annotations
  "patchwork",      # Alternative composition of plots for ggplot2
  "gridExtra",      # Management of multiple graphs
  "ggh4x",          # Advanced extensions of ggplot2
  "pheatmap",       # Heatmap visualization
# Ecology and statistics
  "vegan",          # Multivariate ecological analyses (PERMANOVA, ordinations)
  "rdacca.hp",      # Hierarchical partitioning and variance partitioning
# Phylogeny and networks
  "ape",            # Phylogenetic analysis and tree manipulation
  "ggtree",         # ggplot2 visualization of phylogenetic trees
  "igraph",         # Analysis and visualization of complex networks
  "ggraph",         # ggplot2 visualization of graphs and networks
  "data.tree",      # Hierarchical data representation in tree structures
# Set operations and diagrams
  "venn",           # Venn diagrams
  "ggVennDiagram",  # ggplot2-based Venn diagrams
  "UpSetR"          # Complex set intersections (UpSet plots)
)


missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(missing_cran)) install.packages(missing_cran)
lapply(cran_packages, library, character.only = TRUE)


# Bioconductor Packages

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

bioc_packages <- c(
  "edgeR",                  # Differential expression analysis for RNA-seq
  "limma",                  # Linear models and empirical Bayes for transcriptomics
  "DESeq2",                 # Differential expression analysis for counts
  "apeglm",                 # Shrinkage estimators for log fold changes
  "tximport",               # Import transcript quantification data (Salmon, Kallisto)
  "phyloseq",               # (Optional) Microbiome analysis (for integration)
  "treeio",                 # Read/annotate phylogenetic trees
  "DirichletMultinomial",   # Probabilistic models for microbiome data
  "mia",                    # Multi-omics analysis of microbiome data
  "FEAST",                  # Source contribution analysis of microbial communities
  "EnhancedVolcano",        # Publication-ready volcano plots
  "GO.db",                  # Gene Ontology annotation database
  "SummarizedExperiment"    # Data container for assays + metadata
)

missing_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(missing_bioc))  BiocManager::install(missing_bioc, ask = FALSE)
lapply(bioc_packages, library, character.only = TRUE)

# Github packages

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")


github_packages <- list(
  "jbisanz/qiime2R"             = "qiime2R",           # Import QIIME2 artifacts to R
  "waldronlab/lefser"           = "lefser",            # LEfSe implementation in R for biomarker detection
  "yiluheihei/microbiomeMarker" = "microbiomeMarker",  # Identification and visualization of microbial biomarkers
  "taowenmicro/ggClusterNet"    = "ggClusterNet",      # Co-occurrence network construction
  "zdk123/SpiecEasi"            = "SpiecEasi",         # Ecological networks inference
  "gavinsimpson/ggvegan"        = "ggvegan"            # ggplot2 interface for vegan
)

for (repo in names(github_packages)) {
  pkg <- github_packages[[repo]]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(repo, build_vignettes = FALSE, force = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Conflict resolution: prefer functions from base R and core packages
################################################################################

# edgeR vs SingleCellExperiment
conflicts_prefer(edgeR::cpm)

# base R vs Biostrings/generics
conflicts_prefer(base::intersect)

# dplyr vs stats
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# dplyr vs MASS
conflicts_prefer(dplyr::select)

# dplyr vs randomForest
conflicts_prefer(dplyr::mutate)

# tibble vs base
conflicts_prefer(tibble::rownames_to_column)
conflicts_prefer(tibble::column_to_rownames)

# readr vs utils
conflicts_prefer(readr::read_delim)
conflicts_prefer(readr::write_tsv)
conflicts_prefer(readr::write_csv)

# base R vs Biostrings/generics
conflicts_prefer(base::setdiff)

# base R vs Biostrings/generics/igraph
conflicts_prefer(base::union)
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
# Preferencias de funciones con conflicto
conflicts_prefer(stats::cor)        
conflicts_prefer(base::union)        
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
conflicts_prefer(edgeR::cpm)         
conflicts_prefer(dplyr::filter)     
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::count)     
conflicts_prefer(dplyr::left_join)
conflicts_prefer(dplyr::right_join)
