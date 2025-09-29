###########################################################################################################################
# Script: 01_packages.R
# Description:
#       - Installation of all the packages required for the whole project
#       - Define methods of installation 
#       - Load all required libraries with brief description of their purpose
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###########################################################################################################################

# CRAN Packages

cran_packages <- c(
# Data manipulation and management
  "tidyverse",	    # Collection of packages for data manipulation and visualization
  "dplyr",	    # Data manipulation (select, filter, mutate, summarise)
  "tidyr",          # Data cleaning and structuration
  "tibble",         # Dataframes creation and manipulation
  "stringr",        # Strings manipulation
  "data.table",     # Reading and manipulation of big tables
  "reshape2",       # Wide-long format conversion
  "here",           # File paths management
  "magrittr",       # Pipe operator (%>%) and fluid reading
  "forcats",	    # Factors manipulation and ordination		
# Import and export data
  "readr",          # Fast import of text and tables
  "readxl",         # Read Excel files
  "openxlsx",       # Reading/Writing Excel files without Java
  "writexl",        # Export as Excel file
# Visualization    
  "ggplot2",        # Flexible data visualization
  "ggpubr",         # Utilities and improvements for ggplot2 graphs
  "RColorBrewer",   # Color palettes
  "cowplot",        # Graph composition and anotations
  "patchwork",      # Alternative composition of plots for ggplot2
  "gridExtra",      # Management of high number of graphs
  "ggh4x",          # Advanced extensions of ggplot2
# Microbiome and ecology
  "vegan",          # Alfa and beta diversity analysis, ecological ordinations
  "microeco",       # Complete flow of microbiome analysis (main package)
  "compositions",   # Statistic tools for compositional data (CLR, ILR, etc.)
# Filogeny and networks
  "ape",            # Phylogenetic analysis and tree manipulation
  "ggtree",         # ggplot2 visualization of phylogenetic trees
  "igraph",         # Analysis and visualization of complex networks
  "ggraph",         # ggplot2 visualization of graphs and networks
  "data.tree",      # Hierarchical data representation in tree structures
# Machine learning and prediction
  "randomForest",   # Classification and prediction
  "caret",          # Unified framework for machine learning
  "pROC",           # ROC curves and metrics of predictive models
# Installation of complementary packages
  "BiocManager",     # Get access to Bioconductor to install additional packages
  "remotes"	    # Install packages from GitHub, GitLab, Bitbucket
)

# Install missing CRAN Packages
missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if (length(missing_cran)) install.packages(missing_cran)

# Load CRAN packages
lapply(cran_packages, library, character.only = TRUE)

# Bioconductor Packages:
bioc_packages <- c(
  "DESeq2",			# Diferential Analysis of genes/OTUs
  "treeio",			# Reading and anotation of phylogenetic trees
  "microbiome",			# Microbial communities analysis
  "phyloseq",			# Microbiome analysis based on OTUs/ASVs
  "apeglm",			# Adjusted log-fold changes estimation
  "DirichletMultinomial",	# Probabilistic models for microbiomes
  "mia",			# Multi-omic analysis for microbiomes
  "file2meco",			# Data conversion to microtable format (microeco)
  "FEAST",			# Source and contribution of microbial communities
  "SummarizedExperiment"	# LEfser data structure	
)

install_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if (length(install_bioc_packages)) {
  BiocManager::install(install_bioc_packages, ask = FALSE)
}
# Load Bioconductor packages
lapply(bioc_packages, library, character.only = TRUE)


# Github packages
github_packages <- list(
  "microbiota/amplicon"			= "amplicon",		# R pipeline for amplicon data analysis: It includes data preprocessing, diversity assessment, taxonomy and visualization
  "jbisanz/qiime2R"			= "qiime2R",		# It imports directly qiime2 artifacts to R objects
  "waldronlab/lefser"			= "lefser",		# Lefse implementation in R for biomarkers detection
  "yiluheihei/microbiomeMarker"		= "microbiomeMarker",	# Identification and visualization of microbial biomarkers
  "taowenmicro/ggClusterNet"		= "ggClusterNet",	# Construction and visualization of co-occurrence networks
  "zdk123/SpiecEasi"			= "SpiecEasi"		# Ecological networks inference from abundance data of microbiomes
)

for (repo in names(github_packages)) {
	pkg <- github_packages[[repo]]
	if(!requireNamespace(pkg, quietly = TRUE)){
		remotes::install_github(repo, build_vignettes = FALSE, force = TRUE)
	}
	library(pkg, character.only = TRUE)
}

