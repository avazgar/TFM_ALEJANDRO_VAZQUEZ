###########################################################################################################################
# Script: 02_data_preprocessing.R
# Description:
#       - Import tables and metadata from QIIME2 preprocessing
#       - Data formating
#       - Merge data 
#       - Create phyloseq object for further analysis
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
#########################################################################################################################

# Initial data formatting
# Read taxonomy file and export (.tsv):
tax_table <- read.delim(file.path(raw_path, "taxonomy.tsv"))
# Save tax_table as .xlsx
write.xlsx(tax_table, file.path(processed_path,"tax_table.xlsx"))
# Remove "p__,c__,o__,f__,g__,s__" characters before the taxa name if the taxa file contains them:
taxonomy_clean <- tax_table %>%
separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
	sep = ";\\s*", fill = "right") %>%  # Separate by ';' and space
mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species),
		~str_remove(., "^[a-z]__")))  # Remove prefixes like d__, p__, etc.

# Export the file cleaned:
write.xlsx(taxonomy_clean, file.path(processed_path,"taxonomy_clean.xlsx"), rowNames = FALSE)

# Data loading:
# 1.OTU Table:
otu_table_raw <- read.delim(
file.path(raw_path, "feature-table.tsv"),
skip = 1,
check.names = FALSE
)
otu_table <- otu_table_raw %>%
  column_to_rownames(var = colnames(otu_table_raw)[1])

write.xlsx(otu_table, file.path(raw_path,"otu_table.xlsx"))


# 2.TAXONOMY Table:
tax_table <- readxl::read_excel(file.path(processed_path,"taxonomy_clean.xlsx"))

# 3.METADATA Table:
metadata <- read.delim(file.path(raw_path, "metadata.tsv"))
write.xlsx(metadata, file.path(raw_path,"metadata.xlsx"))

# Preprocess tables format to phyloseq requirements

# ABUNDANCE (OTUs/ASVs): The first column must contain the taxa names (ASVs/OTUs)
# Convert to dataframe and use the first column as rownames:  --> Included as commnent: not always necessary, it depends whether the TAXA IDs appears as rownames or not
# otu_table <- as.data.frame(otu_table)
# rownames(otu_table) <- otu_table[,1] # First column = ASVs IDs
# otu_table <- otu_table[,-1] # Remove the ASVs IDs column (duplicated)

# Convert to phyloseq matrix
otu <- phyloseq::otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)

# TAXONOMY: The first column must contain the ASV/OTU IDs, categories (kingdom, phylum, class..) are displayed along the following columns
tax_table <- as.data.frame(tax_table)
rownames(tax_table) <- tax_table[,1]
tax_table <- tax_table[,-1]

## Convert to phyloseq matrix
tax <- phyloseq::tax_table(as.matrix(tax_table))

# METADATA: The first columns must contain the sampleID
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

## Convert to phyloseq matrix
meta <- phyloseq::sample_data(metadata)

# Phyloseq object (ps) creation from matrix:
ps <- phyloseq(otu, tax, meta)

# Save ps as RDS
saveRDS(ps, file = file.path(processed_path, "boreal_soils.rds"))

