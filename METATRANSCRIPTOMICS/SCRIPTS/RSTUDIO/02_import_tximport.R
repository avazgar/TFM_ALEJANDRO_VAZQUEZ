######################################################################################
# Script: 02_import_tximport.R
# Description:
#   - Import counts and TPM (Salmon) from .tsv exports
#   - Apply initial QC and normalization
#   - Save processed objects for downstream DEA at gene level
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
######################################################################################

# Load base project scripts 
#source("scripts/00_configuration.R") 
#source("scripts/01_packages.R") 

# Read quantification files
counts <- read.delim(file.path(raw_path, "salmon_counts.tsv"),
                     row.names = 1, check.names = FALSE)
tpm <- read.delim(file.path(raw_path, "salmon_TPM.tsv"),
                  row.names = 1, check.names = FALSE)

	# Ensure identical sample order
stopifnot(identical(colnames(counts), colnames(tpm)))

# Metadata alignment
meta <- read.delim(file.path(raw_path, "metadata.tsv"),
                   sep = "\t", row.names = 1, check.names = FALSE)

	# Map SRR IDs to final names
names_final <- setNames(rownames(meta), meta$Original_SampleID)
colnames(counts) <- names_final[colnames(counts)]
colnames(tpm) <- names_final[colnames(tpm)]

	# Align metadata to counts order
meta <- meta[match(colnames(counts), rownames(meta)), , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(meta)))

# Filter genes without expression
counts <- counts[rowSums(counts) > 0, ]

# QC: Histogram of aveLogCPM
counts <- counts[complete.cases(counts), ]
tpm    <- tpm[complete.cases(tpm), ]  # corrected line

ave_cpm <- aveLogCPM(counts)
hist(ave_cpm,
     main = "Genes after filtering",
     xlab = "Average logCPM",
     ylab = "Number of genes",
     xlim = c(-5, 15), breaks = 30,
     col = "lightblue", border = "white")

# Create DGEList 
y <- DGEList(counts = counts, group = meta$SampleZone)
y$samples

# Library sizes
par(mar = c(8, 5, 3, 1))  
heights_millones <- y$samples$lib.size / 1e6

barplot(
  heights_millones,
  names.arg = colnames(y),
  las = 2,
  ylim = c(0, max(heights_millones) * 1.1),
  cex.names = 0.9,
  main = "Library sizes",
  ylab = "Millions of counts"
)
mtext(side = 1, text = "Samples", line = 5)

# Normalization (TMM)
y <- calcNormFactors(y, method = "TMM")

# QC: Boxplot of normalized logCPM
cpms <- edgeR::cpm(y, log = TRUE)
boxplot(cpms,
        xlab = "Samples", ylab = "Log2 CPM", las = 2,
        main = "Distribution of logCPM after TMM normalization")
abline(h = median(cpms), col = "red")

# QC: MDS plot
plotMDS(y, main = "MDS plot - Biological distances")

# Save processed objects
dir.create(file.path(processed_path, "tximport"), recursive = TRUE, showWarnings = FALSE)
saveRDS(counts, file.path(processed_path, "tximport", "salmon_counts_filtered.rds"))
saveRDS(tpm,    file.path(processed_path, "tximport", "salmon_TPM_raw.rds"))
saveRDS(y,      file.path(processed_path, "tximport", "DGEList_gene_level.rds"))
saveRDS(meta,   file.path(processed_path, "tximport", "metadata_aligned.rds"))
