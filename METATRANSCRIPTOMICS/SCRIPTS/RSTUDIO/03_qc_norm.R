######################################################################################
# Script: 03_qc_norm.R
# Description:
#   - Filter genes by expression (filterByExpr)
#   - Define experimental design
#   - Estimate common/trended/tagwise dispersions
#   - Fit NB-GLM quasi-likelihood model
#   - Save processed objects for contrasts and DEA
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Context: Transcriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
######################################################################################

# Load base project scripts
#source(file.path(scripts_path, "00_configuration.R"))
#source(file.path(scripts_path, "01_packages.R"))

# Load objects from previous step
y    <- readRDS(file.path(processed_path, "tximport", "DGEList_gene_level.rds"))
meta <- readRDS(file.path(processed_path, "tximport", "metadata_aligned.rds"))

# Define experimental design
stopifnot("SampleZone" %in% colnames(meta))
GROUP <- factor(meta$SampleZone)
GROUP <- droplevels(GROUP)
y$samples$group <- GROUP

	# Check replicates per condition
print(table(GROUP))

# Filter genes by expression
keep <- filterByExpr(y, group = GROUP, min.count = 10)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize data (TMM)
y <- calcNormFactors(y, method = "TMM")

# Create design matrix
design <- model.matrix(~0 + GROUP)
colnames(design) <- levels(GROUP)
rownames(design) <- rownames(meta)

# Estimate dispersions
y <- estimateDisp(y, design, robust = TRUE)

	# QC: Biological Coefficient of Variation (BCV)
par(mar = c(4,4,2,1))
plotBCV(y, main = "BCV after filtering and TMM normalization")

# Fit NB-GLM (quasi-likelihood)
fit <- glmQLFit(y, design, robust = TRUE)

	# QC: Quasi-likelihood dispersion vs average logCPM
par(mar = c(4,4,2,1))
plotQLDisp(fit, ylab = "Quasi-likelihood dispersion",
           main = "QL dispersion vs AveLogCPM")

# Save objects for contrasts
dir.create(file.path(processed_path, "dds"), recursive = TRUE, showWarnings = FALSE)
saveRDS(y,      file = file.path(processed_path, "dds", "dge_filtered_norm.rds"))
saveRDS(design, file = file.path(processed_path, "dds", "design_matrix.rds"))
saveRDS(fit,    file = file.path(processed_path, "dds", "glmQLFit_object.rds"))
