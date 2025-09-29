######################################################################################
# Script: 06_edger_functional.R
# Description:
#   - Differential Expression Analysis (DEA) at functional level (KO / COG / GO)
#   - Based on functional matrices built from eggNOG-mapper annotations
#   - Workflow:
#       1) Load functional count matrices and metadata
#       2) Filter features with low expression
#       3) Normalize libraries (TMM)
#       4) Estimate dispersion and fit NB-GLM (QL framework)
#       5) Run contrasts between soil horizons (Oxidized, Transition, Unoxidized)
#       6) Export results (full + filtered tables)
#       7) Generate MA plots and Volcano plots for each contrast
# Date: 23/09/2025
# Ownerr: Alejandro Vázquez García
# Context: Transcriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
######################################################################################

# Load environment configuration and packages

source("scripts/00_configuration.R")
source("scripts/01_packages.R")

# Load functional matrix and metadata

KO   <- readRDS(file.path(processed_path, "matrices", "KO_counts.rds"))
COG  <- readRDS(file.path(processed_path, "matrices", "COG_counts.rds"))
GO   <- readRDS(file.path(processed_path, "matrices", "GO_counts.rds"))
meta <- readRDS(file.path(processed_path, "tximport", "metadata_aligned.rds"))

stopifnot(identical(rownames(meta), colnames(KO)))

# Define experimental design

GROUP <- factor(meta$SampleZone, levels = c("Oxidized","Transition","Unoxidized"))
design <- model.matrix(~0 + GROUP)
colnames(design) <- levels(GROUP)
rownames(design) <- rownames(meta)

# Create functional DEA function

run_dea_fun <- function(M, label, fdr_cutoff = 0.05, lfc_cutoff = 1) {
  
  # Filtering features with no expression
  M <- M[rowSums(M) > 0, , drop = FALSE]
  
  # Create edgeR object
  y <- DGEList(counts = M)
  keep <- filterByExpr(y$counts, group = GROUP, min.count = 10)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  
  # Estimate dispersion and fit the model
  y   <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  
  # Define contrasts
  L <- cbind(
    OX_VS_UN = c( 1,  0, -1),
    TR_VS_OX = c(-1,  1,  0),
    TR_VS_UN = c( 0,  1, -1)
  )
  rownames(L) <- colnames(design)
  
  # Create output directories 
  outT <- file.path(results_path, "tables_functional", label); dir.create(outT, TRUE, TRUE)
  outF <- file.path(results_path, "figures_functional", label); dir.create(outF, TRUE, TRUE)
  
  # Run contrast 
  results_list <- list()
  for (cn in colnames(L)) {
    qlf  <- glmQLFTest(fit, contrast = L[, cn])
    full <- topTags(qlf, n = Inf)$table
    full$FeatureID <- rownames(full)
    full <- full[, c("FeatureID","logFC","logCPM","F","PValue","FDR")]
    
    filt <- subset(full, FDR <= fdr_cutoff & abs(logFC) >= lfc_cutoff)
    
    # Save tables
    write.table(full, file.path(outT, paste0("DEA_full_", cn, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(filt, file.path(outT, paste0("DEA_filtered_", cn, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # MA Plot
    png(file.path(outF, paste0("MA_", cn, ".png")), 1200, 900, res = 150)
    plotMD(qlf, status = decideTestsDGE(qlf),
           main = paste(label, "MA:", cn),
           xlab = "Average logCPM")
    abline(h = c(-lfc_cutoff, lfc_cutoff), col = "grey40", lty = 2)
    dev.off()
    
    # Volcano plot
    rownames(full) <- full$FeatureID
    yl <- max(-log10(full$FDR), na.rm = TRUE)
    p <- EnhancedVolcano(
      full, lab = rownames(full),
      x = "logFC", y = "FDR",
      pCutoff = fdr_cutoff, FCcutoff = lfc_cutoff,
      title = paste(label, ":", cn), subtitle = NULL,
      ylim = c(0, yl + 0.5),
      labSize = 3, pointSize = 1.5
    )
    ggsave(file.path(outF, paste0("Volcano_", cn, "_enhanced.png")),
           p, width = 8, height = 6, dpi = 300)
    
    results_list[[cn]] <- list(test = qlf, full = full, filtered = filt)
  }
  
  # Save results
  saveRDS(results_list, file.path(outT, paste0("DEA_results_", label, ".rds")))
  invisible(results_list)
}

# Run functional analysis

res_KO  <- run_dea_fun(KO,  "KO")
# res_COG <- run_dea_fun(COG, "COG")
# res_GO  <- run_dea_fun(GO,  "GO")

message("Functional DEA completed (KO).")
