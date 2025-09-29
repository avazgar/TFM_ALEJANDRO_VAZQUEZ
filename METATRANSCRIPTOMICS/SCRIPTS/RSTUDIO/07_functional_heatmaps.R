######################################################################################
# Script: 07_functional_heatmaps.R
# Description:
#   - Visualization of differential functions at KO level (eggNOG-mapper)
#   - Main objetives of the script:
#       1) Select top-N KO significative functions by FDR and |logFC|
#       2) Create expression heatmaps (logCPM, z-score by row) for visualizing patterns
#       3) Build clustering heatmaps of samples based on differential KO functions
#       4) Relate correlations between samples (Pearson/Spearman) with selected KOs
#   - Notes:
#       · LogCPM normalization is applied to stabilize variances
#       · Rows are scaled to z-scores to facilitate making relationships between functions
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Contexto: Transcriptomic processing  - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
######################################################################################

# Load environment configuration and packages 
source("scripts/00_configuration.R")
source("scripts/01_packages.R")

# Load matrix and metadata

KO   <- readRDS(file.path(processed_path, "matrices", "KO_counts.rds"))
meta <- readRDS(file.path(processed_path, "tximport", "metadata_aligned.rds"))
stopifnot(identical(rownames(meta), colnames(KO)))

# Function selection (top-N by contrast)

dea <- readRDS(file.path(results_path, "tables_functional", "KO", "DEA_results_KO.rds"))
tt  <- dea[["OX_VS_UN"]]$full

N <- 500
top_ids <- tt %>%
  arrange(FDR, desc(abs(logFC))) %>%
  slice_head(n = N) %>%
  pull(FeatureID)

mat <- KO[top_ids, , drop = FALSE]

	# Transform to logCPM and scale by row (z-score)
logcpm <- cpm(mat, log = TRUE, prior.count = 1)
zsc    <- t(scale(t(logcpm)))

ann_col <- data.frame(SampleZone = meta$SampleZone)
rownames(ann_col) <- rownames(meta)

	# Heatmap top-N KO
dir.create(file.path(results_path, "figures_functional", "KO"), TRUE, TRUE)
pheatmap(
  zsc,
  annotation_col = ann_col,
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "KO top (FDR) — z-score de logCPM",
  filename = file.path(results_path, "figures_functional", "KO", "Heatmap_KO_topFDR.png"),
  width = 10, height = 10
)
write.table(rownames(zsc),
            file.path(results_path, "tables_functional", "KO", "KO_ids_in_heatmap.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Heatmap KO top-N created.")

# Union significative DE 

sig_union <- unique(unlist(lapply(dea, \(x) {
  d <- x$full
  rownames(d[d$FDR <= 0.05 & abs(d$logFC) >= 1, , drop = FALSE])
})))

KO_sel <- KO[intersect(sig_union, rownames(KO)), , drop = FALSE]
if (nrow(KO_sel) < 5) {
  # fallback:  top 200 OX_VS_UN
  tt <- dea[["OX_VS_UN"]]$full |> arrange(FDR, desc(abs(logFC)))
  KO_sel <- KO[head(tt$FeatureID, 200), , drop = FALSE]
}

# logCPM + z-score
logcpm <- cpm(KO_sel, log = TRUE, prior.count = 1)
zsc    <- t(scale(t(logcpm)))

ann_col <- data.frame(SampleZone = meta$SampleZone)
rownames(ann_col) <- rownames(meta)

# Samples clustering heatmap 
pheatmap(zsc,
         annotation_col = ann_col,
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",   # opción: "correlation"
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Clustering of samples by KO profiles (DE functions)",
         filename = file.path(results_path, "figures_functional", "KO",
                              "Heatmap_samples_KO_DE_functions.png"),
         width = 10, height = 9)

# Heatmap of correlation between samples

C <- cor(logcpm, method = "pearson")  # opción: "spearman"
ann <- data.frame(SampleZone = meta$SampleZone)
rownames(ann) <- rownames(meta)

pheatmap(C,
         annotation_col = ann, annotation_row = ann,
         clustering_method = "complete",
         main = "Correlation between samples (KO, logCPM)",
         filename = file.path(results_path, "figures_functional", "KO",
                              "Heatmap_samples_KO_correlation.png"),
         width = 8, height = 7)

message("Functional heatmaps completed")
