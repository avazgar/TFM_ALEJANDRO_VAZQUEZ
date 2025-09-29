######################################################################################
# Script: 04_edger_gene_level.R
# Description:
#   - Run pairwise contrasts and DEA at gene level (edgeR QL F-test)
#   - Export full and filtered tables per contrast
#   - MA + Volcano plots (top genes labelled)
#   - Venn diagrams (UP/DOWN) with region exports
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Transcriptomic processing  - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
######################################################################################

# Load base project scripts
#source(file.path(scripts_path, "00_configuration.R"))
#source(file.path(scripts_path, "01_packages.R"))  # edgeR, limma, EnhancedVolcano, ggVennDiagram, ggplot2, dplyr, ggrepel

# Load processed objects (from 03_qc_norm.R)
y       <- readRDS(file.path(processed_path, "dds", "dge_filtered_norm.rds"))       # filtered + TMM-normalized
design  <- readRDS(file.path(processed_path, "dds", "design_matrix.rds"))
fit     <- readRDS(file.path(processed_path, "dds", "glmQLFit_object.rds"))
meta    <- readRDS(file.path(processed_path, "tximport", "metadata_aligned.rds"))

	# Formatting checks
stopifnot(identical(rownames(design), colnames(y)))
stopifnot(identical(rownames(design), rownames(meta)[match(rownames(design), rownames(meta))]))

	# Set factor with explicit order (must match design columns)
GROUP <- factor(meta[colnames(y), "SampleZone"],
                levels = c("Oxidized", "Transition", "Unoxidized"))
stopifnot(identical(colnames(design), levels(GROUP)))

# Define contrasts
	# Columns order in design: "Oxidized","Transition","Unoxidized"
print(colnames(design))

L <- cbind(
  OX_VS_UN = c( 1, 0,-1),
  TR_VS_OX = c(-1, 1, 0),
  TR_VS_UN = c( 0, 1,-1)
)
rownames(L) <- colnames(design)

# Compatibility shim (v.edgeR >= 4.4)
if (!exists("decideTestsDGE")) {
  decideTestsDGE <- function(object, adjust.method = "BH",
                             p.value = 0.05, lfc = 0) {
    tab <- if (is.list(object) && !is.null(object$table)) object$table else object
    if (is.null(tab$PValue) || is.null(tab$logFC)) {
      stop("Object lacks 'PValue' and/or 'logFC'.")
    }
    padj <- p.adjust(tab$PValue, method = adjust.method)
    status <- integer(nrow(tab))
    status[padj <= p.value & tab$logFC >  lfc] <-  1L
    status[padj <= p.value & tab$logFC < -lfc] <- -1L
    status
  }
}

# Define paths and create Input/Output folders
dir.create(results_path,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_path,  showWarnings = FALSE, recursive = TRUE)

dea_results_path  <- file.path(results_path, "DEA_gene_level")
dea_figures_path  <- file.path(figures_path, "DEA_gene_level")
dir.create(dea_results_path, showWarnings = FALSE, recursive = TRUE)
dir.create(dea_figures_path, showWarnings = FALSE, recursive = TRUE)

# Function: run one contrast
run_contrast <- function(fit, contrast_vec, name,
                         fdr_cutoff = 0.05, lfc_cutoff = 1) {
	# QL F-test
  qlf  <- edgeR::glmQLFTest(fit, contrast = contrast_vec)

	# Full table (includes FDR)
  full <- edgeR::topTags(qlf, n = Inf)$table
  full$GeneID <- rownames(full)
  full <- full[, c("GeneID", "logFC", "logCPM", "F", "PValue", "FDR")]
  rownames(full) <- full$GeneID

	# Filtered DEGs (biologically meaningful threshold)
  filt <- subset(full, FDR <= fdr_cutoff & abs(logFC) >= lfc_cutoff)

	# Export tables
  write.table(full,
              file = file.path(dea_results_path, paste0("DEA_full_", name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(filt,
              file = file.path(dea_results_path, paste0("DEA_filtered_", name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

	# MA plot (edgeR base)
  png(file.path(dea_figures_path, paste0("MA_", name, ".png")),
      width = 1200, height = 900, res = 150)
  plotMD(qlf,
         status = decideTestsDGE(qlf),
         main   = paste("MA plot:", name),
         xlab   = "Average logCPM")
  abline(h = c(-lfc_cutoff, lfc_cutoff), col = "grey40", lty = 2)
  dev.off()

	# Volcano (EnhancedVolcano) with top labels
  topGenes <- full |>
    dplyr::arrange(FDR, dplyr::desc(abs(logFC))) |>
    head(10) |>
    dplyr::pull(GeneID)

  vp <- EnhancedVolcano(
    full,
    lab             = rownames(full),
    x               = "logFC",
    y               = "FDR",
    pCutoff         = fdr_cutoff,
    FCcutoff        = lfc_cutoff,
    title           = paste("Volcano:", name),
    subtitle        = NULL,
    selectLab       = topGenes,
    boxedLabels     = TRUE,
    drawConnectors  = TRUE,
    widthConnectors = 0.75,
    colConnectors   = "grey40",
    maxoverlaps     = 45,
    pointSize       = 0.8,
    labSize         = 3.0,
    legendPosition  = "right",
    legendLabSize   = 10,
    legendIconSize  = 3.5,
    colAlpha        = 0.7
  )
  ggplot2::ggsave(
    filename = file.path(dea_figures_path, paste0("volcano_", name, "_enhanced.png")),
    plot = vp, width = 8, height = 6, dpi = 300
  )

  list(test = qlf, full = full, filt = filt, topGenes = topGenes)
}

# Apply run_contrast function: Run all contrasts
results_list <- lapply(colnames(L), function(cn) {
  message("Running contrast: ", cn)
  run_contrast(fit, L[, cn], name = cn)
})
names(results_list) <- colnames(L)

saveRDS(results_list, file = file.path(processed_path, "dds", "DEA_gene_level_results.rds"))

# Export top genes per contrast (CSV)
dir.create(file.path(dea_results_path, "topGenes_tables"), showWarnings = FALSE)
for (cn in names(results_list)) {
  full_df  <- results_list[[cn]]$full
  topGenes <- results_list[[cn]]$topGenes
  top_tab  <- full_df[full_df$GeneID %in% topGenes, c("GeneID","logFC","logCPM","FDR")]
  top_tab  <- top_tab[order(top_tab$FDR, -abs(top_tab$logFC)), ]
  write.csv(top_tab,
            file = file.path(dea_results_path, "topGenes_tables", paste0("topGenes_", cn, ".csv")),
            row.names = FALSE)
}

# Optional: ggplot MA plot with labels (top genes identified in previous volcanos)
plot_MA_topGenes <- function(full_df, topGenes, name, outdir, lfc_cutoff = 1) {
  status <- rep("NotSig", nrow(full_df))
  status[full_df$FDR <= 0.05 & full_df$logFC >=  lfc_cutoff] <- "Up"
  status[full_df$FDR <= 0.05 & full_df$logFC <= -lfc_cutoff] <- "Down"
  full_df$status <- status

  col_map <- c("NotSig" = "black", "Up" = "red", "Down" = "blue")
  top_data <- full_df[full_df$GeneID %in% topGenes, ]

  p <- ggplot(full_df, aes(x = logCPM, y = logFC, color = status)) +
    geom_point(size = 1) +
    scale_color_manual(values = col_map) +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = top_data, aes(label = GeneID),
                             size = 3, box.padding = 0.5, point.padding = 0.3,
                             segment.color = "grey50") +
    labs(title = paste("MA plot:", name), x = "Average logCPM", y = "log-fold-change") +
    theme_minimal(base_size = 14)

  ggsave(file.path(outdir, paste0("MA_", name, "_topGenes.png")),
         plot = p, width = 8, height = 6, dpi = 300)
}

for (cn in names(results_list)) {
  plot_MA_topGenes(results_list[[cn]]$full, results_list[[cn]]$topGenes, cn, dea_figures_path, lfc_cutoff = 1)
}

# Venn (UP/DOWN) + export regions
get_sets <- function(res, fdr = 0.05, lfc = 1) {
  full <- res$full
  up   <- rownames(subset(full, FDR <= fdr & logFC >=  lfc))
  down <- rownames(subset(full, FDR <= fdr & logFC <= -lfc))
  list(up = up, down = down)
}
S_OX_UN <- get_sets(results_list[["OX_VS_UN"]])
S_TR_OX <- get_sets(results_list[["TR_VS_OX"]])
S_TR_UN <- get_sets(results_list[["TR_VS_UN"]])

venn_up <- list(OXvsUN_up = S_OX_UN$up, TRvsOX_up = S_TR_OX$up, TRvsUN_up = S_TR_UN$up)
venn_down <- list(OXvsUN_down = S_OX_UN$down, TRvsOX_down = S_TR_OX$down, TRvsUN_down = S_TR_UN$down)

venn_plot <- function(venn_list, title, high_col = "#E69F00") {
  ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.8, set_size = 3.5) +
    scale_fill_gradient(low = "white", high = high_col, name = "N genes") +
    labs(title = title) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title      = element_text(face = "bold"),
      plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
    )
}

ggsave(file.path(dea_figures_path, "Venn_UP_proportional.png"),
       plot = venn_plot(venn_up, "UP genes: common and exclusive", high_col = "#E69F00"),
       width = 10, height = 6.5, dpi = 300)

ggsave(file.path(dea_figures_path, "Venn_DOWN_proportional.png"),
       plot = venn_plot(venn_down, "DOWN genes: common and exclusive", high_col = "#56B4E9"),
       width = 10, height = 6.5, dpi = 300)

# Export Venn regions & counts
venn_regions_3 <- function(sets_named) {
  stopifnot(length(sets_named) == 3)
  nms <- names(sets_named)
  A <- unique(sets_named[[1]]); B <- unique(sets_named[[2]]); C <- unique(sets_named[[3]])
  A_only <- setdiff(A, union(B, C)); B_only <- setdiff(B, union(A, C)); C_only <- setdiff(C, union(A, B))
  AB_only <- setdiff(intersect(A,B), C); AC_only <- setdiff(intersect(A,C), B); BC_only <- setdiff(intersect(B,C), A)
  ABC <- Reduce(intersect, list(A,B,C))
  setNames(list(sort(A_only), sort(B_only), sort(C_only),
                sort(AB_only), sort(AC_only), sort(BC_only), sort(ABC)),
           c(paste0(nms[1], "_only"), paste0(nms[2], "_only"), paste0(nms[3], "_only"),
             paste0(nms[1], "_", nms[2], "_only"),
             paste0(nms[1], "_", nms[3], "_only"),
             paste0(nms[2], "_", nms[3], "_only"),
             paste0(nms[1], "_", nms[2], "_", nms[3])))
}
export_region_tables_safe <- function(sets_named, prefix) {
  regs <- venn_regions_3(sets_named)
  counts_df <- data.frame(region = names(regs), count = sapply(regs, length))
  write.csv(counts_df, file.path(dea_results_path, paste0(prefix, "_regions_counts.csv")), row.names = FALSE)
  outdir <- file.path(dea_results_path, paste0(prefix, "_regions"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(regs)) {
    fn <- file.path(outdir, paste0(prefix, "_", gsub("[^A-Za-z0-9_]+","_", nm), ".txt"))
    write.table(regs[[nm]], fn, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
save_intersections <- function(sets, filename_prefix) {
  combos <- c(combn(names(sets), 2, simplify = FALSE), combn(names(sets), 3, simplify = FALSE))
  for (cmb in combos) {
    ids <- Reduce(intersect, sets[cmb])
    write.table(ids, file.path(dea_results_path, paste0(filename_prefix, "_", paste(cmb, collapse = "_"), ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
save_intersections(venn_up,   "Venn_UP")
export_region_tables_safe(venn_up,   "Venn_UP")
save_intersections(venn_down, "Venn_DOWN")
export_region_tables_safe(venn_down, "Venn_DOWN")
