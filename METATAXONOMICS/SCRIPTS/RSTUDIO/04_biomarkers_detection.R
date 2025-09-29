###########################################################################################################################
# Script: 04_biomarkers_detection.R
# Description:
#	  - Identification of differential biomarkers through deseq2
#	  - Parwise comparisons between different conditions
#	  - ASV - taxonomy integration
#	  - Volcano plot visualizations
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing IRTA [BIOSISTEMAS] - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###########################################################################################################################

# Load environment variables and packages if needed [added as comments]:
# source("scripts/00_configuration.R")
# source(file.path(scripts_path, "01_packages.R")

# Prepare data
	# Filter valid samples (no NA values) from phyloseq object and define a new deseq phyloseq object
ps_deseq <- subset_samples(ps, !is.na(SampleZone))

	# Extract couns and metadata
otu_counts <- as(otu_table(ps_deseq), "matrix")
if (!taxa_are_rows(ps_deseq)) otu_counts <- t(otu_counts)
sample_info <- as(sample_data(ps_deseq), "data.frame")

	# Extract taxonomy table
tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)
tax[is.na(tax)] <- "Unassigned"
tax$taxon <- paste(tax$Phylum, tax$Genus, sep = ":")

# Create Deseq2 object
	# Create an auxiliary function to build and run Deseq2 with customised reference
dds <- DESeqDataSetFromMatrix(countData = otu_counts,
                              colData = sample_info,
                              design = ~ SampleZone)
dds <- dds[rowSums(counts(dds)) > 50, ]   # filtrar baja abundancia
dds <- DESeq(dds)

# Make parwise comparisons
comparisons <- list(
  Oxidized_vs_Unoxidized   = c("SampleZone", "Oxidized", "Unoxidized"),
  Transition_vs_Unoxidized = c("SampleZone", "Transition", "Unoxidized"),
  Transition_vs_Oxidized   = c("SampleZone", "Transition", "Oxidized")
)

results_list <- lapply(names(comparisons), function(label) {
  res <- results(dds, contrast = comparisons[[label]])
  df <- as.data.frame(res)
  df$ASV <- rownames(df)
  df$comparison <- label
  df <- merge(df, tax[, c("ASV", "taxon")], by = "ASV", all.x = TRUE)
  return(df)
})
names(results_list) <- names(comparisons)

# Save combined results

res_all <- do.call(rbind, results_list)
write.csv(res_all,
          file = file.path(processed_path,"DESeq2_all_pairwise_comparisons.csv"),
          row.names = TRUE)

# Function Volcano Plot for visualization of significant biomarkers
volcano_plot <- function(df, comparison_label, top_n = 5) {
  df$significance <- "Not significant"
  df$significance[df$padj < 0.05 & abs(df$log2FoldChange) > 1] <- "padj < 0.05 & |log2FC| > 1"
  df$significance[df$padj < 0.05 & abs(df$log2FoldChange) <= 1] <- "padj < 0.05 only"
  df$significance[df$padj >= 0.05 & abs(df$log2FoldChange) > 1] <- "|log2FC| > 1 only"
  df$significance <- factor(df$significance,
                            levels = c("padj < 0.05 & |log2FC| > 1",
                                       "padj < 0.05 only",
                                       "|log2FC| > 1 only",
                                       "Not significant"))
  
  top_fc   <- df[order(-abs(df$log2FoldChange)), ][1:top_n, ]
  top_padj <- df[order(df$padj), ][1:top_n, ]
  labels_df <- unique(rbind(top_fc, top_padj))
  
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c(
      "padj < 0.05 & |log2FC| > 1" = "red",
      "padj < 0.05 only" = "blue",
      "|log2FC| > 1 only" = "orange",
      "Not significant" = "grey"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(label = taxon),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Volcano plot -", comparison_label),
         x = "log2 Fold Change", y = "-log10(padj)", color = "Significance") +
    theme(legend.position = "bottom")
  
  return(p)
}

# Generate and save volcanos
for (label in names(results_list)) {
  p <- volcano_plot(results_list[[label]], label)
  ggsave(file.path(figures_path,"volcano", paste0(label, ".png")),
         plot = p, width = 8, height = 6, dpi = 300)
}

