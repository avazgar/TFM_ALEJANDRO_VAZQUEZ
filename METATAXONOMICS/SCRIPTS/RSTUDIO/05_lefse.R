###########################################################################################################################
# Script: 04_biomarkers_detection.R
# Description:
#         - Identification of differential biomarkers through LEfse (lefser package in R)
#         - Application of LEfse in each pairwise comparison
#         - Cleaning and filtering of the results: Minimum taxonomic level --> Class
#         - LDA Barplots generation (png & pdf)
# 	  - Export results in csv format
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing IRTA [BIOSISTEMAS] - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###########################################################################################################################

## ------------------------------------------------------------
## 1) Prepare data: phyloseq -> SummarizedExperiment (relative abundances)
## ------------------------------------------------------------
ps_relab <- transform_sample_counts(ps, function(x) x / sum(x))

otu_rel <- as(otu_table(ps_relab), "matrix")
if (!taxa_are_rows(ps_relab)) otu_rel <- t(otu_rel)

tax <- as.data.frame(tax_table(ps_relab))
rownames(otu_rel) <- apply(tax, 1, function(x) paste(na.omit(x), collapse = "|"))

meta <- as.data.frame(sample_data(ps_relab))

se <- SummarizedExperiment(
  assays  = list(counts = otu_rel),
  colData = meta
)

## ------------------------------------------------------------
## 2) Function: run LEfSe + clean + select top-N
## ------------------------------------------------------------
run_lefse <- function(se_obj, groups,
                      lda_thr = 2, p_thr = 0.05, top_n = 30,
                      class_col = "SampleZone") {
  
  # Subset samples and set class levels/order
  keep <- colData(se_obj)[[class_col]] %in% groups
  se_sub <- se_obj[, keep]
  cls <- droplevels(factor(colData(se_sub)[[class_col]], levels = groups))
  
  # Relative matrix (rows = features, cols = samples)
  mat <- assay(se_sub, "counts")
  mat <- sweep(mat, 2, colSums(mat), "/")
  mat[is.na(mat)] <- 0
  
  # Run LEfSe (matrix + class vector)
  res <- as.data.frame(lefser::lefse(mat, class = cls, lda.threshold = lda_thr))
  
  # -----------------------
  # Normalize column names
  # -----------------------
  if ("features" %in% names(res) && !"feature" %in% names(res)) res <- dplyr::rename(res, feature = features)
  if ("scores"   %in% names(res)) res <- dplyr::rename(res, lda = scores)
  if ("class"    %in% names(res)) res <- dplyr::rename(res, group = class)
  if ("pval"     %in% names(res) && !"pvalue" %in% names(res)) res <- dplyr::rename(res, pvalue = pval)
  
  # -----------------------
  # Filter by minimum taxonomic depth (≥ 3 levels)
  # -----------------------
  feature_std <- gsub(";", "|", as.character(res$feature))
  parts      <- strsplit(feature_std, "\\|")
  tax_depth  <- lengths(parts)
  third_tok  <- vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, "")
  remove_rows <- (tax_depth < 3) | is.na(third_tok) | !grepl("[A-Za-z]", third_tok)
  res <- res[!remove_rows, , drop = FALSE]
  
  # -----------------------
  # Labels
  # -----------------------
  if (!"ASV_label" %in% names(res)) res$ASV_label <- paste0("ASV", seq_len(nrow(res)))
  res$Phylum <- ifelse(grepl("\\|", res$feature),
                       sapply(strsplit(as.character(res$feature), "\\|"), `[`, 1),
                       NA_character_)
  res$feature_label <- paste0(res$ASV_label,
                              ifelse(!is.na(res$Phylum), paste0(" | ", res$Phylum), ""),
                              " | ", gsub(";", "·", gsub("_", " ", res$feature)))
  
  # -----------------------
  # Filter by LDA and p-value (if available)
  # -----------------------
  if ("pvalue" %in% names(res)) {
    res_f <- dplyr::filter(res, abs(lda) >= lda_thr, pvalue <= p_thr)
  } else {
    res_f <- dplyr::filter(res, abs(lda) >= lda_thr)
  }
  
  # Return early if no biomarkers
  if (nrow(res_f) == 0) {
    return(tibble::tibble(feature_label = character(), group = character(),
                          lda = numeric(), lda_signed = numeric()))
  }
  
  # -----------------------
  # Top-N per group (balanced)
  # -----------------------
  res_top <- res_f %>%
    dplyr::group_by(group) %>%
    dplyr::slice_max(order_by = abs(lda), n = ceiling(top_n / 2), with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Signed LDA (left side = second group)
  res_top <- res_top %>%
    dplyr::mutate(
      lda_signed = ifelse(group == groups[1], abs(lda), -abs(lda))
    ) %>%
    dplyr::arrange(lda_signed) %>%
    dplyr::mutate(feature_label = forcats::fct_inorder(feature_label))
  
  return(res_top)
}

## ------------------------------------------------------------
## 3) Plot function
## ------------------------------------------------------------
plot_lefse <- function(res_top, comp_label) {
  pal <- c(Oxidized = "#E67E22", Unoxidized = "#1F77B4", Transition = "#2ECC71")
  if (nrow(res_top) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste("LEfSe –", comp_label, "(no biomarkers at threshold)")))
  }
  lim <- max(abs(res_top$lda_signed)); lim <- ceiling(lim * 1.05)
  
  ggplot(res_top, aes(x = lda_signed, y = feature_label, fill = group)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
    scale_x_continuous(limits = c(-lim, lim), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = pal, name = "Group") +
    labs(title = paste("LEfSe –", comp_label),
         x = "LDA score (log10)", y = NULL) +
    theme_classic(base_size = 16) +
    theme(
      plot.title   = element_text(size = 20, face = "bold"),
      axis.text.y  = element_text(size = 9),
      axis.text.x  = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      plot.margin  = margin(t = 10, r = 20, b = 10, l = 300)
    ) +
    coord_cartesian(clip = "off")
}

## ------------------------------------------------------------
## 4) Comparisons and export
## ------------------------------------------------------------
comparisons <- list(
  "Oxidized_vs_Unoxidized"   = c("Oxidized", "Unoxidized"),
  "Transition_vs_Unoxidized" = c("Transition", "Unoxidized"),
  "Transition_vs_Oxidized"   = c("Transition", "Oxidized")
)

for (label in names(comparisons)) {
  res_top <- run_lefse(se, groups = comparisons[[label]], lda_thr = 3, p_thr = 0.05, top_n = 40)
  p <- plot_lefse(res_top, label)
  
  n_bars <- nrow(res_top)
  height_in <- max(9, n_bars * 0.35)
  
  ggsave(file.path(figures_path, paste0("lefse_biomarkers_barplot_", label, ".png")),
         plot = p, width = 20, height = height_in, units = "in", dpi = 300, limitsize = FALSE)
  ggsave(file.path(figures_path, paste0("lefse_biomarkers_barplot_", label, ".pdf")),
         plot = p, width = 20, height = height_in, units = "in", device = cairo_pdf, limitsize = FALSE)
  
  write.csv(res_top,
            file = file.path(processed_path, paste0("lefse_results_", label, ".csv")),
            row.names = FALSE)
}
