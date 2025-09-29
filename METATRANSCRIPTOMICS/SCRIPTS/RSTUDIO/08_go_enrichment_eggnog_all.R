###############################################################################
# Script: 08_go_enrichment_eggnog_all.R
# Description:
#   - GO (Gene Onthology) enrichment analysis from eggNOG annotations
#   - Use of differential genes (UP/DOWN) of the defined contrasts at gene level
#   - Hypergeometric test (ORA, Over-Representation Analysis) -> DEGs vs GO Universe
#   - Application of FDR correction (Benjamini–Hochberg).
#   - Export tables, individual graphs and summary
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Transcriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###############################################################################

# Load configuration and packages
source("scripts/00_configuration.R")
source("scripts/01_packages.R")

set.seed(123) # reproducibility

# Load eggNOG annotations
ann_path <- file.path(raw_path, "eggnog_annotation.emapper.annotations")
hdr      <- grep("^#query", readLines(ann_path), value = FALSE)
ann      <- read_tsv(ann_path, skip = hdr - 1, show_col_types = FALSE)
names(ann)[names(ann) == "#query"] <- "query"

stopifnot(all(c("query","GOs") %in% names(ann)))

# GOs universe in the dataset
all_gos <- ann %>%
  filter(!is.na(GOs), GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  mutate(GOs = trimws(GOs)) %>%
  filter(GOs != "")

# Load DEA results (DDS) at gene level
dea_gene <- readRDS(file.path(processed_path, "dds", "DEA_gene_level_results.rds"))
contrastes <- names(dea_gene)  # p.ej. c("OX_VS_UN","TR_VS_OX","TR_VS_UN")

# Auxiliary functions
	# Extract GOs associated to any set of genes
extract_gos <- function(gene_ids, ann_df) {
  ann_df %>%
    filter(query %in% gene_ids, !is.na(GOs), GOs != "-") %>%
    separate_rows(GOs, sep = ",") %>%
    mutate(GOs = trimws(GOs)) %>%
    filter(GOs != "") %>%
    count(GOs, sort = TRUE, name = "count")
}

	# Hypergeometric enrichment with FDR adjusting
enrich_calc <- function(sig_gos, all_gos_df) {
  if (nrow(sig_gos) == 0) {
    return(tibble(GOs=character(), count=integer(), all_count=integer(),
                  expected=double(), enrichment=double(), 
                  p_value=double(), padj=double()))
  }
  total_sig <- sum(sig_gos$count)
  total_all <- nrow(all_gos_df)
  sig_gos %>%
    left_join(count(all_gos_df, GOs, name = "all_count"), by = "GOs") %>%
    mutate(
      expected   = all_count / total_all * total_sig,
      enrichment = count / pmax(expected, 1e-12),
      p_value    = phyper(count - 1, all_count, total_all - all_count, total_sig,
                          lower.tail = FALSE),
      padj       = p.adjust(p_value, method = "BH")
    ) %>%
    arrange(padj)
}

	# Enrichment Plot (top-N by FDR)
plot_enrich <- function(enrich_df, title, top_n = 10) {
  if (nrow(enrich_df) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(title, "(sin términos)")))
  }
  top <- enrich_df %>%
    slice_min(padj, n = min(top_n, nrow(enrich_df))) %>%
    mutate(GOs = factor(GOs, levels = rev(GOs)))
  ggplot(top, aes(x = GOs, y = enrichment, fill = -log10(padj))) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "skyblue", high = "red") +
    labs(title = title, x = "GO ID", y = "Enrichment", fill = "-log10(FDR)") +
    theme_minimal(base_size = 12)
}

# Output directories
dir_fig <- file.path(results_path, "figures_functional", "GO_enrichment")
dir_tab <- file.path(results_path, "tables_functional",  "GO_enrichment")
dir.create(dir_fig, TRUE, TRUE); dir.create(dir_tab, TRUE, TRUE)

# Main loop (contrast by direction)
summary_list <- list()

for (ct in contrastes) {
  tt <- dea_gene[[ct]]$full
  up_genes   <- rownames(tt[tt$FDR <= 0.05 & tt$logFC >=  1, , drop = FALSE])
  down_genes <- rownames(tt[tt$FDR <= 0.05 & tt$logFC <= -1, , drop = FALSE])
  
  for (dirn in c("UP","DOWN")) {
    gene_set <- if (dirn == "UP") up_genes else down_genes
    
    sig_gos <- extract_gos(gene_set, ann)
    enrich  <- enrich_calc(sig_gos, all_gos)
    
    # Save table
    out_tsv <- file.path(dir_tab, paste0("GO_enrich_", dirn, "_", ct, ".tsv"))
    write.table(enrich, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Save graph
    p <- plot_enrich(enrich, paste("Top Enriched GO —", dirn, "(", ct, ")"))
    ggsave(file.path(dir_fig, paste0("GO_enrich_", dirn, "_", ct, ".png")),
           p, width = 8, height = 6, dpi = 300)
    
    # resumen global
    summary_list[[paste(ct, dirn, sep="_")]] <- tibble(
      contrast  = ct,
      direction = dirn,
      n_sigGO   = sum(enrich$padj <= 0.05, na.rm = TRUE)
    )
  }
}

# Save global summary
summary_tab <- bind_rows(summary_list)
write.csv(summary_tab, file.path(dir_tab, "GO_enrichment_summary.csv"), row.names = FALSE)

message("GO enrichment completed for all the contrasts (UP/DOWN).")
###############################################################################
