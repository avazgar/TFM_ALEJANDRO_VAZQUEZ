###############################################################################
# Script: 11_kegg_pathway_enrichment_from_kos.R
# Description:
#   - KEGG pathway enrichment using differentially expressed KOs (from functional DEA)
#   - Reads KO-to-pathway mapping from eggNOG-mapper annotations
#   - Performs hypergeometric enrichment test (UP/DOWN separately)
#   - Exports results as TSV tables and barplots (top pathways per contrast)
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Transcriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###############################################################################

source("scripts/00_configuration.R")
source("scripts/01_packages.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# Load eggNOG annotations and extract KO–Pathway pairs
ann_path <- file.path(raw_path, "eggnog_annotation.emapper.annotations")
hdr      <- grep("^#query", readLines(ann_path), value = FALSE)
ann      <- read_tsv(ann_path, skip = hdr - 1, show_col_types = FALSE)
names(ann)[names(ann) == "#query"] <- "query"

	# Unique KO–Pathway mappings
kp <- ann %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "-", !is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(
    KEGG_ko = sub("^ko:", "", trimws(KEGG_ko)),
    KEGG_Pathway = trimws(KEGG_Pathway)
  ) %>%
  filter(KEGG_ko != "", KEGG_Pathway != "") %>%
  distinct(KEGG_ko, KEGG_Pathway)

	# Background universe: all KOs with a pathway assignment
K_all <- unique(kp$KEGG_ko)

# Load KO-level DEA results
dea_KO     <- readRDS(file.path(results_path, "tables_functional", "KO", "DEA_results_KO.rds"))
contrasts  <- names(dea_KO)

# Helper: pathway enrichment with hypergeometric test
path_enrich <- function(K_sig, kp_pairs) {
  K_sig <- base::intersect(K_sig, unique(kp_pairs$KEGG_ko))
  if (length(K_sig) == 0) {
    return(tibble(
      KEGG_Pathway = character(), k_sig = integer(), k_all = integer(),
      expected = double(), enrichment = double(), p_value = double()
    ))
  }
  total_sig <- length(K_sig)
  total_all <- length(unique(kp_pairs$KEGG_ko))
  
  sig_df <- kp_pairs %>% filter(KEGG_ko %in% K_sig) %>% count(KEGG_Pathway, name = "k_sig")
  all_df <- kp_pairs %>% count(KEGG_Pathway, name = "k_all")
  
  sig_df %>%
    right_join(all_df, by = "KEGG_Pathway") %>%
    mutate(k_sig = replace_na(k_sig, 0L)) %>%
    mutate(
      expected   = pmax(k_all / total_all * total_sig, 1e-12),
      enrichment = k_sig / expected,
      p_value    = phyper(k_sig - 1, k_all, total_all - k_all, total_sig, lower.tail = FALSE)
    ) %>%
    arrange(p_value)
}

# Output directories
dir_fig <- file.path(results_path, "figures_functional", "KEGG_enrichment")
dir_tab <- file.path(results_path, "tables_functional", "KEGG_enrichment")
dir.create(dir_fig, TRUE, TRUE); dir.create(dir_tab, TRUE, TRUE)

# Plotting function
plot_path_enrich <- function(df, title, top_n = 15) {
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(title, "(no pathways)")))
  }
  
  top <- df %>%
    slice_min(p_value, n = min(top_n, nrow(df)), with_ties = FALSE) %>%
    mutate(KEGG_Pathway = factor(KEGG_Pathway, levels = rev(KEGG_Pathway)))
  
  ggplot(top, aes(x = KEGG_Pathway, y = enrichment, fill = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "skyblue", high = "red") +
    labs(
      title = title,
      x = "KEGG Pathway (koXXXXX)",
      y = "Enrichment",
      fill = "-log10(p)"
    ) +
    theme_minimal(base_size = 12)
}

# Loop over contrasts × directions
for (ct in contrasts) {
  tt <- dea_KO[[ct]]$full
  K_up   <- tt %>% filter(FDR <= 0.05, logFC >=  1) %>% pull(FeatureID) %>% sub("^ko:", "", .)
  K_down <- tt %>% filter(FDR <= 0.05, logFC <= -1) %>% pull(FeatureID) %>% sub("^ko:", "", .)
  
  for (dirn in c("UP","DOWN")) {
    Ks <- if (dirn == "UP") K_up else K_down
    enr <- path_enrich(Ks, kp)
    
    # Save tables
    write.table(enr, file.path(dir_tab, paste0("KEGG_pathway_", dirn, "_", ct, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Save plots
    p <- plot_path_enrich(enr, paste0("KEGG Pathways — ", dirn, " (", ct, ")"))
    ggsave(file.path(dir_fig, paste0("KEGG_pathway_", dirn, "_", ct, ".png")),
           p, width = 9, height = 6, dpi = 300)
  }
}

message("KEGG pathway enrichment completed.")
###############################################################################
