###############################################################################
# Script: 35_go_enrichment_dotplot_combined.R
# Description:
#   - Generates dotplots for GO enrichment combining UP and DOWN sets
#   - Facets results by ontology (Biological Process / Molecular Function / Cellular Component)
#   - Reads enrichment tables produced by 34_go_enrichment_annotated.R
#   - Selects top-N terms per ontology (based on p-value across UP/DOWN)
#   - Handles missing files and incomplete columns gracefully
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Transcriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###############################################################################

source("scripts/00_configuration.R")
source("scripts/01_packages.R")

# Input/Output directories
dir_tab <- file.path(results_path, "tables_functional", "GO_enrichment_annotated")
dir_fig <- file.path(results_path, "figures_functional", "GO_enrichment_dotplot_combined")
dir.create(dir_fig, TRUE, TRUE)

# Detect contrasts and ontologies 
files <- list.files(dir_tab, pattern = "^GO_enrichment_.*\\.tsv$", full.names = TRUE)
info <- strcapture("^GO_enrichment_(UP|DOWN)_(.+)_([BMC][PF])\\.tsv$",
                   basename(files),
                   proto = list(direction=character(), contrast=character(), ontology=character()))
contrasts   <- unique(info$contrast)
ontologies  <- c("BP","MF","CC")

# Helper: facet plot for a single contrast
make_facet_plot <- function(ct, top_n = 10) {
  
	# Sub-helper: safe reader
  read_enrich <- function(dirn, ont) {
    f <- file.path(dir_tab, paste0("GO_enrichment_", dirn, "_", ct, "_", ont, ".tsv"))
    if (!file.exists(f)) return(tibble())
    df_tmp <- read_tsv(f, show_col_types = FALSE)
	# Ensure required columns are present
    if (!all(c("GO_Name","p_value","enrichment") %in% colnames(df_tmp))) return(tibble())
    df_tmp %>% mutate(Direction = dirn, Ontology = ont)
  }
  
	# Load UP/DOWN for all ontologies
  df <- bind_rows(
    read_enrich("UP","BP"),   read_enrich("DOWN","BP"),
    read_enrich("UP","MF"),   read_enrich("DOWN","MF"),
    read_enrich("UP","CC"),   read_enrich("DOWN","CC")
  )
  
	# If empty: skip
  if (nrow(df) == 0) return(NULL)
  
	# Prepare dataframe
  df <- df %>%
    filter(!is.na(GO_Name), GO_Name != "") %>%
    mutate(logP = -log10(p_value))
  
	# Select top-N per ontology (across UP/DOWN)
  top_terms <- df %>%
    group_by(Ontology) %>%
    slice_min(order_by = p_value, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(Ontology, GO_Name)
  
  df <- df %>% semi_join(top_terms, by = c("Ontology","GO_Name"))
  
	# Ordering terms inside each ontology by global p-value
  term_order <- df %>%
    group_by(Ontology, GO_Name) %>%
    summarise(min_p = min(p_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(Ontology, min_p)
  
  df <- df %>%
    left_join(term_order, by = c("Ontology","GO_Name")) %>%
    group_by(Ontology) %>%
    mutate(GO_Name_ord = factor(GO_Name, levels = rev(unique(GO_Name[order(min_p)])))) %>%
    ungroup()
  
	# Dotplot
  p <- ggplot(df, aes(x = GO_Name_ord, y = enrichment, size = logP, color = Direction)) +
    geom_point(alpha = 0.85) +
    coord_flip() +
    facet_wrap(~ Ontology, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("UP"="#E69F00", "DOWN"="#56B4E9")) +
    labs(
      title = paste0("GO enrichment — ", ct, " (faceted by ontology)"),
      x = "GO term", y = "Enrichment", size = "-log10(p)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  ggsave(file.path(dir_fig, paste0("GO_dotplot_facet_", ct, ".png")),
         p, width = 12, height = 6, dpi = 300)
  p
}

# Run for all contrasts
for (ct in contrasts) {
  make_facet_plot(ct, top_n = 10)
}

message("Facet dotplots per ontology saved in: ", dir_fig)
###############################################################################
