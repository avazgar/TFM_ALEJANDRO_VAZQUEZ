###############################################################################
# Script: 09_go_enrichment_annotated.R
# Description:
#   - Annotates enriched GO IDs with Names + Ontology (BP/MF/CC) using GO.db
#   - Recalculates enrichment separately for Biological Process (BP),
#     Molecular Function (MF), and Cellular Component (CC)
#   - Performs over-representation analysis (hypergeometric test)
#   - Exports tables and barplots for each contrast and direction (UP/DOWN)
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Transcriptomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###############################################################################

source("scripts/00_configuration.R")
source("scripts/01_packages.R")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(AnnotationDbi); library(GO.db)  # required for GO term annotation
})
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

# Load eggNOG annotations (GO universe)
ann_path <- file.path(raw_path, "eggnog_annotation.emapper.annotations")
hdr      <- grep("^#query", readLines(ann_path), value = FALSE)
ann      <- read_tsv(ann_path, skip = hdr - 1, show_col_types = FALSE)
names(ann)[names(ann) == "#query"] <- "query"

stopifnot(all(c("query","GOs") %in% names(ann)))

	# Universe of GOs observed in the dataset
all_gos <- ann %>%
  filter(!is.na(GOs), GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  mutate(GOs = trimws(GOs)) %>%
  filter(GOs != "")

# Build GO dictionary (ID: Name + Ontology)
map_go <- function(go_ids) {
  go_ids <- unique(trimws(go_ids))
  go_ids[!grepl("^GO:", go_ids)] <- paste0("GO:", go_ids[!grepl("^GO:", go_ids)])
  
  res <- AnnotationDbi::select(
    GO.db,
    keys    = go_ids,
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  )
  
  res <- unique(res[!is.na(res$GOID), c("GOID","TERM","ONTOLOGY")])
  tibble::tibble(
    GOs      = res$GOID,
    GO_Name  = res$TERM,
    Ontology = res$ONTOLOGY
  )
}

go_dict <- map_go(all_gos$GOs)

	# Add ontology to the universe
all_gos_onto <- dplyr::left_join(all_gos, go_dict, by = "GOs") %>%
  dplyr::filter(!is.na(Ontology))

# Load DEG results (gene-level)
dea_gene <- readRDS(file.path(processed_path, "dds", "DEA_gene_level_results.rds"))
contrastes <- names(dea_gene)  # e.g. c("OX_VS_UN","TR_VS_OX","TR_VS_UN")

# Helper functions
	# Extract GO counts from gene IDs
extract_gos <- function(gene_ids, ann_df) {
  ann_df %>%
    filter(query %in% gene_ids, !is.na(GOs), GOs != "-") %>%
    separate_rows(GOs, sep = ",") %>%
    mutate(GOs = trimws(GOs)) %>%
    filter(GOs != "") %>%
    count(GOs, sort = TRUE, name = "count")
}

	# Hypergeometric test for enrichment
enrich_calc <- function(sig_gos_counts, universe_gos_df) {
  if (nrow(sig_gos_counts) == 0 || nrow(universe_gos_df) == 0) {
    return(tibble(GOs=character(), count=integer(), all_count=integer(),
                  expected=double(), enrichment=double(), p_value=double()))
  }
  total_sig <- sum(sig_gos_counts$count)
  total_all <- nrow(universe_gos_df)
  sig_gos_counts %>%
    left_join(count(universe_gos_df, GOs, name = "all_count"), by = "GOs") %>%
    mutate(
      expected   = pmax(all_count / total_all * total_sig, 1e-12),
      enrichment = count / expected,
      p_value    = phyper(count - 1, all_count, total_all - all_count, total_sig,
                          lower.tail = FALSE)
    ) %>%
    arrange(p_value)
}

	# Plot top enriched terms
plot_enrich <- function(df_annot, title) {
  if (nrow(df_annot) == 0) 
    return(ggplot() + theme_void() + ggtitle(paste(title, "(no terms)")))
  
  top <- df_annot %>%
    slice_min(p_value, n = min(10, nrow(df_annot))) %>%
    mutate(GO_Label = ifelse(is.na(GO_Name) | GO_Name == "", GOs, GO_Name),
           GO_Label = factor(GO_Label, levels = rev(GO_Label)))
  
  ggplot(top, aes(x = GO_Label, y = enrichment, fill = -log10(p_value))) +
    geom_col() + coord_flip() +
    scale_fill_gradient(low = "skyblue", high = "red") +
    labs(title = title, x = "GO term", y = "Enrichment", fill = "-log10(p)") +
    theme_minimal(base_size = 12)
}

# Output directories
dir_fig <- file.path(results_path, "figures_functional", "GO_enrichment_annotated")
dir_tab <- file.path(results_path, "tables_functional",  "GO_enrichment_annotated")
dir.create(dir_fig, TRUE, TRUE); dir.create(dir_tab, TRUE, TRUE)

# Main loop: contrast × direction × ontology
for (ct in contrastes) {
  tt <- dea_gene[[ct]]$full
  up_genes   <- rownames(tt[tt$FDR <= 0.05 & tt$logFC >=  1, , drop=FALSE])
  down_genes <- rownames(tt[tt$FDR <= 0.05 & tt$logFC <= -1, , drop=FALSE])
  
  for (dirn in c("UP","DOWN")) {
    gene_set <- if (dirn == "UP") up_genes else down_genes
    sig_gos  <- extract_gos(gene_set, ann) %>%
      left_join(go_dict, by = "GOs") %>%
      filter(!is.na(Ontology))
    
    # Recalculate enrichment per ontology category
    for (ont in c("BP","MF","CC")) {
      sig_ont <- sig_gos %>% dplyr::filter(Ontology == ont) %>% dplyr::select(GOs, count)
      uni_ont <- all_gos_onto %>% dplyr::filter(Ontology == ont)
      
      enrich <- enrich_calc(sig_ont, uni_ont) %>%
        left_join(go_dict, by = "GOs") %>%
        relocate(GOs, GO_Name, Ontology)
      
      # Save tables
      out_tsv <- file.path(dir_tab, paste0("GO_enrichment_", dirn, "_", ct, "_", ont, ".tsv"))
      write.table(enrich, out_tsv, sep="\t", quote=FALSE, row.names=FALSE)
      
      # Save plots
      p <- plot_enrich(enrich, paste0("Top GO (", ont, ") — ", dirn, " (", ct, ")"))
      ggsave(file.path(dir_fig, paste0("GO_enrich_", dirn, "_", ct, "_", ont, ".png")),
             p, width = 9, height = 6, dpi = 300)
    }
  }
}

message("Annotated GO enrichment (BP/MF/CC) completed for all contrasts.")
###############################################################################
