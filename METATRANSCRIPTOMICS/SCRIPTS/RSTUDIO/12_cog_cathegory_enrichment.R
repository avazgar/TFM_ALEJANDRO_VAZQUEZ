###############################################################################
# Script: 37_cog_category_enrichment.R
# Description:
#   - Functional enrichment of COG categories (letters A–Z) using DEGs (gene-level)
#   - Expands multi-letter assignments (e.g. "EG" → "E","G")
#   - Applies hypergeometric test vs. background universe of annotated genes
#   - Runs separately for UP/DOWN contrasts
#   - Saves enrichment tables (TSV) and barplots (top categories per contrast)
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

# Load annotations and expand COG letters
ann_path <- file.path(raw_path, "eggnog_annotation.emapper.annotations")
hdr      <- grep("^#query", readLines(ann_path), value = FALSE)
ann      <- read_tsv(ann_path, skip = hdr - 1, show_col_types = FALSE)
names(ann)[names(ann) == "#query"] <- "query"

cog_long <- ann %>%
  filter(!is.na(COG_category), COG_category != "-") %>%
  mutate(COG_category = gsub("\\s", "", COG_category)) %>%
  separate_rows(COG_category, sep = "") %>%
  filter(COG_category != "")

# COG dictionary (letter --> name) 
cog_dict <- tribble(
  ~COG_category, ~COG_name,
  "J","Translation, ribosomal structure and biogenesis",
  "A","RNA processing and modification",
  "K","Transcription",
  "L","Replication, recombination and repair",
  "B","Chromatin structure and dynamics",
  "D","Cell cycle control, cell division, chromosome partitioning",
  "Y","Nuclear structure",
  "V","Defense mechanisms",
  "T","Signal transduction mechanisms",
  "M","Cell wall/membrane/envelope biogenesis",
  "N","Cell motility",
  "Z","Cytoskeleton",
  "W","Extracellular structures",
  "U","Intracellular trafficking, secretion, and vesicular transport",
  "O","Posttranslational modification, protein turnover, chaperones",
  "C","Energy production and conversion",
  "G","Carbohydrate transport and metabolism",
  "E","Amino acid transport and metabolism",
  "F","Nucleotide transport and metabolism",
  "H","Coenzyme transport and metabolism",
  "I","Lipid transport and metabolism",
  "P","Inorganic ion transport and metabolism",
  "Q","Secondary metabolites biosynthesis, transport and catabolism",
  "R","General function prediction only",
  "S","Function unknown"
)

cog_long <- cog_long %>% left_join(cog_dict, by = "COG_category")

	# Background universe: number of genes per category
uni_cog <- cog_long %>% count(COG_category, name = "all_count")

# Load DEG results
dea_gene  <- readRDS(file.path(processed_path, "dds", "DEA_gene_level_results.rds"))
contrasts <- names(dea_gene)

# Enrichment function
enrich_cog <- function(genes_sig, cog_long, uni_cog) {
  sig <- cog_long %>% filter(query %in% genes_sig) %>% count(COG_category, name = "sig_count")
  total_sig <- sum(sig$sig_count, na.rm = TRUE)
  total_all <- sum(uni_cog$all_count, na.rm = TRUE)
  if (total_sig == 0 || total_all == 0) {
    return(tibble(COG_category=character(), sig_count=integer(),
                  all_count=integer(), expected=double(), enrichment=double(),
                  p_value=double()))
  }
  sig %>%
    right_join(uni_cog, by="COG_category") %>%
    mutate(sig_count = replace_na(sig_count, 0L)) %>%
    mutate(
      expected   = pmax(all_count / total_all * total_sig, 1e-12),
      enrichment = sig_count / expected,
      p_value    = phyper(sig_count - 1, all_count, total_all - all_count, total_sig, lower.tail = FALSE)
    ) %>%
    arrange(p_value)
}

# Output directories
dir_fig <- file.path(results_path, "figures_functional", "COG_enrichment")
dir_tab <- file.path(results_path, "tables_functional", "COG_enrichment")
dir.create(dir_fig, TRUE, TRUE); dir.create(dir_tab, TRUE, TRUE)

# Plotting function
plot_cog <- function(df, title, top_n = 15) {
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(title, "(no categories)")))
  }
  
  top <- df %>%
    filter(!is.na(p_value)) %>%
    left_join(cog_dict, by = "COG_category") %>%
    mutate(label = ifelse(is.na(COG_name),
                          COG_category,
                          paste0(COG_category, ": ", COG_name)))
  
  top <- top %>%
    slice_min(order_by = p_value, n = min(top_n, nrow(top)), with_ties = FALSE) %>%
    mutate(label = factor(label, levels = rev(label)))
  
  ggplot(top, aes(x = label, y = enrichment, fill = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "skyblue", high = "red") +
    labs(title = title, x = "COG category", y = "Enrichment", fill = "-log10(p)") +
    theme_minimal(base_size = 12)
}

# Loop over contrasts × directions
for (ct in contrasts) {
  tt <- dea_gene[[ct]]$full
  G_up   <- rownames(tt[tt$FDR <= 0.05 & tt$logFC >=  1, , drop=FALSE])
  G_down <- rownames(tt[tt$FDR <= 0.05 & tt$logFC <= -1, , drop=FALSE])
  
  for (dirn in c("UP","DOWN")) {
    Gs <- if (dirn == "UP") G_up else G_down
    enr <- enrich_cog(Gs, cog_long, uni_cog)
    
    # Save table
    write.table(enr, file.path(dir_tab, paste0("COG_", dirn, "_", ct, ".tsv")),
                sep="\t", quote=FALSE, row.names=FALSE)
    
    # Save plot
    p <- plot_cog(enr, paste0("COG enrichment — ", dirn, " (", ct, ")"))
    ggsave(file.path(dir_fig, paste0("COG_", dirn, "_", ct, ".png")),
           p, width = 10, height = 6, dpi = 300)
  }
}

message("COG category enrichment completed.")
###############################################################################
