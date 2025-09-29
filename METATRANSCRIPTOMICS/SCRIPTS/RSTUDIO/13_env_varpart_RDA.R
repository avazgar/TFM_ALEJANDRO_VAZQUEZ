###############################################################################
# Script: 13_env_PCA_varpart_RDA.R
# Description:
#   - Environmental PCA (standardized), loadings/scores export and scree plot
#   - Hellinger transform of KO/COG/GO (samples × features)
#   - Variation partitioning (vegan::varpart) using PC1 and PC2 as orthogonal drivers
#   - RDA using PCs (PC1+PC2) with envfit overlay of real variables
#   - RDA using original environmental variables (pH, Fe2+, TRS)
#   - Consistent export of tables and figures
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Transcriptomic functional profiling — Environmental drivers - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
###############################################################################


source("scripts/00_configuration.R")
source("scripts/01_packages.R")   # vegan, ggplot2, ggvegan, ggrepel, patchwork, scales, etc.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales)
  library(vegan); library(ggvegan)
})

# ENV data & PCA 

	# Rename if needed
if ("pH_mean" %in% names(meta) && !"pH" %in% names(meta)) {
  names(meta)[names(meta) == "pH_mean"] <- "pH"
}

	# Select variables
stopifnot(all(c("pH","Fe2+","TRS") %in% colnames(meta)))

ENV_raw <- meta[, c("pH","Fe2+","TRS")]
rownames(ENV_raw) <- rownames(meta)

	# Keep only complete cases
keep_env <- stats::complete.cases(ENV_raw)
ENV_raw  <- ENV_raw[keep_env, , drop = FALSE]
meta_env <- meta[keep_env, , drop = FALSE]

	# Standardization and PCA
ENV_std  <- scale(ENV_raw)                         
pca_env  <- prcomp(ENV_std, center = FALSE, scale. = FALSE)

	# Scores (samples) and loadings (variables)
PC_scores   <- as.data.frame(pca_env$x[, 1:2, drop = FALSE])
PC_loadings <- as.data.frame(pca_env$rotation[, 1:2, drop = FALSE])
PC_scores$SampleID <- rownames(PC_scores)

	# Explained variance
imp      <- summary(pca_env)$importance
var_prop <- imp["Proportion of Variance", 1:ncol(pca_env$rotation)]

# Input/Output structure
dir_tab <- file.path(results_path, "tables_env")
dir_fig <- file.path(results_path, "figures_env")
dir.create(dir_tab, TRUE, TRUE); dir.create(dir_fig, TRUE, TRUE)

	# Export correlations and PCA tables
write.table(round(cor(ENV_raw, method = "spearman"), 3),
            file.path(dir_tab, "ENV_spearman_cor.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(round(PC_loadings, 4),
            file.path(dir_tab, "ENV_PCA_loadings.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(round(PC_scores[, c("PC1","PC2")], 4),
            file.path(dir_tab, "ENV_PCA_scores.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

	# Scree plot
var_df <- data.frame(PC = paste0("PC", seq_along(var_prop)),
                     VarExp = as.numeric(var_prop))
p_scree <- ggplot(var_df, aes(x = PC, y = VarExp)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = scales::percent(VarExp, accuracy = 0.1)),
            vjust = -0.4, size = 3.6) +
  labs(title = "Explained variance by PCs", y = "Proportion", x = NULL) +
  theme_minimal(base_size = 13)
ggsave(file.path(dir_fig, "ENV_PCA_scree.png"), p_scree, width = 7.5, height = 4.5, dpi = 300)

	# PCA biplot (samples + variable vectors)
scores_df  <- PC_scores |>
  dplyr::left_join(tibble::rownames_to_column(meta_env, "SampleID"), by = "SampleID")
loads_df   <- tibble::rownames_to_column(PC_loadings, "Variable")
loads_plot <- loads_df; loads_plot[, c("PC1","PC2")] <- loads_plot[, c("PC1","PC2")] * 3

p_biplot <- ggplot(scores_df, aes(x = PC1, y = PC2, color = SampleZone)) +
  geom_point(size = 3.6, alpha = 0.9) +
  geom_text_repel(aes(label = SampleID), size = 3) +
  geom_segment(data = loads_plot,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  geom_text_repel(data = loads_plot,
                  aes(x = PC1, y = PC2, label = Variable),
                  size = 3.4, color = "black") +
  labs(
    title = "Environmental PCA biplot",
    x = sprintf("PC1 (%.1f%%)", 100 * var_prop[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * var_prop[2])
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(dir_fig, "ENV_PCA_biplot.png"), p_biplot, width = 8.5, height = 6.5, dpi = 300)

# Functional matrices (Hellinger)

KO   <- readRDS(file.path(processed_path, "matrices", "KO_counts.rds"))
COG  <- readRDS(file.path(processed_path, "matrices", "COG_counts.rds"))
GO   <- readRDS(file.path(processed_path, "matrices", "GO_counts.rds"))

	# Align with ENV samples
stopifnot(identical(colnames(KO), rownames(meta)))
KO   <- KO[, rownames(meta_env), drop = FALSE]
COG  <- COG[, rownames(meta_env), drop = FALSE]
GO   <- GO[, rownames(meta_env), drop = FALSE]

prep_hellinger <- function(M_counts) decostand(t(M_counts), method = "hellinger")
Y_KO  <- prep_hellinger(KO)
Y_COG <- prep_hellinger(COG)
Y_GO  <- prep_hellinger(GO)

PC_env <- PC_scores[rownames(Y_KO), c("PC1","PC2"), drop = FALSE]

# Variation partitioning (PC1 vs PC2)

run_varpart_pcs <- function(Y, label) {
  vp  <- varpart(Y, ~ PC1, ~ PC2, data = PC_env)
  sink(file.path(dir_tab, sprintf("varpart_PCs_%s.txt", label))); print(vp); sink()
  png(file.path(dir_fig, sprintf("varpart_PCs_%s.png", label)), 1200, 900, res = 140)
  plot(vp, bg = c("#1f77b4", "#ff7f0e"), Xnames = c("PC1","PC2"),
       main = sprintf("Variation partitioning (%s) — PC1 / PC2", label))
  dev.off()
  invisible(vp)
}

vp_KO  <- run_varpart_pcs(Y_KO,  "KO")
vp_COG <- run_varpart_pcs(Y_COG, "COG")
vp_GO  <- run_varpart_pcs(Y_GO,  "GO")

# RDA with PCs (PC1 + PC2)

run_rda_with_pcs <- function(Y, label) {
  mod <- rda(Y ~ PC1 + PC2, data = PC_env)

  aov_glob  <- anova(mod, permutations = 9999)
  aov_terms <- anova(mod, by = "term", permutations = 9999)
  capture.output(aov_glob,  file = file.path(dir_tab, sprintf("RDA_PCs_global_%s.txt", label)))
  capture.output(aov_terms, file = file.path(dir_tab, sprintf("RDA_PCs_terms_%s.txt",  label)))

  	# Base RDA biplot + envfit overlay
  png(file.path(dir_fig, sprintf("RDA_PCs_biplot_%s_base.png", label)), 1200, 900, res = 140)
  plot(mod, scaling = 2, main = sprintf("RDA (Hellinger) — %s ~ PC1 + PC2", label))
  ef <- envfit(mod, ENV_std, permutations = 0, scaling = 2)
  plot(ef, col = "red", add = TRUE)
  dev.off()

 	# ggplot version
  site_scores  <- as.data.frame(scores(mod, display = "sites", scaling = 2))
  site_scores$SampleID <- rownames(site_scores)
  site_scores  <- dplyr::left_join(site_scores,
                                   tibble::rownames_to_column(meta_env, "SampleID"),
                                   by = "SampleID")
  arrow_scores <- as.data.frame(scores(mod, display = "bp", scaling = 2))
  arrow_scores$Axis <- rownames(arrow_scores)

  p <- ggplot(site_scores, aes(x = RDA1, y = RDA2, color = SampleZone)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_segment(data = arrow_scores,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(0.20, "cm")), color = "black") +
    geom_text(data = arrow_scores,
              aes(x = RDA1 * 1.08, y = RDA2 * 1.08, label = Axis),
              size = 3.3, color = "black") +
    labs(title = sprintf("RDA biplot (%s ~ PC1 + PC2)", label)) +
    theme_minimal(base_size = 13)
  ggsave(file.path(dir_fig, sprintf("RDA_PCs_biplot_%s_ggplot.png", label)),
         p, width = 8.5, height = 6.5, dpi = 300)

  invisible(list(model = mod, aov_global = aov_glob, aov_terms = aov_terms))
}

rdaPC_KO  <- run_rda_with_pcs(Y_KO,  "KO")
rdaPC_COG <- run_rda_with_pcs(Y_COG, "COG")
rdaPC_GO  <- run_rda_with_pcs(Y_GO,  "GO")

# RDA with original ENV variables (pH, Fe2+, TRS)

ENV_model <- meta_env[, c("pH","Fe2+","TRS")]
stopifnot(identical(rownames(ENV_model), rownames(Y_KO)))

run_rda_with_env <- function(Y, label) {
  mod <- rda(Y ~ pH + `Fe2+` + TRS, data = ENV_model)

  aov_glob  <- anova(mod, permutations = 9999)
  aov_terms <- anova(mod, by = "term", permutations = 9999)
  capture.output(aov_glob,  file = file.path(dir_tab, sprintf("RDA_ENV_global_%s.txt", label)))
  capture.output(aov_terms, file = file.path(dir_tab, sprintf("RDA_ENV_terms_%s.txt",  label)))

 	# ggvegan autoplot
  p <- autoplot(mod, scaling = 2) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("RDA — %s vs ENV (pH, Fe2+, TRS)", label))
  ggsave(file.path(dir_fig, sprintf("RDA_ENV_%s_autoplot.png", label)),
         plot = p, width = 8, height = 6.5, dpi = 300)

  invisible(list(model = mod, aov_global = aov_glob, aov_terms = aov_terms))
}

rdaENV_KO  <- run_rda_with_env(Y_KO,  "KO")
rdaENV_COG <- run_rda_with_env(Y_COG, "COG")
rdaENV_GO  <- run_rda_with_env(Y_GO,  "GO")

message("ENV PCA, varpart (PC1/PC2), and RDA (PCs & ENV) completed. Outputs in:")
message("  - Tables:  ", dir_tab)
message("  - Figures: ", dir_fig)

###############################################################################
# OPTIONAL: Hierarchical Partitioning with rdacca.hp
# - Useful when many environmental variables are included (>3 predictors)
# - Provides averaged contributions (unique + shared) across all combinations
###############################################################################

suppressPackageStartupMessages({ library(rdacca.hp) })

# Prepare response matrices (Hellinger transformed)
KO_hel  <- decostand(t(KO),  method = "hellinger")
COG_hel <- decostand(t(COG), method = "hellinger")
GO_hel  <- decostand(t(GO),  method = "hellinger")

# Full set of environmental variables (example with pH, Fe2+, TRS)
env_vars <- meta[, c("pH", "Fe2+", "TRS")]

# Function to run rdacca.hp
run_hp <- function(Y, Y_name, env_df, outdir) {
  hp <- rdacca.hp(Y = Y, X = env_df,
                  method = "RDA", type = "adjR2",
                  trace = TRUE, nperm = 9999)
  
  # Export table
  tab_file <- file.path(outdir, paste0("rdacca_hp_", Y_name, ".tsv"))
  write.table(hp$Var.part, tab_file, sep = "\t", quote = FALSE)
  
  # Export figure (barplot of contributions)
  png(file.path(outdir, paste0("rdacca_hp_", Y_name, ".png")), 
      width = 900, height = 700, res = 130)
  barplot(hp$Var.part$Independent, las = 2, col = "skyblue",
          main = paste("Hierarchical partitioning —", Y_name),
          ylab = "Independent contribution (adj R²)")
  dev.off()
  
  invisible(hp)
}

# Example: run for KO, COG, GO
hp_dir <- file.path(results_path, "tables_functional", "ENV_hp")
dir.create(hp_dir, TRUE, TRUE)

hp_KO  <- run_hp(KO_hel,  "KO",  env_vars, hp_dir)
hp_COG <- run_hp(COG_hel, "COG", env_vars, hp_dir)
hp_GO  <- run_hp(GO_hel,  "GO",  env_vars, hp_dir)

message("rdacca.hp completed (optional hierarchical partitioning).")
