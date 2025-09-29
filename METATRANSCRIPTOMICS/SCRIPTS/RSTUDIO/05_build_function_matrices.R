######################################################################################
# Script: 05_build_function_matrices.R
# Description:
# eggNOG-mapper is a functional annotion tool based on eggNOG database. It integrates orthology and annotations in different levels (KO, CO, GO, Pfam, etc.). This tool allows
# map genes or ORFs to known functions, facilitating the functional analysis
#   - Read eggNOG-mapper Leer anotaciones de eggNOG-mapper
#   - Map genes/ORFs -< KO/GO/COG
#   - Sum counts by function and build function matrix by sample
#   - Save matrix for further subsequent analysis (DEA, heatmaps)
# Date: 23/09/25
# Owner: Alejandro Vázquez García
# Context: Transcriptomic processing  - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
######################################################################################

# Inputs: Counts and eggNOG annotations

counts <- readRDS(file.path(processed_path, "tximport", "salmon_counts_filtered.rds"))
meta   <- readRDS(file.path(processed_path, "tximport", "metadata_aligned.rds"))
annotations_file_path <- file.path(raw_path, "eggnog_annotation.emapper.annotations")

# Find header: regex 
header_line <- grep("^#query", readLines(annotations_file_path), value = FALSE)
ann <- read_tsv(annotations_file_path, skip = header_line - 1)
names(ann) <- sub("^#\\s*", "", names(ann))

# Select columns
col_fn <- list(
  KO  = "KEGG_ko",
  COG = "COG_category",
  GO  = "GOs"
)

missing_cols <- setdiff(unlist(col_fn), names(ann))
if (length(missing_cols) > 0) {
  warning("Columns not found in annotations: ", paste(missing_cols, collapse = ", "))
}

# Align genes between counts and annotations

genes_common <- intersect(rownames(counts), ann$query)
counts <- counts[genes_common, , drop = FALSE]
ann    <- dplyr::filter(ann, query %in% genes_common)

# Build function: Matrix by sample

build_matrix <- function(counts_mat, ann_tbl, fun_col, fun_sep = ",") {
  if (!(fun_col %in% names(ann_tbl))) return(matrix(0, nrow = 0, ncol = ncol(counts_mat)))
  
  ann_long <- ann_tbl %>%
    dplyr::select(query, fun = all_of(fun_col)) %>%
    dplyr::mutate(fun = ifelse(is.na(fun), "-", fun)) %>%
    tidyr::separate_rows(fun, sep = fun_sep) %>%
    dplyr::mutate(fun = stringr::str_trim(fun)) %>%
    dplyr::filter(fun != "-")

  if (nrow(ann_long) == 0) return(matrix(0, nrow = 0, ncol = ncol(counts_mat)))

  counts_df <- tibble::rownames_to_column(as.data.frame(counts_mat), "query")

  mat_fun <- ann_long %>%
    dplyr::inner_join(counts_df, by = "query") %>%
    dplyr::select(-query) %>%
    dplyr::group_by(fun) %>%
    dplyr::summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames("fun") %>%
    as.matrix()

  # Order columns according to metadata object
  mat_fun <- mat_fun[, rownames(meta), drop = FALSE]
  return(mat_fun)
}

# Build KO/CG/GO matrix 

KO_counts  <- build_matrix(counts, ann, col_fn$KO,  fun_sep = ",")
COG_counts <- build_matrix(counts, ann, col_fn$COG, fun_sep = "")
GO_counts  <- build_matrix(counts, ann, col_fn$GO,  fun_sep = ",")

message("Resumen funciones: ",
        "KO=", nrow(KO_counts), ", ",
        "COG=", nrow(COG_counts), ", ",
        "GO=", nrow(GO_counts))

# Save matrix

dir.create(file.path(processed_path, "matrices"), showWarnings = FALSE, recursive = TRUE)

saveRDS(KO_counts,  file.path(processed_path, "matrices", "KO_counts.rds"))
saveRDS(COG_counts, file.path(processed_path, "matrices", "COG_counts.rds"))
saveRDS(GO_counts,  file.path(processed_path, "matrices", "GO_counts.rds"))

# Export TSVs
readr::write_tsv(tibble::rownames_to_column(as.data.frame(KO_counts), "Function"),
                 file.path(processed_path, "matrices", "KO_matrix.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(COG_counts), "Function"),
                 file.path(processed_path, "matrices", "COG_matrix.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(GO_counts), "Function"),
                 file.path(processed_path, "matrices", "GO_matrix.tsv"))

message("Functional matrix created and saved in the corresponding directories")
