# TFM_ALEJANDRO_VAZQUEZ
# TFM 2025 — Metataxonomía y Metatranscriptómica
This repository includes the code used in this project. Reproducibility and scalability were considered as priorities during the coding steps. This project is divided in two main complementary branches:
- **METATAXONOMICS_TFM**: metataxonomics pipeline (16S rRNA gene).
- **METATRANSCRIPTOMICS_TFM**: metatranscriptomics pipeline (Differential Expression Analysis, KEGG Orthology).

Each project includes the structure of directories followed as well as a modular pipeline of scripts. 
The code is organized in two main blocks
1) **Bash Preprocessing**: Download, quality control, filtering and trimming, outputs in adequate format.
2) **R analysis**: statistics, visualization and final figures/tables generation.

> Note: Raw and intermediate processing files are not included, only expected paths.

## General structure metataxonomics project:
METATAXONOMICS_TFM/
├── data/
│   ├── raw/                   # Raw data (FASTQ, metadatos) 
│   └── processed/
│       ├── diversity/
│       │   ├── alfa-diversity/   # Alpha diversity results
│       │   └── beta-diversity/   # Beta diversity results (incl. PERMANOVA)
│       └── taxonomy/
│           ├── microeco/         # Microeco files
│           │   └── basic_files/
│           └── Relabunds/        # Relative abundance 
├── scripts/
│   ├── bash/                  # Bash preprocessing
│   └── R/                     # Statistic analysis RStudio
├── results/                   # Final tables and figures
│   ├── figures/
│   └── tables/
└── README.md                  # Pipeline guideline


## General structure metatranscriptomics project:
METATRANSCRIPTOMICS_TFM/
├── data/
│   ├── raw/                     # Raw data
│   │   ├── metadata.tsv         # Metadata
│   │   ├── salmon_counts.tsv    # Salmon Quantification files (counts)
│   │   ├── salmon_TPM.tsv       # Salmon Quantification files (TPM)
│   │   └── eggnog_annotation.emapper.annotations  # EggNOG annotations
│   └── processed/
│       ├── tximport/            # Import salmon results
│       │   ├── DGEList_gene_level.rds
│       │   ├── metadata_aligned.rds
│       │   └── salmon_counts_filtered.rds
│       ├── dds/                 # DESeq2/edgeR objects
│       │   ├── DEA_gene_level_results.rds
│       │   ├── design_matrix.rds
│       │   └── glmQLFit_object.rds
│       └── matrices/            # Functional matrix
│           ├── KO_matrix.tsv
│           ├── COG_matrix.tsv
│           └── GO_matrix.tsv
├── results/
│   └── DEA_gene_level/          # Differential expression analysis (DEA) results
│       ├── DEA_filtered_OX_VS_UN.tsv
│       ├── DEA_filtered_TR_VS_OX.tsv
│       ├── DEA_filtered_TR_VS_UN.tsv
│       ├── DEA_full_OX_VS_UN.tsv
│       ├── DEA_full_TR_VS_OX.tsv
│       ├── DEA_full_TR_VS_UN.tsv
│       └── Venn_DOWN_regions_counts.csv
├── figures/
│   └── DEA_gene_level/          # DEA: main figures
│       ├── MA_OX_VS_UN.png
│       ├── MA_TR_VS_OX.png
│       ├── MA_TR_VS_UN.png
│       ├── volcano_OX_VS_UN_enhanced.png
│       ├── volcano_TR_VS_OX_enhanced.png
│       ├── volcano_TR_VS_UN_enhanced.png
│       └── Venn_UP_proportional.png
├── scripts/
│   ├── bash/                    # BASH Preprocessing
│   └── R/                       # R analysis (DEA, enrichment, visualization)
└── README.md                    # Pipeline guideline

