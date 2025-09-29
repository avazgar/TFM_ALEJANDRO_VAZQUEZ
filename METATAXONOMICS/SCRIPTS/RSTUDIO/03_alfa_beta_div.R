###########################################################################################################################
# Script: 02_data_preprocessing.R
# Description:
#       - Import tables and metadata from QIIME2 preprocessing
#       - Data formating
#       - Merge data 
#       - Create phyloseq object for further analysis
# Date: 23/09/2025
# Owner: Alejandro Vázquez García
# Context: Metataxonomic processing - ANALYSIS AND INTERPRETATION OF RESULTS - RSTUDIO
#########################################################################################################################

############################################ DIVERSITY ANALYSIS PREPROCESSING ###########################################

# Load phyloseq object (ps):
ps <- readRDS(file.path(processed_path, "boreal_soils.rds"))

# Transform to relative abundances for further analysis:
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
# Save phyloseq object of relative abundances:
saveRDS(ps, file = file.path(processed_path, "ps_rel.rds"))


# STUDY DATA DISTRIBUTION 
# Determine sequencing depth distribution
sample_sums(ps) %>% summary()
	#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	#30300   55972  106619   95928  121608  151048		# Specific to each sequencing 

# Determine depth per sample
sample_depths <- sample_sums(ps)

# Samples retained by sequencing depth value:

table(sample_depths >= 30000) #9
table(sample_depths >= 50000) #FALSE 1 TRUE 8
table(sample_depths >= 80000) #FALSE 3 TRUE 6		# Specific to each sequencing

# Rarefaction by minimum reads: Create rarefied phyloseq object

ps_rarefied <- rarefy_even_depth(
  ps,
  sample.size = 30000,
  rngseed = 42,		# Consider a concrete value as seed
  replace = FALSE,
  verbose = TRUE
)

	# `set.seed(42)` was used to initialize repeatable random subsampling.
	# Please record this for your records so others can reproduce.
	# Try `set.seed(42); .Random.seed` for the full vector
	# ...
	# 144OTUs were removed because they are no longer 
	# present in any sample after random subsampling

# DATASET FORMATTING
# Reorder the levels depending on the "Trt" variable <- SampleZone variable regarded as treatment
trt_ordered <- c("Oxidized","Transition","Unoxidized")

# Apply desired order to the factor:
sample_data(ps_rarefied)$SampleZone <- factor(sample_data(ps_rarefied)$SampleZone,
                                       levels = trt_ordered)

# Check if the factor is in the desired order:
levels(sample_data(ps_rarefied)$SampleZone)
#[1] "Oxidized"   "Transition" "Unoxidized"

######################################## TAXONOMY - MICROECO ##############################################

# Create MICROECO object
dataset <- microtable$new(otu_table)
dataset <- microtable$new(otu_table, sample_table = metadata)
dataset <- microtable$new(sample_table = metadata, otu_table, tax_table)
dataset$otu_table
dataset$sample_table
dataset$tax_table
dataset

# Rarefaction by minimum reads 
dataset$sample_sums() %>% range  # 30300 151048
dataset$sample_sums()
	#OX-1   OX-2   OX-3   TR-1   TR-2   TR-3   UN-1   UN-2   UN-3 
	#151048  30300 108124 121608  92786  50854  55972 106619 146037 		# Specific to each sequencing

dataset$rarefy_samples(sample.size=30000)
	#203 features are removed because they are no longer present in any sample after random subsampling ...
	#203 taxa with 0 abundance are removed from the otu_table ...

dataset$sample_sums() %>% range # Check the sequence numbers --> 30000 30000

# Apply tidy_dataset function to trim the dataset -> make the OTU and sample information consistent across all files in the dataset object
dataset$tidy_dataset
dataset

	# > dataset
	# microtable-class object:
	# sample_table have 9 rows and 5 columns
	# otu_table have 1959 rows and 9 columns
	# tax_table have 1959 rows and 8 columns

# Order the levels
dataset$sample_table$SampleZone %<>% factor(., levels = trt_ordered)

# Save filtered otu table, taxonomy and metadata and create the structure of taxonomy directories if they don't exist
output_dir <- file.path(processed_path, "taxonomy", "microeco","basic_files")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
dir.create(file.path(processed_path, "taxonomy", "Relabunds"))

# Save as each basic file as .xlsx
write.xlsx(dataset$otu_table, file.path(output_dir, "otu_table.xlsx"))
write.xlsx(dataset$sample_table, file.path(output_dir, "sample_table.xlsx"))
write.xlsx(dataset$tax_table, file.path(output_dir, "tax_table.xlsx"))

# RELATIVE ABUNDANCES

# Calculate and save the relative abundance at each taxonomic level for each ASV
	# Create cal_abund object
dataset$cal_abund()
	# Save as xlsx
write.xlsx(dataset$taxa_abund$Phylum, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_PHYLUM.xlsx"), rowNames = TRUE)
write.xlsx(dataset$taxa_abund$Class, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_CLASS.xlsx"), rowNames = TRUE)
write.xlsx(dataset$taxa_abund$Order, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_ORDER.xlsx"), rowNames = TRUE)
write.xlsx(dataset$taxa_abund$Family, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_FAMILY.xlsx"), rowNames = TRUE)
write.xlsx(dataset$taxa_abund$Genus, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_GENUS.xlsx"), rowNames = TRUE)


# Transformation and visualization of relative abundances for each taxonomic level grouped
	# Create trans_abund object
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 30) # Top 30 phylums
t1
t1$data_abund
t2 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 30)		# Top 30 classes
t3 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 30)		# Top 30 orders
t4 <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 30)	# Top 30 Families
t5 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 60)		# Top 30 Genus
	# Save as xlsx
write.xlsx(t1$data_abund, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_PHYLUM_group.xlsx"), rowNames = TRUE)
write.xlsx(t2$data_abund, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_CLASS_group.xlsx"), rowNames = TRUE)
write.xlsx(t3$data_abund, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_ORDER_group.xlsx"), rowNames = TRUE)
write.xlsx(t4$data_abund, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_FAMILY_group.xlsx"), rowNames = TRUE)
write.xlsx(t5$data_abund, file.path(processed_path, "taxonomy", "Relabunds", "RELABUND_GENUS_group.xlsx"), rowNames = TRUE)

dir.create(file.path(figures_path, "taxonomy", "microeco"), recursive = TRUE)

# Create the plots of relative abundances at phylum and class levels <- select only top 17 phylums so that the barplot does not appear saturated

	# Barplot & Boxplot Phylum level TOP 17
	# Barplot top 17 phylums:
t1_top <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 17)

g1 <- t1_top$plot_bar(others_color = "grey70", facet = "SampleZone", legend_text_italic = FALSE) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 15),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    strip.text   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 18),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 17),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Phylum", ncol = 1, byrow = TRUE)) +
  ggtitle("16S - Phylum - Relative Abundance (%) by oxidation gradient (Top 17 + Others)")

ggsave(file.path(figures_path, "taxonomy", "microeco","PHYLUM_REL_ABUND_BARPLOT_top17.png"),
       plot = g1, width = 20, height = 8, units = "in", dpi = 300)

	# Boxplot top 17 phylums:

g1.2 <- t1_top$plot_box(group = "SampleZone", xtext_angle = 45) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 15),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    strip.text   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 18),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 17),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Phylum", ncol = 1, byrow = TRUE)) +
  ggtitle("16S - Phylum - Relative Abundance (%) by oxidation gradient (Top 17)")
ggsave(file.path(figures_path, "taxonomy", "microeco","PHYLUM_REL_ABUND_BOXPLOT_top17.jpg"),
       plot = g1.2, width = 20, height = 8, units = "in", dpi = 300)
	# Barplot & Boxplot Class level top 17
	# Barplot top 17 classes:
t2_top <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 17)

g2 <- t2_top$plot_bar(others_color = "grey70", facet = "SampleZone", legend_text_italic = FALSE) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 15),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    strip.text   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 18),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 17),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Class", ncol = 1, byrow = TRUE)) +
  ggtitle("16S rRNA - Class - Relative Abundance (%) by oxidation gradient (Top 17 + Others)")
ggsave(file.path(figures_path, "taxonomy", "microeco","CLASS_REL_ABUND_BARPLOT_top17.png"),
       plot = g2, width = 20, height = 8, units = "in", dpi = 300)

	# Boxplot Top 17 classes:
g2.2 <- t2_top$plot_box(group = "SampleZone", xtext_angle = 45) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 18),
    strip.text   = element_text(size = 15, face = "bold"),
    plot.title   = element_text(size = 30),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Class", ncol = 1, byrow = TRUE)) +
  ggtitle("16S - Class - Relative Abundance (%) por by oxidation gradient (Top 17)")
ggsave(file.path(figures_path, "taxonomy", "microeco","CLASS_REL_ABUND_BOXPLOT_top17.jpg"),
       plot = g2.2, width = 20, height = 8, units = "in", dpi = 300)


	# Barplot & Boxplot Order level top 17
	# Barplot top 17 orders:
t3_top <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 17)

g3 <- t3_top$plot_bar(others_color = "grey70", facet = "SampleZone", legend_text_italic = FALSE) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    strip.text   = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Order", ncol = 2, byrow = TRUE)) +
  ggtitle("16S rRNA - Order - Relative Abundance (%) by oxidation gradient (Top 17 + Others)")

ggsave(file.path(figures_path, "taxonomy", "microeco","ORDER_REL_ABUND_BARPLOT_top17.png"),
       plot = g3, width = 20, height = 8, units = "in", dpi = 300)

	# Boxplot Top 17 orders:
g3.2 <- t3_top$plot_box(group = "SampleZone", xtext_angle = 45) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    strip.text   = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Order", ncol = 2, byrow = TRUE)) +
  ggtitle("16S rRNA - Order - Relative Abundance (%) by oxidation gradient (Top 17)")
ggsave(file.path(figures_path, "taxonomy", "microeco","ORDER_REL_ABUND_BOXPLOT_top17.jpg"),
       plot = g3.2, width = 20, height = 8, units = "in", dpi = 300)


	# Barplot & Boxplot Family level top 17
	# Barplot top 17 families:
t4_top <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 17)

g4 <- t4_top$plot_bar(others_color = "grey70", facet = "SampleZone", legend_text_italic = FALSE) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    strip.text   = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Order", ncol = 2, byrow = TRUE)) +
  ggtitle("16S - Famlily - Relative Abundance (%) por Gradiente de Oxidación (Top 17 + Others)")

ggsave(file.path(figures_path, "taxonomy", "microeco","FAMILY_REL_ABUND_BARPLOT_top17.png"),
       plot = g4, width = 20, height = 8, units = "in", dpi = 300)

	# Boxplot Top 17 families:
g4.2 <- t4_top$plot_box(group = "SampleZone", xtext_angle = 45) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    strip.text   = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Family", ncol = 2, byrow = TRUE)) +
  ggtitle("16S rRNA - Family - Relative Abundance (%) by oxidation gradient (Top 17)")
ggsave(file.path(figures_path, "taxonomy", "microeco","FAMILY_REL_ABUND_BOXPLOT_top17.jpg"),
       plot = g4.2, width = 20, height = 8, units = "in", dpi = 300)


	# Barplot & Boxplot Genus level top 17
	# Barplot top 17 families
t5_top <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 17)

g5 <- t5_top$plot_bar(others_color = "grey70", facet = "SampleZone", legend_text_italic = FALSE) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    strip.text   = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Genus", ncol = 2, byrow = TRUE)) +
  ggtitle("16S rRNA - Genus - Relative Abundance (%) by oxidation gradient (Top 17 + Others)")

ggsave(file.path(figures_path, "taxonomy", "microeco","GENUS_REL_ABUND_BARPLOT_top17.png"),
       plot = g5, width = 20, height = 8, units = "in", dpi = 300)

	# Boxplot Top 17 families:
g5.2 <- t5_top$plot_box(group = "SampleZone", xtext_angle = 45) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    strip.text   = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.size = unit(1.1, "lines")
  ) +
  guides(fill = guide_legend(title = "Genus", ncol = 2, byrow = TRUE)) +
  ggtitle("16S rRNA - Genus - Relative Abundance (%) by oxidation gradient (Top 17)")
ggsave(file.path(figures_path, "taxonomy", "microeco","GENUS_REL_ABUND_BOXPLOT_top17.jpg"),
       plot = g5.2, width = 20, height = 8, units = "in", dpi = 300)

	# Barplot & Boxplot genus level top 17
	# Barplot top 17 genus:
g5 <- t5$plot_bar(others_color = "grey70", facet = "SampleZone", legend_text_italic = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 10, face = "bold"), 
        plot.title = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 8)) +
  ggtitle("16S rRNA - Genus - Relative Abundance (%) by oxidation gradient")

ggsave(file.path(figures_path, "taxonomy", "microeco","GENUS_REL_ABUND_BARPLOT.jpg"),
       plot = g5, width = 20, height = 8, units = "in", dpi = 300)

	# Boxplot top 17 genus:
g5.2 <- t5$plot_box(group = "SampleZone", xtext_angle = 45) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 8)) +
  ggtitle("16S rRNA - Genus - Relative Abundance (%) by oxidation gradient")

ggsave(file.path(figures_path, "taxonomy", "microeco","GENUS_REL_ABUND_BOXPLOT.jpg"),
       plot = g5.2, width = 20, height = 8, units = "in", dpi = 300)


########################################### ALFA-DIVERSITY ############################################

dir.create(file.path(processed_path,"diversity","alfa-diversity"), recursive = TRUE)

# The rarefied phyloseq object (ps_rarefied) is considered to calculate alfa diversity indexes:
	# Calculate alfa diversity indexes
alpha_div <- estimate_richness(ps_rarefied, measures = c("Shannon", "Observed", "Chao1"))

# Data formatting
	# Add metadata
metadata <- as.data.frame(sample_data(ps_rarefied))
alpha_div$SampleZone <- metadata$SampleZone  # "Grouping" variable is applied here

	# Convert to long format. Adequate to manipulate data with ggplot2
alpha_long <- pivot_longer(alpha_div,
                           cols = c("Shannon", "Observed", "Chao1"),
                           names_to = "Index", values_to = "Value")

# Alfa diversity indexes representation
	# Boxplot by index and treatment
library(ggplot2)
p <- ggplot(alpha_long, aes(x = SampleZone, y = Value, fill = SampleZone)) +
  geom_boxplot() +
  facet_wrap(~ Index, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none"
  ) +
  labs(
    title = "Alpha Diversity Indexes by sample zone",
    x = "Sample Zone",
    y = "Diversity Index"
  )

	# Visualize and save the alfa-diversity plot

print(p)

ggsave(file.path(figures_path, "diversity","alfa-diversity","alpha_diversity_indexes_boxplot.png"),
       plot = p, width = 14, height = 8, dpi = 300)

write_xlsx(alpha_long, file.path(processed_path,"diversity", "alfa-diversity", "alpha_diversity_indices_data.xlsx"))


#####################################BETA-DIVERSITY################################################

# Beta-diversity assesment with no rarefaction:
	# References 
	# - McMurdie & Holmes 2014 
	# - Gloor et al. 2017

# Bray-Curtis on relative abundances

dir.create(file.path(processed_path,"diversity","beta-diversity"), recursive = TRUE)
	# Work with phyloseq object of relative abundances initially created (ps_rel)
bray_dist <- phyloseq::distance(ps_rel, method = "bray")

ordination_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)

p_bray <- plot_ordination(ps_rel, ordination_bray, color = "SampleZone") +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(color = SampleZone), type = "t", level = 0.95, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "PCoA - Bray-Curtis (Relative abundances)")
ggsave(file.path(figures_path, "diversity", "beta-diversity", "PCoAs", "PCoA_BrayCurtis.png"),
       plot = p_bray, width = 10, height = 6, dpi = 300)

	# PERMANOVA
meta <- as(sample_data(ps_rel), "data.frame")
permanova_bray <- adonis2(bray_dist ~ SampleZone, data = meta, permutations = 999)
print(permanova_bray)

# CLR + Aitchison (Euclidean)
ps_clr <- microbiome::transform(ps, "clr", zero.replace = 0.5)
clr_dist <- phyloseq::distance(ps_clr, method = "euclidean")

ordination_clr <- ordinate(ps_clr, method = "PCoA", distance = clr_dist)

p_clr <- plot_ordination(ps_clr, ordination_clr, color = "SampleZone") +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "PCoA - Aitchison (CLR + Euclidean)")
ggsave(file.path(figures_path, "diversity","beta-diversity","PCoAs","PCoA_CLR_Aitchison.png"),
       plot = p_clr, width = 10, height = 6, dpi = 300)

	# PERMANOVA
meta_clr <- as(sample_data(ps_clr), "data.frame")
permanova_clr <- adonis2(as.matrix(clr_dist) ~ SampleZone, data = meta_clr, permutations = 999)
print(permanova_clr)

# Hellinger transformation
ps_hell <- microbiome::transform(ps, "hellinger")
hell_dist <- phyloseq::distance(ps_hell, method = "bray")

ordination_hell <- ordinate(ps_hell, method = "PCoA", distance = hell_dist)

p_hell <- plot_ordination(ps_hell, ordination_hell, color = "SampleZone") +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(color = SampleZone), type = "t", level = 0.95, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "PCoA - Bray-Curtis (Hellinger-transformed)")
ggsave(file.path(figures_path,"diversity","beta-diversity","PCoAs","PCoA_Hellinger_BrayCurtis.png"),
       plot = p_hell, width = 10, height = 6, dpi = 300)

	# PERMANOVA
meta_hell <- as(sample_data(ps_hell), "data.frame")
permanova_hell <- adonis2(hell_dist ~ SampleZone, data = meta_hell, permutations = 999)
print(permanova_hell)

# Betadispersion & combined beta-diversity visualization
dispersion <- betadisper(bray_dist, group = meta$SampleZone)
anova_disp <- anova(dispersion)
write.csv(as.data.frame(anova_disp),
          file.path(processed_path,"diversity","beta-diversity","betadispersion_anova.csv"))
	
	# PCoA with ellipses
p1 <- plot_ordination(ps_rel, ordination_bray, color = "SampleZone") +
  stat_ellipse(type = "t", level = 0.95, linetype = 2, linewidth = 1.1) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCoA - Bray-Curtis (95% CI)")
	
	# Dispersion Boxplot
disp_df <- data.frame(DistanceToCentroid = dispersion$distances, Group = dispersion$group)
p2 <- ggplot(disp_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "Dispersion within Groups", x = "Group", y = "Distance to Centroid") +
  theme_minimal() + theme(legend.position = "none")

	# Combined visualization
combined_plot <- plot_grid(p1, p2, labels = c("A", "B"), rel_widths = c(1.4, 1), ncol = 2)
ggsave(file.path(figures_path,"diversity","beta-diversity","beta_diversity_combined.png"),
       plot = combined_plot, width = 12, height = 6, dpi = 300)
