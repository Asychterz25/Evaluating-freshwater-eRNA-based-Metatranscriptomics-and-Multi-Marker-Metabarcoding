#Fluoroprobe Beta-Diversity Ordinations

#Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)
library(cowplot)
library(vegan)

############################

#1)     Load Files

############################

#Load metadata
metadata <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata.tsv", sep="\t", header=TRUE, comment.char="")

#Load Fluoroprobe data
fluoro <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/F_August_10_LEAP2021.tsv", sep="\t", header=TRUE, comment.char="")

######## MT
#Load MT Algal identified through Class
filepath_class_mt <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/MT_class_algae_clean.tsv"
df_mt <- read_tsv(filepath_class_mt)

######## 18S
#Load MT Algal identified through Class
filepath_class_18S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/18S_PR2_class_algal_clean.tsv"
df_18S <- read_tsv(filepath_class_18S)

######## 23S
#Load MT Algal identified through Class
filepath_class_23S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/23S_Phytool_class_algal_clean.tsv"
df_23S <- read_tsv(filepath_class_23S)

############################

#2)     Proccess Files

############################

process_df <- function(df_taxon) {
  
  # 1) Remove the original Taxon column
  df_taxon <- df_taxon %>% select(-Taxon)
  
  # 2) Collapse rows by 'Algae' type (sum numeric columns for each type)
  df_collapsed <- df_taxon %>%
    group_by(Algae) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # 3) Compute relative abundance for all numeric columns
  numeric_cols <- df_collapsed %>% select(where(is.numeric))
  col_sums <- colSums(numeric_cols)
  col_sums[col_sums == 0] <- NA  # avoid division by zero
  
  df_rel <- df_collapsed
  df_rel[ , names(numeric_cols)] <- sweep(numeric_cols, 2, col_sums, FUN = "/") * 100
  df_rel[is.na(df_rel)] <- 0
  
  # 4) Set row names to Algae and remove the column
  df_rel <- as.data.frame(df_rel) 
  rownames(df_rel) <- df_rel$Algae
  df_rel$Algae <- NULL
  
  # 5) Transpose so samples are rows and taxa are columns
  df_transposed <- t(df_rel)
  
  return(df_transposed)
}

#Run function to process files
df_result_MT <- process_df(df_mt)
df_result_18S <- process_df(df_18S)
df_result_23S <- process_df(df_23S)


####Process Fluoroprobe data
fluoro_avg <- fluoro %>%
  group_by(SampleID) %>%                   
  summarise(across(where(is.numeric), 
                   mean, na.rm = TRUE))

#Get fluoroprobe averages
fluoro_avg <- fluoro_avg %>%
  select(SampleID, Green_Algae, Brown_Pigmented_Algae) %>%
  mutate(
    total = Green_Algae + Brown_Pigmented_Algae,
    total = replace(total, total == 0, NA),
    Green_Algae = (Green_Algae / total) * 100,
    Brown_Pigmented_Algae = (Brown_Pigmented_Algae / total) * 100
  ) %>%
  select(SampleID, Green_Algae, Brown_Pigmented_Algae)

# make SampleID the rownames
fluoro_avg <- as.data.frame(fluoro_avg)
rownames(fluoro_avg) <- fluoro_avg$SampleID
fluoro_avg$SampleID <- NULL

##############################

# 3) Bray-Curtis PCoA function

##############################

bray_curtis <- function(df_results, taxonomy) {
  
  # Bray-Curtis distance
  bray_dist <- vegdist(df_results, method = "bray")
  
  # PCoA
  bray_pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
  
  # Variance explained
  bray_var_exp <- round(bray_pcoa$eig / sum(bray_pcoa$eig) * 100, 1)
  
  # PCoA coordinates
  bray_pcoa_df <- as.data.frame(bray_pcoa$points) %>%
    rename(Axis.1 = V1, Axis.2 = V2) %>%
    mutate(
      Method = taxonomy,
      Axis1_Label = paste0("Axis 1 (", bray_var_exp[1], "%)"),
      Axis2_Label = paste0("Axis 2 (", bray_var_exp[2], "%)"),
      SampleID = rownames(.)
    )
  
  return(bray_pcoa_df)
}

# Run Bray-Curtis PCoA
BC_class_MT <- bray_curtis(df_result_MT, "Metatranscriptomics")
BC_class_18S <- bray_curtis(df_result_18S, "18S")
BC_class_23S <- bray_curtis(df_result_23S, "23S")
BC_fluoro <- bray_curtis(fluoro_avg, "Fluoroprobe")

# Combine into one ord_df
ord_df <- bind_rows(BC_class_MT, BC_class_18S, BC_class_23S, BC_fluoro)

##############################

# 4) Prepare metadata

##############################

metadata <- metadata %>%
  mutate(
    nutrient_level = ifelse(grepl("N", treatment), "High", "Moderate"),
    nutrient_level = factor(
      nutrient_level,
      levels = c("Moderate", "High")
    )
  )

# Merge metadata into ord_df
ord_df <- ord_df %>%
  left_join(metadata, by = "SampleID")

# Colors and shapes
nutrient_colors <- c("Moderate" = "#DDCC77", "High" = "#117733")

# Base theme
base_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  )

#Determine global axis limits
x_limits <- range(ord_df$Axis.1, na.rm = TRUE)
y_limits <- range(ord_df$Axis.2, na.rm = TRUE)

#Plot PCoA with global axis limits
bray_curtis_plot <- function(taxonomy) {
  
  df_filtered <- ord_df %>% filter(Method == taxonomy)
  
  ggplot(df_filtered, aes(x = Axis.1, y = Axis.2,
                          color = nutrient_level)) +
    geom_point(size = 3) +
    scale_color_manual(values = nutrient_colors, na.translate = FALSE) +
    labs(
      x = df_filtered$Axis1_Label[1],
      y = df_filtered$Axis2_Label[1],
      color = "Nutrient Level",
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    base_theme
}

##############################

# 5) Create plots

##############################

p1 <- bray_curtis_plot("Metatranscriptomics") + ggtitle("MT")
p2 <- bray_curtis_plot("18S") + ggtitle("18S")
p3 <- bray_curtis_plot("23S") + ggtitle("23S")
p4 <- bray_curtis_plot("Fluoroprobe") + ggtitle("FO")

# Add consistent title styling to base_theme plots
title_theme <- theme(
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
)
p1 <- p1 + title_theme
p2 <- p2 + title_theme
p3 <- p3 + title_theme
p4 <- p4 + title_theme

# Extract legend from one plot, with larger text, positioned at bottom
legend <- get_legend(
  p4 + 
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

#Create row spacer
row_spacer <- ggplot() + theme_void()

#Plot in 2x2
plots_grid <- plot_grid(
  p1, p2,
  row_spacer, row_spacer,
  p3, p4,
  ncol = 2,
  rel_heights = c(1, 0.1, 1),  # middle value controls gap size
  align = "hv",
  axis = "tb"
)

# Combine grid with bottom legend
final_plot <- plot_grid(
  plots_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

# Display final plot
final_plot

#Save
ggsave("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Fluoroprove_BC.png", final_plot, width = 10, height = 9)

##############################

# 6) Mantel tests

#############################

# 1) Subset Fluoroprobe to samples present in each marker
fluoro_MT <- fluoro_avg[rownames(fluoro_avg) %in% rownames(df_result_MT), ]
fluoro_18S <- fluoro_avg[rownames(fluoro_avg) %in% rownames(df_result_18S), ]
fluoro_23S <- fluoro_avg[rownames(fluoro_avg) %in% rownames(df_result_23S), ]

# 2) Subset markers to match Fluoroprobe
df_MT_sub <- df_result_MT[rownames(df_result_MT) %in% rownames(fluoro_MT), ]
df_18S_sub <- df_result_18S[rownames(df_result_18S) %in% rownames(fluoro_18S), ]
df_23S_sub <- df_result_23S[rownames(df_result_23S) %in% rownames(fluoro_23S), ]

# 3) Make sure rows are in the same order
df_MT_sub <- df_MT_sub[order(rownames(df_MT_sub)), ]
fluoro_MT <- fluoro_MT[order(rownames(fluoro_MT)), ]

df_18S_sub <- df_18S_sub[order(rownames(df_18S_sub)), ]
fluoro_18S <- fluoro_18S[order(rownames(fluoro_18S)), ]

df_23S_sub <- df_23S_sub[order(rownames(df_23S_sub)), ]
fluoro_23S <- fluoro_23S[order(rownames(fluoro_23S)), ]

# 4) Compute Bray-Curtis distances
dist_MT <- vegdist(df_MT_sub, method = "bray")
dist_18S <- vegdist(df_18S_sub, method = "bray")
dist_23S <- vegdist(df_23S_sub, method = "bray")

dist_fluoro_MT <- vegdist(fluoro_MT, method = "bray")
dist_fluoro_18S <- vegdist(fluoro_18S, method = "bray")
dist_fluoro_23S <- vegdist(fluoro_23S, method = "bray")

# 5) Run Mantel tests
mantel_MT <- mantel(dist_MT, dist_fluoro_MT, method = "spearman", permutations = 999)
mantel_18S <- mantel(dist_18S, dist_fluoro_18S, method = "spearman", permutations = 999)
mantel_23S <- mantel(dist_23S, dist_fluoro_23S, method = "spearman", permutations = 999)

# 6) Optional: organize results
mantel_results <- data.frame(
  Comparison = c("MT vs Fluoro", "18S vs Fluoro", "23S vs Fluoro"),
  Mantel_r = c(mantel_MT$statistic, mantel_18S$statistic, mantel_23S$statistic),
  p_value = c(mantel_MT$signif, mantel_18S$signif, mantel_23S$signif)
)

mantel_results

##############################

# 7) Betadisper (MT, 18S, 23S, Fluoroprobe)

##############################

# Betadisper function
beta_dispersion <- function(df_results, metadata, treatment_col, marker) {
  
  # Ensure sample order matches
  meta_sub <- metadata %>%
    filter(SampleID %in% rownames(df_results)) %>%
    arrange(match(SampleID, rownames(df_results)))
  
  # Bray-Curtis distance
  bray_dist <- vegdist(df_results, method = "bray")
  
  # Betadisper
  bd <- betadisper(bray_dist, meta_sub[[treatment_col]])
  
  # Permutation test
  bd_perm <- permutest(bd, permutations = 9999)
  
  # Extract distances to centroid
  bd_df <- data.frame(
    SampleID = rownames(df_results),
    Distance_to_Centroid = bd$distances,
    Treatment = meta_sub[[treatment_col]],
    Method = marker
  )
  
  return(list(
    betadisper = bd,
    permutest = bd_perm,
    distances = bd_df
  ))
}

# Run Betadisper for each marker/method
BD_MT <- beta_dispersion(df_result_MT, metadata, treatment_col = "nutrient_level", marker = "Metatranscriptomics")
BD_18S <- beta_dispersion(df_result_18S, metadata, treatment_col = "nutrient_level", marker = "18S")
BD_23S <- beta_dispersion(df_result_23S, metadata, treatment_col = "nutrient_level", marker = "23S")
BD_fluoro <- beta_dispersion(fluoro_avg, metadata, treatment_col = "nutrient_level", marker = "Fluoroprobe")

# Combine and set factor levels/labels
BD_all_raw <- bind_rows(
  BD_MT$distances,
  BD_18S$distances,
  BD_23S$distances,
  BD_fluoro$distances
) %>%
  mutate(
    Method = factor(
      Method,
      levels = c("Metatranscriptomics", "18S", "23S", "Fluoroprobe"),
      labels = c("MT", "18S", "23S", "FM")
    ),
    Treatment = factor(Treatment, levels = c("Moderate", "High"))
  )

nutrient_colors <- c("Moderate" = "#DDCC77", "High" = "#117733")

# Plot beta dispersion as box plot
beta_disp_plot <- ggplot(
  BD_all_raw,
  aes(x = Method, y = Distance_to_Centroid, fill = Treatment)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.size = 0.8,
    alpha = 0.9
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.2,
    alpha = 0.6,
    color = "black",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = nutrient_colors) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "right") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y.left = element_blank(),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    
    axis.text.x.top = element_text(size = 13, face = "bold", margin = margin(b = 8)),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    
    axis.text.y.left = element_blank(),
    axis.text.y.right = element_text(size = 12),
    axis.ticks.y.left = element_blank(),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.size = unit(1.2, "cm")
  ) +
  labs(
    y = "Distance to Centroid (β-dispersion)",
    fill = "Nutrient Enrichment"
  )

beta_disp_plot

ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/BD_BetaDispersion_Fluoro_boxplot.png",
  plot = beta_disp_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
