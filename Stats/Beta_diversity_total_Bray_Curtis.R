########Beta Diversity
# This code with create PcOA for Bray-Curtis based on MT, MM, 18S, COI and 23S at the Class, Family and Genus level
# MT = Metatranscriptomics
# MM = Multi-marker metabarcoding

#Load libraries
library(vegan)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)

############################

#1)     Load Files

############################

#Load metadata
metadata <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata.tsv", sep="\t", header=TRUE, comment.char="")

#Load MT metadata
metadata_mt <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata_mt.tsv", sep="\t", header=TRUE, comment.char="")

######## MT
#Load MT Class
filepath_class_mt <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MT_class_raw.tsv"
df_class_mt <- read_tsv(filepath_class_mt)

#Load MT Family
filepath_family_mt <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MT_family_raw.tsv"
df_family_mt <- read_tsv(filepath_family_mt)

#Load MT Genus
filepath_genus_mt <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MT_genus_raw.tsv"
df_genus_mt <- read_tsv(filepath_genus_mt)

######### MM
#Load MM Class
filepath_class_mm <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MM_class_p_raw.tsv"
df_class_mm <- read_tsv(filepath_class_mm)

#Load MM Family
filepath_family_mm <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MM_family_p_raw.tsv"
df_family_mm <- read_tsv(filepath_family_mm)

#Load MM Genus
filepath_genus_mm <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MM_genus_p_raw.tsv"
df_genus_mm <- read_tsv(filepath_genus_mm)

######### 18S
#Load 18S Class
filepath_class_18S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/18S_PR2_class_raw.tsv"
df_class_18S <- read_tsv(filepath_class_18S)

#Load 18S Family
filepath_family_18S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/18S_PR2_family_raw.tsv"
df_family_18S <- read_tsv(filepath_family_18S)

#Load 18S Genus
filepath_genus_18S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/18S_PR2_genus_raw.tsv"
df_genus_18S <- read_tsv(filepath_genus_18S)

######### COI
#Load COI Class
filepath_class_COI <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/COI_Porter_class_raw.tsv"
df_class_COI <- read_tsv(filepath_class_COI)

#Load COI Family
filepath_family_COI <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/COI_Porter_family_raw.tsv"
df_family_COI <- read_tsv(filepath_family_COI)

#Load COI Genus
filepath_genus_COI <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/COI_Porter_genus_raw.tsv"
df_genus_COI <- read_tsv(filepath_genus_COI)

######### 23S
#Load 23S Class
filepath_class_23S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/23S_Phytool_class_raw.tsv"
df_class_23S <- read_tsv(filepath_class_23S)

#Load 23S Family
filepath_family_23S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/23S_Phytool_family_raw.tsv"
df_family_23S <- read_tsv(filepath_family_23S)

#Load 23S Genus
filepath_genus_23S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/23S_Phytool_genus_raw.tsv"
df_genus_23S <- read_tsv(filepath_genus_23S)

############################

#2)     Process Files

############################

# Function to process files
process_df <- function(df_taxon) {
  # Columns to remove
  cols_to_remove <- c("False_Positives", "Non-Target_Taxa", "Reason", 
                      "Mean", "DB_R1", "DB_R2", "DB_R3", "EB_21", "Marker")
  
  # 1) Remove rows where "Reason" is not blank
  df_clean <- df_taxon[df_taxon$Reason == "" | is.na(df_taxon$Reason), ]
  
  # 2) Remove unwanted columns
  df_clean <- df_clean[, !(names(df_clean) %in% cols_to_remove)]
  
  # 3) Convert to relative abundance (% of column sum)
  taxa <- df_clean$Taxon
  df_clean <- df_clean[, !(names(df_clean) %in% c("Taxon"))]
  
  # Avoid division by zero or NA column sums
  col_sums <- colSums(df_clean, na.rm = TRUE)
  col_sums[col_sums == 0] <- NA  # prevent division by zero
  
  df_rel <- sweep(df_clean, 2, col_sums, FUN = "/") * 100
  df_rel[is.na(df_rel)] <- 0     # replace NaN results with 0
  
  # 4) Set row names to Taxon
  rownames(df_rel) <- taxa
  
  # 5) Transpose (samples as rows, taxa as columns)
  df_transposed <- t(df_rel)
  
  return(df_transposed)
}

#Run function to clean files
df_result_class_MT <- process_df(df_class_mt)
df_result_family_MT <- process_df(df_family_mt)
df_result_genus_MT <- process_df(df_genus_mt)
df_result_class_MM <- process_df(df_class_mm)
df_result_family_MM <- process_df(df_family_mm)
df_result_genus_MM <- process_df(df_genus_mm)
df_result_class_18S <- process_df(df_class_18S)
df_result_family_18S <- process_df(df_family_18S)
df_result_genus_18S <- process_df(df_genus_18S)
df_result_class_COI <- process_df(df_class_COI)
df_result_family_COI <- process_df(df_family_COI)
df_result_genus_COI <- process_df(df_genus_COI)
df_result_class_23S <- process_df(df_class_23S)
df_result_family_23S <- process_df(df_family_23S)
df_result_genus_23S <- process_df(df_genus_23S)

##############################

#3) Calculating Beta Diversity

##############################

# Function to calculate Bray-Curtis diversity and PCoA
bray_curtis <- function(df_results,taxonomy,marker) {
  # 1) Bray-Curtis Distance calculation
  bray_dist <- vegdist(df_results, method = "bray")
  
  # 2) PCoA calculations
  bray_pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
  
  # 3) Calculate percentage variance explained
  bray_var_exp <- round(bray_pcoa$eig / sum(bray_pcoa$eig) * 100, 1)
  
  # 4) Extract PCoA coordinates
  bray_pcoa_df <- as.data.frame(bray_pcoa$points) %>%
    rename(Axis.1 = V1, Axis.2 = V2) %>%
    mutate(
      Method = taxonomy,
      Metric = marker,
      Axis1_Label = paste0("Axis 1 (", bray_var_exp[1], "%)"),
      Axis2_Label = paste0("Axis 2 (", bray_var_exp[2], "%)"),
      SampleID = rownames(.)
    )
  
  # 5) Return the result
  return(bray_pcoa_df)
}

#Run function to calculate Bray Curtis
BC_class_MT <- bray_curtis(df_result_class_MT, "Class", "Metatranscriptomics")
BC_family_MT <- bray_curtis(df_result_family_MT, "Family", "Metatranscriptomics")
BC_genus_MT <- bray_curtis(df_result_genus_MT, "Genus", "Metatranscriptomics")
BC_class_MM <- bray_curtis(df_result_class_MM, "Class", "Multi-Marker")
BC_family_MM <- bray_curtis(df_result_family_MM, "Family", "Multi-Marker")
BC_genus_MM <- bray_curtis(df_result_genus_MM, "Genus", "Multi-Marker")
BC_class_18S <- bray_curtis(df_result_class_18S, "Class", "18S")
BC_family_18S <- bray_curtis(df_result_family_18S, "Family", "18S")
BC_genus_18S <- bray_curtis(df_result_genus_18S, "Genus", "18S")
BC_class_COI <- bray_curtis(df_result_class_COI, "Class", "COI")
BC_family_COI <- bray_curtis(df_result_family_COI, "Family", "COI")
BC_genus_COI <- bray_curtis(df_result_genus_COI, "Genus", "COI")
BC_class_23S <- bray_curtis(df_result_class_23S, "Class", "23S")
BC_family_23S <- bray_curtis(df_result_family_23S, "Family", "23S")
BC_genus_23S <- bray_curtis(df_result_genus_23S, "Genus", "23S")

##############################

#4).     Prepare Plot

##############################

# Combine all ordination results
ord_df <- bind_rows(
  BC_class_MT, BC_family_MT, BC_genus_MT,
  BC_class_MM, BC_family_MM, BC_genus_MM,
  BC_class_18S, BC_family_18S, BC_genus_18S,
  BC_class_COI, BC_family_COI, BC_genus_COI,
  BC_class_23S, BC_family_23S, BC_genus_23S
)

# Merge with metadata so only existing SampleIDs are kept
ord_df <- left_join(ord_df, metadata, by = "SampleID")

#Add temp and nutrient columns to metadata_test
metadata <- metadata %>%
  mutate(
    temp_level = case_when(
      grepl("T", treatment) ~ "Heated",
      TRUE ~ "Ambient"
    ),
    nutrient_level = factor(
      case_when(
        grepl("N", treatment) ~ "High",
        TRUE ~ "Moderate"
      ),
      levels = c("Moderate", "High")
    )
  )

#Merge with ord_df again if metadata_test was changed
ord_df <- left_join(ord_df, metadata, by = c("SampleID" = "SampleID"))

#Define two colors and two shapes
nutrient_colors <- c("Moderate" = "#DDCC77", "High" = "#117733")
temp_shapes <- c("Heated" = 17, "Ambient" = 16)  # circle and triangle

# Base theme for all plots
base_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 7, margin = margin(t = 1)),
    axis.title.y = element_text(size = 7, margin = margin(r = 1)),
    axis.text.x = element_text(size = 6, margin = margin(t = 1)),
    axis.text.y = element_text(size = 6, margin = margin(r = 1)),
    legend.position = "none",
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
  )

##############################

#5).        Plot

##############################

#Calculate global ranges
x_limits <- range(ord_df$Axis.1)
y_limits <- range(ord_df$Axis.2)

#Equalized x and y axis between graphs
axis_range <- max(diff(x_limits), diff(y_limits)) / 2
x_mid <- mean(x_limits)
y_mid <- mean(y_limits)
x_limits <- c(x_mid - axis_range, x_mid + axis_range)
y_limits <- c(y_mid - axis_range, y_mid + axis_range)

# Function to plot Bray-Curtis
bray_curtis_plot <- function(taxonomy, marker) {
  df_filtered <- dplyr::filter(ord_df, Method == taxonomy, Metric == marker)
  
  # Grab axis labels
  x_lab <- unique(df_filtered$Axis1_Label)
  y_lab <- unique(df_filtered$Axis2_Label)
  
  p <- ggplot(
    df_filtered,
    aes(x = Axis.1, y = Axis.2, color = nutrient_level)
  ) +
    geom_point(size = 2.5) +
    scale_color_manual(values = nutrient_colors) +
    scale_x_continuous(
      labels = function(x) ifelse(x >= 0, paste0("\u2007", x), x)
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    labs(x = x_lab, y = y_lab, color = "Nutrient Enrichment") +
    coord_fixed(xlim = x_limits, ylim = y_limits)
  
  return(p)
}

#Run function to calculate Bray Curtis
p1 <- bray_curtis_plot("Class", "Metatranscriptomics")
p2 <- bray_curtis_plot("Family", "Metatranscriptomics")
p3 <- bray_curtis_plot("Genus", "Metatranscriptomics")
p4 <- bray_curtis_plot("Class", "Multi-Marker")
p5 <- bray_curtis_plot("Family", "Multi-Marker")
p6 <- bray_curtis_plot("Genus", "Multi-Marker")
p7 <- bray_curtis_plot("Class", "18S")
p8 <- bray_curtis_plot("Family", "18S")
p9 <- bray_curtis_plot("Genus", "18S")
p10 <- bray_curtis_plot("Class", "COI")
p11 <- bray_curtis_plot("Family", "COI")
p12 <- bray_curtis_plot("Genus", "COI")
p13 <- bray_curtis_plot("Class", "23S")
p14 <- bray_curtis_plot("Family", "23S")
p15 <- bray_curtis_plot("Genus", "23S")


# Get the Nutrient Enrichment Level
legend <- get_legend(
  p4 + theme(legend.position = "bottom", legend.text = element_text(size = 11),
             legend.title = element_text(size = 15))
)

# Combine plots into a 5x3 matrix, with minimal gaps between panels
plots_grid <- plot_grid(
  # Class row
  p1, p4, p7, p10, p13,
  # Family row
  p2, p5, p8, p11, p14,
  # Genus row
  p3, p6, p9, p12, p15,
  ncol = 5,
  align = "hv"
)

# Create left row labels
left_labels <- ggdraw() +
  draw_label("Class",  x = 0.5, y = 0.86, angle = 90, fontface = "bold", size = 14) +
  draw_label("Family", x = 0.5, y = 0.53, angle = 90, fontface = "bold", size = 14) +
  draw_label("Genus",  x = 0.5, y = 0.19, angle = 90, fontface = "bold", size = 14)

# Create top column labels
top_labels <- ggdraw() +
  draw_label("MT",  x = 0.153, y = 0.5, fontface = "bold", size = 14) +
  draw_label("MM",  x = 0.343, y = 0.5, fontface = "bold", size = 14) +
  draw_label("18S", x = 0.534, y = 0.5, fontface = "bold", size = 14) +
  draw_label("COI", x = 0.724, y = 0.5, fontface = "bold", size = 14) +
  draw_label("23S", x = 0.915, y = 0.5, fontface = "bold", size = 14)


#Put Legeend below the plot
plots_with_side_labels <- plot_grid(
  left_labels,
  plots_grid,
  ncol = 2,
  rel_widths = c(0.05, 1)
)

# Combine everything: top labels + (side labels + plots) + legend on its own row at bottom
final_plot <- plot_grid(
  top_labels,
  plots_with_side_labels,
  legend,
  ncol = 1,
  rel_heights = c(0.05, 1, 0.08)
)

# Display final plot
final_plot

##############################

#6).     Save Figure

##############################

# Save figure
ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Total_Bray_Curtis_plots.png",
  plot = final_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


##############################

#7).     PERMANOVA

##############################

#Function PERMANOVA
permanova <- function(df_results, metadata) {
  # Remove empty rows
  df_results <- df_results[rowSums(df_results, na.rm = TRUE) > 0, ]
  
  # Compute Bray-Curtis distance
  bray <- vegdist(df_results, method = "bray")
  
  # Run adonis2
  total_bray <- adonis2(bray ~ phos + temp, data = metadata, permutations = 999)
  print(total_bray)
  
  # Run marginal effects
  variable_bray <- adonis2(bray ~ phos + temp, data = metadata, permutations = 999, by = "margin")
  print(variable_bray)
}

permanova(df_result_class_MT, metadata_mt)
permanova(df_result_family_MT, metadata_mt)
permanova(df_result_genus_MT, metadata_mt)

permanova(df_result_class_MM, metadata_mt)
permanova(df_result_family_MM, metadata_mt)
permanova(df_result_genus_MM, metadata_mt)

permanova(df_result_class_18S, metadata_mt)
permanova(df_result_family_18S, metadata_mt)
permanova(df_result_genus_18S, metadata_mt)

permanova(df_result_class_COI, metadata_mt)
permanova(df_result_family_COI, metadata_mt)
permanova(df_result_genus_COI, metadata_mt)

permanova(df_result_class_23S, metadata_mt)
permanova(df_result_family_23S, metadata_mt)
permanova(df_result_genus_23S, metadata_mt)

##############################

#8).     Mantel Tests

##############################

# Remove G4 from MM, 18S, COI and 23S datasets
df_result_class_MM <- df_result_class_MM[rownames(df_result_class_MM) != "G4", ]
df_result_family_MM <- df_result_family_MM[rownames(df_result_family_MM) != "G4", ]
df_result_genus_MM <- df_result_genus_MM[rownames(df_result_genus_MM) != "G4", ]

df_result_class_18S<- df_result_class_18S[rownames(df_result_class_18S) != "G4", ]
df_result_family_18S <- df_result_family_18S[rownames(df_result_family_18S) != "G4", ]
df_result_genus_18S <- df_result_genus_18S[rownames(df_result_genus_18S) != "G4", ]

df_result_class_COI <- df_result_class_COI[rownames(df_result_class_COI) != "G4", ]
df_result_family_COI <- df_result_family_COI[rownames(df_result_family_COI) != "G4", ]
df_result_genus_COI <- df_result_genus_COI[rownames(df_result_genus_COI) != "G4", ]

df_result_class_23S <- df_result_class_23S[rownames(df_result_class_23S) != "G4", ]
df_result_family_23S <- df_result_family_23S[rownames(df_result_family_23S) != "G4", ]
df_result_genus_23S <- df_result_genus_23S[rownames(df_result_genus_23S) != "G4", ]

# Initialize an empty results table
mantel_summary <- data.frame(
  Level = character(),
  Method = character(),
  Mantel_r = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Function to run Mantel test for multiple distance methods
run_mantel <- function(df_mt, df_mm, level){
  df_mt <- as.matrix(df_mt)
  df_mm <- as.matrix(df_mm)
  
  # Bray-Curtis
  dist_bc_mt <- vegdist(df_mt, method = "bray")
  dist_bc_mm <- vegdist(df_mm, method = "bray")
  mantel_bc <- mantel(dist_bc_mt, dist_bc_mm, method = "spearman", permutations = 999)
  
  # Jaccard (convert to presence/absence)
  dist_j_mt <- vegdist((df_mt > 0) * 1, method = "jaccard")
  dist_j_mm <- vegdist((df_mm > 0) * 1, method = "jaccard")
  mantel_j <- mantel(dist_j_mt, dist_j_mm, method = "spearman", permutations = 999)
  
  # Print
  cat("\nMantel Test for", level, "- Bray-Curtis\n")
  print(mantel_bc)
  cat("\nMantel Test for", level, "- Jaccard\n")
  print(mantel_j)
  
  # Add to summary
  mantel_summary <<- rbind(mantel_summary, 
                           data.frame(Level = level, Method = "Bray-Curtis",
                                      Mantel_r = mantel_bc$statistic,
                                      P_value = mantel_bc$signif),
                           data.frame(Level = level, Method = "Jaccard",
                                      Mantel_r = mantel_j$statistic,
                                      P_value = mantel_j$signif))
}

# Run Mantel tests for all levels
run_mantel(df_result_class_MT, df_result_class_MM, "Class")
run_mantel(df_result_family_MT, df_result_family_MM, "Family")
run_mantel(df_result_genus_MT, df_result_genus_MM, "Genus")

run_mantel(df_result_class_MT, df_result_class_18S, "Class")
run_mantel(df_result_family_MT, df_result_family_18S, "Family")
run_mantel(df_result_genus_MT, df_result_genus_18S, "Genus")

run_mantel(df_result_class_MT, df_result_class_COI, "Class")
run_mantel(df_result_family_MT, df_result_family_COI, "Family")
run_mantel(df_result_genus_MT, df_result_genus_COI, "Genus")

run_mantel(df_result_class_MT, df_result_class_23S, "Class")
run_mantel(df_result_family_MT, df_result_family_23S, "Family")
run_mantel(df_result_genus_MT, df_result_genus_23S, "Genus")

# Show summary table
print(mantel_summary)

write.csv(mantel_summary, "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/mantel_summary.csv", row.names = FALSE)

##############################

#9).     Betadisper

##############################

####Betadisper function
beta_dispersion <- function(df_results, metadata, treatment_col, taxonomy, marker) {
  
  # Ensure sample order matches
  meta_sub <- metadata %>%
    filter(SampleID %in% rownames(df_results)) %>%
    arrange(match(SampleID, rownames(df_results)))
  
  # Bray–Curtis distance
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
    Taxonomy = taxonomy,
    Method = marker
  )
  
  return(list(
    betadisper = bd,
    permutest = bd_perm,
    distances = bd_df
  ))
}

#Run Betadisper
BD_class_MT <- beta_dispersion(df_result_class_MT, metadata_mt, treatment_col = "phos", taxonomy = "Class", marker = "Metatranscriptomics")
BD_family_MT <- beta_dispersion(df_result_family_MT, metadata_mt, treatment_col = "phos", taxonomy = "Family", marker = "Metatranscriptomics")
BD_genus_MT <- beta_dispersion(df_result_genus_MT, metadata_mt, treatment_col = "phos", taxonomy = "Genus", marker = "Metatranscriptomics")

BD_class_MM <- beta_dispersion(df_result_class_MM, metadata_mt, treatment_col = "phos", taxonomy = "Class", marker = "Multi-Marker")
BD_family_MM <- beta_dispersion(df_result_family_MM, metadata_mt, treatment_col = "phos", taxonomy = "Family", marker = "Multi-Marker")
BD_genus_MM <- beta_dispersion(df_result_genus_MM, metadata_mt, treatment_col = "phos", taxonomy = "Genus", marker = "Multi-Marker")

BD_class_18S <- beta_dispersion(df_result_class_18S, metadata, treatment_col = "phos", taxonomy = "Class", marker = "18S")
BD_family_18S <- beta_dispersion(df_result_family_18S, metadata, treatment_col = "phos", taxonomy = "Family", marker = "18S")
BD_genus_18S <- beta_dispersion(df_result_genus_18S, metadata, treatment_col = "phos", taxonomy = "Genus", marker = "18S")

BD_class_COI <- beta_dispersion(df_result_class_COI, metadata, treatment_col = "phos", taxonomy = "Class", marker = "COI")
BD_family_COI <- beta_dispersion(df_result_family_COI, metadata, treatment_col = "phos", taxonomy = "Family", marker = "COI")
BD_genus_COI <- beta_dispersion(df_result_genus_COI, metadata, treatment_col = "phos", taxonomy = "Genus", marker = "COI")

BD_class_23S <- beta_dispersion(df_result_class_23S, metadata, treatment_col = "phos", taxonomy = "Class", marker = "23S")
BD_family_23S <- beta_dispersion(df_result_family_23S, metadata, treatment_col = "phos", taxonomy = "Family", marker = "23S")
BD_genus_23S <- beta_dispersion(df_result_genus_23S, metadata, treatment_col = "phos", taxonomy = "Genus", marker = "23S")


#Summarize results function
summarize_betadisper <- function(bd_obj) {
  bd_obj$distances %>%
    group_by(Treatment, Method, Taxonomy) %>%
    summarise(
      mean_distance = mean(Distance_to_Centroid),
      sd_distance   = sd(Distance_to_Centroid),
      .groups = "drop"
    )
}

##############################

#10).     Plot

##############################

# Method order and labels
method_levels <- c("Metatranscriptomics", "Multi-Marker", "18S", "COI", "23S")
method_labels <- c("MT", "MM", "18S", "COI", "23S")

BD_all <- BD_all %>%
  mutate(
    Treatment = as.character(Treatment),
    Treatment = factor(Treatment,
                       levels = c("40", "280"),
                       labels = c("Moderate", "High"))
  )

nutrient_colors <- c("Moderate" = "#DDCC77", "High" = "#117733")

BD_all_raw <- bind_rows(
  BD_class_MT$distances,
  BD_family_MT$distances,
  BD_genus_MT$distances,
  BD_class_MM$distances,
  BD_family_MM$distances,
  BD_genus_MM$distances,
  BD_class_18S$distances,
  BD_family_18S$distances,
  BD_genus_18S$distances,
  BD_class_COI$distances,
  BD_family_COI$distances,
  BD_genus_COI$distances,
  BD_class_23S$distances,
  BD_family_23S$distances,
  BD_genus_23S$distances
)

BD_all_raw <- BD_all_raw %>%
  mutate(
    Treatment = as.character(Treatment),
    Treatment = factor(
      Treatment,
      levels = c("40", "280"),
      labels = c("Moderate", "High")
    ),
    Method = factor(
      Method,
      levels = c("Metatranscriptomics", "Multi-Marker", "18S", "COI", "23S"),
      labels = c("MT", "MM", "18S", "COI", "23S")
    ),
    Taxonomy = factor(Taxonomy, levels = c("Class", "Family", "Genus"))
  )

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
  facet_wrap(
    ~ Taxonomy,
    ncol = 1,
    scales = "free_y",
    strip.position = "left"
  ) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 90,
      size = 14,
      face = "bold",
      margin = margin(r = 10)
    ),
    
    axis.title.x = element_blank(),
    axis.title.y.left = element_blank(),
    axis.title.y.right = element_text(
      size = 14,
      margin = margin(l = 10)
    ),
    
    axis.text.x.top = element_text(
      size = 13,
      face = "bold",
      margin = margin(b = 8)
    ),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    
    axis.text.y.left = element_blank(),
    axis.text.y.right = element_text(
      size = 12,
      angle = 0
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.spacing = unit(1.2, "lines"),
    
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
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/BD_BetaDispersion_boxplot.png",
  plot = beta_disp_plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

