########Beta Diversity
# This code with create PcOA for Jaccard based on 18S, COI and 23S at the Class, Family and Genus level
# Specifically looking at NovaSeq vs MiSeq

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
metadata_18S <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata_TR.tsv", sep="\t", header=TRUE, comment.char="")
metadata_COI <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata_TR_COI.tsv", sep="\t", header=TRUE, comment.char="")
metadata_23S <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata_TR_23S.tsv", sep="\t", header=TRUE, comment.char="")

#Load 18S ASV class
filepath_18S_c <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/18S_Class.tsv"
df_18S_c <- read_tsv(filepath_18S_c)

#Load 18S ASV family
filepath_18S_f <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/18S_Family.tsv"
df_18S_f <- read_tsv(filepath_18S_f)

#Load 18S ASV class
filepath_18S_g <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/18S_Genus.tsv"
df_18S_g <- read_tsv(filepath_18S_g)


#Load COI ASV
filepath_COI_c <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/COI_Class.tsv"
df_COI_c <- read_tsv(filepath_COI_c)

#Load COI ASV
filepath_COI_f <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/COI_Family.tsv"
df_COI_f <- read_tsv(filepath_COI_f)

#Load COI ASV
filepath_COI_g <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/COI_Genus.tsv"
df_COI_g <- read_tsv(filepath_COI_g)


#Load 23S ASV
filepath_23S_c <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/23S_Class.tsv"
df_23S_c <- read_tsv(filepath_23S_c)

#Load 23S ASV
filepath_23S_f <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/23S_Family.tsv"
df_23S_f <- read_tsv(filepath_23S_f)

#Load 23S ASV
filepath_23S_g <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/All_Replicates/23S_Genus.tsv"
df_23S_g <- read_tsv(filepath_23S_g)

############################

#2)     Process Files

############################

#Function to process files
process_df <- function(df_taxon) {
  # Columns to remove
  cols_to_remove <- c("False_Positives", "Non-Target_Taxa", "Reason", 
                      "Mean", "DB_R1", "DB_R2", "DB_R3", "EB_21", "Confidence")
  
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
df_result_class_18S <- process_df(df_18S_c)
df_result_family_18S <- process_df(df_18S_f)
df_result_genus_18S <- process_df(df_18S_g)
df_result_class_COI <- process_df(df_COI_c)
df_result_family_COI <- process_df(df_COI_f)
df_result_genus_COI <- process_df(df_COI_g)
df_result_class_23S <- process_df(df_23S_c)
df_result_family_23S <- process_df(df_23S_f)
df_result_genus_23S <- process_df(df_23S_g)

##############################

#3) Calculating Beta Diversity

##############################

# Function to calculate presence-absence Jaccard diversity and PCoA
jacc <- function(df_results, taxonomy, marker) {
  # 1) Convert to numeric presence/absence
  df_pa <- as.data.frame(df_results > 0) * 1  # converts TRUE/FALSE to 1/0
  
  # 2) Jaccard distance
  jacc_dist <- vegdist(df_pa, method = "jaccard", binary = TRUE)
  
  # 3) PCoA
  jacc_pcoa <- cmdscale(jacc_dist, k = 2, eig = TRUE)
  
  # 4) Variance explained
  jacc_var_exp <- round(jacc_pcoa$eig / sum(jacc_pcoa$eig) * 100, 1)
  
  # 5) Extract coordinates
  jacc_pcoa_df <- as.data.frame(jacc_pcoa$points) %>%
    rename(Axis.1 = V1, Axis.2 = V2) %>%
    mutate(
      Method = taxonomy,
      Metric = marker,
      Axis1_Label = paste0("Axis 1 (", jacc_var_exp[1], "%)"),
      Axis2_Label = paste0("Axis 2 (", jacc_var_exp[2], "%)"),
      SampleID = rownames(.)
    )
  
  return(jacc_pcoa_df)
}

#Run function to calculate Jaccard
J_class_18S <- jacc(df_result_class_18S, "Class", "18S")
J_family_18S <- jacc(df_result_family_18S, "Family", "18S")
J_genus_18S <- jacc(df_result_genus_18S, "Genus", "18S")
J_class_COI <- jacc(df_result_class_COI, "Class", "COI")
J_family_COI <- jacc(df_result_family_COI, "Family", "COI")
J_genus_COI <- jacc(df_result_genus_COI, "Genus", "COI")
J_class_23S <- jacc(df_result_class_23S, "Class", "23S")
J_family_23S <- jacc(df_result_family_23S, "Family", "23S")
J_genus_23S <- jacc(df_result_genus_23S, "Genus", "23S")

##############################

#4).     Prepare Plot

##############################

# Combine all ordination results
ord_df <- bind_rows(
  J_class_18S, J_family_18S, J_genus_18S,
  J_class_COI, J_family_COI, J_genus_COI,
  J_class_23S, J_family_23S, J_genus_23S
)

# Merge with metadata so only existing SampleIDs are kept
ord_df <- left_join(ord_df, metadata_COI, by = "SampleID")

#Color Map
color_map <- c(
  "F2" = "#E69F00",
  "F6" = "#56B4E9",
  "F7" = "#009E73",
  "G2" = "#F0E442",
  "G4" = "#0072B2",
  "G5" = "#D55E00",
  "G7" = "#CC79A7",
  "J1" = "#117733",
  "J3" = "#999933",
  "J5" = "#DDCC77",
  "J8" = "#661100",
  "K5" = "#999999",
  "blank" = "#000000"
)

#Merge with ord_df again if metadata_test was changed
ord_df <- left_join(ord_df, metadata_COI, by = c("SampleID" = "SampleID"))

# Base theme for all plots
base_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"  # show legend
  )

# Function to plot Jaccard PCoA with custom color mapping
jacc_curtis_plot <- function(df_ord, metadata, taxonomy, marker, color_map = NULL) {
  
  # Filter ordination results
  df_filtered <- df_ord %>% 
    filter(Method == taxonomy, Metric == marker)
  
  # Merge with metadata
  df_filtered <- left_join(df_filtered, metadata, by = "SampleID")
  
  # Check if Sample column exists
  if(!"Sample" %in% colnames(df_filtered)) {
    stop("Column 'Sample' not found in metadata or join failed.")
  }
  
  # Use provided color_map if given, otherwise generate automatically
  if(is.null(color_map)){
    all_samples <- unique(df_filtered$Sample)
    n_samples <- length(all_samples)
    
    if(n_samples >= 3 & n_samples <= 12){
      colors <- RColorBrewer::brewer.pal(n_samples, "Set3")
    } else {
      colors <- grDevices::rainbow(n_samples)
    }
    color_map <- setNames(colors, all_samples)
  }
  
  # Grab axis labels
  x_lab <- unique(df_filtered$Axis1_Label)
  y_lab <- unique(df_filtered$Axis2_Label)
  
  # Plot
  p <- ggplot(df_filtered, aes(x = Axis.1, y = Axis.2,
                               color = Sample, shape = SampleType)) +
    geom_point(size = 3) +
    scale_color_manual(values = color_map) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right"
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      color = "Sample",
      shape = "Sequencing Type"
    )
  
  return(p)
}

#Run function to calculate Jaccard
p1 <- jacc_curtis_plot(ord_df, metadata_18S, "Class", "18S", color_map)
p2 <- jacc_curtis_plot(ord_df, metadata_18S, "Family", "18S", color_map)
p3 <- jacc_curtis_plot(ord_df, metadata_18S, "Genus", "18S", color_map)
p4 <- jacc_curtis_plot(ord_df, metadata_COI, "Class", "COI", color_map)
p5 <- jacc_curtis_plot(ord_df, metadata_COI, "Family", "COI", color_map)
p6 <- jacc_curtis_plot(ord_df, metadata_COI, "Genus", "COI", color_map)
p7 <- jacc_curtis_plot(ord_df, metadata_23S, "Class", "23S", color_map)
p8 <- jacc_curtis_plot(ord_df, metadata_23S, "Family", "23S", color_map)
p9 <- jacc_curtis_plot(ord_df, metadata_23S, "Genus", "23S", color_map)

# Extract legend from p5, with larger title and text
legend <- get_legend(
  p5 + theme(
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 11)
  )
)

# Remove legends from individual plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")

# Combine plots into a 3x3 matrix: columns = marker (18S, COI, 23S), rows = taxonomic level
plots_grid <- plot_grid(
  p1, p4, p7,   # Class row: 18S, COI, 23S
  p2, p5, p8,   # Family row: 18S, COI, 23S
  p3, p6, p9,   # Genus row: 18S, COI, 23S
  ncol = 3,
  align = "hv"
)

# Top labels now mark marker (columns)
top_labels <- ggdraw() +
  draw_label("18S", x = 0.2, y = 0.7, fontface = "bold", size = 14) +
  draw_label("COI", x = 0.48, y = 0.7, fontface = "bold", size = 14) +
  draw_label("23S", x = 0.76, y = 0.7, fontface = "bold", size = 14)

# Left labels now mark taxonomic level (rows)
left_labels <- ggdraw() +
  draw_label("Class", x = 0.5, y = 0.855, angle = 90, fontface = "bold", size = 14) +
  draw_label("Family", x = 0.5, y = 0.53, angle = 90, fontface = "bold", size = 14) +
  draw_label("Genus", x = 0.5, y = 0.2, angle = 90, fontface = "bold", size = 14)

# Combine everything: labels + plots + shared legend
final_plot <- plot_grid(
  top_labels,
  plot_grid(
    left_labels,
    plots_grid,
    legend,
    ncol = 3,
    rel_widths = c(0.05, 1, 0.15)
  ),
  ncol = 1,
  rel_heights = c(0.07, 1)
)

# Display final plot
final_plot

##############################

#5).     Save Figure

##############################

# Save the figure
ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/MiSeq_vs_NovaSeq_plot.png",
  plot = final_plot,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)
