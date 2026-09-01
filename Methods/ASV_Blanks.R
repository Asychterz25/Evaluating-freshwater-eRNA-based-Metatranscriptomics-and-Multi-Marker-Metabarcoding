############## Blank ASVs
#Creates Jaccard PCoA graphs of ASVs of 18S, COI and 23S
#Done to show 1) Blanks significantly different than samples, 2) No difference between MiSeq and NovaSeq

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

#Load 18S ASV
filepath_18S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/ASV/18S_ASV.tsv"
df_18S <- read_tsv(filepath_18S)

#Load COI ASV
filepath_COI <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/ASV/COI_ASV.tsv"
df_COI <- read_tsv(filepath_COI)

#Load 23S ASV
filepath_23S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/ASV/23S_ASV.tsv"
df_23S <- read_tsv(filepath_23S)

############################

#2)     Process Files

############################

process_df <- function(df_taxon) {
  # Store ASV IDs separately
  taxa <- df_taxon$ASV_ID
  
  # Keep only numeric columns for calculations
  df_numeric <- df_taxon[, !(names(df_taxon) %in% c("ASV_ID", "Taxon"))]
  
  # Ensure numeric type
  df_numeric <- as.data.frame(lapply(df_numeric, as.numeric))
  
  # Avoid division by zero
  col_sums <- colSums(df_numeric, na.rm = TRUE)
  col_sums[col_sums == 0] <- NA
  
  # Convert to relative abundance
  df_rel <- sweep(df_numeric, 2, col_sums, FUN = "/") * 100
  df_rel[is.na(df_rel)] <- 0
  
  # Assign ASV IDs as row names
  rownames(df_rel) <- taxa
  
  # Transpose (samples as rows, ASVs as columns)
  df_transposed <- t(df_rel)
  
  return(df_transposed)
}

#Run function to clean files
df_result_18S <- process_df(df_18S)
df_result_COI <- process_df(df_COI)
df_result_23S <- process_df(df_23S)

##############################

#3) Calculating Beta Diversity

##############################

# Function to calculate presence-absence Jaccard diversity and PCoA
jacc <- function(df_results) {
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
      Axis1_Label = paste0("Axis 1 (", jacc_var_exp[1], "%)"),
      Axis2_Label = paste0("Axis 2 (", jacc_var_exp[2], "%)"),
      SampleID = rownames(.)
    )
  
  return(jacc_pcoa_df)
}

#Run function to calculate Jaccard
J_18S <- jacc(df_result_18S)
J_COI <- jacc(df_result_COI)
J_23S <- jacc(df_result_23S)

##############################

#4).     Prepare Plot

##############################

#Color Map
color_map <- c(
  "blank" = "#000000",
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
  "K5" = "#999999"
)


# General plotting function for Jaccard PCoA with custom colors
jacc_plot <- function(ord_df, metadata, marker_name, color_map = NULL) {
  
  # Merge ordination results with metadata
  df_merged <- ord_df %>%
    left_join(metadata, by = "SampleID")
  
  # Extract axis labels
  x_lab <- unique(df_merged$Axis1_Label)
  y_lab <- unique(df_merged$Axis2_Label)
  
  # Base plot
  p <- ggplot(df_merged, aes(x = Axis.1, y = Axis.2, 
                             color = Sample, shape = SampleType)) +
    geom_point(size = 3) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "right"
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      color = "Sample",
      shape = "Sequencing Type"
    )
  
  # Apply custom colors if provided
  if (!is.null(color_map)) {
    p <- p + scale_color_manual(values = color_map)
  }
  
  return(p)
}

#Run plotting function
p18S <- jacc_plot(J_18S, metadata_18S, "18S", color_map)
pCOI <- jacc_plot(J_COI, metadata_COI, "COI", color_map)
p23S <- jacc_plot(J_23S, metadata_23S, "23S", color_map)

# Display plots
p18S
pCOI
p23S


##############################

#5).     Merge Figure

##############################

# Save the figure
ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plots/18S_total_Jaccard.png",
  plot = p18S,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

##############################

#6).     PERMANOVA

##############################

# PERMANOVA
permanova <- function(df_results, metadata) {
  # Remove empty rows
  df_results <- df_results[rowSums(df_results, na.rm = TRUE) > 0, ]
  
  # Compute Jaccard distance (presence/absence)
  jacc <- vegdist(df_results, method = "jacc")
  
  # Run PERMANOVA for SampleType
  sampletype_jacc <- adonis2(jacc ~ SampleType, data = metadata, permutations = 999)
  print(sampletype_jacc)
}

permanova(df_result_18S, metadata_18S)
permanova(df_result_COI, metadata_COI)
permanova(df_result_23S, metadata_23S)

##############################

#7).     PERMANOVA within each nutrient treatment

##############################

# PERMANOVA within each treatment
permanova_by_treatment <- function(df_results, metadata, treatment_col = "phos") {
  
  # Remove empty rows
  df_results <- df_results[rowSums(df_results, na.rm = TRUE) > 0, ]
  
  # Check that SampleIDs match
  metadata <- metadata[metadata$SampleID %in% rownames(df_results), ]
  metadata <- metadata[match(rownames(df_results), metadata$SampleID), ]
  
  # Get all unique treatments
  treatments <- unique(metadata[[treatment_col]])
  
  # Loop over treatments
  for (tr in treatments) {
    cat("\n===== PERMANOVA for", treatment_col, "=", tr, "=====\n")
    
    # Subset metadata and ASV table
    metadata_sub <- metadata[metadata[[treatment_col]] == tr, ]
    df_sub <- df_results[rownames(df_results) %in% metadata_sub$SampleID, ]
    
    # Align rows
    metadata_sub <- metadata_sub[match(rownames(df_sub), metadata_sub$SampleID), ]
    
    # Compute Jaccard distance
    jacc <- vegdist(df_sub, method = "jacc")
    
    # Run PERMANOVA for SampleType
    result <- adonis2(jacc ~ SampleType, data = metadata_sub, permutations = 999)
    print(result)
  }
}

permanova_by_treatment(df_result_18S, metadata_18S)
permanova_by_treatment(df_result_COI, metadata_COI)
permanova_by_treatment(df_result_23S, metadata_23S)

##############################

# 8).     PERMDISP within each nutrient treatment

##############################

# PERMDISP within each treatment
permdisp_by_treatment <- function(df_results, metadata, treatment_col = "phos") {
  
  # Remove empty rows
  df_results <- df_results[rowSums(df_results, na.rm = TRUE) > 0, ]
  
  # Ensure SampleIDs match
  metadata <- metadata[metadata$SampleID %in% rownames(df_results), ]
  metadata <- metadata[match(rownames(df_results), metadata$SampleID), ]
  
  # Get all unique treatments
  treatments <- unique(metadata[[treatment_col]])
  
  # Loop over treatments
  for (tr in treatments) {
    cat("\n===== PERMDISP for", treatment_col, "=", tr, "=====\n")
    
    # Subset metadata and ASV table
    metadata_sub <- metadata[metadata[[treatment_col]] == tr, ]
    df_sub <- df_results[rownames(df_results) %in% metadata_sub$SampleID, ]
    
    # Align rows
    metadata_sub <- metadata_sub[match(rownames(df_sub), metadata_sub$SampleID), ]
    
    # Compute Jaccard distance
    jacc <- vegdist(df_sub, method = "jacc")
    
    # Compute multivariate dispersion (PERMDISP)
    disp <- betadisper(jacc, metadata_sub$SampleType)
    
    # Run permutation test
    result <- permutest(disp, permutations = 999)
    
    # Print results
    print(result)
    
    # Optional: show mean distances to centroid
    cat("\nMean distances to centroid:\n")
    print(tapply(disp$distances, metadata_sub$SampleType, mean))
  }
}

# Run for each marker
permdisp_by_treatment(df_result_18S, metadata_18S)
permdisp_by_treatment(df_result_COI, metadata_COI)
permdisp_by_treatment(df_result_23S, metadata_23S)

