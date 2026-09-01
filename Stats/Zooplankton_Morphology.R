##################### Comparisons with MOology

#This Code will compare Zooplankton counts with MT, MM, COI at the family level

#Load libraries
library(ggVennDiagram)
library(readr)
library(dplyr)
library(vegan)

############################

#1)     Load Files

############################
#Load metadata
metadata <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata.tsv", sep="\t", header=TRUE, comment.char="")

#Load MT metadata
metadata_mt <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata_mt.tsv", sep="\t", header=TRUE, comment.char="")

######## MT
#Load MT Metazoans identified through Family
filepath_mt <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/MT_metazoan.tsv"
df_mt <- read_tsv(filepath_mt)

######## MM
#Load MM Metazoans identified through Family
filepath_mm <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/MM_metazoan.tsv"
df_mm <- read_tsv(filepath_mm)

######## COI
#Load COI Metazoans identified through Family
filepath_COI <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/COI_metazoan.tsv"
df_COI <- read_tsv(filepath_COI)

######## MO
#Load MOology Metazoans identified at Family
filepath_MO <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/Morphology_family.tsv"
df_MO <- read_tsv(filepath_MO)


############################

#2)     Process Files for Venn Diagram

############################

#Function to process files
process_df <- function(df_taxon) {
  # Columns to remove
  cols_to_remove <- c("False_Positives", "Non-Target_Taxa", "Reason", "Mean", "DB_R1", "DB_R2", "DB_R3", "EB_21", "G4")
  
  # 1) Remove rows where "Reason" is not blank
  df_clean <- df_taxon[df_taxon$Reason == "" | is.na(df_taxon$Reason), ]
  
  # 2) Remove unwanted columns
  df_clean <- df_clean[, !(names(df_clean) %in% cols_to_remove)]
  
  # 3) Save Taxon column separately
  taxa <- df_clean$Taxon
  
  # 4) Compute relative abundance
  df_clean_numeric <- df_clean[, !(names(df_clean) %in% c("Taxon"))]
  col_sums <- colSums(df_clean_numeric, na.rm = TRUE)
  col_sums[col_sums == 0] <- NA
  
  df_rel <- sweep(df_clean_numeric, 2, col_sums, FUN = "/") * 100
  df_rel[is.na(df_rel)] <- 0
  
  # 5) Add back the Taxon column
  df_rel$Taxon <- taxa
  
  return(df_rel)
}

#Run function to clean files
df_result_venn_MT <- process_df(df_mt)
df_result_venn_MM <- process_df(df_mm)
df_result_venn_MO <- process_df(df_MO)


################################################# Venn Diagram ########################################################

############################

#3)  Calculate shared Taxonomy (MT, MM, MO) and graph

############################

#Function to calculate the unique taxon from MT, 18S (total), COI, 23S(algae)
unique_taxon <- function(MT, MM, MO, plot_name) {
  taxon_col <- function(df) {
    col <- grep("taxon", names(df), ignore.case = TRUE, value = TRUE)
    if (length(col) == 0) stop("No column named 'taxon' found in one of the inputs.")
    unique(trimws(as.character(df[[col[1]]])))
  }
  
  sets <- list(
    "MT" = taxon_col(MT),
    "MM" = taxon_col(MM),
    "MO" = taxon_col(MO)
  )
  
  sets <- sets[sapply(sets, length) > 0]
  
  p <- ggVennDiagram(sets, label_alpha = 0, label = "count", category.names = names(sets)) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme_void() +
    theme(legend.position = "none")
  
  # assign the plot to the name given in plot_name
  assign(plot_name, p, envir = .GlobalEnv)
}

venn_diagram_MO <- unique_taxon(df_result_venn_MT, df_result_venn_MM, df_result_venn_MO, "MOology Venn Diagram")

venn_diagram_MO

#Read Taxon
taxon_col <- function(df) {
  col <- grep("taxon", names(df), ignore.case = TRUE, value = TRUE)
  unique(trimws(as.character(df[[col[1]]])))
}

mt_taxa   <- taxon_col(df_result_venn_MT)
mm_taxa   <- taxon_col(df_result_venn_MM)
MO_taxa <- taxon_col(df_result_venn_MO)

shared_all <- Reduce(intersect, list(mt_taxa, mm_taxa, MO_taxa))
shared_mt_mm <- intersect(mt_taxa, mm_taxa)
shared_mt_MO <- intersect(mt_taxa, MO_taxa)
shared_mm_MO <- intersect(mm_taxa, MO_taxa)
unique_mt    <- setdiff(mt_taxa, union(mm_taxa, MO_taxa))
unique_mm    <- setdiff(mm_taxa, union(mt_taxa, MO_taxa))
unique_MO <- setdiff(MO_taxa, union(mt_taxa, mm_taxa))


shared_all
shared_mt_mm
shared_mt_MO
shared_mm_MO

unique_mt
unique_mm
unique_MO


############################

#4)  Process Files for Beta Diversity 

############################

# Function to process files
process_df_BC <- function(df_taxon) {
  # Columns to remove
  cols_to_remove <- c("False_Positives", "Non-Target_Taxa", "Reason", 
                      "Mean", "DB_R1", "DB_R2", "DB_R3", "EB_21")
  
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

#Remove column G4 from MOology for comparison with MT
df_MO_MT <- df_MO[, !names(df_MO) %in% "G4"]

#Run function to clean files
df_result_BC_MT <- process_df_BC(df_mt)
df_result_BC_MM <- process_df_BC(df_mm)
df_result_BC_COI <- process_df_BC(df_COI)
df_result_BC_MO <- process_df_BC(df_MO)
df_result_BC_MO_MT <- process_df_BC(df_MO_MT)


# Function to calculate Bray-Curtis diversity and PCoA (marker-only)
bray_curtis <- function(df_results, marker) {
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
      Marker = marker,
      Axis1_Label = paste0("Axis 1 (", bray_var_exp[1], "%)"),
      Axis2_Label = paste0("Axis 2 (", bray_var_exp[2], "%)"),
      SampleID = rownames(.)
    )
  
  # 5) Return the result
  return(bray_pcoa_df)
}

# Run function to calculate Bray-Curtis
BC_MT    <- bray_curtis(df_result_BC_MT, "Metatranscriptomics")
BC_MM    <- bray_curtis(df_result_BC_MM, "Multi-Marker")
BC_COI   <- bray_curtis(df_result_BC_COI, "COI")
BC_MO <- bray_curtis(df_result_BC_MO, "MOology")


############################

#5)  Plot Beta Diversity

############################

#Change High to Heated in the temp column
metadata_new <- metadata %>%
  mutate(temp = recode(temp,
                       "High" = "Heated",
                       "Ambient" = "Ambient"))

# Assuming you've run the marker-only bray_curtis() function already:
ord_df <- bind_rows(
  BC_MT,
  BC_MM,
  BC_COI,
  BC_MO
)

# Merge metadata
ord_df <- left_join(ord_df, metadata_new, by = "SampleID") %>%
  rename(
    temp_level = temp,
    nutrient_level = phos
  )

ord_df <- ord_df %>%
  mutate(
    nutrient_level = factor(nutrient_level, levels = c(40, 280), labels = c("Moderate", "High"))
  )

# Define aesthetics
nutrient_colors <- c("Moderate" = "#DDCC77", "High" = "#117733")
temp_shapes <- c("Heated" = 17, "Ambient" = 16) 

# Base theme
base_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

#Determine global axis limits
x_limits <- range(ord_df$Axis.1, na.rm = TRUE)
y_limits <- range(ord_df$Axis.2, na.rm = TRUE)

#Function to plot Bray-Curtis PCoA by marker
bray_curtis_plot <- function(marker, title_name) {
  df_filtered <- dplyr::filter(ord_df, Marker == marker)
  
  # Grab labels directly from data
  x_lab <- unique(df_filtered$Axis1_Label)
  y_lab <- unique(df_filtered$Axis2_Label)
  
  ggplot(df_filtered, aes(x = Axis.1, y = Axis.2,
                          color = nutrient_level, shape = temp_level)) +
    geom_point(size = 3) +
    scale_color_manual(values = nutrient_colors) +
    scale_shape_manual(values = temp_shapes) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    base_theme +
    labs(
      x = x_lab,
      y = NULL,  # remove y-axis label
      color = "Nutrient Level",
      shape = "Temperature",
      title = title_name
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
}

########################################

# 6) Generate plots for each marker

########################################

p1 <- bray_curtis_plot("Metatranscriptomics", "MT")
p2 <- bray_curtis_plot("Multi-Marker", "MM")
p3 <- bray_curtis_plot("COI", "COI")
p4 <- bray_curtis_plot("MOology", "MO")

# Extract legend with larger text
legend <- get_legend(
  p1 + theme(
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.box = "horizontal",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 11)
  ) +
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    )
)

########################################

# 7) Combine in a 1×4 row

########################################

# Combine Beta Diversity plots (1 row, 4 cols)
plots_grid <- plot_grid(p1, p2, p3, p4, ncol = 4, align = "hv")

# Combine with legend on the right
final_plot <- plot_grid(
  plots_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.15) 
)
final_plot

############################

#8)  Combine Beta Diversity and Venn Diagram

############################

# Venn diagram plot is `venn_diagram_MO`
# Beta Diversity plots are `final_plot` (from previous step)
combined_plot <- plot_grid(
  venn_diagram_MO,
  final_plot,
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("A)", "B)"),
  label_size = 18
)

combined_plot

ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Morphology_Venn_Beta.png",
  plot = combined_plot,      
  width = 12, 
  height = 8,
  dpi = 300   
)


############################

#9)  Mantel Tests

############################

# Bray-Curtis distance matrices
dist_MT    <- vegdist(df_result_BC_MT, method = "bray")
dist_MM    <- vegdist(df_result_BC_MM, method = "bray")
dist_COI   <- vegdist(df_result_BC_COI, method = "bray")
dist_MO <- vegdist(df_result_BC_MO, method = "bray")
dist_MO_MT <- vegdist(df_result_BC_MO_MT, method = "bray")

# Morphology vs Multi-Marker
mantel_MM <- mantel(dist_MO, dist_MM, method = "spearman", permutations = 999)
mantel_MM

# Morphology vs COI
mantel_COI <- mantel(dist_MO, dist_COI, method = "spearman", permutations = 999)
mantel_COI

# Morphology vs Metatranscriptomics
mantel_MT <- mantel(dist_MO_MT, dist_MT, method = "spearman", permutations = 999)
mantel_MT

mantel_results <- data.frame(
  Comparison = c("MO vs MT", "MO vs MM", "MO vs COI"),
  Mantel_r   = c(mantel_MT$statistic, mantel_MM$statistic, mantel_COI$statistic),
  p_value    = c(mantel_MT$signif, mantel_MM$signif, mantel_COI$signif)
)

mantel_results

