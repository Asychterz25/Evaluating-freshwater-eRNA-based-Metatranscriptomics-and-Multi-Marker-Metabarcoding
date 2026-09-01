######## Venn Diagrams
# This code with create venn diagrams showing which taxa are unique and shared between methods
# MT = Metatranscriptomics
# MM = Multi-marker metabarcoding
#Removes G4 from calculations

#Load libraries
library(vegan)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggVennDiagram)
library(ggforce)

############################

#1)     Load Files

############################
#Load MT metadata (only used in Spearman Rank)
metadata_mt <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata_MT.tsv", sep="\t", header=TRUE, comment.char="")

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

#Function to process files
process_df <- function(df_taxon) {
  # Columns to remove
  cols_to_remove <- c("False_Positives", "Non-Target_Taxa", "Reason", "Mean", "DB_R1", "DB_R2", "DB_R3", "EB_21", "G4", "Marker")
  
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

################################################# Total Venn Diagram ########################################################

############################

#3)  Calculate shared Taxonomy (MT, 18S, COI, 23S) and graph

############################

#Function to calculate the unique taxon from MT, 18S (total), COI, 23S(algae)
unique_taxon <- function(MT, total, COI, algae, plot_name) {
  taxon_col <- function(df) {
    col <- grep("taxon", names(df), ignore.case = TRUE, value = TRUE)
    if (length(col) == 0) stop("No column named 'taxon' found in one of the inputs.")
    unique(trimws(as.character(df[[col[1]]])))
  }
  
  sets <- list(
    "MT" = taxon_col(MT),
    "18S" = taxon_col(total),
    "COI" = taxon_col(COI),
    "23S" = taxon_col(algae)
  )
  
  sets <- sets[sapply(sets, length) > 0]
  
  p <- ggVennDiagram(sets, label_alpha = 0, label = "count", category.names = names(sets)) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme_void() +
    theme(legend.position = "none")
  
  # assign the plot to the name given in plot_name
  assign(plot_name, p, envir = .GlobalEnv)
}

class <- unique_taxon(df_result_class_MT, df_result_class_18S, df_result_class_COI, df_result_class_23S, "class_p")
family <- unique_taxon(df_result_family_MT, df_result_family_18S, df_result_family_COI, df_result_family_23S, "family_p")
genus <- unique_taxon(df_result_genus_MT, df_result_genus_18S, df_result_genus_COI, df_result_genus_23S, "genus_p")

############################

#4)  Combine graphs for Venn Diagram

############################

# Combine plots into a 4x3 matrix
plots_grid <- plot_grid(
  class_p, family_p, genus_p,
  ncol = 3,
  align = "hv",
  scale = 1
)

# Create top column labels
top_labels <- ggdraw() +
  draw_label("Class", x = 0.17, y = 0.21, fontface = "bold", size = 14) +
  draw_label("Family", x = 0.5, y = 0.21, fontface = "bold", size = 14) +
  draw_label("Genus", x = 0.83, y = 0.21, fontface = "bold", size = 14)

# Combine everything: labels + plots + legend
final_plot <- plot_grid(
  top_labels,
  plots_grid,
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
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Venn_diagram_supplementary_total.png",
  plot = final_plot,
  width = 12,
  height = 4,
  dpi = 300,
  bg = "white"
)


######################################## MT vs MM Venn diagram + Spearman ##################################################

############################

#3)  Calculate shared Taxonomy (MT, MM) and graph

############################

#Function to calculate the unique taxon from MT and MM and create a venn diagram
unique_taxon_ggplot <- function(MT, MM, plot_name = "venn_plot") {
  
  taxon_col <- function(df) {
    col <- grep("taxon", names(df), ignore.case = TRUE, value = TRUE)
    if (length(col) == 0) stop("No column named 'taxon' found in one of the inputs.")
    unique(trimws(as.character(df[[col[1]]])))
  }
  
  set_MT <- taxon_col(MT)
  set_MM <- taxon_col(MM)
  intersect_taxa <- intersect(set_MT, set_MM)
  
  n_MT <- length(set_MT)
  n_MM <- length(set_MM)
  n_intersect <- length(intersect_taxa)
  
  n_MT_only <- length(setdiff(set_MT, set_MM))
  n_MM_only <- length(setdiff(set_MM, set_MT))
  n_intersect <- length(intersect(set_MT, set_MM))
  
  # Radii
  r_MT <- sqrt(n_MT / pi)
  r_MM <- sqrt(n_MM / pi)
  
  # Distance between circle centers for proportional overlap
  d <- r_MT + r_MM - sqrt(n_intersect / pi)
  
  df_circles <- data.frame(
    x = c(0, d),
    y = c(0, 0),
    r = c(r_MT, r_MM),
    label = c("MT", "MM"),
    fill = c("#4981BF", "#F4FAFE")
  )
  
  # Calculate true intersection center
  x_overlap <- (r_MT^2 - r_MM^2 + d^2) / (2*d)
  
  # Plot
  p <- ggplot() +
    geom_circle(data = df_circles, aes(x0 = x, y0 = y, r = r, fill = fill),
                alpha = 0.5, color = "black") +
    # Labels above circles
    geom_text(data = df_circles,
              aes(x = x, y = y + r + 0.1 * max(r), label = label),
              fontface = "bold") +
    # Counts inside circles
    geom_text(aes(x = 0, y = 0, label = n_MT_only), fontface = "bold") +
    geom_text(aes(x = d, y = 0, label = n_MM_only), fontface = "bold") +
    geom_text(aes(x = x_overlap, y = 0, label = n_intersect), fontface = "bold") +
    coord_fixed() +
    theme_void() +
    scale_fill_identity()
  
  assign(plot_name, p, envir = .GlobalEnv)
  return(p)
}

class2 <- unique_taxon_ggplot(df_result_class_MT, df_result_class_MM, "class2_p")
family2 <- unique_taxon_ggplot(df_result_family_MT, df_result_family_MM, "family2_p")
genus2 <- unique_taxon_ggplot(df_result_genus_MT, df_result_genus_MM, "genus2_p")

############################

#4)  Spearman rank correlation

############################

# Spearman correlation function
spearman_rho_plot <- function(MT_df, MM_df, metadata, taxonomic_level = "Level", plot_name = NULL) {
  
  # Ensure taxa match
  common_taxa <- intersect(MT_df$Taxon, MM_df$Taxon)
  MT_df <- MT_df %>% filter(Taxon %in% common_taxa) %>% arrange(Taxon)
  MM_df <- MM_df %>% filter(Taxon %in% common_taxa) %>% arrange(Taxon)
  
  # Remove Taxon column for calculation
  MT_matrix <- MT_df %>% select(-Taxon)
  MM_matrix <- MM_df %>% select(-Taxon)
  
  # Compute Spearman correlation per pond
  pond_names <- colnames(MT_matrix)
  spearman_results <- data.frame(
    SampleID = pond_names,
    Spearman_rho = NA,
    p_value = NA
  )
  
  for (i in seq_along(pond_names)) {
    pond <- pond_names[i]
    x <- MT_matrix[[pond]]
    y <- MM_matrix[[pond]]
    
    if(length(unique(x)) > 1 && length(unique(y)) > 1) {
      cor_test <- cor.test(x, y, method = "spearman")
      spearman_results$Spearman_rho[i] <- cor_test$estimate
      spearman_results$p_value[i] <- cor_test$p.value
    } else {
      spearman_results$Spearman_rho[i] <- NA
      spearman_results$p_value[i] <- NA
    }
  }
  
  # Merge metadata
  spearman_results <- spearman_results %>%
    left_join(metadata, by = "SampleID") %>%
    mutate(
      Nutrient = ifelse(phos == 40, "Moderate", "High"),
      Temperature = case_when(
        temp == "High"    ~ "Heated",
        temp == "Ambient" ~ "Ambient",
        TRUE ~ NA_character_
      ),
      treatment = factor(treatment, levels = c("tn", "Tn", "tN", "TN"))
    ) %>%
    mutate(
      Nutrient = factor(Nutrient, levels = c("Moderate", "High"))
    )
  
  # Define colors and shapes
  nutrient_colors <- c("Moderate" = "#DDCC77", "High" = "#117733")
  temp_shapes <- c("Heated" = 17, "Ambient" = 16)
  
  # Plot Spearman rho per treatment
  p <- ggplot(spearman_results, aes(x = treatment, y = Spearman_rho,
                                    color = Nutrient, shape = Temperature)) +
    geom_point(size = 4, position = position_jitter(width = 0.05, height = 0)) +
    geom_hline(yintercept = c(1, 0.5, 0, -0.5),
               color = "grey30", linewidth = 0.1) +
    scale_color_manual(values = nutrient_colors) +
    scale_shape_manual(values = temp_shapes) +
    ylim(-1, 1) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = "Treatment", y = "Spearman rho") +
    theme_classic(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      
      # Axis lines: all same thickness and grey
      axis.line = element_line(color = "grey30", linewidth = 0.8),
      axis.line.x.top = element_line(color = "grey30", linewidth = 0.8),
      axis.line.y.right = element_line(color = "grey30", linewidth = 0.8),
      
      # Axis ticks
      axis.ticks = element_line(color = "grey30"),
      axis.ticks.length = unit(0.15, "cm"),
      
      # Axis text and title (black, larger)
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(color = "black", size = 17),
      
      plot.title = element_blank(),
      
      # Legend
      legend.text = element_text(size = 9.2),
      legend.title = element_text(size = 10.7),
      legend.spacing.y = unit(-0.2, "cm")
    ) +
    guides(
      color = guide_legend(title = "Nutrient Enrichment"),
      shape = guide_legend(title = "Temperature")
    )
  
  if(!is.null(plot_name)) assign(plot_name, p, envir = .GlobalEnv)
  
  return(list(correlation_table = spearman_results, plot = p))
}

#Using function to create Spearman scatterplots for class, family and genus
result_class <- spearman_rho_plot(df_result_class_MT, df_result_class_MM,
                                  metadata = metadata_mt,
                                  taxonomic_level = "Class",
                                  plot_name = "class2p")
result_family <- spearman_rho_plot(df_result_family_MT, df_result_family_MM,
                                  metadata = metadata_mt,
                                  taxonomic_level = "Family",
                                  plot_name = "family2p")
result_genus <- spearman_rho_plot(df_result_genus_MT, df_result_genus_MM,
                                  metadata = metadata_mt,
                                  taxonomic_level = "Genus",
                                  plot_name = "genus2p")


############################

#5)  Combine Venn Diagram and Spearman

############################

#Create column titles
title_class  <- ggdraw() + 
  draw_label("Class", fontface = "bold", size = 20, x = 0.5, y = 0.5, hjust = 0.5)

title_family <- ggdraw() + 
  draw_label("Family", fontface = "bold", size = 20, x = 0.5, y = 0.5, hjust = 0.5)

title_genus  <- ggdraw() + 
  draw_label("Genus", fontface = "bold", size = 20, x = 0.5, y = 0.5, hjust = 0.5)


titles <- plot_grid(title_class, title_family, title_genus, ncol = 3)

#Extract legend from one of the Spearman plots
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}
legend <- get_legend(class2p)

# Remove legends from each Spearman plot
class2p_noleg  <- class2p  + theme(legend.position = "none",
                                   axis.title.x = element_text(color = "white"))

family2p_noleg <- family2p + theme(legend.position = "none",
                                   axis.title.y = element_blank(),
                                   axis.text.y  = element_blank(),
                                   axis.ticks.y = element_blank())

genus2p_noleg  <- genus2p  + theme(legend.position = "none",
                                   axis.title.y = element_blank(),
                                   axis.text.y  = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.x = element_text(color = "white"))

#Create top row of Venn Diagrams and bottom row of Spearman plots
rowA <- plot_grid(class2_p, family2_p, genus2_p, ncol = 3, labels = NULL)

#Add legend box to the Class Spearman plot 
legend_box <- ggdraw() + draw_plot(legend) + theme_void() + 
  theme( plot.background = element_rect(color = "black", fill = NA, linewidth = 1), 
         panel.background = element_rect(fill = NA, color = NA)) 

# Overlay legend on RIGHT side of genus2p plot 
genus2p_with_legend <- ggdraw() + draw_plot(genus2p_noleg) +
  draw_plot(legend_box, x = 0.50, y = 0.697, width = 0.40, height = 0.5)

#Combine bottom row
rowB <- plot_grid(class2p_noleg, family2p_noleg, genus2p_with_legend, ncol = 3)

#Combine everything
final_plot <- plot_grid(
  # A) Venn diagrams
  plot_grid(rowA, ncol = 1, labels = c("A)"), label_size = 18,
            label_fontface = "bold", label_x = 0, label_y = 1),
  # B) Spearman
  plot_grid(rowB, ncol = 1, labels = c("B)"), label_size = 18,
            label_fontface = "bold", label_x = 0, label_y = 1.10),
  ncol = 1,
  rel_heights = c(1, 1)
)

#Add titles above the graph
final_plot_with_titles <- plot_grid(
  titles,       
  final_plot,
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1)
)

#Display
final_plot_with_titles

############################

#6)  Save plot

############################

# Save the figure
ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Venn_diagaram_w_Spearman_MT_MM.png",
  plot = final_plot_with_titles,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

