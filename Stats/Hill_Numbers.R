######## Alpha Diversity 
# This code with calculate Alpha Diversity (Hill Numbers) based on MT, MM, 18S, COI and 23S at the Class, Family and Genus level
# MT = Metatranscriptomics
# MM = Multi-marker metabarcoding

#Load Libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(cowplot)

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

############################

#3)     Hill Numbers

############################

# Put all your cleaned tables into a named list
df_list <- list(
  class_MT = df_result_class_MT,
  family_MT = df_result_family_MT,
  genus_MT = df_result_genus_MT,
  
  class_MM = df_result_class_MM,
  family_MM = df_result_family_MM,
  genus_MM = df_result_genus_MM,
  
  class_18S = df_result_class_18S,
  family_18S = df_result_family_18S,
  genus_18S = df_result_genus_18S,
  
  class_COI = df_result_class_COI,
  family_COI = df_result_family_COI,
  genus_COI = df_result_genus_COI,
  
  class_23S = df_result_class_23S,
  family_23S = df_result_family_23S,
  genus_23S = df_result_genus_23S
)

# Empty master Hill number summary table
hill_master <- data.frame(
  Dataset = character(),
  Level = character(),
  q0 = numeric(),
  q1 = numeric(),
  q2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop and compute q=0,1,2 for each dataset
hill_master <- data.frame()

for (df_name in names(df_list)) {
  df <- df_list[[df_name]]
  
  df_mat <- df[, colnames(df) != "Taxon"]
  comm <- t(df_mat)
  comm <- as.matrix(comm)
  mode(comm) <- "numeric"
  
  q0 <- specnumber(comm)
  q1 <- exp(diversity(comm, index="shannon"))
  q2 <- 1 / diversity(comm, index="simpson")
  
  parts <- strsplit(df_name, "_")[[1]]
  level <- parts[1]
  dataset <- parts[2]
  
  hill_master <- rbind(hill_master,
                       data.frame(
                         Dataset = dataset,
                         Level = level,
                         Sample = rownames(comm),
                         q0 = q0,
                         q1 = q1,
                         q2 = q2
                       ))
}

write.csv(hill_master, "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/hill_diversity.csv", row.names = FALSE)

############################

#4)     Plot

############################

#Merge Metatdata
hill_master_mt <- hill_master %>%
  filter(Dataset == "MT") %>%
  left_join(metadata_mt %>% select(SampleID, treatment), 
            by = c("Sample" = "SampleID"))

hill_master_mm <- hill_master %>%
  filter(Dataset == "MM") %>%
  left_join(metadata %>% select(SampleID, treatment), 
            by = c("Sample" = "SampleID"))

hill_master_18S <- hill_master %>%
  filter(Dataset == "18S") %>%
  left_join(metadata %>% select(SampleID, treatment), by = c("Sample" = "SampleID"))

hill_master_COI <- hill_master %>%
  filter(Dataset == "COI") %>%
  left_join(metadata %>% select(SampleID, treatment), by = c("Sample" = "SampleID"))

hill_master_23S <- hill_master %>%
  filter(Dataset == "23S") %>%
  left_join(metadata %>% select(SampleID, treatment), by = c("Sample" = "SampleID"))

# Combine all into one master table
hill_master_all <- bind_rows(
  hill_master_mt,
  hill_master_mm,
  hill_master_18S,
  hill_master_COI,
  hill_master_23S
)

# Reorder Dataset factor
hill_master_all$Dataset <- factor(
  hill_master_all$Dataset,
  levels = c("MT", "MM", "18S", "COI", "23S")
)

#Calculate Mean + SD
hill_summary <- hill_master_all %>%
  group_by(Dataset, Level, treatment) %>%
  summarise(
    mean_q1 = mean(q1, na.rm=TRUE),
    sd_q1 = sd(q1, na.rm=TRUE),
    .groups="drop"
  )

#Barplot function
plot_q1_level <- function(df, level_name) {
  df_plot <- df %>% filter(Level == level_name)
  
  ggplot(df_plot, aes(x=Dataset, y=mean_q1, fill=treatment)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
    geom_errorbar(aes(ymin=mean_q1 - sd_q1, ymax=mean_q1 + sd_q1),
                  position=position_dodge(width=0.8), width=0.2) +
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
    labs(y="Hill q1 (Shannon Diversity) Index", x=NULL) +
    theme_classic() +  # removes grid lines
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle=0, hjust=0.5),
      axis.line = element_line(),
      axis.ticks = element_line()
    )
}

#Generate plots
p_class <- plot_q1_level(hill_summary, "class")
p_family <- plot_q1_level(hill_summary, "family")
p_genus <- plot_q1_level(hill_summary, "genus")

#Extract Legend
legend <- get_legend(p_class)

# Remove legends from individual plots
p_class <- p_class + theme(legend.position="none")
p_family <- p_family + theme(legend.position="none")
p_genus <- p_genus + theme(legend.position="none")

# Stack plots without legend
plots_stacked <- plot_grid(
  p_class, p_family, p_genus,
  ncol = 1,
  align = "v",
  axis = "lr"
)

# Add A/B/C labels with space on the left
final_plot <- ggdraw() +
  draw_plot(plots_stacked, x = 0.05, width = 0.95) +   # shift plots right a bit
  draw_plot_label(
    label = c("A)","B)","C)"),
    x = 0.05,                     # now this is safely visible
    y = c(1, 0.68, 0.34),      # top/middle/bottom
    hjust = 1,
    size = 14
  )

# Add shared legend on the right
final_plot_with_legend <- plot_grid(final_plot, legend, ncol=2, rel_widths=c(1,0.15))

final_plot_with_legend

# Save
ggsave("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/q1_barplots.png", final_plot_with_legend,
       width=8, height=10, dpi=300)

############################

#5)     Wilcoxon 

############################

# Define all combinations
methods <- c("MM", "18S", "COI", "23S")
datasets <- c("MT", "MM", "18S", "COI", "23S")
levels <- c("class", "family", "genus")
hill_qs <- c("q0", "q1", "q2")

# Make sure 'phos' column exists in hill_master_all
hill_master_all <- hill_master_all %>%
  left_join(metadata %>% select(SampleID, phos),
            by = c("Sample" = "SampleID")) %>%
  mutate(phos = ifelse(Dataset == "MT",
                       metadata_mt$phos[match(Sample, metadata_mt$SampleID)],
                       phos))

# Initialize empty results table
results_all <- data.frame(
  Level = character(),
  Hill = character(),
  TestType = character(),
  Comparison = character(),
  Method = character(),
  p_value = numeric(),
  adj_p = numeric(),
  stringsAsFactors = FALSE
)

# Loop over taxonomic levels and Hill numbers
for(lv in levels){
  for(hq in hill_qs){
    
    # --- MT vs other methods ---
    for(m in methods){
      mt_data <- hill_master_all %>% 
        filter(Dataset == "MT", Level == lv) %>% 
        pull(all_of(hq))
      
      other_data <- hill_master_all %>% 
        filter(Dataset == m, Level == lv) %>% 
        pull(all_of(hq))
      
      test <- wilcox.test(mt_data, other_data)
      
      results_all <- rbind(
        results_all,
        data.frame(
          Level = lv,
          Hill = hq,
          TestType = "MT_vs_other",
          Comparison = paste("MT vs", m),
          Method = m,
          p_value = test$p.value,
          adj_p = NA
        )
      )
    }
    
    # --- Nutrient effect within each method ---
    for(ds in datasets){
      df <- hill_master_all %>% filter(Dataset == ds, Level == lv)
      
      test <- wilcox.test(df[[hq]] ~ df$phos)
      
      results_all <- rbind(
        results_all,
        data.frame(
          Level = lv,
          Hill = hq,
          TestType = "Eutrophication",
          Comparison = "High vs Moderate",
          Method = ds,
          p_value = test$p.value,
          adj_p = NA
        )
      )
    }
  }
}

# Adjust p-values for multiple testing separately by TestType
results_all <- results_all %>%
  group_by(TestType) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# View results
head(results_all)

# Save to CSV
write.csv(results_all, "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/hill_alpha_stats.csv",
          row.names = FALSE)
