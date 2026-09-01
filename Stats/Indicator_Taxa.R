######## Co-Occurence 
# This code will run Simper to determine the taxa driving the variation
# This code will also create Co-Occurence networks
# MT = Metatranscriptomics
# MM = Multi-marker metabarcoding

#Load libraries
library(vegan)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

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

############################

#3)     Simper

############################

#Simper Anlaysis function
simper_analysis <- function(df_taxa, metadata, level_name, marker_name, top_n = 10) {
  
  # 1) Perform simper analysis
  simper_x <- simper(df_taxa, metadata$phos, permutations = 999, parallel = 1)
  
  # 2) Obtain summary
  simper_list <- summary(simper_x, ordered = TRUE,
                         digits = max(3, getOption("digits") - 3))
  
  # 3) Extract top N taxa and add taxa + Contri%
  do.call(rbind, lapply(names(simper_list), function(contrast_name) {
    df <- simper_list[[contrast_name]]
    
    # Convert row names to a column called 'taxa'
    df$taxa <- rownames(df)
    rownames(df) <- NULL
    
    # Compute contribution percentage
    df$Contri_percent <- (df$average / sum(df$average)) * 100
    
    # Reorder and subset top N
    df <- df[order(-df$average), ]
    df <- head(df, top_n)
    
    df$contrast <- contrast_name
    df$taxa_level <- level_name
    df$marker <- marker_name
    df
  }))
}

s_c_MT <- simper_analysis(df_result_class_MT, metadata_mt, "class", "MT", 10)
s_f_MT <- simper_analysis(df_result_family_MT, metadata_mt, "family", "MT", 10)
s_g_MT <- simper_analysis(df_result_genus_MT, metadata_mt, "genus", "MT", 10)

s_c_MM <- simper_analysis(df_result_class_MM, metadata, "class", "MM", 10)
s_f_MM <- simper_analysis(df_result_family_MM, metadata, "family", "MM", 10)
s_g_MM <- simper_analysis(df_result_genus_MM, metadata, "genus", "MM", 10)

s_c_18S <- simper_analysis(df_result_class_18S, metadata, "class", "18S", 10)
s_f_18S <- simper_analysis(df_result_family_18S, metadata, "family", "18S", 10)
s_g_18S <- simper_analysis(df_result_genus_18S, metadata, "genus", "18S", 10)

s_c_COI <- simper_analysis(df_result_class_COI, metadata, "class", "COI", 10)
s_f_COI <- simper_analysis(df_result_family_COI, metadata, "family", "COI", 10)
s_g_COI <- simper_analysis(df_result_genus_COI, metadata, "genus", "COI", 10)

s_c_23S <- simper_analysis(df_result_class_23S, metadata, "class", "23S", 10)
s_f_23S <- simper_analysis(df_result_family_23S, metadata, "family", "23S", 10)
s_g_23S <- simper_analysis(df_result_genus_23S, metadata, "genus", "23S", 10)


# Combine all levels into one data frame
total_simper_top10 <- rbind(s_c_MT, s_f_MT, s_g_MT, s_c_MM, s_f_MM, s_g_MM, s_c_18S, s_f_18S, s_g_18S, s_c_COI, s_f_COI, s_g_COI, s_c_23S, s_f_23S, s_g_23S)

# View first few rows
head(total_simper_top10)

# Save as CSV
write.csv(total_simper_top10, "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/top10_simper.csv", row.names = FALSE)
