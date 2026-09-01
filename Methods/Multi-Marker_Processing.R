########### Multi-Marker processing #########
#Combines 18S, COI and 23S at Class, Family and Genus Level

library(dplyr)
library(readr)
library(stringr)
library(tidyverse)

############################

#1)     Set Directories

############################

input_dir  <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/"
output_dir <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MM_Marker_Partition/"

dir.create(output_dir, showWarnings = FALSE)

############################

#2)     Process Files

############################

process_taxa_file <- function(filepath, marker) {
  
  df <- read_tsv(filepath, show_col_types = FALSE)
  
  # 2. Remove unwanted columns
  df <- df %>%
    select(-c(DB_R1, DB_R2, DB_R3, EB_21, Mean))
  
  # 3. Remove rows where Reason contains any text (non-empty)
  df <- df %>%
    filter(is.na(Reason) | Reason == "")
  
  # 4. Remove metadata columns
  df <- df %>%
    select(-False_Positives, -"Non-Target_Taxa", -Reason)
  
  # Identify sample columns
  sample_cols <- setdiff(names(df), "Taxon")
  
  # 5. Remove rows where all sample abundances sum to zero
  df <- df %>%
    rowwise() %>%
    mutate(total_abundance = sum(c_across(all_of(sample_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    filter(total_abundance > 0) %>%
    select(-total_abundance)
  
  # 6. Convert to relative abundance (sum = 100 per sample)
  df_ra <- df
  
  df_ra[sample_cols] <- apply(df[sample_cols], 2, function(x) {
    if (sum(x) == 0) {
      return(rep(0, length(x)))   # avoid division by zero
    } else {
      return((x / sum(x)) * 100)
    }
  })
  
  # Add marker column
  df_ra$Marker <- marker
  
  return(df_ra)
}

############################

#3)     List groups by taxonomic level

############################

files_class <- c(
  "18S_PR2_class_raw.tsv"      = "18S",
  "COI_Porter_class_raw.tsv"   = "COI",
  "23S_Phytool_class_raw.tsv"  = "23S"
)

files_family <- c(
  "18S_PR2_family_raw.tsv"      = "18S",
  "COI_Porter_family_raw.tsv"   = "COI",
  "23S_Phytool_family_raw.tsv"  = "23S"
)

files_genus <- c(
  "18S_PR2_genus_raw.tsv"      = "18S",
  "COI_Porter_genus_raw.tsv"   = "COI",
  "23S_Phytool_genus_raw.tsv"  = "23S"
)

#Process and combine at each level
process_group <- function(file_list) {
  bind_rows(
    lapply(names(file_list), function(fname) {
      process_taxa_file(
        filepath = file.path(input_dir, fname),
        marker = file_list[[fname]]
      )
    })
  )
}

combined_class  <- process_group(files_class)
combined_family <- process_group(files_family)
combined_genus  <- process_group(files_genus)

############################

#4)     Save outputs

############################

write_tsv(combined_class,
          file.path(output_dir, "MM_class_p_raw.tsv"))

write_tsv(combined_family,
          file.path(output_dir, "MM_family_p_raw.tsv"))

write_tsv(combined_genus,
          file.path(output_dir, "MM_genus_p_raw.tsv"))



############################

#5)     Marker Partition

############################

# Directory containing your combined files
dir <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/final_taxonomic_lists/MM_Marker_Partition/"

# File names
files <- c("MM_class_p_raw.tsv",
           "MM_family_p_raw.tsv",
           "MM_genus_p_raw.tsv")

# Function to normalize abundances so each sample sums to 100
normalize_file <- function(filename) {
  
  path <- file.path(dir, filename)
  
  # Read file
  df <- read_tsv(path, show_col_types = FALSE)
  
  # Identify non-abundance columns
  # (Taxon + Marker columns stay as-is)
  non_abund_cols <- c("Taxon", "Marker")
  
  # All other columns assumed to be abundances
  abund_cols <- setdiff(colnames(df), non_abund_cols)
  
  # Convert abundance columns to numeric safely
  df[abund_cols] <- lapply(df[abund_cols], as.numeric)
  
  # Normalize so each sample sums to 100
  df_norm <- df %>%
    mutate(across(all_of(abund_cols), ~ .x / sum(.x) * 100))
  
  # Save over original file
  write_tsv(df_norm, path)
  
  message(paste("✓ Normalized and saved:", filename))
}


# Run normalization on all 3 files
lapply(files, normalize_file)

