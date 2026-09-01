######## Mantel Tests
# This code will run Mantel tests to compare MT and MM to environmental variables
# Will also look at 18S, COI and 23S but the graphs are only MT and MM
# MT = Metatranscriptomics
# MM = Multi-marker metabarcoding
#Removes G4 from calculations

#Load libraries
library(vegan)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(linkET)

############################

#1)     Load Files

############################

#Load MT metadata
metadata_mt <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Metadata/LEAP2021_metadata.tsv", sep="\t", header=TRUE, comment.char="")

#Load environmental data
env <- readr::read_tsv(
  "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/physiochemical/Environmental_variables.tsv",
  locale = readr::locale(encoding = "UTF-8")
)

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
  col_sums[col_sums == 0] <- NA
  
  df_rel <- sweep(df_clean, 2, col_sums, FUN = "/") * 100
  df_rel[is.na(df_rel)] <- 0  
  
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

#3)     Mantel Test

############################

# Function to run Mantel tests for a single dataset against all environmental variables
run_mantel_analysis <- function(df_taxon, env_data, method_name, tax_level, distance_method) {
  
  # Ensure sample IDs match between datasets
  # Get sample IDs from the taxon data (rownames)
  taxon_samples <- rownames(df_taxon)
  
  # Match samples between taxon data and environmental data
  common_samples <- intersect(taxon_samples, env_data$SampleID)
  
  if(length(common_samples) < 3) {
    warning(paste("Not enough common samples for", method_name, tax_level, distance_method))
    return(NULL)
  }
  
  # Subset and order both datasets to match
  df_taxon_matched <- df_taxon[common_samples, ]
  env_matched <- env_data[match(common_samples, env_data$SampleID), ]
  
  # Remove Taxon column if it exists in the taxon data
  if("Taxon" %in% colnames(df_taxon_matched)) {
    df_taxon_matched <- df_taxon_matched[, !colnames(df_taxon_matched) %in% "Taxon"]
  }
  
  # Calculate community distance matrix
  if(distance_method == "Bray-Curtis") {
    comm_dist <- vegdist(df_taxon_matched, method = "bray")
  } else if(distance_method == "Jaccard") {
    # Convert to presence/absence for Jaccard
    df_pa <- (df_taxon_matched > 0) * 1
    comm_dist <- vegdist(df_pa, method = "jaccard")
  }
  
  # Get environmental variable columns (exclude SampleID)
  env_vars <- setdiff(colnames(env_matched), "SampleID")
  
  # Initialize results list
  results_list <- list()
  
  # Run Mantel test for each environmental variable
  for(var in env_vars) {
    # Create distance matrix for single environmental variable
    env_single <- as.data.frame(env_matched[[var]])
    colnames(env_single) <- var
    rownames(env_single) <- common_samples
    
    # Skip if variable has missing values
    if(any(is.na(env_single[[var]]))) {
      warning(paste("Skipping", var, "due to missing values"))
      next
    }
    
    # Calculate Euclidean distance for environmental variable
    env_dist <- dist(env_single, method = "euclidean")
    
    # Run Mantel test
    mantel_result <- mantel(comm_dist, env_dist, method = "pearson", permutations = 999)
    
    # Store results
    results_list[[var]] <- data.frame(
      Method_Name = method_name,
      Taxonomic_Level = tax_level,
      Distance_Method = distance_method,
      Environmental_Variable = var,
      Mantel_R = mantel_result$statistic,
      Mantel_P = mantel_result$signif,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results for this dataset
  if(length(results_list) > 0) {
    return(bind_rows(results_list))
  } else {
    return(NULL)
  }
}

# Read environmental data
env <- readr::read_tsv(
  "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/physiochemical/Environmental_variables.tsv",
  locale = readr::locale(encoding = "UTF-8")
)

# Create list of all datasets with their metadata
datasets <- list(
  list(data = df_result_class_MT, method = "Metatranscriptomics", level = "Class"),
  list(data = df_result_family_MT, method = "Metatranscriptomics", level = "Family"),
  list(data = df_result_genus_MT, method = "Metatranscriptomics", level = "Genus"),
  list(data = df_result_class_MM, method = "Multi-marker", level = "Class"),
  list(data = df_result_family_MM, method = "Multi-marker", level = "Family"),
  list(data = df_result_genus_MM, method = "Multi-marker", level = "Genus"),
  list(data = df_result_class_18S, method = "18S", level = "Class"),
  list(data = df_result_family_18S, method = "18S", level = "Family"),
  list(data = df_result_genus_18S, method = "18S", level = "Genus"),
  list(data = df_result_class_COI, method = "COI", level = "Class"),
  list(data = df_result_family_COI, method = "COI", level = "Family"),
  list(data = df_result_genus_COI, method = "COI", level = "Genus"),
  list(data = df_result_class_23S, method = "23S", level = "Class"),
  list(data = df_result_family_23S, method = "23S", level = "Family"),
  list(data = df_result_genus_23S, method = "23S", level = "Genus")
)

# Run analysis for all datasets and both distance methods
all_results <- list()

for(i in seq_along(datasets)) {
  ds <- datasets[[i]]
  
  cat("Processing:", ds$method, "-", ds$level, "\n")
  
  # Bray-Curtis
  result_bc <- run_mantel_analysis(
    df_taxon = ds$data,
    env_data = env,
    method_name = ds$method,
    tax_level = ds$level,
    distance_method = "Bray-Curtis"
  )
  
  if(!is.null(result_bc)) {
    all_results[[length(all_results) + 1]] <- result_bc
  }
  
  # Jaccard
  result_j <- run_mantel_analysis(
    df_taxon = ds$data,
    env_data = env,
    method_name = ds$method,
    tax_level = ds$level,
    distance_method = "Jaccard"
  )
  
  if(!is.null(result_j)) {
    all_results[[length(all_results) + 1]] <- result_j
  }
}

# Combine all results
final_results <- bind_rows(all_results)

# Reshape to wide format for easier reading
final_results_wide <- final_results %>%
  pivot_wider(
    names_from = Environmental_Variable,
    values_from = c(Mantel_R, Mantel_P),
    names_glue = "{Environmental_Variable}_{.value}"
  )

# Display results
print(final_results_wide)

# Save results
write.csv(final_results, 
          "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/mantel_results_long_format.csv", 
          row.names = FALSE)

write.csv(final_results_wide, 
          "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/mantel_results_wide_format.csv", 
          row.names = FALSE)

cat("- mantel_results_long_format.csv (one row per test)\n")
cat("- mantel_results_wide_format.csv (environmental variables as columns)\n")

############################

#4)     Mantel Test Plot

############################

# Read environmental data
env <- readr::read_tsv(
  "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/physiochemical/Environmental_variables.tsv",
  locale = readr::locale(encoding = "UTF-8")
)

# Read Mantel results
results <- read.csv("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/mantel_results_long_format.csv")

# Function to prepare mantel data for linkET format
prepare_mantel_data <- function(mantel_results, method_filter, taxonomic_level) {
  
  # Filter for specific method and taxonomic level
  mantel_subset <- mantel_results %>%
    filter(Method_Name == method_filter,
           Taxonomic_Level == taxonomic_level)
  
  # Create spec column including method and distance type
  mantel_subset <- mantel_subset %>%
    mutate(spec = case_when(
      Distance_Method == "Bray-Curtis" & method_filter == "Metatranscriptomics" ~ "MT Abundance",
      Distance_Method == "Jaccard" & method_filter == "Metatranscriptomics" ~ "MT Diversity",
      Distance_Method == "Bray-Curtis" & method_filter == "Multi-marker" ~ "MM Abundance",
      Distance_Method == "Jaccard" & method_filter == "Multi-marker" ~ "MM Diversity"
    ))
  
  # Rename columns to match linkET requirements
  mantel_subset <- mantel_subset %>%
    select(spec, env = Environmental_Variable, 
           r = Mantel_R, p = Mantel_P)
  
  # Create categories for p-values
  mantel_subset <- mantel_subset %>%
    mutate(pd = case_when(
      p <= 0.001 ~ "p ≤ 0.001",
      p <= 0.01 ~ "p ≤ 0.01",
      p <= 0.05 ~ "p ≤ 0.05",
      TRUE ~ "p > 0.05"
    ))
  
  # Create categories for r-values (absolute value for sizing)
  mantel_subset <- mantel_subset %>%
    mutate(rd = case_when(
      abs(r) >= 0.5 ~ "> 0.5",
      abs(r) >= 0.25 ~ "0.25-0.5",
      TRUE ~ "< 0.25"
    ))
  
  # Set factor levels
  mantel_subset$pd <- factor(mantel_subset$pd,
                             levels = c("p ≤ 0.001", "p ≤ 0.01", 
                                        "p ≤ 0.05", "p > 0.05"))
  
  mantel_subset$rd <- factor(mantel_subset$rd,
                             levels = c("> 0.5", "0.25-0.5", "< 0.25"))
  
  return(mantel_subset)
}

# Function to create the plot
create_linkET_plot <- function(env_data, mantel_data, title_text) {
  
  # Prepare environmental data (remove SampleID)
  env_numeric <- env_data %>%
    select(-SampleID) %>%
    select_if(is.numeric)
  
  # Calculate correlation matrix
  env_corr <- linkET::correlate(env_numeric, method = "pearson")
  
  # Define color palette for p-values
  p_colors <- c("p ≤ 0.001" = "#D7191C",
                "p ≤ 0.01" = "#FDAE61",
                "p ≤ 0.05" = "#FEE08B",
                "p > 0.05" = "#CCCCCC")
  
  # Create the plot
  p <- qcorrplot(env_corr, type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(data = mantel_data,
                aes(colour = pd, size = rd),
                curvature = nice_curvature()) +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(11, "RdBu"),
                         limits = c(-1, 1)) +
    scale_size_manual(values = c("> 0.5" = 2,
                                 "0.25-0.5" = 1,
                                 "< 0.25" = 0.5)) +
    scale_colour_manual(values = p_colors) +
    guides(size = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"),
                               order = 2),
           colour = guide_legend(title = "Mantel's p",
                                 override.aes = list(size = 3),
                                 order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3)) +
    labs(title = title_text) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

# Prepare mantel data for MT and MM
mantel_MT <- prepare_mantel_data(results, "Metatranscriptomics", "Family")
mantel_MM <- prepare_mantel_data(results, "Multi-marker", "Family")

mantel_combined <- bind_rows(mantel_MT, mantel_MM)

# Create plots
plot_family <- create_linkET_plot(
  env_data = env,
  mantel_data = mantel_combined,
  title_text = ""
)

# Display plots
print(plot_family)

# Save the figure
ggsave(
  filename = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Mantel_MT_MM_Family.png",
  plot = plot_family,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)
