##################### Comparisons with Fluoroprobe

#This Code will compare FLuoroprobe counts with MT, MM, 18S and 23S

#Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)
library(cowplot)

############################

#1)     Load Files

############################

#Load Fluoroprobe data
fluoro <- read.table("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/F_August_10_LEAP2021.tsv", sep="\t", header=TRUE, comment.char="")

######## MT
#Load MT Algal identified through Class
filepath_class_mt <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/MT_class_raw_algae.tsv"
df_mt <- read_tsv(filepath_class_mt)

######## 18S
#Load MT Algal identified through Class
filepath_class_18S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/18S_PR2_class_algal.tsv"
df_18S <- read_tsv(filepath_class_18S)

######## 23S
#Load MT Algal identified through Class
filepath_class_23S <- "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Fluoro_Morph/23S_Phytool_class_algal.tsv"
df_23S <- read_tsv(filepath_class_23S)

############################

#2)     Proccess Files

############################

#Function to process files
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
  
  return(df_rel)
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

############################

#3)     Linear Regression

############################

# List your molecular datasets for convenience
molecular_datasets <- list(
  MT = df_result_MT,
  `18S` = df_result_18S,
  `23S` = df_result_23S
)

# Algae types to compare
algae_pairs <- list(
  list(fluoro_col = "Green_Algae", molecular_col = "Green_Algae"),
  list(fluoro_col = "Brown_Pigmented_Algae", molecular_col = "Brown_Pigmented_Algae")
)

# Function to extract and merge relevant data
prepare_data <- function(fluoro_df, mol_df, algae_type) {
  mol_df_long <- mol_df %>%
    filter(Algae == algae_type$molecular_col) %>%
    pivot_longer(-Algae, names_to = "SampleID", values_to = "Molecular_value")
  
  merged_df <- fluoro_df %>%
    select(SampleID, !!algae_type$fluoro_col) %>%
    rename(Fluoro_value = !!algae_type$fluoro_col) %>%
    inner_join(mol_df_long, by = "SampleID")
  
  return(merged_df)
}

# Run regressions for all datasets and algae types
results <- expand.grid(
  Method = names(molecular_datasets),
  Algae = c("Green_Algae", "Brown_Pigmented_Algae"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Data = pmap(list(Method, Algae), function(Method, Algae) {
      mol_df <- molecular_datasets[[Method]]
      algae_type <- list(fluoro_col = ifelse(Algae == "Brown_Pigmented_Algae", "Brown_Pigmented_Algae", "Green_Algae"),
                         molecular_col = Algae)
      prepare_data(fluoro_avg, mol_df, algae_type)
    }),
    Model = map(Data, ~ lm(Fluoro_value ~ Molecular_value, data = .x)),
    Summary = map(Model, broom::tidy),
    R2 = map_dbl(Model, ~ summary(.x)$r.squared)
  )

# View regression summaries
results %>%
  select(Method, Algae, R2, Summary) %>%
  print()

# Create a flat summary table of regression results
results_table <- results %>%
  select(Method, Algae, R2, Summary) %>%
  unnest(Summary) %>%
  select(Method, Algae, term, estimate, std.error, statistic, p.value, R2)

# Print as a nice table
print(results_table)

write.csv(results_table, "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Fluoro_linear_regression_results.csv", row.names = FALSE)


############################

#4)     Plot

############################

plots <- results %>%
  mutate(
    Plot = map2(Data, paste(Method, Algae, sep = "_"), function(df, title) {
      ggplot(df, aes(y = Molecular_value, x = Fluoro_value)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", se = FALSE, color = "black") +
        coord_fixed(ratio = 1) +   # ← makes each panel square
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_blank(), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10)
        ) +
        labs(
          x = "Biomass (µg/L)",
          y = "Normalized read count"
        )
    })
  )

# Plot each regression
plots <- results %>%
  mutate(
    Plot = pmap(list(Data, Method, Algae, R2), function(df, Method, Algae, r2) {
      
      r2_label <- paste0("R² = ", round(r2, 2))
      
      base_plot <- ggplot(df, aes(y = Molecular_value, x = Fluoro_value)) +
        geom_jitter(width = 0.1, height = 0.1, size = 2) +
        geom_smooth(method = "lm", se = FALSE, color = "black") +
        annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2,
                 label = r2_label, size = 4) +
        theme_classic(base_size = 12) +
        labs(
          x = NULL,     # remove x title
          y = NULL
        ) +
        theme(
          axis.title.x = element_blank()   # ⬅ ONLY remove X title, keep numbers
        )
      
      # Special rule: ONLY MT keeps y-axis text + label
      if (Method == "MT") {
        base_plot <- base_plot +
          labs(y = "Normalized read count") +
          theme(
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 20),
            axis.ticks.y = element_line()
          )
      } else {
        base_plot <- base_plot +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
      }
      
      return(base_plot)
    })
  )

#Create a grid of only GREEN algae and only MT / 18S / 23S
green_plots <- plots %>%
  filter(
    Algae == "Green_Algae",
    Method != "MM"
  ) %>%
  arrange(factor(Method, levels = c("MT", "18S", "23S"))) %>%
  pull(Plot)

# Add method labels above each column
labels_top <- ggdraw() +
  draw_label("MT",  x = 0.2, y = 0.5, size = 20, fontface = "bold") +
  draw_label("18S", x = 0.53,  y = 0.5, size = 20, fontface = "bold") +
  draw_label("23S", x = 0.87, y = 0.5, size = 20, fontface = "bold")

# Single row of plots
green_row <- plot_grid(
  plotlist = green_plots,
  ncol = 3,
  align = "hv"
)

# Combine labels + row
final_plot <- plot_grid(
  labels_top,
  green_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Add bottom X label
x_label <- ggdraw() + 
  draw_label(
    expression("Fluorometry estimates (" * mu * "gL"^{-1} * ")"),
    size = 20,
    fontface = "bold"
  )

final_plot_with_xlabel <- plot_grid(
  final_plot,
  x_label,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

final_plot_with_xlabel

############################

#5)     Save Plot

############################

ggsave("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/Fluoroprobe_Comparison.png", final_plot_with_xlabel, width = 14, height = 8, dpi = 300)

