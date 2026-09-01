library(tidyverse)

# 1. Load table
timeline <- read_tsv("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/Taxa_Lists/Timeline.tsv")

# 2. Convert to long format
timeline_long <- timeline %>%
  pivot_longer(cols = -Day, names_to = "Variable", values_to = "Value") %>%
  filter(tolower(Value) == "yes")

# 3. Reorder variables so eRNA is last
timeline_long$Variable <- factor(
  timeline_long$Variable,
  levels = c(setdiff(unique(timeline_long$Variable), "eRNA"), "eRNA")
)

# 4. Plot
ggplot(timeline_long, aes(x = Day, y = Variable)) +
  # black dots
  geom_point(size = 3, color = "black") +
  
  # horizontal reference lines for each variable
  geom_hline(
    yintercept = seq_along(unique(timeline_long$Variable)),
    color = "black", linewidth = 0.2, alpha = 0.5
  ) +
  
  # vertical event lines (stop before labels)
  geom_segment(aes(x = 14, xend = 14, y = 0.5, yend = length(unique(timeline_long$Variable)) + 0.1),
               color = "black", linewidth = 0.8) +
  geom_segment(aes(x = 35, xend = 35, y = 0.5, yend = length(unique(timeline_long$Variable)) + 0.1),
               color = "black", linewidth = 0.8) +
  
  # horizontal text labels (slightly above plot, visible now)
  annotate("text", x = 14, y = length(unique(timeline_long$Variable)) + 0.3,
           label = "Nutrients Added", angle = 0, vjust = 0, hjust = 0.5,
           size = 4.2, fontface = "bold") +
  annotate("text", x = 35, y = length(unique(timeline_long$Variable)) + 0.3,
           label = "Heaters Added", angle = 0, vjust = 0, hjust = 0.5,
           size = 4.2, fontface = "bold") +
  
  # extend x-axis to include day 0
  scale_x_continuous(
    breaks = seq(0, max(timeline$Day), by = 5),
    limits = c(0, max(timeline$Day))
  ) +
  
  # compact vertical spacing
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.15))) +
  
  # axis labels + base theme
  labs(x = "Day of Experiment", y = NULL) +
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(t = 80, r = 20, b = 10, l = 10),
    clip = "off"
  )

# 5. Save
ggsave(
  "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_1/figures/Final_Plot_MM_Marker_Partition/timeline_plot.png",
  width = 6,
  height = 5,
  dpi = 300
)



