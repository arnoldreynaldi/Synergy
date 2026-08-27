library(dplyr)
library(ggplot2)

# Create an empty list to store median dataframes
median_list <- list()

for (i in seq_len(nrow(plot_combinations))) {
  
  meas      <- plot_combinations$Measurement_type[i]
  cell_type <- plot_combinations$Cell_Type[i]
  
  sub_data <- plot_data %>%
    filter(Measurement_type == meas, Cell_Type == cell_type)
  
  if (nrow(sub_data) == 0 || sum(!is.na(sub_data$Value)) == 0) next
  
  unit_label <- unique(sub_data$Unit)[1]
  
  # For count data, keep only positive values so log scale works
  if (unit_label == "number") {
    sub_data <- sub_data %>% filter(Value > 0)
  }
  
  # Compute median per Vaccine and Day
  med_summary <- sub_data %>%
    group_by(Vaccine, Day) %>%
    summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")
  
  # Save this summary with metadata for later
  med_summary <- med_summary %>%
    mutate(Measurement_type = meas, Cell_Type = cell_type, Unit = unit_label)
  median_list[[length(median_list) + 1]] <- med_summary
  
  # Plot only the median lines
  p <- ggplot(med_summary, aes(x = Day, y = Median, color = Vaccine)) +
    geom_line(size = 1) +
    geom_point(size = 2) +   # optional, remove if you want only lines
    labs(
      title = paste0(meas, "  [", unit_label, "]"),
      x = "Day",
      y = "Median Value",
      color = "Vaccine"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
  
  # Log scale for counts
  if (unit_label == "number") {
    p <- p + scale_y_log10()
  }
  
  # Save PNG
  safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
  file_name <- paste0(cell_type, "_", safe_meas, "_median.png")
  file_path <- file.path(figure_folder, file_name)
  
  ggsave(
    filename = file_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Combine all median dataframes into one
all_medians <- bind_rows(median_list)

# Save as CSV
write.csv(all_medians, "median_values.csv", row.names = FALSE)

# View the first few rows
head(all_medians)



library(dplyr)

# Assuming all_medians is already created from the previous loop
normalized_medians <- all_medians %>%
  group_by(Measurement_type, Cell_Type, Vaccine) %>%
  mutate(
    baseline_Day = min(Day),
    baseline_Median = Median[Day == baseline_Day][1],
    Normalized = Median / baseline_Median
  ) %>%
  ungroup() %>%
  filter(!is.na(Normalized), is.finite(Normalized))



library(ggplot2)

# Create output folder for normalized plots
norm_folder <- "normalized_median_plots"
dir.create(norm_folder, showWarnings = FALSE)

# Loop over combinations again using all_medians or normalized_medians
for (i in seq_len(nrow(plot_combinations))) {
  
  meas      <- plot_combinations$Measurement_type[i]
  cell_type <- plot_combinations$Cell_Type[i]
  
  sub_norm <- normalized_medians %>%
    filter(Measurement_type == meas, Cell_Type == cell_type)
  
  if (nrow(sub_norm) == 0) next
  
  unit_label <- unique(sub_norm$Unit)[1]
  
  p <- ggplot(sub_norm, aes(x = Day, y = Normalized, color = Vaccine)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0(meas, "  [", unit_label, "]  (Normalized)"),
      x = "Day",
      y = "Normalized median (fold change)",
      color = "Vaccine"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    ) +
    scale_y_log10()   # remove if you prefer linear
  
  safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
  file_name <- paste0(cell_type, "_", safe_meas, "_normalized_median.png")
  ggsave(file.path(norm_folder, file_name), p, width = 8, height = 6, dpi = 300)
}


write.csv(normalized_medians, "normalized_medians.csv", row.names = FALSE)