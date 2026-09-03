process_summary <- function(stat = c("median", "mean")) {
  stat <- match.arg(stat)
  stat_fun <- if (stat == "median") median else mean
  stat_label <- tools::toTitleCase(stat)
  
  # Folders
  plot_dir <- file.path(figure_folder, paste0("Normalised_", stat))
  table_dir <- file.path(table_folder, paste0("Normalised_", stat))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Vaccine order and colours
  vaccine_levels <- c("TM", "TMd21", "SOL", "IC")
  formulation_cols <- c(
    TM    = "#2F75B5",
    TMd21 = "#D99A2B",
    SOL   = "#4A9E7D",
    IC    = "#B86691"
  )
  
  # Compute summary for each combination
  summary_list <- list()
  for (i in seq_len(nrow(plot_combinations))) {
    meas <- plot_combinations$Measurement_type[i]
    cell_type <- plot_combinations$Cell_Type[i]
    
    sub_data <- plot_data %>%
      filter(Measurement_type == meas, Cell_Type == cell_type)
    
    if (nrow(sub_data) == 0 || sum(!is.na(sub_data$Value)) == 0) next
    
    unit_label <- unique(sub_data$Unit)[1]
    if (unit_label == "number") sub_data <- filter(sub_data, Value > 0)
    
    sum_df <- sub_data %>%
      group_by(Vaccine, Day) %>%
      summarise(Value = stat_fun(Value, na.rm = TRUE), .groups = "drop") %>%
      mutate(Measurement_type = meas, Cell_Type = cell_type, Unit = unit_label)
    
    summary_list[[length(summary_list) + 1]] <- sum_df
  }
  
  all_summaries <- bind_rows(summary_list)
  
  # Normalize by baseline (minimum Day)
  normalised <- all_summaries %>%
    group_by(Measurement_type, Cell_Type, Vaccine) %>%
    mutate(
      baseline_Day = min(Day),
      baseline_Value = Value[Day == baseline_Day][1],
      normalised = Value / baseline_Value
    ) %>%
    ungroup() %>%
    filter(!is.na(normalised), is.finite(normalised))
  
  # Plot and save
  for (i in seq_len(nrow(plot_combinations))) {
    meas <- plot_combinations$Measurement_type[i]
    cell_type <- plot_combinations$Cell_Type[i]
    
    sub_norm <- normalised %>%
      filter(Measurement_type == meas, Cell_Type == cell_type) %>%
      mutate(Vaccine = factor(Vaccine, levels = vaccine_levels))
    
    if (nrow(sub_norm) == 0) next
    
    unit_label <- unique(sub_norm$Unit)[1]
    
    # Map unit to display label for title
    unit_label_display <- switch(
      unit_label,
      "number" = "Number of Cells",
      "percentage" = "%",
      unit_label
    )
    
    p <- ggplot(sub_norm, aes(x = Day, y = normalised, color = Vaccine)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.5) +
      labs(
        title = paste0(meas, " [", unit_label_display, "] (normalised ", stat_label, ")"),
        x = "Day",
        y = paste("normalised", stat_label, "(fold change)"),
        color = "Vaccine"
      ) +
      scale_color_manual(
        values = formulation_cols,
        breaks = vaccine_levels
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right"
      )
    
    if (unit_label == "number") p <- p + scale_y_log10()
    
    safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
    file_name <- paste0(cell_type, "_", safe_meas, "_normalised_", stat, ".png")
    ggsave(file.path(plot_dir, file_name), p, width = 8, height = 6, dpi = 300)
  }
  
  # Save CSV
  write.csv(normalised, file.path(table_dir, paste0("normalised_", stat, "s.csv")), row.names = FALSE)
}

# Run for both median and mean
process_summary("median")
process_summary("mean")