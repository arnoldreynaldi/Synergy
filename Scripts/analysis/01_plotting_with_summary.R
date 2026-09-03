plot_data <- combined_data %>%
  mutate(
    Day = as.numeric(Day),
    Rep_Num = as.character(Rep_Num),
    Vaccine = as.character(Vaccine)
  ) %>%
  filter(
    !is.na(Day),
    !is.na(Rep_Num),
    Rep_Num != "NA",
    Rep_Num != "NaN",
    !is.na(Vaccine),
    Vaccine != "NA",
    Vaccine != "NaN"
  )

# Vaccine order and colours
vaccine_levels <- c("TM", "TMd21", "SOL", "IC")

formulation_cols <- c(
  TM    = "#2F75B5",
  TMd21 = "#D99A2B",
  SOL   = "#4A9E7D",
  IC    = "#B86691"
)

# Plot combinations
plot_combinations <- plot_data %>%
  distinct(Measurement_type, Cell_Type) %>%
  arrange(
    factor(Cell_Type, levels = c("CD4", "CD8")),
    Measurement_type
  )

# Calculate sensible widths using the smallest interval between days
day_values <- sort(unique(plot_data$Day))

if (length(day_values) > 1) {
  min_day_gap <- min(diff(day_values))
} else {
  min_day_gap <- 14
}

# These values are in units of days
dodge_width <- 0.4 * min_day_gap
box_width <- 0.3 * min_day_gap
jitter_width <- 0.08 * min_day_gap


# Function to create and save raw-data plots
save_raw_plots <- function(stat = c("median", "mean")) {
  
  stat <- match.arg(stat)
  
  stat_fun <- if (stat == "median") {
    median
  } else {
    mean
  }
  
  out_folder <- file.path(
    figure_folder,
    paste0("Raw_data_", stat)
  )
  
  dir.create(
    out_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for (i in seq_len(nrow(plot_combinations))) {
    
    meas <- plot_combinations$Measurement_type[i]
    cell_type <- plot_combinations$Cell_Type[i]
    
    sub_data <- plot_data %>%
      filter(
        Measurement_type == meas,
        Cell_Type == cell_type
      ) %>%
      mutate(
        Vaccine = factor(
          Vaccine,
          levels = vaccine_levels
        ),
        Rep_Num = factor(Rep_Num)
      )
    
    # Skip empty plots
    if (
      nrow(sub_data) == 0 ||
      sum(!is.na(sub_data$Value)) == 0
    ) {
      next
    }
    
    # Identify measurement unit
    unit_values <- unique(na.omit(sub_data$Unit))
    
    if (length(unit_values) == 0) {
      unit_label <- "value"
    } else {
      unit_label <- as.character(unit_values[1])
    }
    
    # Remove zero values before log transformation
    if (unit_label == "number") {
      sub_data <- sub_data %>%
        filter(!is.na(Value), Value > 0)
    }
    
    if (
      nrow(sub_data) == 0 ||
      sum(!is.na(sub_data$Value)) == 0
    ) {
      next
    }
    
    # Y-axis label
    y_label <- switch(
      unit_label,
      "number" = "Number of Cells",
      "percentage" = "%",
      "Value"
    )
    
    # Plot
    p <- ggplot(
      sub_data,
      aes(
        x = Day,
        y = Value,
        color = Vaccine,
        shape = Rep_Num
      )
    ) +
      geom_boxplot(
        aes(
          group = interaction(Day, Vaccine),
          fill = Vaccine
        ),
        width = box_width,
        alpha = 0.25,
        outlier.shape = NA,
        linewidth = 0.45,
        position = position_dodge(
          width = dodge_width
        ),
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      geom_point(
        position = position_jitterdodge(
          jitter.width = jitter_width,
          jitter.height = 0,
          dodge.width = dodge_width
        ),
        size = 2,
        alpha = 0.6,
        na.rm = TRUE
      ) +
      stat_summary(
        aes(group = Vaccine),
        fun = stat_fun,
        geom = "line",
        linewidth = 0.8,
        position = position_dodge(
          width = dodge_width
        ),
        na.rm = TRUE
      ) +
      stat_summary(
        aes(group = Vaccine),
        fun = stat_fun,
        geom = "point",
        size = 2,
        position = position_dodge(
          width = dodge_width
        ),
        na.rm = TRUE
      ) +
      labs(
        title = paste0(meas, " [", unit_label, "]"),
        x = "Day",
        y = y_label,
        color = "Vaccine",
        shape = "Replicate"
      ) +
      scale_x_continuous(
        breaks = sort(unique(sub_data$Day)),
        expand = expansion(
          mult = c(0.03, 0.05)
        )
      ) +
      scale_color_manual(
        values = formulation_cols,
        breaks = vaccine_levels,
        na.translate = FALSE
      ) +
      scale_fill_manual(
        values = formulation_cols,
        breaks = vaccine_levels,
        na.translate = FALSE
      ) +
      scale_shape_discrete(
        na.translate = FALSE
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )
    
    # Log scale for cell-number measurements
    if (unit_label == "number") {
      p <- p +
        scale_y_log10()
    }
    
    # Safe filename
    safe_meas <- gsub(
      "[^A-Za-z0-9]",
      "_",
      meas
    )
    
    file_name <- paste0(
      cell_type,
      "_",
      safe_meas,
      ".png"
    )
    
    ggsave(
      filename = file.path(out_folder, file_name),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}


# Generate median plots
save_raw_plots("median")

# Generate mean plots
save_raw_plots("mean")