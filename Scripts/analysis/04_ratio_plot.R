# ==================================================================
# Ratio plots 
# ==================================================================
# ------------------------------------------------------------------
# 0. Defaults if not already in the environment
# ------------------------------------------------------------------
if (!exists("vaccine_levels")) {
  vaccine_levels <- c("TM", "TMd21", "SOL", "IC")
}
if (!exists("formulation_cols")) {
  formulation_cols <- c(
    TM    = "#2F75B5",
    TMd21 = "#D99A2B",
    SOL   = "#4A9E7D",
    IC    = "#B86691"
  )
}

# ------------------------------------------------------------------
# 1. Bootstrap 95% CI for the mean
# ------------------------------------------------------------------
boot_mean_ci <- function(x, n_boot = 5000, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(data.frame(y = NA_real_, ymin = NA_real_, ymax = NA_real_))
  if (length(x) == 1)
    return(data.frame(y = x, ymin = NA_real_, ymax = NA_real_))
  boots <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  ci <- quantile(boots, c((1 - conf) / 2, 1 - (1 - conf) / 2),
                 na.rm = TRUE)
  data.frame(y    = mean(x),
             ymin = unname(ci[1]),
             ymax = unname(ci[2]))
}

# ------------------------------------------------------------------
# 2. Generic ratio plotter (linear + log10)
# ------------------------------------------------------------------
plot_ratio <- function(numerator_label,
                       denominator_label,
                       numerator_short,
                       denominator_short,
                       cell_type         = "CD4",
                       unit              = "percentage",
                       file_stub         = NULL,
                       n_boot            = 5000) {
  
  ratio_data <- combined_data %>%
    mutate(
      Day       = as.numeric(Day),
      Vaccine   = trimws(as.character(Vaccine)),
      Exp_Date  = trimws(as.character(Exp_Date)),
      Exp_label = trimws(as.character(Exp_label)),
      Exp_label = if_else(
        is.na(Exp_label) | Exp_label == "",
        paste0("unknown_", row_number()),
        Exp_label
      ),
      Mouse_ID  = paste(Exp_Date, Exp_label, sep = "_"),
      Measurement_type = trimws(Measurement_type),
      Measurement_type = gsub("[–—]", "-", Measurement_type)
    ) %>%
    filter(
      Cell_Type == cell_type,
      Unit      == unit,
      Measurement_type %in% c(numerator_label, denominator_label),
      !is.na(Day),
      !is.na(Vaccine), !Vaccine %in% c("NA", "NaN")
    ) %>%
    mutate(
      Pop = if_else(Measurement_type == numerator_label, "NUM", "DEN")
    ) %>%
    select(Day, Vaccine, Mouse_ID, Pop, Value) %>%
    pivot_wider(
      id_cols     = c(Day, Vaccine, Mouse_ID),
      names_from  = Pop,
      values_from = Value,
      values_fn   = mean
    )
  
  if (!all(c("NUM", "DEN") %in% names(ratio_data))) {
    warning("SKIP ", numerator_label, " / ", denominator_label,
            ": one side of the ratio not found for ", cell_type)
    return(invisible(NULL))
  }
  
  ratio_data <- ratio_data %>%
    filter(!is.na(NUM), !is.na(DEN), DEN > 0) %>%
    mutate(
      Ratio   = NUM / DEN,
      Vaccine = factor(Vaccine, levels = vaccine_levels)
    )
  
  if (nrow(ratio_data) == 0) {
    warning("SKIP ", numerator_label, " / ", denominator_label,
            ": no valid pairs after filtering")
    return(invisible(NULL))
  }
  
  message("OK  ", cell_type, "  ",
          numerator_short, "/", denominator_short,
          "  (", nrow(ratio_data), " rows)")
  
  day_values_local <- sort(unique(ratio_data$Day))
  min_gap <- if (length(day_values_local) > 1) {
    min(diff(day_values_local))
  } else {
    14
  }
  dodge_w  <- 0.4  * min_gap
  jitter_w <- 0.15 * min_gap
  ci_w     <- 0.18 * min_gap
  
  # ---- Title uses FULL marker names; y-axis is just "Ratio" ----
  title_txt <- paste0(
    cell_type, ": ",
    numerator_label,
    "  /  ",
    denominator_label
  )
  y_txt <- "Ratio"
  
  # ---- Build a plot for a given y-scale ----
  build_plot <- function(dat, scale_type = c("linear", "log10")) {
    scale_type <- match.arg(scale_type)
    
    p <- ggplot(
      dat,
      aes(x = Day, y = Ratio, color = Vaccine)
    ) +
      geom_point(
        position = position_jitterdodge(
          jitter.width  = jitter_w,
          jitter.height = 0,
          dodge.width   = dodge_w
        ),
        size  = 1.8,
        alpha = 0.6,
        na.rm = TRUE
      ) +
      stat_summary(
        aes(group = Vaccine),
        fun.data    = function(x) boot_mean_ci(x, n_boot = n_boot),
        geom        = "errorbar",
        width       = ci_w,
        linewidth   = 1.6,
        position    = position_dodge(width = dodge_w),
        na.rm       = TRUE,
        show.legend = FALSE
      ) +
      stat_summary(
        aes(group = Vaccine),
        fun      = mean,
        geom     = "point",
        size     = 2.5,
        position = position_dodge(width = dodge_w),
        na.rm    = TRUE
      ) +
      stat_summary(
        aes(group = Vaccine),
        fun       = mean,
        geom      = "line",
        linewidth = 1.8,
        position  = position_dodge(width = dodge_w),
        na.rm     = TRUE
      ) +
      labs(
        title = if (scale_type == "log10")
          paste0(title_txt, "  (log10)")
        else
          title_txt,
        x     = "Day",
        y     = if (scale_type == "log10") "Ratio (log10 scale)" else y_txt,
        color = "Vaccine"
      ) +
      scale_x_continuous(
        breaks = sort(unique(dat$Day)),
        expand = expansion(mult = c(0.03, 0.05))
      ) +
      scale_color_manual(
        values       = formulation_cols,
        breaks       = vaccine_levels,
        na.translate = FALSE
      ) +
      theme_bw(base_size = 16) +
      theme(
        plot.title         = element_text(size = 15, face = "bold"),
        axis.title         = element_text(size = 18, face = "bold"),
        axis.text          = element_text(size = 15),
        legend.title       = element_text(size = 16),
        legend.text        = element_text(size = 14),
        legend.position    = "right",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank()
      )
    
    if (scale_type == "log10") {
      p <- p + scale_y_log10(
        breaks = scales::breaks_log(n = 6),
        labels = scales::label_log(base = 10),
        expand = expansion(mult = c(0.05, 0.05))
      )
    } else {
      p <- p + scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.05))
      )
    }
    p
  }
  
  # ---- File names ----
  if (is.null(file_stub)) {
    safe <- function(s) gsub("[^A-Za-z0-9]+", "_", s)
    file_stub <- paste0(
      cell_type, "_", safe(numerator_short),
      "_over_", safe(denominator_short), "_ratio"
    )
  }
  
  out_dir <- file.path(figure_folder, "Ratio")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- 1) Linear version ----
  p_lin <- build_plot(ratio_data, "linear")
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, ".png")),
    plot     = p_lin,
    width    = 10,
    height   = 6.5,
    dpi      = 300
  )
  
  # ---- 2) log10 version (requires strictly positive ratios) ----
  ratio_data_log <- ratio_data %>%
    filter(is.finite(Ratio), Ratio > 0)
  
  n_dropped <- nrow(ratio_data) - nrow(ratio_data_log)
  if (n_dropped > 0) {
    message("   log10: dropped ", n_dropped,
            " non-positive ratio value(s)")
  }
  
  p_log <- NULL
  if (nrow(ratio_data_log) > 0) {
    p_log <- build_plot(ratio_data_log, "log10")
    ggsave(
      filename = file.path(out_dir, paste0(file_stub, "_log10.png")),
      plot     = p_log,
      width    = 10,
      height   = 6.5,
      dpi      = 300
    )
  } else {
    warning("SKIP log10 for ", numerator_short, "/", denominator_short,
            ": no positive ratio values")
  }
  
  # ---- Summary table (linear-scale summary) ----
  ratio_summary <- ratio_data %>%
    group_by(Day, Vaccine) %>%
    summarise(
      n       = sum(!is.na(Ratio)),
      mean    = mean(Ratio, na.rm = TRUE),
      sd      = sd(Ratio,   na.rm = TRUE),
      ci_low  = boot_mean_ci(Ratio, n_boot = n_boot)$ymin,
      ci_high = boot_mean_ci(Ratio, n_boot = n_boot)$ymax,
      .groups = "drop"
    )
  
  write.csv(
    ratio_summary,
    file.path(out_dir, paste0(file_stub, "_summary.csv")),
    row.names = FALSE
  )
  
  invisible(list(plot       = p_lin,
                 plot_log10 = p_log,
                 data       = ratio_data,
                 data_log10 = ratio_data_log,
                 summary    = ratio_summary))
}


# ==================================================================
# 3. Run all the ratios
# ==================================================================

# ------------------------------------------------------------------
# CD4
# ------------------------------------------------------------------

plot_ratio(
  numerator_label   = "CD44+ CD62L- TEM (% tetramer)",
  denominator_label = "CD44+ CD62L+ TCM (% tetramer)",
  numerator_short   = "TEM",
  denominator_short = "TCM",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_TEM_over_TCM_ratio"
)

plot_ratio(
  numerator_label   = "CXCR3- CCR6- Th2 (% tetramer)",
  denominator_label = "CXCR3+ CCR6- Th1 (% tetramer)",
  numerator_short   = "Th2",
  denominator_short = "Th1",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_Th2_over_Th1_ratio"
)

plot_ratio(
  numerator_label   = "CXCR3+ CCR6- Th1 (% tetramer)",
  denominator_label = "CXCR3- CCR6+ Th17 (% tetramer)",
  numerator_short   = "Th1",
  denominator_short = "Th17",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_Th1_over_Th17_ratio"
)

plot_ratio(
  numerator_label   = "CXCR3- CCR6- Th2 (% tetramer)",
  denominator_label = "CXCR3- CCR6+ Th17 (% tetramer)",
  numerator_short   = "Th2",
  denominator_short = "Th17",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_Th2_over_Th17_ratio"
)

plot_ratio(
  numerator_label   = "PD-1+ CXCR5+ BCL6+ TFH (% tetramer)",
  denominator_label = "Ly6C+ FOXP3- Teff (% tetramer)",
  numerator_short   = "TFH",
  denominator_short = "Teff",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_TFH_over_Teff_ratio"
)

plot_ratio(
  numerator_label   = "CXCR3+ CCR6- Th1 (% tetramer)",
  denominator_label = "PD-1+ CXCR5+ BCL6+ TFH (% tetramer)",
  numerator_short   = "Th1",
  denominator_short = "TFH",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_Th1_over_TFH_ratio"
)

plot_ratio(
  numerator_label   = "Ly6C+ FOXP3- Teff (% tetramer)",
  denominator_label = "Tetramer+  (% CD4)",
  numerator_short   = "Teff",
  denominator_short = "Tet",
  cell_type         = "CD4", unit = "percentage",
  file_stub         = "CD4_Teff_over_Tetramer_ratio"
)

# ------------------------------------------------------------------
# CD8
# ------------------------------------------------------------------

plot_ratio(
  numerator_label   = "CD127+ KLRG1- memory (% tetramer)",
  denominator_label = "KLRG1+ CD127- effector (% tetramer)",
  numerator_short   = "Memory",
  denominator_short = "Effector",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_Memory_over_Effector_ratio"
)

plot_ratio(
  numerator_label   = "CD127+ KLRG1- TCF1+ SLAMF6+ Tsl (% tetramer)",
  denominator_label = "KLRG1+ CD127- effector (% tetramer)",
  numerator_short   = "Tsl",
  denominator_short = "TermEff",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_Tsl_over_TerminalEffector_ratio"
)

plot_ratio(
  numerator_label   = "CD62L+ CX3CR1- TCM (% memory)",
  denominator_label = "CD62L- CX3CR1hi TEM (% memory)",
  numerator_short   = "TCM",
  denominator_short = "TEM",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_TCM_over_TEM_ratio"
)

plot_ratio(
  numerator_label   = "CD127+ TCF1+ SLAMF6+ CD62L+ Tscm (% tetramer)",
  denominator_label = "CD62L- CX3CR1hi TEM (% tetramer)",
  numerator_short   = "Tscm",
  denominator_short = "TEM",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_Tscm_over_TEM_ratio"
)

plot_ratio(
  numerator_label   = "CD127+ KLRG1- TCF1+ SLAMF6+ Tsl (% tetramer)",
  denominator_label = "CD127+ TCF1+ SLAMF6+ CD62L+ Tscm (% tetramer)",
  numerator_short   = "Tsl",
  denominator_short = "Tscm",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_Tsl_over_Tscm_ratio"
)

plot_ratio(
  numerator_label   = "PD-1+ CXCR5+ TFH like (% tetramer)",
  denominator_label = "CD127+ KLRG1- memory (% tetramer)",
  numerator_short   = "TFHlike",
  denominator_short = "Memory",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_TFHlike_over_Memory_ratio"
)

plot_ratio(
  numerator_label   = "CD127+ KLRG1- TCF1+ SLAMF6+ Tsl (% tetramer)",
  denominator_label = "Tetramer+ (% of CD8)",
  numerator_short   = "Tsl",
  denominator_short = "Tet",
  cell_type         = "CD8", unit = "percentage",
  file_stub         = "CD8_Tsl_over_Tetramer_ratio"
)