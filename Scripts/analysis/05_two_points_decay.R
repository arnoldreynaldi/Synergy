# =========================================================================
# Helpers
# =========================================================================
sanitize <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

big_text_theme <- theme_bw(base_size = 16) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(size = 17, face = "bold", margin = margin(b = 6)),
    strip.text       = element_text(size = 15, face = "bold", lineheight = 1.1),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(size = 13),
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 13),
    plot.margin      = margin(5, 5, 5, 5)
  )

halflife_from_slope <- function(slope) -log10(2) / slope

group_two_point_slope <- function(dat, phase = c("early", "late")) {
  phase <- match.arg(phase)
  d <- sort(unique(dat$Day))
  if (length(d) < 2) return(NA_real_)
  idx <- if (phase == "early") c(1, 2) else c(length(d) - 1, length(d))
  d1 <- d[idx[1]]; d2 <- d[idx[2]]
  y1 <- mean(dat$log10_Value[dat$Day == d1], na.rm = TRUE)
  y2 <- mean(dat$log10_Value[dat$Day == d2], na.rm = TRUE)
  (y2 - y1) / (d2 - d1)
}

# Normalise emmeans column names
emm_clean <- function(emm) {
  d  <- as.data.frame(emm, stringsAsFactors = FALSE)
  nm <- names(d)
  nm <- sub("^asymp\\.LCL$", "lower.CL", nm, ignore.case = TRUE)
  nm <- sub("^asymp\\.UCL$", "upper.CL", nm, ignore.case = TRUE)
  nm <- sub("^LCL$",         "lower.CL", nm, ignore.case = TRUE)
  nm <- sub("^UCL$",         "upper.CL", nm, ignore.case = TRUE)
  nm <- sub("^conf\\.low$",  "lower.CL", nm, ignore.case = TRUE)
  nm <- sub("^conf\\.high$", "upper.CL", nm, ignore.case = TRUE)
  nm <- sub("^std\\.error$", "SE",       nm, ignore.case = TRUE)
  nm <- sub("^std_error$",   "SE",       nm, ignore.case = TRUE)
  nm <- sub("^p_value$",     "p.value",  nm, ignore.case = TRUE)
  nm <- sub("^pvalue$",      "p.value",  nm, ignore.case = TRUE)
  names(d) <- nm
  d
}

# =========================================================================
# Per-phase emmeans
# =========================================================================
emm_one_phase <- function(dat, phase_label) {
  if (is.null(dat) || nrow(dat) < 4) return(NULL)
  dat$Vaccine <- droplevels(factor(dat$Vaccine))
  if (nlevels(dat$Vaccine) < 2) return(NULL)
  n_per_vax <- table(dat$Vaccine)
  if (any(n_per_vax < 2)) return(NULL)
  
  fit <- tryCatch(lm(log10_Value ~ Vaccine * Day, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  if (!is.finite(summary(fit)$sigma) || summary(fit)$sigma < 1e-10) return(NULL)
  if (any(is.na(coef(fit)))) return(NULL)
  
  emm_slopes <- tryCatch(
    emtrends(fit, ~ Vaccine, var = "Day"),
    error = function(e) NULL
  )
  if (is.null(emm_slopes)) return(NULL)
  
  slope_raw <- emm_clean(emm_slopes)
  needed    <- c("Vaccine", "Day.trend", "SE", "lower.CL", "upper.CL")
  if (!all(needed %in% names(slope_raw))) return(NULL)
  
  # Half-life + CI.
  #   halflife_from_slope() is INCREASING in slope, so for a slope CI entirely
  #   on one side of 0:  HL_lower = f(Slope_lower),  HL_upper = f(Slope_upper).
  #   If the slope CI straddles 0 -> half-life CI is not an interval; suppress it.
  slope_df <- slope_raw %>%
    mutate(
      Slope          = Day.trend,
      Slope_SE       = SE,
      Slope_Lower    = lower.CL,
      Slope_Upper    = upper.CL,
      Phase          = phase_label,
      HalfLife       = ifelse(abs(Slope) < 1e-10, NA_real_,
                              halflife_from_slope(Slope)),
      .ci_ok         = Slope_Lower < 0 & Slope_Upper < 0 |
        Slope_Lower >= 0 & Slope_Upper >= 0,
      .ci_ok         = .ci_ok & (Slope_Lower < 0 | Slope_Upper < 0),
      HalfLife_Lower = ifelse(.ci_ok, halflife_from_slope(Slope_Lower), NA_real_),
      HalfLife_Upper = ifelse(.ci_ok, halflife_from_slope(Slope_Upper), NA_real_)
    ) %>%
    select(Vaccine, Phase, Slope, Slope_SE, Slope_Lower, Slope_Upper,
           HalfLife, HalfLife_Lower, HalfLife_Upper)
  
  # ---- pairwise p-values ----
  pairs_obj <- tryCatch(pairs(emm_slopes, adjust = "tukey"),
                        error = function(e) NULL)
  if (is.null(pairs_obj)) {
    return(list(slopes = slope_df, pairs = NULL, model = fit))
  }
  
  pairs_df <- tryCatch({
    t <- emm_clean(as.data.frame(pairs_obj, stringsAsFactors = FALSE))
    ccol <- intersect(c("contrast", "term"), names(t))
    if (length(ccol) == 0) stop("no contrast/term column")
    ccol <- ccol[1]
    p_col <- intersect(c("p.value", "p_value"), names(t))
    p_col <- if (length(p_col)) p_col[1] else NA_character_
    data.frame(
      Phase    = phase_label,
      Contrast = t[[ccol]],
      p_value  = if (!is.na(p_col)) t[[p_col]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
  
  if (is.null(pairs_df) || nrow(pairs_df) == 0) {
    return(list(slopes = slope_df, pairs = NULL, model = fit))
  }
  
  list(slopes = slope_df, pairs = pairs_df, model = fit)
}

# =========================================================================
# Half-life plot with pairwise p-value brackets.
# p < 0.05 -> bold red label.  Vaccine colour scheme preserved on points/CI.
# =========================================================================
make_halflife_pvalue_plot <- function(per_vaccine_hl, emm_pairs) {
  
  vaccine_order <- levels(per_vaccine_hl$Vaccine)
  if (is.null(vaccine_order))
    vaccine_order <- sort(unique(per_vaccine_hl$Vaccine))
  vx <- setNames(seq_along(vaccine_order), vaccine_order)
  
  per_vaccine_hl <- per_vaccine_hl %>%
    mutate(Phase   = as.character(Phase),
           Vaccine = factor(Vaccine, levels = vaccine_order)) %>%
    filter(is.finite(HalfLife))
  
  phase_top <- per_vaccine_hl %>%
    group_by(Phase) %>%
    summarise(
      data_top = suppressWarnings(max(c(HalfLife, HalfLife_Upper), na.rm = TRUE)),
      data_bot = suppressWarnings(min(c(HalfLife, HalfLife_Lower), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(pad = pmax((data_top - data_bot) * 0.1, 1e-6))
  
  if (!is.null(emm_pairs) && nrow(emm_pairs) > 0) {
    br <- emm_pairs %>%
      mutate(Phase = as.character(Phase)) %>%
      separate(Contrast, into = c("g1", "g2"),
               sep = "\\s*-\\s*", remove = FALSE) %>%
      mutate(g1 = trimws(g1), g2 = trimws(g2),
             x1 = vx[g1], x2 = vx[g2],
             xmin = pmin(x1, x2), xmax = pmax(x1, x2),
             p_label = ifelse(is.na(p_value), "", sprintf("%.3f", p_value))) %>%
      filter(!is.na(xmin), !is.na(xmax)) %>%
      left_join(phase_top, by = "Phase") %>%
      group_by(Phase) %>%
      arrange(xmin, xmax, .by_group = TRUE) %>%
      mutate(bracket_i = row_number(),
             y = data_top + pad * (1 + bracket_i)) %>%
      ungroup()
  } else {
    br <- NULL
  }
  
  br_sig <- if (!is.null(br) && nrow(br) > 0) {
    br %>% filter(!is.na(p_value), p_value < 0.05)
  } else NULL
  
  headroom <- phase_top %>%
    mutate(y = data_top + pad * (if (is.null(br) || nrow(br) == 0) 1.5
                                 else max(br$bracket_i[br$Phase == Phase]) + 2)) %>%
    transmute(Phase,
              Vaccine = factor(vaccine_order[1], levels = vaccine_order),
              y)
  
  p <- ggplot() +
    geom_blank(data = headroom, aes(x = Vaccine, y = y)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(data = per_vaccine_hl,
               aes(x = Vaccine, y = HalfLife, color = Vaccine), size = 3) +
    geom_errorbar(data = per_vaccine_hl,
                  aes(x = Vaccine, ymin = HalfLife_Lower,
                      ymax = HalfLife_Upper, color = Vaccine),
                  width = 0.2, na.rm = TRUE) +
    facet_wrap(~ Phase, scales = "free_y") +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Half-life (days) with 95% CI\nand adjusted pairwise p-values",
         y = "Half-life (days)", x = NULL) +
    big_text_theme
  
  if (!is.null(br) && nrow(br) > 0) {
    p <- p +
      geom_segment(data = br,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, color = "grey30", linewidth = 0.4) +
      # Base layer: all p-labels in plain black
      geom_text(data = br,
                aes(x = (xmin + xmax) / 2, y = y, label = p_label),
                inherit.aes = FALSE, vjust = -0.3, size = 3.2,
                color = "black", fontface = "plain") +
      # Overlay: significant p-labels in bold red on top
      geom_text(data = br_sig,
                aes(x = (xmin + xmax) / 2, y = y, label = p_label),
                inherit.aes = FALSE, vjust = -0.3, size = 3.2,
                color = "red", fontface = "bold")
  }
  
  p
}

# =========================================================================
# Main analysis function
# =========================================================================
run_analysis <- function(data, cell_type, measurement_type) {
  
  subset_df <- data %>%
    filter(Cell_Type == cell_type,
           Measurement_type == measurement_type,
           Unit == "number") %>%
    mutate(log10_Value = log10(Value)) %>%
    filter(is.finite(log10_Value)) %>%
    filter(!Mouse_number %in% subject_to_exclude)
  
  if (nrow(subset_df) == 0) return(NULL)
  
  subset_df$Vaccine <- factor(subset_df$Vaccine, levels = vaccine_levels)
  vaccines <- levels(droplevels(subset_df$Vaccine))
  if (length(vaccines) < 2) return(NULL)
  
  group_decay_two_point <- subset_df %>%
    group_by(Vaccine) %>%
    group_modify(~data.frame(
      early_slope = group_two_point_slope(.x, "early"),
      late_slope  = group_two_point_slope(.x, "late")
    )) %>%
    ungroup() %>%
    mutate(
      early_halflife = halflife_from_slope(early_slope),
      late_halflife  = halflife_from_slope(late_slope)
    )
  
  early_data <- subset_df %>%
    group_by(Vaccine) %>%
    filter(Day %in% sort(unique(Day))[1:2]) %>%
    ungroup()
  
  late_data <- subset_df %>%
    group_by(Vaccine) %>%
    filter(Day %in% sort(unique(Day))[
      (length(unique(Day)) - 1):length(unique(Day))]) %>%
    ungroup()
  
  early_emm <- emm_one_phase(early_data, "Early")
  late_emm  <- emm_one_phase(late_data,  "Late")
  
  if (is.null(early_emm) && is.null(late_emm)) return(NULL)
  
  per_vaccine_decay <- dplyr::bind_rows(
    if (!is.null(early_emm)) early_emm$slopes else NULL,
    if (!is.null(late_emm))  late_emm$slopes  else NULL
  )
  
  emm_pairs <- dplyr::bind_rows(
    if (!is.null(early_emm) && !is.null(early_emm$pairs)) early_emm$pairs else NULL,
    if (!is.null(late_emm)  && !is.null(late_emm$pairs))  late_emm$pairs  else NULL
  )
  
  if (nrow(per_vaccine_decay) == 0) return(NULL)
  
  anchor_days <- subset_df %>%
    group_by(Vaccine) %>%
    summarise(
      d1  = sort(unique(Day))[1],
      d2  = sort(unique(Day))[2],
      dn1 = sort(unique(Day))[length(unique(Day)) - 1],
      dn  = sort(unique(Day))[length(unique(Day))],
      .groups = "drop"
    )
  
  anchor_means <- subset_df %>%
    left_join(anchor_days, by = "Vaccine") %>%
    filter(Day %in% c(d1, d2, dn1, dn)) %>%
    mutate(Phase = ifelse(Day %in% c(d1, d2), "Early", "Late")) %>%
    group_by(Vaccine, Phase, Day) %>%
    summarise(mean_log10 = mean(log10_Value, na.rm = TRUE), .groups = "drop")
  
  slope_lookup <- per_vaccine_decay %>%
    select(Vaccine, Phase, Slope)
  
  segments_df <- anchor_means %>%
    left_join(slope_lookup, by = c("Vaccine", "Phase")) %>%
    group_by(Vaccine, Phase) %>%
    arrange(Day) %>%
    summarise(
      x1 = Day[1],
      x2 = Day[2],
      y1 = mean_log10[1],
      y2 = mean_log10[1] + Slope[1] * (Day[2] - Day[1]),
      .groups = "drop"
    )
  
  raw_anchors <- subset_df %>%
    left_join(anchor_days, by = "Vaccine") %>%
    filter(Day %in% c(d1, d2, dn1, dn)) %>%
    mutate(Phase = ifelse(Day %in% c(d1, d2), "Early", "Late"))
  
  plot_group_slopes <- ggplot(
    tidyr::pivot_longer(group_decay_two_point,
                        cols = c(early_slope, late_slope),
                        names_to = "Phase", values_to = "Slope") %>%
      mutate(Phase = ifelse(Phase == "early_slope", "Early", "Late")),
    aes(x = Vaccine, y = Slope, color = Vaccine)) +
    geom_point(size = 3) +
    facet_wrap(~Phase) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Group-level two-point\ndecay slopes",
         y = "Slope (log10 per day)", x = NULL) +
    big_text_theme
  
  plot_group_halflife <- ggplot(
    tidyr::pivot_longer(group_decay_two_point,
                        cols = c(early_halflife, late_halflife),
                        names_to = "Phase", values_to = "HalfLife") %>%
      mutate(Phase = ifelse(Phase == "early_halflife", "Early", "Late")),
    aes(x = Vaccine, y = HalfLife, color = Vaccine)) +
    geom_point(size = 3) +
    facet_wrap(~Phase) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Group-level two-point half-life\n(days, log10 scale)",
         y = "Half-life (days)", x = NULL) +
    big_text_theme
  
  plot_halflife_pvalue <- make_halflife_pvalue_plot(per_vaccine_decay, emm_pairs)
  
  plot_anchor_segments <- ggplot() +
    geom_point(data = raw_anchors,
               aes(x = Day, y = log10_Value, color = Vaccine),
               alpha = 0.6, size = 1.6) +
    geom_segment(data = segments_df,
                 aes(x = x1, xend = x2, y = y1, yend = y2, color = Vaccine),
                 linewidth = 1.3) +
    facet_wrap(~Phase, scales = "free_x") +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Datapoints with\nfitted slope",
         y = "log10(Cell count)", x = "Day") +
    big_text_theme
  
  plot_base  <- file.path(figure_folder, paste0(cell_type, "_analysis"),
                          sanitize(measurement_type))
  table_base <- file.path(table_folder,  paste0(cell_type, "_analysis"),
                          sanitize(measurement_type))
  dir.create(plot_base,  recursive = TRUE, showWarnings = FALSE)
  dir.create(table_base, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(plot_base, "1_Group_slopes_two_point.png"),
         plot_group_slopes,     width = 5.5, height = 3.2, dpi = 300)
  ggsave(file.path(plot_base, "2_Group_halflife_two_point.png"),
         plot_group_halflife,   width = 5.5, height = 3.2, dpi = 300)
  ggsave(file.path(plot_base, "3_Halflife_with_CI_and_pvalues.png"),
         plot_halflife_pvalue,  width = 7,   height = 4.2, dpi = 300)
  ggsave(file.path(plot_base, "5_AnchorPoints_emm_slope.png"),
         plot_anchor_segments,  width = 6,   height = 3.2, dpi = 300)
  
  write.csv(group_decay_two_point,
            file.path(table_base, "group_decay_two_point.csv"),
            row.names = FALSE)
  write.csv(per_vaccine_decay,
            file.path(table_base, "emmeans_slopes_CI.csv"),
            row.names = FALSE)
  if (!is.null(emm_pairs) && nrow(emm_pairs) > 0) {
    write.csv(emm_pairs,
              file.path(table_base, "emmeans_pairwise_pvalues.csv"),
              row.names = FALSE)
  }
  
  list(
    group_decay_two_point = group_decay_two_point,
    per_vaccine_decay     = per_vaccine_decay,
    emm_pairs             = emm_pairs
  )
}

# =========================================================================
# Main loop
# =========================================================================
combos <- combined_data %>%
  filter(Unit == "number", Cell_Type %in% c("CD4", "CD8")) %>%
  distinct(Cell_Type, Measurement_type) %>%
  arrange(Cell_Type, Measurement_type)

results_list <- list()
for (i in seq_len(nrow(combos))) {
  ct <- combos$Cell_Type[i]
  mt <- combos$Measurement_type[i]
  res <- tryCatch(
    run_analysis(combined_data, ct, mt),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    results_list[[paste(ct, mt, sep = " | ")]] <- res
  }
}