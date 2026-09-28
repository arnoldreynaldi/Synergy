# =========================================================================
# Libraries
# =========================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
subject_to_exclude <- c(subject_to_exclude, "11171", "10236")
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
  if (is.null(dat) || nrow(dat) < 4) {
    message("  [", phase_label, "] skipped: < 4 rows")
    return(NULL)
  }
  dat$Vaccine <- droplevels(factor(dat$Vaccine))
  if (nlevels(dat$Vaccine) < 2) {
    message("  [", phase_label, "] skipped: < 2 vaccine levels")
    return(NULL)
  }
  n_per_vax <- table(dat$Vaccine)
  if (any(n_per_vax < 2)) {
    message("  [", phase_label, "] skipped: some vaccine < 2 obs")
    return(NULL)
  }
  fit <- tryCatch(lm(log10_Value ~ Vaccine * Day, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) {
    message("  [", phase_label, "] skipped: lm failed")
    return(NULL)
  }
  if (isTRUE(summary(fit)$sigma == 0)) {
    message("  [", phase_label, "] skipped: zero residual variance")
    return(NULL)
  }
  if (any(is.na(coef(fit)))) {
    message("  [", phase_label, "] skipped: NA coefficients")
    return(NULL)
  }
  
  emm_slopes <- tryCatch(
    emtrends(fit, ~ Vaccine, var = "Day"),
    error = function(e) NULL
  )
  if (is.null(emm_slopes)) {
    message("  [", phase_label, "] skipped: emtrends failed")
    return(NULL)
  }
  
  slope_raw <- emm_clean(emm_slopes)
  needed    <- c("Vaccine", "Day.trend", "SE", "lower.CL", "upper.CL")
  if (!all(needed %in% names(slope_raw))) {
    message("  [", phase_label, "] skipped: emtrends columns ",
            paste(names(slope_raw), collapse = ", "))
    return(NULL)
  }
  
  slope_df <- slope_raw %>%
    mutate(
      Slope          = Day.trend,
      Slope_SE       = SE,
      Slope_Lower    = lower.CL,
      Slope_Upper    = upper.CL,
      Phase          = phase_label,
      HalfLife       = halflife_from_slope(Slope),
      HalfLife_Lower = ifelse(Slope_Lower < 0,
                              halflife_from_slope(Slope_Lower), NA_real_),
      HalfLife_Upper = ifelse(Slope_Upper < 0,
                              halflife_from_slope(Slope_Upper), NA_real_)
    ) %>%
    select(Vaccine, Phase, Slope, Slope_SE, Slope_Lower, Slope_Upper,
           HalfLife, HalfLife_Lower, HalfLife_Upper)
  
  # ---- pairwise p-values: defensive column lookup ----
  pairs_obj <- tryCatch(
    pairs(emm_slopes, adjust = "tukey"),
    error = function(e) {
      message("  [", phase_label, "] pairs() failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(pairs_obj)) {
    return(list(slopes = slope_df, pairs = NULL, model = fit))
  }
  
  pairs_df <- tryCatch({
    t <- emm_clean(as.data.frame(pairs_obj, stringsAsFactors = FALSE))
    ccol <- intersect(c("contrast", "term"), names(t))
    if (length(ccol) == 0) stop("no contrast/term column; have: ",
                                paste(names(t), collapse = ", "))
    ccol <- ccol[1]
    se_col <- intersect(c("SE", "std.error"), names(t))
    se_col <- if (length(se_col)) se_col[1] else NA_character_
    p_col  <- intersect(c("p.value", "p_value"), names(t))
    p_col  <- if (length(p_col)) p_col[1] else NA_character_
    has_lc <- "lower.CL" %in% names(t)
    has_uc <- "upper.CL" %in% names(t)
    n <- nrow(t)
    data.frame(
      Phase      = rep(phase_label, n),
      Contrast   = t[[ccol]],
      Diff_Slope = t$estimate,
      Diff_SE    = if (!is.na(se_col)) t[[se_col]] else rep(NA_real_, n),
      Diff_Lower = if (has_lc) t$lower.CL else rep(NA_real_, n),
      Diff_Upper = if (has_uc) t$upper.CL else rep(NA_real_, n),
      p_value    = if (!is.na(p_col)) t[[p_col]] else rep(NA_real_, n),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("  [", phase_label, "] pairs extraction failed: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(pairs_df) || nrow(pairs_df) == 0) {
    return(list(slopes = slope_df, pairs = NULL, model = fit))
  }
  
  list(slopes = slope_df, pairs = pairs_df, model = fit)
}

# =========================================================================
# Half-life plot with pairwise p-value brackets
# (significant p-values, p < 0.05, are printed in BOLD)
# =========================================================================
make_halflife_pvalue_plot <- function(per_vaccine_hl, emm_pairs) {
  
  vaccine_order <- levels(per_vaccine_hl$Vaccine)
  if (is.null(vaccine_order)) vaccine_order <- sort(unique(per_vaccine_hl$Vaccine))
  vx <- setNames(seq_along(vaccine_order), vaccine_order)
  
  per_vaccine_hl <- per_vaccine_hl %>%
    mutate(Phase   = as.character(Phase),
           Vaccine = factor(Vaccine, levels = vaccine_order))
  
  phase_top <- per_vaccine_hl %>%
    group_by(Phase) %>%
    summarise(data_top = max(HalfLife_Upper_plot, na.rm = TRUE),
              .groups = "drop")
  
  if (!is.null(emm_pairs) && nrow(emm_pairs) > 0) {
    br <- emm_pairs %>%
      mutate(Phase = as.character(Phase)) %>%
      separate(Contrast, into = c("g1", "g2"), sep = "\\s*-\\s*", remove = FALSE) %>%
      mutate(g1 = trimws(g1), g2 = trimws(g2),
             x1 = vx[g1], x2 = vx[g2],
             xmin = pmin(x1, x2), xmax = pmax(x1, x2),
             p_label = sprintf("%.3f", p_value),
             p_bold  = ifelse(p_value < 0.05, "bold", "plain")) %>%
      filter(!is.na(xmin), !is.na(xmax)) %>%
      left_join(phase_top, by = "Phase") %>%
      group_by(Phase) %>%
      arrange(xmin, xmax, .by_group = TRUE) %>%
      mutate(bracket_i = row_number(),
             y = data_top * (1 + 0.13 * bracket_i)) %>%
      ungroup()
  } else {
    br <- NULL
  }
  
  if (!is.null(br) && nrow(br) > 0) {
    ytops <- br %>%
      group_by(Phase) %>%
      summarise(ytop = max(y) * 1.12, .groups = "drop")
  } else {
    ytops <- phase_top %>% mutate(ytop = data_top * 1.15)
  }
  
  dummy <- data.frame(
    Phase   = ytops$Phase,
    Vaccine = factor(vaccine_order[1], levels = vaccine_order),
    y       = ytops$ytop
  )
  
  p <- ggplot() +
    geom_blank(data = dummy, aes(x = Vaccine, y = y)) +
    geom_point(data = per_vaccine_hl,
               aes(x = Vaccine, y = HalfLife, color = Vaccine), size = 3) +
    geom_errorbar(data = per_vaccine_hl,
                  aes(x = Vaccine, ymin = HalfLife_Lower,
                      ymax = HalfLife_Upper_plot, color = Vaccine),
                  width = 0.2) +
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
      geom_text(data = br,
                aes(x = (xmin + xmax) / 2, y = y, label = p_label,
                    fontface = p_bold),
                inherit.aes = FALSE, vjust = -0.3, size = 3.2)
  }
  
  p
}

# =========================================================================
# Main analysis function
# =========================================================================
run_analysis <- function(data, cell_type, measurement_type) {
  
  message(">>> ", cell_type, " / ", measurement_type)
  
  subset_df <- data %>%
    filter(Cell_Type == cell_type,
           Measurement_type == measurement_type,
           Unit == "number") %>%
    mutate(log10_Value = log10(Value))
  
  subset_df <- subset_df %>% filter(is.finite(log10_Value))
  subset_df <- subset_df %>% filter(!Mouse_number %in% subject_to_exclude)
  
  n_per_vaccine <- subset_df %>% group_by(Vaccine) %>% summarise(n = n())
  if (any(n_per_vaccine$n < 1)) return(NULL)
  
  subset_df$Vaccine <- factor(subset_df$Vaccine, levels = vaccine_levels)
  vaccines <- levels(droplevels(subset_df$Vaccine))
  if (length(vaccines) < 2) return(NULL)
  
  # ---- 1. Group-level two-point decay rates -------------------------------
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
  
  # ---- 2. Per-vaccine first-2 / last-2 days -------------------------------
  early_data <- subset_df %>%
    group_by(Vaccine) %>%
    filter(Day %in% sort(unique(Day))[1:2]) %>%
    ungroup()
  
  late_data <- subset_df %>%
    group_by(Vaccine) %>%
    filter(Day %in% sort(unique(Day))[
      (length(unique(Day)) - 1):length(unique(Day))]) %>%
    ungroup()
  
  # ---- 3. emmeans per phase ----------------------------------------------
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
  
  # ---- 4. Anchor segments -------------------------------------------------
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
  
  # ---- 5. Plots -----------------------------------------------------------
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
  
  per_vaccine_hl <- per_vaccine_decay %>%
    mutate(HalfLife_Upper_plot = ifelse(is.na(HalfLife_Upper) & HalfLife > 0,
                                        300, HalfLife_Upper))
  
  plot_halflife_pvalue <- make_halflife_pvalue_plot(per_vaccine_hl, emm_pairs)
  
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
  
  # ---- 6. Save ------------------------------------------------------------
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
  
  return(list(
    group_decay_two_point = group_decay_two_point,
    per_vaccine_decay     = per_vaccine_decay,
    emm_pairs             = emm_pairs
  ))
}

# =========================================================================
# Main loop
# =========================================================================
combos <- combined_data %>%
  filter(Unit == "number", Cell_Type %in% c("CD4", "CD8")) %>%
  distinct(Cell_Type, Measurement_type) %>%
  arrange(Cell_Type, Measurement_type)

results_list <- list()
for (i in 1:nrow(combos)) {
  ct <- combos$Cell_Type[i]
  mt <- combos$Measurement_type[i]
  res <- tryCatch(
    run_analysis(combined_data, ct, mt),
    error = function(e) {
      message("Skipped: ", ct, " / ", mt, " -- ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(res)) {
    results_list[[paste(ct, mt, sep = " | ")]] <- res
  }
}