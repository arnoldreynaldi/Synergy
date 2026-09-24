# =========================================================================
# Libraries
# =========================================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# Assumes these globals exist in your environment:
#   combined_data, subject_to_exclude, vaccine_levels, formulation_cols,
#   figure_folder, table_folder

# =========================================================================
# Helpers
# =========================================================================
sanitize <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

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

generate_partitions <- function(n) {
  if (n == 1) return(list(1L))
  result <- list()
  prev <- generate_partitions(n - 1)
  for (p in prev) {
    max_cluster <- max(p)
    for (c in 1:max_cluster) result <- c(result, list(c(p, c)))
    result <- c(result, list(c(p, max_cluster + 1L)))
  }
  result
}

format_group <- function(cl, vaccines) {
  if (is.null(cl)) return("")
  groups <- split(vaccines, cl)
  paste(sapply(groups, paste, collapse = " & "), collapse = " | ")
}

fit_slope_partition <- function(dat, cl, vaccines) {
  dat$Vaccine <- factor(dat$Vaccine, levels = vaccines)
  slope_cluster <- factor(cl[match(dat$Vaccine, vaccines)])
  for (g in levels(slope_cluster)) {
    dat[[paste0("Dg", g)]] <- ifelse(slope_cluster == g, dat$Day, 0)
  }
  fml <- as.formula(paste("log10_Value ~ Vaccine +",
                          paste(paste0("Dg", levels(slope_cluster)),
                                collapse = " + ")))
  mod <- tryCatch(lm(fml, data = dat), error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  list(aic = AIC(mod), k = length(coef(mod)), cl = cl, model = mod)
}

slope_partition_search <- function(dat, phase_label, vaccines) {
  if (is.null(dat) || nrow(dat) < 4 || length(vaccines) < 2) return(NULL)
  parts <- generate_partitions(length(vaccines))
  fits  <- lapply(parts, function(p) fit_slope_partition(dat, p, vaccines))
  keep  <- !vapply(fits, is.null, logical(1))
  fits  <- fits[keep]
  if (length(fits) == 0) return(NULL)
  
  aics <- vapply(fits, function(x) x$aic, numeric(1))
  ks   <- vapply(fits, function(x) x$k,   numeric(1))
  ord  <- order(aics)
  fits <- fits[ord]; aics <- aics[ord]; ks <- ks[ord]
  
  best_aic <- aics[1]
  delta    <- aics - best_aic
  eligible <- which(delta < 10)
  min_k    <- min(ks[eligible])
  cand     <- eligible[ks[eligible] == min_k]
  parsim_i <- cand[which.min(aics[cand])]
  
  tbl <- data.frame(
    Phase        = phase_label,
    Groups       = vapply(fits, function(x) format_group(x$cl, vaccines), character(1)),
    n_slope_grps = vapply(fits, function(x) length(unique(x$cl)), integer(1)),
    k            = ks,
    AIC          = aics,
    Delta_AIC    = delta,
    Is_Best      = seq_along(fits) == 1,
    Is_Parsim    = seq_along(fits) == parsim_i,
    row.names    = NULL
  )
  list(table = tbl, best = fits[[1]], parsimonious = fits[[parsim_i]])
}

halflife_with_ci <- function(dat, cl, vaccines, phase_label, model_label) {
  dat$Vaccine <- factor(dat$Vaccine, levels = vaccines)
  slope_cluster <- factor(cl[match(dat$Vaccine, vaccines)])
  for (g in levels(slope_cluster)) {
    dat[[paste0("Dg", g)]] <- ifelse(slope_cluster == g, dat$Day, 0)
  }
  fml <- as.formula(paste("log10_Value ~ Vaccine +",
                          paste(paste0("Dg", levels(slope_cluster)),
                                collapse = " + ")))
  mod <- lm(fml, data = dat)
  cf  <- coef(mod)
  vc  <- vcov(mod)
  
  do.call(rbind, lapply(seq_along(vaccines), function(i) {
    v  <- vaccines[i]
    g  <- as.character(cl[i])
    nm <- paste0("Dg", g)
    est <- cf[[nm]]
    se  <- sqrt(vc[nm, nm])
    lo  <- est - 1.96 * se
    hi  <- est + 1.96 * se
    data.frame(
      Phase          = phase_label,
      Model          = model_label,
      Vaccine        = v,
      Slope          = est,
      Slope_Lower    = lo,
      Slope_Upper    = hi,
      HalfLife       = halflife_from_slope(est),
      HalfLife_Lower = if (!is.na(lo) && lo < 0) halflife_from_slope(lo) else NA_real_,
      HalfLife_Upper = if (!is.na(hi) && hi < 0) halflife_from_slope(hi) else NA_real_,
      row.names      = NULL
    )
  }))
}

# =========================================================================
# Main analysis function
# =========================================================================
run_analysis <- function(data, cell_type, measurement_type) {
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
  if (length(vaccines) < 1) return(NULL)
  
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
  
  # ---- 3. AIC search: separately for early and late -----------------------
  early_search <- slope_partition_search(early_data, "Early (first 2 timepoints)", vaccines)
  late_search  <- slope_partition_search(late_data,  "Late  (last 2 timepoints)",  vaccines)
  
  aic_decay <- rbind(
    if (!is.null(early_search)) early_search$table else NULL,
    if (!is.null(late_search))  late_search$table  else NULL
  )
  
  # ---- 4. Per-vaccine slope / half-life / 95% CI --------------------------
  all_diff_cl <- seq_along(vaccines)
  
  per_vaccine_decay <- rbind(
    if (!is.null(early_search))
      halflife_with_ci(early_data, early_search$best$cl, vaccines,
                       "Early (first 2 timepoints)",
                       model_label = "Best AIC partition") else NULL,
    if (!is.null(late_search))
      halflife_with_ci(late_data, late_search$best$cl, vaccines,
                       "Late  (last 2 timepoints)",
                       model_label = "Best AIC partition") else NULL,
    
    if (!is.null(early_search))
      halflife_with_ci(early_data, all_diff_cl, vaccines,
                       "Early (first 2 timepoints)",
                       model_label = "All 4 different") else NULL,
    if (!is.null(late_search))
      halflife_with_ci(late_data, all_diff_cl, vaccines,
                       "Late  (last 2 timepoints)",
                       model_label = "All 4 different") else NULL
  )
  
  # =========================================================================
  # 4b. Segment: line between the two anchor timepoints,
  #     slope from the best-AIC model, anchored at mean of first timepoint
  # =========================================================================
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
    filter(Model == "Best AIC partition") %>%
    mutate(Phase = ifelse(grepl("^Early", Phase), "Early", "Late")) %>%
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
    ) %>%
    mutate(Phase = ifelse(Phase == "Early",
                          "Early (first 2 timepoints)",
                          "Late  (last 2 timepoints)"))
  
  raw_anchors <- subset_df %>%
    left_join(anchor_days, by = "Vaccine") %>%
    filter(Day %in% c(d1, d2, dn1, dn)) %>%
    mutate(Phase = ifelse(Day %in% c(d1, d2), "Early", "Late")) %>%
    mutate(Phase = ifelse(Phase == "Early",
                          "Early (first 2 timepoints)",
                          "Late  (last 2 timepoints)"))
  
  # =========================================================================
  # 5. Plots
  # =========================================================================
  
  plot_group_slopes <- ggplot(
    tidyr::pivot_longer(group_decay_two_point,
                        cols = c(early_slope, late_slope),
                        names_to = "Phase", values_to = "Slope") %>%
      mutate(Phase = ifelse(Phase == "early_slope",
                            "Early (first 2 timepoints)",
                            "Late  (last 2 timepoints)")),
    aes(x = Vaccine, y = Slope, color = Vaccine)) +
    geom_point(size = 3) +
    facet_wrap(~Phase) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Group-level two-point decay slopes",
         y = "Slope (log10 per day)", x = NULL) +
    theme_bw() + theme(legend.position = "none")
  
  plot_group_halflife <- ggplot(
    tidyr::pivot_longer(group_decay_two_point,
                        cols = c(early_halflife, late_halflife),
                        names_to = "Phase", values_to = "HalfLife") %>%
      mutate(Phase = ifelse(Phase == "early_halflife",
                            "Early (first 2 timepoints)",
                            "Late  (last 2 timepoints)")),
    aes(x = Vaccine, y = HalfLife, color = Vaccine)) +
    geom_point(size = 3) +
    facet_wrap(~Phase) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Group-level two-point half-life (days, log10 scale)",
         y = "Half-life (days)", x = NULL) +
    theme_bw() + theme(legend.position = "none")
  
  # ---- Half-life with CI; free_y so early and late get their own scale ----
  make_hl_plot <- function(df, title, group_lookup = NULL) {
    if (!is.null(group_lookup)) {
      df <- df %>%
        left_join(group_lookup, by = "Phase") %>%
        mutate(Phase_label = paste0(Phase, "\nslope grouping: ", Groups))
      p <- ggplot(df, aes(x = Vaccine, y = HalfLife, color = Vaccine)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = HalfLife_Lower, ymax = HalfLife_Upper),
                      width = 0.2) +
        facet_wrap(~Phase_label, scales = "free_y")
    } else {
      p <- ggplot(df, aes(x = Vaccine, y = HalfLife, color = Vaccine)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = HalfLife_Lower, ymax = HalfLife_Upper),
                      width = 0.2) +
        facet_wrap(~Phase, scales = "free_y")
    }
    p +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = title, y = "Half-life (days)", x = NULL) +
      theme_bw() + theme(legend.position = "none")
  }
  
  best_group_lookup <- aic_decay %>%
    filter(Is_Best) %>%
    select(Phase, Groups)
  
  plot_halflife_ci_best <- make_hl_plot(
    per_vaccine_decay %>% filter(Model == "Best AIC partition"),
    "Half-life (days) with 95% CI - from best-AIC model",
    group_lookup = best_group_lookup
  )
  
  plot_halflife_ci_alldiff <- make_hl_plot(
    per_vaccine_decay %>% filter(Model == "All 4 different"),
    "Half-life (days) with 95% CI - all 4 groups different"
  )
  
  # ---- Anchor points (no jitter) + best-AIC slope segment -----------------
  plot_anchor_segments <- ggplot() +
    geom_point(data = raw_anchors,
               aes(x = Day, y = log10_Value, color = Vaccine),
               alpha = 0.6, size = 1.6) +
    geom_segment(data = segments_df,
                 aes(x = x1, xend = x2, y = y1, yend = y2, color = Vaccine),
                 linewidth = 1.3) +
    facet_wrap(~Phase, scales = "free_x") +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Datapoints with best-AIC fitted slope",
         y = "log10(Cell count)", x = "Day") +
    theme_bw() + theme(legend.position = "none")
  
  # =========================================================================
  # 6. Save outputs
  # =========================================================================
  plot_base  <- file.path(figure_folder, paste0(cell_type, "_analysis"),
                          sanitize(measurement_type))
  table_base <- file.path(table_folder,  paste0(cell_type, "_analysis"),
                          sanitize(measurement_type))
  dir.create(plot_base,  recursive = TRUE, showWarnings = FALSE)
  dir.create(table_base, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(plot_base, "1_Group_slopes_two_point.png"),
         plot_group_slopes,        width = 7, height = 4)
  ggsave(file.path(plot_base, "2_Group_halflife_two_point.png"),
         plot_group_halflife,      width = 7, height = 4)
  ggsave(file.path(plot_base, "3_Halflife_with_CI_bestAIC.png"),
         plot_halflife_ci_best,    width = 8, height = 4)
  ggsave(file.path(plot_base, "4_Halflife_with_CI_all4diff.png"),
         plot_halflife_ci_alldiff, width = 7, height = 4)
  ggsave(file.path(plot_base, "5_AnchorPoints_bestAICslope.png"),
         plot_anchor_segments,     width = 8, height = 4)
  
  write.csv(group_decay_two_point,
            file.path(table_base, "group_decay_two_point.csv"), row.names = FALSE)
  write.csv(aic_decay,
            file.path(table_base, "AIC_decay_partitions.csv"),  row.names = FALSE)
  write.csv(per_vaccine_decay,
            file.path(table_base, "decay_per_vaccine_with_CI.csv"),
            row.names = FALSE)
  write.csv(
    per_vaccine_decay %>%
      select(Phase, Model, Vaccine, Slope, Slope_Lower, Slope_Upper,
             HalfLife, HalfLife_Lower, HalfLife_Upper),
    file.path(table_base, "halflife_CI_table.csv"),
    row.names = FALSE
  )
  
  return(aic_decay)
}

# =========================================================================
# Main loop
# =========================================================================
combos <- combined_data %>%
  filter(Unit == "number", Cell_Type %in% c("CD4", "CD8")) %>%
  distinct(Cell_Type, Measurement_type) %>%
  arrange(Cell_Type, Measurement_type)

for (i in 1:nrow(combos)) {
  ct <- combos$Cell_Type[i]
  mt <- combos$Measurement_type[i]
  tryCatch(
    run_analysis(combined_data, ct, mt),
    error = function(e) {
      # silently skip, or message(conditionMessage(e))
    }
  )
}