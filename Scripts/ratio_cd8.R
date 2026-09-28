library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------------
# 1. Build ratio data (day 8 and day 14, CD8, number unit)
#    NA Mouse_number rows are kept but made unique per row.
# ------------------------------------------------------------------
ratio_raw <- combined_data %>%
  filter(
    Cell_Type == "CD8",
    Unit == "number",
    Day %in% c(8, 14),
    Measurement_type %in% c(
      "Number of CD127+ KLRG1- memory",
      "Number of KLRG1+ CD127- effector"
    )
  ) %>%
  mutate(
    Cell = ifelse(grepl("memory", Measurement_type),
                  "Memory", "Effector"),
    Mouse_number = ifelse(is.na(Mouse_number),
                          paste0("NA_", row_number()),
                          as.character(Mouse_number))
  ) %>%
  select(Mouse_number, Vaccine, Day, Cell, Value) %>%
  pivot_wider(names_from = Cell, values_from = Value) %>%
  mutate(
    Ratio       = Memory / Effector,
    log10_Ratio = log10(Ratio)
  ) %>%
  filter(is.finite(log10_Ratio)) %>%
  filter(!Mouse_number %in% subject_to_exclude)

# Sanity checks
ratio_raw %>% count(Vaccine, Day) %>% print()
ratio_raw %>% count(Mouse_number, Vaccine, Day) %>% filter(n > 1) %>% print()

# ------------------------------------------------------------------
# 2. AIC partition search on the ratio slopes (15 partitions)
# ------------------------------------------------------------------
fit_ratio_slope_partition <- function(dat, cl, vaccines) {
  dat$Vaccine <- factor(dat$Vaccine, levels = vaccines)
  slope_cluster <- factor(cl[match(dat$Vaccine, vaccines)])
  for (g in levels(slope_cluster)) {
    dat[[paste0("Dg", g)]] <- ifelse(slope_cluster == g, dat$Day, 0)
  }
  fml <- as.formula(paste("log10_Ratio ~ Vaccine +",
                          paste(paste0("Dg", levels(slope_cluster)),
                                collapse = " + ")))
  mod <- tryCatch(lm(fml, data = dat), error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  list(aic = AIC(mod), k = length(coef(mod)), cl = cl, model = mod)
}

ratio_aic_sweep <- function(dat) {
  vaccines <- levels(droplevels(factor(dat$Vaccine)))
  parts <- generate_partitions(length(vaccines))       # Bell(4) = 15
  fits  <- lapply(parts, function(p) fit_ratio_slope_partition(dat, p, vaccines))
  keep  <- !vapply(fits, is.null, logical(1))
  fits  <- fits[keep]
  aics  <- vapply(fits, function(x) x$aic, numeric(1))
  ks    <- vapply(fits, function(x) x$k,   numeric(1))
  ord   <- order(aics)
  fits  <- fits[ord]; aics <- aics[ord]; ks <- ks[ord]
  
  best_aic <- aics[1]
  delta    <- aics - best_aic
  eligible <- which(delta < 10)
  min_k    <- min(ks[eligible])
  cand     <- eligible[ks[eligible] == min_k]
  parsim_i <- cand[which.min(aics[cand])]
  
  tbl <- data.frame(
    Groups       = vapply(fits, function(x) format_group(x$cl, vaccines), character(1)),
    n_slope_grps = vapply(fits, function(x) length(unique(x$cl)), integer(1)),
    k            = ks,
    AIC          = aics,
    Delta_AIC    = delta,
    Is_Best      = seq_along(fits) == 1,
    Is_Parsim    = seq_along(fits) == parsim_i
  )
  list(table = tbl, best = fits[[1]], parsimonious = fits[[parsim_i]],
       vaccines = vaccines)
}

ratio_sweep <- ratio_aic_sweep(ratio_raw)
ratio_sweep$table

# ------------------------------------------------------------------
# 3. Per-vaccine slope, half-life, CI for log10(M/E)
# ------------------------------------------------------------------
ratio_slopes_ci <- function(dat, cl, vaccines, model_label) {
  dat$Vaccine <- factor(dat$Vaccine, levels = vaccines)
  slope_cluster <- factor(cl[match(dat$Vaccine, vaccines)])
  for (g in levels(slope_cluster)) {
    dat[[paste0("Dg", g)]] <- ifelse(slope_cluster == g, dat$Day, 0)
  }
  fml <- as.formula(paste("log10_Ratio ~ Vaccine +",
                          paste(paste0("Dg", levels(slope_cluster)),
                                collapse = " + ")))
  mod <- lm(fml, data = dat)
  cf  <- coef(mod); vc <- vcov(mod)
  
  do.call(rbind, lapply(seq_along(vaccines), function(i) {
    v  <- vaccines[i]
    g  <- as.character(cl[i])
    nm <- paste0("Dg", g)
    est <- cf[[nm]]
    se  <- sqrt(vc[nm, nm])
    lo  <- est - 1.96 * se
    hi  <- est + 1.96 * se
    data.frame(
      Model          = model_label,
      Vaccine        = v,
      Slope          = est,
      Slope_Lower    = lo,
      Slope_Upper    = hi,
      # Half-life of the ratio (doubling if slope>0, halving if slope<0)
      Ratio_HalfLife = -log10(2) / est,
      Ratio_HL_Lower = if (lo < 0) -log10(2) / lo else NA_real_,
      Ratio_HL_Upper = if (hi < 0) -log10(2) / hi else NA_real_
    )
  }))
}

all_diff_cl <- seq_along(ratio_sweep$vaccines)

ratio_slopes <- rbind(
  ratio_slopes_ci(ratio_raw, ratio_sweep$best$cl, ratio_sweep$vaccines,
                  "Best AIC partition"),
  ratio_slopes_ci(ratio_raw, all_diff_cl,         ratio_sweep$vaccines,
                  "All 4 different")
)
ratio_slopes

# ------------------------------------------------------------------
# 4. Plots
# ------------------------------------------------------------------
ggplot(ratio_raw, aes(x = Day, y = log10_Ratio, color = Vaccine)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
  scale_color_manual(values = formulation_cols, drop = FALSE) +
  labs(title = "log10(Memory / Effector), day 8 vs 14",
       y = "log10(Memory/Effector)", x = "Day") +
  big_text_theme +
  theme(legend.position = "right")

ggplot(ratio_slopes %>% filter(Model == "Best AIC partition"),
       aes(x = Vaccine, y = Slope, color = Vaccine)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Slope_Lower, ymax = Slope_Upper), width = 0.2) +
  scale_color_manual(values = formulation_cols, drop = FALSE) +
  labs(title = "Slope of log10(M/E) per day\n(best-AIC partition)",
       y = "Slope (log10 ratio per day)", x = NULL) +
  big_text_theme

# ------------------------------------------------------------------
# 5. Direct interaction test (independent of AIC partition)
# ------------------------------------------------------------------
m_int <- lm(log10_Ratio ~ Vaccine * Day, data = ratio_raw)
anova(m_int)



# ------------------------------------------------------------------
# 3. Per-vaccine slope + CI: best-AIC and parsimonious
# ------------------------------------------------------------------
ratio_slopes_ci <- function(dat, cl, vaccines, model_label) {
  dat$Vaccine <- factor(dat$Vaccine, levels = vaccines)
  slope_cluster <- factor(cl[match(dat$Vaccine, vaccines)])
  for (g in levels(slope_cluster)) {
    dat[[paste0("Dg", g)]] <- ifelse(slope_cluster == g, dat$Day, 0)
  }
  fml <- as.formula(paste("log10_Ratio ~ Vaccine +",
                          paste(paste0("Dg", levels(slope_cluster)),
                                collapse = " + ")))
  mod <- lm(fml, data = dat)
  cf  <- coef(mod); vc <- vcov(mod)
  
  do.call(rbind, lapply(seq_along(vaccines), function(i) {
    v  <- vaccines[i]
    g  <- as.character(cl[i])
    nm <- paste0("Dg", g)
    est <- cf[[nm]]
    se  <- sqrt(vc[nm, nm])
    lo  <- est - 1.96 * se
    hi  <- est + 1.96 * se
    data.frame(
      Model          = model_label,
      Vaccine        = v,
      Slope          = est,
      Slope_Lower    = lo,
      Slope_Upper    = hi,
      Ratio_HalfLife = -log10(2) / est,
      Ratio_HL_Lower = if (lo < 0) -log10(2) / lo else NA_real_,
      Ratio_HL_Upper = if (hi < 0) -log10(2) / hi else NA_real_
    )
  }))
}

all_diff_cl <- seq_along(ratio_sweep$vaccines)

ratio_slopes <- rbind(
  ratio_slopes_ci(ratio_raw, ratio_sweep$best$cl,         ratio_sweep$vaccines,
                  "Best AIC partition"),
  ratio_slopes_ci(ratio_raw, ratio_sweep$parsimonious$cl, ratio_sweep$vaccines,
                  "Parsimonious partition"),
  ratio_slopes_ci(ratio_raw, all_diff_cl,                 ratio_sweep$vaccines,
                  "All 4 different")
)
ratio_slopes

# ------------------------------------------------------------------
# 4. Combined faceted plot: Best AIC vs Parsimonious
# ------------------------------------------------------------------
plot_labels <- ratio_sweep$table %>%
  filter(Is_Best | Is_Parsim) %>%
  mutate(
    Model = ifelse(Is_Best, "Best AIC partition", "Parsimonious partition"),
    label = paste0(Model, "\n", Groups)
  ) %>%
  select(Model, label) %>%
  distinct()

plot_df <- ratio_slopes %>%
  filter(Model %in% c("Best AIC partition", "Parsimonious partition")) %>%
  left_join(plot_labels, by = "Model") %>%
  mutate(Model = factor(Model,
                        levels = c("Best AIC partition",
                                   "Parsimonious partition")))

plot_both_slopes <- ggplot(
  plot_df,
  aes(x = Vaccine, y = Slope, color = Vaccine)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Slope_Lower, ymax = Slope_Upper), width = 0.2) +
  facet_wrap(~ label, ncol = 2) +
  scale_color_manual(values = formulation_cols, drop = FALSE) +
  labs(
    title = "Slope of log10(M/E) per day",
    y = "Slope (log10 ratio per day)",
    x = NULL
  ) +
  big_text_theme

plot_both_slopes

# ------------------------------------------------------------------
# 5. Same, but as ratio doubling time (days)
# ------------------------------------------------------------------
plot_both_halflife <- ggplot(
  plot_df,
  aes(x = Vaccine, y = Ratio_HalfLife, color = Vaccine)
) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Ratio_HL_Lower, ymax = Ratio_HL_Upper), width = 0.2) +
  facet_wrap(~ label, ncol = 2) +
  scale_color_manual(values = formulation_cols, drop = FALSE) +
  labs(
    title = "Doubling time of M/E ratio (days)",
    y = "Doubling time (days)",
    x = NULL
  ) +
  big_text_theme

plot_both_halflife

# ------------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------------
ggsave(file.path(figure_folder, "CD8_analysis",
                 "Number_of_CD127_KLRG1_memory_vs_effector",
                 "MEmodel_bestAIC_vs_parsimonious_slope.png"),
       plot_both_slopes, width = 7.5, height = 3.6, dpi = 300)

ggsave(file.path(figure_folder, "CD8_analysis",
                 "Number_of_CD127_KLRG1_memory_vs_effector",
                 "MEmodel_bestAIC_vs_parsimonious_doubling.png"),
       plot_both_halflife, width = 7.5, height = 3.6, dpi = 300)