library(ggplot2)
library(dplyr)
library(tidyr)

theme_set(theme_bw())

vaccine_levels <- c("TM", "TMd21", "SOL", "IC")
formulation_cols <- c(
  TM    = "#2F75B5",
  TMd21 = "#D99A2B",
  SOL   = "#4A9E7D",
  IC    = "#B86691"
)

# ------------------------------------------------------------------
# Read the data — note sep = "\t"
# ------------------------------------------------------------------
df <- read.table(text = "
Vaccine	Day	Value	Measurement_type
IC	8	20.35	CD44+ CD62L+ TCM (% tetramer)
IC	14	33.08333333	CD44+ CD62L+ TCM (% tetramer)
IC	28	42.63333333	CD44+ CD62L+ TCM (% tetramer)
IC	60	46.66666667	CD44+ CD62L+ TCM (% tetramer)
SOL	8	12.7	CD44+ CD62L+ TCM (% tetramer)
SOL	14	31.71666667	CD44+ CD62L+ TCM (% tetramer)
SOL	28	39.76666667	CD44+ CD62L+ TCM (% tetramer)
SOL	60	41.35833333	CD44+ CD62L+ TCM (% tetramer)
TM	8	10.52	CD44+ CD62L+ TCM (% tetramer)
TM	14	16.045	CD44+ CD62L+ TCM (% tetramer)
TM	28	29.18333333	CD44+ CD62L+ TCM (% tetramer)
TM	60	33.99166667	CD44+ CD62L+ TCM (% tetramer)
TMd21	8	14.54	CD44+ CD62L+ TCM (% tetramer)
TMd21	14	18.79583333	CD44+ CD62L+ TCM (% tetramer)
TMd21	28	27.01666667	CD44+ CD62L+ TCM (% tetramer)
TMd21	60	32.11666667	CD44+ CD62L+ TCM (% tetramer)
IC	8	79.16666667	CD44+ CD62L- TEM (% tetramer)
IC	14	67.05833333	CD44+ CD62L- TEM (% tetramer)
IC	28	57.36666667	CD44+ CD62L- TEM (% tetramer)
IC	60	52.65833333	CD44+ CD62L- TEM (% tetramer)
SOL	8	86.61666667	CD44+ CD62L- TEM (% tetramer)
SOL	14	68.29166667	CD44+ CD62L- TEM (% tetramer)
SOL	28	60.23333333	CD44+ CD62L- TEM (% tetramer)
SOL	60	58.08333333	CD44+ CD62L- TEM (% tetramer)
TM	8	88.86666667	CD44+ CD62L- TEM (% tetramer)
TM	14	84.00833333	CD44+ CD62L- TEM (% tetramer)
TM	28	70.66666667	CD44+ CD62L- TEM (% tetramer)
TM	60	65.00833333	CD44+ CD62L- TEM (% tetramer)
TMd21	8	84.75	CD44+ CD62L- TEM (% tetramer)
TMd21	14	81.28333333	CD44+ CD62L- TEM (% tetramer)
TMd21	28	72.8	CD44+ CD62L- TEM (% tetramer)
TMd21	60	67.68333333	CD44+ CD62L- TEM (% tetramer)
", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# Reshape: one row per Vaccine-Day with TCM and TEM columns
# ------------------------------------------------------------------
df_wide <- df %>%
  mutate(Measurement = ifelse(grepl("TCM", Measurement_type), "TCM", "TEM")) %>%
  select(Vaccine, Day, Measurement, Value) %>%
  pivot_wider(names_from = Measurement, values_from = Value) %>%
  mutate(
    Ratio   = TEM / TCM,
    Vaccine = factor(Vaccine, levels = vaccine_levels)
  )

# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------
p <- ggplot(df_wide, aes(x = Day, y = Ratio, colour = Vaccine, group = Vaccine)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = formulation_cols, drop = FALSE) +
  labs(
    x      = "Day",
    y      = "Ratio (CD62L– TEM / CD62L+ TCM)",
    title  = "Conversion of CD62L– to CD62L+ over time",
    colour = "Vaccine"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p)

# Optional: log10 scale on y-axis
# p + scale_y_log10() + labs(y = "log10 (Ratio TEM / TCM)")