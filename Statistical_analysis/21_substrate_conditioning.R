# =============================================================================
# Script 21: Substrate conditioning — further evidence
# Two analyses:
#   A) Within-group dose-response: does mScarlet expression (proxy for
#      viral transduction) predict enwrapment reduction within each group?
#      If C6ST1 conditions the substrate, C6ST1_ADAMTS15 animals with higher
#      expression should show stronger effects than equivalent ADAMTS15 animals.
#   B) Spatial gradient comparison: does C6ST1_ADAMTS15 show a broader
#      spatial footprint than ADAMTS15? Compare slopes, recovery distances,
#      and area under the enwrapment-vs-distance curve.
# Inputs:
#   14_scarlet_expression/scarlet_expression_animal.csv
#   15_spatial_gradient/animal_gradient.csv
#   15_spatial_gradient/cell_level_distances.csv
#   19_ipsi_contra_ratio/animal_ratios.csv
# Output: results/21_substrate_conditioning/
# =============================================================================

library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------
BASE    <- "/path/to/results"
OUT_DIR <- file.path(BASE, "21_substrate_conditioning")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TREATMENT_ORDER <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")

PALETTE_2 <- c(ADAMTS15 = "#8e44ad", C6ST1_ADAMTS15 = "#27ae60")
PALETTE   <- c(mScarlet       = "#c0392b", ADAMTS4        = "#2980b9",
               ADAMTS4_MD     = "#1abc9c", ADAMTS15       = "#8e44ad",
               C6ST1          = "#e67e22", C6ST1_ADAMTS15 = "#27ae60")

# =============================================================================
# PART A: Within-group dose-response
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part A: Within-group dose-response\n")
cat("══════════════════════════════════════════════════════\n")

# Load expression data
expr <- fread(file.path(BASE, "14_scarlet_expression", "scarlet_expression_animal.csv"))

# Load enwrapment outcome (log ipsi/contra ratio from script 19)
ratios <- fread(file.path(BASE, "19_ipsi_contra_ratio", "animal_ratios.csv"))
ratios <- ratios[, .(mouse_id, treatment, log_ratio, ratio)]

# Join — match on mouse_id

expr <- expr[mouse_id != ""]   # Exlude animal if needed (injection failure)

merged <- merge(ratios, expr, by = "mouse_id", suffixes = c("", "_expr"))
merged[, treatment := treatment]   # use corrected treatment from ratios

cat(sprintf("Merged: %d animals\n", nrow(merged)))

# Restrict to treated groups (not mScarlet — no protease to show dose-response)
treated <- merged[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15",
                                    "ADAMTS4", "ADAMTS4_MD", "C6ST1")]

# Correlate log_ratio (more negative = stronger effect) with Core expression
# Use n_scarlet_core as the expression proxy (cell count)
doseresponse <- treated[, {
  if (.N >= 3) {
    ct <- cor.test(n_scarlet_core, log_ratio, method = "pearson")
    .(r       = ct$estimate,
      p       = ct$p.value,
      n       = .N,
      direction = ifelse(ct$estimate < 0,
                         "more expression → stronger effect",
                         "more expression → weaker effect"))
  } else {
    .(r = NA_real_, p = NA_real_, n = .N, direction = NA_character_)
  }
}, by = treatment]

cat("\nWithin-group correlation: n_scarlet_core vs log_ratio\n")
cat("(negative r = more expression → lower ratio = stronger PNN reduction)\n")
print(doseresponse)
fwrite(doseresponse, file.path(OUT_DIR, "within_group_doseresponse.csv"))

# Also save the merged table for plotting
fwrite(merged, file.path(OUT_DIR, "expression_enwrapment_merged.csv"))

# Fig A: scatter plots for ADAMTS15 and C6ST1_ADAMTS15
two_tx <- treated[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15")]
two_tx[, treatment := factor(treatment, levels = c("ADAMTS15", "C6ST1_ADAMTS15"))]

pA <- ggplot(two_tx, aes(x = n_scarlet_core, y = log_ratio, colour = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(label = mouse_id), vjust = -0.8, size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = PALETTE_2) +
  facet_wrap(~ treatment, scales = "free_x") +
  labs(
    title    = "Within-group dose-response: expression vs enwrapment reduction",
    subtitle = "More negative log ratio = stronger PNN reduction",
    x        = "mScarlet+ cell count in Core zone (expression proxy)",
    y        = "log(ipsi/contra enwrapment ratio)",
    colour   = NULL,
    caption  = "Pearson r; shaded = 95% CI"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "fig_doseresponse.pdf"), pA, width = 8, height = 4.5)

# =============================================================================
# PART B: Spatial gradient comparison
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part B: Spatial gradient comparison\n")
cat("══════════════════════════════════════════════════════\n")

# Load per-animal slopes
slopes <- fread(file.path(BASE, "15_spatial_gradient", "animal_gradient.csv"))

# Apply corrections
slopes[, treatment := gsub("C6ST1_ADAMTS4", "C6ST1_ADAMTS15", treatment)]
slopes <- slopes[mouse_id != "C6ST1_ADAMTS4_4"]

cat("\nPer-animal slopes (beta_dist = enwrapment increase per mm from centre):\n")
print(slopes[order(treatment), .(mouse_id, treatment,
                                  slope = round(beta_dist, 4),
                                  r     = round(r, 3),
                                  p     = round(p, 4),
                                  n_cells)])

# Compare slopes between ADAMTS15 and C6ST1_ADAMTS15
two_slopes <- slopes[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15")]
slope_test <- t.test(beta_dist ~ treatment, data = two_slopes,
                     var.equal = FALSE)
cat(sprintf("\nSlope comparison (ADAMTS15 vs C6ST1_ADAMTS15):\n"))
cat(sprintf("  ADAMTS15       mean slope = %.4f\n",
            mean(two_slopes[treatment == "ADAMTS15", beta_dist])))
cat(sprintf("  C6ST1_ADAMTS15 mean slope = %.4f\n",
            mean(two_slopes[treatment == "C6ST1_ADAMTS15", beta_dist])))
cat(sprintf("  t = %.3f, df = %.1f, p = %.4f\n",
            slope_test$statistic, slope_test$parameter, slope_test$p.value))

slope_summary <- data.table(
  treatment   = c("ADAMTS15", "C6ST1_ADAMTS15"),
  mean_slope  = c(mean(two_slopes[treatment == "ADAMTS15",       beta_dist]),
                  mean(two_slopes[treatment == "C6ST1_ADAMTS15", beta_dist])),
  sd_slope    = c(sd(two_slopes[treatment == "ADAMTS15",       beta_dist]),
                  sd(two_slopes[treatment == "C6ST1_ADAMTS15", beta_dist])),
  n           = c(nrow(two_slopes[treatment == "ADAMTS15"]),
                  nrow(two_slopes[treatment == "C6ST1_ADAMTS15"])),
  t_stat      = slope_test$statistic,
  p_slope     = slope_test$p.value
)
fwrite(slope_summary, file.path(OUT_DIR, "slope_comparison.csv"))

# Load cell-level distances for AUC and recovery analysis
cat("\nLoading cell_level_distances.csv ...\n")
cells <- fread(file.path(BASE, "15_spatial_gradient", "cell_level_distances.csv"))

# Apply corrections
cells[, treatment := gsub("C6ST1_ADAMTS4", "C6ST1_ADAMTS15", treatment)]
cells <- cells[mouse_id != "C6ST1_ADAMTS4_4"]
cells[mouse_id == "mScarlet_4",   treatment := "ADAMTS4_MD"]
cells[mouse_id == "ADAMTS4_MD_4", treatment := "mScarlet"]

cat(sprintf("Cell-level rows: %d\n", nrow(cells)))

# Bin distances and compute mean enwrapment per bin per treatment
cells[, dist_bin := round(dist_mm / 0.1) * 0.1]   # 0.1 mm bins
binned <- cells[, .(
  frac_enwrapped = mean(enwrapped, na.rm = TRUE),
  n_cells        = .N
), by = .(treatment, dist_bin)]

# Only keep bins with >= 20 cells
binned <- binned[n_cells >= 20]
binned[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

# Compute mScarlet baseline at each bin
mscarlet_baseline <- binned[treatment == "mScarlet",
                             .(dist_bin, baseline = frac_enwrapped)]

# Join baseline and compute recovery fraction
binned <- merge(binned, mscarlet_baseline, by = "dist_bin", all.x = TRUE)
binned[, pct_of_baseline := frac_enwrapped / baseline]

# Find recovery distance: first bin where each treatment reaches >= 90% baseline
recovery <- binned[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15") &
                   !is.na(pct_of_baseline)][
  order(treatment, dist_bin)][
  pct_of_baseline >= 0.90, .SD[1], by = treatment]

cat("\nRecovery to 90% of mScarlet baseline:\n")
if (nrow(recovery) > 0) {
  print(recovery[, .(treatment, dist_bin, frac_enwrapped,
                     pct_of_baseline = round(pct_of_baseline, 3))])
} else {
  cat("  Neither treatment fully recovers to 90% baseline within observed range\n")
}
fwrite(recovery, file.path(OUT_DIR, "recovery_distance.csv"))

# AUC comparison: integrate enwrapment deficit vs mScarlet over 0–1 mm
auc_dat <- binned[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15", "mScarlet") &
                  dist_bin <= 1.0 & !is.na(baseline)]

auc_summary <- auc_dat[treatment != "mScarlet", .(
  auc_deficit = sum((baseline - frac_enwrapped) * 0.1, na.rm = TRUE),
  n_bins      = .N
), by = treatment]

cat("\nAUC of enwrapment deficit vs mScarlet (0–1 mm):\n")
cat("(larger = greater total PNN reduction across the spatial extent)\n")
print(auc_summary)
fwrite(auc_summary, file.path(OUT_DIR, "auc_deficit.csv"))

fwrite(binned, file.path(OUT_DIR, "binned_gradient_with_baseline.csv"))

# Fig B1: Spatial gradient curves ADAMTS15 vs C6ST1_ADAMTS15 vs mScarlet
plot_tx <- c("mScarlet", "ADAMTS15", "C6ST1_ADAMTS15")
plot_dat <- binned[treatment %in% plot_tx & dist_bin <= 1.5]
plot_dat[, treatment := factor(treatment, levels = plot_tx)]

pB1 <- ggplot(plot_dat, aes(x = dist_bin, y = frac_enwrapped,
                              colour = treatment)) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.15,
              linewidth = 0.9) +
  scale_colour_manual(values = c(mScarlet       = "#c0392b",
                                  ADAMTS15       = "#8e44ad",
                                  C6ST1_ADAMTS15 = "#27ae60")) +
  labs(
    title   = "Spatial enwrapment profile: 0–1.5 mm from injection centre",
    x       = "Distance from injection centre (mm)",
    y       = "Fraction enwrapped",
    colour  = NULL,
    caption = "LOESS smooth; shaded = 95% CI"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_spatial_profiles.pdf"), pB1, width = 7, height = 5)

# Fig B2: Per-animal slope dotplot
slopes_plot <- slopes[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15")]
slopes_plot[, treatment := factor(treatment,
                                   levels = c("ADAMTS15", "C6ST1_ADAMTS15"))]

pB2 <- ggplot(slopes_plot, aes(x = treatment, y = beta_dist,
                                 colour = treatment)) +
  geom_hline(yintercept = mean(slopes[treatment == "mScarlet", beta_dist]),
             linetype = "dashed", colour = "#c0392b", linewidth = 0.8) +
  annotate("text", x = 0.6, y = mean(slopes[treatment == "mScarlet", beta_dist]),
           label = "mScarlet mean", colour = "#c0392b", hjust = 0,
           vjust = -0.5, size = 3.2) +
  geom_jitter(width = 0.08, size = 4, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.3, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE_2) +
  labs(
    title   = "Per-animal spatial gradient slopes",
    subtitle= sprintf("t-test p = %.4f", slope_test$p.value),
    x       = NULL,
    y       = "Slope (frac enwrapped per mm from injection centre)",
    colour  = NULL,
    caption = "Steeper slope = faster recovery with distance = more localised effect"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "fig_slope_comparison.pdf"), pB2, width = 5, height = 5)

# =============================================================================
# Summary
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 21 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  within_group_doseresponse.csv\n")
cat("  expression_enwrapment_merged.csv\n")
cat("  slope_comparison.csv\n")
cat("  recovery_distance.csv\n")
cat("  auc_deficit.csv\n")
cat("  binned_gradient_with_baseline.csv\n")
cat("  fig_doseresponse.pdf\n")
cat("  fig_spatial_profiles.pdf\n")
cat("  fig_slope_comparison.pdf\n")
cat("══════════════════════════════════════════════════════\n")
