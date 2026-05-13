# =============================================================================
# Script 22: ADAMTS4 vs ADAMTS4_MD — monocular deprivation effect
# Mirrors script 21 structure but tests RQ4:
#   Does MD shift the spatial profile or dose-response of ADAMTS4?
# Three analyses:
#   A) Direct pairwise contrast in ipsi/contra ratio model
#   B) Within-group dose-response: expression vs enwrapment per group
#   C) Spatial gradient comparison: slopes, AUC, recovery distance
#
# Expected outcome: genuine null — confirms RQ4 rigorously

# Output: results/22_adamts4_md_effect/
# =============================================================================

library(data.table)
library(emmeans)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------
BASE    <- "/path/to/results"
OUT_DIR <- file.path(BASE, "22_adamts4_md_effect")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TREATMENT_ORDER <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")

PALETTE_2 <- c(ADAMTS4 = "#2980b9", ADAMTS4_MD = "#1abc9c")

# =============================================================================
# PART A: Direct pairwise contrast in ratio model
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part A: Pairwise ratio contrast — ADAMTS4 vs ADAMTS4_MD\n")
cat("══════════════════════════════════════════════════════\n")

animal_dat <- fread(file.path(BASE, "19_ipsi_contra_ratio", "animal_ratios.csv"))
animal_dat[, treatment := factor(treatment, levels = TREATMENT_ORDER)]


cat("\nADAMTS4 and ADAMTS4_MD animals:\n")
print(animal_dat[treatment %in% c("ADAMTS4", "ADAMTS4_MD"),
                 .(mouse_id, treatment, ratio = round(ratio, 3),
                   log_ratio = round(log_ratio, 3), is_swapped)])

# Full model (all animals)
m_log  <- lm(log_ratio ~ treatment, data = animal_dat)
em_log <- emmeans(m_log, ~ treatment)
ct_log <- as.data.table(contrast(em_log, method = "pairwise", adjust = "none"))
ct_log[, sig := ifelse(p.value < 0.05, "*", "")]

key <- ct_log[contrast %in% c("ADAMTS4 - ADAMTS4_MD", "ADAMTS4_MD - ADAMTS4")]
cat("\nDirect contrast: ADAMTS4 vs ADAMTS4_MD (log ratio, all animals)\n")
print(key[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

# Sensitivity: exclude swapped animal
animal_nosw <- animal_dat[(!is_swapped)]
m_nosw  <- lm(log_ratio ~ treatment, data = animal_nosw)
em_nosw <- emmeans(m_nosw, ~ treatment)
ct_nosw <- as.data.table(contrast(em_nosw, method = "pairwise", adjust = "none"))
ct_nosw[, sig := ifelse(p.value < 0.05, "*", "")]

key_nosw <- ct_nosw[contrast %in% c("ADAMTS4 - ADAMTS4_MD", "ADAMTS4_MD - ADAMTS4")]
cat("\nDirect contrast: ADAMTS4 vs ADAMTS4_MD (log ratio, swapped animal excluded)\n")
print(key_nosw[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

# Group means both ways
gmeans_all <- animal_dat[treatment %in% c("ADAMTS4", "ADAMTS4_MD"), .(
  mean_ratio     = mean(ratio),
  sd_ratio       = sd(ratio),
  mean_log_ratio = mean(log_ratio),
  n              = .N
), by = treatment]

gmeans_nosw <- animal_nosw[treatment %in% c("ADAMTS4", "ADAMTS4_MD"), .(
  mean_ratio     = mean(ratio),
  sd_ratio       = sd(ratio),
  mean_log_ratio = mean(log_ratio),
  n              = .N
), by = treatment]

cat("\nGroup means (all animals):\n");  print(gmeans_all)
cat("\nGroup means (swapped animal excluded):\n"); print(gmeans_nosw)

fwrite(ct_log,      file.path(OUT_DIR, "pairwise_log_contrasts.csv"))
fwrite(ct_nosw,     file.path(OUT_DIR, "pairwise_log_contrasts_noswaped.csv"))
fwrite(gmeans_all,  file.path(OUT_DIR, "group_means_all.csv"))
fwrite(gmeans_nosw, file.path(OUT_DIR, "group_means_noswapped.csv"))

# =============================================================================
# PART B: Within-group dose-response
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part B: Within-group dose-response\n")
cat("══════════════════════════════════════════════════════\n")

expr <- fread(file.path(BASE, "14_scarlet_expression", "scarlet_expression_animal.csv"))

# Remap IDs to match ratios file
ratios <- fread(file.path(BASE, "19_ipsi_contra_ratio", "animal_ratios.csv"))
ratios <- ratios[, .(mouse_id, treatment, log_ratio, ratio)]
ratios[, mouse_id_for_join := gsub("C6ST1_ADAMTS4_", "C6ST1_ADAMTS15_", mouse_id)]

merged <- merge(ratios,
                expr[, .(mouse_id, n_scarlet_core, mean_scarlet_core)],
                by.x = "mouse_id_for_join", by.y = "mouse_id")

cat("\nADAMTS4 dose-response:\n")
d4 <- merged[treatment == "ADAMTS4"]
print(d4[, .(mouse_id, log_ratio = round(log_ratio,3),
             n_scarlet_core, mean_scarlet_core = round(mean_scarlet_core,0))])
ct4 <- cor.test(d4$n_scarlet_core, d4$log_ratio)
cat(sprintf("r = %.3f, p = %.3f, n = %d — %s\n",
            ct4$estimate, ct4$p.value, nrow(d4),
            ifelse(ct4$estimate < 0,
                   "more expression -> stronger effect",
                   "more expression -> weaker effect")))

cat("\nADAMTS4_MD dose-response (all animals including swapped):\n")
dmd <- merged[treatment == "ADAMTS4_MD"]
print(dmd[, .(mouse_id, log_ratio = round(log_ratio,3),
              n_scarlet_core, mean_scarlet_core = round(mean_scarlet_core,0),
              swapped = mouse_id %in% c("mScarlet_4","ADAMTS4_MD_4"))])
ctmd <- cor.test(dmd$n_scarlet_core, dmd$log_ratio)
cat(sprintf("r = %.3f, p = %.3f, n = %d — %s\n",
            ctmd$estimate, ctmd$p.value, nrow(dmd),
            ifelse(ctmd$estimate < 0,
                   "more expression -> stronger effect",
                   "more expression -> weaker effect")))

cat("\nADAMTS4_MD dose-response (swapped animal excluded):\n")
dmd_nosw <- dmd[!mouse_id %in% c("mScarlet_4")]
if (nrow(dmd_nosw) >= 3) {
  ctmd_nosw <- cor.test(dmd_nosw$n_scarlet_core, dmd_nosw$log_ratio)
  cat(sprintf("r = %.3f, p = %.3f, n = %d — %s\n",
              ctmd_nosw$estimate, ctmd_nosw$p.value, nrow(dmd_nosw),
              ifelse(ctmd_nosw$estimate < 0,
                     "more expression -> stronger effect",
                     "more expression -> weaker effect")))
} else {
  cat("Insufficient animals after exclusion\n")
}

doseresponse <- data.table(
  treatment = c("ADAMTS4", "ADAMTS4_MD (all)", "ADAMTS4_MD (no swap)"),
  r         = c(ct4$estimate,
                ctmd$estimate,
                ifelse(nrow(dmd_nosw) >= 3, ctmd_nosw$estimate, NA)),
  p         = c(ct4$p.value,
                ctmd$p.value,
                ifelse(nrow(dmd_nosw) >= 3, ctmd_nosw$p.value, NA)),
  n         = c(nrow(d4), nrow(dmd), nrow(dmd_nosw))
)
doseresponse[, direction := ifelse(r < 0,
                                    "more expression -> stronger effect",
                                    "more expression -> weaker effect")]
fwrite(doseresponse, file.path(OUT_DIR, "within_group_doseresponse.csv"))

# =============================================================================
# PART C: Spatial gradient comparison
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part C: Spatial gradient comparison\n")
cat("══════════════════════════════════════════════════════\n")

slopes <- fread(file.path(BASE, "15_spatial_gradient", "animal_gradient.csv"))
slopes[, treatment := gsub("C6ST1_ADAMTS4", "C6ST1_ADAMTS15", treatment)]
slopes <- slopes[mouse_id != "C6ST1_ADAMTS4_4"]

cat("\nPer-animal slopes — ADAMTS4 and ADAMTS4_MD:\n")
print(slopes[treatment %in% c("ADAMTS4", "ADAMTS4_MD"),
             .(mouse_id, treatment,
               slope = round(beta_dist, 4),
               r     = round(r, 3),
               p     = round(p, 4),
               n_cells,
               swapped = mouse_id == "mScarlet_4")])

# Slope t-test: all animals
two_slopes <- slopes[treatment %in% c("ADAMTS4", "ADAMTS4_MD")]
slope_test <- t.test(beta_dist ~ treatment, data = two_slopes, var.equal = FALSE)
cat(sprintf("\nSlope comparison (all): t=%.3f, p=%.4f\n",
            slope_test$statistic, slope_test$p.value))

# Slope t-test: exclude swapped
two_slopes_nosw <- two_slopes[mouse_id != "mScarlet_4"]
if (nrow(two_slopes_nosw[treatment=="ADAMTS4_MD"]) >= 2) {
  slope_nosw <- t.test(beta_dist ~ treatment, data = two_slopes_nosw, var.equal = FALSE)
  cat(sprintf("Slope comparison (no swap): t=%.3f, p=%.4f\n",
              slope_nosw$statistic, slope_nosw$p.value))
}

slope_summary <- two_slopes[, .(
  mean_slope = mean(beta_dist),
  sd_slope   = sd(beta_dist),
  n          = .N
), by = treatment]
cat("\nSlope summary:\n"); print(slope_summary)
fwrite(slope_summary, file.path(OUT_DIR, "slope_comparison.csv"))

# Load cell-level distances
cells <- fread(file.path(BASE, "15_spatial_gradient", "cell_level_distances.csv"))
cells[, treatment := gsub("C6ST1_ADAMTS4", "C6ST1_ADAMTS15", treatment)]
cells <- cells[mouse_id != "C6ST1_ADAMTS4_4"]
cells[mouse_id == "mScarlet_4",   treatment := "ADAMTS4_MD"]
cells[mouse_id == "ADAMTS4_MD_4", treatment := "mScarlet"]

# Bin distances
cells[, dist_bin := round(dist_mm / 0.1) * 0.1]
binned <- cells[, .(
  frac_enwrapped = mean(enwrapped, na.rm = TRUE),
  n_cells        = .N
), by = .(treatment, dist_bin)]
binned <- binned[n_cells >= 20]
binned[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

# mScarlet baseline
mscarlet_baseline <- binned[treatment == "mScarlet",
                             .(dist_bin, baseline = frac_enwrapped)]
binned <- merge(binned, mscarlet_baseline, by = "dist_bin", all.x = TRUE)
binned[, pct_of_baseline := frac_enwrapped / baseline]

# Recovery distance
recovery <- binned[treatment %in% c("ADAMTS4", "ADAMTS4_MD") &
                   !is.na(pct_of_baseline)][
  order(treatment, dist_bin)][
  pct_of_baseline >= 0.90, .SD[1], by = treatment]
cat("\nRecovery to 90% of mScarlet baseline:\n")
if (nrow(recovery) > 0) {
  print(recovery[, .(treatment, dist_bin,
                     pct_of_baseline = round(pct_of_baseline, 3))])
} else {
  cat("  Neither treatment recovers to 90% baseline within observed range\n")
}
fwrite(recovery, file.path(OUT_DIR, "recovery_distance.csv"))

# AUC deficit 0-1 mm
auc_dat <- binned[treatment %in% c("ADAMTS4", "ADAMTS4_MD", "mScarlet") &
                  dist_bin <= 1.0 & !is.na(baseline)]
auc_summary <- auc_dat[treatment != "mScarlet", .(
  auc_deficit = sum((baseline - frac_enwrapped) * 0.1, na.rm = TRUE),
  n_bins      = .N
), by = treatment]
cat("\nAUC of enwrapment deficit vs mScarlet (0-1 mm):\n")
print(auc_summary)
fwrite(auc_summary, file.path(OUT_DIR, "auc_deficit.csv"))
fwrite(binned,      file.path(OUT_DIR, "binned_gradient_with_baseline.csv"))

# =============================================================================
# Figures
# =============================================================================

# Fig A: Ratio dotplot with swapped animal flagged
plot_two <- animal_dat[treatment %in% c("ADAMTS4", "ADAMTS4_MD")]
plot_two[, treatment := factor(treatment, levels = c("ADAMTS4", "ADAMTS4_MD"))]

pA <- ggplot(plot_two, aes(x = treatment, y = ratio, colour = treatment)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_jitter(aes(shape = is_swapped), width = 0.08, size = 4, alpha = 0.9) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8),
                     labels = c("FALSE" = "Normal", "TRUE" = "Swapped animal")) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.3, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE_2) +
  scale_y_continuous(limits = c(0.3, 1.1)) +
  labs(
    title   = "Ipsi/contra ratio: ADAMTS4 vs ADAMTS4_MD",
    x       = NULL,
    y       = "Ipsi / contra enwrapment ratio",
    colour  = NULL,
    shape   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_ratio_dotplot.pdf"), pA, width = 5, height = 5)

# Fig B: Dose-response scatter
plot_dr <- merged[treatment %in% c("ADAMTS4", "ADAMTS4_MD")]
plot_dr[, treatment := factor(treatment, levels = c("ADAMTS4", "ADAMTS4_MD"))]
plot_dr[, is_swapped := mouse_id == "mScarlet_4"]

pB <- ggplot(plot_dr, aes(x = n_scarlet_core, y = log_ratio, colour = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  geom_point(aes(shape = is_swapped), size = 3.5) +
  geom_text(aes(label = mouse_id), vjust = -0.8, size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = PALETTE_2) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8),
                     labels = c("FALSE" = "Normal", "TRUE" = "Swapped animal")) +
  facet_wrap(~ treatment, scales = "free_x") +
  labs(
    title   = "Within-group dose-response: ADAMTS4 vs ADAMTS4_MD",
    x       = "mScarlet+ cell count in Core zone",
    y       = "log(ipsi/contra enwrapment ratio)",
    colour  = NULL, shape = NULL,
    caption = "* = swapped animal (correctly attributed to ADAMTS4_MD)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_doseresponse.pdf"), pB, width = 8, height = 4.5)

# Fig C: Spatial profiles
plot_tx  <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD")
plot_dat <- binned[treatment %in% plot_tx & dist_bin <= 1.5]
plot_dat[, dist_bin  := as.numeric(dist_bin)]
plot_dat[, treatment := factor(treatment, levels = plot_tx)]

pC <- ggplot(plot_dat, aes(x = dist_bin, y = frac_enwrapped, colour = treatment)) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.15,
              linewidth = 0.9) +
  scale_colour_manual(values = c(mScarlet   = "#c0392b",
                                  ADAMTS4    = "#2980b9",
                                  ADAMTS4_MD = "#1abc9c")) +
  labs(
    title   = "Spatial enwrapment profile: ADAMTS4 vs ADAMTS4_MD",
    x       = "Distance from injection centre (mm)",
    y       = "Fraction enwrapped",
    colour  = NULL,
    caption = "LOESS smooth; shaded = 95% CI"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_spatial_profiles.pdf"), pC, width = 7, height = 5)

# Fig D: Per-animal slope dotplot
slopes_plot <- slopes[treatment %in% c("ADAMTS4", "ADAMTS4_MD")]
slopes_plot[, treatment  := factor(treatment, levels = c("ADAMTS4", "ADAMTS4_MD"))]
slopes_plot[, is_swapped := mouse_id == "mScarlet_4"]

pD <- ggplot(slopes_plot, aes(x = treatment, y = beta_dist, colour = treatment)) +
  geom_hline(yintercept = mean(slopes[treatment == "mScarlet", beta_dist]),
             linetype = "dashed", colour = "#c0392b", linewidth = 0.8) +
  annotate("text", x = 0.6,
           y = mean(slopes[treatment == "mScarlet", beta_dist]),
           label = "mScarlet mean", colour = "#c0392b",
           hjust = 0, vjust = -0.5, size = 3.2) +
  geom_jitter(aes(shape = is_swapped), width = 0.08, size = 4, alpha = 0.9) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.3, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE_2) +
  labs(
    title   = "Per-animal spatial gradient slopes: ADAMTS4 vs ADAMTS4_MD",
    x       = NULL,
    y       = "Slope (frac enwrapped per mm)",
    colour  = NULL, shape  = NULL,
    caption = sprintf("t-test p = %.4f (all animals)", slope_test$p.value)
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_slope_comparison.pdf"), pD, width = 5, height = 5)

# =============================================================================
# Summary
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 22 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  pairwise_log_contrasts.csv\n")
cat("  pairwise_log_contrasts_noswaped.csv\n")
cat("  group_means_all.csv\n")
cat("  group_means_noswapped.csv\n")
cat("  within_group_doseresponse.csv\n")
cat("  slope_comparison.csv\n")
cat("  recovery_distance.csv\n")
cat("  auc_deficit.csv\n")
cat("  binned_gradient_with_baseline.csv\n")
cat("  fig_ratio_dotplot.pdf\n")
cat("  fig_doseresponse.pdf\n")
cat("  fig_spatial_profiles.pdf\n")
cat("  fig_slope_comparison.pdf\n")
cat("══════════════════════════════════════════════════════\n")
