# =============================================================================
# Script 20: ADAMTS15 vs C6ST1_ADAMTS15 — direct comparison
# Three analyses:
#   A) Direct pairwise contrast in ipsi/contra ratio model (script 19 data)
#   B) Bayesian posterior: P(C6ST1_ADAMTS15 effect > ADAMTS15 effect)
#      using saved brms M1 model from script 13
#   C) Effect size summary across all metrics
# Output: results/20_adamts15_vs_c6st1_adamts15/
# =============================================================================

library(data.table)
library(emmeans)
library(ggplot2)
library(brms)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------
RATIO_DIR <- "/path/to/results/19_ipsi_contra_ratio"
BRMS_DIR  <- "/path/to/results/13_bayesian_lmm"
OUT_DIR   <- "/path/to/results/20_adamts15_vs_c6st1_adamts15"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TREATMENT_ORDER <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")

PALETTE_2 <- c(
  ADAMTS15       = "#8e44ad",
  C6ST1_ADAMTS15 = "#27ae60"
)

# -----------------------------------------------------------------------------
# A. Direct pairwise contrast: ADAMTS15 vs C6ST1_ADAMTS15 in ratio model
# -----------------------------------------------------------------------------
cat("══════════════════════════════════════════════════════\n")
cat("Part A: Pairwise ratio contrast\n")
cat("══════════════════════════════════════════════════════\n")

animal_dat <- fread(file.path(RATIO_DIR, "animal_ratios.csv"))
animal_dat[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

# Log ratio model — all pairwise, no correction (single pre-specified contrast)
m_log  <- lm(log_ratio ~ treatment, data = animal_dat)
em_log <- emmeans(m_log, ~ treatment)
ct_log <- as.data.table(contrast(em_log, method = "pairwise", adjust = "none"))
ct_log[, sig := ifelse(p.value < 0.05, "*", "")]

key_log <- ct_log[contrast %in% c("ADAMTS15 - C6ST1_ADAMTS15",
                                   "C6ST1_ADAMTS15 - ADAMTS15")]
cat("\nDirect contrast: ADAMTS15 vs C6ST1_ADAMTS15 (log ratio)\n")
print(key_log[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

# Raw ratio model
m_raw  <- lm(ratio ~ treatment, data = animal_dat)
em_raw <- emmeans(m_raw, ~ treatment)
ct_raw <- as.data.table(contrast(em_raw, method = "pairwise", adjust = "none"))
ct_raw[, sig := ifelse(p.value < 0.05, "*", "")]

key_raw <- ct_raw[contrast %in% c("ADAMTS15 - C6ST1_ADAMTS15",
                                   "C6ST1_ADAMTS15 - ADAMTS15")]
cat("\nDirect contrast: ADAMTS15 vs C6ST1_ADAMTS15 (raw ratio)\n")
print(key_raw[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

# Group means
group_means <- animal_dat[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15"), .(
  mean_ratio     = mean(ratio),
  sd_ratio       = sd(ratio),
  mean_log_ratio = mean(log_ratio),
  n              = .N
), by = treatment]
cat("\nGroup means:\n")
print(group_means)

fwrite(ct_log,      file.path(OUT_DIR, "pairwise_log_contrasts.csv"))
fwrite(ct_raw,      file.path(OUT_DIR, "pairwise_ratio_contrasts.csv"))
fwrite(group_means, file.path(OUT_DIR, "group_means_ratio.csv"))

# -----------------------------------------------------------------------------
# B. Bayesian posterior: P(C6ST1_ADAMTS15 effect > ADAMTS15 effect)
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Part B: Bayesian posterior comparison\n")
cat("══════════════════════════════════════════════════════\n")

# Use lowercase filename as it appears on disk
m1_path <- file.path(BRMS_DIR, "m1_enwrap_weak.rds")
cat(sprintf("Loading model: %s\n", m1_path))
m1 <- readRDS(m1_path)

# Extract posterior draws
draws <- as.data.table(as_draws_df(m1))
cat(sprintf("Posterior draws: %d\n", nrow(draws)))

# Show all treatment-related columns so we can identify exact names
tx_cols <- grep("treatment", names(draws), value = TRUE)
cat("Treatment columns found:\n")
print(tx_cols)

# Identify the correct column names
col_a15   <- grep("treatmentADAMTS15$", names(draws), value = TRUE)
col_c6a15 <- grep("treatmentC6ST1_ADAMTS15", names(draws), value = TRUE)

if (length(col_a15) == 0 || length(col_c6a15) == 0) {
  cat("WARNING: Could not auto-detect column names.\n")
  cat("Please inspect 'tx_cols' above and set manually.\n")
  # Fallback: try positional match
  col_a15   <- tx_cols[grepl("ADAMTS15", tx_cols) & !grepl("C6ST1", tx_cols)][1]
  col_c6a15 <- tx_cols[grepl("C6ST1_ADAMTS15", tx_cols)][1]
}

col_a15   <- col_a15[1]
col_c6a15 <- col_c6a15[1]
cat(sprintf("\nUsing:\n  ADAMTS15:       %s\n  C6ST1_ADAMTS15: %s\n",
            col_a15, col_c6a15))

# Posterior probability that C6ST1_ADAMTS15 has larger negative effect
draws[, diff := get(col_c6a15) - get(col_a15)]  # negative = C6ST1_ADAMTS15 stronger

p_c6a15_stronger <- mean(draws$diff < 0)
cat(sprintf("\nP(C6ST1_ADAMTS15 effect > ADAMTS15 effect) = %.3f\n",
            p_c6a15_stronger))
cat(sprintf("P(effects equal or ADAMTS15 stronger)      = %.3f\n",
            1 - p_c6a15_stronger))

# Posterior summaries
post_summary <- data.table(
  parameter = c("ADAMTS15", "C6ST1_ADAMTS15", "difference (C6-A15)"),
  mean      = c(mean(draws[[col_a15]]),
                mean(draws[[col_c6a15]]),
                mean(draws$diff)),
  sd        = c(sd(draws[[col_a15]]),
                sd(draws[[col_c6a15]]),
                sd(draws$diff)),
  ci_89_lo  = c(quantile(draws[[col_a15]],  0.055),
                quantile(draws[[col_c6a15]], 0.055),
                quantile(draws$diff, 0.055)),
  ci_89_hi  = c(quantile(draws[[col_a15]],  0.945),
                quantile(draws[[col_c6a15]], 0.945),
                quantile(draws$diff, 0.945)),
  p_negative = c(mean(draws[[col_a15]]  < 0),
                 mean(draws[[col_c6a15]] < 0),
                 mean(draws$diff < 0))
)
cat("\nPosterior summary:\n")
print(post_summary)
fwrite(post_summary, file.path(OUT_DIR, "posterior_comparison.csv"))

# -----------------------------------------------------------------------------
# C. Effect size summary across all metrics
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Part C: Effect size comparison across metrics\n")
cat("══════════════════════════════════════════════════════\n")

metrics <- data.table(
  metric = c(
    "Enwrapment (M1 Bayesian p_neg)",
    "Enwrapment (ipsi/contra mean ratio)",
    "Spatial gradient: enwrapment at injection centre",
    "Spatial gradient: slope ratio vs mScarlet",
    "Composite score (M2 Bayesian p_neg)",
    "Hemisphere lateralisation (M3 Bayesian p_neg)",
    "Off-target specificity (M6 Bayesian p_neg)"
  ),
  ADAMTS15 = c(0.938, 0.674, 0.184, 3.1, 0.719, 0.999, 1.000),
  C6ST1_ADAMTS15 = c(0.951, 0.449, 0.107, 2.7, 0.725, 0.997, 0.999),
  higher_is_stronger = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
)

metrics[, C6ST1_ADAMTS15_stronger := fifelse(
  higher_is_stronger,
  C6ST1_ADAMTS15 > ADAMTS15,
  C6ST1_ADAMTS15 < ADAMTS15
)]

cat("\nEffect size comparison:\n")
print(metrics[, .(metric, ADAMTS15, C6ST1_ADAMTS15, C6ST1_ADAMTS15_stronger)])
cat(sprintf("\nC6ST1_ADAMTS15 numerically stronger in %d / %d metrics\n",
            sum(metrics$C6ST1_ADAMTS15_stronger), nrow(metrics)))

fwrite(metrics, file.path(OUT_DIR, "effect_size_comparison.csv"))

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

# Fig A: Posterior density overlay
plot_post <- rbind(
  data.table(value = draws[[col_a15]],   treatment = "ADAMTS15"),
  data.table(value = draws[[col_c6a15]], treatment = "C6ST1_ADAMTS15")
)

p1 <- ggplot(plot_post, aes(x = value, fill = treatment, colour = treatment)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = PALETTE_2) +
  scale_colour_manual(values = PALETTE_2) +
  labs(
    title    = "Posterior distributions: enwrapment effect vs mScarlet",
    subtitle = sprintf("P(C6ST1_ADAMTS15 stronger than ADAMTS15) = %.3f",
                       p_c6a15_stronger),
    x        = "Posterior beta (frac enwrapped vs mScarlet)",
    y        = "Density",
    fill     = NULL, colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_posterior_overlap.pdf"),
       p1, width = 7, height = 5)

# Fig B: Posterior density of the difference
p2 <- ggplot(draws, aes(x = diff)) +
  geom_density(fill = "#2c3e50", alpha = 0.6, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  annotate("text",
           x     = quantile(draws$diff, 0.02),
           y     = Inf,
           label = sprintf("P(diff < 0) = %.3f", p_c6a15_stronger),
           hjust = 0, vjust = 1.5, size = 4) +
  labs(
    title    = "Posterior: C6ST1_ADAMTS15 effect minus ADAMTS15 effect",
    subtitle = "Negative = C6ST1_ADAMTS15 has larger reduction",
    x        = "Difference in posterior beta",
    y        = "Density"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "fig_posterior_difference.pdf"),
       p2, width = 7, height = 5)

# Fig C: Dot plot of ratios, two groups only
two_groups <- animal_dat[treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15")]
two_groups[, treatment := factor(treatment,
                                  levels = c("ADAMTS15", "C6ST1_ADAMTS15"))]

p3 <- ggplot(two_groups, aes(x = treatment, y = ratio, colour = treatment)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.08, size = 4, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.3, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE_2) +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(
    title  = "Ipsi/contra ratio: ADAMTS15 vs C6ST1_ADAMTS15",
    x      = NULL,
    y      = "Ipsi / contra enwrapment ratio",
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "fig_ratio_two_groups.pdf"),
       p3, width = 4.5, height = 5)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 20 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  pairwise_log_contrasts.csv\n")
cat("  pairwise_ratio_contrasts.csv\n")
cat("  group_means_ratio.csv\n")
cat("  posterior_comparison.csv\n")
cat("  effect_size_comparison.csv\n")
cat("  fig_posterior_overlap.pdf\n")
cat("  fig_posterior_difference.pdf\n")
cat("  fig_ratio_two_groups.pdf\n")
cat("══════════════════════════════════════════════════════\n")