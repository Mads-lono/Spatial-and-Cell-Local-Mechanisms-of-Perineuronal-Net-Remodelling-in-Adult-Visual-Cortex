# =============================================================================
# Script 23: Bayesian analysis of C6ST1 + ADAMTS15 additive effect (RQ3)
# Research question: Does combining C6ST1 with ADAMTS15 enhance PNN
#   degradation beyond ADAMTS15 alone, consistent with substrate conditioning?
#
# Three analyses in order of complexity:
#
#   LM:  log_ratio ~ treatment (ref = ADAMTS15)
#        Single pre-specified contrast: C6ST1_ADAMTS15 - ADAMTS15
#        Direct frequentist analog of script 04's MD contrast.
#        No correction applied — single pre-specified test.
#
#   M7:  log_ratio ~ treatment_ref15 (brms, ref = ADAMTS15)
#        Bayesian version of LM above.
#        C6ST1_ADAMTS15 coefficient = posterior of enhancement over ADAMTS15.
#
#   M8:  log_ratio ~ ADAMTS15_present * C6ST1_present (brms)
#        Factorial ANCOVA encoding the biological hypothesis directly.
#        Interaction term = enhancement when both components present.
#        Groups: ADAMTS4 (neither), ADAMTS15 (A15 only),
#                C6ST1 (C6ST1 only), C6ST1_ADAMTS15 (both)
#        P(interaction < 0) = P(combination > sum of individual contributions)
#
# All three models use log(ipsi/contra ratio) as outcome — removes between-
# animal baseline variance that swamps effects in absolute enwrapment.
#
# Input:
#   19_ipsi_contra_ratio/animal_ratios.csv
# Output: results/23_bayesian_additive/
# =============================================================================

library(data.table)
library(emmeans)
library(brms)
library(posterior)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Paths & constants
# -----------------------------------------------------------------------------
BASE    <- "/path/to/results"
OUT_DIR <- file.path(BASE, "23_bayesian_additive")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

RATIO_PATH <- file.path(BASE, "19_ipsi_contra_ratio", "animal_ratios.csv")

TREATMENT_ORDER <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")

PALETTE <- c(
  mScarlet       = "#c0392b",
  ADAMTS4        = "#2980b9",
  ADAMTS4_MD     = "#1abc9c",
  ADAMTS15       = "#8e44ad",
  C6ST1          = "#e67e22",
  C6ST1_ADAMTS15 = "#27ae60"
)

# brms settings
MCMC_ITER    <- 4000L
MCMC_WARMUP  <- 1000L
MCMC_CHAINS  <- 4L
MCMC_CORES   <- 4L
SEED         <- 42L

# M8-specific settings (needs higher adapt_delta for small n + interaction)
M8_ITER      <- 8000L
M8_WARMUP    <- 2000L
M8_ADAPT     <- 0.99

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------
cat("Loading animal_ratios.csv...\n")
ratios <- fread(RATIO_PATH)
ratios[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("Animals loaded: %d\n", nrow(ratios)))
print(ratios[order(treatment),
             .(mouse_id, treatment,
               ratio     = round(ratio, 3),
               log_ratio = round(log_ratio, 3))])

# Relevel for LM and M7: ADAMTS15 as reference
# → C6ST1_ADAMTS15 coefficient directly = enhancement over ADAMTS15
ratios[, treatment_ref15 := relevel(treatment, ref = "ADAMTS15")]

# M8 data: four groups relevant to factorial design
# ADAMTS4 = neither component (reference)
# ADAMTS15 = A15 only; C6ST1 = C6ST1 only; C6ST1_ADAMTS15 = both
m8_groups <- c("ADAMTS4", "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")
m8_dat    <- ratios[treatment %in% m8_groups]
m8_dat[, ADAMTS15_present := as.integer(treatment %in% c("ADAMTS15", "C6ST1_ADAMTS15"))]
m8_dat[, C6ST1_present    := as.integer(treatment %in% c("C6ST1",    "C6ST1_ADAMTS15"))]

cat(sprintf("\nM8 data: %d animals across 4 groups\n", nrow(m8_dat)))
print(m8_dat[order(treatment),
             .(mouse_id, treatment, log_ratio = round(log_ratio, 3),
               ADAMTS15_present, C6ST1_present)])

# =============================================================================
# LM: log_ratio ~ treatment (ref = ADAMTS15)
# Single pre-specified contrast — no correction applied
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("LM: Pre-specified contrast — C6ST1_ADAMTS15 vs ADAMTS15\n")
cat("══════════════════════════════════════════════════════\n")

m_lm <- lm(log_ratio ~ treatment_ref15, data = ratios)
cat("\nModel summary (ref = ADAMTS15):\n")
print(summary(m_lm))

em_lm <- emmeans(m_lm, ~ treatment_ref15)
ct_lm <- contrast(em_lm,
                  method = list("C6ST1_ADAMTS15 - ADAMTS15" = c(0, 0, 0, -1, 0, 1)),
                  adjust = "none")
ct_lm_dt <- as.data.table(summary(ct_lm))
ct_lm_dt[, sig := ifelse(p.value < 0.05, "*", "")]
cat("\nPre-specified contrast (no correction — single test):\n")
print(ct_lm_dt[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

fwrite(ct_lm_dt, file.path(OUT_DIR, "lm_direct_contrast.csv"))

coef_dt <- as.data.table(summary(m_lm)$coefficients, keep.rownames = "term")
setnames(coef_dt, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
                  c("estimate", "SE", "t.ratio", "p.value"))
coef_dt[, sig := fcase(
  p.value < 0.001, "***",
  p.value < 0.01,  "**",
  p.value < 0.05,  "*",
  p.value < 0.10,  ".",
  default          = ""
)]
fwrite(coef_dt, file.path(OUT_DIR, "lm_full_coefficients.csv"))

# =============================================================================
# M7: log_ratio ~ treatment (ref = ADAMTS15), Bayesian
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("M7: Bayesian direct contrast — C6ST1_ADAMTS15 vs ADAMTS15\n")
cat("══════════════════════════════════════════════════════\n")

prior_m7 <- c(
  prior(normal(0, 0.5), class = Intercept),
  prior(normal(0, 0.2), class = b),
  prior(exponential(1), class = sigma)
)

m7_cache <- file.path(OUT_DIR, "m7_direct_contrast.rds")
if (file.exists(m7_cache)) {
  cat("Loading cached M7...\n")
  m7 <- readRDS(m7_cache)
} else {
  cat("Fitting M7...\n")
  m7 <- brm(
    log_ratio ~ treatment_ref15,
    data    = ratios,
    prior   = prior_m7,
    iter    = MCMC_ITER,
    warmup  = MCMC_WARMUP,
    chains  = MCMC_CHAINS,
    cores   = MCMC_CORES,
    seed    = SEED,
    control = list(adapt_delta = 0.95),
    backend = "rstan",
    file    = m7_cache
  )
}

cat("\nM7 summary:\n")
print(summary(m7))

draws_m7      <- as_draws_df(m7)
col_c6a15_m7  <- grep("treatment_ref15.*C6ST1_ADAMTS15",
                       names(draws_m7), value = TRUE)[1]
cat(sprintf("Enhancement column: %s\n", col_c6a15_m7))

enhancement_draws <- draws_m7[[col_c6a15_m7]]
p_enhancement     <- mean(enhancement_draws < 0)

m7_summary <- data.table(
  parameter      = "C6ST1_ADAMTS15 vs ADAMTS15 (log_ratio)",
  post_mean      = round(mean(enhancement_draws), 4),
  post_sd        = round(sd(enhancement_draws), 4),
  ci_89_lo       = round(quantile(enhancement_draws, 0.055), 4),
  ci_89_hi       = round(quantile(enhancement_draws, 0.945), 4),
  p_negative     = round(p_enhancement, 4),
  interpretation = fcase(
    p_enhancement > 0.95, "Strong evidence for enhancement",
    p_enhancement > 0.80, "Moderate evidence for enhancement",
    default              = "Weak/no evidence for enhancement"
  )
)

cat("\n── M7 result ──\n")
print(m7_summary)
cat(sprintf("\nP(C6ST1_ADAMTS15 reduces log_ratio more than ADAMTS15) = %.3f\n",
            p_enhancement))

fwrite(m7_summary, file.path(OUT_DIR, "m7_direct_contrast_result.csv"))

# =============================================================================
# M8: log_ratio ~ ADAMTS15_present * C6ST1_present (factorial ANCOVA)
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("M8: Factorial ANCOVA — substrate conditioning test\n")
cat("══════════════════════════════════════════════════════\n")
cat("Interaction term = extra reduction when both components present\n")
cat("H1: interaction < 0 (combination > sum of individual contributions)\n\n")

# Informative priors grounded in observed data ranges
prior_m8 <- c(
  prior(normal(-0.25, 0.15), class = Intercept),
  prior(normal(-0.15, 0.15), class = b, coef = "ADAMTS15_present"),
  prior(normal(0,     0.10), class = b, coef = "C6ST1_present"),
  prior(normal(-0.20, 0.15), class = b, coef = "ADAMTS15_present:C6ST1_present"),
  prior(exponential(2),      class = sigma)
)

m8_cache <- file.path(OUT_DIR, "m8_factorial_ancova.rds")
if (file.exists(m8_cache)) {
  cat("Loading cached M8...\n")
  m8 <- readRDS(m8_cache)
} else {
  cat("Fitting M8...\n")
  m8 <- brm(
    log_ratio ~ ADAMTS15_present * C6ST1_present,
    data    = m8_dat,
    prior   = prior_m8,
    iter    = M8_ITER,
    warmup  = M8_WARMUP,
    chains  = MCMC_CHAINS,
    cores   = MCMC_CORES,
    seed    = SEED,
    control = list(adapt_delta = M8_ADAPT, max_treedepth = 15),
    backend = "rstan",
    file    = m8_cache
  )
}

cat("\nM8 summary:\n")
print(summary(m8))

draws_m8    <- as_draws_df(m8)
inter_draws <- draws_m8[["b_ADAMTS15_present:C6ST1_present"]]
a15_draws   <- draws_m8[["b_ADAMTS15_present"]]
c6_draws    <- draws_m8[["b_C6ST1_present"]]

pred_neither  <- draws_m8$b_Intercept
pred_a15_only <- pred_neither + a15_draws
pred_c6_only  <- pred_neither + c6_draws
pred_both     <- pred_neither + a15_draws + c6_draws + inter_draws

p_inter_neg <- mean(inter_draws < 0)
p_c6_zero   <- mean(abs(c6_draws) < 0.05)

cat(sprintf("\nP(interaction < 0) = %.3f\n", p_inter_neg))
cat(sprintf("P(C6ST1 alone ~ 0) = %.3f\n", p_c6_zero))
cat(sprintf("Interaction: mean=%.3f, 89%% CI [%.3f, %.3f]\n",
            mean(inter_draws),
            quantile(inter_draws, 0.055),
            quantile(inter_draws, 0.945)))

if (p_inter_neg > 0.95) {
  cat("→ Strong evidence: combination > sum of individual contributions\n")
} else if (p_inter_neg > 0.80) {
  cat("→ Moderate evidence: combination > sum of individual contributions\n")
} else {
  cat("→ Weak evidence for super-additive interaction\n")
}

# Parameter summary table
m8_result <- data.table(
  parameter = c("ADAMTS4 (intercept)", "ADAMTS15_present",
                "C6ST1_present",       "Interaction (A15 x C6ST1)"),
  post_mean = round(c(mean(pred_neither), mean(a15_draws),
                      mean(c6_draws),     mean(inter_draws)), 4),
  post_sd   = round(c(sd(pred_neither),   sd(a15_draws),
                      sd(c6_draws),        sd(inter_draws)), 4),
  ci_89_lo  = round(c(quantile(pred_neither,  0.055),
                      quantile(a15_draws,     0.055),
                      quantile(c6_draws,      0.055),
                      quantile(inter_draws,   0.055)), 4),
  ci_89_hi  = round(c(quantile(pred_neither,  0.945),
                      quantile(a15_draws,     0.945),
                      quantile(c6_draws,      0.945),
                      quantile(inter_draws,   0.945)), 4)
)

cat("\n── M8 parameter summary ──\n")
print(m8_result)
fwrite(m8_result, file.path(OUT_DIR, "m8_factorial_result.csv"))

# Scalar summary for easy reporting
m8_scalars <- data.table(
  metric = c(
    "P(interaction < 0)",
    "P(C6ST1 alone ~ 0, |effect| < 0.05)",
    "Interaction mean",
    "Interaction 89% CI lo",
    "Interaction 89% CI hi",
    "ADAMTS15 effect mean",
    "C6ST1 effect mean",
    "Predicted log_ratio: ADAMTS4",
    "Predicted log_ratio: ADAMTS15",
    "Predicted log_ratio: C6ST1",
    "Predicted log_ratio: C6ST1_ADAMTS15"
  ),
  value = round(c(
    p_inter_neg, p_c6_zero,
    mean(inter_draws),
    quantile(inter_draws, 0.055),
    quantile(inter_draws, 0.945),
    mean(a15_draws), mean(c6_draws),
    mean(pred_neither), mean(pred_a15_only),
    mean(pred_c6_only), mean(pred_both)
  ), 4)
)
fwrite(m8_scalars, file.path(OUT_DIR, "m8_scalar_summary.csv"))

# =============================================================================
# Figures
# =============================================================================
cat("\nGenerating figures...\n")

# Fig 1: LM coefficient plot — all treatments vs ADAMTS15
coef_plot_dt <- coef_dt[term != "(Intercept)"]
coef_plot_dt[, treatment := gsub("treatment_ref15", "", term)]
coef_plot_dt[, treatment := factor(treatment, levels = setdiff(TREATMENT_ORDER, "ADAMTS15"))]
coef_plot_dt[, colour := ifelse(p.value < 0.05, "significant", "ns")]

p0 <- ggplot(coef_plot_dt,
             aes(x = estimate, y = treatment, colour = colour)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * SE,
                     xmax = estimate + 1.96 * SE),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(significant = "#e74c3c", ns = "grey50")) +
  labs(
    title    = "LM: Treatment effects vs ADAMTS15 (ref)",
    subtitle = "Outcome: log(ipsi/contra enwrapment ratio)",
    x        = "Coefficient (vs ADAMTS15)",
    y        = NULL,
    colour   = NULL,
    caption  = "Error bars = 95% CI; single pre-specified contrast, no correction"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "fig_lm_coefficients.pdf"), p0, width = 7, height = 5)

# Fig 2: M7 posterior density — enhancement of C6ST1_ADAMTS15 over ADAMTS15
p1 <- ggplot(data.frame(x = enhancement_draws), aes(x = x)) +
  geom_density(fill = "#27ae60", alpha = 0.5, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  annotate("text",
           x     = quantile(enhancement_draws, 0.05),
           y     = Inf,
           label = sprintf("P(enhancement) = %.3f", p_enhancement),
           hjust = 0, vjust = 1.5, size = 4) +
  labs(
    title    = "M7: Posterior — C6ST1_ADAMTS15 vs ADAMTS15 (log ratio)",
    subtitle = "Negative = C6ST1_ADAMTS15 reduces log ratio more than ADAMTS15 alone",
    x        = "Posterior beta (C6ST1_ADAMTS15 − ADAMTS15)",
    y        = "Density",
    caption  = "brms: log_ratio ~ treatment(ref=ADAMTS15)"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "fig_m7_posterior.pdf"), p1, width = 7, height = 5)

# Fig 3: M8 interaction posterior
p2 <- ggplot(data.frame(x = inter_draws), aes(x = x)) +
  geom_density(fill = "#2c3e50", alpha = 0.5, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  annotate("text",
           x     = quantile(inter_draws, 0.05),
           y     = Inf,
           label = sprintf("P(interaction < 0) = %.3f", p_inter_neg),
           hjust = 0, vjust = 1.5, size = 4) +
  labs(
    title    = "M8: Posterior — ADAMTS15 x C6ST1 interaction",
    subtitle = "Negative = combination exceeds sum of individual contributions\nC6ST1-alone effect near zero by design",
    x        = "Posterior: ADAMTS15_present x C6ST1_present interaction",
    y        = "Density",
    caption  = "brms: log_ratio ~ ADAMTS15_present * C6ST1_present"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "fig_m8_interaction.pdf"), p2, width = 7, height = 5)

# Fig 4: M8 predicted group posteriors with observed data
pred_long <- rbind(
  data.table(group = "ADAMTS4\n(neither)",     value = as.numeric(pred_neither)),
  data.table(group = "ADAMTS15\n(A15 only)",   value = as.numeric(pred_a15_only)),
  data.table(group = "C6ST1\n(C6ST1 only)",    value = as.numeric(pred_c6_only)),
  data.table(group = "C6ST1_ADAMTS15\n(both)", value = as.numeric(pred_both))
)
group_levels <- c("ADAMTS4\n(neither)", "ADAMTS15\n(A15 only)",
                  "C6ST1\n(C6ST1 only)", "C6ST1_ADAMTS15\n(both)")
pred_long[, group := factor(group, levels = group_levels)]

obs_long <- m8_dat[, .(
  group = fcase(
    treatment == "ADAMTS4",        "ADAMTS4\n(neither)",
    treatment == "ADAMTS15",       "ADAMTS15\n(A15 only)",
    treatment == "C6ST1",          "C6ST1\n(C6ST1 only)",
    treatment == "C6ST1_ADAMTS15", "C6ST1_ADAMTS15\n(both)"
  ),
  value = log_ratio
)]
obs_long[, group := factor(group, levels = group_levels)]

p3 <- ggplot(pred_long, aes(x = group, y = value)) +
  geom_violin(fill = "steelblue", alpha = 0.35, linewidth = 0.5) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.3, linewidth = 0.8, colour = "navy") +
  geom_jitter(data = obs_long, aes(x = group, y = value),
              width = 0.07, size = 3.5, colour = "black", alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title    = "M8: Posterior predicted vs observed log(ipsi/contra) by group",
    subtitle = "Violin = posterior predictive; points = observed animals; crossbar = posterior mean",
    x        = NULL,
    y        = "log(ipsi/contra enwrapment ratio)",
    caption  = "brms: log_ratio ~ ADAMTS15_present * C6ST1_present"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 10))

ggsave(file.path(OUT_DIR, "fig_m8_group_posteriors.pdf"), p3, width = 7, height = 5)

# =============================================================================
# Summary
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 23 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("\nFiles written:\n")
cat("  lm_direct_contrast.csv\n")
cat("  lm_full_coefficients.csv\n")
cat("  m7_direct_contrast_result.csv\n")
cat("  m8_factorial_result.csv\n")
cat("  m8_scalar_summary.csv\n")
cat("  fig_lm_coefficients.pdf\n")
cat("  fig_m7_posterior.pdf\n")
cat("  fig_m8_interaction.pdf\n")
cat("  fig_m8_group_posteriors.pdf\n")
cat("══════════════════════════════════════════════════════\n")