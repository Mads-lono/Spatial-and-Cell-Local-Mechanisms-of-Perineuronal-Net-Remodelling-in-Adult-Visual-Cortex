# =============================================================================
# 01_primary_LMM.R
# =============================================================================
# INPUT:  results/00_enwrapment/slice_enwrapment.csv
#         results/00_enwrapment/animal_enwrapment.csv
# OUTPUT: results/01_primary_LMM/
#   lmm_results.csv            — LMM on ratio (primary)
#   lmm_log_results.csv        — LMM on log(ratio) (sensitivity)
#   residual_check.pdf/png     — QQ + residual plots for both models
#   fig1_forest.pdf/png        — forest plot, both models side by side
#   fig2_animal_dots.pdf/png   — raw dots + LMM-predicted means and 95% CI
#
# MODELS
#   Primary:     ratio      ~ treatment + (1 | mouse_id)
#   Sensitivity: log(ratio) ~ treatment + (1 | mouse_id)
#     β on log scale → exponentiate to get multiplicative fold-change
#     β = 0 on log scale ≡ ratio = 1.0 (no change)
#
# Both models use the reference group specified by TREAT_ORDER[1],
# REML = TRUE, and Satterthwaite degrees of freedom.
# Singular fit is reported explicitly.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(patchwork)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS_DIR <- "/path/to/results"
ENW_DIR     <- file.path(RESULTS_DIR, "00_enwrapment")
OUT_DIR     <- file.path(RESULTS_DIR, "01_primary_LMM")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PALETTE <- c(
  ADAMTS4        = "#4e9af1",
  ADAMTS4_MD     = "#f17c4e",
  ADAMTS15       = "#4ef196",
  C6ST1          = "#c44ef1",
  C6ST1_ADAMTS15 = "#f1c44e",
  mScarlet       = "#888888"
)
TREAT_ORDER <- c("mScarlet", "C6ST1", "ADAMTS4", "ADAMTS4_MD",
                 "ADAMTS15", "C6ST1_ADAMTS15")

# =============================================================================
# 1. LOAD
# =============================================================================
slice_df  <- fread(file.path(ENW_DIR, "slice_enwrapment.csv"))
animal_df <- fread(file.path(ENW_DIR, "animal_enwrapment.csv"))

slice_df[,  treatment := factor(treatment, levels = TREAT_ORDER)]
animal_df[, treatment := factor(treatment, levels = TREAT_ORDER)]

# Log transform requires ratio > 0; report and remove any zero/negative rows
n_bad <- sum(slice_df$ratio <= 0, na.rm = TRUE)
if (n_bad > 0)
  cat(sprintf("  Removing %d slices with ratio <= 0 before log model.\n", n_bad))
slice_df_log <- slice_df[ratio > 0]
slice_df_log[, log_ratio := log(ratio)]

cat(sprintf("Slices (ratio model): %d  |  Slices (log model): %d  |  Animals: %d\n",
            nrow(slice_df), nrow(slice_df_log), nrow(animal_df)))

# =============================================================================
# 2. FIT BOTH MODELS
# =============================================================================
cat("\n── Fitting models ───────────────────────────────────────────────────────\n")

fit_ratio <- lmer(
  ratio ~ treatment + (1 | mouse_id),
  data    = slice_df,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

fit_log <- lmer(
  log_ratio ~ treatment + (1 | mouse_id),
  data    = slice_df_log,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

# ── Singular fit check ────────────────────────────────────────────────────────
check_singular <- function(fit, name) {
  if (isSingular(fit)) {
    cat(sprintf("\n  *** SINGULAR FIT: %s ***\n", name))
    cat("  Between-animal variance estimated at zero.\n")
    cat("  Model is collapsing to a standard linear model.\n")
    cat("  p-values remain valid but random effect is uninformative.\n\n")
  } else {
    cat(sprintf("  %s: no singularity.\n", name))
  }
}
check_singular(fit_ratio, "ratio model")
check_singular(fit_log,   "log(ratio) model")

# ── Helper: extract fixed effects ─────────────────────────────────────────────
extract_coefs <- function(fit) {
  coefs <- as.data.frame(coef(summary(fit, ddf = "Satterthwaite")))
  coefs$term <- rownames(coefs)
  setDT(coefs)
  res <- coefs[term != "(Intercept)"]
  setnames(res,
    c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)"),
    c("beta",     "se",         "df", "t",       "p"))
  res[, treatment := gsub("^treatment", "", term)]
  res[, ci_lo := beta - qt(0.975, df) * se]
  res[, ci_hi := beta + qt(0.975, df) * se]
  res[, sig := fcase(
    p < 0.001, "***",
    p < 0.01,  "**",
    p < 0.05,  "*",
    p < 0.10,  ".",
    default =  "ns"
  )]
  res[, .(treatment, beta, se, ci_lo, ci_hi, t, df, p, sig)]
}

res_ratio <- extract_coefs(fit_ratio)
res_log   <- extract_coefs(fit_log)
res_log[, fold_change := exp(beta)]
res_log[, fc_ci_lo    := exp(ci_lo)]
res_log[, fc_ci_hi    := exp(ci_hi)]
res_log[, pct_change  := (fold_change - 1) * 100]
res_log[, pct_ci_lo   := (fc_ci_lo   - 1) * 100]
res_log[, pct_ci_hi   := (fc_ci_hi   - 1) * 100]

fwrite(res_ratio, file.path(OUT_DIR, "lmm_results.csv"))
fwrite(res_log,   file.path(OUT_DIR, "lmm_log_results.csv"))

cat("\n=== Primary: ratio model ===\n")
print(res_ratio)
cat("\n=== Sensitivity: log(ratio) model ===\n")
print(res_log[, .(treatment, beta, p, sig, fold_change, pct_change, pct_ci_lo, pct_ci_hi)])

# ── Variance components ───────────────────────────────────────────────────────
# ICC = between-animal variance / total variance.
# ICC > 0.5: animal differences dominate.
# ICC < 0.1: very consistent across animals.
cat("\n── Variance components ──────────────────────────────────────────────────\n")
for (nm in c("ratio", "log")) {
  fit <- if (nm == "ratio") fit_ratio else fit_log
  vc  <- as.data.frame(VarCorr(fit))
  between <- vc[vc$grp == "mouse_id", "vcov"]
  total   <- sum(vc$vcov)
  icc     <- between / total
  interp  <- if (icc > 0.5) "animal differences dominate"
             else if (icc > 0.2) "moderate between-animal variance"
             else if (icc > 0.1) "low-moderate between-animal variance"
             else "very consistent across animals"
  cat(sprintf("  ICC %-12s: %.3f  (%.0f%% between-animal) — %s\n",
              nm, icc, icc * 100, interp))
}

# =============================================================================
# 3. RESIDUAL DIAGNOSTICS
# =============================================================================
cat("\n── Saving residual diagnostics ──────────────────────────────────────────\n")

diag_panel <- function(fit, title) {
  df <- data.frame(
    fitted    = fitted(fit),
    std_resid = as.numeric(scale(resid(fit, type = "pearson")))
  )
  p_qq <- ggplot(df, aes(sample = std_resid)) +
    stat_qq(size = 1.2, alpha = 0.4) +
    stat_qq_line(colour = "red", linewidth = 0.7) +
    labs(title = paste("QQ —", title),
         x = "Theoretical", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_rv <- ggplot(df, aes(x = fitted, y = std_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1, alpha = 0.35) +
    geom_smooth(method = "loess", se = FALSE, colour = "red",
                linewidth = 0.6, formula = y ~ x) +
    labs(title = paste("Resid vs fitted —", title),
         x = "Fitted", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_qq + p_rv
}

p_diag <- diag_panel(fit_ratio, "ratio") /
          diag_panel(fit_log,   "log(ratio)")

ggsave(file.path(OUT_DIR, "residual_check.pdf"), p_diag,
       width = 9, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "residual_check.png"), p_diag,
       width = 9, height = 8, dpi = 150)
cat("  Saved residual_check — inspect QQ plots before choosing primary model.\n")

# =============================================================================
# 4. FIGURE 1 — Forest plots, both models side by side
# =============================================================================
PLOT_ORDER <- c("C6ST1_ADAMTS15", "ADAMTS15", "ADAMTS4_MD", "ADAMTS4", "C6ST1")

make_forest <- function(res, x_label, title) {
  d <- copy(res[treatment %in% PLOT_ORDER])
  d[, treatment := factor(treatment, levels = PLOT_ORDER)]
  x_range <- diff(range(c(d$ci_lo, d$ci_hi)))
  ggplot(d, aes(x = beta, y = treatment, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, y = treatment),
                  orientation = "y", width = 0.25, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(x = ci_hi + x_range * 0.05, label = sig),
              hjust = 0, size = 4, fontface = "bold", colour = "black") +
    scale_colour_manual(values = PALETTE, guide = "none") +
    scale_x_continuous(name   = x_label,
                       expand = expansion(mult = c(0.05, 0.15))) +
    scale_y_discrete(name = NULL) +
    labs(title    = title,
         subtitle = "vs reference  |  Satterthwaite df  |  95% CI") +
    theme_bw(base_size = 10) +
    theme(plot.title       = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank())
}

p1 <- make_forest(res_ratio, expression(beta ~ "(ratio scale)"),
                  "ratio ~ treatment + (1|mouse_id)") |
      make_forest(res_log,   expression(beta ~ "(log ratio scale)"),
                  "log(ratio) ~ treatment + (1|mouse_id)") +
  plot_annotation(
    title    = "Treatment effects on PV enwrapment — whole hemisphere",
    subtitle = "Left: untransformed  |  Right: log-transformed",
    theme    = theme(plot.title    = element_text(face = "bold"),
                     plot.subtitle = element_text(colour = "grey40"))
  )

# =============================================================================
# 5. FIGURE 2 — Raw dots + LMM-predicted means and 95% CI
# =============================================================================
new_data <- data.table(
  treatment = factor(TREAT_ORDER, levels = TREAT_ORDER),
  mouse_id  = NA_character_
)
pred        <- predict(fit_ratio, newdata = new_data, re.form = NA, se.fit = TRUE)
med_df      <- median(res_ratio$df)
new_data[, pred_mean := pred$fit]
new_data[, pred_lo   := pred$fit - qt(0.975, med_df) * pred$se.fit]
new_data[, pred_hi   := pred$fit + qt(0.975, med_df) * pred$se.fit]

p2 <- ggplot(animal_df, aes(x = treatment, y = ratio, colour = treatment)) +
  geom_hline(yintercept = 1, linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.85) +
  geom_errorbar(data = new_data,
                aes(y = pred_mean, ymin = pred_lo, ymax = pred_hi),
                colour = "black", width = 0.3, linewidth = 1.1) +
  geom_point(data = new_data, aes(y = pred_mean),
             colour = "black", size = 4.5, shape = 18) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Enwrapment ratio (ipsi / contra)") +
  labs(
    title    = "Animal-level enwrapment ratios",
    subtitle = "Points = individual animals  |  Diamond + bar = LMM-predicted mean ± 95% CI"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 30, hjust = 1),
        plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# =============================================================================
# 6. SAVE
# =============================================================================
for (ext in c("pdf", "png")) {
  dpi <- 300
  ggsave(file.path(OUT_DIR, paste0("fig1_forest.",      ext)), p1,
         width = 11, height = 4, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig2_animal_dots.", ext)), p2,
         width = 7,  height = 5, dpi = dpi)
}

cat(sprintf("\n── All outputs → %s\n", OUT_DIR))
