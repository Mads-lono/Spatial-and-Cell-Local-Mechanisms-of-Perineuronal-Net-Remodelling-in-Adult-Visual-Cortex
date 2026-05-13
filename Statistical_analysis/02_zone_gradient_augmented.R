# =============================================================================
# 02_zone_gradient.R
# =============================================================================
# INPUT:  results/00_enwrapment/zone_enwrapment.csv
# OUTPUT: results/02_zone_gradient/
#   zone_interaction_anova.csv   — Type III ANOVA: treatment*zone interaction
#   zone_contrasts_raw.csv       — emmeans zone contrasts, raw scale, BH-FDR
#   zone_contrasts_logit.csv     — emmeans zone contrasts, logit scale, BH-FDR
#   residual_check.pdf/png
#   fig3_zone_forest.pdf/png
#   fig4_zone_dots.pdf/png
#
# MODELS (single model across all zones — pools variance)
#
#   Primary:
#     frac_enwrapped ~ treatment * zone + (1 | mouse_id)
#
#   Sensitivity:
#     logit(frac_enwrapped + offset) ~ treatment * zone + (1 | mouse_id)
#     Continuity offset retains all records including zero fractions.
#
# FDR: BH correction within each treatment across the 3 zone contrasts.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  if (!requireNamespace("emmeans", quietly = TRUE))
    install.packages("emmeans", repos = "https://cloud.r-project.org")
  library(emmeans)
  library(ggplot2)
  library(patchwork)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS_DIR <- "/path/to/results"
ENW_DIR     <- file.path(RESULTS_DIR, "00_enwrapment")
OUT_DIR     <- file.path(RESULTS_DIR, "02_zone_gradient")
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
ZONE_ORDER  <- c("Core", "Penumbra", "Outside")
PLOT_ORDER  <- c("C6ST1_ADAMTS15", "ADAMTS15", "ADAMTS4_MD", "ADAMTS4", "C6ST1")

# =============================================================================
# 1. LOAD & TRANSFORM
# =============================================================================
zone_df <- fread(file.path(ENW_DIR, "zone_enwrapment.csv"))
zone_df[, treatment := factor(treatment, levels = TREAT_ORDER)]
zone_df[, zone      := factor(zone,      levels = ZONE_ORDER)]

# Continuity offset for logit — retains all records including zero fractions
LOGIT_OFFSET <- 0.01
zone_df[, logit_fe := log(
  (frac_enwrapped + LOGIT_OFFSET) / (1 - frac_enwrapped + LOGIT_OFFSET)
)]

n_zero <- sum(zone_df$frac_enwrapped <= 0)
cat(sprintf("Records: %d  |  zero fractions handled via logit offset (%.2f): %d kept\n",
            nrow(zone_df), LOGIT_OFFSET, n_zero))

# =============================================================================
# 2. FIT INTERACTION MODELS
# =============================================================================
cat("\n── Fitting treatment × zone interaction models ───────────────────────────\n")

fit_raw <- lmer(
  frac_enwrapped ~ treatment * zone + (1 | mouse_id),
  data    = zone_df,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

fit_logit <- lmer(
  logit_fe ~ treatment * zone + (1 | mouse_id),
  data    = zone_df,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

# ── Singular fit check ────────────────────────────────────────────────────────
for (nm in c("raw", "logit")) {
  fit <- if (nm == "raw") fit_raw else fit_logit
  if (isSingular(fit)) {
    cat(sprintf("  *** SINGULAR FIT [%s] — between-animal variance → 0 ***\n", nm))
  } else {
    cat(sprintf("  %s model: no singularity.\n", nm))
  }
}

# ── Variance components ───────────────────────────────────────────────────────
cat("\n── Variance components ──────────────────────────────────────────────────\n")
for (nm in c("raw", "logit")) {
  fit     <- if (nm == "raw") fit_raw else fit_logit
  vc      <- as.data.frame(VarCorr(fit))
  between <- vc[vc$grp == "mouse_id", "vcov"]
  icc     <- between / sum(vc$vcov)
  interp  <- if (icc > 0.5) "animal differences dominate"
             else if (icc > 0.2) "moderate between-animal variance"
             else if (icc > 0.1) "low-moderate between-animal variance"
             else "very consistent across animals"
  cat(sprintf("  ICC %-8s: %.3f (%.0f%%) — %s\n", nm, icc, icc * 100, interp))
}

# =============================================================================
# 3. TYPE III ANOVA — interaction significance
# =============================================================================
cat("\n── Type III ANOVA (Satterthwaite) ───────────────────────────────────────\n")
cat("  Key row: treatment:zone — is the effect spatially graded?\n\n")
anova_raw   <- as.data.frame(anova(fit_raw,   ddf = "Satterthwaite"))
anova_logit <- as.data.frame(anova(fit_logit, ddf = "Satterthwaite"))

anova_raw$term   <- rownames(anova_raw)
anova_logit$term <- rownames(anova_logit)

cat("Raw model:\n");   print(anova_raw[,   c("term","NumDF","DenDF","F value","Pr(>F)")])
cat("\nLogit model:\n"); print(anova_logit[, c("term","NumDF","DenDF","F value","Pr(>F)")])

fwrite(setDT(anova_raw)[, model := "raw"],
       file.path(OUT_DIR, "zone_interaction_anova.csv"))

# =============================================================================
# 4. EMMEANS — zone-specific contrasts vs reference
# =============================================================================
extract_zone_contrasts <- function(fit, model_name) {
  emm  <- emmeans(fit, ~ treatment | zone)
  cont <- contrast(emm, method = "trt.vs.ctrl", ref = "mScarlet")
  res  <- as.data.table(summary(cont, infer = TRUE))
  setnames(res,
    c("estimate", "SE", "df", "lower.CL", "upper.CL", "t.ratio", "p.value"),
    c("beta",     "se", "df", "ci_lo",    "ci_hi",    "t",       "p"),
    skip_absent = TRUE)

  # BH-FDR within each treatment across zones
  res[, p_adj := p.adjust(p, method = "BH"), by = contrast]
  res[, sig := fcase(
    p_adj < 0.001, "***",
    p_adj < 0.01,  "**",
    p_adj < 0.05,  "*",
    p_adj < 0.10,  ".",
    default =      "ns"
  )]
  res[, treatment := gsub(" - mScarlet", "", contrast)]
  res[, model := model_name]
  res[, .(model, zone, treatment, beta, se, ci_lo, ci_hi, t, df, p, p_adj, sig)]
}

cat("\n── Extracting emmeans contrasts ─────────────────────────────────────────\n")
res_raw   <- extract_zone_contrasts(fit_raw,   "raw")
res_logit <- extract_zone_contrasts(fit_logit, "logit")

# For logit model add % change (back-transform via plogis)
intercepts_logit <- as.data.table(
  summary(emmeans(fit_logit, ~ zone, at = list(treatment = "mScarlet")))
)[, .(zone, mScarlet_logit = emmean)]
res_logit <- merge(res_logit, intercepts_logit, by = "zone")
res_logit[, pct_change := (plogis(mScarlet_logit + beta) /
                            plogis(mScarlet_logit) - 1) * 100]
res_logit[, mScarlet_logit := NULL]

fwrite(res_raw,   file.path(OUT_DIR, "zone_contrasts_raw.csv"))
fwrite(res_logit, file.path(OUT_DIR, "zone_contrasts_logit.csv"))

cat("\n=== Zone contrasts — raw scale (BH-FDR per treatment) ===\n")
print(res_raw[order(zone, treatment),
              .(zone, treatment, beta, ci_lo, ci_hi, p, p_adj, sig)])

cat("\n=== Zone contrasts — logit scale (% change from reference baseline) ===\n")
print(res_logit[order(zone, treatment),
                .(zone, treatment, beta, p, p_adj, sig, pct_change)])

cat("\n=== Agreement raw vs logit ===\n")
agree <- merge(
  res_raw[,   .(zone, treatment, sig_raw   = sig)],
  res_logit[, .(zone, treatment, sig_logit = sig)],
  by = c("zone", "treatment"))
agree[, agree := sig_raw == sig_logit]
print(agree[order(zone, treatment)])

# =============================================================================
# 5. RESIDUAL DIAGNOSTICS
# =============================================================================
diag_panel <- function(fit, title) {
  df <- data.frame(
    fitted    = fitted(fit),
    std_resid = as.numeric(scale(resid(fit, type = "pearson")))
  )
  p_qq <- ggplot(df, aes(sample = std_resid)) +
    stat_qq(size = 0.8, alpha = 0.3) +
    stat_qq_line(colour = "red", linewidth = 0.6) +
    labs(title = paste("QQ —", title), x = "Theoretical", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_rv <- ggplot(df, aes(x = fitted, y = std_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.7, alpha = 0.25) +
    geom_smooth(method = "loess", se = FALSE, colour = "red",
                linewidth = 0.5, formula = y ~ x) +
    labs(title = paste("Resid —", title), x = "Fitted", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_qq + p_rv
}

p_diag <- diag_panel(fit_raw,   "raw interaction model") /
          diag_panel(fit_logit, "logit interaction model") +
  plot_annotation(
    title = "Residual diagnostics — treatment × zone interaction models",
    theme = theme(plot.title = element_text(face = "bold"))
  )

ggsave(file.path(OUT_DIR, "residual_check.pdf"), p_diag,
       width = 9, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "residual_check.png"), p_diag,
       width = 9, height = 8, dpi = 150)

# =============================================================================
# 6. FIGURE 3 — Forest plot (emmeans contrasts, both models)
# =============================================================================
make_forest_emm <- function(res, x_label, model_name) {
  d <- copy(res[treatment %in% PLOT_ORDER])
  d[, treatment := factor(treatment, levels = PLOT_ORDER)]
  d[, zone      := factor(zone, levels = ZONE_ORDER)]

  d[, x_offset := {
    r <- diff(range(c(ci_lo, ci_hi)))
    rep(r * 0.06, .N)
  }, by = zone]

  ggplot(d, aes(x = beta, y = treatment, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, y = treatment),
                  orientation = "y", width = 0.25, linewidth = 0.7) +
    geom_point(size = 3) +
    geom_text(aes(x = ci_hi + x_offset, label = sig),
              hjust = 0, size = 3.5, fontface = "bold", colour = "black") +
    facet_wrap(~ zone, ncol = 3, scales = "free_x") +
    scale_colour_manual(values = PALETTE, guide = "none") +
    scale_x_continuous(name   = x_label,
                       expand = expansion(mult = c(0.05, 0.18))) +
    scale_y_discrete(name = NULL) +
    labs(title    = model_name,
         subtitle = "emmeans contrasts vs reference  |  BH-FDR per treatment  |  95% CI") +
    theme_bw(base_size = 9) +
    theme(plot.title       = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold"))
}

p3 <- make_forest_emm(
        res_raw,
        expression(beta ~ "(frac enwrapped, raw scale)"),
        "frac_enwrapped ~ treatment * zone + (1|mouse_id)  [primary]") /
      make_forest_emm(
        res_logit,
        expression(beta ~ "(logit scale)"),
        "logit(frac_enwrapped) ~ treatment * zone + (1|mouse_id)  [sensitivity]") +
  plot_annotation(
    title    = "Zone-stratified treatment effects on PV enwrapment",
    subtitle = paste0("Single interaction model — pooled variance across all zones\n",
                      "Top: raw  |  Bottom: logit  |  Consistent results = robust finding"),
    theme    = theme(plot.title    = element_text(face = "bold"),
                     plot.subtitle = element_text(colour = "grey40", size = 8))
  )

# =============================================================================
# 7. FIGURE 4 — Animal dots + emmeans predicted CI
# =============================================================================
animal_zone <- zone_df[, .(frac_enwrapped = mean(frac_enwrapped)),
                         by = .(treatment, mouse_id, zone)]
animal_zone[, treatment := factor(treatment, levels = TREAT_ORDER)]
animal_zone[, zone      := factor(zone, levels = ZONE_ORDER)]

emm_pred <- as.data.table(summary(emmeans(fit_raw, ~ treatment | zone)))
setnames(emm_pred,
  c("emmean", "lower.CL", "upper.CL"),
  c("pred_mean", "pred_lo",  "pred_hi"),
  skip_absent = TRUE)
emm_pred[, treatment := factor(treatment, levels = TREAT_ORDER)]
emm_pred[, zone      := factor(zone, levels = ZONE_ORDER)]

p4 <- ggplot(animal_zone,
             aes(x = treatment, y = frac_enwrapped, colour = treatment)) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
  geom_errorbar(data = emm_pred,
                aes(y = pred_mean, ymin = pred_lo, ymax = pred_hi),
                colour = "black", width = 0.3, linewidth = 1) +
  geom_point(data = emm_pred, aes(y = pred_mean),
             colour = "black", size = 4, shape = 18) +
  facet_wrap(~ zone, ncol = 3) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Fraction PV cells enwrapped (ipsilateral)") +
  labs(
    title    = "Animal-level enwrapment per zone",
    subtitle = "Points = individual animals  |  Diamond + bar = emmeans predicted mean ± 95% CI"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x     = element_text(angle = 35, hjust = 1),
        plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"))

# =============================================================================
# 8. SAVE
# =============================================================================
for (ext in c("pdf", "png")) {
  dpi <- 300
  ggsave(file.path(OUT_DIR, paste0("fig3_zone_forest.", ext)), p3,
         width = 11, height = 7,   dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig4_zone_dots.",   ext)), p4,
         width = 10, height = 4.5, dpi = dpi)
}

cat(sprintf("\n── All outputs → %s\n", OUT_DIR))
cat("   zone_interaction_anova.csv  — F-test for treatment:zone interaction\n")
cat("   zone_contrasts_raw.csv      — emmeans, raw scale\n")
cat("   zone_contrasts_logit.csv    — emmeans, logit scale (+ pct_change)\n")
cat("   residual_check.pdf\n")
cat("   fig3_zone_forest.pdf\n")
cat("   fig4_zone_dots.pdf\n")
