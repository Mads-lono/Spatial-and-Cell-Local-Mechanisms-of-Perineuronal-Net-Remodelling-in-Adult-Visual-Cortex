# =============================================================================
# 06_MD_enhancement.R
# =============================================================================
# INPUT:  results/00_enwrapment/slice_enwrapment.csv
#         results/00_enwrapment/zone_enwrapment.csv
#         results/00_enwrapment/atlas_enwrapment.csv
#         results/01_primary_LMM/sensitivity_primary.csv  (for Fig 11 context)
# OUTPUT: results/06_MD_enhancement/
#   md_global.csv            — whole-hemisphere comparison
#   md_zone.csv              — zone-stratified: treatment × zone interaction
#   md_layer.csv             — layer-stratified: treatment × layer interaction
#   md_zone_anova.csv        — ANOVA for zone interaction F-test
#   md_layer_anova.csv       — ANOVA for layer interaction F-test
#   fig11_md_global.pdf/png
#   fig12_md_zone.pdf/png
#   fig13_md_layer.pdf/png
#   residual_check_md.pdf/png
#
# QUESTION: Does activity deprivation (e.g. monocular deprivation) provide
#   additive PNN disruption on top of enzymatic cleavage?
#   The two groups compared are set in MD_GROUPS below (first = reference).


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
OUT_DIR     <- file.path(RESULTS_DIR, "06_MD_enhancement")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Groups compared ───────────────────────────────────────────────────────────
# First entry is the reference group; second is the group with activity
# deprivation. Adjust to match your treatment labels.
MD_GROUPS <- c("ADAMTS4", "ADAMTS4_MD")   # ADAMTS4 = enzyme only (reference)

MD_PALETTE <- c(
  ADAMTS4    = "#4e9af1",
  ADAMTS4_MD = "#f17c4e"
)

ZONE_ORDER        <- c("Core", "Penumbra", "Outside")
LAYER_ORDER_MODEL <- c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6")
LOGIT_OFFSET      <- 0.01

cat("── Script 06: Activity-deprivation enhancement ──────────────────────────────\n")
cat(sprintf("   Reference:   %s\n",   MD_GROUPS[1]))
cat(sprintf("   Comparison:  %s\n",   MD_GROUPS[2]))
cat("   β < 0: comparison group shows GREATER enwrapment reduction\n")
cat("   β > 0: reference group shows greater reduction\n\n")

# =============================================================================
# 1. LOAD & SUBSET
# =============================================================================
slice_df <- fread(file.path(ENW_DIR, "slice_enwrapment.csv"))
zone_df  <- fread(file.path(ENW_DIR, "zone_enwrapment.csv"))
atlas_df <- fread(file.path(ENW_DIR, "atlas_enwrapment.csv"))

md_slice <- slice_df[treatment %in% MD_GROUPS]
md_zone  <- zone_df[treatment  %in% MD_GROUPS]
md_atlas <- atlas_df[treatment %in% MD_GROUPS]

md_slice[, treatment := factor(treatment, levels = MD_GROUPS)]
md_zone[,  treatment := factor(treatment, levels = MD_GROUPS)]
md_atlas[, treatment := factor(treatment, levels = MD_GROUPS)]

n_zero_slice <- sum(md_slice$ratio <= 0)
md_slice_log <- md_slice[ratio > 0]
md_slice_log[, log_ratio := log(ratio)]

md_zone[, logit_fe := log(
  (frac_enwrapped + LOGIT_OFFSET) / (1 - frac_enwrapped + LOGIT_OFFSET)
)]
md_zone[, zone := factor(zone, levels = ZONE_ORDER)]

n_zero_atlas <- sum(md_atlas$ratio <= 0)
md_atlas_log <- md_atlas[ratio > 0]
md_atlas_log[, log_ratio   := log(ratio)]
md_atlas_log[, layer_model := fcase(
  layer %in% c("layer 6a", "layer 6b"), "layer 6",
  default = as.character(layer)
)]
md_atlas_log[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]

cat(sprintf("Slice records: %d (log model: %d, %d zero excluded)\n",
            nrow(md_slice), nrow(md_slice_log), n_zero_slice))
cat(sprintf("Zone records:  %d (logit offset=%.2f)\n",
            nrow(md_zone), LOGIT_OFFSET))
cat(sprintf("Atlas records: %d (log model: %d, %d zero excluded)\n\n",
            nrow(md_atlas), nrow(md_atlas_log), n_zero_atlas))

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
fit_md <- function(formula, data, label) {
  # REML=FALSE: required for 2-group LMMs with small n per group
  fit <- lmer(formula, data = data, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa"))
  if (isSingular(fit)) {
    cat(sprintf("  SINGULAR [%s]: between-animal variance → 0. ", label))
    cat("p-values valid; CIs conservative.\n")
  } else {
    cat(sprintf("  %s: no singularity.\n", label))
  }
  fit
}

report_icc <- function(fit, label) {
  vc      <- as.data.frame(VarCorr(fit))
  between <- vc[vc$grp == "mouse_id", "vcov"]
  icc     <- between / sum(vc$vcov)
  interp  <- if (icc > 0.5) "animal differences dominate"
             else if (icc > 0.2) "moderate between-animal variance"
             else if (icc > 0.1) "low-moderate"
             else "very consistent across animals (singular likely)"
  cat(sprintf("  ICC [%-18s]: %.3f (%.0f%%) — %s\n",
              label, icc, icc * 100, interp))
}

extract_md <- function(fit, strat_var = NULL) {
  ref   <- MD_GROUPS[1]
  comp  <- MD_GROUPS[2]
  label <- paste(comp, "-", ref)

  if (is.null(strat_var)) {
    coefs <- as.data.frame(coef(summary(fit, ddf = "Satterthwaite")))
    coefs$term <- rownames(coefs)
    setDT(coefs)
    res <- coefs[term != "(Intercept)"]
    setnames(res,
      c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)"),
      c("beta",     "se",         "df", "t",       "p"))
    res[, contrast := label]
    res[, ci_lo    := beta - qt(0.975, df) * se]
    res[, ci_hi    := beta + qt(0.975, df) * se]
  } else {
    f   <- as.formula(paste("~ treatment |", strat_var))
    emm <- emmeans(fit, f)
    con <- contrast(emm, method = "revpairwise")
    res <- as.data.table(summary(con, infer = TRUE))
    old <- intersect(c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"),
                     names(res))
    new <- c("beta","se","df","ci_lo","ci_hi","t","p")[
             match(old, c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"))]
    setnames(res, old, new)
    res[, contrast := label]
  }

  res[, fold_change := exp(beta)]
  res[, pct_change  := (fold_change - 1) * 100]

  if (!is.null(strat_var) && nrow(res) > 1)
    res[, p_adj := p.adjust(p, method = "BH")]
  else
    res[, p_adj := p]

  res[, sig := fcase(
    p_adj < 0.001, "***",
    p_adj < 0.01,  "**",
    p_adj < 0.05,  "*",
    p_adj < 0.10,  ".",
    default =      "ns"
  )]
  res
}

# =============================================================================
# 3. MODEL A — Global
# =============================================================================
cat("\n── Model A: global — log(ratio) ~ treatment + (1|mouse_id) ─────────────────\n\n")

fit_global <- fit_md(
  log_ratio ~ treatment + (1 | mouse_id),
  md_slice_log, "global"
)
report_icc(fit_global, "global")

res_global <- extract_md(fit_global)
res_global[, model := "global"]
fwrite(res_global, file.path(OUT_DIR, "md_global.csv"))

cat(sprintf(
  "\n=== Global result ===\n  %s vs %s: β=%.3f  pct=%.1f%%  95%%CI[%.3f, %.3f]  p=%.4f  %s\n",
  MD_GROUPS[2], MD_GROUPS[1],
  res_global$beta, res_global$pct_change,
  res_global$ci_lo, res_global$ci_hi,
  res_global$p, res_global$sig
))
cat("  β < 0 → comparison group shows greater reduction\n")

# =============================================================================
# 4. MODEL B — Zone
# =============================================================================
cat("\n── Model B: zone — logit(frac) ~ treatment * zone + (1|mouse_id) ───────────\n\n")

fit_zone <- fit_md(
  logit_fe ~ treatment * zone + (1 | mouse_id),
  md_zone, "zone"
)
report_icc(fit_zone, "zone")

cat("\n── Type III ANOVA — zone model ──────────────────────────────────────────────\n")
anova_zone <- as.data.frame(anova(fit_zone, ddf = "Satterthwaite"))
anova_zone$term <- rownames(anova_zone)
print(anova_zone[, c("term", "NumDF", "DenDF", "F value", "Pr(>F)")])
fwrite(setDT(anova_zone), file.path(OUT_DIR, "md_zone_anova.csv"))

zone_ix_p <- anova_zone[anova_zone$term == "treatment:zone", "Pr(>F)"]
if (length(zone_ix_p) > 0) {
  if (zone_ix_p < 0.05)
    cat(sprintf("\n  *** treatment:zone p=%.4f — SIGNIFICANT ***\n", zone_ix_p))
  else
    cat(sprintf("\n  treatment:zone p=%.4f — not significant.\n", zone_ix_p))
}

res_zone <- extract_md(fit_zone, "zone")
res_zone[, model := "zone"]
fwrite(res_zone, file.path(OUT_DIR, "md_zone.csv"))

cat("\n=== Zone contrasts ===\n")
res_zone[, zone := factor(zone, levels = ZONE_ORDER)]
print(res_zone[order(zone), .(zone, beta, pct_change, ci_lo, ci_hi, p, p_adj, sig)])

# Spatial gradient check
core_b <- res_zone[zone == "Core",    beta]
out_b  <- res_zone[zone == "Outside", beta]
if (length(core_b) > 0 && length(out_b) > 0) {
  cat(sprintf("\n  Core β = %.3f  |  Outside β = %.3f\n", core_b, out_b))
  if (core_b < 0 && out_b > 0)
    cat("  Injection-centred spatial gradient: comparison group stronger at core, reference group stronger distally.\n")
  else if (core_b < 0 && out_b < 0)
    cat("  Comparison group stronger across all zones — no spatial gradient.\n")
  else
    cat("  Mixed spatial pattern.\n")
}

# ── df fragility check ────────────────────────────────────────────────────────
fragile_zone <- res_zone[df < 5]
if (nrow(fragile_zone) > 0) {
  cat("\n  *** LOW df CONTRASTS (df < 5) — treat as indicative only ***\n")
  print(fragile_zone[, .(zone, beta, pct_change, p, p_adj, sig, df)])
}

# =============================================================================
# 5. MODEL C — Layer
# =============================================================================
cat("\n── Model C: layer — log(ratio) ~ treatment * layer_model + (1|mouse_id) ────\n\n")

fit_layer <- fit_md(
  log_ratio ~ treatment * layer_model + (1 | mouse_id),
  md_atlas_log, "layer"
)
report_icc(fit_layer, "layer")

cat("\n── Type III ANOVA — layer model ─────────────────────────────────────────────\n")
anova_layer <- as.data.frame(anova(fit_layer, ddf = "Satterthwaite"))
anova_layer$term <- rownames(anova_layer)
print(anova_layer[, c("term", "NumDF", "DenDF", "F value", "Pr(>F)")])
fwrite(setDT(anova_layer), file.path(OUT_DIR, "md_layer_anova.csv"))

layer_ix_p <- anova_layer[anova_layer$term == "treatment:layer_model", "Pr(>F)"]
if (length(layer_ix_p) > 0) {
  if (layer_ix_p < 0.05)
    cat(sprintf("\n  *** treatment:layer_model p=%.4f — SIGNIFICANT ***\n", layer_ix_p))
  else
    cat(sprintf("\n  treatment:layer_model p=%.4f — not significant.\n", layer_ix_p))
}

res_layer <- extract_md(fit_layer, "layer_model")
res_layer[, model := "layer"]
fwrite(res_layer, file.path(OUT_DIR, "md_layer.csv"))

cat("\n=== Layer contrasts ===\n")
res_layer[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]
print(res_layer[order(layer_model),
                .(layer_model, beta, pct_change, ci_lo, ci_hi, p, p_adj, sig, df)])

fragile_layer <- res_layer[df < 5]
if (nrow(fragile_layer) > 0) {
  cat("\n  *** LOW df CONTRASTS (df < 5) — treat as indicative only ***\n")
  print(fragile_layer[, .(layer_model, beta, pct_change, p, p_adj, sig, df)])
}

# =============================================================================
# 6. RESIDUAL DIAGNOSTICS
# =============================================================================
diag_panel <- function(fit, title) {
  df_d <- data.frame(
    fitted    = fitted(fit),
    std_resid = as.numeric(scale(resid(fit, type = "pearson")))
  )
  p_qq <- ggplot(df_d, aes(sample = std_resid)) +
    stat_qq(size = 1, alpha = 0.5) + stat_qq_line(colour = "red", linewidth = 0.6) +
    labs(title = paste("QQ —", title), x = "Theoretical", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_rv <- ggplot(df_d, aes(x = fitted, y = std_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.9, alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, colour = "red",
                linewidth = 0.5, formula = y ~ x) +
    labs(title = paste("Resid —", title), x = "Fitted", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_qq + p_rv
}

p_diag <- diag_panel(fit_global, "global") /
          diag_panel(fit_zone,   "zone") /
          diag_panel(fit_layer,  "layer") +
  plot_annotation(title = "Residual diagnostics — activity-deprivation models",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(OUT_DIR, "residual_check_md.pdf"), p_diag, width=9, height=11, dpi=300)
ggsave(file.path(OUT_DIR, "residual_check_md.png"), p_diag, width=9, height=11, dpi=150)

# =============================================================================
# 7. FIGURE 11 — Global effect sizes in context
# =============================================================================
# Load individual vs-reference effect sizes from script 01 primary analysis
sens_csv <- file.path(RESULTS_DIR, "05_sensitivity", "sensitivity_primary.csv")
if (file.exists(sens_csv)) {
  script01_ref <- fread(sens_csv)[treatment %in% MD_GROUPS,
                                  .(treatment, beta, ci_lo, ci_hi, pct_change)]
  script01_ref[, source := "vs reference (script 01)"]
  script01_ref[, treatment := factor(treatment, levels = MD_GROUPS)]
} else {
  cat("  NOTE: sensitivity_primary.csv not found — Fig 11 context bars omitted.\n")
  script01_ref <- NULL
}

delta_row <- data.table(
  treatment  = paste0(MD_GROUPS[2], "\nvs ", MD_GROUPS[1]),
  beta       = res_global$beta,
  ci_lo      = res_global$ci_lo,
  ci_hi      = res_global$ci_hi,
  pct_change = res_global$pct_change,
  source     = paste0("Activity contribution (script 06)\np=",
                      round(res_global$p, 3), " ", res_global$sig)
)

if (!is.null(script01_ref)) {
  p11 <- ggplot(script01_ref,
                aes(x = beta, y = treatment, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                  orientation = "y", width = 0.2, linewidth = 0.9) +
    geom_point(size = 4) +
    geom_errorbar(data = delta_row,
                  aes(y = treatment, xmin = ci_lo, xmax = ci_hi, x = beta),
                  colour = "black", width = 0.2, linewidth = 0.9) +
    geom_point(data = delta_row, aes(y = treatment, x = beta),
               colour = "black", size = 4, shape = 18) +
    geom_text(aes(x = ci_hi + 0.03, label = sprintf("%.1f%%", pct_change)),
              hjust = 0, size = 3.2, colour = "grey30") +
    geom_text(data = delta_row,
              aes(x = ci_hi + 0.03, y = treatment,
                  label = sprintf("%.1f%% %s", pct_change, res_global$sig)),
              hjust = 0, size = 3.2, colour = "black") +
    scale_colour_manual(values = MD_PALETTE, guide = "none") +
    scale_x_continuous(name   = expression(beta ~ "(log ratio scale)"),
                       expand = expansion(mult = c(0.05, 0.30))) +
    scale_y_discrete(name = NULL) +
    labs(title    = "Activity contribution: global effect sizes",
         subtitle = paste0(
           "Coloured = individual effects vs reference (script 01)  |  ",
           "Black = comparison − reference (this script)\n",
           "β < 0 = greater enwrapment reduction.")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
} else {
  p11 <- ggplot(delta_row, aes(x = beta, y = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.2) +
    geom_point(size = 4) +
    labs(title = "Activity contribution: direct comparison",
         x = expression(beta ~ "(log ratio scale)"), y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
}

# =============================================================================
# 8. FIGURES 12 & 13 — Zone and layer profiles
# =============================================================================
animal_zone <- md_zone[, .(frac_enwrapped = mean(frac_enwrapped)),
                         by = .(treatment, mouse_id, zone)]
animal_zone[, treatment := factor(treatment, levels = MD_GROUPS)]
animal_zone[, zone      := factor(zone, levels = ZONE_ORDER)]

emm_zone <- as.data.table(summary(emmeans(fit_zone, ~ treatment | zone)))
setnames(emm_zone, c("emmean","lower.CL","upper.CL"),
                   c("pred_logit","pred_lo","pred_hi"), skip_absent = TRUE)
emm_zone[, pred_mean := plogis(pred_logit) - LOGIT_OFFSET]
emm_zone[, pred_lo_p := plogis(pred_lo)    - LOGIT_OFFSET]
emm_zone[, pred_hi_p := plogis(pred_hi)    - LOGIT_OFFSET]
emm_zone[, treatment := factor(treatment, levels = MD_GROUPS)]
emm_zone[, zone      := factor(zone, levels = ZONE_ORDER)]

res_zone_lab <- copy(res_zone)
res_zone_lab[, label := sprintf("Δ=%.0f%% %s", pct_change, sig)]

p12 <- ggplot(animal_zone,
              aes(x = treatment, y = frac_enwrapped, colour = treatment)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.85) +
  geom_errorbar(data = emm_zone,
                aes(y = pred_mean, ymin = pred_lo_p, ymax = pred_hi_p),
                colour = "black", width = 0.25, linewidth = 1) +
  geom_point(data = emm_zone, aes(y = pred_mean),
             colour = "black", size = 4.5, shape = 18) +
  geom_text(data = res_zone_lab,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.4, hjust = 0.5, size = 3.2,
            colour = "black", inherit.aes = FALSE) +
  facet_wrap(~ zone, ncol = 3) +
  scale_colour_manual(values = MD_PALETTE, guide = "none") +
  scale_x_discrete(name = NULL,
                   labels = c(MD_GROUPS[1], paste0(MD_GROUPS[2]))) +
  scale_y_continuous(name = "Fraction PV cells enwrapped (ipsi)") +
  labs(title    = "Zone-stratified activity contribution",
       subtitle = "Δ = % change of comparison group vs reference.") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title  = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"))

animal_layer <- md_atlas_log[, .(log_ratio = mean(log_ratio)),
                               by = .(treatment, mouse_id, layer_model)]
animal_layer[, treatment  := factor(treatment, levels = MD_GROUPS)]
animal_layer[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]

emm_layer <- as.data.table(summary(emmeans(fit_layer, ~ treatment | layer_model)))
setnames(emm_layer, c("emmean","lower.CL","upper.CL"),
                    c("pred_mean","pred_lo","pred_hi"), skip_absent = TRUE)
emm_layer[, treatment  := factor(treatment, levels = MD_GROUPS)]
emm_layer[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]

res_layer_lab <- copy(res_layer)
res_layer_lab[, label := sprintf("Δ=%.0f%% %s", pct_change, sig)]

p13 <- ggplot(animal_layer,
              aes(x = log_ratio, y = layer_model, colour = treatment)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_jitter(height = 0.1, size = 2.5, alpha = 0.8) +
  geom_errorbar(data = emm_layer,
                aes(x = pred_mean, xmin = pred_lo, xmax = pred_hi),
                orientation = "y", width = 0.3, linewidth = 1,
                position = position_dodge(width = 0.5)) +
  geom_point(data = emm_layer, aes(x = pred_mean),
             size = 4, shape = 18,
             position = position_dodge(width = 0.5)) +
  geom_text(data = res_layer_lab,
            aes(x = Inf, y = layer_model, label = label),
            hjust = 1.05, size = 3, colour = "black",
            inherit.aes = FALSE) +
  scale_colour_manual(values = MD_PALETTE, name = NULL) +
  scale_x_continuous(name   = "log(enwrapment ratio)",
                     expand = expansion(mult = c(0.05, 0.22))) +
  scale_y_discrete(name = "Layer") +
  labs(title    = "Layer-stratified activity contribution",
       subtitle = "Points = individual animals  |  Diamond + bar = LMM mean ± 95% CI") +
  theme_bw(base_size = 10) +
  theme(plot.title       = element_text(face = "bold"),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

# =============================================================================
# 9. SAVE ALL
# =============================================================================
for (ext in c("pdf", "png")) {
  dpi <- if (ext == "pdf") 300 else 150
  ggsave(file.path(OUT_DIR, paste0("fig11_md_global.", ext)),  p11,
         width = 8,  height = 4,   dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig12_md_zone.",   ext)),  p12,
         width = 10, height = 4.5, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig13_md_layer.",  ext)),  p13,
         width = 8,  height = 5,   dpi = dpi)
}

cat(sprintf("\n── All outputs → %s\n", OUT_DIR))
cat("   md_global.csv          — whole-hemisphere comparison\n")
cat("   md_zone.csv            — zone contrasts (ipsi-only)\n")
cat("   md_zone_anova.csv      — zone interaction F-test\n")
cat("   md_layer.csv           — layer contrasts\n")
cat("   md_layer_anova.csv     — layer interaction F-test\n")
cat("   fig11_md_global        — effect sizes in context\n")
cat("   fig12_md_zone          — zone profile\n")
cat("   fig13_md_layer         — layer profile\n")
