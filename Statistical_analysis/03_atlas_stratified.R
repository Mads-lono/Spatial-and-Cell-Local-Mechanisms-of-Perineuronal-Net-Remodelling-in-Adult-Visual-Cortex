# =============================================================================
# 03_atlas_stratified.R
# =============================================================================
# INPUT:  results/00_enwrapment/atlas_enwrapment.csv
# OUTPUT: results/03_atlas_stratified/
#   layer_interaction_anova.csv    — Type III ANOVA: treatment*layer
#   layer_contrasts.csv            — emmeans contrasts vs reference, BH-FDR
#   area_interaction_anova.csv     — Type III ANOVA: treatment*area_type
#   area_contrasts.csv             — emmeans contrasts vs reference, BH-FDR
#   fig5_layer_forest.pdf/png
#   fig6_area_forest.pdf/png
#   fig7_heatmap.pdf/png           — exploratory: treatment × area × layer
#   residual_check.pdf/png
#
# OUTCOME: log(ratio) — consistent with script 01.
#   ratio = frac_enwrapped_ipsi / frac_enwrapped_contra per animal per area×layer.
#   Rows with ratio <= 0 are excluded (log undefined).
#   Log scale: β = 0 → no change; β = -0.3 → ~26% reduction.
#
# MODEL 1 — Layer stratification (primary):
#   log(ratio) ~ treatment * layer + (1 | mouse_id)
#   Layer 6a and 6b are merged into "layer 6" to ensure sufficient cell counts
#   across groups. Original layers are retained for the exploratory heatmap.
#   emmeans contrasts vs reference within each layer, BH-FDR per treatment.
#
# MODEL 2 — Area type (primary vs secondary visual cortex):
#   log(ratio) ~ treatment * area_type + (1 | mouse_id)
#   Sparse areas (< 2 animals in any group) are excluded from this model.
#   emmeans contrasts vs reference within each area_type, BH-FDR per treatment.
#
# HEATMAP (exploratory, no per-cell statistics):
#   Full treatment × area × layer grid. Too sparse for per-cell LMM in many
#   cells. Shows raw mean log(ratio) per cell for orientation only.
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
OUT_DIR     <- file.path(RESULTS_DIR, "03_atlas_stratified")
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

LAYER_ORDER       <- c("layer 1", "layer 2/3", "layer 4",
                        "layer 5", "layer 6a",  "layer 6b")   # heatmap only
LAYER_ORDER_MODEL <- c("layer 1", "layer 2/3", "layer 4",
                        "layer 5", "layer 6")                  # primary model

# Adjust area lists to match your atlas; sparse areas are excluded from models
PRIMARY_AREA   <- "VISp"        # primary target area short name
SECONDARY_AREAS <- c("VISal", "VISam", "VISl", "VISpm")
EXCLUDED_AREAS  <- c("VISpl")   # too sparse for model — edit as needed
AREA_ORDER      <- c(PRIMARY_AREA, SECONDARY_AREAS)

PLOT_ORDER <- c("C6ST1_ADAMTS15", "ADAMTS15", "ADAMTS4_MD", "ADAMTS4", "C6ST1")

# =============================================================================
# 1. LOAD & PREPARE
# =============================================================================
atlas <- fread(file.path(ENW_DIR, "atlas_enwrapment.csv"))
atlas[, treatment := factor(treatment, levels = TREAT_ORDER)]
atlas[, layer     := factor(layer,     levels = LAYER_ORDER)]

n_zero <- sum(atlas$ratio <= 0, na.rm = TRUE)
cat(sprintf("Atlas records: %d  |  ratio <= 0 excluded: %d\n",
            nrow(atlas), n_zero))
if (n_zero > 0)
  print(atlas[ratio <= 0, .(treatment, mouse_id, area_short, layer, ratio)])

atlas_log <- atlas[ratio > 0]
atlas_log[, log_ratio := log(ratio)]

# Area type — primary vs secondary; excluded areas dropped from area model
atlas_log[, area_type := fcase(
  area_short == PRIMARY_AREA,              "primary",
  area_short %in% EXCLUDED_AREAS,          NA_character_,
  area_short %in% SECONDARY_AREAS,         "secondary",
  default = NA_character_
)]
atlas_area <- atlas_log[!is.na(area_type)]
atlas_area[, area_type := factor(area_type, levels = c("secondary", "primary"))]

# Merge thin layers into "layer 6" for the primary model.
# Any layer with fewer than the minimum animals per group in the reference
# treatment should be merged with an adjacent layer.
atlas_log[, layer_model := fcase(
  layer %in% c("layer 6a", "layer 6b"), "layer 6",
  default = as.character(layer)
)]
atlas_log[,  layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]
atlas_area[, layer_model := fcase(
  layer %in% c("layer 6a", "layer 6b"), "layer 6",
  default = as.character(layer)
)]
atlas_area[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]

cat(sprintf("Layer model: %d records  |  Area model: %d records\n",
            nrow(atlas_log), nrow(atlas_area)))

# Flag any group × layer cell with fewer than 3 animals — contrasts will be fragile
sparse_check <- atlas_log[, .(n = uniqueN(mouse_id)),
                            by = .(treatment, layer_model)]
sparse_cells <- sparse_check[n < 3]
if (nrow(sparse_cells) > 0) {
  cat("  NOTE: sparse group × layer cells (n < 3 animals) — contrasts will have low df:\n")
  print(sparse_cells)
}

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
check_singular <- function(fit, label) {
  if (isSingular(fit))
    cat(sprintf("  *** SINGULAR FIT [%s] ***\n", label))
  else
    cat(sprintf("  %s: no singularity.\n", label))
}

report_icc <- function(fit, label) {
  vc      <- as.data.frame(VarCorr(fit))
  between <- vc[vc$grp == "mouse_id", "vcov"]
  icc     <- between / sum(vc$vcov)
  interp  <- if (icc > 0.5) "animal differences dominate"
             else if (icc > 0.2) "moderate between-animal variance"
             else if (icc > 0.1) "low-moderate between-animal variance"
             else "very consistent across animals"
  cat(sprintf("  ICC %-18s: %.3f (%.0f%%) — %s\n",
              label, icc, icc * 100, interp))
}

extract_contrasts <- function(fit, strat_var, model_label) {
  f   <- as.formula(paste("~ treatment |", strat_var))
  emm <- emmeans(fit, f)
  con <- contrast(emm, method = "trt.vs.ctrl", ref = "mScarlet")
  res <- as.data.table(summary(con, infer = TRUE))

  old <- intersect(c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"),
                   names(res))
  new <- c("beta","se","df","ci_lo","ci_hi","t","p")[
           match(old, c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"))]
  setnames(res, old, new)

  con_col <- if ("contrast" %in% names(res)) "contrast" else names(res)[1]
  res[, treatment := gsub(" - mScarlet", "", get(con_col))]

  res[, p_adj := p.adjust(p, method = "BH"), by = treatment]
  res[, sig := fcase(
    p_adj < 0.001, "***",
    p_adj < 0.01,  "**",
    p_adj < 0.05,  "*",
    p_adj < 0.10,  ".",
    default =      "ns"
  )]
  res[, fold_change := exp(beta)]
  res[, pct_change  := (fold_change - 1) * 100]
  res[, model := model_label]

  cols <- c("model", strat_var, "treatment",
            "beta", "se", "ci_lo", "ci_hi", "t", "df",
            "p", "p_adj", "sig", "fold_change", "pct_change")
  res[, .SD, .SDcols = intersect(cols, names(res))]
}

# =============================================================================
# 3. MODEL 1 — LAYER
# =============================================================================
cat("\n── Model 1: log(ratio) ~ treatment * layer + (1|mouse_id) ──────────────\n")

fit_layer <- lmer(
  log_ratio ~ treatment * layer_model + (1 | mouse_id),
  data    = atlas_log,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
check_singular(fit_layer, "layer model")
report_icc(fit_layer,     "layer model")

cat("\n── Type III ANOVA — layer model ─────────────────────────────────────────\n")
anova_layer <- as.data.frame(anova(fit_layer, ddf = "Satterthwaite"))
anova_layer$term <- rownames(anova_layer)
print(anova_layer[, c("term", "NumDF", "DenDF", "F value", "Pr(>F)")])
fwrite(setDT(anova_layer), file.path(OUT_DIR, "layer_interaction_anova.csv"))

layer_ix_p <- anova_layer[anova_layer$term == "treatment:layer_model", "Pr(>F)"]
if (length(layer_ix_p) > 0) {
  if (layer_ix_p < 0.05) {
    cat(sprintf("\n  *** treatment:layer_model p=%.4f SIGNIFICANT ***\n", layer_ix_p))
    cat("  Layer contrasts are PRIMARY. Report as headline result.\n")
  } else {
    cat(sprintf("\n  treatment:layer_model p=%.4f — not significant.\n", layer_ix_p))
    cat("  Report treatment MAIN EFFECT only. Layer contrasts are EXPLORATORY.\n")
  }
}

cat("\n── emmeans contrasts vs reference per layer ──────────────────────────────\n")
res_layer <- extract_contrasts(fit_layer, "layer_model", "layer_model")
fwrite(res_layer, file.path(OUT_DIR, "layer_contrasts.csv"))

fragile <- res_layer[df < 3]
if (nrow(fragile) > 0) {
  cat("\n  *** FRAGILE CONTRASTS (df < 3) — do not report as primary results ***\n")
  print(fragile[, .(layer_model, treatment, df, p, sig)])
}

print(res_layer[order(layer_model, treatment),
                .(layer_model, treatment, beta, pct_change, p, p_adj, sig)])

# =============================================================================
# 4. MODEL 2 — AREA TYPE
# =============================================================================
cat("\n── Model 2: log(ratio) ~ treatment * area_type + (1|mouse_id) ───────────\n")

fit_area <- lmer(
  log_ratio ~ treatment * area_type + (1 | mouse_id),
  data    = atlas_area,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
check_singular(fit_area, "area_type model")
report_icc(fit_area,     "area_type model")

cat("\n── Type III ANOVA — area model ──────────────────────────────────────────\n")
anova_area <- as.data.frame(anova(fit_area, ddf = "Satterthwaite"))
anova_area$term <- rownames(anova_area)
print(anova_area[, c("term", "NumDF", "DenDF", "F value", "Pr(>F)")])
fwrite(setDT(anova_area), file.path(OUT_DIR, "area_interaction_anova.csv"))

area_ix_p <- anova_area[anova_area$term == "treatment:area_type", "Pr(>F)"]
if (length(area_ix_p) > 0) {
  if (area_ix_p < 0.05) {
    cat(sprintf("\n  *** treatment:area_type p=%.4f — primary and secondary areas differ ***\n", area_ix_p))
  } else {
    cat(sprintf("\n  treatment:area_type p=%.4f — no primary vs secondary difference.\n", area_ix_p))
    cat("  Area contrasts are EXPLORATORY.\n")
  }
}

cat("\n── emmeans contrasts vs reference per area type ───────────────────────────\n")
res_area <- extract_contrasts(fit_area, "area_type", "area_model")
fwrite(res_area, file.path(OUT_DIR, "area_contrasts.csv"))
print(res_area[order(area_type, treatment),
               .(area_type, treatment, beta, pct_change, p, p_adj, sig)])

# =============================================================================
# 5. RESIDUAL DIAGNOSTICS
# =============================================================================
diag_panel <- function(fit, title) {
  df <- data.frame(
    fitted    = fitted(fit),
    std_resid = as.numeric(scale(resid(fit, type = "pearson")))
  )
  p_qq <- ggplot(df, aes(sample = std_resid)) +
    stat_qq(size = 0.8, alpha = 0.4) +
    stat_qq_line(colour = "red", linewidth = 0.6) +
    labs(title = paste("QQ —", title), x = "Theoretical", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_rv <- ggplot(df, aes(x = fitted, y = std_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1, alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, colour = "red",
                linewidth = 0.5, formula = y ~ x) +
    labs(title = paste("Resid —", title), x = "Fitted", y = "Std. residuals") +
    theme_bw(base_size = 9)
  p_qq + p_rv
}

p_diag <- diag_panel(fit_layer, "layer model") /
          diag_panel(fit_area,  "area_type model") +
  plot_annotation(title = "Residual diagnostics — atlas models",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(OUT_DIR, "residual_check.pdf"), p_diag, width=9, height=8, dpi=300)
ggsave(file.path(OUT_DIR, "residual_check.png"), p_diag, width=9, height=8, dpi=150)

# =============================================================================
# 6–8. FIGURES
# =============================================================================
make_forest <- function(res, strat_var, strat_levels, x_label, title, subtitle) {
  d <- copy(res[treatment %in% PLOT_ORDER])
  d[, treatment := factor(treatment, levels = PLOT_ORDER)]
  d[[strat_var]] <- factor(d[[strat_var]], levels = strat_levels)

  d[, x_offset := {
    r <- diff(range(c(ci_lo, ci_hi), na.rm = TRUE))
    rep(r * 0.06, .N)
  }, by = strat_var]

  ggplot(d, aes(x = beta, y = treatment, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, y = treatment),
                  orientation = "y", width = 0.25, linewidth = 0.7) +
    geom_point(size = 3) +
    geom_text(aes(x = ci_hi + x_offset, label = sig),
              hjust = 0, size = 3.5, fontface = "bold", colour = "black") +
    facet_wrap(as.formula(paste("~", strat_var)),
               ncol = length(strat_levels), scales = "free_x") +
    scale_colour_manual(values = PALETTE, guide = "none") +
    scale_x_continuous(name = x_label,
                       expand = expansion(mult = c(0.05, 0.18))) +
    scale_y_discrete(name = NULL) +
    labs(title = title, subtitle = subtitle) +
    theme_bw(base_size = 9) +
    theme(plot.title       = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold", size = 8))
}

p5 <- make_forest(
  res_layer, "layer_model", LAYER_ORDER_MODEL,
  expression(beta ~ "(log ratio scale)"),
  "Layer-stratified treatment effects  |  log(ratio) ~ treatment * layer + (1|mouse_id)",
  "emmeans vs reference  |  BH-FDR per treatment across layers  |  95% CI"
)

p6 <- make_forest(
  res_area, "area_type", c("primary", "secondary"),
  expression(beta ~ "(log ratio scale)"),
  "Primary vs secondary visual areas  |  log(ratio) ~ treatment * area_type + (1|mouse_id)",
  "emmeans vs reference  |  BH-FDR per treatment"
)

hmap_data <- atlas_log[area_short %in% AREA_ORDER,
  .(mean_log_ratio = mean(log_ratio, na.rm = TRUE),
    n_animals = uniqueN(mouse_id)),
  by = .(treatment, area_short, layer)]
hmap_data[, treatment := factor(treatment, levels = TREAT_ORDER)]
hmap_data[, area_short := factor(area_short, levels = AREA_ORDER)]
hmap_data[, layer      := factor(layer,      levels = LAYER_ORDER)]

p7 <- ggplot(hmap_data, aes(x = area_short, y = layer, fill = mean_log_ratio)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", mean_log_ratio)),
            size = 2.2, colour = "grey20") +
  facet_wrap(~ treatment, ncol = 3) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d73027", midpoint = 0,
    name = "mean\nlog(ratio)", limits = c(-1.2, 1.2), oob = scales::squish
  ) +
  scale_x_discrete(name = "Visual area") +
  scale_y_discrete(name = "Layer") +
  labs(
    title    = "Mean log(ratio) per treatment × area × layer  [exploratory]",
    subtitle = "Blue = enwrapment reduced vs contra  |  Red = increased  |  White = no change\nNo per-cell statistics — too sparse for LMM in many cells"
  ) +
  theme_bw(base_size = 9) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(colour = "grey40", size = 7),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_text(angle = 35, hjust = 1),
        legend.position  = "right")

for (ext in c("pdf", "png")) {
  dpi <- 300
  ggsave(file.path(OUT_DIR, paste0("fig5_layer_forest.", ext)), p5,
         width = 13, height = 4.5, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig6_area_forest.",  ext)), p6,
         width = 7,  height = 4.5, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig7_heatmap.",      ext)), p7,
         width = 11, height = 7,   dpi = dpi)
}

cat(sprintf("\n── All outputs → %s\n", OUT_DIR))
cat("   layer_interaction_anova.csv  — F-test: does layer modulate treatment effect?\n")
cat("   layer_contrasts.csv          — emmeans per layer vs reference\n")
cat("   area_interaction_anova.csv   — F-test: primary vs secondary areas\n")
cat("   area_contrasts.csv           — emmeans per area_type vs reference\n")
cat("   residual_check.pdf\n")
cat("   fig5_layer_forest.pdf        — layer panels\n")
cat("   fig6_area_forest.pdf         — area type panels\n")
cat("   fig7_heatmap.pdf             — exploratory area × layer grid\n")
