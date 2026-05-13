# =============================================================================
# 05_sensitivity.R
# =============================================================================
# Repeats the primary model (script 01), zone model (script 02), and layer
# model (script 03) with OUTLIER excluded from all datasets.  The excluded
# animal is specified in the configuration block below.
#
# INPUT:  results/00_enwrapment/slice_enwrapment.csv
#         results/00_enwrapment/zone_enwrapment.csv
#         results/00_enwrapment/atlas_enwrapment.csv
#         results/01_primary_LMM/lmm_log_results.csv
#         results/02_zone_gradient/zone_contrasts_raw.csv
#         results/03_atlas_stratified/layer_contrasts.csv
# OUTPUT: results/05_sensitivity/
#   sensitivity_primary.csv    — script 01 model, outlier excluded
#   sensitivity_zone.csv       — script 02 zone contrasts, outlier excluded
#   sensitivity_layer.csv      — script 03 layer contrasts, outlier excluded
#   comparison_primary.csv     — primary vs sensitivity, side by side
#   comparison_zone.csv
#   comparison_layer.csv
#   fig11_sensitivity.pdf/png  — β comparison primary vs sensitivity
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
OUT_DIR     <- file.path(RESULTS_DIR, "05_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Sensitivity exclusion ─────────────────────────────────────────────────────
# Set OUTLIER to the animal ID to exclude (e.g. confirmed injection failure).
# Set to NULL to run without any exclusion.
OUTLIER <- NULL   # e.g. "GROUP_ANIMAL_4"

PALETTE <- c(
  ADAMTS4        = "#4e9af1",
  ADAMTS4_MD     = "#f17c4e",
  ADAMTS15       = "#4ef196",
  C6ST1          = "#c44ef1",
  C6ST1_ADAMTS15 = "#f1c44e",
  mScarlet       = "#888888"
)
TREAT_ORDER       <- c("mScarlet","C6ST1","ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1_ADAMTS15")
ZONE_ORDER        <- c("Core","Penumbra","Outside")
LAYER_ORDER_MODEL <- c("layer 1","layer 2/3","layer 4","layer 5","layer 6")
PLOT_ORDER        <- c("C6ST1_ADAMTS15","ADAMTS15","ADAMTS4_MD","ADAMTS4","C6ST1")
LOGIT_OFFSET      <- 0.01

# =============================================================================
# 1. LOAD, EXCLUDE OUTLIER, TRANSFORM
# =============================================================================
slice_df <- fread(file.path(ENW_DIR, "slice_enwrapment.csv"))
zone_df  <- fread(file.path(ENW_DIR, "zone_enwrapment.csv"))
atlas_df <- fread(file.path(ENW_DIR, "atlas_enwrapment.csv"))

for (df in list(slice_df, zone_df, atlas_df))
  df[, treatment := factor(treatment, levels = TREAT_ORDER)]

if (!is.null(OUTLIER)) {
  n_removed <- c(
    slice = nrow(slice_df[mouse_id == OUTLIER]),
    zone  = nrow(zone_df[mouse_id  == OUTLIER]),
    atlas = nrow(atlas_df[mouse_id == OUTLIER])
  )
  cat(sprintf("Excluding %s:\n", OUTLIER))
  cat(sprintf("  slice: %d rows  |  zone: %d rows  |  atlas: %d rows\n",
              n_removed["slice"], n_removed["zone"], n_removed["atlas"]))
  slice_df <- slice_df[mouse_id != OUTLIER]
  zone_df  <- zone_df[mouse_id  != OUTLIER]
  atlas_df <- atlas_df[mouse_id != OUTLIER]
} else {
  cat("No outlier exclusion specified — running on full dataset.\n")
}

sens_slice <- slice_df
sens_zone  <- zone_df
sens_atlas <- atlas_df

sens_slice[ratio > 0, log_ratio := log(ratio)]
sens_slice_log <- sens_slice[ratio > 0]

sens_zone[, logit_fe := log(
  (frac_enwrapped + LOGIT_OFFSET) / (1 - frac_enwrapped + LOGIT_OFFSET)
)]
sens_zone[, zone := factor(zone, levels = ZONE_ORDER)]

sens_atlas[ratio > 0, log_ratio := log(ratio)]
sens_atlas_log <- sens_atlas[ratio > 0]
sens_atlas_log[, layer_model := fcase(
  layer %in% c("layer 6a","layer 6b"), "layer 6",
  default = as.character(layer)
)]
sens_atlas_log[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
check_singular <- function(fit, label) {
  if (isSingular(fit))
    cat(sprintf("  *** SINGULAR [%s] ***\n", label))
  else
    cat(sprintf("  %s: no singularity.\n", label))
}

extract_coefs <- function(fit, strat_var = NULL, model_label = "") {
  if (is.null(strat_var)) {
    coefs <- as.data.frame(coef(summary(fit, ddf = "Satterthwaite")))
    coefs$term <- rownames(coefs)
    setDT(coefs)
    res <- coefs[term != "(Intercept)"]
    setnames(res,
      c("Estimate","Std. Error","df","t value","Pr(>|t|)"),
      c("beta","se","df","t","p"))
    res[, treatment := gsub("^treatment","",term)]
    res[, ci_lo     := beta - qt(0.975, df) * se]
    res[, ci_hi     := beta + qt(0.975, df) * se]
  } else {
    f   <- as.formula(paste("~ treatment |", strat_var))
    emm <- emmeans(fit, f)
    con <- contrast(emm, method = "trt.vs.ctrl", ref = "mScarlet")
    res <- as.data.table(summary(con, infer = TRUE))
    old <- intersect(c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"), names(res))
    new <- c("beta","se","df","ci_lo","ci_hi","t","p")[
             match(old, c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"))]
    setnames(res, old, new)
    res[, treatment := gsub(" - mScarlet","", contrast)]
  }
  res[, fold_change := exp(beta)]
  res[, pct_change  := (fold_change - 1) * 100]
  by_col <- if (!is.null(strat_var)) strat_var else "treatment"
  res[, p_adj := p.adjust(p, method = "BH"), by = eval(by_col)]
  res[, sig := fcase(
    p_adj < 0.001, "***",
    p_adj < 0.01,  "**",
    p_adj < 0.05,  "*",
    p_adj < 0.10,  ".",
    default =      "ns"
  )]
  res[, model := model_label]
  res
}

# =============================================================================
# 3. SENSITIVITY PRIMARY
# =============================================================================
cat("── Sensitivity primary: log(ratio) ~ treatment + (1|mouse_id) ───────────\n")
fit_prim <- lmer(
  log_ratio ~ treatment + (1 | mouse_id),
  data    = sens_slice_log,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
check_singular(fit_prim, "primary sensitivity")

res_prim <- extract_coefs(fit_prim, strat_var = NULL, model_label = "sensitivity")
res_prim[, p_adj := p.adjust(p, method = "BH")]
res_prim[, sig := fcase(
  p_adj < 0.001, "***", p_adj < 0.01, "**",
  p_adj < 0.05, "*",   p_adj < 0.10, ".",
  default = "ns"
)]

fwrite(res_prim[, .(treatment, beta, se, ci_lo, ci_hi, t, df, p, p_adj, sig,
                    fold_change, pct_change, model)],
       file.path(OUT_DIR, "sensitivity_primary.csv"))

cat("\n=== Primary sensitivity results ===\n")
print(res_prim[order(p), .(treatment, beta, pct_change, p, p_adj, sig)])

# =============================================================================
# 4. SENSITIVITY ZONE
# =============================================================================
cat("\n── Sensitivity zone: frac ~ treatment * zone + (1|mouse_id) ─────────────\n")
fit_zone <- lmer(
  logit_fe ~ treatment * zone + (1 | mouse_id),
  data    = sens_zone,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
check_singular(fit_zone, "zone sensitivity")

res_zone <- extract_coefs(fit_zone, strat_var = "zone", model_label = "sensitivity")
fwrite(res_zone[, .(zone, treatment, beta, pct_change, ci_lo, ci_hi, p, p_adj, sig,
                    fold_change, model)],
       file.path(OUT_DIR, "sensitivity_zone.csv"))

cat("\n=== Zone sensitivity results ===\n")
print(res_zone[order(zone, p), .(zone, treatment, beta, pct_change, p, p_adj, sig)])

# =============================================================================
# 5. SENSITIVITY LAYER
# =============================================================================
cat("\n── Sensitivity layer: log(ratio) ~ treatment * layer + (1|mouse_id) ─────\n")
fit_layer <- lmer(
  log_ratio ~ treatment * layer_model + (1 | mouse_id),
  data    = sens_atlas_log,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
check_singular(fit_layer, "layer sensitivity")

res_layer <- extract_coefs(fit_layer, strat_var = "layer_model",
                           model_label = "sensitivity")
fwrite(res_layer[, .(layer_model, treatment, beta, pct_change, ci_lo, ci_hi,
                     p, p_adj, sig, fold_change, model)],
       file.path(OUT_DIR, "sensitivity_layer.csv"))

cat("\n=== Layer sensitivity results ===\n")
res_layer[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]
print(res_layer[order(layer_model, p),
                .(layer_model, treatment, beta, pct_change, p, p_adj, sig)])

# =============================================================================
# 6. COMPARISON TABLES — primary vs sensitivity
# =============================================================================
cat("\n── Loading primary results for comparison ────────────────────────────────\n")

primary_prim  <- fread(file.path(RESULTS_DIR, "01_primary_LMM",    "lmm_log_results.csv"))
primary_zone  <- fread(file.path(RESULTS_DIR, "02_zone_gradient",  "zone_contrasts_raw.csv"))
primary_layer <- fread(file.path(RESULTS_DIR, "03_atlas_stratified","layer_contrasts.csv"))

make_comparison <- function(prim_dt, sens_dt, join_cols) {
  prim_dt <- copy(prim_dt)
  sens_dt <- copy(sens_dt)

  for (dt in list(prim_dt, sens_dt)) {
    if (!"p_adj" %in% names(dt) && "p" %in% names(dt)) {
      dt[, p_adj := p]
      cat("  NOTE: p_adj missing — using raw p as fallback.\n")
    }
    if (!"sig" %in% names(dt) && "p_adj" %in% names(dt)) {
      dt[, sig := fcase(
        p_adj < 0.001, "***", p_adj < 0.01, "**",
        p_adj < 0.05,  "*",   p_adj < 0.10, ".",
        default = "ns")]
    }
    if (!"pct_change" %in% names(dt) && "beta" %in% names(dt))
      dt[, pct_change := (exp(beta) - 1) * 100]
  }

  base_cols <- c("beta", "pct_change", "p", "p_adj", "sig")
  prim_cols <- c(join_cols, intersect(base_cols, names(prim_dt)))
  sens_cols <- c(join_cols, intersect(base_cols, names(sens_dt)))

  comp <- merge(
    prim_dt[, prim_cols, with = FALSE],
    sens_dt[, sens_cols, with = FALSE],
    by = join_cols, suffixes = c("_conservative", "_sensitivity")
  )

  comp[, robust := fcase(
    sig_conservative == "ns" & sig_sensitivity == "ns",  "ns→ns (stable null)",
    sig_conservative != "ns" & sig_sensitivity != "ns",  "sig→sig (robust)",
    sig_conservative != "ns" & sig_sensitivity == "ns",  "conservative sig, sensitivity ns — flag",
    sig_conservative == "ns" & sig_sensitivity != "ns",  "unmasked by exclusion",
    default = "changed"
  )]
  comp[, beta_shift := beta_sensitivity - beta_conservative]
  comp
}

comp_prim  <- make_comparison(primary_prim,  res_prim,  "treatment")
comp_zone  <- make_comparison(primary_zone,  res_zone,  c("zone","treatment"))
comp_layer <- make_comparison(primary_layer, res_layer, c("layer_model","treatment"))

fwrite(comp_prim,  file.path(OUT_DIR, "comparison_primary.csv"))
fwrite(comp_zone,  file.path(OUT_DIR, "comparison_zone.csv"))
fwrite(comp_layer, file.path(OUT_DIR, "comparison_layer.csv"))

cat("\n=== ROBUSTNESS SUMMARY ===\n\n")
cat("── Primary model ────────────────────────────────────────────────────────\n")
print(comp_prim[order(p_sensitivity),
                .(treatment, beta_conservative, sig_conservative,
                  beta_sensitivity, sig_sensitivity, beta_shift, robust)])

cat("\n── Zone model ───────────────────────────────────────────────────────────\n")
print(comp_zone[order(zone, p_sensitivity),
                .(zone, treatment, beta_conservative, sig_conservative,
                  beta_sensitivity, sig_sensitivity, beta_shift, robust)])

cat("\n── Layer model ──────────────────────────────────────────────────────────\n")
comp_layer[, layer_model := factor(layer_model, levels = LAYER_ORDER_MODEL)]
print(comp_layer[order(layer_model, p_sensitivity),
                 .(layer_model, treatment, beta_conservative, sig_conservative,
                   beta_sensitivity, sig_sensitivity, beta_shift, robust)])

# Flag any significance changes
changes <- rbindlist(list(
  comp_prim[robust %like% "unmasked|flag",
            .(context = "primary", stratum = treatment, robust)],
  comp_zone[robust %like% "unmasked|flag",
            .(context = paste("zone:", zone), stratum = treatment, robust)],
  comp_layer[robust %like% "unmasked|flag",
             .(context = paste("layer:", layer_model), stratum = treatment, robust)]
), fill = TRUE)

if (nrow(changes) > 0) {
  cat("\n── Significance changes (sensitivity vs conservative) ───────────────────\n")
  print(changes[, .(context, stratum, robust)])
} else {
  cat("\n✓ All findings consistent across conservative and sensitivity analyses.\n")
}

# =============================================================================
# 7. FIGURE 11 — β comparison
# =============================================================================
make_comp_plot <- function(comp_dt, stratum_col, facet_col = NULL, title) {
  d <- copy(comp_dt[treatment %in% PLOT_ORDER])
  d[, treatment := factor(treatment, levels = PLOT_ORDER)]

  long <- rbindlist(list(
    d[, .(treatment, stratum = get(stratum_col), beta = beta_sensitivity,
          model = "Sensitivity", sig = sig_sensitivity)],
    d[, .(treatment, stratum = get(stratum_col), beta = beta_conservative,
          model = "Conservative", sig = sig_conservative)]
  ))
  long[, model := factor(model, levels = c("Conservative", "Sensitivity"))]

  p <- ggplot(long, aes(x = beta, y = treatment,
                        colour = treatment, alpha = model, shape = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3.5, position = position_dodge(width = 0.6)) +
    geom_text(data = d[robust %like% "unmasked"],
              aes(x = beta_sensitivity, y = treatment, label = "!"),
              hjust = -0.8, size = 5, fontface = "bold", colour = "red",
              inherit.aes = FALSE) +
    scale_colour_manual(values = PALETTE, guide = "none") +
    scale_alpha_manual(values = c("Conservative" = 0.4, "Sensitivity" = 1.0),
                       name = "Analysis") +
    scale_shape_manual(values = c("Conservative" = 16, "Sensitivity" = 18),
                       name = "Analysis") +
    labs(title = title, x = expression(beta ~ "(log scale)"), y = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  if (!is.null(facet_col))
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), scales = "free_x")
  p
}

pA <- make_comp_plot(comp_prim,  "treatment",  NULL,          "Global")
pB <- make_comp_plot(comp_zone,  "treatment",  "zone",        "Zone")
pC <- make_comp_plot(comp_layer, "treatment",  "layer_model", "Layer")

p11 <- (pA / pB / pC) +
  plot_layout(heights = c(1, 1, 2)) +
  plot_annotation(
    title    = "Sensitivity analysis: impact of outlier exclusion",
    subtitle = "Solid diamonds = sensitivity; faded circles = conservative (full dataset)\nRed '!' = effects unmasked by exclusion.",
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

for (ext in c("pdf","png")) {
  dpi <- 300
  ggsave(file.path(OUT_DIR, paste0("fig11_sensitivity.", ext)), p11,
         width = 12, height = 14, dpi = dpi)
}

cat(sprintf("\n── All outputs → %s\n", OUT_DIR))
cat("   sensitivity_primary/zone/layer.csv — outlier-excluded results\n")
cat("   comparison_primary/zone/layer.csv  — side-by-side comparisons\n")
cat("   fig11_sensitivity.pdf              — β comparison with robustness flags\n")
