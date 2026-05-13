# =============================================================================
# 08_off_target_survey.R
# =============================================================================
# RQ8: (A) Is the protease effect specific to visual cortex?
#       (B) How far from the injection site do effects travel?
#
# INPUT:  analysis_results_with_zones/cells_with_zones.csv  — raw cell data
#         results/00_enwrapment/atlas_enwrapment.csv         — visual cortex
# OUTPUT: results/08_off_target_survey/
#   off_target_atlas_enwrapment.csv   — Part 0: computed off-target enwrapment
#   system_lmm_results.csv            — Part A: per-system LMMs (both FDRs)
#   fdr_comparison.csv                — Part A: discrepancies between FDR families
#   visual_vs_other_interaction.csv   — Part B: specificity contrasts
#   off_target_layer_results.csv      — Part C: layer × treatment off-target
#   spread_per_area_betas.csv         — Part D: per-area β for distance plot
#   fig17_system_heatmap.pdf/png
#   fig18_specificity_contrast.pdf/png
#   fig19_distance_gradient.pdf/png
#   fig20_off_target_layer.pdf/png    (only if Part C runs)
#   residual_check_rq8.pdf/png
#
# ANALYSIS PLAN
# ─────────────────────────────────────────────────────────────────────────────
# Part 0 — Compute off-target enwrapment (mirrors script 00 section 6)
#   Identical nn2 nearest-neighbour computation; MIN_PV_ANIMAL = 30 (stricter
#   than visual cortex, lower PV density expected off-target). Excludes visual
#   cortex. Pools slices per animal × area × layer.
#
# Part A — Per-system safety survey
#   log(ratio) ~ treatment + (1|mouse_id), per functional system.
#   FDR reported with BOTH per-system BH and global BH across all
#   system × treatment combinations; discrepancies flagged in output.
#
# Part B — Specificity: is visual cortex more affected than off-target?
#   log(ratio) ~ treatment * is_visual + (1|mouse_id), all areas combined.
#   Primary question. Significant treatment:is_visual = effect localised.
#
# Part C — Layer × treatment in off-target systems (mirrors script 03)
#   log(ratio) ~ treatment * layer_model + (1|mouse_id), per system.
#   Only run where ≥ 3 layers are represented with ≥ 2 animals per cell.
#
# Part D — Spatial spread: how far does the protease travel?
#   Distance tiers defined from Allen CCF topology (VISp = tier 0).
#   Per-area LMM → β per area × treatment. Weighted linear trend:
#   β ~ dist_tier (weights = 1/SE²). Negative slope = effect decays
#   with distance from injection site (localised).
#
# DISTANCE TIERS (Allen CCF topological adjacency, VISp = injection target):
#   0  Primary visual area (VISp)
#   1  Higher visual areas (VISam, VISl, VISal, VISpm, VISpl)
#   2  Directly adjacent — Retrosplenial, Posterior parietal (RSP, PTLp)
#   3  One step further — Auditory, posterior Somatosensory, Subiculum
#   4  Two steps — SS-anterior, Motor, Hippocampus proper
#   5  Distant — Prefrontal/Association, Subcortical
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(RANN)
  library(lme4)
  library(lmerTest)
  if (!requireNamespace("emmeans", quietly = TRUE))
    install.packages("emmeans", repos = "https://cloud.r-project.org")
  library(emmeans)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV   <- "/path/to/analysis_results/cells_with_zones.csv"
RESULTS_DIR <- "/path/to/results"
ENW_DIR     <- file.path(RESULTS_DIR, "00_enwrapment")
OUT_DIR     <- file.path(RESULTS_DIR, "08_off_target_survey")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants (identical to scripts 00–07) ────────────────────────────────────
TREATMENTS      <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS
EXCLUDE_PRIMARY <- "C6ST1_ADAMTS15_4"
COLOC_THRESH    <- 30L     # px — same as script 00
MIN_PV_ANIMAL   <- 30L     # stricter than script 00 (20); lower PV density off-target
LAYER_ORDER     <- c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6")

PALETTE <- c(
  mScarlet       = "#888888",
  ADAMTS4        = "#4e9af1",
  ADAMTS4_MD     = "#f17c4e",
  ADAMTS15       = "#4ef196",
  C6ST1          = "#c44ef1",
  C6ST1_ADAMTS15 = "#f1c44e"
)

SYSTEM_ORDER <- c(
  "Visual (target)",
  "Retrosplenial / Parietal",
  "Auditory",
  "Somatosensory",
  "Motor",
  "Hippocampus",
  "Prefrontal / Association",
  "Other cortex",
  "Subcortical"
)

cat("── Script 08: Off-target survey (RQ8) ──────────────────────────────────────\n")
cat("   Parts: (0) compute off-target enwrapment  (A) system LMMs\n")
cat("          (B) specificity interaction         (C) off-target layers\n")
cat("          (D) spatial spread / distance gradient\n\n")

# ── Helpers ───────────────────────────────────────────────────────────────────
compute_frac_enwrapped <- function(pv_xy, pnn_xy) {
  if (nrow(pv_xy) == 0L || nrow(pnn_xy) == 0L) return(NA_real_)
  nn <- nn2(data = pnn_xy, query = pv_xy, k = 1L)
  mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
}

fit_lmm <- function(formula, data, label, reml = TRUE) {
  fit <- lmer(formula, data = data, REML = reml,
              control = lmerControl(optimizer = "bobyqa"))
  if (isSingular(fit)) cat(sprintf("  SINGULAR [%s]\n", label))
  else                  cat(sprintf("  OK [%s]\n", label))
  fit
}

extract_contrasts <- function(fit, strat_var = NULL, model_label = "") {
  if (is.null(strat_var)) {
    coefs <- as.data.frame(coef(summary(fit, ddf = "Satterthwaite")))
    coefs$term <- rownames(coefs)
    setDT(coefs)
    res <- coefs[term != "(Intercept)"]
    setnames(res,
      c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)"),
      c("beta",     "se",         "df", "t",       "p"))
    res[, treatment := gsub("^treatment", "", term)]
  } else {
    f   <- as.formula(paste("~ treatment |", strat_var))
    emm <- emmeans(fit, f)
    con <- contrast(emm, method = "trt.vs.ctrl", ref = "mScarlet")
    res <- as.data.table(summary(con, infer = TRUE))
    old <- intersect(c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"),
                     names(res))
    new <- c("beta","se","df","ci_lo","ci_hi","t","p")[
             match(old, c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"))]
    setnames(res, old, new)
    res[, treatment := gsub(" - mScarlet", "", contrast)]
  }
  res[, ci_lo       := beta - qt(0.975, df) * se]
  res[, ci_hi       := beta + qt(0.975, df) * se]
  res[, fold_change := exp(beta)]
  res[, pct_change  := (fold_change - 1) * 100]
  res[, model       := model_label]
  res
}

add_fdr <- function(dt, p_col = "p", by_col = NULL, suffix = "_adj") {
  col_sig  <- paste0("sig",   suffix)
  col_padj <- paste0("p_adj", suffix)
  if (!is.null(by_col)) {
    dt[, (col_padj) := p.adjust(get(p_col), method = "BH"), by = by_col]
  } else {
    dt[, (col_padj) := p.adjust(get(p_col), method = "BH")]
  }
  dt[, (col_sig) := fcase(
    get(col_padj) < 0.001, "***",
    get(col_padj) < 0.01,  "**",
    get(col_padj) < 0.05,  "*",
    get(col_padj) < 0.10,  ".",
    default = "ns"
  )]
  dt
}

assign_system <- function(area_vec) {
  sys <- rep("Subcortical", length(area_vec))
  sys[grepl("relimbic|nfralimbic|ingulate|rbitofrontal|ssociation cortex|Frontal|Prefrontal",
             area_vec, ignore.case = TRUE)] <- "Prefrontal / Association"
  sys[grepl("ippocamp|CA[123]|Dentate|ubiculum|ntorhinal",
             area_vec, ignore.case = TRUE)] <- "Hippocampus"
  sys[grepl("Motor|omotor",       area_vec, ignore.case = TRUE)] <- "Motor"
  sys[grepl("omatosensory",       area_vec, ignore.case = TRUE)] <- "Somatosensory"
  sys[grepl("uditory",            area_vec, ignore.case = TRUE)] <- "Auditory"
  sys[grepl("etrosplenial|osterior parietal|PTLp",
             area_vec, ignore.case = TRUE)] <- "Retrosplenial / Parietal"
  sys[sys == "Subcortical" &
      grepl("cortex|Parietal|Temporal|Insular|Ecto|Perirhinal",
            area_vec, ignore.case = TRUE)] <- "Other cortex"
  sys[grepl("isual", area_vec, ignore.case = TRUE)] <- "Visual (target)"
  sys
}

assign_dist_tier <- function(area_vec) {
  tier <- rep(5L, length(area_vec))
  tier[grepl("relimbic|nfralimbic|ingulate|rbitofrontal|ssociation|Frontal",
              area_vec, ignore.case = TRUE)] <- 5L
  tier[grepl("Motor|omotor|omatosensory|ippocamp|CA[123]|Dentate",
              area_vec, ignore.case = TRUE)] <- 4L
  tier[grepl("uditory|ubiculum|ntorhinal",
              area_vec, ignore.case = TRUE)] <- 3L
  tier[grepl("etrosplenial|osterior parietal|PTLp",
              area_vec, ignore.case = TRUE)] <- 2L
  tier[grepl("isual", area_vec, ignore.case = TRUE) &
       !grepl("Primary visual", area_vec)]   <- 1L
  tier[grepl("Primary visual", area_vec)]    <- 0L
  tier
}

# =============================================================================
# PART 0 — COMPUTE OFF-TARGET ENWRAPMENT
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART 0: Computing off-target enwrapment (mirrors script 00, section 6)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

LOAD_COLS <- c("mouse_id", "slice_id", "cell_type", "hemisphere",
               "x_hires", "y_hires", "brain_area", "treatment")

cat("Loading cells_with_zones.csv...\n")
t0 <- proc.time()
dt <- fread(CELLS_CSV, select = LOAD_COLS)
cat(sprintf("  Loaded %s rows in %.1fs\n",
    format(nrow(dt), big.mark = ","), (proc.time() - t0)[["elapsed"]]))

# Exclude visual cortex, keep off-target only
dt <- dt[!str_detect(brain_area, regex("isual", ignore_case = TRUE))]
dt <- dt[treatment %in% TREATMENTS]
dt <- dt[mouse_id  != EXCLUDE_PRIMARY]
cat(sprintf("  After off-target filter: %s rows\n",
    format(nrow(dt), big.mark = ",")))

# Parse area and layer (same format as script 00: "Area name, layer X")
dt[, c("area", "layer_raw") := {
  parts <- str_split_fixed(brain_area, ", ", 2)
  list(str_trim(parts[, 1]), str_to_lower(str_trim(parts[, 2])))
}]
dt[, layer_model := fcase(
  layer_raw %in% c("layer 6a", "layer 6b"), "layer 6",
  layer_raw == "",                            NA_character_,
  default = layer_raw
)]
dt[, layer_model := factor(layer_model, levels = LAYER_ORDER)]

pv  <- dt[cell_type == "PV"]
pnn <- dt[cell_type == "PNN"]
cat(sprintf("  PV: %s  |  PNN: %s\n",
    format(nrow(pv), big.mark = ","), format(nrow(pnn), big.mark = ",")))

# Atlas-level loop: pools all slices per animal × area × layer
atlas_keys_ot <- unique(pv[!is.na(layer_model), .(treatment, mouse_id, area, layer_model)])
cat(sprintf("\n  Processing %d animal × area × layer strata...\n", nrow(atlas_keys_ot)))

t1 <- proc.time()
ot_records <- vector("list", nrow(atlas_keys_ot))

for (i in seq_len(nrow(atlas_keys_ot))) {
  key   <- atlas_keys_ot[i]
  pv_k  <- pv[mouse_id == key$mouse_id & area == key$area &
               layer_model == key$layer_model]
  pnn_k <- pnn[mouse_id == key$mouse_id & area == key$area &
                layer_model == key$layer_model]

  rec <- list(treatment = key$treatment, mouse_id = key$mouse_id,
              area = key$area, layer = as.character(key$layer_model))
  ok  <- TRUE

  for (hemi in c("left", "right")) {
    label <- if (hemi == "left") "ipsi" else "contra"
    pv_h  <- pv_k[hemisphere == hemi, .(x_hires, y_hires)]
    pnn_h <- pnn_k[hemisphere == hemi, .(x_hires, y_hires)]
    n_pv  <- nrow(pv_h)
    if (n_pv < MIN_PV_ANIMAL || nrow(pnn_h) == 0L) { ok <- FALSE; break }
    rec[[paste0("n_pv_",  label)]] <- n_pv
    rec[[paste0("frac_",  label)]] <- compute_frac_enwrapped(pv_h, pnn_h)
  }

  if (!ok) { ot_records[[i]] <- NULL; next }

  frac_i <- rec[["frac_ipsi"]]
  frac_c <- rec[["frac_contra"]]
  rec[["ratio"]] <- if (!is.na(frac_i) && !is.na(frac_c) && frac_c > 0)
                      frac_i / frac_c else NA_real_
  ot_records[[i]] <- rec

  if (i %% 1000L == 0L)
    cat(sprintf("  %d / %d (%.0fs)\n", i, nrow(atlas_keys_ot),
        (proc.time() - t1)[["elapsed"]]))
}

ot_df <- rbindlist(Filter(Negate(is.null), ot_records))
ot_df <- ot_df[!is.na(ratio)]
cat(sprintf("  Valid off-target records: %d\n", nrow(ot_df)))
fwrite(ot_df, file.path(OUT_DIR, "off_target_atlas_enwrapment.csv"))
cat("  Saved: off_target_atlas_enwrapment.csv\n\n")
rm(dt, pv, pnn, ot_records); gc()

# ── Annotate ──────────────────────────────────────────────────────────────────
ot_df[, system    := factor(assign_system(area), levels = SYSTEM_ORDER)]
ot_df[, dist_tier := assign_dist_tier(area)]
ot_df[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
ot_df[, log_ratio := log(ratio)]
ot_df <- ot_df[is.finite(log_ratio)]
ot_df[, is_visual  := FALSE]
ot_df[, layer_model := factor(layer, levels = LAYER_ORDER)]

# ── Load visual cortex atlas (Parts B & D need combined dataset) ──────────────
vis_df <- fread(file.path(ENW_DIR, "atlas_enwrapment.csv"))
vis_df <- vis_df[mouse_id != EXCLUDE_PRIMARY & ratio > 0]
vis_df[, log_ratio   := log(ratio)]
vis_df[, is_visual   := TRUE]
vis_df[, system      := "Visual (target)"]
vis_df[, dist_tier   := assign_dist_tier(area)]
vis_df[, treatment   := factor(treatment, levels = TREATMENT_ORDER)]
vis_df[, layer_model := fcase(
  layer %in% c("layer 6a", "layer 6b"), "layer 6",
  default = as.character(layer)
)]
vis_df[, layer_model := factor(layer_model, levels = LAYER_ORDER)]

COMMON_COLS <- c("treatment", "mouse_id", "area", "layer_model",
                 "log_ratio", "system", "is_visual", "dist_tier")
both_df <- rbind(
  vis_df[, ..COMMON_COLS],
  ot_df[,  ..COMMON_COLS],
  fill = TRUE
)
both_df[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("Combined dataset: %d rows (%d visual + %d off-target)\n\n",
    nrow(both_df), nrow(vis_df), nrow(ot_df)))

# Report system counts
sys_counts <- ot_df[, .(n_areas = uniqueN(area), n_records = .N), by = system]
cat("── Off-target system composition ─────────────────────────────────────────────\n")
print(sys_counts[order(match(system, SYSTEM_ORDER))])
cat("\n")

# =============================================================================
# PART A — PER-SYSTEM LMMs  (per-system BH + global BH, flag discrepancies)
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART A: Per-system safety survey\n")
cat("        log(ratio) ~ treatment + (1|mouse_id), BH per-system AND global\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

systems_to_test <- setdiff(SYSTEM_ORDER, "Visual (target)")

all_sys_res <- rbindlist(lapply(systems_to_test, function(sys) {
  d <- ot_df[system == sys]
  n_per_group <- d[, uniqueN(mouse_id), by = treatment]
  min_anim    <- if (nrow(n_per_group) == 0L) 0L else min(n_per_group$V1, na.rm = TRUE)
  if (nrow(d) < 20 || min_anim < 2) {
    cat(sprintf("  SKIP [%s]: %g rows, min %g animals/group\n",
        sys, nrow(d), min_anim))
    return(NULL)
  }
  cat(sprintf("── %s (%d rows, %d areas) ──\n", sys, nrow(d), uniqueN(d$area)))
  fit <- fit_lmm(log_ratio ~ treatment + (1 | mouse_id), d, sys)
  res <- extract_contrasts(fit, model_label = sys)
  res[, system := sys]
  anova_sys <- as.data.table(anova(fit, ddf = "Satterthwaite"), keep.rownames = "term")
  res[, treatment_F_p := anova_sys[term == "treatment", `Pr(>F)`]]
  cat(sprintf("  Treatment F: p = %.4f\n\n", res$treatment_F_p[1]))
  res
}), fill = TRUE)

# Both FDR families
all_sys_res <- add_fdr(all_sys_res, "p", by_col = "system", suffix = "_local")
all_sys_res <- add_fdr(all_sys_res, "p", by_col = NULL,     suffix = "_global")

all_sys_res[, fdr_discrepancy := fcase(
  sig_local != "ns" & sig_global == "ns",
    "local-only (report with caveat)",
  sig_global != "ns" & sig_local == "ns",
    "global-only (unexpected — flag)",
  default = "consistent"
)]

fwrite(all_sys_res[, .(system, treatment, beta, se, ci_lo, ci_hi, pct_change,
                        p, p_adj_local, sig_local, p_adj_global, sig_global,
                        fdr_discrepancy, treatment_F_p)],
       file.path(OUT_DIR, "system_lmm_results.csv"))

discr <- all_sys_res[fdr_discrepancy != "consistent"]
fwrite(discr[, .(system, treatment, p, p_adj_local, sig_local,
                  p_adj_global, sig_global, fdr_discrepancy)],
       file.path(OUT_DIR, "fdr_comparison.csv"))

cat("── Part A summary ──────────────────────────────────────────────────────────\n")
cat(sprintf("  Local FDR hits:  %d\n", nrow(all_sys_res[sig_local  != "ns"])))
cat(sprintf("  Global FDR hits: %d\n", nrow(all_sys_res[sig_global != "ns"])))
if (nrow(discr) > 0) {
  cat("  *** FDR discrepancies:\n")
  print(discr[, .(system, treatment, p, p_adj_local, p_adj_global, fdr_discrepancy)])
} else {
  cat("  ✓ Local and global FDR consistent on all calls.\n")
}
cat("\n")

# =============================================================================
# PART B — SPECIFICITY INTERACTION (PRIMARY QUESTION)
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART B: Specificity — treatment × is_visual  (PRIMARY)\n")
cat("        log(ratio) ~ treatment * is_visual + (1|mouse_id), all areas\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

fit_spec <- fit_lmm(log_ratio ~ treatment * is_visual + (1 | mouse_id),
                    both_df, "specificity")
anova_spec <- as.data.table(anova(fit_spec, ddf = "Satterthwaite"),
                              keep.rownames = "term")
cat("\n  Type III ANOVA:\n")
print(anova_spec[, .(term, NumDF, DenDF, `F value`, `Pr(>F)`)])
ix_p <- anova_spec[term == "treatment:is_visual", `Pr(>F)`]
cat(sprintf("\n  treatment:is_visual F p = %.4f  %s\n\n", ix_p,
    ifelse(ix_p < 0.05, "✓ Specific to visual cortex",
                         "NS — not significantly stronger in visual cortex")))

emm_spec <- emmeans(fit_spec, ~ treatment | is_visual)
con_spec  <- contrast(emm_spec, method = "trt.vs.ctrl", ref = "mScarlet")
res_spec  <- as.data.table(summary(con_spec, infer = TRUE))
old <- intersect(c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"), names(res_spec))
new <- c("beta","se","df","ci_lo","ci_hi","t","p")[
       match(old, c("estimate","SE","df","lower.CL","upper.CL","t.ratio","p.value"))]
setnames(res_spec, old, new)
res_spec <- add_fdr(res_spec, "p", by_col = "is_visual", suffix = "")
res_spec[, pct_change      := (exp(beta) - 1) * 100]
res_spec[, region          := ifelse(is_visual, "Visual (target)", "Off-target (all)")]
res_spec[, treatment_label := gsub(" - mScarlet", "", contrast)]
fwrite(res_spec, file.path(OUT_DIR, "visual_vs_other_interaction.csv"))

cat("  Stratified contrasts:\n")
print(res_spec[order(-is_visual, contrast),
               .(region, treatment_label, beta, pct_change, p, p_adj, sig)])
cat("\n")

# =============================================================================
# PART C — LAYER × TREATMENT IN OFF-TARGET SYSTEMS
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART C: Layer × treatment in off-target systems (mirrors script 03)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

layer_res_list <- lapply(systems_to_test, function(sys) {
  d <- ot_df[system == sys & !is.na(layer_model)]
  n_layers <- uniqueN(d$layer_model)
  n_min_dt <- d[, uniqueN(mouse_id), by = .(treatment, layer_model)]
  n_min    <- if (nrow(n_min_dt) == 0L) 0L else min(n_min_dt$V1, na.rm = TRUE)
  if (n_layers < 3 || nrow(d) < 50 || n_min < 2) {
    cat(sprintf("  SKIP layer [%s]: %g layers, min %g animals/cell\n",
        sys, n_layers, n_min))
    return(NULL)
  }
  cat(sprintf("── Layer model: %s (%d rows, %d layers) ──\n",
      sys, nrow(d), n_layers))
  fit <- fit_lmm(log_ratio ~ treatment * layer_model + (1 | mouse_id), d,
                 paste(sys, "layer"))
  anova_l <- as.data.table(anova(fit, ddf = "Satterthwaite"), keep.rownames = "term")
  ix_p_l  <- anova_l[term == "treatment:layer_model", `Pr(>F)`]
  cat(sprintf("  treatment:layer_model F: p = %.4f\n\n", ix_p_l))
  res <- extract_contrasts(fit, strat_var = "layer_model",
                           model_label = paste(sys, "layer"))
  res[, system        := sys]
  res[, interaction_p := ix_p_l]
  res <- add_fdr(res, "p", by_col = "layer_model", suffix = "")
  res
})
layer_res <- rbindlist(Filter(Negate(is.null), layer_res_list), fill = TRUE)

if (nrow(layer_res) > 0) {
  layer_res[, layer_model := factor(layer_model, levels = LAYER_ORDER)]
  fwrite(layer_res[, .(system, layer_model, treatment, beta, se, ci_lo, ci_hi,
                        pct_change, p, p_adj, sig, interaction_p)],
         file.path(OUT_DIR, "off_target_layer_results.csv"))
  sig_layer <- layer_res[p_adj < 0.05]
  cat(sprintf("  Saved: off_target_layer_results.csv (%d rows)\n", nrow(layer_res)))
  if (nrow(sig_layer) > 0) {
    cat("  *** Significant off-target layer contrasts:\n")
    print(sig_layer[, .(system, layer_model, treatment, beta, pct_change, p_adj, sig)])
  } else {
    cat("  ✓ No significant layer-specific off-target effects.\n")
  }
} else {
  cat("  No systems had sufficient data for layer models.\n")
}
cat("\n")

# =============================================================================
# PART D — SPATIAL SPREAD: distance gradient
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART D: Spatial spread — β vs distance tier from VISp\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

# Per-area LMMs (pool layers within each area)
all_areas   <- unique(both_df$area)
area_res_list <- lapply(all_areas, function(ar) {
  d <- both_df[area == ar]
  n_per_tx <- d[, uniqueN(mouse_id), by = treatment]
  if (nrow(d) < 15 || sum(n_per_tx$V1 >= 2) < 2) return(NULL)
  fit <- tryCatch(
    fit_lmm(log_ratio ~ treatment + (1 | mouse_id), d, ar),
    error = function(e) {
      cat(sprintf("  ERROR [%s]: %s\n", ar, e$message)); NULL
    }
  )
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, model_label = ar)
  res[, area      := ar]
  res[, dist_tier := assign_dist_tier(ar)[1]]
  res[, system    := assign_system(ar)[1]]
  res
})
area_res <- rbindlist(Filter(Negate(is.null), area_res_list), fill = TRUE)
area_res  <- add_fdr(area_res, "p", by_col = "treatment", suffix = "")
area_res[, treatment := factor(treatment, levels = TREATMENT_ORDER[-1])]
cat(sprintf("  Per-area β computed: %d areas × %d treatments\n",
    uniqueN(area_res$area), uniqueN(area_res$treatment)))
fwrite(area_res[, .(area, dist_tier, system, treatment,
                     beta, se, ci_lo, ci_hi, pct_change, p, p_adj, sig)],
       file.path(OUT_DIR, "spread_per_area_betas.csv"))

# Weighted linear trend: β ~ dist_tier per treatment (weights = 1/SE²)
cat("\n  Weighted linear trend (β ~ dist_tier):\n")
cat(sprintf("  %-20s  slope    SE      p       interpretation\n", "treatment"))
cat(sprintf("  %s\n", strrep("─", 75)))

trend_res <- rbindlist(lapply(levels(area_res$treatment), function(tx) {
  d <- area_res[treatment == tx & is.finite(beta) & is.finite(se) & se > 0]
  if (nrow(d) < 4) return(NULL)
  w  <- 1 / d$se^2
  s  <- summary(lm(beta ~ dist_tier, data = d, weights = w))$coefficients
  slope   <- s["dist_tier", "Estimate"]
  slope_se <- s["dist_tier", "Std. Error"]
  slope_p  <- s["dist_tier", "Pr(>|t|)"]
  interp   <- if (slope_p < 0.05 && slope < 0) "decays with distance ✓"
               else if (slope_p < 0.05 && slope > 0) "INCREASES with distance (concern)"
               else "no significant gradient"
  cat(sprintf("  %-20s  %+.4f  %.4f  %.4f  %s\n",
      tx, slope, slope_se, slope_p, interp))
  data.table(treatment = tx, slope, slope_se, slope_p, interpretation = interp)
}))
cat("\n")

# =============================================================================
# RESIDUAL DIAGNOSTICS
# =============================================================================
diag_panel <- function(fit, title) {
  df_d <- data.frame(
    fitted    = fitted(fit),
    std_resid = as.numeric(scale(resid(fit, type = "pearson")))
  )
  pq <- ggplot(df_d, aes(sample = std_resid)) +
    stat_qq(size = 0.7, alpha = 0.4) + stat_qq_line(colour = "red", linewidth = 0.6) +
    labs(title = paste("QQ —", title), x = "Theoretical", y = "Std. resid.") +
    theme_bw(base_size = 9)
  pr <- ggplot(df_d, aes(x = fitted, y = std_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.5, alpha = 0.3) +
    geom_smooth(method = "loess", se = FALSE, colour = "red",
                linewidth = 0.5, formula = y ~ x) +
    labs(title = paste("Resid —", title), x = "Fitted", y = "Std. resid.") +
    theme_bw(base_size = 9)
  pq + pr
}
p_diag <- diag_panel(fit_spec, "B: specificity") +
  plot_annotation(title = "Residual diagnostics — script 08",
                  theme = theme(plot.title = element_text(face = "bold")))
for (ext in c("pdf","png"))
  ggsave(file.path(OUT_DIR, paste0("residual_check_rq8.", ext)), p_diag,
         width = 9, height = 4, dpi = ifelse(ext == "pdf", 300, 150))

# =============================================================================
# FIGURE 17 — System × Treatment heatmap
# =============================================================================
# Visual cortex row: inject script 01 primary betas for context
vis_ref <- data.table(
  system     = "Visual (target)",
  treatment  = factor(c("ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1","C6ST1_ADAMTS15"),
                      levels = TREATMENT_ORDER[-1]),
  beta       = c(-0.186, -0.256, -0.327, -0.027, -0.699),
  sig_local  = c(".",    "*",    "*",    "ns",   "***"),
  sig_global = c(".",    "*",    "*",    "ns",   "***")
)
plot_heat <- rbind(
  vis_ref,
  all_sys_res[treatment %in% TREATMENT_ORDER[-1],
    .(system, treatment = factor(treatment, levels = TREATMENT_ORDER[-1]),
      beta, sig_local, sig_global)],
  fill = TRUE
)
plot_heat[, system    := factor(system,    levels = rev(SYSTEM_ORDER))]
plot_heat[, treatment := factor(treatment, levels = rev(TREATMENT_ORDER[-1]))]

BETA_LIM <- 1.0
plot_heat[, fill_val := pmax(-BETA_LIM, pmin(BETA_LIM, beta))]
# Label: local sig; dagger if global disagrees
plot_heat[, label := fcase(
  sig_local %in% c("***","**","*") & sig_local != sig_global, paste0(sig_local, "†"),
  default = sig_local
)]
plot_heat[label == "ns", label := ""]

p17 <- ggplot(plot_heat, aes(x = treatment, y = system, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.2, vjust = 0.8, colour = "grey10") +
  geom_hline(yintercept = length(SYSTEM_ORDER) - 0.5,
             colour = "black", linewidth = 1.0) +
  scale_fill_gradient2(
    low = "#D55E00", mid = "white", high = "#0072B2",
    midpoint = 0, limits = c(-BETA_LIM, BETA_LIM), oob = squish,
    name = expression(beta ~ "(log)")
  ) +
  scale_x_discrete(name = NULL) + scale_y_discrete(name = NULL) +
  labs(
    title    = "RQ8 — Part A: PNN enwrapment across functional systems",
    subtitle = paste0(
      "β vs reference (log scale). Orange = reduction. Stars = local BH-FDR.\n",
      "† = significant only under local (not global) FDR.\n",
      "Visual cortex row from script 01 primary analysis (reference)."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"),
        panel.grid = element_blank(), legend.position = "right")

# =============================================================================
# FIGURE 18 — Specificity forest plot
# =============================================================================
res_spec_plot <- copy(res_spec)
res_spec_plot[, treatment := factor(treatment_label, levels = TREATMENT_ORDER[-1])]
res_spec_plot[, region    := factor(region,
    levels = c("Visual (target)", "Off-target (all)"))]

p18 <- ggplot(res_spec_plot,
              aes(x = beta, y = treatment,
                  colour = treatment, shape = region, alpha = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.25, linewidth = 0.9, position = position_dodge(0.55)) +
  geom_point(size = 4, position = position_dodge(0.55)) +
  geom_text(aes(x = ci_hi + 0.04, label = sig), hjust = 0, size = 3.5,
            position = position_dodge(0.55), show.legend = FALSE) +
  scale_colour_manual(values = PALETTE[-1], guide = "none") +
  scale_shape_manual(values = c("Visual (target)" = 18, "Off-target (all)" = 16), name = NULL) +
  scale_alpha_manual(values = c("Visual (target)" = 1.0, "Off-target (all)" = 0.55), name = NULL) +
  scale_x_continuous(name = expression(beta ~ "(log scale)"),
                     expand = expansion(mult = c(0.05, 0.25))) +
  scale_y_discrete(name = NULL) +
  labs(
    title    = "RQ8 — Part B: Specificity test (interaction model)",
    subtitle = sprintf(
      "treatment:is_visual F p = %.4f  %s\nDiamond = visual cortex  |  Circle = off-target (all areas pooled)",
      ix_p, ifelse(ix_p < 0.05, "✓ Specific", "NS"))
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.minor = element_blank())

# =============================================================================
# FIGURE 19 — Distance gradient
# =============================================================================
TIER_LABELS <- c("0" = "VISp\n(target)", "1" = "Higher\nvisual",
                  "2" = "RSP /\nPTLp",   "3" = "AUD /\nSS-post",
                  "4" = "SS / Motor\n/ Hipp", "5" = "Prefrontal\n/ Subcort")

# Weighted trend lines per treatment
trend_lines <- rbindlist(lapply(levels(area_res$treatment), function(tx) {
  d <- area_res[treatment == tx & is.finite(beta) & is.finite(se) & se > 0]
  if (nrow(d) < 4) return(NULL)
  w  <- 1 / d$se^2
  lm_fit <- lm(beta ~ dist_tier, data = d, weights = w)
  tiers  <- data.table(dist_tier = 0:5)
  pred   <- predict(lm_fit, newdata = tiers, se.fit = TRUE)
  tiers[, `:=`(treatment = tx, pred = pred$fit,
               pred_lo = pred$fit - 1.96 * pred$se.fit,
               pred_hi = pred$fit + 1.96 * pred$se.fit)]
  tiers
}), fill = TRUE)
if (nrow(trend_lines) > 0)
  trend_lines[, treatment := factor(treatment, levels = TREATMENT_ORDER[-1])]

plot_spread <- copy(area_res)
set.seed(42)
plot_spread[, tier_jitter := dist_tier + runif(.N, -0.18, 0.18)]

p19 <- ggplot(plot_spread, aes(x = dist_tier, y = beta, colour = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 1.5, linetype = "dotted", colour = "grey60") +
  {if (nrow(trend_lines) > 0)
     geom_ribbon(data = trend_lines,
                 aes(x = dist_tier, ymin = pred_lo, ymax = pred_hi,
                     fill = treatment, y = NULL),
                 alpha = 0.15, inherit.aes = FALSE)} +
  {if (nrow(trend_lines) > 0)
     geom_line(data = trend_lines,
               aes(x = dist_tier, y = pred, colour = treatment),
               linewidth = 0.9, inherit.aes = FALSE)} +
  geom_point(aes(x = tier_jitter), size = 2.5, alpha = 0.75) +
  scale_colour_manual(values = PALETTE[-1], name = NULL) +
  scale_fill_manual(  values = PALETTE[-1], name = NULL, guide = "none") +
  scale_x_continuous(breaks = 0:5, labels = TIER_LABELS,
                     name = "Anatomical distance from injection site (VISp)") +
  scale_y_continuous(name = expression(beta ~ "(log) vs mScarlet")) +
  facet_wrap(~ treatment, ncol = 3) +
  labs(
    title    = "RQ8 — Part D: Spatial spread from injection site",
    subtitle = paste0(
      "Each point = one brain area. Line = weighted linear trend (1/SE²).\n",
      "Negative slope = effect decays with distance (localised).\n",
      "Dotted line: right of line = off-target (tiers 2–5)."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7))

# =============================================================================
# FIGURE 20 — Off-target layer profiles (Part C, if any ran)
# =============================================================================
if (nrow(layer_res) > 0) {
  plot_lr <- layer_res[treatment %in% TREATMENT_ORDER[-1]]
  plot_lr[, treatment  := factor(treatment,  levels = TREATMENT_ORDER[-1])]
  plot_lr[, layer_model := factor(layer_model, levels = LAYER_ORDER)]
  p20 <- ggplot(plot_lr, aes(x = layer_model, y = beta,
                              colour = treatment, group = treatment)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(aes(shape = (p_adj < 0.05)), size = 3) +
    scale_colour_manual(values = PALETTE[-1], name = NULL) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 19),
                       name = NULL, labels = c("ns", "p_adj < 0.05")) +
    facet_wrap(~ system, ncol = 2) +
    labs(title    = "RQ8 — Part C: Layer profiles in off-target systems",
         subtitle = "β vs mScarlet. Filled = p_adj < 0.05 (BH per system × layer).",
         x = "Layer", y = expression(beta ~ "(log)")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"), legend.position = "right")
  for (ext in c("pdf","png"))
    ggsave(file.path(OUT_DIR, paste0("fig20_off_target_layer.", ext)), p20,
           width = 12, height = 6, dpi = ifelse(ext == "pdf", 300, 150))
}

# =============================================================================
# SAVE FIGURES 17–19
# =============================================================================
for (ext in c("pdf","png")) {
  dpi <- ifelse(ext == "pdf", 300, 150)
  ggsave(file.path(OUT_DIR, paste0("fig17_system_heatmap.",       ext)), p17,
         width = 11, height = 7,  dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig18_specificity_contrast.", ext)), p18,
         width = 9,  height = 5,  dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig19_distance_gradient.",    ext)), p19,
         width = 12, height = 8,  dpi = dpi)
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("RQ8 SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Part A — Off-target hits (local FDR):   %d\n",
    nrow(all_sys_res[sig_local  != "ns"])))
cat(sprintf("  Part A — Off-target hits (global FDR):  %d\n",
    nrow(all_sys_res[sig_global != "ns"])))
cat(sprintf("  Part A — FDR discrepancies:             %d\n", nrow(discr)))
cat(sprintf("  Part B — treatment:is_visual F p = %.4f  %s\n", ix_p,
    ifelse(ix_p < 0.05, "✓ Specific", "NS")))
cat(sprintf("  Part C — Off-target layer models ran:   %d systems\n",
    uniqueN(layer_res$system)))
if (!is.null(trend_res) && nrow(trend_res) > 0) {
  cat("\n  Part D — Distance gradient slopes:\n")
  print(trend_res[, .(treatment, slope, slope_p, interpretation)])
}

cat(sprintf("\n── Outputs → %s\n", OUT_DIR))
cat("   off_target_atlas_enwrapment.csv   system_lmm_results.csv\n")
cat("   fdr_comparison.csv                visual_vs_other_interaction.csv\n")
cat("   off_target_layer_results.csv      spread_per_area_betas.csv\n")
cat("   fig17–20 (pdf + png)              residual_check_rq8\n")