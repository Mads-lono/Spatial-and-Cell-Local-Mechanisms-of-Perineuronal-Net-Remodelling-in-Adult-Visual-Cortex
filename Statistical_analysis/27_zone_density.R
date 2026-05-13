# =============================================================================
# Script 27: Zone-level cell density LMM
# =============================================================================
# PURPOSE:
#   Extends script 11 by computing cell density at the zone level
#   (Core, Penumbra, Outside) rather than hemisphere level only.
#
#   Script 11 was restricted to zone = "Total" because the main merged_dataset
#   only contains hemisphere-level area. This script uses merged_dataset_zones,
#   which contains areaMm2 per brain_area × zone × hemisphere, and joins it to
#   per-zone cell counts derived from cells_with_zones.csv.
#
#   Density = n_cells / areaMm2, where areaMm2 is the total tissue area in
#   that zone across all sections contributing to that animal × region record.
#   This is consistent with how the merged dataset aggregates across slices.
#
# Model structure (mirrors script 11):
#   Part A — Primary zone model: log(density) ~ treatment + (1|animal_id)
#             Fitted separately per zone × staining, visual cortex only,
#             ipsilateral hemisphere. BH-FDR correction within each stratum.
#
#   Part B — Zone × treatment interaction:
#             log(density) ~ treatment * zone + (1|animal_id)
#             Tests whether the treatment effect on density differs across
#             zones (Core vs Penumbra vs Outside). LRT vs additive model.
#
#   Part C — PV density specificity control (zone-level):
#             Same models as Parts A and B applied to PV staining.
#             PV density should not change if enzymatic effects are specific
#             to PNN components.
#
#   Part D — Planned contrasts for RQ3 and RQ4 (pre-specified, no correction):
#             C6ST1_ADAMTS15 vs ADAMTS15, ADAMTS4_MD vs ADAMTS4
#             Fitted per zone, WFA and PV.
#
# Outputs: results/27_zone_density/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(patchwork)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
MERGED_ZONES_CSV <- "/path/to/merged_datasets/merged_dataset_zones.csv"
CELLS_CSV        <- "/path/to/analysis_results/cells_with_zones.csv"
RESULTS_DIR      <- "/path/to/results"
OUT_DIR          <- file.path(RESULTS_DIR, "27_zone_density")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENT_ORDER  <- c("mScarlet","ADAMTS4","ADAMTS4_MD",
                      "ADAMTS15","C6ST1","C6ST1_ADAMTS15")
ZONE_ORDER       <- c("Core","Penumbra","Outside")
STAIN_ORDER      <- c("WFA","PV")
EXCLUDE          <- ""    #Exclude animal as needed (e.g., injection failure)
VISUAL_RE        <- regex("isual", ignore_case = TRUE)
MIN_CELLS        <- 5L      # minimum cells per animal × region × zone
MIN_AREA         <- 0.001   # minimum area in mm² (filters degenerate zones)
ALPHA            <- 0.05

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

cat("Script 27: Zone-level cell density LMM\n")
cat("============================================================\n\n")

# ── Helpers ───────────────────────────────────────────────────────────────────

fit_lmm <- function(formula, data, label, reml = TRUE) {
  fit <- tryCatch(
    lmer(formula, data = data, REML = reml,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) {
      cat(sprintf("  ERROR [%s]: %s\n", label, e$message)); NULL
    }
  )
  if (!is.null(fit)) {
    if (isSingular(fit)) cat(sprintf("  SINGULAR [%s]\n", label))
    else                  cat(sprintf("  OK        [%s]\n", label))
  }
  fit
}

extract_contrasts <- function(fit, model_label = "") {
  if (is.null(fit)) return(NULL)
  coefs <- as.data.frame(coef(summary(fit, ddf = "Satterthwaite")))
  coefs$term <- rownames(coefs)
  setDT(coefs)
  res <- coefs[term != "(Intercept)"]
  old <- intersect(names(res), c("Estimate","Std. Error","df","t value","Pr(>|t|)"))
  setnames(res, old, c("beta","se","df","t","p")[seq_along(old)])
  res[, model := model_label]
  res
}

add_fdr <- function(dt, suffix = "_adj", by_col = NULL) {
  if (is.null(by_col)) {
    dt[, paste0("p", suffix) := p.adjust(p, method = "BH")]
  } else {
    dt[, paste0("p", suffix) := p.adjust(p, method = "BH"), by = by_col]
  }
  dt[, paste0("sig", suffix) := fcase(
    get(paste0("p", suffix)) < 0.001, "***",
    get(paste0("p", suffix)) < 0.01,  "**",
    get(paste0("p", suffix)) < 0.05,  "*",
    get(paste0("p", suffix)) < 0.10,  ".",
    default = "ns"
  )]
  dt
}

pct_change <- function(beta) (exp(beta) - 1) * 100

# =============================================================================
# STEP 1: Load area data from merged_dataset_zones
# =============================================================================
cat("Loading merged_dataset_zones for area data ...\n")

mz <- fread(MERGED_ZONES_CSV)
mz <- fix_names(mz)
mz <- mz[animal_id != EXCLUDE]
mz <- mz[zone %in% ZONE_ORDER]

# Keep one row per animal × brain_area × zone × hemisphere × staining
# areaMm2 is already the total area across all sections (n_slices tells how many)
# Use the WFA rows for area (PV and WFA share the same tissue area per zone)
area_dt <- unique(mz[staining == "WFA",
                      .(animal_id, treatment, brain_area, zone, hemisphere,
                        areaMm2, n_slices)])
area_dt <- area_dt[areaMm2 > MIN_AREA]

cat(sprintf("  Area records (WFA, zone-level): %d\n", nrow(area_dt)))
cat(sprintf("  Animals: %d\n\n", uniqueN(area_dt$animal_id)))

# =============================================================================
# STEP 2: Count cells per animal × brain_area × zone × hemisphere × staining
# =============================================================================
cat("Loading cells_with_zones.csv for cell counts ...\n")

cells <- fread(CELLS_CSV,
               select = c("mouse_id","slice_id","cell_type","hemisphere",
                          "x_hires","y_hires","brain_area","zone","treatment"))

# Rename mouse_id → animal_id for join consistency
setnames(cells, "mouse_id", "animal_id")
setnames(cells, "cell_type", "staining_raw")

# Apply corrections
cells <- fix_names(cells)
cells <- cells[animal_id != EXCLUDE]
cells <- cells[str_detect(brain_area, VISUAL_RE)]
cells <- cells[zone %in% ZONE_ORDER]

# Map cell_type to staining label: PNN → WFA, PV → PV
cells[, staining := fcase(
  staining_raw == "PNN", "WFA",
  staining_raw == "PV",  "PV",
  default = NA_character_
)]
cells <- cells[!is.na(staining)]

# Hemisphere: left = ipsilateral (injected side)
cells[, hemisphere := fcase(
  hemisphere == "left",  "left",
  hemisphere == "right", "right",
  default = NA_character_
)]
cells <- cells[!is.na(hemisphere)]

# Count cells per animal × brain_area × zone × hemisphere × staining
cell_counts <- cells[, .(n_cells = .N),
                      by = .(animal_id, treatment, brain_area, zone,
                             hemisphere, staining)]

cat(sprintf("  Cell count records: %d\n\n", nrow(cell_counts)))

# =============================================================================
# STEP 3: Join counts to area and compute density
# =============================================================================
cat("Computing zone-level density ...\n")

# Join on animal_id × brain_area × zone × hemisphere
# Use staining-specific join so WFA cell counts get WFA area and PV get PV area
# (area is the same for both but the join structure mirrors the data model)
density_dt <- merge(
  cell_counts,
  area_dt[, .(animal_id, brain_area, zone, hemisphere, areaMm2)],
  by   = c("animal_id","brain_area","zone","hemisphere"),
  all.x = TRUE
)

# Exclude records with no area match or insufficient cells/area
density_dt <- density_dt[!is.na(areaMm2) & areaMm2 > MIN_AREA]
density_dt <- density_dt[n_cells >= MIN_CELLS]

# Compute density
density_dt[, density := n_cells / areaMm2]
density_dt[, log_density := log(density)]

# Apply treatment factor
density_dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
density_dt[, zone      := factor(zone,      levels = ZONE_ORDER)]
density_dt[, staining  := factor(staining,  levels = STAIN_ORDER)]

# Restrict to visual cortex, ipsilateral hemisphere for primary models
density_vis_ipsi <- density_dt[
  str_detect(brain_area, VISUAL_RE) & hemisphere == "left"]

cat(sprintf("Zone-level density records (visual cortex, ipsilateral): %d\n",
            nrow(density_vis_ipsi)))
print(density_vis_ipsi[, .N, by = .(staining, zone)])
cat("\n")

fwrite(density_dt, file.path(OUT_DIR, "zone_density_all.csv"))

# =============================================================================
# PART A — Primary zone-stratified models: per zone per staining
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part A: Zone-stratified density models (visual cortex, ipsilateral)\n")
cat("        log(density) ~ treatment + (1|animal_id), per zone × staining\n")
cat("══════════════════════════════════════════════════════\n\n")

zone_res_a <- rbindlist(lapply(STAIN_ORDER, function(st) {
  rbindlist(lapply(ZONE_ORDER, function(z) {
    d <- density_vis_ipsi[staining == st & zone == z]
    if (uniqueN(d$animal_id) < 3) return(NULL)
    fit <- fit_lmm(log_density ~ treatment + (1 | animal_id), d,
                   sprintf("A_%s_%s", st, z))
    if (is.null(fit)) return(NULL)
    res <- extract_contrasts(fit, sprintf("A_%s_%s", st, z))
    res[, `:=`(staining = st, zone = z, pct_change = pct_change(beta))]
    res
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(zone_res_a) && nrow(zone_res_a) > 0) {
  zone_res_a <- add_fdr(zone_res_a, suffix = "_adj",
                         by_col = c("staining","zone"))
  fwrite(zone_res_a, file.path(OUT_DIR, "partA_zone_density_contrasts.csv"))
  cat("\nSignificant zone-level density results (FDR-corrected):\n")
  sig_a <- zone_res_a[sig_adj != "ns"]
  if (nrow(sig_a) > 0) {
    print(sig_a[, .(staining, zone, term,
                    beta = round(beta,3), pct_change = round(pct_change,1),
                    p = round(p,4), p_adj = round(p_adj,4), sig_adj)])
  } else {
    cat("  No significant hits after FDR correction\n")
  }
}

# =============================================================================
# PART B — Zone × treatment interaction LMM
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part B: Zone × treatment interaction\n")
cat("        log(density) ~ treatment * zone + (1|animal_id)\n")
cat("══════════════════════════════════════════════════════\n\n")

zone_lrt_b <- rbindlist(lapply(STAIN_ORDER, function(st) {
  d <- density_vis_ipsi[staining == st]
  if (uniqueN(d$animal_id) < 3) return(NULL)
  cat(sprintf("  %s staining:\n", st))

  fit_inter <- fit_lmm(log_density ~ treatment * zone + (1 | animal_id),
                        d, sprintf("B_inter_%s", st), reml = FALSE)
  fit_add   <- fit_lmm(log_density ~ treatment + zone + (1 | animal_id),
                        d, sprintf("B_add_%s",   st), reml = FALSE)
  if (is.null(fit_inter) || is.null(fit_add)) return(NULL)

  lrt <- anova(fit_add, fit_inter)
  cat(sprintf("    LRT p = %.4f\n", lrt[2, "Pr(>Chisq)"]))

  lrt_dt <- as.data.table(lrt, keep.rownames = "model")
  lrt_dt[, staining := st]
  lrt_dt
}), fill = TRUE)

if (!is.null(zone_lrt_b) && nrow(zone_lrt_b) > 0) {
  fwrite(zone_lrt_b, file.path(OUT_DIR, "partB_zone_interaction_lrt.csv"))
  cat("\nLRT summary:\n")
  print(zone_lrt_b[, .(staining, model, AIC, Chisq = round(Chisq,3),
                         Df, `Pr(>Chisq)` = round(`Pr(>Chisq)`,4))])
}

# Zone-specific emmeans from interaction model (WFA)
wfa_data <- density_vis_ipsi[staining == "WFA"]
if (uniqueN(wfa_data$animal_id) >= 3) {
  fit_wfa_inter <- fit_lmm(
    log_density ~ treatment * zone + (1 | animal_id),
    wfa_data, "B_emm_WFA", reml = TRUE)

  if (!is.null(fit_wfa_inter)) {
    em_b <- emmeans(fit_wfa_inter, ~ treatment | zone)
    ct_b <- as.data.table(
      contrast(em_b, method = "trt.vs.ctrl", ref = "mScarlet", adjust = "BH"))
    if ("estimate" %in% names(ct_b)) setnames(ct_b, "estimate", "beta")
    if ("p.value"  %in% names(ct_b)) setnames(ct_b, "p.value",  "p")
    ct_b[, `:=`(pct_change = pct_change(beta),
                sig = fcase(p < 0.001,"***", p < 0.01,"**",
                             p < 0.05,"*",   p < 0.10,".",
                             default = "ns"))]
    cat("\nWFA density contrasts per zone (BH-corrected within zone):\n")
    print(ct_b[sig != "ns",
               .(zone, contrast, beta = round(beta,3),
                 pct_change = round(pct_change,1),
                 p = round(p,4), sig)])
    fwrite(ct_b, file.path(OUT_DIR, "partB_wfa_zone_contrasts.csv"))
  }
}

# =============================================================================
# PART C — PV density specificity control (zone-level)
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part C: PV density specificity control (zone-level)\n")
cat("══════════════════════════════════════════════════════\n\n")

pv_data <- density_vis_ipsi[staining == "PV"]

pv_zone_res <- rbindlist(lapply(ZONE_ORDER, function(z) {
  d <- pv_data[zone == z]
  if (uniqueN(d$animal_id) < 3) return(NULL)
  fit <- fit_lmm(log_density ~ treatment + (1 | animal_id),
                  d, sprintf("C_PV_%s", z))
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, sprintf("C_PV_%s", z))
  res[, `:=`(staining = "PV", zone = z, pct_change = pct_change(beta))]
  res
}), fill = TRUE)

if (!is.null(pv_zone_res) && nrow(pv_zone_res) > 0) {
  pv_zone_res <- add_fdr(pv_zone_res, suffix = "_adj", by_col = "zone")
  fwrite(pv_zone_res, file.path(OUT_DIR, "partC_pv_zone_contrasts.csv"))
  cat("PV density zone contrasts (FDR-corrected):\n")
  sig_pv <- pv_zone_res[sig_adj != "ns"]
  if (nrow(sig_pv) > 0) {
    print(sig_pv[, .(zone, term, beta = round(beta,3),
                      pct_change = round(pct_change,1),
                      p = round(p,4), p_adj = round(p_adj,4), sig_adj)])
  } else {
    cat("  No significant PV density changes — consistent with PNN specificity\n")
  }
}

# =============================================================================
# PART D — Planned contrasts (RQ3 and RQ4), per zone, no correction
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part D: Planned contrasts — RQ3 and RQ4, per zone\n")
cat("══════════════════════════════════════════════════════\n\n")

planned_d <- rbindlist(lapply(STAIN_ORDER, function(st) {
  rbindlist(lapply(ZONE_ORDER, function(z) {
    d <- density_vis_ipsi[staining == st & zone == z]
    if (uniqueN(d$animal_id) < 3) return(NULL)
    fit <- fit_lmm(log_density ~ treatment + (1 | animal_id),
                    d, sprintf("D_%s_%s", st, z))
    if (is.null(fit)) return(NULL)
    em  <- emmeans(fit, ~ treatment)
    # RQ3: C6ST1_ADAMTS15 vs ADAMTS15
    ct_rq3 <- tryCatch(
      as.data.table(contrast(em, method = list(
        "C6ST1_ADAMTS15 - ADAMTS15" = c(0,0,0,-1,0,1)),
        adjust = "none")),
      error = function(e) NULL)
    # RQ4: ADAMTS4_MD vs ADAMTS4
    ct_rq4 <- tryCatch(
      as.data.table(contrast(em, method = list(
        "ADAMTS4_MD - ADAMTS4" = c(0,-1,1,0,0,0)),
        adjust = "none")),
      error = function(e) NULL)
    res <- rbind(
      if (!is.null(ct_rq3)) ct_rq3[, rq := "RQ3"],
      if (!is.null(ct_rq4)) ct_rq4[, rq := "RQ4"],
      fill = TRUE
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    if ("estimate" %in% names(res)) setnames(res, "estimate", "beta")
    if ("p.value"  %in% names(res)) setnames(res, "p.value",  "p")
    res[, `:=`(staining = st, zone = z, pct_change = pct_change(beta),
               sig = fcase(p < 0.001,"***", p < 0.01,"**",
                            p < 0.05,"*",   p < 0.10,".",
                            default = "ns"))]
    res
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(planned_d) && nrow(planned_d) > 0) {
  fwrite(planned_d, file.path(OUT_DIR, "partD_planned_contrasts.csv"))
  cat("Planned contrasts (no correction — pre-specified):\n")
  print(planned_d[, .(rq, staining, zone, contrast,
                       beta = round(beta,3), pct_change = round(pct_change,1),
                       p = round(p,4), sig)])
}


# =============================================================================
# PART E — Factorial density model (RQ3): mirrors M8 from script 23
# =============================================================================
cat("\n══════════════════════════════════════════════════════\n")
cat("Part E: Factorial density model (RQ3)\n")
cat("        log(density) ~ ADAMTS15_present * C6ST1_present + (1|animal_id)\n")
cat("        Core zone only | ADAMTS4 as reference (neither component)\n")
cat("        P(interaction < 0) = density analog of M8\n")
cat("══════════════════════════════════════════════════════\n\n")

# Four-group factorial subset — Core zone, WFA, ipsilateral, visual cortex
fact_groups <- c("ADAMTS4","ADAMTS15","C6ST1","C6ST1_ADAMTS15")
fact_e <- density_vis_ipsi[
  staining == "WFA" & zone == "Core" &
  treatment %in% fact_groups]

fact_e[, ADAMTS15_present := as.integer(
  treatment %in% c("ADAMTS15","C6ST1_ADAMTS15"))]
fact_e[, C6ST1_present    := as.integer(
  treatment %in% c("C6ST1","C6ST1_ADAMTS15"))]

cat("Factorial groups (Core zone, WFA, n per animal):\n")
print(fact_e[, .N, by = .(treatment, ADAMTS15_present, C6ST1_present)])
cat("\n")

# Frequentist: interaction model vs additive model (LRT)
fit_e_inter <- fit_lmm(
  log_density ~ ADAMTS15_present * C6ST1_present + (1 | animal_id),
  fact_e, "E_inter", reml = FALSE)

fit_e_add <- fit_lmm(
  log_density ~ ADAMTS15_present + C6ST1_present + (1 | animal_id),
  fact_e, "E_add", reml = FALSE)

if (!is.null(fit_e_inter) && !is.null(fit_e_add)) {
  lrt_e <- anova(fit_e_add, fit_e_inter)
  cat("LRT: interaction vs additive (density, Core zone):\n")
  print(lrt_e)
  fwrite(as.data.table(lrt_e, keep.rownames = "model"),
         file.path(OUT_DIR, "partE_factorial_lrt.csv"))

  # Extract interaction coefficient from REML fit
  fit_e_reml <- fit_lmm(
    log_density ~ ADAMTS15_present * C6ST1_present + (1 | animal_id),
    fact_e, "E_inter_reml", reml = TRUE)

  if (!is.null(fit_e_reml)) {
    coef_e <- as.data.frame(coef(summary(fit_e_reml, ddf = "Satterthwaite")))
    coef_e$term <- rownames(coef_e)
    setDT(coef_e)
    setnames(coef_e,
             intersect(c("Estimate","Std. Error","df","t value","Pr(>|t|)"),
                       names(coef_e)),
             c("beta","se","df","t","p")[
               seq_along(intersect(
                 c("Estimate","Std. Error","df","t value","Pr(>|t|)"),
                 names(coef_e)))])
    coef_e[, pct_change := pct_change(beta)]
    coef_e[, sig := fcase(p < 0.001,"***", p < 0.01,"**",
                           p < 0.05,"*",   p < 0.10,".",
                           default = "ns")]

    cat("\nFactorial model coefficients:\n")
    print(coef_e[, .(term, beta = round(beta,3), se = round(se,3),
                      df = round(df,1), t = round(t,3),
                      p = round(p,4), pct_change = round(pct_change,1), sig)])
    fwrite(coef_e, file.path(OUT_DIR, "partE_factorial_coefficients.csv"))

    inter_row <- coef_e[str_detect(term, ":")]
    if (nrow(inter_row) > 0) {
      cat(sprintf("\nInteraction (ADAMTS15 × C6ST1): β = %.3f, p = %.4f, %s\n",
                  inter_row$beta, inter_row$p, inter_row$sig))
      cat(sprintf("Interpretation: combination produces %.1f%% %s than\n",
                  abs(inter_row$pct_change),
                  ifelse(inter_row$beta < 0, "more depletion","less depletion")))
      cat(sprintf("  expected from additive individual contributions\n"))
    }
  }
}

# Bayesian version: P(interaction < 0)
cat("\nFitting Bayesian factorial model (brms) ...\n")
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(posterior))

priors_e <- c(
  prior(normal(0, 1), class = b),
  prior(normal(0, 1), class = sigma)
)

# Aggregate to animal level (one value per animal) to mirror M8
animal_e <- fact_e[, .(
  log_density = mean(log_density),
  n_records   = .N
), by = .(animal_id, treatment, ADAMTS15_present, C6ST1_present)]

m_e_bayes <- brm(
  log_density ~ ADAMTS15_present * C6ST1_present,
  data    = animal_e,
  prior   = priors_e,
  chains  = 4L, iter = 4000L, warmup = 2000L,
  cores   = 4L, seed  = 2025L,
  control = list(adapt_delta = 0.95),
  file    = file.path(OUT_DIR, "m_bayes_factorial_density"),
  silent  = 2L
)

draws_e     <- as_draws_df(m_e_bayes)
inter_col_e <- grep("ADAMTS15_present:C6ST1_present",
                     names(draws_e), value = TRUE)[1]

if (!is.na(inter_col_e)) {
  inter_e_draws <- draws_e[[inter_col_e]]
  p_inter_neg_e <- mean(inter_e_draws < 0)

  cat(sprintf("\nBayesian factorial density model:\n"))
  cat(sprintf("  Interaction mean = %.4f, 89%% CI [%.4f, %.4f]\n",
              mean(inter_e_draws),
              quantile(inter_e_draws, 0.055),
              quantile(inter_e_draws, 0.945)))
  cat(sprintf("  P(interaction < 0) = %.4f\n", p_inter_neg_e))
  cat(sprintf("  Interpretation: %s\n",
              ifelse(p_inter_neg_e > 0.95,
                     "Strong evidence: combination reduces density more than sum of parts",
                     ifelse(p_inter_neg_e > 0.80,
                            "Moderate evidence for super-additive density reduction",
                            "Weak evidence for interaction on density"))))

  bayes_e_summary <- data.table(
    outcome          = "WFA density (Core zone)",
    interaction_mean = round(mean(inter_e_draws), 4),
    ci_89_lo         = round(quantile(inter_e_draws, 0.055), 4),
    ci_89_hi         = round(quantile(inter_e_draws, 0.945), 4),
    p_inter_neg      = round(p_inter_neg_e, 4),
    interpretation   = ifelse(p_inter_neg_e > 0.80,
                              "super-additive density reduction",
                              "no super-additive effect on density")
  )
  fwrite(bayes_e_summary, file.path(OUT_DIR, "partE_bayesian_summary.csv"))
  cat("  Saved partE_bayesian_summary.csv\n\n")
} else {
  cat("  Could not identify interaction column in Bayesian draws.\n\n")
}

# =============================================================================
# PART F — Core/Outside density ratio (RQ3): within-animal spatial contrast
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part F: Core/Outside density ratio (RQ3)\n")
cat("        log(Core density / Outside density) per animal\n")
cat("        Pre-specified contrast: C6ST1_ADAMTS15 vs ADAMTS15 (no correction)\n")
cat("══════════════════════════════════════════════════════\n\n")

# Compute per-animal mean density in Core and Outside (WFA, visual cortex, ipsi)
core_dens <- density_vis_ipsi[staining == "WFA" & zone == "Core", .(
  density_core = mean(density)
), by = .(animal_id, treatment)]

out_dens <- density_vis_ipsi[staining == "WFA" & zone == "Outside", .(
  density_outside = mean(density)
), by = .(animal_id, treatment)]

ratio_f <- merge(core_dens, out_dens, by = c("animal_id","treatment"))
ratio_f <- ratio_f[density_core > 0 & density_outside > 0]
ratio_f[, log_ratio := log(density_core / density_outside)]
ratio_f[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat("Core/Outside log density ratio per animal:\n")
print(ratio_f[order(treatment),
               .(animal_id, treatment,
                 density_core    = round(density_core, 1),
                 density_outside = round(density_outside, 1),
                 log_ratio       = round(log_ratio, 3))])

fwrite(ratio_f, file.path(OUT_DIR, "partF_core_outside_ratio.csv"))

# Primary LMM: log_ratio ~ treatment
fit_f <- lm(log_ratio ~ treatment, data = ratio_f)

em_f   <- emmeans(fit_f, ~ treatment)

# Pre-specified contrast: C6ST1_ADAMTS15 vs ADAMTS15 (no correction)
ct_f_rq3 <- tryCatch(
  as.data.table(contrast(em_f,
    method = list("C6ST1_ADAMTS15 - ADAMTS15" = c(0,0,0,-1,0,1)),
    adjust = "none")),
  error = function(e) NULL)

# Pre-specified contrast: ADAMTS4_MD vs ADAMTS4 (no correction)
ct_f_rq4 <- tryCatch(
  as.data.table(contrast(em_f,
    method = list("ADAMTS4_MD - ADAMTS4" = c(0,-1,1,0,0,0)),
    adjust = "none")),
  error = function(e) NULL)

planned_f <- rbind(
  if (!is.null(ct_f_rq3)) ct_f_rq3[, rq := "RQ3"],
  if (!is.null(ct_f_rq4)) ct_f_rq4[, rq := "RQ4"],
  fill = TRUE
)

if (!is.null(planned_f) && nrow(planned_f) > 0) {
  if ("estimate" %in% names(planned_f)) setnames(planned_f, "estimate", "beta")
  if ("p.value"  %in% names(planned_f)) setnames(planned_f, "p.value",  "p")
  planned_f[, pct_change := pct_change(beta)]
  planned_f[, sig := fcase(p < 0.001,"***", p < 0.01,"**",
                             p < 0.05,"*",   p < 0.10,".",
                             default = "ns")]
  cat("\nPre-specified contrasts on Core/Outside log density ratio:\n")
  print(planned_f[, .(rq, contrast, beta = round(beta,3),
                       SE = round(SE,3), df = round(df,1),
                       t.ratio = round(t.ratio,3),
                       p = round(p,4), sig)])
  fwrite(planned_f, file.path(OUT_DIR, "partF_planned_contrasts.csv"))
}

# Group means for reporting
means_f <- ratio_f[, .(
  mean_ratio = mean(log_ratio),
  sd_ratio   = sd(log_ratio),
  n          = .N
), by = treatment]
cat("\nGroup means (Core/Outside log density ratio):\n")
print(means_f[order(factor(treatment, levels = TREATMENT_ORDER)),
               .(treatment,
                 mean_ratio = round(mean_ratio, 3),
                 sd_ratio   = round(sd_ratio, 3), n)])
fwrite(means_f, file.path(OUT_DIR, "partF_group_means.csv"))

# Figure F: dot plot of Core/Outside ratio per animal
p_f <- ggplot(ratio_f,
              aes(x = treatment, y = log_ratio, colour = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.5) +
  geom_jitter(width = 0.1, size = 3.5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.35, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  scale_x_discrete(limits = TREATMENT_ORDER) +
  labs(
    title    = "Part F: Core/Outside WFA density ratio per animal",
    subtitle = "log ratio < 0 = fewer PNN cells in Core than Outside (injection-site depletion)",
    x        = NULL,
    y        = "log(Core density / Outside density)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig_partF_core_outside_ratio.pdf"),
       p_f, width = 8, height = 5)
cat("  Saved fig_partF_core_outside_ratio.pdf\n\n")


# =============================================================================
# FIGURES
# =============================================================================
cat("\nGenerating figures ...\n")

# Figure 1: WFA density % change per zone — forest plot
if (!is.null(zone_res_a) && nrow(zone_res_a) > 0) {
  wfa_a <- zone_res_a[staining == "WFA" & str_detect(term, "treatment")]
  wfa_a[, tx := str_remove(term, "treatment")]
  wfa_a[, tx := factor(tx, levels = rev(TREATMENT_ORDER[-1]))]
  wfa_a[, zone := factor(zone, levels = ZONE_ORDER)]

  p1 <- ggplot(wfa_a, aes(x = pct_change, y = tx, colour = tx)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = pct_change - 1.96*se*100,
                      xmax = pct_change + 1.96*se*100),
                  width = 0.25, linewidth = 0.7) +
    geom_point(size = 3) +
    geom_text(aes(x = pct_change + sign(pct_change)*5,
                  label = sig_adj),
              hjust = 0.5, size = 3.5, fontface = "bold",
              colour = "black") +
    facet_wrap(~ zone, ncol = 3) +
    scale_colour_manual(values = PALETTE[-1], guide = "none") +
    labs(
      title    = "WFA cell density: % change from mScarlet per zone (ipsilateral V1)",
      subtitle = "Error bars = 95% CI | * FDR < 0.05",
      x        = "% change in cells/mm²",
      y        = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(OUT_DIR, "fig1_wfa_density_per_zone.pdf"),
         p1, width = 12, height = 4)
  cat("  Saved fig1_wfa_density_per_zone.pdf\n")
}

# Figure 2: WFA vs PV density in Core zone side by side
core_both <- density_vis_ipsi[zone == "Core", .(
  mean_density = mean(density),
  sd_density   = sd(density),
  n            = .N
), by = .(treatment, staining)]
core_both[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
core_both[, staining  := factor(staining,  levels = STAIN_ORDER)]

p2 <- ggplot(core_both, aes(x = treatment, y = mean_density,
                               fill = treatment)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_density - sd_density,
                    ymax = mean_density + sd_density),
                width = 0.25) +
  facet_wrap(~ staining, scales = "free_y") +
  scale_fill_manual(values = PALETTE, guide = "none") +
  labs(
    title = "Mean cell density in Core zone by treatment — WFA vs PV",
    x     = NULL,
    y     = "Cells per mm²"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig2_core_density_wfa_vs_pv.pdf"),
       p2, width = 9, height = 4)
cat("  Saved fig2_core_density_wfa_vs_pv.pdf\n")

# Figure 3: Density across all three zones — line plot per treatment (WFA)
zone_means <- density_vis_ipsi[staining == "WFA", .(
  mean_density = mean(density)
), by = .(treatment, zone)]
zone_means[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
zone_means[, zone      := factor(zone,      levels = ZONE_ORDER)]

p3 <- ggplot(zone_means,
             aes(x = zone, y = mean_density,
                 colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_colour_manual(values = PALETTE) +
  labs(
    title   = "WFA cell density across zones (ipsilateral V1)",
    x       = "Zone",
    y       = "Mean cells per mm²",
    colour  = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig3_wfa_density_zone_profile.pdf"),
       p3, width = 8, height = 5)
cat("  Saved fig3_wfa_density_zone_profile.pdf\n\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Script 27 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  zone_density_all.csv\n")
cat("  partA_zone_density_contrasts.csv\n")
cat("  partB_zone_interaction_lrt.csv\n")
cat("  partB_wfa_zone_contrasts.csv\n")
cat("  partC_pv_zone_contrasts.csv\n")
cat("  partD_planned_contrasts.csv\n")
cat("  fig1_wfa_density_per_zone.pdf\n")
cat("  fig2_core_density_wfa_vs_pv.pdf\n")
cat("  fig3_wfa_density_zone_profile.pdf\n")
cat("  partE_factorial_lrt.csv\n")
cat("  partE_factorial_coefficients.csv\n")
cat("  partE_bayesian_summary.csv\n")
cat("  partF_core_outside_ratio.csv\n")
cat("  partF_planned_contrasts.csv\n")
cat("  partF_group_means.csv\n")
cat("  fig_partF_core_outside_ratio.pdf\n")
cat("══════════════════════════════════════════════════════\n")
