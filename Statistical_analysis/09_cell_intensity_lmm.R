# =============================================================================
# 09_cell_intensity_lmm.R
# Cell-level WFA and PV fluorescence intensity LMM
#
# Outcome: normalized_btm20 (cell fluorescence intensity normalised against
#   the 20th percentile of the Outside zone background, per cell)
#
# Research questions addressed:
#   RQ1 — Does any construct reduce WFA/PV cell intensity in V1?
#   RQ2 — Is any effect spatially graded across injection zones?
#   RQ3 — Does C6ST-1 augment the effect of ADAMTS-15? (planned contrast)
#   RQ4 — Does monocular deprivation augment ADAMTS-4? (planned contrast)
#   RQ5 — Do effects differ across cortical layers?
#   RQ6 — Is any effect lateralised to the injected hemisphere?
#
# Model structure:
#   Part A — Primary model: intensity ~ treatment + (1|mouse_id)
#             Pooled across zones, ipsilateral hemisphere, visual cortex only
#   Part B — Zone-stratified models: run separately for Core, Penumbra, Outside
#             Tests spatial gradient (RQ2)
#   Part C — Zone × treatment interaction: intensity ~ treatment * zone + (1|mouse_id)
#             Formally tests whether treatment effect differs across zones (RQ2)
#   Part D — Layer-stratified models: run separately per cortical layer (RQ5)
#   Part E — Hemisphere comparison: intensity ~ treatment * hemisphere + (1|mouse_id)
#             Tests lateralisation of effect (RQ6)
#   Part F — Planned contrasts: C6ST1_ADAMTS15 vs ADAMTS15 (RQ3)
#                               ADAMTS4_MD vs ADAMTS4 (RQ4)
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(lme4)
  library(lmerTest)
  if (!requireNamespace("emmeans", quietly = TRUE))
    install.packages("emmeans", repos = "https://cloud.r-project.org")
  library(emmeans)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV   <- "/path/to/analysis_results/cells_with_zones.csv"
RESULTS_DIR <- "/path/to/results"
OUT_DIR     <- file.path(RESULTS_DIR, "09_cell_intensity_lmm")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENTS      <- c("mScarlet","ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1","C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS
EXCLUDE_PRIMARY <- ""    # mouse_id to exclude from primary model (e.g. injection failure)
ZONE_ORDER      <- c("Core","Penumbra","Outside")
LAYER_ORDER     <- c("layer 1","layer 2/3","layer 4","layer 5","layer 6")
MIN_CELLS       <- 10L   # minimum cells per animal per zone/layer to include
ALPHA           <- 0.05

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

cat("── Script 09: Cell-level intensity LMM ─────────────────────────────────────\n")
cat("   Parts: (A) primary  (B) zone-stratified  (C) zone interaction\n")
cat("          (D) layer-stratified  (E) hemisphere  (F) planned contrasts\n\n")

# ── Helpers ───────────────────────────────────────────────────────────────────
fit_lmm <- function(formula, data, label, reml = TRUE) {
  fit <- tryCatch(
    lmer(formula, data = data, REML = reml,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { cat(sprintf("  ERROR [%s]: %s\n", label, e$message)); NULL }
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
  # handle both t-based and z-based column names from emmeans
  old_names <- intersect(names(res), c("Estimate","Std. Error","df","t value","Pr(>|t|)"))
  setnames(res, old_names, c("beta","se","df","t","p")[seq_along(old_names)])
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

# ── Load and correct data ──────────────────────────────────────────────────────
cat("Loading cells_with_zones.csv...\n")
dt <- fread(CELLS_CSV)



# Exclude injection failure
dt <- dt[mouse_id != EXCLUDE_PRIMARY]

# Restrict to ipsilateral hemisphere and visual cortex
# Visual cortex areas contain "visual" (case-insensitive) in brain_area
dt_ipsi <- dt[hemisphere == "left" &
               str_detect(tolower(brain_area), "visual")]

# Extract layer from brain_area (e.g. "Primary visual area, layer 2/3" -> "layer 2/3")
dt_ipsi[, layer := str_extract(brain_area, "layer [0-9/]+")]
dt_ipsi <- dt_ipsi[!is.na(layer)]
dt_ipsi[, layer := factor(layer, levels = LAYER_ORDER)]

# Set factor levels
dt_ipsi[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
dt_ipsi[, zone      := factor(zone,      levels = ZONE_ORDER)]

# Aggregate to section level: mean normalized_btm20 per mouse × slice × zone × cell_type
section_dt <- dt_ipsi[!is.na(normalized_btm20),
  .(intensity    = mean(normalized_btm20, na.rm = TRUE),
    log_intensity = mean(log(normalized_btm20[normalized_btm20 > 0]), na.rm = TRUE),
    n_cells       = .N),
  by = .(mouse_id, treatment, slice_id, zone, cell_type, layer)]

# Apply minimum cell threshold
section_dt <- section_dt[n_cells >= MIN_CELLS]

cat(sprintf("Sections after filtering: %d\n", nrow(section_dt)))
cat(sprintf("Animals: %d | Treatments: %s\n\n",
    uniqueN(section_dt$mouse_id),
    paste(sort(unique(section_dt$treatment)), collapse = ", ")))

# =============================================================================
# PART A — Primary model: pooled across zones, ipsilateral, visual cortex
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART A: Primary intensity model (RQ1)\n")
cat("        log_intensity ~ treatment + (1|mouse_id), per cell_type\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

primary_res <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  d <- section_dt[cell_type == ct]
  cat(sprintf("── Cell type: %s (%d sections, %d animals) ──\n",
      ct, nrow(d), uniqueN(d$mouse_id)))
  fit <- fit_lmm(log_intensity ~ treatment + (1 | mouse_id), d,
                 sprintf("Primary_%s", ct))
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, model_label = sprintf("primary_%s", ct))
  res[, `:=`(cell_type = ct,
             pct_change = pct_change(beta))]
  res
}), fill = TRUE)

if (!is.null(primary_res) && nrow(primary_res) > 0) {
  primary_res <- add_fdr(primary_res, suffix = "_adj", by_col = "cell_type")
  fwrite(primary_res, file.path(OUT_DIR, "primary_lmm_results.csv"))
  cat("\n── Part A results ───────────────────────────────────────────────────────────\n")
  print(primary_res[, .(cell_type, term, beta, pct_change, p, p_adj, sig_adj)])
}

# =============================================================================
# PART B — Zone-stratified models (RQ2)
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART B: Zone-stratified intensity models (RQ2)\n")
cat("        log_intensity ~ treatment + (1|mouse_id), per zone × cell_type\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

zone_res <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  rbindlist(lapply(ZONE_ORDER, function(z) {
    d <- section_dt[cell_type == ct & zone == z]
    if (uniqueN(d$mouse_id) < 3 || nrow(d) < MIN_CELLS) {
      cat(sprintf("  SKIP [%s, %s]: too few observations\n", ct, z))
      return(NULL)
    }
    cat(sprintf("── %s | %s (%d sections) ──\n", ct, z, nrow(d)))
    fit <- fit_lmm(log_intensity ~ treatment + (1 | mouse_id), d,
                   sprintf("%s_%s", ct, z))
    if (is.null(fit)) return(NULL)
    res <- extract_contrasts(fit, model_label = sprintf("zone_%s_%s", ct, z))
    res[, `:=`(cell_type = ct, zone = z, pct_change = pct_change(beta))]
    res
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(zone_res) && nrow(zone_res) > 0) {
  zone_res <- add_fdr(zone_res, suffix = "_adj", by_col = c("cell_type","zone"))
  fwrite(zone_res, file.path(OUT_DIR, "zone_lmm_results.csv"))
  cat("\n── Part B results (significant or trending) ─────────────────────────────────\n")
  print(zone_res[sig_adj != "ns", .(cell_type, zone, term, beta, pct_change, p, p_adj, sig_adj)])
}

# =============================================================================
# PART C — Zone × treatment interaction (RQ2)
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART C: Zone × treatment interaction (RQ2)\n")
cat("        log_intensity ~ treatment * zone + (1|mouse_id), per cell_type\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

interaction_res <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  d <- section_dt[cell_type == ct]
  cat(sprintf("── Cell type: %s ──\n", ct))
  fit <- fit_lmm(log_intensity ~ treatment * zone + (1 | mouse_id), d,
                 sprintf("Interaction_%s", ct))
  if (is.null(fit)) return(NULL)
  # ANOVA-style test of interaction term
  aov_res <- as.data.frame(anova(fit))
  aov_res$term <- rownames(aov_res)
  setDT(aov_res)
  aov_res[, cell_type := ct]
  aov_res
}), fill = TRUE)

if (!is.null(interaction_res) && nrow(interaction_res) > 0) {
  fwrite(interaction_res, file.path(OUT_DIR, "zone_interaction_anova.csv"))
  cat("\n── Part C: Interaction ANOVA ────────────────────────────────────────────────\n")
  print(interaction_res[, .(cell_type, term, `F value`, `Pr(>F)`)])
}

# =============================================================================
# PART D — Layer-stratified models (RQ5)
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART D: Layer-stratified intensity models (RQ5)\n")
cat("        log_intensity ~ treatment + (1|mouse_id), per layer × cell_type\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

layer_res <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  rbindlist(lapply(LAYER_ORDER, function(ly) {
    d <- section_dt[cell_type == ct & layer == ly]
    if (uniqueN(d$mouse_id) < 3 || nrow(d) < MIN_CELLS) {
      cat(sprintf("  SKIP [%s, %s]: too few observations\n", ct, ly))
      return(NULL)
    }
    cat(sprintf("── %s | %s (%d sections) ──\n", ct, ly, nrow(d)))
    fit <- fit_lmm(log_intensity ~ treatment + (1 | mouse_id), d,
                   sprintf("%s_%s", ct, ly))
    if (is.null(fit)) return(NULL)
    res <- extract_contrasts(fit, model_label = sprintf("layer_%s_%s", ct, ly))
    res[, `:=`(cell_type = ct, layer = ly, pct_change = pct_change(beta))]
    res
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(layer_res) && nrow(layer_res) > 0) {
  layer_res <- add_fdr(layer_res, suffix = "_adj", by_col = c("cell_type","layer"))
  fwrite(layer_res, file.path(OUT_DIR, "layer_lmm_results.csv"))
  cat("\n── Part D results (significant or trending) ─────────────────────────────────\n")
  print(layer_res[sig_adj != "ns", .(cell_type, layer, term, beta, pct_change, p, p_adj, sig_adj)])
}

# =============================================================================
# PART E — Hemisphere comparison (RQ6)
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART E: Hemisphere comparison (RQ6)\n")
cat("        log_intensity ~ treatment * hemisphere + (1|mouse_id), per cell_type\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

# For hemisphere comparison use full dataset (both hemispheres)
dt_both <- dt[str_detect(tolower(brain_area), "visual")]
dt_both[, mouse_id  := str_replace(mouse_id, "C6ST1_ADAMTS4_", "C6ST1_ADAMTS15_")]
dt_both[mouse_id == "mScarlet_4",   treatment := "ADAMTS4_MD"]
dt_both[mouse_id == "ADAMTS4_MD_4", treatment := "mScarlet"]
dt_both <- dt_both[mouse_id != EXCLUDE_PRIMARY]
dt_both[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
dt_both[, hemisphere := factor(hemisphere, levels = c("right","left"),
                                labels = c("contralateral","ipsilateral"))]

section_both <- dt_both[!is.na(normalized_btm20),
  .(log_intensity = mean(log(normalized_btm20[normalized_btm20 > 0]), na.rm = TRUE),
    n_cells = .N),
  by = .(mouse_id, treatment, slice_id, hemisphere, cell_type)]
section_both <- section_both[n_cells >= MIN_CELLS]

hemi_anova <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  d <- section_both[cell_type == ct]
  cat(sprintf("── Cell type: %s ──\n", ct))
  fit <- fit_lmm(log_intensity ~ treatment * hemisphere + (1 | mouse_id), d,
                 sprintf("Hemisphere_%s", ct))
  if (is.null(fit)) return(NULL)
  aov_res <- as.data.frame(anova(fit))
  aov_res$term <- rownames(aov_res)
  setDT(aov_res)
  aov_res[, cell_type := ct]
  aov_res
}), fill = TRUE)

# Ipsilateral contrasts
hemi_ipsi <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  d <- section_both[cell_type == ct & hemisphere == "ipsilateral"]
  fit <- fit_lmm(log_intensity ~ treatment + (1 | mouse_id), d,
                 sprintf("Ipsi_%s", ct))
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, sprintf("ipsi_%s", ct))
  res[, `:=`(cell_type = ct, hemisphere = "ipsilateral",
             pct_change = pct_change(beta))]
  res
}), fill = TRUE)

# Contralateral contrasts
hemi_contra <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  d <- section_both[cell_type == ct & hemisphere == "contralateral"]
  fit <- fit_lmm(log_intensity ~ treatment + (1 | mouse_id), d,
                 sprintf("Contra_%s", ct))
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, sprintf("contra_%s", ct))
  res[, `:=`(cell_type = ct, hemisphere = "contralateral",
             pct_change = pct_change(beta))]
  res
}), fill = TRUE)

if (nrow(hemi_anova) > 0)  fwrite(hemi_anova,  file.path(OUT_DIR, "hemisphere_interaction_anova.csv"))
if (nrow(hemi_ipsi) > 0) {
  hemi_ipsi  <- add_fdr(hemi_ipsi,  suffix = "_adj", by_col = "cell_type")
  fwrite(hemi_ipsi,  file.path(OUT_DIR, "hemisphere_ipsi_contrasts.csv"))
}
if (nrow(hemi_contra) > 0) {
  hemi_contra <- add_fdr(hemi_contra, suffix = "_adj", by_col = "cell_type")
  fwrite(hemi_contra, file.path(OUT_DIR, "hemisphere_contra_contrasts.csv"))
}

cat("\n── Part E: Hemisphere interaction ANOVA ─────────────────────────────────────\n")
print(hemi_anova[, .(cell_type, term, `F value`, `Pr(>F)`)])

# =============================================================================
# PART F — Planned contrasts: RQ3 and RQ4
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART F: Planned contrasts\n")
cat("        RQ3: C6ST1_ADAMTS15 vs ADAMTS15 (synergy)\n")
cat("        RQ4: ADAMTS4_MD vs ADAMTS4 (monocular deprivation)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

planned_res <- rbindlist(lapply(c("PNN","PV"), function(ct) {
  rbindlist(lapply(ZONE_ORDER, function(z) {
    d <- section_dt[cell_type == ct & zone == z]
    if (uniqueN(d$mouse_id) < 3) return(NULL)
    fit <- fit_lmm(log_intensity ~ treatment + (1 | mouse_id), d,
                   sprintf("Planned_%s_%s", ct, z))
    if (is.null(fit)) return(NULL)
    emm <- emmeans(fit, ~ treatment)
    # RQ3: synergy
    synergy <- as.data.table(contrast(emm,
      list("C6ST1_ADAMTS15 - ADAMTS15" = c(
        mScarlet=0, ADAMTS4=0, ADAMTS4_MD=0, ADAMTS15=-1,
        C6ST1=0, C6ST1_ADAMTS15=1))))
    synergy[, `:=`(contrast_label = "C6ST1_ADAMTS15 vs ADAMTS15",
                   rq = "RQ3", cell_type = ct, zone = z)]
    # RQ4: MD
    md <- as.data.table(contrast(emm,
      list("ADAMTS4_MD - ADAMTS4" = c(
        mScarlet=0, ADAMTS4=-1, ADAMTS4_MD=1, ADAMTS15=0,
        C6ST1=0, C6ST1_ADAMTS15=0))))
    md[, `:=`(contrast_label = "ADAMTS4_MD vs ADAMTS4",
              rq = "RQ4", cell_type = ct, zone = z)]
    rbind(synergy, md, fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(planned_res) && nrow(planned_res) > 0) {
  # Standardise column names from emmeans output
  if ("estimate" %in% names(planned_res)) setnames(planned_res, "estimate", "beta")
  if ("p.value"  %in% names(planned_res)) setnames(planned_res, "p.value",  "p")
  planned_res[, pct_change := pct_change(beta)]
  planned_res <- add_fdr(planned_res, suffix = "_adj", by_col = c("rq","cell_type"))
  fwrite(planned_res, file.path(OUT_DIR, "planned_contrasts.csv"))
  cat("\n── Part F results ───────────────────────────────────────────────────────────\n")
  print(planned_res[, .(rq, cell_type, zone, contrast_label, beta, pct_change, p, p_adj, sig_adj)])
}


# =============================================================================
# Script 09 — Part G: Focused layer contrast (RQ5 addendum)
# =============================================================================
# Biologically motivated single contrast:
#   Superficial (layers 1 + 2/3) vs Deep (layers 4 + 5 + 6)
#
# Rationale: PV cells and PNN-PV interactions are densest in layers 2/3 and 4.
# If ADAMTS15 acts specifically on PNNs enwrapping active PV cells, the effect
# should be stronger in these layers than in layer 1 (sparse PV, few PNNs)
# or layers 5/6 (sparser PNN coverage). This is a single pre-specified contrast
# — no BH correction applied.
#
# Input:  cells_with_zones.csv (same as rest of script 09)
# Appends outputs to: results/09_cell_intensity_lmm/
# =============================================================================

library(data.table)
library(lme4)
library(lmerTest)
library(emmeans)

# -----------------------------------------------------------------------------
# Paths — match script 09 constants
# -----------------------------------------------------------------------------
CELLS_PATH <- "/path/to/analysis_results/cells_with_zones.csv"
OUT_DIR    <- "/path/to/results/09_cell_intensity_lmm"

TREATMENT_ORDER <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")

MIN_PV_SECTION <- 5L

# -----------------------------------------------------------------------------
# Load and correct data (same corrections as script 09)
# -----------------------------------------------------------------------------
cat("Loading cells_with_zones.csv...\n")
cells <- fread(CELLS_PATH)

cells[treatment == "C6ST1_ADAMTS4",  treatment := "C6ST1_ADAMTS15"]
cells[mouse_id  == "mScarlet_4",     treatment := "ADAMTS4_MD"]
cells[mouse_id  == "ADAMTS4_MD_4",   treatment := "mScarlet"]
cells <- cells[mouse_id != "C6ST1_ADAMTS15_4"]

cells[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

# Extract layer from brain_area — matches script 09 logic
cells[, layer := fcase(
  brain_area %like% "layer 1$",   "L1",
  brain_area %like% "layer 2",    "L2/3",
  brain_area %like% "layer 4$",   "L4",
  brain_area %like% "layer 5$",   "L5",
  brain_area %like% "layer 6",    "L6",
  default                          = NA_character_
)]
cells <- cells[!is.na(layer)]

# Assign layer group: superficial vs deep
cells[, layer_group := fcase(
  layer %in% c("L1", "L2/3"), "Superficial",
  layer %in% c("L4", "L5", "L6"), "Deep",
  default = NA_character_
)]
cells <- cells[!is.na(layer_group)]
cells[, layer_group := factor(layer_group, levels = c("Superficial", "Deep"))]

# Restrict to ipsilateral, visual cortex, PV cells
cells_pv <- cells[
  hemisphere == "left" &
  brain_area %like% "(?i)visual" &
  cell_type  == "PV"
]

cat(sprintf("PV cells after filtering: %d\n", nrow(cells_pv)))

# Aggregate to section level (mouse_id × slice_id × layer_group)
sec <- cells_pv[, .(
  mean_pv   = mean(normalized_btm20, na.rm = TRUE),
  n_pv      = .N,
  treatment = treatment[1]
), by = .(mouse_id, slice_id, layer_group)]

sec <- sec[n_pv >= MIN_PV_SECTION]
cat(sprintf("Sections after min-cell filter: %d\n", nrow(sec)))

# -----------------------------------------------------------------------------
# Part G1: Does PV intensity differ between superficial and deep layers?
# Model: mean_pv ~ treatment * layer_group + (1|mouse_id)
# Primary interest: treatment:layer_group interaction
# Secondary: treatment contrasts within each layer group
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Part G1: PV intensity ~ treatment * layer_group\n")
cat("══════════════════════════════════════════════════════\n")

m_interact <- tryCatch(
  lmer(mean_pv ~ treatment * layer_group + (1|mouse_id),
       data = sec, REML = FALSE),
  error = function(e) { cat("Interaction model failed:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(m_interact)) {
  if (isSingular(m_interact)) cat("  Singular fit — interpret with caution\n")

  # LRT: does adding layer_group interaction improve fit?
  m_main <- lmer(mean_pv ~ treatment + layer_group + (1|mouse_id),
                 data = sec, REML = FALSE)
  lrt <- anova(m_main, m_interact)
  cat("\nLRT — interaction vs main effects:\n")
  print(lrt)

  fwrite(as.data.table(lrt), file.path(OUT_DIR, "partG_layer_interaction_lrt.csv"))
}

# Part G2: Superficial vs Deep contrast per treatment
# Using interaction model so each treatment gets its own estimate
m_interact_em <- emmeans(m_interact, ~ layer_group | treatment)
ct_layer <- as.data.table(
  contrast(m_interact_em, method = "pairwise", adjust = "none")
)
ct_layer[, sig := ifelse(p.value < 0.05, "*", "")]
cat("\nSuperficial vs Deep per treatment (interaction model, no correction):\n")
print(ct_layer[, .(treatment, contrast, estimate, SE, df, t.ratio, p.value, sig)])

fwrite(ct_layer, file.path(OUT_DIR, "partG_superficial_vs_deep.csv"))

# -----------------------------------------------------------------------------
# Part G3: Key pre-specified contrast for RQ5
# Does ADAMTS15 show a stronger superficial-vs-deep gradient than mScarlet?
# This tests whether the protease effect is layer-specific, not just global.
# Model: mean_pv ~ treatment * layer_group + (1|mouse_id)
# Contrast: (ADAMTS15 Superficial - ADAMTS15 Deep) vs
#           (mScarlet  Superficial - mScarlet  Deep)
# = difference in layer gradients between ADAMTS15 and control
# -----------------------------------------------------------------------------
cat("\n── Part G3: ADAMTS15 layer gradient vs mScarlet layer gradient ──\n")

if (!is.null(m_interact)) {
  em_full <- emmeans(m_interact, ~ treatment * layer_group)

  # Difference-in-differences: does ADAMTS15 show steeper superficial-vs-deep
  # gradient than mScarlet?
  # Positive = ADAMTS15 superficial effect relatively larger
  did_contrast <- contrast(em_full,
    method = list(
      "ADAMTS15 gradient vs mScarlet gradient" = c(
        # order: mScarlet×Superficial, mScarlet×Deep,
        #        ADAMTS4×Sup, ADAMTS4×Deep, ADAMTS4_MD×Sup, ADAMTS4_MD×Deep,
        #        ADAMTS15×Sup, ADAMTS15×Deep,
        #        C6ST1×Sup, C6ST1×Deep, C6ST1_ADAMTS15×Sup, C6ST1_ADAMTS15×Deep
        # mScarlet: +1 Sup, -1 Deep (baseline gradient)
        # ADAMTS15: -1 Sup, +1 Deep (treatment gradient, flipped for difference)
        1, -1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0
      ),
      "C6ST1_ADAMTS15 gradient vs mScarlet gradient" = c(
        1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1
      )
    ),
    adjust = "none"
  )

  did_dt <- as.data.table(summary(did_contrast))
  did_dt[, sig := ifelse(p.value < 0.05, "*", "")]
  cat("\nDifference-in-differences: layer gradient contrast\n")
  cat("(Negative = treatment shows shallower Sup-Deep gradient = effect in both layers)\n")
  cat("(Positive = treatment shows steeper Sup-Deep gradient = concentrated in one layer)\n")
  print(did_dt[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

  fwrite(did_dt, file.path(OUT_DIR, "partG_layer_gradient_did.csv"))
}

# -----------------------------------------------------------------------------
# Part G4: Animal-level summary by layer group for plotting
# -----------------------------------------------------------------------------
animal_layer <- sec[, .(
  mean_pv   = mean(mean_pv, na.rm = TRUE),
  n_sec     = .N,
  treatment = treatment[1]
), by = .(mouse_id, layer_group)]

fwrite(animal_layer, file.path(OUT_DIR, "partG_animal_layer_summary.csv"))

# Quick descriptive print
cat("\nMean PV intensity by treatment × layer group:\n")
print(animal_layer[, .(mean_pv = round(mean(mean_pv), 3)),
                   by = .(treatment, layer_group)][order(treatment, layer_group)])

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 09 Part G complete\n")
cat(sprintf("Outputs appended to: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  partG_layer_interaction_lrt.csv\n")
cat("  partG_superficial_vs_deep.csv\n")
cat("  partG_layer_gradient_did.csv\n")
cat("  partG_animal_layer_summary.csv\n")
cat("══════════════════════════════════════════════════════\n")


# =============================================================================
# FIGURES
# =============================================================================
cat("\n── Generating figures ───────────────────────────────────────────────────────\n")

# Figure 1: Primary model — % change from mScarlet, PNN and PV side by side
if (exists("primary_res") && !is.null(primary_res) && nrow(primary_res) > 0) {
  pd <- primary_res[str_detect(term, "treatment")]
  pd[, treatment := str_remove(term, "treatment")]
  pd[, treatment := factor(treatment, levels = TREATMENT_ORDER[-1])]

  p1 <- ggplot(pd, aes(x = treatment, y = pct_change, fill = treatment)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pct_change - se*100, ymax = pct_change + se*100),
                  width = 0.25) +
    geom_text(aes(label = sig_adj,
                  y = ifelse(pct_change >= 0,
                             pct_change + se*100 + 1,
                             pct_change - se*100 - 1)),
              size = 4, fontface = "bold") +
    facet_wrap(~cell_type) +
    scale_fill_manual(values = PALETTE[-1]) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(title = "Cell intensity: % change from mScarlet (primary model)",
         x = NULL, y = "% change in normalised intensity") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(OUT_DIR, "fig1_primary_pct_change.pdf"),
         p1, width = 8, height = 4)
  cat("  Saved fig1_primary_pct_change.pdf\n")
}

# Figure 2: Zone gradient — mean intensity by zone per treatment, PNN only
zone_means <- section_dt[cell_type == "PNN",
  .(mean_int = mean(log_intensity, na.rm = TRUE),
    se_int   = sd(log_intensity, na.rm = TRUE) / sqrt(.N)),
  by = .(treatment, zone)]
zone_means[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
zone_means[, zone := factor(zone, levels = ZONE_ORDER)]

p2 <- ggplot(zone_means, aes(x = zone, y = mean_int,
                              colour = treatment, group = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_int - se_int, ymax = mean_int + se_int),
                width = 0.15) +
  scale_colour_manual(values = PALETTE) +
  labs(title = "PNN cell intensity by injection zone (mean ± SEM)",
       x = "Zone", y = "Mean log normalised intensity",
       colour = "Treatment") +
  theme_classic(base_size = 11)

ggsave(file.path(OUT_DIR, "fig2_zone_gradient.pdf"),
       p2, width = 7, height = 4)
cat("  Saved fig2_zone_gradient.pdf\n")

# Figure 3: Layer profile — heatmap of % change per layer × treatment, PNN
if (exists("layer_res") && !is.null(layer_res) && nrow(layer_res) > 0) {
  ld <- layer_res[cell_type == "PNN" & str_detect(term, "treatment")]
  ld[, treatment := str_remove(term, "treatment")]
  ld[, treatment := factor(treatment, levels = TREATMENT_ORDER[-1])]
  ld[, layer := factor(layer, levels = LAYER_ORDER)]

  p3 <- ggplot(ld, aes(x = layer, y = treatment, fill = pct_change)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_adj), size = 3.5) +
    scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#E84C4C",
                         midpoint = 0, name = "% change") +
    labs(title = "PNN intensity: % change from mScarlet by cortical layer",
         x = "Layer", y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(OUT_DIR, "fig3_layer_heatmap.pdf"),
         p3, width = 7, height = 4)
  cat("  Saved fig3_layer_heatmap.pdf\n")
}

cat("\n── Script 09 complete ───────────────────────────────────────────────────────\n")
cat(sprintf("All outputs saved to: %s\n", OUT_DIR))
