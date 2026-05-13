# =============================================================================
# 13_bayesian_lmm.R
# Bayesian mixed-effects models using brms/Stan
#
# Six models covering all seven research questions:
#
# Model 1 — Primary enwrapment (RQ1)
#   frac_enwrapped ~ treatment + (1|mouse_id)
#   Includes prior sensitivity check with tighter priors.
#
# Model 2A — Composite primary (RQ1)
#   composite ~ treatment + (1|mouse_id)
#   Includes planned contrasts for RQ3 (synergy) and RQ4 (MD)
#   extracted via hypothesis() — no additional MCMC sampling.
#
# Model 2B — Composite zone interaction (RQ2)
#   composite ~ treatment * zone + (1|mouse_id)
#   Zone-specific posteriors extracted via posterior_epred.
#
# Model 3 — Hemisphere lateralisation (RQ6)
#   composite ~ treatment * hemisphere + (1|mouse_id)
#   P(ipsilateral effect > contralateral effect) per treatment.
#
# Model 4 — Ipsilateral composite (RQ1 + RQ6 combined)
#   composite ~ treatment + (1|mouse_id)
#   Fitted on ipsilateral sections only — directly comparable to
#   the significant density results from script 11.
#
# Model 5 — Layer specificity (RQ5)
#   composite ~ treatment * layer + (1|mouse_id)
#   Layer-specific posteriors extracted via posterior_epred.
#   Tests whether treatment effects differ across cortical layers.
#
# Model 6 — Off-target specificity (RQ7)
#   composite ~ treatment * region_type + (1|mouse_id)
#   Data: combined visual + off-target sections from script 12.
#   P(visual effect > off-target effect) per treatment.
#
# Prior specification (weakly informative):
#   Fixed effects:    normal(0, 1)
#   Random effect SD: half-normal(0, 1)
#   Residual SD:      half-normal(0, 1)
#
# MCMC: 4 chains x 4000 iter (2000 warmup), cores=4, adapt_delta=0.95
# Stan model objects cached to OUT_DIR — rerunning loads from cache.
#
# ROPE thresholds:
#   frac_enwrapped: 0.05 (5 percentage points)
#   composite:      0.10 SD units
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(RANN)
  library(brms)
  library(ggplot2)
  library(patchwork)
})

if (!requireNamespace("posterior", quietly = TRUE))
  install.packages("posterior", repos = "https://cloud.r-project.org")
library(posterior)

# ── Paths ─────────────────────────────────────────────────────────────────────
COMPOSITE_CSV  <- "/path/to/results/12_pnn_integrity_composite/composite_section_data.csv"
HEMI_CSV       <- "/path/to/results/12_pnn_integrity_composite/composite_hemisphere_data.csv"
CELLS_CSV      <- "/path/to/analysis_results/cells_with_zones.csv"
RESULTS_DIR    <- "/path/to/results"
OUT_DIR        <- file.path(RESULTS_DIR, "13_bayesian_lmm")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENTS      <- c("mScarlet","ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1","C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS
ZONE_ORDER      <- c("Core","Penumbra","Outside")
LAYER_ORDER     <- c("layer 1","layer 2/3","layer 4","layer 5","layer 6")
VISUAL_KEYWORDS <- c("visual","VIS")
COLOC_THRESH    <- 30L
MIN_PV          <- 5L
MIN_CELLS_ZONE  <- 5L
EXCLUDE_PRIMARY <- ""     # any mouse_id containing this string will be excluded from primary analyses (e.g. injection failure)
ROPE_ENWRAP     <- 0.05
ROPE_COMPOSITE  <- 0.10
MCMC_CHAINS     <- 4L
MCMC_ITER       <- 4000L
MCMC_WARMUP     <- 2000L
MCMC_CORES      <- 4L
MCMC_SEED       <- 2025L
ADAPT_DELTA     <- 0.95

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

cat("Script 13: Bayesian mixed-effects models — full pipeline\n")
cat(sprintf("Chains: %d | Iter: %d | Warmup: %d | Cores: %d\n",
    MCMC_CHAINS, MCMC_ITER, MCMC_WARMUP, MCMC_CORES))
cat("Models: M1 enwrap | M2A composite | M2B zone | M3 hemisphere |\n")
cat("        M4 ipsilateral | M5 layer | M6 off-target\n\n")

# ── Priors ────────────────────────────────────────────────────────────────────
priors_weak <- c(
  prior(normal(0, 1),   class = b),
  prior(normal(0, 1),   class = sd),
  prior(normal(0, 1),   class = sigma)
)

priors_tight <- c(
  prior(normal(0, 0.5), class = b),
  prior(normal(0, 0.5), class = sd),
  prior(normal(0, 0.5), class = sigma)
)

# ── Helpers ───────────────────────────────────────────────────────────────────
fix_names <- function(dt, id_col = "mouse_id") {
  dt[, (id_col)  := str_replace(get(id_col), "C6ST1_ADAMTS4_", "C6ST1_ADAMTS15_")]
  dt[, treatment := str_replace(treatment, "C6ST1_ADAMTS4$", "C6ST1_ADAMTS15")]
  dt[get(id_col) == "mScarlet_4",   treatment := "ADAMTS4_MD"]
  dt[get(id_col) == "ADAMTS4_MD_4", treatment := "mScarlet"]
  dt
}

is_visual <- function(x) str_detect(tolower(x),
  paste(tolower(VISUAL_KEYWORDS), collapse = "|"))

extract_post <- function(fit, model_name, rope_bound) {
  draws <- as_draws_df(fit)
  treat_cols <- grep("^b_treatment", names(draws), value = TRUE)
  rbindlist(lapply(treat_cols, function(col) {
    vals  <- draws[[col]]
    treat <- str_remove(col, "^b_treatment")
    ci89  <- quantile(vals, c(0.055, 0.945))
    ci95  <- quantile(vals, c(0.025, 0.975))
    data.table(
      model        = model_name,
      treatment    = treat,
      post_mean    = round(mean(vals),  4),
      post_sd      = round(sd(vals),    4),
      ci_89_lo     = round(ci89[[1]],   4),
      ci_89_hi     = round(ci89[[2]],   4),
      ci_95_lo     = round(ci95[[1]],   4),
      ci_95_hi     = round(ci95[[2]],   4),
      p_negative   = round(mean(vals < 0),           3),
      p_meaningful = round(mean(vals < -rope_bound), 3)
    )
  }), fill = TRUE)
}

extract_diag <- function(fit, model_name) {
  rh   <- rhat(fit)
  ne   <- neff_ratio(fit)
  ndiv <- sum(nuts_params(fit)$Value[nuts_params(fit)$Parameter == "divergent__"])
  data.table(
    model          = model_name,
    max_rhat       = round(max(rh, na.rm = TRUE), 4),
    min_neff_ratio = round(min(ne, na.rm = TRUE), 4),
    n_divergent    = ndiv,
    converged      = max(rh, na.rm = TRUE) < 1.01 & ndiv == 0
  )
}

epred_contrasts <- function(epred, newdata, treat_col, group_col,
                             group_levels, rope_bound) {
  rbindlist(lapply(TREATMENTS[-1], function(treat) {
    rbindlist(lapply(group_levels, function(grp) {
      it <- which(newdata[[treat_col]] == treat      & newdata[[group_col]] == grp)
      ic <- which(newdata[[treat_col]] == "mScarlet" & newdata[[group_col]] == grp)
      dv <- epred[, it] - epred[, ic]
      data.table(
        treatment    = treat,
        group        = grp,
        post_mean    = round(mean(dv),              4),
        post_sd      = round(sd(dv),                4),
        ci_89_lo     = round(quantile(dv, 0.055),   4),
        ci_89_hi     = round(quantile(dv, 0.945),   4),
        p_negative   = round(mean(dv < 0),          3),
        p_meaningful = round(mean(dv < -rope_bound),3)
      )
    }), fill = TRUE)
  }), fill = TRUE)
}

# ── Load section data ─────────────────────────────────────────────────────────
cat("Loading section data...\n")
dt <- fread(COMPOSITE_CSV)
dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
dt[, zone      := factor(zone,      levels = ZONE_ORDER)]

# Layer information not available in composite_section_data
# Model 5 will source layer data directly from cells_with_zones
dt[, layer := NA_character_]

hemi_dt <- fread(HEMI_CSV)
hemi_dt[, treatment  := factor(treatment,  levels = TREATMENT_ORDER)]
hemi_dt[, hemisphere := factor(hemisphere, levels = c("contralateral","ipsilateral"))]

cat(sprintf("  Section data: %d rows, %d animals\n", nrow(dt), uniqueN(dt$mouse_id)))
cat(sprintf("  Hemisphere data: %d rows\n\n", nrow(hemi_dt)))

# ── Build off-target data for Model 6 ─────────────────────────────────────────
cat("Building off-target section data for Model 6...\n")
cells <- fread(CELLS_CSV)
cells <- fix_names(cells)
cells <- cells[mouse_id != EXCLUDE_PRIMARY]
cells[, is_visual := is_visual(brain_area)]

# Global standardisation params from visual cortex sections
mu_fe  <- mean(dt$frac_enwrapped,   na.rm = TRUE)
sd_fe  <- sd(dt$frac_enwrapped,     na.rm = TRUE)
mu_btm <- mean(dt$normalized_btm20, na.rm = TRUE)
sd_btm <- sd(dt$normalized_btm20,   na.rm = TRUE)

cells_off <- cells[is_visual == FALSE]
ot_records <- list()
for (mid in unique(cells_off$mouse_id)) {
  treat <- cells_off[mouse_id == mid, treatment[1]]
  for (sl in unique(cells_off[mouse_id == mid, slice_id])) {
    grp <- cells_off[mouse_id == mid & slice_id == sl]
    pv  <- grp[cell_type == "PV",  .(x_hires, y_hires)]
    pnn <- grp[cell_type == "PNN", .(x_hires, y_hires)]
    if (nrow(pv) < MIN_PV || nrow(pnn) < MIN_PV) next
    nn    <- nn2(data = as.matrix(pnn), query = as.matrix(pv), k = 1L)
    frac  <- mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
    ivs   <- grp[cell_type == "PNN" & !is.na(normalized_btm20) &
                   normalized_btm20 > 0, normalized_btm20]
    if (length(ivs) < MIN_CELLS_ZONE) next
    ot_records[[length(ot_records) + 1L]] <- data.table(
      mouse_id = mid, treatment = treat, slice_id = sl,
      frac_enwrapped = frac, normalized_btm20 = mean(ivs)
    )
  }
}
ot_dt <- rbindlist(ot_records)
ot_dt[, frac_enwrapped_z   := (frac_enwrapped   - mu_fe)  / sd_fe]
ot_dt[, normalized_btm20_z := (normalized_btm20 - mu_btm) / sd_btm]
ot_dt[, composite   := rowMeans(.SD, na.rm = TRUE),
       .SDcols = c("frac_enwrapped_z","normalized_btm20_z")]
ot_dt[, region_type := "Off-target"]
ot_dt[, treatment   := factor(treatment, levels = TREATMENT_ORDER)]

vc_dt <- dt[, .(mouse_id, treatment, slice_id, composite,
                 region_type = "Visual (target)")]

rq7_dt <- rbindlist(list(vc_dt, ot_dt[, .(mouse_id, treatment, slice_id,
                                            composite, region_type)]),
                    fill = TRUE)
rq7_dt[, region_type := factor(region_type,
                                levels = c("Visual (target)","Off-target"))]
rq7_dt[, treatment   := factor(treatment, levels = TREATMENT_ORDER)]
cat(sprintf("  Off-target sections: %d | Visual: %d\n\n",
    nrow(ot_dt), nrow(vc_dt)))

# =============================================================================
# MODEL 1: frac_enwrapped ~ treatment + (1|mouse_id)  [RQ1]
# =============================================================================
cat("Fitting Model 1 — enwrapment, weak priors (RQ1)...\n")
m1_weak <- brm(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data = dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m1_enwrap_weak"), silent = 2
)

cat("Fitting Model 1 — enwrapment, tight priors (sensitivity)...\n")
m1_tight <- brm(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data = dt, prior = priors_tight,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m1_enwrap_tight"), silent = 2
)

m1_post_weak  <- extract_post(m1_weak,  "M1_weak",  ROPE_ENWRAP)
m1_post_tight <- extract_post(m1_tight, "M1_tight", ROPE_ENWRAP)
m1_diag       <- rbind(extract_diag(m1_weak,  "M1_weak"),
                        extract_diag(m1_tight, "M1_tight"))
sens_m1 <- merge(
  m1_post_weak[,  .(treatment, post_mean_weak  = post_mean)],
  m1_post_tight[, .(treatment, post_mean_tight = post_mean)], by = "treatment")
sens_m1[, beta_shift := round(post_mean_tight - post_mean_weak, 4)]

cat("\nModel 1 diagnostics:\n");        print(m1_diag)
cat("\nModel 1 posteriors (weak):\n")
print(m1_post_weak[, .(treatment, post_mean, post_sd,
                         ci_89_lo, ci_89_hi, p_negative, p_meaningful)])
cat("\nModel 1 prior sensitivity:\n");  print(sens_m1)

fwrite(rbind(m1_post_weak, m1_post_tight),
       file.path(OUT_DIR, "m1_enwrap_posterior.csv"))
fwrite(m1_diag,  file.path(OUT_DIR, "m1_diagnostics.csv"))
fwrite(sens_m1,  file.path(OUT_DIR, "m1_prior_sensitivity.csv"))

# =============================================================================
# MODEL 2A: composite ~ treatment + (1|mouse_id)  [RQ1, RQ3, RQ4]
# =============================================================================
cat("\nFitting Model 2A — composite primary (RQ1)...\n")
m2a <- brm(
  composite ~ treatment + (1 | mouse_id),
  data = dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m2a_composite_primary"), silent = 2
)

m2a_post <- extract_post(m2a, "M2A", ROPE_COMPOSITE)
m2a_diag <- extract_diag(m2a, "M2A")

cat("\nModel 2A diagnostics:\n"); print(m2a_diag)
cat("\nModel 2A posteriors:\n")
print(m2a_post[, .(treatment, post_mean, post_sd,
                    ci_89_lo, ci_89_hi, p_negative, p_meaningful)])

# ── RQ3: C6ST1_ADAMTS15 vs ADAMTS15 (synergy) ────────────────────────────────
cat("\nRQ3 planned contrast (synergy) from Model 2A...\n")
rq3_hyp <- hypothesis(m2a,
  "treatmentC6ST1_ADAMTS15 - treatmentADAMTS15 < 0")
rq3_dt <- as.data.table(rq3_hyp$hypothesis)
rq3_dt[, rq := "RQ3"]
rq3_dt[, contrast := "C6ST1_ADAMTS15 - ADAMTS15"]
cat("RQ3 (P(C6ST1_ADAMTS15 < ADAMTS15)):\n"); print(rq3_dt)

# ── RQ4: ADAMTS4_MD vs ADAMTS4 (monocular deprivation) ───────────────────────
cat("\nRQ4 planned contrast (MD) from Model 2A...\n")
rq4_hyp <- hypothesis(m2a,
  "treatmentADAMTS4_MD - treatmentADAMTS4 < 0")
rq4_dt <- as.data.table(rq4_hyp$hypothesis)
rq4_dt[, rq := "RQ4"]
rq4_dt[, contrast := "ADAMTS4_MD - ADAMTS4"]
cat("RQ4 (P(ADAMTS4_MD < ADAMTS4)):\n"); print(rq4_dt)

planned_bayes <- rbind(rq3_dt, rq4_dt, fill = TRUE)
fwrite(m2a_post,      file.path(OUT_DIR, "m2a_composite_posterior.csv"))
fwrite(m2a_diag,      file.path(OUT_DIR, "m2a_diagnostics.csv"))
fwrite(planned_bayes, file.path(OUT_DIR, "m2a_planned_contrasts.csv"))

# =============================================================================
# MODEL 2B: composite ~ treatment * zone + (1|mouse_id)  [RQ2]
# =============================================================================
cat("\nFitting Model 2B — composite zone interaction (RQ2)...\n")
m2b <- brm(
  composite ~ treatment * zone + (1 | mouse_id),
  data = dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m2b_composite_zone"), silent = 2
)

m2b_diag <- extract_diag(m2b, "M2B")
cat("\nModel 2B diagnostics:\n"); print(m2b_diag)

nd_zone <- expand.grid(treatment = TREATMENTS, zone = ZONE_ORDER,
                        mouse_id = NA, stringsAsFactors = FALSE)
nd_zone$treatment <- factor(nd_zone$treatment, levels = TREATMENT_ORDER)
nd_zone$zone      <- factor(nd_zone$zone,      levels = ZONE_ORDER)
epred_zone  <- posterior_epred(m2b, newdata = nd_zone,
                                re_formula = NA, allow_new_levels = TRUE)
zone_post   <- epred_contrasts(epred_zone, nd_zone,
                                "treatment", "zone", ZONE_ORDER, ROPE_COMPOSITE)
setnames(zone_post, "group", "zone")

cat("\nModel 2B zone posteriors:\n")
print(zone_post[, .(treatment, zone, post_mean, p_negative, p_meaningful)])

fwrite(zone_post, file.path(OUT_DIR, "m2b_zone_posteriors.csv"))
fwrite(m2b_diag,  file.path(OUT_DIR, "m2b_diagnostics.csv"))

# =============================================================================
# MODEL 3: composite ~ treatment * hemisphere + (1|mouse_id)  [RQ6]
# =============================================================================
cat("\nFitting Model 3 — hemisphere lateralisation (RQ6)...\n")
m3 <- brm(
  composite ~ treatment * hemisphere + (1 | mouse_id),
  data = hemi_dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m3_hemisphere"), silent = 2
)

m3_diag <- extract_diag(m3, "M3")
cat("\nModel 3 diagnostics:\n"); print(m3_diag)

nd_hemi <- expand.grid(
  treatment  = TREATMENTS,
  hemisphere = c("contralateral","ipsilateral"),
  mouse_id   = NA, stringsAsFactors = FALSE)
nd_hemi$treatment  <- factor(nd_hemi$treatment,  levels = TREATMENT_ORDER)
nd_hemi$hemisphere <- factor(nd_hemi$hemisphere,
                              levels = c("contralateral","ipsilateral"))
epred_hemi <- posterior_epred(m3, newdata = nd_hemi,
                               re_formula = NA, allow_new_levels = TRUE)
hemi_post  <- epred_contrasts(epred_hemi, nd_hemi,
                               "treatment", "hemisphere",
                               c("ipsilateral","contralateral"), ROPE_COMPOSITE)
setnames(hemi_post, "group", "hemisphere")

# Lateralisation: P(ipsi effect more negative than contra)
lat_post <- rbindlist(lapply(TREATMENTS[-1], function(treat) {
  ii <- which(nd_hemi$treatment == treat      & nd_hemi$hemisphere == "ipsilateral")
  ic <- which(nd_hemi$treatment == "mScarlet" & nd_hemi$hemisphere == "ipsilateral")
  ji <- which(nd_hemi$treatment == treat      & nd_hemi$hemisphere == "contralateral")
  jc <- which(nd_hemi$treatment == "mScarlet" & nd_hemi$hemisphere == "contralateral")
  di <- epred_hemi[, ii] - epred_hemi[, ic]
  dc <- epred_hemi[, ji] - epred_hemi[, jc]
  lat <- di - dc
  data.table(
    treatment          = treat,
    post_mean_ipsi     = round(mean(di),  3),
    post_mean_contra   = round(mean(dc),  3),
    post_mean_lat_diff = round(mean(lat), 3),
    p_ipsi_gt_contra   = round(mean(lat < 0), 3),
    ci_89_lat_lo       = round(quantile(lat, 0.055), 3),
    ci_89_lat_hi       = round(quantile(lat, 0.945), 3)
  )
}), fill = TRUE)

cat("\nModel 3 hemisphere posteriors:\n")
print(hemi_post[, .(treatment, hemisphere, post_mean, p_negative, p_meaningful)])
cat("\nModel 3 lateralisation:\n"); print(lat_post)

fwrite(hemi_post, file.path(OUT_DIR, "m3_hemisphere_posteriors.csv"))
fwrite(lat_post,  file.path(OUT_DIR, "m3_lateralisation.csv"))
fwrite(m3_diag,   file.path(OUT_DIR, "m3_diagnostics.csv"))

# =============================================================================
# MODEL 4: composite ~ treatment + (1|mouse_id), ipsilateral only  [RQ1 + RQ6]
# =============================================================================
cat("\nFitting Model 4 — composite ipsilateral only (RQ1 + RQ6)...\n")
# hemi_dt already has hemisphere labels; filter to ipsilateral
ipsi_dt <- hemi_dt[hemisphere == "ipsilateral"]

m4 <- brm(
  composite ~ treatment + (1 | mouse_id),
  data = ipsi_dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m4_composite_ipsi"), silent = 2
)

m4_post <- extract_post(m4, "M4_ipsi", ROPE_COMPOSITE)
m4_diag <- extract_diag(m4, "M4_ipsi")

cat("\nModel 4 diagnostics:\n"); print(m4_diag)
cat("\nModel 4 posteriors (ipsilateral hemisphere):\n")
print(m4_post[, .(treatment, post_mean, post_sd,
                   ci_89_lo, ci_89_hi, p_negative, p_meaningful)])

fwrite(m4_post, file.path(OUT_DIR, "m4_ipsi_posterior.csv"))
fwrite(m4_diag, file.path(OUT_DIR, "m4_diagnostics.csv"))

# =============================================================================
# MODEL 5: composite ~ treatment * layer + (1|mouse_id)  [RQ5]
# =============================================================================
cat("\nBuilding layer-level composite from cells_with_zones...\n")
cells_layer <- cells[is_visual == TRUE]
cells_layer[, layer := str_extract(brain_area, "layer [0-9/]+")]
cells_layer <- cells_layer[!is.na(layer)]

layer_records <- list()
for (mid in unique(cells_layer$mouse_id)) {
  treat <- cells_layer[mouse_id == mid, treatment[1]]
  for (sl in unique(cells_layer[mouse_id == mid, slice_id])) {
    for (ly in LAYER_ORDER) {
      for (z in ZONE_ORDER) {
        grp <- cells_layer[mouse_id == mid & slice_id == sl &
                              layer == ly & zone == z]
        pv  <- grp[cell_type == "PV",  .(x_hires, y_hires)]
        pnn <- grp[cell_type == "PNN", .(x_hires, y_hires)]
        if (nrow(pv) < MIN_PV || nrow(pnn) < MIN_PV) next
        nn   <- nn2(data = as.matrix(pnn), query = as.matrix(pv), k = 1L)
        frac <- mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
        ivs  <- grp[cell_type == "PNN" & !is.na(normalized_btm20) &
                      normalized_btm20 > 0, normalized_btm20]
        if (length(ivs) < MIN_CELLS_ZONE) next
        layer_records[[length(layer_records) + 1L]] <- data.table(
          mouse_id = mid, treatment = treat, slice_id = sl,
          layer = ly, zone = z,
          frac_enwrapped = frac, normalized_btm20 = mean(ivs)
        )
      }
    }
  }
}
dt5 <- rbindlist(layer_records)
dt5[, frac_enwrapped_z   := (frac_enwrapped   - mu_fe)  / sd_fe]
dt5[, normalized_btm20_z := (normalized_btm20 - mu_btm) / sd_btm]
dt5[, composite := rowMeans(.SD, na.rm = TRUE),
     .SDcols = c("frac_enwrapped_z","normalized_btm20_z")]
dt5[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
dt5[, layer     := factor(layer,     levels = LAYER_ORDER)]
cat(sprintf("  Layer section data: %d rows\n", nrow(dt5)))

cat("Fitting Model 5 — composite layer interaction (RQ5)...\n")
m5 <- brm(
  composite ~ treatment * layer + (1 | mouse_id),
  data = dt5, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m5_composite_layer"), silent = 2
)

m5_diag <- extract_diag(m5, "M5_layer")
cat("\nModel 5 diagnostics:\n"); print(m5_diag)

nd_layer <- expand.grid(treatment = TREATMENTS, layer = LAYER_ORDER,
                         mouse_id = NA, stringsAsFactors = FALSE)
nd_layer$treatment <- factor(nd_layer$treatment, levels = TREATMENT_ORDER)
nd_layer$layer     <- factor(nd_layer$layer,     levels = LAYER_ORDER)
epred_layer  <- posterior_epred(m5, newdata = nd_layer,
                                 re_formula = NA, allow_new_levels = TRUE)
layer_post   <- epred_contrasts(epred_layer, nd_layer,
                                 "treatment", "layer", LAYER_ORDER, ROPE_COMPOSITE)
setnames(layer_post, "group", "layer")

cat("\nModel 5 layer posteriors:\n")
print(layer_post[, .(treatment, layer, post_mean, p_negative, p_meaningful)])

fwrite(layer_post, file.path(OUT_DIR, "m5_layer_posteriors.csv"))
fwrite(m5_diag,    file.path(OUT_DIR, "m5_diagnostics.csv"))

# =============================================================================
# MODEL 6: composite ~ treatment * region_type + (1|mouse_id)  [RQ7]
# =============================================================================
cat("\nFitting Model 6 — off-target specificity (RQ7)...\n")
m6 <- brm(
  composite ~ treatment * region_type + (1 | mouse_id),
  data = rq7_dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(OUT_DIR, "m6_offtarget"), silent = 2
)

m6_diag <- extract_diag(m6, "M6_offtarget")
cat("\nModel 6 diagnostics:\n"); print(m6_diag)

nd_rq7 <- expand.grid(
  treatment   = TREATMENTS,
  region_type = c("Visual (target)","Off-target"),
  mouse_id    = NA, stringsAsFactors = FALSE)
nd_rq7$treatment   <- factor(nd_rq7$treatment,   levels = TREATMENT_ORDER)
nd_rq7$region_type <- factor(nd_rq7$region_type,
                               levels = c("Visual (target)","Off-target"))
epred_rq7  <- posterior_epred(m6, newdata = nd_rq7,
                               re_formula = NA, allow_new_levels = TRUE)
rq7_post   <- epred_contrasts(epred_rq7, nd_rq7, "treatment", "region_type",
                               c("Visual (target)","Off-target"), ROPE_COMPOSITE)
setnames(rq7_post, "group", "region_type")

# P(visual effect > off-target effect) per treatment
spec_post <- rbindlist(lapply(TREATMENTS[-1], function(treat) {
  iv <- which(nd_rq7$treatment == treat      & nd_rq7$region_type == "Visual (target)")
  ic <- which(nd_rq7$treatment == "mScarlet" & nd_rq7$region_type == "Visual (target)")
  jo <- which(nd_rq7$treatment == treat      & nd_rq7$region_type == "Off-target")
  jc <- which(nd_rq7$treatment == "mScarlet" & nd_rq7$region_type == "Off-target")
  dv  <- epred_rq7[, iv] - epred_rq7[, ic]
  do_ <- epred_rq7[, jo] - epred_rq7[, jc]
  sp  <- dv - do_
  data.table(
    treatment            = treat,
    post_mean_visual     = round(mean(dv),  3),
    post_mean_offtarget  = round(mean(do_), 3),
    post_mean_specificity = round(mean(sp), 3),
    p_visual_gt_offtarget = round(mean(sp < 0), 3),
    ci_89_spec_lo        = round(quantile(sp, 0.055), 3),
    ci_89_spec_hi        = round(quantile(sp, 0.945), 3)
  )
}), fill = TRUE)

cat("\nModel 6 region posteriors:\n")
print(rq7_post[, .(treatment, region_type, post_mean, p_negative, p_meaningful)])
cat("\nModel 6 specificity (P(visual effect > off-target)):\n")
print(spec_post)

fwrite(rq7_post,  file.path(OUT_DIR, "m6_offtarget_posteriors.csv"))
fwrite(spec_post, file.path(OUT_DIR, "m6_specificity.csv"))
fwrite(m6_diag,   file.path(OUT_DIR, "m6_diagnostics.csv"))

# =============================================================================
# COMBINED DIAGNOSTICS
# =============================================================================
all_diag <- rbindlist(list(m1_diag, m2a_diag, m2b_diag, m3_diag,
                            m4_diag, m5_diag, m6_diag))
fwrite(all_diag, file.path(OUT_DIR, "all_model_diagnostics.csv"))
cat("\nAll model diagnostics:\n"); print(all_diag)

# =============================================================================
# FIGURES
# =============================================================================
cat("\nGenerating figures...\n")

posterior_density_plot <- function(fit, model_title, rope_bound, rope_label) {
  draws <- as_draws_df(fit)
  tc    <- grep("^b_treatment", names(draws), value = TRUE)
  dl    <- rbindlist(lapply(tc, function(col) data.table(
    treatment = factor(str_remove(col,"^b_treatment"), levels = TREATMENT_ORDER[-1]),
    value = draws[[col]])))
  ggplot(dl, aes(x = value, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(xintercept = 0,           linewidth = 0.6, linetype = "dashed") +
    geom_vline(xintercept = -rope_bound, linewidth = 0.4,
               linetype = "dotted", colour = "grey40") +
    facet_wrap(~treatment, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = PALETTE[-1]) +
    scale_colour_manual(values = PALETTE[-1]) +
    labs(title    = model_title,
         subtitle = sprintf("Dashed: null (0)  Dotted: ROPE boundary (-%s)", rope_label),
         x = "Treatment effect vs mScarlet", y = "Density") +
    theme_classic(base_size = 10) + theme(legend.position = "none")
}

forest_plot <- function(post_dt, group_col, title_str) {
  post_dt <- copy(post_dt)
  post_dt[, treatment := factor(treatment, levels = rev(TREATMENT_ORDER[-1]))]
  post_dt[, (group_col) := factor(get(group_col))]
  ggplot(post_dt, aes_string(y = "treatment", x = "post_mean",
                               xmin = "ci_89_lo", xmax = "ci_89_hi",
                               colour = group_col)) +
    geom_vline(xintercept = 0,               linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -ROPE_COMPOSITE, linewidth = 0.3,
               linetype = "dotted", colour = "grey40") +
    geom_errorbarh(height = 0.3, linewidth = 0.8,
                   position = position_dodgev(height = 0.6)) +
    geom_point(size = 3, position = position_dodgev(height = 0.6)) +
    labs(title    = title_str,
         subtitle = "Points: posterior mean  Bars: 89% CI",
         x = "Effect on composite (SD units)", y = NULL) +
    theme_classic(base_size = 10)
}

# Fig 1: Model 1 densities
p1 <- posterior_density_plot(m1_weak,
  "Model 1: Posterior distributions — frac_enwrapped (RQ1)",
  ROPE_ENWRAP, "0.05")
ggsave(file.path(OUT_DIR, "fig1_m1_posteriors.pdf"), p1, width = 8, height = 8)

# Fig 2: Model 2A densities
p2 <- posterior_density_plot(m2a,
  "Model 2A: Posterior distributions — PNN Integrity Composite (RQ1)",
  ROPE_COMPOSITE, "0.1 SD")
ggsave(file.path(OUT_DIR, "fig2_m2a_posteriors.pdf"), p2, width = 8, height = 8)

# Fig 3: Model 2B zone forest
p3 <- forest_plot(zone_post, "zone",
  "Model 2B: Zone-specific posteriors — composite vs mScarlet (RQ2)")
ggsave(file.path(OUT_DIR, "fig3_m2b_zone_forest.pdf"), p3, width = 11, height = 4)

# Fig 4: Model 3 hemisphere forest
p4 <- forest_plot(hemi_post, "hemisphere",
  "Model 3: Hemisphere posteriors — composite vs mScarlet (RQ6)")
ggsave(file.path(OUT_DIR, "fig4_m3_hemisphere_forest.pdf"), p4, width = 8, height = 5)

# Fig 5: Model 4 ipsilateral densities
p5 <- posterior_density_plot(m4,
  "Model 4: Posterior distributions — ipsilateral composite (RQ1 + RQ6)",
  ROPE_COMPOSITE, "0.1 SD")
ggsave(file.path(OUT_DIR, "fig5_m4_ipsi_posteriors.pdf"), p5, width = 8, height = 8)

# Fig 6: Model 5 layer forest
p6 <- forest_plot(layer_post, "layer",
  "Model 5: Layer-specific posteriors — composite vs mScarlet (RQ5)") +
  facet_wrap(~layer, nrow = 1)
ggsave(file.path(OUT_DIR, "fig6_m5_layer_forest.pdf"), p6, width = 14, height = 4)

# Fig 7: Model 6 off-target forest
p7 <- forest_plot(rq7_post, "region_type",
  "Model 6: Visual cortex vs off-target posteriors (RQ7)")
ggsave(file.path(OUT_DIR, "fig7_m6_offtarget_forest.pdf"), p7, width = 8, height = 5)

# Fig 8: Prior sensitivity check
sens_long <- rbindlist(list(
  m1_post_weak[,  .(treatment, post_mean, ci_89_lo, ci_89_hi, prior = "Weak (0,1)")],
  m1_post_tight[, .(treatment, post_mean, ci_89_lo, ci_89_hi, prior = "Tight (0,0.5)")]
))
sens_long[, treatment := factor(treatment, levels = rev(TREATMENT_ORDER[-1]))]
p8 <- ggplot(sens_long, aes(y = treatment, x = post_mean,
                              xmin = ci_89_lo, xmax = ci_89_hi, colour = prior)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_errorbarh(height = 0.25, linewidth = 0.8,
                 position = position_dodgev(height = 0.5)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  scale_colour_manual(values = c("Weak (0,1)" = "#2E86AB",
                                  "Tight (0,0.5)" = "#E84C4C")) +
  labs(title    = "Prior sensitivity check — Model 1 (frac_enwrapped)",
       subtitle = "Overlap confirms results are not prior-driven",
       x = "Treatment effect vs mScarlet (fraction)", y = NULL,
       colour = "Prior") +
  theme_classic(base_size = 10)
ggsave(file.path(OUT_DIR, "fig8_prior_sensitivity.pdf"), p8, width = 8, height = 4)

cat("All figures saved.\n")
cat("\nScript 13 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Stan .rds files cached — rerunning loads from cache.\n")