# =============================================================================
# 12_pnn_integrity_composite.R
# PNN Integrity Composite Score — construction, validation, and LMM analysis
#
# The composite PNN integrity score (PNN-IS) combines two metrics available
# at the section level:
#
#   frac_enwrapped   — fraction of PV cells with adjacent PNN (r=0.75 with density)
#   normalized_btm20 — cell-level WFA fluorescence intensity  (r=0.56 with enwrap)
#
# avgPxIntensity and density are incorporated in the validation table (Part G)
# but not in the primary LMM because they are only available at region/animal
# level — using them in the LMM would collapse the data to one observation per
# animal, making the random intercept unidentifiable.
#
# diffFluo is excluded: near-zero correlation with frac_enwrapped (r=-0.08)
# and negatively correlated with normalized_btm20 (r=-0.62).
#
# Construction:
#   Each metric is standardised to mean=0, SD=1 across all sections.
#   The composite is the mean of the two standardised scores.
#   Higher composite = more intact PNNs; lower = more degraded.
#   Each animal contributes ~25 sections, ensuring the random intercept
#   for mouse_id is identifiable in all models.
#
# Research questions addressed:
#   RQ1 — Does any construct reduce overall PNN integrity in V1?
#   RQ2 — Is any effect spatially graded across injection zones?
#   RQ3 — Does C6ST-1 augment the effect of ADAMTS-15? (planned contrast)
#   RQ4 — Does monocular deprivation augment ADAMTS-4? (planned contrast)
#   RQ6 — Is any effect lateralised to the injected hemisphere?
#   RQ7 — Are effects confined to V1 or detectable in surrounding regions?
#
# Model structure:
#   Part A — Primary: composite ~ treatment + (1|mouse_id), V1 ipsilateral
#   Part B — Zone-stratified models (RQ2)
#   Part C — Zone x treatment interaction (RQ2)
#   Part D — Hemisphere comparison using lm() (RQ6)
#             (one value per animal per hemisphere after aggregation)
#   Part E — Planned contrasts: RQ3 (synergy) and RQ4 (MD)
#   Part F — Off-target specificity (RQ7)
#   Part G — Validation: all four metrics at animal level
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
  library(patchwork)
  library(scales)
})

# Paths
CELLS_CSV        <- "/path/to/analysis_results/cells_with_zones.csv"
MERGED_ZONES_CSV <- "/path/to/merged_datasets/merged_dataset_zones.csv"
MERGED_MAIN_CSV  <- "/path/to/merged_datasets/merged_dataset.csv"
RESULTS_DIR      <- "/path/to/results"
OUT_DIR          <- file.path(RESULTS_DIR, "12_pnn_integrity_composite")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Constants
TREATMENTS      <- c("mScarlet","ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1","C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS
EXCLUDE_PRIMARY <- "C6ST1_ADAMTS15_4"
ZONE_ORDER      <- c("Core","Penumbra","Outside")
COLOC_THRESH    <- 30L
MIN_PV          <- 5L
MIN_CELLS_ZONE  <- 5L
VISUAL_KEYWORDS <- c("visual","VIS")
ALPHA           <- 0.05

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

cat("Script 12: PNN Integrity Composite Score\n")
cat("Section-level composite: frac_enwrapped + normalized_btm20\n\n")

# Helpers

is_visual <- function(x) str_detect(tolower(x),
  paste(tolower(VISUAL_KEYWORDS), collapse = "|"))

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

extract_lm_contrasts <- function(fit, model_label = "") {
  if (is.null(fit)) return(NULL)
  coefs <- as.data.frame(summary(fit)$coefficients)
  coefs$term <- rownames(coefs)
  setDT(coefs)
  res <- coefs[term != "(Intercept)"]
  old <- intersect(names(res), c("Estimate","Std. Error","t value","Pr(>|t|)"))
  setnames(res, old, c("beta","se","t","p")[seq_along(old)])
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

emm_options(pbkrtest.limit = 20000, lmerTest.limit = 20000)

# =============================================================================
# CONSTRUCT SECTION-LEVEL COMPOSITE
# =============================================================================
cat("Loading cells_with_zones...\n")
cells <- fread(CELLS_CSV)
cells <- fix_names(cells)
cells <- cells[mouse_id != EXCLUDE_PRIMARY]
cells[, is_visual := is_visual(brain_area)]

# frac_enwrapped per section x zone
cat("Computing frac_enwrapped per section x zone...\n")
enwrap_records <- list()
for (mid in unique(cells$mouse_id)) {
  treat <- cells[mouse_id == mid, treatment[1]]
  for (sl in unique(cells[mouse_id == mid, slice_id])) {
    for (z in ZONE_ORDER) {
      grp <- cells[mouse_id == mid & slice_id == sl & zone == z & is_visual == TRUE]
      pv  <- grp[cell_type == "PV",  .(x_hires, y_hires)]
      pnn <- grp[cell_type == "PNN", .(x_hires, y_hires)]
      if (nrow(pv) < MIN_PV || nrow(pnn) < MIN_PV) next
      nn   <- nn2(data = as.matrix(pnn), query = as.matrix(pv), k = 1L)
      frac <- mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
      enwrap_records[[length(enwrap_records) + 1L]] <- data.table(
        mouse_id = mid, treatment = treat, slice_id = sl,
        zone = z, frac_enwrapped = frac, n_pv = nrow(pv)
      )
    }
  }
}
enwrap_dt <- rbindlist(enwrap_records)
cat(sprintf("  frac_enwrapped: %d section x zone records\n", nrow(enwrap_dt)))

# normalized_btm20 per section x zone
cat("Computing normalized_btm20 per section x zone...\n")
intensity_dt <- cells[cell_type == "PNN" & is_visual == TRUE &
                        !is.na(normalized_btm20) & normalized_btm20 > 0,
  .(normalized_btm20 = mean(normalized_btm20, na.rm = TRUE),
    n_cells = .N),
  by = .(mouse_id, treatment, slice_id, zone)]
intensity_dt <- intensity_dt[n_cells >= MIN_CELLS_ZONE]
cat(sprintf("  normalized_btm20: %d section x zone records\n", nrow(intensity_dt)))

# Merge and standardise
section_dt <- merge(enwrap_dt, intensity_dt,
                    by = c("mouse_id","treatment","slice_id","zone"),
                    all = FALSE)

mu_fe  <- mean(section_dt$frac_enwrapped,   na.rm = TRUE)
sd_fe  <- sd(section_dt$frac_enwrapped,     na.rm = TRUE)
mu_btm <- mean(section_dt$normalized_btm20, na.rm = TRUE)
sd_btm <- sd(section_dt$normalized_btm20,   na.rm = TRUE)

section_dt[, frac_enwrapped_z   := (frac_enwrapped   - mu_fe)  / sd_fe]
section_dt[, normalized_btm20_z := (normalized_btm20 - mu_btm) / sd_btm]
section_dt[, composite := rowMeans(.SD, na.rm = TRUE),
            .SDcols = c("frac_enwrapped_z","normalized_btm20_z")]
section_dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
section_dt[, zone      := factor(zone,      levels = ZONE_ORDER)]

cat(sprintf("\nSection-level composite: %d records, %d animals, ~%d sections/animal\n",
    nrow(section_dt), uniqueN(section_dt$mouse_id),
    round(nrow(section_dt) / uniqueN(section_dt$mouse_id))))

fwrite(section_dt, file.path(OUT_DIR, "composite_section_data.csv"))

# =============================================================================
# PART A: Primary model (RQ1)
# =============================================================================
cat("\nPART A: Primary composite model (RQ1)\n\n")

fit_primary <- fit_lmm(composite ~ treatment + (1 | mouse_id),
                        section_dt, "Primary_composite")
primary_res <- extract_contrasts(fit_primary, "primary_composite")
if (!is.null(primary_res)) {
  primary_res <- add_fdr(primary_res, suffix = "_adj")
  fwrite(primary_res, file.path(OUT_DIR, "primary_lmm_results.csv"))
  cat("Part A results:\n")
  print(primary_res[, .(term, beta, se, df, t, p, p_adj, sig_adj)])
}

# =============================================================================
# PART B: Zone-stratified models (RQ2)
# =============================================================================
cat("\nPART B: Zone-stratified composite models (RQ2)\n\n")

zone_res <- rbindlist(lapply(ZONE_ORDER, function(z) {
  d <- section_dt[zone == z]
  cat(sprintf("Zone: %s (%d sections, %d animals)\n",
      z, nrow(d), uniqueN(d$mouse_id)))
  fit <- fit_lmm(composite ~ treatment + (1 | mouse_id), d, sprintf("Zone_%s", z))
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, sprintf("zone_%s", z))
  res[, zone := z]
  res
}), fill = TRUE)

if (!is.null(zone_res) && nrow(zone_res) > 0) {
  zone_res <- add_fdr(zone_res, suffix = "_adj", by_col = "zone")
  fwrite(zone_res, file.path(OUT_DIR, "zone_lmm_results.csv"))
  sig_b <- zone_res[sig_adj != "ns"]
  cat("Part B significant/trending results:\n")
  if (nrow(sig_b) > 0) print(sig_b[, .(zone, term, beta, se, p, p_adj, sig_adj)])
  else cat("  No significant hits after FDR correction\n")
}

# =============================================================================
# PART C: Zone x treatment interaction (RQ2)
# =============================================================================
cat("\nPART C: Zone x treatment interaction (RQ2)\n\n")

fit_interaction <- fit_lmm(composite ~ treatment * zone + (1 | mouse_id),
                            section_dt, "Interaction_composite")
if (!is.null(fit_interaction)) {
  aov_res <- as.data.frame(anova(fit_interaction))
  aov_res$term <- rownames(aov_res)
  setDT(aov_res)
  fwrite(aov_res, file.path(OUT_DIR, "zone_interaction_anova.csv"))
  cat("Part C interaction ANOVA:\n")
  print(aov_res[, .(term, `F value`, `Pr(>F)`)])
}

# =============================================================================
# PART D: Hemisphere comparison using lm() (RQ6)
# =============================================================================
cat("\nPART D: Hemisphere comparison (RQ6)\n")
cat("lm() used: one composite per animal per hemisphere after aggregation\n\n")

# Compute composite for both hemispheres from cells
hemi_records <- list()
for (mid in unique(cells$mouse_id)) {
  treat <- cells[mouse_id == mid, treatment[1]]
  for (hemi in c("left","right")) {
    for (sl in unique(cells[mouse_id == mid & hemisphere == hemi, slice_id])) {
      for (z in ZONE_ORDER) {
        grp <- cells[mouse_id == mid & slice_id == sl &
                       hemisphere == hemi & zone == z & is_visual == TRUE]
        pv  <- grp[cell_type == "PV",  .(x_hires, y_hires)]
        pnn <- grp[cell_type == "PNN", .(x_hires, y_hires)]
        if (nrow(pv) < MIN_PV || nrow(pnn) < MIN_PV) next
        nn   <- nn2(data = as.matrix(pnn), query = as.matrix(pv), k = 1L)
        frac <- mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
        intens_vals <- grp[cell_type == "PNN" & !is.na(normalized_btm20) &
                             normalized_btm20 > 0, normalized_btm20]
        if (length(intens_vals) < MIN_CELLS_ZONE) next
        intens <- mean(intens_vals)
        hemi_records[[length(hemi_records) + 1L]] <- data.table(
          mouse_id = mid, treatment = treat, slice_id = sl,
          hemisphere = hemi, zone = z,
          frac_enwrapped = frac, normalized_btm20 = intens
        )
      }
    }
  }
}
hemi_section_dt <- rbindlist(hemi_records)
hemi_section_dt[, frac_enwrapped_z   := (frac_enwrapped   - mu_fe)  / sd_fe]
hemi_section_dt[, normalized_btm20_z := (normalized_btm20 - mu_btm) / sd_btm]
hemi_section_dt[, composite := rowMeans(.SD, na.rm = TRUE),
                 .SDcols = c("frac_enwrapped_z","normalized_btm20_z")]

# Aggregate to animal x hemisphere
hemi_animal <- hemi_section_dt[,
  .(composite = mean(composite, na.rm = TRUE)),
  by = .(mouse_id, treatment, hemisphere)]
hemi_animal[, treatment  := factor(treatment, levels = TREATMENT_ORDER)]
hemi_animal[, hemisphere := factor(hemisphere,
                                    levels = c("right","left"),
                                    labels = c("contralateral","ipsilateral"))]
fwrite(hemi_animal, file.path(OUT_DIR, "composite_hemisphere_data.csv"))

# lm() for ipsilateral and contralateral separately
fit_ipsi   <- lm(composite ~ treatment,
                  data = hemi_animal[hemisphere == "ipsilateral"])
fit_contra <- lm(composite ~ treatment,
                  data = hemi_animal[hemisphere == "contralateral"])

ipsi_res   <- extract_lm_contrasts(fit_ipsi,   "ipsi_composite")
contra_res <- extract_lm_contrasts(fit_contra, "contra_composite")

if (!is.null(ipsi_res)) {
  ipsi_res[, hemisphere := "ipsilateral"]
  ipsi_res <- add_fdr(ipsi_res, suffix = "_adj")
  fwrite(ipsi_res, file.path(OUT_DIR, "hemisphere_ipsi_contrasts.csv"))
  cat("Part D ipsilateral contrasts:\n")
  print(ipsi_res[, .(term, beta, se, p, p_adj, sig_adj)])
}
if (!is.null(contra_res)) {
  contra_res[, hemisphere := "contralateral"]
  contra_res <- add_fdr(contra_res, suffix = "_adj")
  fwrite(contra_res, file.path(OUT_DIR, "hemisphere_contra_contrasts.csv"))
}

# =============================================================================
# PART E: Planned contrasts RQ3 and RQ4
# =============================================================================
cat("\nPART E: Planned contrasts (RQ3 synergy, RQ4 MD)\n\n")

planned_res <- rbindlist(lapply(ZONE_ORDER, function(z) {
  d <- section_dt[zone == z]
  cat(sprintf("Zone: %s (%d sections)\n", z, nrow(d)))
  fit <- fit_lmm(composite ~ treatment + (1 | mouse_id), d,
                 sprintf("Planned_%s", z))
  if (is.null(fit)) return(NULL)
  emm <- emmeans(fit, ~ treatment)
  synergy <- as.data.table(contrast(emm,
    list("C6ST1_ADAMTS15 - ADAMTS15" = c(
      mScarlet=0, ADAMTS4=0, ADAMTS4_MD=0, ADAMTS15=-1,
      C6ST1=0, C6ST1_ADAMTS15=1))))
  synergy[, `:=`(contrast_label = "C6ST1_ADAMTS15 vs ADAMTS15",
                 rq = "RQ3", zone = z)]
  md <- as.data.table(contrast(emm,
    list("ADAMTS4_MD - ADAMTS4" = c(
      mScarlet=0, ADAMTS4=-1, ADAMTS4_MD=1, ADAMTS15=0,
      C6ST1=0, C6ST1_ADAMTS15=0))))
  md[, `:=`(contrast_label = "ADAMTS4_MD vs ADAMTS4",
            rq = "RQ4", zone = z)]
  rbind(synergy, md, fill = TRUE)
}), fill = TRUE)

if (!is.null(planned_res) && nrow(planned_res) > 0) {
  if ("estimate" %in% names(planned_res)) setnames(planned_res, "estimate", "beta")
  if ("p.value"  %in% names(planned_res)) setnames(planned_res, "p.value",  "p")
  planned_res <- add_fdr(planned_res, suffix = "_adj", by_col = "rq")
  fwrite(planned_res, file.path(OUT_DIR, "planned_contrasts.csv"))
  cat("Part E results:\n")
  print(planned_res[, .(rq, zone, contrast_label, beta, p, p_adj, sig_adj)])
}

# =============================================================================
# PART F: Off-target specificity (RQ7)
# =============================================================================
cat("\nPART F: Off-target specificity (RQ7)\n\n")

cells_off <- cells[is_visual == FALSE]
ot_records <- list()
for (mid in unique(cells_off$mouse_id)) {
  treat <- cells_off[mouse_id == mid, treatment[1]]
  for (sl in unique(cells_off[mouse_id == mid, slice_id])) {
    grp <- cells_off[mouse_id == mid & slice_id == sl]
    pv  <- grp[cell_type == "PV",  .(x_hires, y_hires)]
    pnn <- grp[cell_type == "PNN", .(x_hires, y_hires)]
    if (nrow(pv) < MIN_PV || nrow(pnn) < MIN_PV) next
    nn   <- nn2(data = as.matrix(pnn), query = as.matrix(pv), k = 1L)
    frac <- mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
    intens_vals <- grp[cell_type == "PNN" & !is.na(normalized_btm20) &
                         normalized_btm20 > 0, normalized_btm20]
    if (length(intens_vals) < MIN_CELLS_ZONE) next
    ot_records[[length(ot_records) + 1L]] <- data.table(
      mouse_id = mid, treatment = treat, slice_id = sl,
      frac_enwrapped = frac, normalized_btm20 = mean(intens_vals)
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

vc_dt <- section_dt[, .(composite = mean(composite, na.rm = TRUE)),
                     by = .(mouse_id, treatment, slice_id)]
vc_dt[, region_type := "Visual (target)"]

combined_rq7 <- rbindlist(
  list(vc_dt, ot_dt[, .(mouse_id, treatment, slice_id, composite, region_type)]),
  fill = TRUE)
combined_rq7[, region_type := factor(region_type,
                                      levels = c("Visual (target)","Off-target"))]
combined_rq7[, treatment   := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("Off-target sections: %d | Visual sections: %d\n",
    nrow(ot_dt), nrow(vc_dt)))

fit_rq7 <- fit_lmm(composite ~ treatment * region_type + (1 | mouse_id),
                    combined_rq7, "RQ7_interaction")
if (!is.null(fit_rq7)) {
  rq7_anova <- as.data.frame(anova(fit_rq7))
  rq7_anova$term <- rownames(rq7_anova)
  setDT(rq7_anova)
  fwrite(rq7_anova, file.path(OUT_DIR, "offtarget_interaction_anova.csv"))
  cat("Part F interaction ANOVA:\n")
  print(rq7_anova[, .(term, `F value`, `Pr(>F)`)])
}

rq7_contrasts <- rbindlist(lapply(c("Visual (target)","Off-target"), function(reg) {
  d <- combined_rq7[region_type == reg]
  fit <- fit_lmm(composite ~ treatment + (1 | mouse_id), d,
                 sprintf("RQ7_%s", str_replace_all(reg, "[ ()]", "_")))
  if (is.null(fit)) return(NULL)
  res <- extract_contrasts(fit, sprintf("rq7_%s", reg))
  res[, region := reg]
  res
}), fill = TRUE)

if (nrow(rq7_contrasts) > 0) {
  rq7_contrasts <- add_fdr(rq7_contrasts, suffix = "_adj", by_col = "region")
  fwrite(rq7_contrasts, file.path(OUT_DIR, "offtarget_contrasts.csv"))
  cat("Part F visual cortex contrasts (significant only):\n")
  vis <- rq7_contrasts[region == "Visual (target)" & sig_adj != "ns"]
  if (nrow(vis) > 0) print(vis[, .(region, term, beta, p, p_adj, sig_adj)])
  else cat("  No significant visual cortex hits\n")
  cat("Part F off-target contrasts (significant only):\n")
  off <- rq7_contrasts[region == "Off-target" & sig_adj != "ns"]
  if (nrow(off) > 0) print(off[, .(region, term, beta, p, p_adj, sig_adj)])
  else cat("  No significant off-target hits - spatial confinement confirmed\n")
}

# =============================================================================
# PART G: Validation - all four metrics at animal level
# =============================================================================
cat("\nPART G: Validation - all four metrics at animal level\n\n")

mz <- fread(MERGED_ZONES_CSV)
setnames(mz, "animal_id", "mouse_id")
mz <- fix_names(mz)
mz <- mz[mouse_id != EXCLUDE_PRIMARY]
mz[, is_visual := is_visual(brain_area)]
avgpx_animal <- mz[staining == "WFA" & is_visual == TRUE & !is.na(avgPxIntensity),
  .(avgPxIntensity = mean(avgPxIntensity, na.rm = TRUE)),
  by = .(mouse_id, treatment)]

mm <- fread(MERGED_MAIN_CSV)
setnames(mm, "animal_id", "mouse_id")
mm <- fix_names(mm)
mm <- mm[mouse_id != EXCLUDE_PRIMARY & zone == "Total"]
mm[, is_visual := is_visual(brain_area)]
density_animal <- mm[staining == "WFA" & is_visual == TRUE & !is.na(density),
  .(density = mean(density, na.rm = TRUE)),
  by = .(mouse_id, treatment)]

composite_animal <- section_dt[,
  .(frac_enwrapped   = mean(frac_enwrapped, na.rm = TRUE),
    normalized_btm20 = mean(normalized_btm20, na.rm = TRUE),
    composite        = mean(composite, na.rm = TRUE)),
  by = .(mouse_id, treatment)]

validation_dt <- merge(composite_animal, avgpx_animal,
                        by = c("mouse_id","treatment"), all.x = TRUE)
validation_dt <- merge(validation_dt, density_animal,
                        by = c("mouse_id","treatment"), all.x = TRUE)
validation_dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
fwrite(validation_dt, file.path(OUT_DIR, "composite_validation_animal.csv"))

cat("Animal-level validation table:\n")
print(validation_dt[,
  .(frac_enwrapped   = round(mean(frac_enwrapped,   na.rm=TRUE), 3),
    normalized_btm20 = round(mean(normalized_btm20, na.rm=TRUE), 3),
    avgPxIntensity   = round(mean(avgPxIntensity,   na.rm=TRUE), 3),
    density          = round(mean(density,           na.rm=TRUE), 3),
    composite        = round(mean(composite,         na.rm=TRUE), 3)),
  by = treatment][order(match(treatment, TREATMENT_ORDER))])

# =============================================================================
# FIGURES
# =============================================================================
cat("\nGenerating figures...\n")

# Figure 1: Primary composite bar chart with animal dots
anim_means <- composite_animal[,
  .(mean_comp = mean(composite, na.rm=TRUE),
    se_comp   = sd(composite,   na.rm=TRUE) / sqrt(.N)),
  by = treatment]
anim_means[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

p1 <- ggplot(anim_means, aes(x = treatment, y = mean_comp, fill = treatment)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_comp - se_comp, ymax = mean_comp + se_comp),
                width = 0.25, linewidth = 0.8) +
  geom_point(data = composite_animal,
             aes(x = treatment, y = composite),
             shape = 21, fill = "white", size = 2.5, stroke = 0.8,
             position = position_jitter(width = 0.08, seed = 42)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = PALETTE) +
  labs(title = "PNN Integrity Composite Score",
       subtitle = "Mean +/- SEM; dots = individual animals",
       x = NULL, y = "Composite score (SD units)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

# Figure 2: Zone gradient
zone_anim <- section_dt[,
  .(mean_comp = mean(composite, na.rm=TRUE),
    se_comp   = sd(composite,   na.rm=TRUE) / sqrt(.N)),
  by = .(treatment, zone)]
zone_anim[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
zone_anim[, zone      := factor(zone,      levels = ZONE_ORDER)]

p2 <- ggplot(zone_anim, aes(x = zone, y = mean_comp,
                              colour = treatment, group = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_comp - se_comp, ymax = mean_comp + se_comp),
                width = 0.15) +
  scale_colour_manual(values = PALETTE) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  labs(title = "PNN Integrity Composite - injection zone gradient",
       x = "Zone", y = "Composite score (SD units)",
       colour = "Treatment") +
  theme_classic(base_size = 11)

# Figure 3: All four metrics as % change from mScarlet
val_long <- melt(validation_dt,
                  id.vars = c("mouse_id","treatment"),
                  measure.vars = c("frac_enwrapped","normalized_btm20",
                                   "avgPxIntensity","density","composite"),
                  variable.name = "metric", value.name = "value")
val_long[, metric := factor(metric,
  levels = c("frac_enwrapped","normalized_btm20","avgPxIntensity","density","composite"),
  labels = c("Enwrapment","Cell\nintensity","Region\nintensity","Density","Composite"))]
val_long[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

ms_means <- val_long[treatment == "mScarlet",
  .(ms_mean = mean(value, na.rm=TRUE)), by = metric]
val_long <- merge(val_long, ms_means, by = "metric")
val_long[, pct_change := (value / ms_mean - 1) * 100]

profile_means <- val_long[,
  .(mean_pct = mean(pct_change, na.rm=TRUE),
    se_pct   = sd(pct_change,   na.rm=TRUE) / sqrt(.N)),
  by = .(treatment, metric)]

p3 <- ggplot(profile_means[treatment != "mScarlet"],
             aes(x = treatment, y = mean_pct, fill = treatment)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
                width = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = PALETTE) +
  labs(title = "Individual metric contributions (% change from mScarlet)",
       x = NULL, y = "% change") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

combined_fig <- p1 / p2 / p3 + plot_layout(heights = c(1.2, 1, 1))
ggsave(file.path(OUT_DIR, "fig_composite_validation.pdf"),
       combined_fig, width = 12, height = 13)
cat("Saved fig_composite_validation.pdf\n")

cat("\nScript 12 complete\n")
cat(sprintf("All outputs saved to: %s\n", OUT_DIR))