# =============================================================================
# Script 28: Continuous viral load and expression heterogeneity analyses
# =============================================================================
# Two analyses motivated by the 50:50 dual-construct delivery in the
# C6ST1_ADAMTS15 group, where each transduced cell expresses both enzymes
# but at stochastically varying ratios.
#
#   Analysis 1 — Continuous local viral load score (RQ3 + RQ1)
#     Replaces the binary 30px virus+/virus- flag from script 24 with a
#     continuous predictor: the distance-weighted sum of mScarlet intensities
#     of all C3 cells within a search radius (LOAD_RADIUS px) of each PV cell.
#     This captures how much enzymatic activity a PV cell is likely exposed to
#     from its immediate neighbourhood, proportional to both the number of
#     nearby transduced cells and their expression level.
#
#     GLMM: enwrapped ~ log(local_viral_load + 1) * treatment
#                       + (1|mouse_id/slice_id)
#     Restricted to Core zone, ipsilateral, visual cortex.
#
#     Key test: is the dose-response gradient (more local viral load → lower
#     enwrapment) steeper in C6ST1_ADAMTS15 than in single-enzyme groups?
#     Under the paracrine model, C6ST1_ADAMTS15 cells benefit from both
#     sulfation sensitisation and protease secretion simultaneously, so local
#     expression level should predict enwrapment reduction more strongly in
#     the combination group than in either enzyme alone.
#
#     Also fits per-treatment logistic regressions and reports slope comparison.
#
#   Analysis 2 — Expression heterogeneity as section-level predictor (RQ3)
#     Under the paracrine model, sections with high variance in mScarlet
#     intensity (many cells at C6ST1-dominant and ADAMTS15-dominant extremes)
#     should show more enwrapment reduction than sections with low variance
#     (uniform mid-level expression). High variance means the chance that a
#     C6ST1-expressing cell neighbours an ADAMTS15-expressing cell is higher.
#
#     For C6ST1_ADAMTS15 animals: correlate per-section mScarlet intensity CV
#     (coefficient of variation) with Core-zone enwrapment fraction.
#     The prediction is negative: more heterogeneous expression → more depletion.
#
#     Comparison: same correlation in single-enzyme groups should be absent or
#     weaker, because variance in a single enzyme does not confer paracrine
#     benefit in the same way.
#
# Inputs:
#   results/24_virus_cell_enwrapment/cell_virus_tags.csv  — PV coords + tags
#   analysis_results/cells_with_zones.csv                 — PNN coords
#   /media/.../Counts_C3-*.csv                            — C3 coords + intensity
#   results/00_enwrapment/slice_enwrapment.csv             — section enwrapment
#
# Outputs: results/28_viral_load_heterogeneity/
#
# NOTE ON INTERPRETATION:
#   Both constructs use mScarlet as the reporter via P2A. The mScarlet signal
#   therefore reflects total viral load, not the ratio of C6ST1 to ADAMTS15
#   expression in any given cell. Local viral load score conflates both enzyme
#   types. All cell-level inferences in C6ST1_ADAMTS15 animals are about total
#   local expression rather than either enzyme specifically. This is stated as
#   a design limitation in the methods.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(RANN)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(patchwork)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV   <- "/path/to/analysis_results/cells_with_zones.csv"
C3_DIR      <- "/path/to/c3_results"
TAGS_CSV    <- "/path/to/results/24_virus_cell_enwrapment/cell_virus_tags.csv"
SLICE_ENW   <- "/path/to/results/00_enwrapment/slice_enwrapment.csv"
RESULTS_DIR <- "/path/to/results"
OUT_DIR     <- file.path(RESULTS_DIR, "28_viral_load_heterogeneity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
COLOC_THRESH    <- 30L
LOAD_RADIUS     <- 150L     # px radius for local viral load score
EXCLUDE         <- ""       #Exclude animal (e.g., injection failure)
TREATMENT_ORDER <- c("mScarlet","ADAMTS4","ADAMTS4_MD",
                     "ADAMTS15","C6ST1","C6ST1_ADAMTS15")
VISUAL_RE       <- regex("isual", ignore_case = TRUE)
MIN_PV_SECTION  <- 5L
MIN_PV_GLMM     <- 20L     # minimum cells per treatment for GLMM slope

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)


cat("Script 28: Continuous viral load and expression heterogeneity\n")
cat("============================================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# ── PV cell coordinates with virus tags ──────────────────────────────────────
cat("Loading cell_virus_tags.csv ...\n")
tags <- fread(TAGS_CSV)
tags <- fix_names(tags)
tags <- tags[mouse_id != EXCLUDE]
tags[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

# ── PNN coordinates for enwrapment computation ────────────────────────────────
cat("Loading PNN coordinates ...\n")
pnn_cols <- c("mouse_id","slice_id","cell_type","hemisphere",
              "x_hires","y_hires","brain_area","zone","treatment")
cwz  <- fread(CELLS_CSV, select = pnn_cols)
cwz  <- fix_names(cwz)
cwz  <- cwz[mouse_id != EXCLUDE & str_detect(brain_area, VISUAL_RE)]
cwz  <- cwz[hemisphere == "left"]
pnn  <- cwz[cell_type == "PNN"]

# ── Recompute enwrapment ──────────────────────────────────────────────────────
cat("Computing enwrapment per PV cell ...\n")
tags[, enwrapped := FALSE]
ew_keys <- unique(tags[, .(mouse_id, slice_id)])
for (i in seq_len(nrow(ew_keys))) {
  mid <- ew_keys$mouse_id[i]; sid <- ew_keys$slice_id[i]
  pv_sl  <- tags[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  pnn_sl <- pnn[mouse_id  == mid & slice_id == sid, .(x_hires, y_hires)]
  if (nrow(pv_sl) == 0L || nrow(pnn_sl) == 0L) next
  nn  <- nn2(data = as.matrix(pnn_sl), query = as.matrix(pv_sl), k = 1L)
  idx <- which(tags$mouse_id == mid & tags$slice_id == sid)
  tags[idx, enwrapped := nn$nn.dists[, 1L] <= COLOC_THRESH]
}
tags[, enwrapped_int := as.integer(enwrapped)]

# ── C3 coordinates WITH intensity ─────────────────────────────────────────────
cat("Loading C3 coordinates and intensities ...\n")

parse_c3_name <- function(fname) {
  m <- regmatches(fname,
        regexec("Counts_C3-(.+)_(\\d+)_(s\\d+)\\.csv", fname))[[1]]
  if (length(m) != 4L) return(NULL)
  list(treatment_raw = m[2], animal_num = m[3], slice_id = m[4])
}

c3_files <- list.files(C3_DIR, pattern = "^Counts_C3-.*\\.csv$",
                       full.names = TRUE)
c3_list  <- vector("list", length(c3_files))

for (i in seq_along(c3_files)) {
  meta <- parse_c3_name(basename(c3_files[i]))
  if (is.null(meta)) next
  # Load coordinates AND intensity
  df <- tryCatch(
    fread(c3_files[i],
          select = c("Global_X","Global_Y","Mean_Intensity")),
    error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) next
  if (!all(c("Global_X","Global_Y","Mean_Intensity") %in% names(df))) next
  df[, `:=`(treatment_raw = meta$treatment_raw,
             animal_num    = meta$animal_num,
             slice_id      = meta$slice_id)]
  setnames(df, c("Global_X","Global_Y","Mean_Intensity"),
               c("x_c3","y_c3","intensity_c3"))
  c3_list[[i]] <- df
}

c3 <- rbindlist(c3_list, fill = TRUE)
c3[, mouse_id := paste0(
  str_replace(treatment_raw, "C6ST1_ADAMTS4$", "C6ST1_ADAMTS15"),
  "_", animal_num)]
c3[mouse_id == "mScarlet_4",   treatment_raw := "ADAMTS4_MD"]
c3[mouse_id == "ADAMTS4_MD_4", treatment_raw := "mScarlet"]
c3 <- c3[mouse_id != EXCLUDE]
cat(sprintf("  %s C3 cells with intensity loaded\n\n",
            format(nrow(c3), big.mark = ",")))

# =============================================================================
# ANALYSIS 1 — Continuous local viral load score
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Analysis 1: Continuous local viral load score\n")
cat(sprintf("  Search radius: %d px (~%.0f µm)\n",
            LOAD_RADIUS, LOAD_RADIUS * 0.325))
cat("  Score = sum(intensity_c3 / distance) for all C3 cells within radius\n")
cat("══════════════════════════════════════════════════════\n\n")

# Restrict to Core zone, ipsilateral, visual cortex
core_tags <- tags[zone == "Core"]

cat("Computing local viral load scores per PV cell ...\n")
core_tags[, local_viral_load := 0]

slice_keys_1 <- unique(core_tags[, .(mouse_id, slice_id)])

for (i in seq_len(nrow(slice_keys_1))) {
  mid <- slice_keys_1$mouse_id[i]
  sid <- slice_keys_1$slice_id[i]

  pv_sl <- core_tags[mouse_id == mid & slice_id == sid,
                      .(x_hires, y_hires)]
  c3_sl <- c3[mouse_id == mid & slice_id == sid,
               .(x_c3, y_c3, intensity_c3)]

  if (nrow(pv_sl) == 0L || nrow(c3_sl) == 0L) next

  # For each PV cell: distance-weighted intensity sum within LOAD_RADIUS
  # Use nn2 with k = min(nrow(c3_sl), 50) to get nearest C3 neighbours
  k_use <- min(nrow(c3_sl), 50L)
  nn    <- nn2(data  = as.matrix(c3_sl[, .(x_c3, y_c3)]),
               query = as.matrix(pv_sl),
               k     = k_use)

  idx <- which(core_tags$mouse_id == mid & core_tags$slice_id == sid)

  scores <- numeric(nrow(pv_sl))
  for (j in seq_len(nrow(pv_sl))) {
    dists    <- nn$nn.dists[j, ]
    in_range <- dists <= LOAD_RADIUS & dists > 0
    if (!any(in_range)) next
    nn_idx      <- nn$nn.idx[j, in_range]
    nn_dists    <- dists[in_range]
    nn_intens   <- c3_sl$intensity_c3[nn_idx]
    # Distance-weighted sum: closer cells contribute more
    scores[j] <- sum(nn_intens / nn_dists, na.rm = TRUE)
  }
  core_tags[idx, local_viral_load := scores]

  if (i %% 200L == 0L)
    cat(sprintf("  %d / %d slices done\n", i, nrow(slice_keys_1)))
}

core_tags[, log_viral_load := log(local_viral_load + 1)]
cat(sprintf("\n  %d%% of Core PV cells have non-zero local load\n",
            round(100 * mean(core_tags$local_viral_load > 0))))

fwrite(core_tags[, .(mouse_id, slice_id, treatment, zone,
                      x_hires, y_hires, enwrapped, enwrapped_int,
                      virus_positive, local_viral_load, log_viral_load)],
       file.path(OUT_DIR, "pv_cells_viral_load.csv"))

# ── Per-treatment GLMM: enwrapped ~ log_viral_load + (1|mouse_id/slice_id) ───
cat("\nFitting per-treatment logistic regressions (viral load → enwrapment) ...\n")

per_tx_slopes <- rbindlist(lapply(TREATMENT_ORDER, function(tx) {
  d <- core_tags[treatment == tx & local_viral_load > 0]
  if (nrow(d) < MIN_PV_GLMM) return(NULL)
  d[, slice_in_mouse := paste0(mouse_id, "_", slice_id)]
  fit <- tryCatch(
    glmer(enwrapped_int ~ log_viral_load + (1 | mouse_id/slice_id),
          data    = d,
          family  = binomial(link = "logit"),
          control = glmerControl(optimizer  = "bobyqa",
                                 optCtrl    = list(maxfun = 2e5))),
    error   = function(e) NULL,
    warning = function(w) {
      tryCatch(
        glmer(enwrapped_int ~ log_viral_load + (1 | mouse_id),
              data    = d,
              family  = binomial(link = "logit"),
              control = glmerControl(optimizer  = "bobyqa",
                                     optCtrl    = list(maxfun = 2e5))),
        error = function(e2) NULL)
    }
  )
  if (is.null(fit)) return(NULL)
  coef_fit <- coef(summary(fit))
  if (!"log_viral_load" %in% rownames(coef_fit)) return(NULL)
  data.table(
    treatment  = tx,
    slope      = coef_fit["log_viral_load","Estimate"],
    SE         = coef_fit["log_viral_load","Std. Error"],
    z          = coef_fit["log_viral_load","z value"],
    p          = coef_fit["log_viral_load","Pr(>|z|)"],
    n_cells    = nrow(d),
    OR_per_unit = exp(coef_fit["log_viral_load","Estimate"])
  )
}), fill = TRUE)

if (!is.null(per_tx_slopes) && nrow(per_tx_slopes) > 0) {
  per_tx_slopes[, sig := fcase(p < 0.001,"***", p < 0.01,"**",
                                p < 0.05,"*",   p < 0.10,".",
                                default = "ns")]
  per_tx_slopes[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
  cat("\nViral load → enwrapment slope per treatment (Core zone):\n")
  cat("  Negative slope = more local viral load → lower enwrapment probability\n\n")
  print(per_tx_slopes[order(treatment),
                       .(treatment, slope = round(slope,4),
                         SE = round(SE,4), z = round(z,3),
                         p = round(p,4), OR_per_unit = round(OR_per_unit,3),
                         sig, n_cells)])
  fwrite(per_tx_slopes, file.path(OUT_DIR, "analysis1_per_tx_slopes.csv"))
}

# ── Interaction GLMM: does slope differ between C6ST1_ADAMTS15 and ADAMTS15? ─
cat("\nGLMM: enwrapped ~ log_viral_load * treatment + (1|mouse_id)\n")
cat("  Restricted to ADAMTS15 and C6ST1_ADAMTS15\n\n")

rq3_load <- core_tags[treatment %in% c("ADAMTS15","C6ST1_ADAMTS15") &
                        local_viral_load > 0]
rq3_load[, treatment := factor(treatment,
                                levels = c("ADAMTS15","C6ST1_ADAMTS15"))]

fit_inter_1 <- tryCatch(
  glmer(enwrapped_int ~ log_viral_load * treatment + (1 | mouse_id),
        data    = rq3_load,
        family  = binomial(link = "logit"),
        control = glmerControl(optimizer  = "bobyqa",
                               optCtrl    = list(maxfun = 2e5))),
  error   = function(e) { cat("  GLMM failed:", e$message, "\n"); NULL },
  warning = function(w) {
    cat("  Warning:", w$message, "\n")
    tryCatch(
      glmer(enwrapped_int ~ log_viral_load * treatment + (1 | mouse_id),
            data    = rq3_load,
            family  = binomial(link = "logit"),
            control = glmerControl(optimizer  = "bobyqa",
                                   optCtrl    = list(maxfun = 2e5))),
      error = function(e2) NULL)
  }
)

if (!is.null(fit_inter_1)) {
  coef_1 <- as.data.table(coef(summary(fit_inter_1)), keep.rownames = "term")
  setnames(coef_1,
           intersect(c("Estimate","Std. Error","z value","Pr(>|z|)"),
                     names(coef_1)),
           c("estimate","SE","z","p")[
             seq_along(intersect(c("Estimate","Std. Error","z value","Pr(>|z|)"),
                                  names(coef_1)))])
  coef_1[, OR := exp(estimate)]
  coef_1[, sig := fcase(p < 0.001,"***", p < 0.01,"**",
                         p < 0.05,"*",   p < 0.10,".",
                         default = "ns")]
  cat("Interaction model: ADAMTS15 vs C6ST1_ADAMTS15 dose-response:\n")
  print(coef_1[, .(term, estimate = round(estimate,4),
                    SE = round(SE,4), z = round(z,3),
                    p = round(p,4), OR = round(OR,3), sig)])
  fwrite(coef_1, file.path(OUT_DIR, "analysis1_interaction_glmm.csv"))

  inter_row_1 <- coef_1[str_detect(term, ":")]
  if (nrow(inter_row_1) > 0) {
    cat(sprintf(
      "\nInteraction term (viral load × C6ST1_ADAMTS15): OR=%.3f, p=%.4f, %s\n",
      inter_row_1$OR, inter_row_1$p, inter_row_1$sig))
    cat(sprintf("  OR > 1: C6ST1_ADAMTS15 dose-response WEAKER than ADAMTS15\n"))
    cat(sprintf("  OR < 1: C6ST1_ADAMTS15 dose-response STEEPER than ADAMTS15\n"))
  }
}

# ── Figure 1: LOESS curves of enwrapment vs viral load per treatment ──────────
# Bin viral load for plotting clarity
core_tags[, load_bin := cut(log_viral_load,
                              breaks = quantile(log_viral_load[log_viral_load > 0],
                                                probs = seq(0, 1, 0.1),
                                                na.rm = TRUE),
                              include.lowest = TRUE)]

plot_load <- core_tags[local_viral_load > 0 & !is.na(load_bin), .(
  frac_enwrapped = mean(enwrapped),
  load_mid       = median(log_viral_load),
  n              = .N
), by = .(treatment, load_bin)]
plot_load <- plot_load[n >= 20]
plot_load[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

p1 <- ggplot(core_tags[local_viral_load > 0],
             aes(x = log_viral_load, y = enwrapped_int,
                 colour = treatment)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, alpha = 0.15, linewidth = 0.9) +
  scale_colour_manual(values = PALETTE) +
  facet_wrap(~ treatment, nrow = 2) +
  labs(
    title    = "Analysis 1: Enwrapment probability vs local viral load (Core zone)",
    subtitle = "Steeper negative slope = stronger dose-response",
    x        = "log(distance-weighted mScarlet intensity + 1)",
    y        = "P(enwrapped)",
    colour   = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig1_viral_load_dose_response.pdf"),
       p1, width = 12, height = 6)
cat("\n  Saved fig1_viral_load_dose_response.pdf\n\n")

# Forest plot of slopes
if (!is.null(per_tx_slopes) && nrow(per_tx_slopes) > 0) {
  per_tx_slopes[, tx := factor(treatment, levels = rev(TREATMENT_ORDER))]
  p1b <- ggplot(per_tx_slopes,
                aes(x = slope, y = tx, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = slope - 1.96*SE, xmax = slope + 1.96*SE),
                  width = 0.25, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(x = slope + sign(slope)*0.05,
                  label = sig),
              hjust = 0.5, size = 3.5, fontface = "bold", colour = "black") +
    scale_colour_manual(values = PALETTE, guide = "none") +
    labs(
      title    = "Viral load → enwrapment slope per treatment (logit scale)",
      subtitle = "Negative = more local expression → lower enwrapment",
      x        = "Slope (log viral load coefficient)",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(OUT_DIR, "fig1b_slope_forest.pdf"),
         p1b, width = 8, height = 4)
  cat("  Saved fig1b_slope_forest.pdf\n\n")
}

# =============================================================================
# ANALYSIS 2 — Expression heterogeneity as section-level predictor
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Analysis 2: Expression heterogeneity and enwrapment\n")
cat("  Hypothesis: sections with high mScarlet CV show more PNN depletion\n")
cat("  in C6ST1_ADAMTS15 (paracrine model) but not in single-enzyme groups\n")
cat("══════════════════════════════════════════════════════\n\n")

# Per-section C3 intensity statistics
cat("Computing per-section mScarlet intensity statistics ...\n")

c3_section <- c3[, .(
  n_c3          = .N,
  mean_intens   = mean(intensity_c3, na.rm = TRUE),
  sd_intens     = sd(intensity_c3,   na.rm = TRUE),
  median_intens = median(intensity_c3, na.rm = TRUE),
  total_intens  = sum(intensity_c3,  na.rm = TRUE)
), by = .(mouse_id, slice_id)]

# CV = sd / mean (only meaningful when mean > 0 and n >= 3)
c3_section[, cv_intens := fifelse(
  mean_intens > 0 & n_c3 >= 3L,
  sd_intens / mean_intens,
  NA_real_)]

# Section-level enwrapment from script 00
cat("Loading slice_enwrapment.csv ...\n")
slice_enw <- fread(SLICE_ENW)
slice_enw <- fix_names(slice_enw)
slice_enw <- slice_enw[mouse_id != EXCLUDE]
# Use ipsilateral enwrapment fraction
slice_enw_ipsi <- slice_enw[, .(mouse_id, slice_id, treatment,
                                  frac_ipsi, ratio)]

# Join heterogeneity stats to enwrapment
het_dt <- merge(c3_section, slice_enw_ipsi,
                by = c("mouse_id","slice_id"))
het_dt <- het_dt[!is.na(cv_intens) & !is.na(frac_ipsi)]
het_dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("  Sections with valid CV and enwrapment: %d\n\n",
            nrow(het_dt)))

fwrite(het_dt, file.path(OUT_DIR, "analysis2_section_heterogeneity.csv"))

# ── Per-treatment Pearson correlation: CV vs enwrapment ──────────────────────
cat("Pearson correlation: mScarlet CV vs ipsilateral enwrapment per treatment:\n")
cat("  (negative r = more heterogeneous expression → lower enwrapment)\n\n")

cor_results <- rbindlist(lapply(TREATMENT_ORDER, function(tx) {
  d <- het_dt[treatment == tx]
  if (nrow(d) < 5L) return(NULL)
  ct <- cor.test(d$cv_intens, d$frac_ipsi, method = "pearson")
  data.table(
    treatment = tx,
    r         = round(ct$estimate, 4),
    p         = round(ct$p.value, 4),
    n         = nrow(d),
    direction = ifelse(ct$estimate < 0,
                       "more heterogeneous → lower enwrapment",
                       "more heterogeneous → higher enwrapment")
  )
}), fill = TRUE)

print(cor_results[, .(treatment, r, p, n, direction)])
fwrite(cor_results, file.path(OUT_DIR, "analysis2_cv_correlations.csv"))

# ── Test whether C6ST1_ADAMTS15 correlation differs from single-enzyme groups ─
cat("\nFisher z-test: is C6ST1_ADAMTS15 correlation stronger than ADAMTS15?\n")

r_c6a15 <- cor_results[treatment == "C6ST1_ADAMTS15", r]
r_a15   <- cor_results[treatment == "ADAMTS15", r]
n_c6a15 <- cor_results[treatment == "C6ST1_ADAMTS15", n]
n_a15   <- cor_results[treatment == "ADAMTS15", n]

if (!is.null(r_c6a15) && !is.null(r_a15) &&
    !is.na(r_c6a15) && !is.na(r_a15)) {
  z1    <- 0.5 * log((1 + r_c6a15) / (1 - r_c6a15))
  z2    <- 0.5 * log((1 + r_a15)   / (1 - r_a15))
  se_diff <- sqrt(1/(n_c6a15 - 3) + 1/(n_a15 - 3))
  z_diff  <- (z1 - z2) / se_diff
  p_diff  <- 2 * pnorm(-abs(z_diff))
  cat(sprintf("  C6ST1_ADAMTS15: r=%.3f (n=%d)\n", r_c6a15, n_c6a15))
  cat(sprintf("  ADAMTS15:       r=%.3f (n=%d)\n", r_a15,   n_a15))
  cat(sprintf("  Fisher z-test: z=%.3f, p=%.4f\n", z_diff, p_diff))
  cat(sprintf("  %s\n\n",
              ifelse(r_c6a15 < r_a15 && p_diff < 0.10,
                     "C6ST1_ADAMTS15 shows stronger negative correlation — supports paracrine model",
                     "No significant difference in correlation strength")))

  fisher_summary <- data.table(
    comparison       = "C6ST1_ADAMTS15 vs ADAMTS15 (CV correlation)",
    r_c6a15          = round(r_c6a15, 4),
    r_a15            = round(r_a15,   4),
    z_stat           = round(z_diff,  4),
    p_value          = round(p_diff,  4),
    interpretation   = ifelse(r_c6a15 < r_a15 && p_diff < 0.10,
                              "stronger paracrine gradient in combination",
                              "no significant difference")
  )
  fwrite(fisher_summary, file.path(OUT_DIR, "analysis2_fisher_test.csv"))
}

# ── Figure 2: scatter plots CV vs enwrapment per treatment ────────────────────
p2 <- ggplot(het_dt[treatment %in% c("ADAMTS15","C6ST1_ADAMTS15",
                                       "C6ST1","mScarlet")],
             aes(x = cv_intens, y = frac_ipsi, colour = treatment)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text(
    data = cor_results[treatment %in% c("ADAMTS15","C6ST1_ADAMTS15",
                                         "C6ST1","mScarlet")],
    aes(x = Inf, y = Inf,
        label = sprintf("r=%.3f, p=%.3f", r, p),
        colour = treatment),
    hjust = 1.1, vjust = 1.3, size = 3.2, inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  facet_wrap(~ treatment, nrow = 2) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  labs(
    title    = "Analysis 2: mScarlet expression CV vs enwrapment per section",
    subtitle = paste("Negative slope in C6ST1_ADAMTS15 = more expression heterogeneity",
                     "→ more PNN depletion (paracrine prediction)"),
    x        = "CV of mScarlet intensity per section",
    y        = "Ipsilateral enwrapment fraction"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig2_cv_vs_enwrapment.pdf"),
       p2, width = 10, height = 7)
cat("  Saved fig2_cv_vs_enwrapment.pdf\n\n")

# ── Figure 3: CV distribution per treatment ───────────────────────────────────
p3 <- ggplot(het_dt,
             aes(x = treatment, y = cv_intens,
                 colour = treatment, fill = treatment)) +
  geom_violin(alpha = 0.3, linewidth = 0.6) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.3, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  scale_fill_manual(values   = PALETTE, guide = "none") +
  scale_x_discrete(limits = TREATMENT_ORDER) +
  labs(
    title    = "mScarlet intensity CV per section by treatment",
    subtitle = "Higher CV = more variable expression across transduced cells",
    x        = NULL,
    y        = "Coefficient of variation (SD / mean)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig3_cv_distribution.pdf"),
       p3, width = 8, height = 5)
cat("  Saved fig3_cv_distribution.pdf\n\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Script 28 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  pv_cells_viral_load.csv\n")
cat("  analysis1_per_tx_slopes.csv\n")
cat("  analysis1_interaction_glmm.csv\n")
cat("  analysis2_section_heterogeneity.csv\n")
cat("  analysis2_cv_correlations.csv\n")
cat("  analysis2_fisher_test.csv\n")
cat("  fig1_viral_load_dose_response.pdf\n")
cat("  fig1b_slope_forest.pdf\n")
cat("  fig2_cv_vs_enwrapment.pdf\n")
cat("  fig3_cv_distribution.pdf\n")
cat("══════════════════════════════════════════════════════\n")
