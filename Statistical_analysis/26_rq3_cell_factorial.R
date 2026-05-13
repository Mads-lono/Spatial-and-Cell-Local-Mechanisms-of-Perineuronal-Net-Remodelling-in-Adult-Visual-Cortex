# =============================================================================
# Script 26: Cell-level RQ3 factorial analysis — substrate conditioning
# =============================================================================
# Research question:
#   Does co-expression of C6ST-1 with ADAMTS-15 produce a greater reduction
#   in PNN enwrapment at the cell level than ADAMTS-15 alone, and does any
#   enhancement reflect a local proximity effect consistent with substrate
#   conditioning?
#
# Animal-level evidence for RQ3 is provided by scripts 20, 21, and 23
# (M7, M8). This script brings the same logic down to the cell level using
# the virus+ tagging from script 24. It provides three complementary tests:
#
#   Part A — Three-way factorial GLMM (primary cell-level RQ3 test)
#             enwrapped ~ ADAMTS15_present * C6ST1_present * virus_positive
#                         + (1|mouse_id/slice_id)
#             Groups: ADAMTS4=neither (reference), ADAMTS15, C6ST1, C6ST1_ADAMTS15
#             The three-way interaction is the key quantity:
#             P(combination produces a larger virus+/virus- gap than the sum
#               of individual gaps) = cell-level analog of M8.
#             Also reports: Bayesian version with brms for posterior inference.
#
#   Part B — Direct virus+ comparison: C6ST1_ADAMTS15 vs ADAMTS15 (RQ3)
#             Are PV cells immediately adjacent to transduced cells more
#             depleted in C6ST1_ADAMTS15 than in ADAMTS15?
#             Section-level LMM on virus+ cells only:
#             frac_enwrapped ~ treatment + (1|mouse_id)
#             Pre-specified contrast: C6ST1_ADAMTS15 vs ADAMTS15, no correction.
#             Also: full distance-bin profile for virus+ cells, two groups only.
#
#   Part C — Diffuse enhancement: C6ST1_ADAMTS15 vs ADAMTS15, virus- cells
#             Does C6ST-1 co-expression also enhance the diffuse extracellular
#             effect of ADAMTS15 on non-adjacent cells?
#             Section-level LMM on virus- cells only, same contrast.
#             If positive: C6ST1 sensitises the broader ECM, not only at
#             cell-proximal contacts.
#             If null: enhancement is purely local.
#
# Inputs:
#   results/24_virus_cell_enwrapment/cell_virus_tags.csv  — per-cell virus tags
#   analysis_results/cells_with_zones.csv                 — PNN coords for nn2
#   /media/.../Counts_C3-*.csv                            — C3 coords (Part B dist)
#
# Outputs: results/26_rq3_cell_factorial/
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
  library(brms)
  library(posterior)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV   <- "/path/to/analysis_results/cells_with_zones.csv"
C3_DIR      <- "/path/to/c3_results"
TAGS_CSV    <- "/path/to/results/24_virus_cell_enwrapment/cell_virus_tags.csv"
RESULTS_DIR <- "/path/to/results"
OUT_DIR     <- file.path(RESULTS_DIR, "26_rq3_cell_factorial")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
COLOC_THRESH    <- 30L
MIN_PV_SECTION  <- 5L
MIN_PV_BIN      <- 10L
PX_TO_MM        <- 0.325 / 1000
EXCLUDE         <- ""    # Exclude mouse_id if needed (e.g. injection failure)
TREATMENT_ORDER <- c("mScarlet","ADAMTS4","ADAMTS4_MD",
                     "ADAMTS15","C6ST1","C6ST1_ADAMTS15")

# Factorial groups only (mirrors M8 in script 23)
FACTORIAL_GROUPS <- c("ADAMTS4","ADAMTS15","C6ST1","C6ST1_ADAMTS15")
VISUAL_RE        <- regex("isual", ignore_case = TRUE)

DIST_BREAKS <- c(0, 30, 100, 300, 600, Inf)
DIST_LABELS <- c("0–30px","30–100px","100–300px","300–600px","600px+")

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)
PALETTE_RQ3 <- c(ADAMTS15 = "#4ef196", C6ST1_ADAMTS15 = "#f1c44e")


cat("Script 26: Cell-level RQ3 factorial analysis\n")
cat("============================================================\n\n")

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

# ── Load virus tags and recompute enwrapment (mirrors script 25) ──────────────
cat("Loading cell_virus_tags.csv ...\n")
tags <- fread(TAGS_CSV)
tags <- fix_names(tags)
tags <- tags[mouse_id != EXCLUDE]
tags[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
cat(sprintf("  %s PV cells loaded\n", format(nrow(tags), big.mark = ",")))

cat("Loading PNN coordinates from cells_with_zones.csv ...\n")
pnn_load_cols <- c("mouse_id","slice_id","cell_type","hemisphere",
                   "x_hires","y_hires","brain_area","treatment")
cwz <- fread(CELLS_CSV, select = pnn_load_cols)
cwz <- fix_names(cwz)
cwz <- cwz[mouse_id != EXCLUDE]
cwz <- cwz[str_detect(brain_area, VISUAL_RE)]
cwz <- cwz[hemisphere == "left"]
pnn_coords <- cwz[cell_type == "PNN"]
cat(sprintf("  %s PNN cells loaded\n\n", format(nrow(pnn_coords), big.mark = ",")))

cat("Computing enwrapment per PV cell via nn2 ...\n")
tags[, enwrapped := FALSE]
enwrap_keys <- unique(tags[, .(mouse_id, slice_id)])

for (i in seq_len(nrow(enwrap_keys))) {
  mid <- enwrap_keys$mouse_id[i]
  sid <- enwrap_keys$slice_id[i]
  pv_sl  <- tags[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  pnn_sl <- pnn_coords[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  if (nrow(pv_sl) == 0L || nrow(pnn_sl) == 0L) next
  nn   <- nn2(data = as.matrix(pnn_sl), query = as.matrix(pv_sl), k = 1L)
  flag <- nn$nn.dists[, 1L] <= COLOC_THRESH
  idx  <- which(tags$mouse_id == mid & tags$slice_id == sid)
  tags[idx, enwrapped := flag]
}

tags[, enwrapped_int := as.integer(enwrapped)]
cat(sprintf("  Overall enwrapment: %.1f%%\n\n", 100 * mean(tags$enwrapped)))

# ── Encode factorial indicators (mirrors M8 in script 23) ────────────────────
fact <- tags[treatment %in% FACTORIAL_GROUPS & zone == "Core"]
fact[, ADAMTS15_present := as.integer(treatment %in% c("ADAMTS15","C6ST1_ADAMTS15"))]
fact[, C6ST1_present    := as.integer(treatment %in% c("C6ST1","C6ST1_ADAMTS15"))]
fact[, virus_factor     := factor(virus_positive,
                                   levels = c(FALSE, TRUE),
                                   labels = c("virus-","virus+"))]

cat(sprintf("Factorial cell-level dataset (Core zone, 4 groups): %s cells\n\n",
            format(nrow(fact), big.mark = ",")))
print(fact[, .N, by = .(treatment, virus_factor)])

# ── Load C3 coords for distance analysis in Part B ────────────────────────────
cat("\nLoading C3 coordinates ...\n")

parse_c3_name <- function(fname) {
  m <- regmatches(fname,
        regexec("Counts_C3-(.+)_(\\d+)_(s\\d+)\\.csv", fname))[[1]]
  if (length(m) != 4L) return(NULL)
  list(treatment_raw = m[2], animal_num = m[3], slice_id = m[4])
}

c3_files <- list.files(C3_DIR, pattern = "^Counts_C3-.*\\.csv$", full.names = TRUE)
c3_list  <- vector("list", length(c3_files))

for (i in seq_along(c3_files)) {
  meta <- parse_c3_name(basename(c3_files[i]))
  if (is.null(meta)) next
  df <- tryCatch(fread(c3_files[i], select = c("Global_X","Global_Y")),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) next
  if (!all(c("Global_X","Global_Y") %in% names(df))) next
  df[, `:=`(treatment_raw = meta$treatment_raw,
             animal_num    = meta$animal_num,
             slice_id      = meta$slice_id)]
  setnames(df, c("Global_X","Global_Y"), c("x_c3","y_c3"))
  c3_list[[i]] <- df
}

c3 <- rbindlist(c3_list, fill = TRUE)
c3[, mouse_id := paste0(
  str_replace(treatment_raw, "C6ST1_ADAMTS4$", "C6ST1_ADAMTS15"),
  "_", animal_num)]
c3[mouse_id == "mScarlet_4",   treatment_raw := "ADAMTS4_MD"]
c3[mouse_id == "ADAMTS4_MD_4", treatment_raw := "mScarlet"]
c3 <- c3[mouse_id != EXCLUDE]
cat(sprintf("  %s C3 cells loaded\n\n", format(nrow(c3), big.mark = ",")))

# =============================================================================
# PART A — Three-way factorial GLMM + Bayesian version
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part A: Three-way factorial GLMM\n")
cat("══════════════════════════════════════════════════════\n\n")
cat("  enwrapped ~ ADAMTS15_present * C6ST1_present * virus_factor\n")
cat("             + (1|mouse_id/slice_id)\n\n")
cat("  Key: three-way interaction = substrate conditioning at cell level\n")
cat("  Positive interaction (on logit scale, virus+ direction) means\n")
cat("  the combination produces a LARGER virus+/virus- gap than expected\n")
cat("  from adding individual C6ST1 and ADAMTS15 gaps.\n\n")

# ── Frequentist GLMM ─────────────────────────────────────────────────────────
fit_3way <- tryCatch(
  glmer(enwrapped_int ~
          ADAMTS15_present * C6ST1_present * virus_factor +
          (1 | mouse_id/slice_id),
        data    = fact,
        family  = binomial(link = "logit"),
        control = glmerControl(optimizer  = "bobyqa",
                               optCtrl    = list(maxfun = 2e5))),
  error   = function(e) { cat("  Nested GLMM failed:", e$message, "\n"); NULL },
  warning = function(w) {
    cat("  Nested GLMM warning:", w$message,
        "\n  Retrying with (1|mouse_id) ...\n")
    tryCatch(
      glmer(enwrapped_int ~
              ADAMTS15_present * C6ST1_present * virus_factor +
              (1 | mouse_id),
            data    = fact,
            family  = binomial(link = "logit"),
            control = glmerControl(optimizer  = "bobyqa",
                                   optCtrl    = list(maxfun = 2e5))),
      error = function(e2) { cat("  Fallback also failed.\n"); NULL })
  }
)

if (!is.null(fit_3way)) {
  coef_a <- as.data.table(coef(summary(fit_3way)), keep.rownames = "term")
  setnames(coef_a,
           intersect(c("Estimate","Std. Error","z value","Pr(>|z|)"), names(coef_a)),
           intersect(c("estimate","SE","z","p"), c("estimate","SE","z","p")),
           skip_absent = TRUE)
  coef_a[, OR  := exp(estimate)]
  coef_a[, sig := fcase(p < 0.001,"***", p < 0.01,"**",
                         p < 0.05,"*", p < 0.10,".",
                         default = "ns")]

  cat("Three-way GLMM — key terms:\n")
  key_terms <- coef_a[str_detect(term, "virus|C6ST1|ADAMTS15")]
  print(key_terms[, .(term, estimate = round(estimate,4),
                       SE = round(SE,4), z = round(z,3),
                       p = round(p,5), OR = round(OR,3), sig)])

  fwrite(coef_a, file.path(OUT_DIR, "partA_glmm_coefficients.csv"))
  cat("  Saved partA_glmm_coefficients.csv\n\n")

  # ── LRT: does three-way interaction improve fit? ──────────────────────────
  cat("LRT: three-way interaction vs two-way additive model\n")
  fit_2way <- tryCatch(
    update(fit_3way, . ~ ADAMTS15_present * C6ST1_present +
             ADAMTS15_present * virus_factor +
             C6ST1_present * virus_factor +
             (1 | mouse_id)),
    error = function(e) NULL
  )
  if (!is.null(fit_2way)) {
    lrt_a <- anova(fit_2way, fit_3way)
    cat("\nLRT result:\n")
    print(lrt_a)
    fwrite(as.data.table(lrt_a, keep.rownames = "model"),
           file.path(OUT_DIR, "partA_lrt_threeway.csv"))
  }
} else {
  cat("  GLMM failed — skipping frequentist Part A.\n\n")
}

# ── Bayesian three-way model ──────────────────────────────────────────────────
cat("\nFitting Bayesian three-way model (brms) ...\n")
cat("  Weakly informative priors; 4 chains × 4000 iterations, warmup 2000\n\n")

# Section-level aggregation for Bayesian model (more tractable than cell-level)
# Aggregate per mouse × slice × ADAMTS15_present × C6ST1_present × virus_factor
sec_fact <- fact[, .(
  frac_enwrapped = mean(enwrapped),
  n_pv           = .N
), by = .(mouse_id, treatment, slice_id,
          ADAMTS15_present, C6ST1_present, virus_factor)]
sec_fact <- sec_fact[n_pv >= MIN_PV_SECTION]
sec_fact[, virus_num := as.integer(virus_factor == "virus+")]

cat(sprintf("  Section-level records for Bayes model: %d\n\n",
            nrow(sec_fact)))

priors_a <- c(
  prior(normal(0, 1), class = b),
  prior(normal(0, 1), class = sd),
  prior(normal(0, 1), class = sigma)
)

m_bayes_3way <- brm(
  frac_enwrapped ~
    ADAMTS15_present * C6ST1_present * virus_num +
    (1 | mouse_id),
  data    = sec_fact,
  prior   = priors_a,
  chains  = 4L, iter  = 4000L, warmup = 2000L,
  cores   = 4L, seed  = 2025L,
  control = list(adapt_delta = 0.95),
  file    = file.path(OUT_DIR, "m_bayes_3way"),
  silent  = 2L
)

draws_3w <- as_draws_df(m_bayes_3way)

# Three-way interaction term
inter3_col <- grep("ADAMTS15_present:C6ST1_present:virus_num",
                   names(draws_3w), value = TRUE)[1]

if (!is.na(inter3_col)) {
  inter3_draws <- draws_3w[[inter3_col]]
  # Negative three-way interaction means: adding virus_num (going from virus-
  # to virus+) produces a LARGER reduction when both enzymes are present
  # than when summing individual enzyme-specific gaps.
  p_synergy_local <- mean(inter3_draws < 0)

  cat(sprintf("Three-way interaction (Bayesian):\n"))
  cat(sprintf("  mean = %.4f, 89%% CI [%.4f, %.4f]\n",
              mean(inter3_draws),
              quantile(inter3_draws, 0.055),
              quantile(inter3_draws, 0.945)))
  cat(sprintf("  P(three-way interaction < 0) = %.4f\n", p_synergy_local))
  cat(sprintf("  Interpretation: %s\n\n",
              ifelse(p_synergy_local > 0.95,
                     "Strong evidence for enhanced local effect in combination",
                     ifelse(p_synergy_local > 0.80,
                            "Moderate evidence for enhanced local effect",
                            "Weak evidence for enhanced local effect"))))

  bayes_3w_summary <- data.table(
    parameter        = c("ADAMTS15_present","C6ST1_present","virus_num",
                         "A15×C6ST1","A15×virus","C6ST1×virus",
                         "A15×C6ST1×virus (key)"),
    col              = c(
      grep("^b_ADAMTS15_present$",       names(draws_3w), value=TRUE)[1],
      grep("^b_C6ST1_present$",          names(draws_3w), value=TRUE)[1],
      grep("^b_virus_num$",              names(draws_3w), value=TRUE)[1],
      grep("ADAMTS15_present:C6ST1_present$", names(draws_3w), value=TRUE)[1],
      grep("ADAMTS15_present:virus_num$",     names(draws_3w), value=TRUE)[1],
      grep("C6ST1_present:virus_num$",        names(draws_3w), value=TRUE)[1],
      inter3_col
    )
  )

  bayes_3w_summary[, `:=`(
    post_mean = sapply(col, function(cn)
      if (!is.na(cn)) round(mean(draws_3w[[cn]]), 4) else NA_real_),
    ci_89_lo  = sapply(col, function(cn)
      if (!is.na(cn)) round(quantile(draws_3w[[cn]], 0.055), 4) else NA_real_),
    ci_89_hi  = sapply(col, function(cn)
      if (!is.na(cn)) round(quantile(draws_3w[[cn]], 0.945), 4) else NA_real_),
    p_negative = sapply(col, function(cn)
      if (!is.na(cn)) round(mean(draws_3w[[cn]] < 0), 4) else NA_real_)
  )]
  bayes_3w_summary[, col := NULL]

  cat("Bayesian three-way model — posterior summaries:\n")
  print(bayes_3w_summary)
  fwrite(bayes_3w_summary, file.path(OUT_DIR, "partA_bayesian_summary.csv"))
} else {
  cat("  Could not identify three-way interaction column in draws.\n\n")
}

# ── Figure A: posterior of three-way interaction ──────────────────────────────
if (exists("inter3_draws")) {
  pA_post <- ggplot(data.frame(x = inter3_draws), aes(x = x)) +
    geom_density(fill = "#f1c44e", alpha = 0.6, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    annotate("text",
             x     = quantile(inter3_draws, 0.05),
             y     = Inf,
             label = sprintf("P(interaction < 0) = %.4f\n(negative = combination amplifies local effect)",
                             p_synergy_local),
             hjust = 0, vjust = 1.4, size = 3.5) +
    labs(
      title    = "Part A: Posterior — ADAMTS15 × C6ST1 × virus proximity interaction",
      subtitle = "Negative = combination produces larger local PNN depletion than sum of parts",
      x        = "Posterior: three-way interaction coefficient",
      y        = "Density"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(OUT_DIR, "fig_partA_threeway_posterior.pdf"),
         pA_post, width = 8, height = 5)
  cat("  Saved fig_partA_threeway_posterior.pdf\n\n")
}

# =============================================================================
# PART B — Direct virus+ comparison: C6ST1_ADAMTS15 vs ADAMTS15
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part B: Virus+ cells — C6ST1_ADAMTS15 vs ADAMTS15\n")
cat("══════════════════════════════════════════════════════\n\n")
cat("  Pre-specified contrast, no correction (mirrors script 20/21 approach)\n\n")

# Section-level aggregation: virus+ cells, Core zone, two groups
sec_pos_rq3 <- tags[
  virus_positive == TRUE &
  zone == "Core" &
  treatment %in% c("ADAMTS15","C6ST1_ADAMTS15"), .(
  n_pv           = .N,
  frac_enwrapped = mean(enwrapped)
), by = .(mouse_id, treatment, slice_id)]
sec_pos_rq3 <- sec_pos_rq3[n_pv >= MIN_PV_SECTION]
sec_pos_rq3[, treatment := factor(treatment,
                                   levels = c("ADAMTS15","C6ST1_ADAMTS15"))]

cat(sprintf("Virus+ sections per group (Core, >= %d cells):\n", MIN_PV_SECTION))
print(sec_pos_rq3[, .N, by = treatment])

fit_b <- lmer(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data    = sec_pos_rq3,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

if (isSingular(fit_b)) cat("  *** SINGULAR FIT: Part B model ***\n")

em_b  <- emmeans(fit_b, ~ treatment)
ct_b  <- as.data.table(contrast(em_b, method = "pairwise", adjust = "none"))
ct_b[, sig := fcase(p.value < 0.001,"***", p.value < 0.01,"**",
                     p.value < 0.05,"*", p.value < 0.10,".",
                     default = "ns")]

cat("\nPre-specified contrast: C6ST1_ADAMTS15 vs ADAMTS15 (virus+ cells):\n")
print(ct_b[, .(contrast, estimate = round(estimate,4),
                SE = round(SE,4), df = round(df,1),
                t.ratio = round(t.ratio,3), p.value = round(p.value,4), sig)])

fwrite(ct_b, file.path(OUT_DIR, "partB_viruspos_contrast.csv"))

# ── Distance-bin profile for virus+ cells, two groups + mScarlet baseline ────
cat("\nComputing distance-to-C3 for Part B distance profile ...\n")

rq3_tags <- tags[treatment %in% c("ADAMTS15","C6ST1_ADAMTS15") & zone == "Core"]
rq3_tags[, dist_to_virus_px := NA_real_]
slice_keys_b <- unique(rq3_tags[, .(mouse_id, slice_id)])

for (i in seq_len(nrow(slice_keys_b))) {
  mid <- slice_keys_b$mouse_id[i]
  sid <- slice_keys_b$slice_id[i]
  pv_sl <- rq3_tags[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  c3_sl <- c3[mouse_id == mid & slice_id == sid, .(x_c3, y_c3)]
  if (nrow(pv_sl) == 0L || nrow(c3_sl) == 0L) next
  nn  <- nn2(data  = as.matrix(c3_sl),
             query = as.matrix(pv_sl), k = 1L)
  idx <- which(rq3_tags$mouse_id == mid & rq3_tags$slice_id == sid)
  rq3_tags[idx, dist_to_virus_px := nn$nn.dists[, 1L]]
}

rq3_tags[, dist_bin := cut(dist_to_virus_px,
                             breaks = DIST_BREAKS,
                             labels = DIST_LABELS,
                             right  = FALSE,
                             include.lowest = TRUE)]

binned_b <- rq3_tags[!is.na(dist_bin) & virus_positive == TRUE, .(
  frac_enwrapped = mean(enwrapped),
  n_cells        = .N
), by = .(treatment, dist_bin)]
binned_b <- binned_b[n_cells >= MIN_PV_BIN]
binned_b[, treatment := factor(treatment,
                                levels = c("ADAMTS15","C6ST1_ADAMTS15"))]
binned_b[, dist_bin  := factor(dist_bin, levels = DIST_LABELS)]

cat("\nVirus+ enwrapment by distance bin — ADAMTS15 vs C6ST1_ADAMTS15:\n")
print(dcast(binned_b, dist_bin ~ treatment, value.var = "frac_enwrapped"))

fwrite(binned_b, file.path(OUT_DIR, "partB_viruspos_distance_bins.csv"))

# ── Figure B ─────────────────────────────────────────────────────────────────
pB_bin <- ggplot(binned_b,
                 aes(x = dist_bin, y = frac_enwrapped,
                     colour = treatment, group = treatment)) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(size = n_cells)) +
  scale_colour_manual(values = PALETTE_RQ3) +
  scale_size_continuous(range = c(1.5, 4), name = "n cells") +
  labs(
    title    = "Part B: Virus+ cell enwrapment vs distance — ADAMTS15 vs C6ST1_ADAMTS15",
    subtitle = "Lower at 0-30px in C6ST1_ADAMTS15 = enhanced local substrate conditioning",
    x        = "Distance to nearest transduced cell",
    y        = "Fraction enwrapped (virus+ cells only)",
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 30, hjust = 1),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig_partB_viruspos_distance.pdf"),
       pB_bin, width = 8, height = 5)
cat("  Saved fig_partB_viruspos_distance.pdf\n\n")

# =============================================================================
# PART C — Diffuse enhancement: C6ST1_ADAMTS15 vs ADAMTS15, virus- cells
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part C: Virus- cells — C6ST1_ADAMTS15 vs ADAMTS15\n")
cat("══════════════════════════════════════════════════════\n\n")
cat("  Tests whether C6ST1 co-expression enhances the DIFFUSE extracellular\n")
cat("  effect of ADAMTS15 on non-adjacent cells.\n\n")

sec_neg_rq3 <- tags[
  virus_positive == FALSE &
  zone == "Core" &
  treatment %in% c("ADAMTS15","C6ST1_ADAMTS15"), .(
  n_pv           = .N,
  frac_enwrapped = mean(enwrapped)
), by = .(mouse_id, treatment, slice_id)]
sec_neg_rq3 <- sec_neg_rq3[n_pv >= MIN_PV_SECTION]
sec_neg_rq3[, treatment := factor(treatment,
                                   levels = c("ADAMTS15","C6ST1_ADAMTS15"))]

cat(sprintf("Virus- sections per group (Core, >= %d cells):\n", MIN_PV_SECTION))
print(sec_neg_rq3[, .N, by = treatment])

fit_c <- lmer(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data    = sec_neg_rq3,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

if (isSingular(fit_c)) cat("  *** SINGULAR FIT: Part C model ***\n")

em_c  <- emmeans(fit_c, ~ treatment)
ct_c  <- as.data.table(contrast(em_c, method = "pairwise", adjust = "none"))
ct_c[, sig := fcase(p.value < 0.001,"***", p.value < 0.01,"**",
                     p.value < 0.05,"*", p.value < 0.10,".",
                     default = "ns")]

cat("\nPre-specified contrast: C6ST1_ADAMTS15 vs ADAMTS15 (virus- cells):\n")
print(ct_c[, .(contrast, estimate = round(estimate,4),
                SE = round(SE,4), df = round(df,1),
                t.ratio = round(t.ratio,3), p.value = round(p.value,4), sig)])

fwrite(ct_c, file.path(OUT_DIR, "partC_virusneg_contrast.csv"))

# ── Bayesian: P(C6ST1_ADAMTS15 virus- < ADAMTS15 virus-) ─────────────────────
cat("\nFitting Bayesian model for virus- comparison (Part C) ...\n")

priors_c <- c(
  prior(normal(0, 1), class = b),
  prior(normal(0, 1), class = sd),
  prior(normal(0, 1), class = sigma)
)

all_neg_rq3 <- rbind(sec_neg_rq3, sec_pos_rq3, fill = TRUE,
                      use.names = TRUE)
all_neg_rq3[, virus_factor := fifelse(
  !is.na(frac_enwrapped) & mouse_id %in% sec_neg_rq3$mouse_id,
  "virus-","virus+")]

m_bayes_c <- brm(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data    = sec_neg_rq3,
  prior   = priors_c,
  chains  = 4L, iter = 4000L, warmup = 2000L,
  cores   = 4L, seed = 2025L,
  control = list(adapt_delta = 0.95),
  file    = file.path(OUT_DIR, "m_bayes_virusneg"),
  silent  = 2L
)

draws_c <- as_draws_df(m_bayes_c)
col_c6a15 <- grep("treatmentC6ST1_ADAMTS15", names(draws_c), value = TRUE)[1]

if (!is.na(col_c6a15)) {
  inter_c_draws <- draws_c[[col_c6a15]]
  p_c_neg       <- mean(inter_c_draws < 0)
  cat(sprintf("\nBayesian: P(C6ST1_ADAMTS15 virus- < ADAMTS15 virus-) = %.4f\n",
              p_c_neg))
  cat(sprintf("  mean = %.4f, 89%% CI [%.4f, %.4f]\n",
              mean(inter_c_draws),
              quantile(inter_c_draws, 0.055),
              quantile(inter_c_draws, 0.945)))
  cat(sprintf("  Interpretation: %s\n\n",
              ifelse(p_c_neg > 0.80,
                     "C6ST1 co-expression enhances diffuse ADAMTS15 effect",
                     "No evidence for diffuse enhancement by C6ST1")))

  bayes_c_summary <- data.table(
    comparison     = "C6ST1_ADAMTS15 vs ADAMTS15 (virus- cells)",
    post_mean      = round(mean(inter_c_draws), 4),
    ci_89_lo       = round(quantile(inter_c_draws, 0.055), 4),
    ci_89_hi       = round(quantile(inter_c_draws, 0.945), 4),
    p_negative     = round(p_c_neg, 4),
    interpretation = ifelse(p_c_neg > 0.80,
                            "diffuse enhancement present",
                            "diffuse enhancement absent")
  )
  fwrite(bayes_c_summary, file.path(OUT_DIR, "partC_bayesian_summary.csv"))
}

# =============================================================================
# COMBINED FIGURE — Parts B and C side by side
# =============================================================================
cat("Saving combined RQ3 summary figure ...\n")

# Animal-level dot plot: frac_enwrapped per animal × virus status × group
animal_rq3 <- tags[
  treatment %in% c("ADAMTS15","C6ST1_ADAMTS15") & zone == "Core", .(
  frac_enwrapped = mean(enwrapped)
), by = .(mouse_id, treatment, virus_positive)]

animal_rq3[, virus_label := fifelse(virus_positive, "virus+","virus-")]
animal_rq3[, treatment   := factor(treatment,
                                    levels = c("ADAMTS15","C6ST1_ADAMTS15"))]
animal_rq3[, group_label := paste0(treatment,"\n",virus_label)]

pBC_dot <- ggplot(animal_rq3,
                  aes(x = group_label, y = frac_enwrapped,
                      colour = treatment)) +
  geom_jitter(width = 0.08, size = 3.5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.35, colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = PALETTE_RQ3, guide = "none") +
  scale_x_discrete(limits = c("ADAMTS15\nvirus-","ADAMTS15\nvirus+",
                               "C6ST1_ADAMTS15\nvirus-",
                               "C6ST1_ADAMTS15\nvirus+")) +
  geom_hline(yintercept = 0, colour = "grey80") +
  labs(
    title    = "RQ3 cell-level: enwrapment by virus status — ADAMTS15 vs C6ST1_ADAMTS15",
    subtitle = "Parts B (virus+) and C (virus-) side by side | Core zone | crossbar = animal mean",
    x        = NULL,
    y        = "Fraction enwrapped"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(size = 9),
        panel.grid.minor = element_blank())

p_rq3_combined <- pBC_dot / pB_bin +
  plot_annotation(
    title    = "Script 26: Cell-level RQ3 summary",
    subtitle = "Top: animal means by virus status | Bottom: distance profile (virus+ cells only)"
  )

ggsave(file.path(OUT_DIR, "fig_rq3_combined.pdf"),
       p_rq3_combined, width = 10, height = 10)
cat("  Saved fig_rq3_combined.pdf\n\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Script 26 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  partA_glmm_coefficients.csv\n")
cat("  partA_lrt_threeway.csv\n")
cat("  partA_bayesian_summary.csv\n")
cat("  partB_viruspos_contrast.csv\n")
cat("  partB_viruspos_distance_bins.csv\n")
cat("  partC_virusneg_contrast.csv\n")
cat("  partC_bayesian_summary.csv\n")
cat("  fig_partA_threeway_posterior.pdf\n")
cat("  fig_partB_viruspos_distance.pdf\n")
cat("  fig_rq3_combined.pdf\n")
cat("══════════════════════════════════════════════════════\n")
