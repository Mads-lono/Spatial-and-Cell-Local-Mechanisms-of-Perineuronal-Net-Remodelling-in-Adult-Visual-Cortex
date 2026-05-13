# =============================================================================
# Script 24: Cell-Level Transduction Analysis
# =============================================================================
# Research question (addendum to RQ1 and RQ2):
#   Do PV cells that are directly adjacent to transduced (mScarlet+) cells
#   show lower PNN enwrapment than non-adjacent PV cells within the same
#   section and animal?
#
# This is a within-animal, within-section comparison that provides
# cell-level evidence for the treatment effects established at group level
# in scripts 01/19. The mScarlet control group is the critical internal
# negative control: virus+ PV cells in mScarlet animals carry reporter only,
# so any enwrapment difference in that group cannot be enzymatic in origin.
#
# Four analyses:
#   Part A — Triple colocalization summary (descriptive)
#             Per animal × zone: proportions of PV cells in each (virus, PNN)
#             status category. Straightforward and maximally interpretable.
#
#   Part B — Section-level LMM (primary frequentist test)
#             frac_enwrapped ~ virus_positive * treatment + (1|mouse_id)
#             Key test: treatment × virus_positive interaction.
#             Is the virus+/virus- enwrapment gap larger in enzyme groups
#             than in mScarlet controls?
#
#   Part C — Cell-level GLMM (logistic)
#             enwrapped ~ virus_positive * treatment + (1|mouse_id/slice_id)
#             Uses all cell-level information without aggregation.
#             Reports OR for virus_positive per treatment group.
#
#   Part D — Bayesian section-level model
#             brms version of Part B.
#             Reports P(virus+ enwrapment < virus- enwrapment) per treatment.
#
# Inputs:
#   cells_with_zones.csv                        — PV/PNN cell coordinates
#   /media/.../Counts_C3-*.csv                  — mScarlet cell coordinates
#
# Outputs: results/24_virus_cell_enwrapment/
#   cell_virus_tags.csv             — per PV cell: virus_positive flag
#   section_virus_enwrapment.csv    — section-level aggregation (Part A/B)
#   partA_triple_colocalization.csv — descriptive proportions
#   partB_lmm_results.csv           — section LMM contrasts
#   partB_interaction_lrt.csv       — LRT for interaction term
#   partC_glmm_results.csv          — cell-level GLMM odds ratios
#   partD_bayesian_results.csv      — posterior summaries
#   fig_virus_enwrapment_*.pdf      — figures
#
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
RESULTS_DIR <- "/path/to/results"
OUT_DIR     <- file.path(RESULTS_DIR, "24_virus_cell_enwrapment")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
COLOC_THRESH    <- 30L      # pixels — same threshold throughout pipeline
VIRUS_THRESH    <- 30L      # pixels — PV is virus+ if nearest C3 cell ≤ this
MIN_PV_SLICE    <- 5L       # min PV cells per slice × hemisphere
MIN_PV_GROUP    <- 3L       # min PV cells per virus-status group per section
EXCLUDE         <- ""       # exclude injection failure animal
TREATMENT_ORDER <- c("mScarlet","ADAMTS4","ADAMTS4_MD",
                     "ADAMTS15","C6ST1","C6ST1_ADAMTS15")
VISUAL_RE       <- regex("isual", ignore_case = TRUE)

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

cat("Script 24: Cell-level transduction analysis\n")
cat("============================================================\n\n")


# =============================================================================
# STEP 1: Load cells_with_zones — PV and PNN coordinates
# =============================================================================
cat("Loading cells_with_zones.csv ...\n")

LOAD_COLS <- c("mouse_id","slice_id","cell_type","hemisphere",
               "x_hires","y_hires","brain_area","zone","treatment")

dt <- fread(CELLS_CSV, select = LOAD_COLS)
dt <- fix_names(dt)
dt <- dt[mouse_id != EXCLUDE]
dt <- dt[str_detect(brain_area, VISUAL_RE)]
dt <- dt[treatment %in% TREATMENT_ORDER]
dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

# Ipsilateral hemisphere only (left = injected side)
dt_ipsi <- dt[hemisphere == "left"]

pv  <- dt_ipsi[cell_type == "PV"]
pnn <- dt_ipsi[cell_type == "PNN"]

cat(sprintf("  PV cells (ipsi, visual): %s\n", format(nrow(pv),  big.mark = ",")))
cat(sprintf("  PNN cells (ipsi, visual): %s\n\n", format(nrow(pnn), big.mark = ",")))

# =============================================================================
# STEP 2: Load C3 (mScarlet) cell coordinates from raw Counts_C3 files
# =============================================================================
cat("Loading C3 mScarlet cell coordinates from Counts_C3 files ...\n")

# Filename format: Counts_C3-{treatment}_{animal_num}_{final_slice}.csv
# final_slice (e.g. s001) = slice_id in cells_with_zones

parse_c3_name <- function(fname) {
  # Greedy match on treatment to handle underscores (e.g. C6ST1_ADAMTS15)
  m <- regmatches(fname, regexec("Counts_C3-(.+)_(\\d+)_(s\\d+)\\.csv", fname))[[1]]
  if (length(m) != 4L) return(NULL)
  list(treatment_raw = m[2], animal_num = m[3], slice_id = m[4])
}

c3_files <- list.files(C3_DIR, pattern = "^Counts_C3-.*\\.csv$",
                       full.names = TRUE)
cat(sprintf("  Found %d Counts_C3 files\n", length(c3_files)))

c3_list <- vector("list", length(c3_files))

for (i in seq_along(c3_files)) {
  meta <- parse_c3_name(basename(c3_files[i]))
  if (is.null(meta)) next

  df <- tryCatch(fread(c3_files[i], select = c("Global_X", "Global_Y")),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) next
  if (!all(c("Global_X","Global_Y") %in% names(df))) next

  df[, `:=`(
    treatment_raw = meta$treatment_raw,
    animal_num    = meta$animal_num,
    slice_id      = meta$slice_id
  )]
  setnames(df, c("Global_X","Global_Y"), c("x_hires","y_hires"))
  c3_list[[i]] <- df
}

c3 <- rbindlist(c3_list, fill = TRUE)
c3[, mouse_id := paste0(
  str_replace(treatment_raw, "C6ST1_ADAMTS4$", "C6ST1_ADAMTS15"),
  "_", animal_num)]

# Apply standard corrections
c3[mouse_id == "mScarlet_4",   treatment_raw := "ADAMTS4_MD"]
c3[mouse_id == "ADAMTS4_MD_4", treatment_raw := "mScarlet"]
c3 <- c3[mouse_id != EXCLUDE]
c3[, treatment := treatment_raw]
c3[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("  mScarlet cells loaded: %s across %d animals\n\n",
            format(nrow(c3), big.mark = ","),
            uniqueN(c3$mouse_id)))

# =============================================================================
# STEP 3: Tag each PV cell as virus_positive / virus_negative
# =============================================================================
cat("Tagging PV cells as virus+ / virus- via nn2 (COLOC_THRESH =",
    VIRUS_THRESH, "px) ...\n")

# Join keys: mouse_id × slice_id
slice_keys <- unique(pv[, .(mouse_id, slice_id)])

pv[, virus_positive := FALSE]

n_tagged <- 0L

for (i in seq_len(nrow(slice_keys))) {
  mid <- slice_keys$mouse_id[i]
  sid <- slice_keys$slice_id[i]

  pv_sl <- pv[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  c3_sl <- c3[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]

  if (nrow(pv_sl) == 0L || nrow(c3_sl) == 0L) next

  nn   <- nn2(data  = as.matrix(c3_sl),
              query = as.matrix(pv_sl),
              k     = 1L)
  flag <- nn$nn.dists[, 1L] <= VIRUS_THRESH

  # Write back into pv
  idx <- which(pv$mouse_id == mid & pv$slice_id == sid)
  pv[idx, virus_positive := flag]
  n_tagged <- n_tagged + sum(flag)
}

pv_tagged <- pv[, .(mouse_id, slice_id, treatment, zone,
                    x_hires, y_hires, virus_positive)]

cat(sprintf("  Tagged %s virus+ PV cells (%.1f%% of all ipsi-visual PV cells)\n\n",
            format(n_tagged, big.mark = ","),
            100 * n_tagged / nrow(pv)))

fwrite(pv_tagged, file.path(OUT_DIR, "cell_virus_tags.csv"))

# =============================================================================
# STEP 4: Compute enwrapment per PV cell using same nn2 logic as script 00
# =============================================================================
cat("Computing per-cell enwrapment ...\n")

pv_tagged[, enwrapped := FALSE]

enwrap_keys <- unique(pv_tagged[, .(mouse_id, slice_id)])

for (i in seq_len(nrow(enwrap_keys))) {
  mid <- enwrap_keys$mouse_id[i]
  sid <- enwrap_keys$slice_id[i]

  pv_sl  <- pv_tagged[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  pnn_sl <- pnn[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]

  if (nrow(pv_sl) == 0L || nrow(pnn_sl) == 0L) next

  nn   <- nn2(data  = as.matrix(pnn_sl),
              query = as.matrix(pv_sl),
              k     = 1L)
  flag <- nn$nn.dists[, 1L] <= COLOC_THRESH

  idx <- which(pv_tagged$mouse_id == mid & pv_tagged$slice_id == sid)
  pv_tagged[idx, enwrapped := flag]
}

cat(sprintf("  Overall ipsi enwrapment: %.1f%%\n\n",
            100 * mean(pv_tagged$enwrapped)))

# =============================================================================
# STEP 5: Aggregate to section level
# =============================================================================
# For each mouse × slice × virus_status: compute frac_enwrapped
# Require at least MIN_PV_GROUP cells per group per section

section_dt <- pv_tagged[, .(
  n_pv         = .N,
  frac_enwrapped = mean(enwrapped)
), by = .(mouse_id, treatment, slice_id, zone, virus_positive)]

section_dt <- section_dt[n_pv >= MIN_PV_GROUP]
section_dt[, virus_label := fifelse(virus_positive, "virus+", "virus-")]

fwrite(section_dt, file.path(OUT_DIR, "section_virus_enwrapment.csv"))

cat(sprintf("Section-level records (min %d cells/group): %d\n\n",
            MIN_PV_GROUP, nrow(section_dt)))

# =============================================================================
# PART A — Triple colocalization summary
# =============================================================================
cat("── Part A: Triple colocalization summary ─────────────────────────────────\n")

triple <- pv_tagged[, .(
  n_total          = .N,
  n_viruspos_enwrap  = sum(virus_positive  & enwrapped),
  n_viruspos_noenwrap= sum(virus_positive  & !enwrapped),
  n_virusneg_enwrap  = sum(!virus_positive & enwrapped),
  n_virusneg_noenwrap= sum(!virus_positive & !enwrapped)
), by = .(mouse_id, treatment, zone)]

triple[, `:=`(
  pct_viruspos_enwrap   = 100 * n_viruspos_enwrap   / n_total,
  pct_viruspos_noenwrap = 100 * n_viruspos_noenwrap / n_total,
  pct_virusneg_enwrap   = 100 * n_virusneg_enwrap   / n_total,
  pct_virusneg_noenwrap = 100 * n_virusneg_noenwrap / n_total
)]

cat("Triple colocalization summary (% of all ipsi-visual PV cells per animal):\n")
print(triple[zone == "Core",
             .(mouse_id, treatment,
               pct_viruspos_enwrap   = round(pct_viruspos_enwrap, 1),
               pct_viruspos_noenwrap = round(pct_viruspos_noenwrap, 1),
               pct_virusneg_enwrap   = round(pct_virusneg_enwrap, 1),
               pct_virusneg_noenwrap = round(pct_virusneg_noenwrap, 1))])

fwrite(triple, file.path(OUT_DIR, "partA_triple_colocalization.csv"))

# Figure A: stacked bar per treatment × virus status (Core zone)
plot_triple <- melt(
  triple[zone == "Core"],
  id.vars       = c("mouse_id","treatment"),
  measure.vars  = c("pct_viruspos_enwrap","pct_viruspos_noenwrap",
                    "pct_virusneg_enwrap","pct_virusneg_noenwrap"),
  variable.name = "category",
  value.name    = "pct"
)
plot_triple[, category := factor(category,
  levels = c("pct_viruspos_enwrap","pct_viruspos_noenwrap",
             "pct_virusneg_enwrap","pct_virusneg_noenwrap"),
  labels = c("Virus+ / enwrapped","Virus+ / not enwrapped",
             "Virus- / enwrapped","Virus- / not enwrapped"))]
plot_triple[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

pA <- ggplot(plot_triple, aes(x = mouse_id, y = pct, fill = category)) +
  geom_col(position = "stack") +
  facet_wrap(~ treatment, scales = "free_x", nrow = 2) +
  scale_fill_manual(
    values = c("Virus+ / enwrapped"     = "#2c7bb6",
               "Virus+ / not enwrapped" = "#abd9e9",
               "Virus- / enwrapped"     = "#d7191c",
               "Virus- / not enwrapped" = "#fdae61")) +
  labs(title   = "Triple colocalization: PV × virus × PNN (Core zone, ipsilateral)",
       x       = NULL,
       y       = "% of all ipsi PV cells",
       fill    = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_partA_triple_colocalization.pdf"),
       pA, width = 12, height = 6)
cat("  Saved fig_partA_triple_colocalization.pdf\n\n")

# =============================================================================
# PART B — Section-level LMM: treatment × virus_positive interaction
# =============================================================================
cat("── Part B: Section-level LMM ─────────────────────────────────────────────\n")

# Use Core zone — strongest expected signal
# Fit interaction model and compare to additive model via LRT

core_sec <- section_dt[zone == "Core"]
core_sec[, treatment    := factor(treatment, levels = TREATMENT_ORDER)]
core_sec[, virus_factor := factor(virus_positive, levels = c(FALSE, TRUE),
                                   labels = c("virus-","virus+"))]

cat("Core zone sections (virus+ and virus-):\n")
print(core_sec[, .N, by = .(treatment, virus_factor)])

# ── Full interaction model ────────────────────────────────────────────────────
fit_interact <- tryCatch(
  lmer(frac_enwrapped ~ treatment * virus_factor + (1 | mouse_id),
       data    = core_sec,
       REML    = FALSE,
       control = lmerControl(optimizer = "bobyqa")),
  error = function(e) { cat("  Interaction model failed:", e$message, "\n"); NULL }
)

# ── Additive model for LRT ───────────────────────────────────────────────────
fit_additive <- tryCatch(
  lmer(frac_enwrapped ~ treatment + virus_factor + (1 | mouse_id),
       data    = core_sec,
       REML    = FALSE,
       control = lmerControl(optimizer = "bobyqa")),
  error = function(e) { cat("  Additive model failed:", e$message, "\n"); NULL }
)

if (!is.null(fit_interact) && !is.null(fit_additive)) {
  lrt <- anova(fit_additive, fit_interact)
  cat("\nLRT: interaction vs additive model\n")
  print(lrt)
  fwrite(as.data.table(lrt, keep.rownames = "model"),
         file.path(OUT_DIR, "partB_interaction_lrt.csv"))

  # Re-fit with REML for estimates
  fit_interact_reml <- lmer(
    frac_enwrapped ~ treatment * virus_factor + (1 | mouse_id),
    data    = core_sec,
    REML    = TRUE,
    control = lmerControl(optimizer = "bobyqa"))

  # ── Contrasts: virus+/virus- gap per treatment ────────────────────────────
  em_b  <- emmeans(fit_interact_reml, ~ virus_factor | treatment)
  ct_b  <- as.data.table(contrast(em_b, method = "revpairwise", adjust = "none"))
  ct_b[, sig := fcase(p.value < 0.001,"***", p.value < 0.01,"**",
                       p.value < 0.05,"*", p.value < 0.10,".",
                       default = "ns")]

  cat("\nVirus+ vs virus- contrast per treatment (Core zone):\n")
  print(ct_b[, .(treatment, contrast, estimate = round(estimate,3),
                  SE = round(SE,3), df = round(df,1),
                  t.ratio = round(t.ratio,3), p.value = round(p.value,4), sig)])

  fwrite(ct_b, file.path(OUT_DIR, "partB_lmm_results.csv"))

  # ── Figure B: dot + CI forest ─────────────────────────────────────────────
  ct_b[, treatment := factor(treatment, levels = rev(TREATMENT_ORDER))]

  pB <- ggplot(ct_b, aes(x = estimate, y = treatment, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = estimate - 1.96*SE, xmax = estimate + 1.96*SE),
                  width = 0.25, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(x = estimate + abs(estimate)*0.15 + 0.02, label = sig),
              hjust = 0, size = 4, fontface = "bold", colour = "black") +
    scale_colour_manual(values = PALETTE, guide = "none") +
    labs(title    = "Virus+ vs virus- enwrapment gap per treatment (Core, section LMM)",
         subtitle = "Negative = virus+ cells less enwrapped than virus- cells",
         x        = "Estimated difference (virus+ − virus-)",
         y        = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(OUT_DIR, "fig_partB_lmm_contrasts.pdf"),
         pB, width = 8, height = 5)
  cat("  Saved fig_partB_lmm_contrasts.pdf\n\n")
} else {
  cat("  Skipping LMM figures due to model fitting failure.\n\n")
}

# =============================================================================
# PART C — Cell-level GLMM (logistic)
# =============================================================================
cat("── Part C: Cell-level GLMM (logistic) ────────────────────────────────────\n")
cat("  NOTE: this model uses all individual PV cells as observations.\n")
cat("  Runtime may be long — nested random effects (1|mouse_id/slice_id).\n\n")

# Restrict to Core zone for consistency with Part B
core_cells <- pv_tagged[zone == "Core"]
core_cells[, treatment    := factor(treatment, levels = TREATMENT_ORDER)]
core_cells[, virus_factor := factor(virus_positive, levels = c(FALSE, TRUE),
                                     labels = c("virus-","virus+"))]
core_cells[, enwrapped_int := as.integer(enwrapped)]
# Nested grouping
core_cells[, slice_in_mouse := paste0(mouse_id,"_",slice_id)]

cat(sprintf("  Cell-level observations: %s\n\n",
            format(nrow(core_cells), big.mark = ",")))

fit_glmm <- tryCatch(
  glmer(enwrapped_int ~ treatment * virus_factor + (1 | mouse_id/slice_id),
        data   = core_cells,
        family = binomial(link = "logit"),
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl   = list(maxfun = 2e5))),
  error   = function(e) { cat("  GLMM failed:", e$message, "\n"); NULL },
  warning = function(w) {
    cat("  GLMM warning:", w$message, "\n")
    tryCatch(
      glmer(enwrapped_int ~ treatment * virus_factor + (1 | mouse_id),
            data   = core_cells,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl   = list(maxfun = 2e5))),
      error = function(e2) NULL)
  }
)

if (!is.null(fit_glmm)) {
  em_c   <- emmeans(fit_glmm, ~ virus_factor | treatment, type = "response")
  ct_c   <- as.data.table(contrast(em_c, method = "revpairwise", adjust = "none"))

  # Add OR from coefficient table
  coef_c <- as.data.table(coef(summary(fit_glmm)), keep.rownames = "term")
  setnames(coef_c, c("Estimate","Std. Error","z value","Pr(>|z|)"),
                   c("log_or","se","z","p"), skip_absent = TRUE)
  coef_c[, OR := exp(log_or)]

  ct_c[, sig := fcase(p.value < 0.001,"***", p.value < 0.01,"**",
                       p.value < 0.05,"*", p.value < 0.10,".",
                       default = "ns")]

  cat("Cell-level GLMM — virus+/virus- OR per treatment:\n")
  print(ct_c[, .(treatment, contrast, odds.ratio = round(odds.ratio,3),
                  SE = round(SE,3), z.ratio = round(z.ratio,3),
                  p.value = round(p.value,4), sig)])

  fwrite(ct_c,   file.path(OUT_DIR, "partC_glmm_results.csv"))
  fwrite(coef_c, file.path(OUT_DIR, "partC_glmm_coefficients.csv"))

  cat("  Saved partC_glmm_results.csv\n\n")
} else {
  cat("  GLMM could not be fitted. Skipping Part C output.\n\n")
}

# =============================================================================
# PART D — Bayesian section-level model
# =============================================================================
cat("── Part D: Bayesian section-level model ──────────────────────────────────\n")

priors_d <- c(
  prior(normal(0, 1),  class = b),
  prior(normal(0, 1),  class = sd),
  prior(normal(0, 1),  class = sigma)
)

m_bayes <- brm(
  frac_enwrapped ~ treatment * virus_factor + (1 | mouse_id),
  data    = core_sec,
  prior   = priors_d,
  chains  = 4L, iter = 4000L, warmup = 2000L,
  cores   = 4L, seed = 2025L,
  control = list(adapt_delta = 0.95),
  file    = file.path(OUT_DIR, "m_virus_bayes"),
  silent  = 2L
)

# ── Posterior predictions: virus+ vs virus- per treatment ─────────────────────
nd_d <- expand.grid(
  treatment   = TREATMENT_ORDER,
  virus_factor = c("virus-","virus+"),
  mouse_id     = NA,
  stringsAsFactors = FALSE
)
nd_d$treatment    <- factor(nd_d$treatment,   levels = TREATMENT_ORDER)
nd_d$virus_factor <- factor(nd_d$virus_factor, levels = c("virus-","virus+"))

epred_d <- posterior_epred(m_bayes, newdata = nd_d,
                            re_formula = NA, allow_new_levels = TRUE)

# For each treatment: P(virus+ enwrapment < virus- enwrapment)
bayes_results <- rbindlist(lapply(TREATMENT_ORDER, function(tx) {
  i_pos <- which(nd_d$treatment == tx & nd_d$virus_factor == "virus+")
  i_neg <- which(nd_d$treatment == tx & nd_d$virus_factor == "virus-")
  diff  <- epred_d[, i_pos] - epred_d[, i_neg]   # positive = virus+ > virus-
  data.table(
    treatment            = tx,
    post_mean_viruspos   = round(mean(epred_d[, i_pos]), 4),
    post_mean_virusneg   = round(mean(epred_d[, i_neg]), 4),
    post_mean_diff       = round(mean(diff), 4),
    ci_89_lo             = round(quantile(diff, 0.055), 4),
    ci_89_hi             = round(quantile(diff, 0.945), 4),
    # P that virus+ cells are LESS enwrapped (negative difference)
    p_viruspos_lower     = round(mean(diff < 0), 4)
  )
}))

cat("\nBayesian posterior: P(virus+ enwrapment < virus- enwrapment) per treatment:\n")
print(bayes_results[, .(treatment, post_mean_viruspos, post_mean_virusneg,
                         post_mean_diff, p_viruspos_lower)])

fwrite(bayes_results, file.path(OUT_DIR, "partD_bayesian_results.csv"))

# ── Figure D: posterior distribution of the gap per treatment ─────────────────
post_long <- rbindlist(lapply(TREATMENT_ORDER, function(tx) {
  i_pos <- which(nd_d$treatment == tx & nd_d$virus_factor == "virus+")
  i_neg <- which(nd_d$treatment == tx & nd_d$virus_factor == "virus-")
  data.table(
    treatment = tx,
    diff      = as.numeric(epred_d[, i_pos] - epred_d[, i_neg])
  )
}))
post_long[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

pD <- ggplot(post_long, aes(x = diff, fill = treatment, colour = treatment)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ treatment, scales = "free_y", ncol = 2) +
  scale_fill_manual(values   = PALETTE, guide = "none") +
  scale_colour_manual(values = PALETTE, guide = "none") +
  geom_text(
    data = bayes_results,
    aes(x = -Inf, y = Inf,
        label = sprintf("P(virus+ < virus-) = %.3f", p_viruspos_lower),
        colour = treatment),
    hjust = -0.05, vjust = 1.3, size = 3, inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  labs(
    title    = "Bayesian posteriors: virus+ minus virus- enwrapment (Core zone)",
    subtitle = "Negative = virus+ cells less enwrapped than virus- cells in same section",
    x        = "Posterior: frac_enwrapped(virus+) − frac_enwrapped(virus-)",
    y        = "Density"
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(OUT_DIR, "fig_partD_bayesian_posteriors.pdf"),
       pD, width = 10, height = 8)
cat("  Saved fig_partD_bayesian_posteriors.pdf\n\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Script 24 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  cell_virus_tags.csv\n")
cat("  section_virus_enwrapment.csv\n")
cat("  partA_triple_colocalization.csv\n")
cat("  partB_interaction_lrt.csv\n")
cat("  partB_lmm_results.csv\n")
cat("  partC_glmm_results.csv\n")
cat("  partC_glmm_coefficients.csv\n")
cat("  partD_bayesian_results.csv\n")
cat("  fig_partA_triple_colocalization.pdf\n")
cat("  fig_partB_lmm_contrasts.pdf\n")
cat("  fig_partD_bayesian_posteriors.pdf\n")
cat("══════════════════════════════════════════════════════\n")
