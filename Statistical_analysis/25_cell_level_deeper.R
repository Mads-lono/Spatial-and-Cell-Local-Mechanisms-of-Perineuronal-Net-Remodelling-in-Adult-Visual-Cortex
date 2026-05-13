# =============================================================================
# Script 25: Deeper cell-level analyses
# =============================================================================
# Three analyses motivated by script 24 results, each addressing existing RQs
# at the cell level:
#
#   Part A — Diffuse secretion test (RQ1 + RQ2)
#             Compares enwrapment of virus-NEGATIVE PV cells across treatment
#             groups against mScarlet virus-negative cells.
#             Logic: if ADAMTS15 acts by diffusing through the ECM after
#             secretion, virus- cells in ADAMTS15 sections should be as
#             depleted as virus+ cells. If C6ST1 acts locally, its virus-
#             cells should look like mScarlet.
#             Model: frac_enwrapped ~ treatment + (1|mouse_id)
#             Fitted separately on virus- and virus+ subsets; results compared.
#
#   Part B — Continuous distance-to-nearest-transduced-cell analysis (RQ2)
#             Replaces the binary 30px threshold from script 24 with a
#             continuous distance variable: Euclidean distance from each PV
#             cell to the nearest C3 (mScarlet) cell in the same slice.
#             Bins: 0-30px, 30-100px, 100-300px, 300-600px, 600px+
#             Plots mean enwrapment per bin per treatment.
#             Also fits a GLMM: enwrapped ~ log(dist_to_virus+1) * treatment
#             + (1|mouse_id/slice_id) to test whether the distance gradient
#             differs between enzyme groups and mScarlet.
#             Prediction: C6ST1 shows sharp local drop-off (local action);
#             ADAMTS15 shows flat profile (diffuse extracellular action).
#
#   Part C — ADAMTS4 paradox: direct test of the sampling artefact hypothesis
#             (RQ1)
#             Script 24 showed virus+ cells MORE enwrapped in ADAMTS4 sections.
#             Two competing explanations:
#               (i)  Sampling artefact: Core zone near bolus has denser PNNs
#                    → test: compare absolute enwrapment of virus- cells in
#                      ADAMTS4 sections vs mScarlet virus- cells
#               (ii) True biology: ADAMTS4 has no local effect
#                    → consistent with the null group-level result in script 01
#             Model: same section-level LMM as Part A, virus- subset only,
#             focused on ADAMTS4 vs mScarlet.
#             Also: within-ADAMTS4 correlation between virus+ cell density
#             (proxy for injection concentration) and local enwrapment.
#
# Inputs:
#   results/24_virus_cell_enwrapment/cell_virus_tags.csv  — per-cell tags
#   analysis_results/cells_with_zones.csv                 — PNN coords for nn2
#   /media/.../Counts_C3-*.csv                            — C3 coords
#
# Outputs: results/25_cell_level_deeper/
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
})

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV    <- "/path/to/analysis_results/cells_with_zones.csv"
C3_DIR       <- "/path/to/c3_results"
TAGS_CSV     <- "/path/to/results/24_virus_cell_enwrapment/cell_virus_tags.csv"
RESULTS_DIR  <- "/path/to/results"
OUT_DIR      <- file.path(RESULTS_DIR, "25_cell_level_deeper")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
COLOC_THRESH    <- 30L
MIN_PV_SECTION  <- 5L
MIN_PV_BIN      <- 10L      # minimum cells per distance bin for plotting
PX_TO_MM        <- 0.325 / 1000  # microns per pixel → mm
EXCLUDE         <- ""       # mouse_id to exclude if injection failed
TREATMENT_ORDER <- c("mScarlet","ADAMTS4","ADAMTS4_MD",
                     "ADAMTS15","C6ST1","C6ST1_ADAMTS15")
VISUAL_RE       <- regex("isual", ignore_case = TRUE)

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


cat("Script 25: Deeper cell-level analyses\n")
cat("============================================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# ── Cell virus tags from script 24 ────────────────────────────────────────────
# NOTE: enwrapment is recomputed here via nn2 against PNN coordinates from
# cells_with_zones.csv, consistent with the logic in script 24 Step 4.

cat("Loading cell_virus_tags.csv ...\n")
tags <- fread(TAGS_CSV)
tags <- fix_names(tags)
tags <- tags[mouse_id != EXCLUDE]
tags[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
cat(sprintf("  %s PV cells loaded\n", format(nrow(tags), big.mark = ",")))

# ── Load PNN coordinates to compute enwrapment ───────────────────────────────
cat("Loading PNN coordinates from cells_with_zones.csv ...\n")
pnn_load_cols <- c("mouse_id","slice_id","cell_type","hemisphere",
                   "x_hires","y_hires","brain_area","treatment")
cwz <- fread(CELLS_CSV, select = pnn_load_cols)
cwz <- fix_names(cwz)
cwz <- cwz[mouse_id != EXCLUDE]
cwz <- cwz[str_detect(brain_area, VISUAL_RE)]
cwz <- cwz[hemisphere == "left"]   # ipsilateral only
pnn_coords <- cwz[cell_type == "PNN"]
cat(sprintf("  %s PNN cells loaded\n\n", format(nrow(pnn_coords), big.mark = ",")))

# ── Compute enwrapped per PV cell via nn2 ─────────────────────────────────────
cat("Computing enwrapment for all PV cells in tags ...\n")
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

  idx <- which(tags$mouse_id == mid & tags$slice_id == sid)
  tags[idx, enwrapped := flag]
}

tags[, enwrapped_int := as.integer(enwrapped)]
cat(sprintf("  Overall enwrapment: %.1f%%\n\n",
            100 * mean(tags$enwrapped)))

# ── C3 coordinates (for Part B distance computation) ─────────────────────────
cat("Loading C3 coordinates ...\n")

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
c3[, treatment := str_replace(treatment_raw, "C6ST1_ADAMTS4$", "C6ST1_ADAMTS15")]
cat(sprintf("  %s C3 cells loaded\n\n", format(nrow(c3), big.mark = ",")))

# =============================================================================
# PART A — Diffuse secretion test
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part A: Diffuse secretion test\n")
cat("══════════════════════════════════════════════════════\n\n")

# Aggregate to section level separately for virus+ and virus- subsets
agg_virus <- function(dt, virus_flag) {
  sub <- dt[virus_positive == virus_flag & zone == "Core"]
  sec <- sub[, .(
    n_pv           = .N,
    frac_enwrapped = mean(enwrapped)
  ), by = .(mouse_id, treatment, slice_id)]
  sec[n_pv >= MIN_PV_SECTION]
}

sec_virusneg <- agg_virus(tags, FALSE)
sec_viruspos <- agg_virus(tags, TRUE)

cat(sprintf("Virus- sections (Core, >= %d cells): %d\n",
            MIN_PV_SECTION, nrow(sec_virusneg)))
cat(sprintf("Virus+ sections (Core, >= %d cells): %d\n\n",
            MIN_PV_SECTION, nrow(sec_viruspos)))

# ── LMM on virus- cells only ──────────────────────────────────────────────────
cat("LMM on virus-NEGATIVE PV cells only (Core zone):\n")
cat("  frac_enwrapped ~ treatment + (1|mouse_id)\n\n")

fit_virusneg <- lmer(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data    = sec_virusneg,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

if (isSingular(fit_virusneg))
  cat("  *** SINGULAR FIT: virus- model ***\n\n")

em_neg  <- emmeans(fit_virusneg, ~ treatment)
ct_neg  <- as.data.table(contrast(em_neg, method = "trt.vs.ctrl",
                                   ref = "mScarlet", adjust = "BH"))
ct_neg[, sig := fcase(p.value < 0.001,"***", p.value < 0.01,"**",
                       p.value < 0.05,"*", p.value < 0.10,".",
                       default = "ns")]
ct_neg[, subset := "virus-"]

cat("Contrasts vs mScarlet (virus- cells):\n")
print(ct_neg[, .(treatment = str_remove(contrast," - mScarlet"),
                  estimate = round(estimate,3), SE = round(SE,3),
                  df = round(df,1), t.ratio = round(t.ratio,3),
                  p.value = round(p.value,4), sig)])

# ── LMM on virus+ cells only ──────────────────────────────────────────────────
cat("\nLMM on virus-POSITIVE PV cells only (Core zone):\n")
cat("  frac_enwrapped ~ treatment + (1|mouse_id)\n\n")

fit_viruspos <- lmer(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data    = sec_viruspos,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

if (isSingular(fit_viruspos))
  cat("  *** SINGULAR FIT: virus+ model ***\n\n")

em_pos  <- emmeans(fit_viruspos, ~ treatment)
ct_pos  <- as.data.table(contrast(em_pos, method = "trt.vs.ctrl",
                                   ref = "mScarlet", adjust = "BH"))
ct_pos[, sig := fcase(p.value < 0.001,"***", p.value < 0.01,"**",
                       p.value < 0.05,"*", p.value < 0.10,".",
                       default = "ns")]
ct_pos[, subset := "virus+"]

cat("Contrasts vs mScarlet (virus+ cells):\n")
print(ct_pos[, .(treatment = str_remove(contrast," - mScarlet"),
                  estimate = round(estimate,3), SE = round(SE,3),
                  df = round(df,1), t.ratio = round(t.ratio,3),
                  p.value = round(p.value,4), sig)])

# ── Save and plot ─────────────────────────────────────────────────────────────
ct_all_a <- rbind(ct_neg, ct_pos)
ct_all_a[, treatment := factor(
  str_remove(contrast," - mScarlet"), levels = TREATMENT_ORDER[-1])]

fwrite(ct_all_a, file.path(OUT_DIR, "partA_diffuse_secretion_contrasts.csv"))

pA <- ggplot(ct_all_a, aes(x = estimate, y = treatment,
                             colour = treatment, shape = subset)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = estimate - 1.96*SE, xmax = estimate + 1.96*SE),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = PALETTE[-1], guide = "none") +
  scale_shape_manual(values = c("virus-" = 16, "virus+" = 17),
                     name = "PV cell subset") +
  facet_wrap(~ subset, ncol = 2) +
  labs(
    title    = "Part A: Enwrapment vs mScarlet — virus+ and virus- cells separately (Core)",
    subtitle = paste("Diffusion test: if ADAMTS15 secretes into ECM, virus- cells should",
                     "also be depleted (left panel)"),
    x        = "Estimated difference vs mScarlet (frac enwrapped)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

ggsave(file.path(OUT_DIR, "fig_partA_diffuse_secretion.pdf"),
       pA, width = 11, height = 5)
cat("\n  Saved fig_partA_diffuse_secretion.pdf\n\n")

# =============================================================================
# PART B — Continuous distance-to-nearest-transduced-cell
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part B: Continuous distance-to-virus analysis\n")
cat("══════════════════════════════════════════════════════\n\n")

cat("Computing distance from each PV cell to nearest C3 cell per slice ...\n")

# Work on Core-zone PV cells only (richest signal)
core_tags <- tags[zone == "Core"]
slice_keys_b <- unique(core_tags[, .(mouse_id, slice_id)])

core_tags[, dist_to_virus_px := NA_real_]

for (i in seq_len(nrow(slice_keys_b))) {
  mid <- slice_keys_b$mouse_id[i]
  sid <- slice_keys_b$slice_id[i]

  pv_sl <- core_tags[mouse_id == mid & slice_id == sid, .(x_hires, y_hires)]
  c3_sl <- c3[mouse_id == mid & slice_id == sid, .(x_c3, y_c3)]

  if (nrow(pv_sl) == 0L || nrow(c3_sl) == 0L) next

  nn  <- nn2(data  = as.matrix(c3_sl),
             query = as.matrix(pv_sl),
             k     = 1L)

  idx <- which(core_tags$mouse_id == mid & core_tags$slice_id == sid)
  core_tags[idx, dist_to_virus_px := nn$nn.dists[, 1L]]
}

core_tags[, dist_to_virus_mm := dist_to_virus_px * PX_TO_MM]
core_tags[, dist_bin := cut(dist_to_virus_px,
                             breaks = DIST_BREAKS,
                             labels = DIST_LABELS,
                             right  = FALSE,
                             include.lowest = TRUE)]

cat(sprintf("  %s cells with distance computed\n",
            format(sum(!is.na(core_tags$dist_to_virus_px)), big.mark=",")))

# ── Binned enwrapment per treatment ──────────────────────────────────────────
binned_b <- core_tags[!is.na(dist_bin), .(
  frac_enwrapped = mean(enwrapped),
  n_cells        = .N
), by = .(treatment, dist_bin)]

binned_b <- binned_b[n_cells >= MIN_PV_BIN]
binned_b[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
binned_b[, dist_bin  := factor(dist_bin,  levels = DIST_LABELS)]

fwrite(binned_b, file.path(OUT_DIR, "partB_binned_enwrapment.csv"))

cat("\nBinned enwrapment by distance to nearest C3 cell (Core zone):\n")
print(dcast(binned_b, dist_bin ~ treatment, value.var = "frac_enwrapped",
            fun.aggregate = mean))

# ── GLMM: distance gradient per treatment ─────────────────────────────────────
cat("\nFitting GLMM: enwrapped ~ log(dist+1) * treatment + (1|mouse_id/slice_id)\n")

core_tags[, log_dist := log(dist_to_virus_px + 1)]
core_tags_valid <- core_tags[!is.na(log_dist)]

fit_dist_glmm <- tryCatch(
  glmer(enwrapped_int ~ log_dist * treatment + (1 | mouse_id/slice_id),
        data    = core_tags_valid,
        family  = binomial(link = "logit"),
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl   = list(maxfun = 2e5))),
  error   = function(e) { cat("  Nested GLMM failed:", e$message, "\n"); NULL },
  warning = function(w) {
    cat("  Nested GLMM warning:", w$message, "\n  Retrying with (1|mouse_id) ...\n")
    tryCatch(
      glmer(enwrapped_int ~ log_dist * treatment + (1 | mouse_id),
            data    = core_tags_valid,
            family  = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl   = list(maxfun = 2e5))),
      error = function(e2) NULL)
  }
)

if (!is.null(fit_dist_glmm)) {
  coef_b <- as.data.table(coef(summary(fit_dist_glmm)), keep.rownames = "term")
  setnames(coef_b,
           intersect(c("Estimate","Std. Error","z value","Pr(>|z|)"), names(coef_b)),
           intersect(c("estimate","SE","z","p"), c("estimate","SE","z","p")),
           skip_absent = TRUE)

  # The key terms: log_dist:treatment — does distance-decay slope differ from mScarlet?
  interaction_terms <- coef_b[str_detect(term, "log_dist:treatment")]
  cat("\nDistance × treatment interaction terms (vs mScarlet reference):\n")
  print(interaction_terms[, .(
    term,
    estimate = round(estimate, 4),
    SE       = round(SE, 4),
    z        = round(z, 3),
    p        = round(p, 4)
  )])

  fwrite(coef_b, file.path(OUT_DIR, "partB_distance_glmm_coefficients.csv"))
  cat("  Saved partB_distance_glmm_coefficients.csv\n")
} else {
  cat("  GLMM could not be fitted. Binned plot still saved.\n")
}

# ── Figure B: line plot of enwrapment vs distance bin per treatment ───────────
pB <- ggplot(binned_b,
             aes(x = dist_bin, y = frac_enwrapped,
                 colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n_cells)) +
  scale_colour_manual(values = PALETTE) +
  scale_size_continuous(range = c(1, 4), name = "n cells") +
  labs(
    title    = "Part B: Enwrapment vs distance to nearest transduced cell (Core zone)",
    subtitle = paste("Sharp drop at 0-30px = local action.",
                     "Flat profile = diffuse extracellular action."),
    x        = "Distance to nearest mScarlet+ cell",
    y        = "Fraction enwrapped",
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig_partB_distance_enwrapment.pdf"),
       pB, width = 9, height = 5)
cat("\n  Saved fig_partB_distance_enwrapment.pdf\n\n")

# =============================================================================
# PART C — ADAMTS4 paradox: sampling artefact vs true null
# =============================================================================
cat("══════════════════════════════════════════════════════\n")
cat("Part C: ADAMTS4 paradox\n")
cat("══════════════════════════════════════════════════════\n\n")

# ── Test 1: are virus- cells in ADAMTS4 sections depleted vs mScarlet? ────────
# Use the virus- section-level LMM from Part A, focused on ADAMTS4 vs mScarlet

cat("Test 1: Virus- cells in ADAMTS4 vs mScarlet (from Part A model)\n")
ct_adamts4_neg <- ct_neg[str_detect(contrast, "ADAMTS4 -") &
                          !str_detect(contrast, "MD")]
cat(sprintf("  ADAMTS4 vs mScarlet (virus- cells): estimate=%.3f, p=%.4f, sig=%s\n\n",
            ct_adamts4_neg$estimate,
            ct_adamts4_neg$p.value,
            ct_adamts4_neg$sig))

# ── Test 2: within-ADAMTS4, does virus+ density predict local enwrapment? ──────
cat("Test 2: Within ADAMTS4 — does local virus+ density predict enwrapment?\n")
cat("  If high local transduction → lower enwrapment: supports biological effect.\n")
cat("  If high local transduction → higher enwrapment: supports sampling artefact.\n\n")

# Per section in ADAMTS4: n_virus+ cells, mean enwrapment of virus- cells
adamts4_sections <- tags[treatment == "ADAMTS4" & zone == "Core",
                          .(n_viruspos  = sum(virus_positive),
                            frac_virusneg_enwrapped = mean(enwrapped[!virus_positive]),
                            n_virusneg   = sum(!virus_positive)),
                          by = .(mouse_id, slice_id)]
adamts4_sections <- adamts4_sections[n_virusneg >= MIN_PV_SECTION]

cor_test_c <- cor.test(adamts4_sections$n_viruspos,
                        adamts4_sections$frac_virusneg_enwrapped,
                        method = "pearson")
cat(sprintf("  Correlation (n_virus+ vs virus- enwrapment in ADAMTS4 Core sections):\n"))
cat(sprintf("  r = %.3f, p = %.4f, n = %d\n",
            cor_test_c$estimate,
            cor_test_c$p.value,
            nrow(adamts4_sections)))
cat(sprintf("  Direction: %s\n\n",
            ifelse(cor_test_c$estimate < 0,
                   "more virus+ → LOWER enwrapment (biological effect)",
                   "more virus+ → HIGHER enwrapment (sampling artefact)")))

fwrite(adamts4_sections,
       file.path(OUT_DIR, "partC_adamts4_section_data.csv"))

# ── Test 3: compare virus+ enwrapment across distance bins for ADAMTS4 ─────────
cat("Test 3: Does ADAMTS4 virus+ enwrapment show any distance gradient?\n")

adamts4_binned <- core_tags[treatment == "ADAMTS4" & !is.na(dist_bin), .(
  frac_enwrapped_viruspos = mean(enwrapped[virus_positive],  na.rm = TRUE),
  frac_enwrapped_virusneg = mean(enwrapped[!virus_positive], na.rm = TRUE),
  n_viruspos = sum(virus_positive),
  n_virusneg = sum(!virus_positive)
), by = dist_bin]

adamts4_binned <- adamts4_binned[n_viruspos >= 5L & n_virusneg >= 5L]
adamts4_binned[, dist_bin := factor(dist_bin, levels = DIST_LABELS)]

cat("ADAMTS4 enwrapment by distance bin:\n")
print(adamts4_binned[order(dist_bin),
                      .(dist_bin,
                        virus_pos = round(frac_enwrapped_viruspos, 3),
                        virus_neg = round(frac_enwrapped_virusneg, 3),
                        n_viruspos, n_virusneg)])

fwrite(adamts4_binned,
       file.path(OUT_DIR, "partC_adamts4_distance_bins.csv"))

# ── Figure C: scatter + correlation for Test 2, bin plot for Test 3 ──────────
pC1 <- ggplot(adamts4_sections,
              aes(x = n_viruspos, y = frac_virusneg_enwrapped)) +
  geom_smooth(method = "lm", se = TRUE, colour = PALETTE["ADAMTS4"],
              alpha = 0.2, linewidth = 0.8) +
  geom_point(colour = PALETTE["ADAMTS4"], size = 3, alpha = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.3f\np = %.4f\nn = %d sections",
                           cor_test_c$estimate,
                           cor_test_c$p.value,
                           nrow(adamts4_sections)),
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(
    title    = "ADAMTS4 (Core): virus+ cell density vs virus- enwrapment",
    subtitle = "Negative r = more transduction → more depletion (biological); Positive = artefact",
    x        = "Number of virus+ cells per section",
    y        = "Mean enwrapment of virus- cells"
  ) +
  theme_bw(base_size = 11)

# Compare ADAMTS4 virus+ vs ADAMTS15 virus+ across distance bins
compare_binned <- core_tags[treatment %in% c("ADAMTS4","ADAMTS15","mScarlet",
                                               "C6ST1","C6ST1_ADAMTS15") &
                              virus_positive == TRUE &
                              !is.na(dist_bin), .(
  frac_enwrapped = mean(enwrapped),
  n_cells        = .N
), by = .(treatment, dist_bin)]
compare_binned <- compare_binned[n_cells >= MIN_PV_BIN]
compare_binned[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
compare_binned[, dist_bin  := factor(dist_bin,  levels = DIST_LABELS)]

pC2 <- ggplot(compare_binned,
              aes(x = dist_bin, y = frac_enwrapped,
                  colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = PALETTE) +
  labs(
    title    = "Virus+ cells only: enwrapment vs distance (Core zone)",
    subtitle = "Tests whether ADAMTS4 virus+ cells show any local depletion gradient",
    x        = "Distance to nearest mScarlet+ cell",
    y        = "Fraction enwrapped (virus+ cells only)",
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 30, hjust = 1),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

pC <- pC1 / pC2

ggsave(file.path(OUT_DIR, "fig_partC_adamts4_paradox.pdf"),
       pC, width = 9, height = 9)
cat("\n  Saved fig_partC_adamts4_paradox.pdf\n\n")

# =============================================================================
# COMBINED SUMMARY FIGURE
# =============================================================================
cat("Saving combined summary figure ...\n")

# Side-by-side: Part A (virus- contrasts) and Part B (distance curves)
pA_neg <- ggplot(ct_all_a[subset == "virus-"],
                 aes(x = estimate,
                     y = factor(str_remove(contrast," - mScarlet"),
                                levels = rev(TREATMENT_ORDER[-1])),
                     colour = factor(str_remove(contrast," - mScarlet"),
                                     levels = TREATMENT_ORDER))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = estimate - 1.96*SE, xmax = estimate + 1.96*SE),
                width = 0.3, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  labs(title    = "A: Virus- cells vs mScarlet (diffusion test)",
       subtitle = "Depletion of non-transduced cells = diffuse ECM action",
       x        = "Δ frac enwrapped vs mScarlet virus-",
       y        = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_combined <- pA_neg | pB
ggsave(file.path(OUT_DIR, "fig_summary_combined.pdf"),
       p_combined, width = 14, height = 5)
cat("  Saved fig_summary_combined.pdf\n\n")

# =============================================================================
# SAVE ALL CSVs
# =============================================================================
fwrite(ct_all_a,        file.path(OUT_DIR, "partA_contrasts_both_subsets.csv"))
fwrite(binned_b,        file.path(OUT_DIR, "partB_binned_enwrapment.csv"))
fwrite(adamts4_sections,file.path(OUT_DIR, "partC_adamts4_section_data.csv"))
fwrite(adamts4_binned,  file.path(OUT_DIR, "partC_adamts4_distance_bins.csv"))

cat("══════════════════════════════════════════════════════\n")
cat("Script 25 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  partA_contrasts_both_subsets.csv\n")
cat("  partA_diffuse_secretion_contrasts.csv\n")
cat("  partB_binned_enwrapment.csv\n")
cat("  partB_distance_glmm_coefficients.csv\n")
cat("  partC_adamts4_section_data.csv\n")
cat("  partC_adamts4_distance_bins.csv\n")
cat("  fig_partA_diffuse_secretion.pdf\n")
cat("  fig_partB_distance_enwrapment.pdf\n")
cat("  fig_partC_adamts4_paradox.pdf\n")
cat("  fig_summary_combined.pdf\n")
cat("══════════════════════════════════════════════════════\n")