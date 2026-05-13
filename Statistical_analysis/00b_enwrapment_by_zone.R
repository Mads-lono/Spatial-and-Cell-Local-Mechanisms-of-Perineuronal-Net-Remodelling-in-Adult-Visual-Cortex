# =============================================================================
# 00b_enwrapment_by_zone.R
# =============================================================================
# Extends script 00 by stratifying the nearest-neighbour enwrapment
# computation across injection zones (Core / Penumbra / Outside) rather
# than hemispheres.
#
# INPUT:  cells_with_zones.csv
# OUTPUT: results/00b_enwrapment_by_zone/
#   zone_enwrapment.csv           — per-animal × area × zone enwrapment fractions
#   zone_lmm_results.csv          — LMM: log(Core/Outside) ~ treatment + (1|mouse_id)
#   zone_fraction_results.csv     — LMM per zone: frac ~ treatment + (1|mouse_id)
#   fig_zone_enwrapment_ratio.pdf/png
#   fig_zone_enwrapment_profile.pdf/png
#   residual_check_00b.pdf/png


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

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV   <- "/path/to/cells_with_zones.csv"
RESULTS_DIR <- "/path/to/results"
OUT_DIR     <- file.path(RESULTS_DIR, "00b_enwrapment_by_zone")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENTS      <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS

# Animal to exclude (e.g. confirmed injection failure). Set to NULL to include all.
EXCLUDE_ANIMAL <- NULL   # e.g. "GROUP_ANIMAL_4"

COLOC_THRESH <- 30L
MIN_PV_ZONE  <- 10L   # lower than script 00 (20) — zones subdivide tissue
ZONE_ORDER   <- c("Core", "Penumbra", "Outside")

# Brain area filter
VIS_AREA_REGEX <- "isual"
IPSI_HEMI      <- "left"

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

cat("── Script 00b: Enwrapment by injection zone ─────────────────────────────────\n")
cat("   Parts: (0) compute  (A) Core/Outside ratio LMM\n")
cat("          (B) per-zone fraction LMMs  (C) zone × treatment interaction\n\n")

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
  else                  cat(sprintf("  OK        [%s]\n", label))
  fit
}

extract_contrasts <- function(fit, model_label = "") {
  coefs <- as.data.frame(coef(summary(fit, ddf = "Satterthwaite")))
  coefs$term <- rownames(coefs)
  setDT(coefs)
  res <- coefs[term != "(Intercept)"]
  setnames(res,
    c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)"),
    c("beta",     "se",         "df", "t",       "p"))
  res[, treatment := gsub("^treatment", "", term)]
  res[, ci_lo       := beta - qt(0.975, df) * se]
  res[, ci_hi       := beta + qt(0.975, df) * se]
  res[, fold_change := exp(beta)]
  res[, pct_change  := (fold_change - 1) * 100]
  res[, model       := model_label]
  res
}

add_fdr <- function(dt, p_col = "p", by_col = NULL, suffix = "_adj") {
  col_padj <- paste0("p_adj", suffix)
  col_sig  <- paste0("sig",   suffix)
  if (!is.null(by_col))
    dt[, (col_padj) := p.adjust(get(p_col), method = "BH"), by = by_col]
  else
    dt[, (col_padj) := p.adjust(get(p_col), method = "BH")]
  dt[, (col_sig) := fcase(
    get(col_padj) < 0.001, "***",
    get(col_padj) < 0.01,  "**",
    get(col_padj) < 0.05,  "*",
    get(col_padj) < 0.10,  ".",
    default = "ns"
  )]
  dt
}

# =============================================================================
# PART 0 — COMPUTE ZONE-LEVEL ENWRAPMENT
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART 0: Computing enwrapment per injection zone\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

LOAD_COLS <- c("mouse_id", "slice_id", "cell_type", "hemisphere",
               "x_hires", "y_hires", "brain_area", "zone", "treatment")

cat("Loading cells_with_zones.csv...\n")
t0 <- proc.time()
dt <- fread(CELLS_CSV, select = LOAD_COLS)
cat(sprintf("  Loaded %s rows in %.1fs\n",
    format(nrow(dt), big.mark = ","), (proc.time() - t0)[["elapsed"]]))

dt <- dt[str_detect(brain_area, regex(VIS_AREA_REGEX, ignore_case = TRUE))]
dt <- dt[treatment %in% TREATMENTS]
if (!is.null(EXCLUDE_ANIMAL))
  dt <- dt[mouse_id != EXCLUDE_ANIMAL]
dt <- dt[zone %in% ZONE_ORDER]
dt <- dt[hemisphere == IPSI_HEMI]
cat(sprintf("  After filters: %s rows\n", format(nrow(dt), big.mark = ",")))

dt[, area := str_trim(str_split_fixed(brain_area, ", ", 2)[, 1])]

pv  <- dt[cell_type == "PV"]
pnn <- dt[cell_type == "PNN"]
cat(sprintf("  PV: %s  |  PNN: %s\n\n",
    format(nrow(pv),  big.mark = ","),
    format(nrow(pnn), big.mark = ",")))

zone_keys    <- unique(pv[, .(treatment, mouse_id, area, zone)])
cat(sprintf("  Processing %d animal × area × zone strata...\n", nrow(zone_keys)))

t1 <- proc.time()
zone_records <- vector("list", nrow(zone_keys))

for (i in seq_len(nrow(zone_keys))) {
  key   <- zone_keys[i]
  pv_k  <- pv[ mouse_id == key$mouse_id & area == key$area & zone == key$zone,
                .(x_hires, y_hires)]
  pnn_k <- pnn[mouse_id == key$mouse_id & area == key$area & zone == key$zone,
                .(x_hires, y_hires)]

  n_pv <- nrow(pv_k)
  if (n_pv < MIN_PV_ZONE || nrow(pnn_k) == 0L) next

  frac <- compute_frac_enwrapped(pv_k, pnn_k)
  if (is.na(frac)) next

  zone_records[[i]] <- list(
    treatment      = key$treatment,
    mouse_id       = key$mouse_id,
    area           = key$area,
    zone           = key$zone,
    n_pv           = n_pv,
    n_pnn          = nrow(pnn_k),
    frac_enwrapped = frac
  )

  if (i %% 500L == 0L)
    cat(sprintf("  %d / %d (%.0fs)\n", i, nrow(zone_keys),
        (proc.time() - t1)[["elapsed"]]))
}

zone_df <- rbindlist(Filter(Negate(is.null), zone_records))
cat(sprintf("\n  Valid zone records: %d across %d areas\n",
    nrow(zone_df), uniqueN(zone_df$area)))

zone_df[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
zone_df[, zone      := factor(zone,      levels = ZONE_ORDER)]

fwrite(zone_df, file.path(OUT_DIR, "zone_enwrapment.csv"))
cat("  Saved: zone_enwrapment.csv\n\n")
rm(dt, pv, pnn, zone_records); gc()

# Animal-level means for descriptive summary and profile plot
animal_zone <- zone_df[, .(
  frac_enwrapped = mean(frac_enwrapped, na.rm = TRUE),
  n_areas        = .N
), by = .(treatment, mouse_id, zone)]

cat("── Animal-level zone means (descriptive) ────────────────────────────────────\n")
print(dcast(animal_zone, treatment ~ zone, value.var = "frac_enwrapped",
            fun.aggregate = mean)[order(match(treatment, TREATMENT_ORDER))])
cat("\n")

# =============================================================================
# PART A — CORE / OUTSIDE RATIO LMM
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART A: Core/Outside ratio LMM\n")
cat("        log(frac_core/frac_outside) ~ treatment + (1|mouse_id)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

core_out <- merge(
  zone_df[zone == "Core",    .(treatment, mouse_id, area, frac_core    = frac_enwrapped)],
  zone_df[zone == "Outside", .(treatment, mouse_id, area, frac_outside = frac_enwrapped)],
  by = c("treatment", "mouse_id", "area")
)
core_out <- core_out[frac_outside > 0]
core_out[, ratio     := frac_core / frac_outside]
core_out[, log_ratio := log(ratio)]
core_out <- core_out[is.finite(log_ratio)]
core_out[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("  Core/Outside pairs: %d area-level rows, %d animals\n\n",
    nrow(core_out), uniqueN(core_out$mouse_id)))

desc_ratio <- core_out[, .(
  mean_ratio = mean(exp(log_ratio)),
  sd         = sd(exp(log_ratio)),
  n_areas    = .N
), by = treatment]
cat("── Mean Core/Outside enwrapment ratio by treatment ─────────────────────────\n")
print(desc_ratio[order(match(treatment, TREATMENT_ORDER))])
cat("\n")

fit_ratio <- fit_lmm(
  log_ratio ~ treatment + (1 | mouse_id),
  core_out, "Core/Outside ratio"
)

res_ratio <- extract_contrasts(fit_ratio, model_label = "Core/Outside")
res_ratio <- add_fdr(res_ratio, suffix = "_adj")

cat("\n── Part A results ───────────────────────────────────────────────────────────\n")
print(res_ratio[, .(treatment, beta, se, ci_lo, ci_hi, pct_change, p, p_adj_adj, sig_adj)])

# =============================================================================
# PART B — PER-ZONE FRACTION LMMs
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART B: Per-zone enwrapment fraction LMMs\n")
cat("        frac_enwrapped ~ treatment + (1|mouse_id), per zone\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

zone_res <- rbindlist(lapply(ZONE_ORDER, function(z) {
  d <- zone_df[zone == z]
  if (uniqueN(d$mouse_id) < 3 || nrow(d) < 10) {
    cat(sprintf("  SKIP [%s]: too few observations\n", z))
    return(NULL)
  }
  cat(sprintf("── Zone: %s (%d area-level rows, %d animals) ──\n",
      z, nrow(d), uniqueN(d$mouse_id)))
  fit <- fit_lmm(frac_enwrapped ~ treatment + (1 | mouse_id), d, z)
  res <- extract_contrasts(fit, model_label = z)
  res[, zone := z]
  res
}), fill = TRUE)

if (nrow(zone_res) > 0) {
  zone_res[, zone := factor(zone, levels = ZONE_ORDER)]
  zone_res <- add_fdr(zone_res, suffix = "_adj", by_col = "zone")

  cat("\n── Part B results (significant only) ───────────────────────────────────────\n")
  sig_b <- zone_res[sig_adj != "ns"]
  if (nrow(sig_b) > 0)
    print(sig_b[, .(zone, treatment, beta, pct_change, p, p_adj_adj, sig_adj)])
  else
    cat("  No significant hits after BH-FDR correction per zone.\n")
}

# =============================================================================
# PART C — ZONE × TREATMENT INTERACTION
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("PART C: Zone × treatment interaction\n")
cat("        frac_enwrapped ~ treatment * zone + (1|mouse_id)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

zone_df[, zone := factor(zone, levels = ZONE_ORDER)]

fit_ix <- fit_lmm(
  frac_enwrapped ~ treatment * zone + (1 | mouse_id),
  zone_df, "treatment*zone interaction"
)

anova_ix   <- anova(fit_ix, ddf = "Satterthwaite")
ix_p_treat <- anova_ix["treatment",      "Pr(>F)"]
ix_p_zone  <- anova_ix["zone",           "Pr(>F)"]
ix_p_inter <- anova_ix["treatment:zone", "Pr(>F)"]

cat(sprintf("\n  treatment F p        = %.4f\n",  ix_p_treat))
cat(sprintf("  zone F p             = %.4f\n",  ix_p_zone))
cat(sprintf("  treatment:zone F p   = %.4f  %s\n", ix_p_inter,
    ifelse(ix_p_inter < 0.05, "✓ Zone-specific effect", "NS")))

emm_ix <- emmeans(fit_ix, ~ treatment | zone)
con_ix  <- contrast(emm_ix, method = "trt.vs.ctrl", ref = "mScarlet")
res_ix  <- as.data.table(summary(con_ix, infer = TRUE))
setnames(res_ix,
  c("estimate","SE","lower.CL","upper.CL","t.ratio","p.value"),
  c("beta",    "se","ci_lo",   "ci_hi",   "t",      "p"))
res_ix[, treatment := gsub(" - mScarlet", "", contrast)]
res_ix[, zone      := factor(zone, levels = ZONE_ORDER)]
res_ix <- add_fdr(res_ix, suffix = "_adj", by_col = "zone")
res_ix[, pct_change := (exp(beta) - 1) * 100]

cat("\n── Part C: treatment × zone contrasts (vs reference) ────────────────────────\n")
print(res_ix[, .(zone, treatment, beta, pct_change, p, p_adj_adj, sig_adj)][
  order(zone, match(treatment, TREATMENT_ORDER[-1]))])

# =============================================================================
# SAVE RESULTS
# =============================================================================
all_results <- rbindlist(list(
  if (exists("res_ratio")) res_ratio[, zone := "Core/Outside (ratio)"] else NULL,
  if (nrow(zone_res) > 0) zone_res else NULL
), fill = TRUE)

fwrite(all_results, file.path(OUT_DIR, "zone_lmm_results.csv"))
fwrite(res_ix,      file.path(OUT_DIR, "zone_interaction_results.csv"))
cat(sprintf("\n  Saved: zone_lmm_results.csv, zone_interaction_results.csv\n"))

# =============================================================================
# FIGURE 1 — Forest plot: Core/Outside enwrapment ratio
# =============================================================================
plot_ratio <- copy(res_ratio)
plot_ratio[, treatment := factor(treatment, levels = rev(TREATMENT_ORDER[-1]))]

p1 <- ggplot(plot_ratio, aes(x = pct_change, y = treatment, colour = treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = (exp(ci_lo)-1)*100, xmax = (exp(ci_hi)-1)*100),
                width = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(x = (exp(ci_hi)-1)*100 + 2, label = sig_adj),
            hjust = 0, size = 4, show.legend = FALSE) +
  scale_colour_manual(values = PALETTE[-1], guide = "none") +
  scale_x_continuous(name = "% change in Core/Outside enwrapment ratio vs reference",
                     expand = expansion(mult = c(0.05, 0.20))) +
  scale_y_discrete(name = NULL) +
  labs(
    title    = "Core/Outside PV enwrapment ratio",
    subtitle = "Negative = reduced enwrapment at injection core relative to distal tissue.\nBH-FDR corrected."
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

# =============================================================================
# FIGURE 2 — Zone profile: mean enwrapment fraction per zone × treatment
# =============================================================================
profile_means <- animal_zone[, .(
  mean_frac = mean(frac_enwrapped, na.rm = TRUE),
  se_frac   = sd(frac_enwrapped, na.rm = TRUE) / sqrt(.N)
), by = .(treatment, zone)]
profile_means[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
profile_means[, zone      := factor(zone,      levels = ZONE_ORDER)]

p2 <- ggplot(profile_means,
             aes(x = zone, y = mean_frac, colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_frac - se_frac, ymax = mean_frac + se_frac),
                width = 0.15, linewidth = 0.6) +
  scale_colour_manual(values = PALETTE, name = NULL) +
  scale_x_discrete(name = "Injection zone") +
  scale_y_continuous(name = "Mean fraction of PV cells enwrapped",
                     labels = percent_format(accuracy = 1)) +
  labs(
    title    = "PV enwrapment fraction across injection zones",
    subtitle = "Mean ± SEM across animals."
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right",
        panel.grid.minor = element_blank())

# =============================================================================
# FIGURE 3 — Heatmap of Part C contrasts
# =============================================================================
if (nrow(res_ix) > 0) {
  res_ix[, treatment_f := factor(treatment, levels = TREATMENT_ORDER[-1])]
  p3 <- ggplot(res_ix, aes(x = zone, y = treatment_f, fill = pct_change)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig_adj != "ns", sig_adj, "")),
              size = 4, fontface = "bold") +
    scale_fill_gradient2(
      low = "#d73027", mid = "white", high = "#4575b4",
      midpoint = 0, name = "% change\nvs reference"
    ) +
    scale_x_discrete(name = "Injection zone") +
    scale_y_discrete(name = NULL) +
    labs(
      title    = "Zone × treatment interaction",
      subtitle = sprintf(
        "treatment:zone F p = %.4f  %s\nStars = BH-FDR < 0.05 per zone.",
        ix_p_inter,
        ifelse(ix_p_inter < 0.05, "✓ Zone-specific effect", "NS"))
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())
}

# =============================================================================
# RESIDUAL CHECK
# =============================================================================
resid_data <- data.frame(
  fitted   = fitted(fit_ratio),
  residual = residuals(fit_ratio)
)
p_resid <- ggplot(resid_data, aes(x = fitted, y = residual)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, colour = "red", linewidth = 0.7) +
  labs(title = "Residual check — ratio LMM (Part A)",
       x = "Fitted values", y = "Residuals") +
  theme_bw(base_size = 10)

# =============================================================================
# SAVE FIGURES
# =============================================================================
for (ext in c("pdf", "png")) {
  dpi <- ifelse(ext == "pdf", 300, 150)
  ggsave(file.path(OUT_DIR, paste0("fig_zone_enwrapment_ratio.",   ext)),
         p1, width = 9, height = 5, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("fig_zone_enwrapment_profile.", ext)),
         p2, width = 9, height = 5, dpi = dpi)
  if (nrow(res_ix) > 0)
    ggsave(file.path(OUT_DIR, paste0("fig_zone_enwrapment_heatmap.", ext)),
           p3, width = 7, height = 5, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("residual_check_00b.", ext)),
         p_resid, width = 6, height = 4, dpi = dpi)
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("SCRIPT 00b SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Part A — Core/Outside ratio hits (FDR): %d\n",
    nrow(res_ratio[sig_adj != "ns"])))
cat(sprintf("  Part B — Per-zone hits (FDR, any zone):  %d\n",
    if (nrow(zone_res) > 0) nrow(zone_res[sig_adj != "ns"]) else 0))
cat(sprintf("  Part C — treatment:zone interaction F p = %.4f  %s\n",
    ix_p_inter, ifelse(ix_p_inter < 0.05, "✓ Zone-specific", "NS")))
cat(sprintf("\n── Outputs → %s\n", OUT_DIR))
cat("   zone_enwrapment.csv           zone_lmm_results.csv\n")
cat("   zone_interaction_results.csv  fig_zone_* (pdf + png)\n")
