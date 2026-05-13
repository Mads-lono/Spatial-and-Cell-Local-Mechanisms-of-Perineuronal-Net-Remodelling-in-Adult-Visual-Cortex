# =============================================================================
# Script 17: PV cell intensity LMM
# Question: Does PV fluorescence intensity change with treatment?
# Input:  cells_with_zones.csv (cell_type == "PV", normalized_btm20)
# Model:  section_mean_pv ~ treatment + (1|mouse_id)
# Also:   zone-stratified version (Core / Penumbra / Outside separately)
# Output: results/17_pv_intensity/
# =============================================================================

library(data.table)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Paths & constants
# -----------------------------------------------------------------------------
CELLS_PATH  <- "/path/to/analysis_results/cells_with_zones.csv"
OUT_DIR     <- "/path/to/results/17_pv_intensity"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TREATMENT_ORDER <- c("mScarlet", "ADAMTS4", "ADAMTS4_MD",
                     "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15")

PALETTE <- c(
  mScarlet       = "#c0392b",
  ADAMTS4        = "#2980b9",
  ADAMTS4_MD     = "#1abc9c",
  ADAMTS15       = "#8e44ad",
  C6ST1          = "#e67e22",
  C6ST1_ADAMTS15 = "#27ae60"
)

MIN_PV_SECTION <- 5L   # minimum PV cells per section to include

# -----------------------------------------------------------------------------
# 1. Load & apply data corrections
# -----------------------------------------------------------------------------
cat("Loading cells_with_zones.csv ...\n")
cells <- fread(CELLS_PATH)

# Keep PV cells only
cells <- cells[cell_type == "PV"]
cat(sprintf("PV cells loaded: %d\n", nrow(cells)))

# Fix naming error
cells[treatment == "C6ST1_ADAMTS4", treatment := "C6ST1_ADAMTS15"]


# Exclude injection failure
cells <- cells[mouse_id != ""]  #replace with actual mouse_id if you need to exclude a specific animal(e.g., injection failure case)

# Factor with reference level = mScarlet
cells[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat(sprintf("After corrections: %d cells, %d animals\n",
            nrow(cells), uniqueN(cells$mouse_id)))

# -----------------------------------------------------------------------------
# 2. Aggregate to section level (avoid pseudoreplication)
# Per section: mean PV intensity, cell count, zone
# -----------------------------------------------------------------------------
# section = mouse_id x slice_id x zone x hemisphere
sec <- cells[, .(
  mean_pv       = mean(normalized_btm20, na.rm = TRUE),
  n_pv          = .N,
  treatment     = treatment[1]
), by = .(mouse_id, slice_id, zone, hemisphere)]

# Apply minimum cell count filter
sec <- sec[n_pv >= MIN_PV_SECTION]
cat(sprintf("Sections after min-cell filter (%d): %d\n", MIN_PV_SECTION, nrow(sec)))

# -----------------------------------------------------------------------------
# 3. Primary LMM — all zones pooled
# -----------------------------------------------------------------------------
cat("\n── Primary LMM: mean_pv ~ treatment + (1|mouse_id) ──\n")

m_primary <- lmer(mean_pv ~ treatment + (1|mouse_id), data = sec, REML = FALSE)
print(summary(m_primary))

# Pairwise contrasts vs mScarlet
em_primary <- emmeans(m_primary, ~ treatment)
contrasts_primary <- as.data.table(
  contrast(em_primary, method = "trt.vs.ctrl", ref = "mScarlet", adjust = "BH")
)
contrasts_primary[, sig := ifelse(p.value < 0.05, "*", "")]
cat("\nPrimary contrasts (BH-adjusted):\n")
print(contrasts_primary[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

fwrite(contrasts_primary,
       file.path(OUT_DIR, "primary_contrasts.csv"))

# -----------------------------------------------------------------------------
# 4. Zone-stratified LMMs
# -----------------------------------------------------------------------------
zones <- c("Core", "Penumbra", "Outside")
zone_results <- list()

for (z in zones) {
  d <- sec[zone == z]
  if (nrow(d) < 10 || uniqueN(d$mouse_id) < 4) {
    cat(sprintf("\nZone %s: insufficient data, skipping\n", z))
    next
  }
  cat(sprintf("\n── Zone LMM: %s ──\n", z))
  m <- tryCatch(
    lmer(mean_pv ~ treatment + (1|mouse_id), data = d, REML = FALSE),
    error = function(e) { cat("  Model failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) next

  em <- emmeans(m, ~ treatment)
  ct <- as.data.table(
    contrast(em, method = "trt.vs.ctrl", ref = "mScarlet", adjust = "BH")
  )
  ct[, zone := z]
  ct[, sig  := ifelse(p.value < 0.05, "*", "")]
  cat(sprintf("  Sections: %d | Animals: %d\n", nrow(d), uniqueN(d$mouse_id)))
  print(ct[, .(contrast, zone, estimate, SE, t.ratio, p.value, sig)])
  zone_results[[z]] <- ct
}

zone_dt <- rbindlist(zone_results)
fwrite(zone_dt, file.path(OUT_DIR, "zone_contrasts.csv"))

# -----------------------------------------------------------------------------
# 5. Animal-level summary (for plotting & reporting)
# -----------------------------------------------------------------------------
animal_summary <- sec[, .(
  mean_pv  = mean(mean_pv, na.rm = TRUE),
  sd_pv    = sd(mean_pv, na.rm = TRUE),
  n_sec    = .N,
  treatment = treatment[1]
), by = mouse_id]

fwrite(animal_summary, file.path(OUT_DIR, "animal_summary.csv"))

# Also zone-stratified animal summary
animal_zone_summary <- sec[, .(
  mean_pv  = mean(mean_pv, na.rm = TRUE),
  n_sec    = .N,
  treatment = treatment[1]
), by = .(mouse_id, zone)]

fwrite(animal_zone_summary, file.path(OUT_DIR, "animal_zone_summary.csv"))

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------

# 6a: Boxplot — all zones pooled
animal_summary[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

p1 <- ggplot(animal_summary, aes(x = treatment, y = mean_pv, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = PALETTE) +
  labs(
    title    = "PV cell fluorescence intensity by treatment",
    subtitle = "Section-level means, all zones pooled",
    x        = NULL,
    y        = "Mean PV intensity (normalized_btm20)",
    caption  = "LMM: mean_pv ~ treatment + (1|mouse_id); BH-FDR"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "fig_pv_intensity_primary.pdf"),
       p1, width = 7, height = 5)

# 6b: Faceted boxplot by zone
if (nrow(zone_dt) > 0) {
  animal_zone_summary[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
  animal_zone_summary[, zone := factor(zone, levels = zones)]

  p2 <- ggplot(animal_zone_summary,
               aes(x = treatment, y = mean_pv, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.9) +
    scale_fill_manual(values = PALETTE) +
    facet_wrap(~ zone) +
    labs(
      title    = "PV fluorescence intensity by treatment and zone",
      x        = NULL,
      y        = "Mean PV intensity (normalized_btm20)",
      caption  = "Zone-stratified LMMs; BH-FDR within each zone"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(OUT_DIR, "fig_pv_intensity_by_zone.pdf"),
         p2, width = 10, height = 5)
}

# -----------------------------------------------------------------------------
# 7. Summary printout
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 17 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  primary_contrasts.csv\n")
cat("  zone_contrasts.csv\n")
cat("  animal_summary.csv\n")
cat("  animal_zone_summary.csv\n")
cat("  fig_pv_intensity_primary.pdf\n")
cat("  fig_pv_intensity_by_zone.pdf\n")
cat("══════════════════════════════════════════════════════\n")
