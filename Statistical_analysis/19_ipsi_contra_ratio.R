# =============================================================================
# Script 19: Ipsilateral / contralateral enwrapment ratio
# Question: Does treatment reduce ipsilateral enwrapment relative to the
#           within-animal contralateral internal control?
# Input:  slice_enwrapment.csv (frac_ipsi, frac_contra per slice per animal)
# Model:  ratio ~ treatment  (simple lm, one value per animal)
# Output: results/19_ipsi_contra_ratio/
# =============================================================================

library(data.table)
library(emmeans)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Paths & constants
# -----------------------------------------------------------------------------
SLICE_PATH <- "/path/to/results/00_enwrapment/slice_enwrapment.csv"
OUT_DIR    <- "/path/to/results/19_ipsi_contra_ratio"
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

MIN_PV_SLICE <- 5L   # minimum PV cells per hemisphere per slice

# -----------------------------------------------------------------------------
# 1. Load & apply data corrections
# -----------------------------------------------------------------------------
cat("Loading slice_enwrapment.csv ...\n")
dat <- fread(SLICE_PATH)

# Fix naming error (old label in this file)
dat[treatment == "C6ST1_ADAMTS4", treatment := "C6ST1_ADAMTS15"]


# Exclude injection failure — use ORIGINAL mouse_id label as it appears

dat <- dat[mouse_id != ""]   # Exclude slice with missing mouse_id (injection failure)

cat(sprintf("Slices loaded: %d across %d animals\n",
            nrow(dat), uniqueN(dat$mouse_id)))

# -----------------------------------------------------------------------------
# 2. Filter slices with insufficient cells, then aggregate to animal level
# Weighted mean: weight each slice by n_pv_ipsi / n_pv_contra
# -----------------------------------------------------------------------------
dat <- dat[n_pv_ipsi >= MIN_PV_SLICE & n_pv_contra >= MIN_PV_SLICE]
cat(sprintf("Slices after min-cell filter: %d\n", nrow(dat)))

animal_dat <- dat[, .(
  mean_ipsi   = weighted.mean(frac_ipsi,   w = n_pv_ipsi,   na.rm = TRUE),
  mean_contra = weighted.mean(frac_contra, w = n_pv_contra,  na.rm = TRUE),
  n_slices    = .N,
  treatment   = treatment[1]
), by = mouse_id]

# Compute ratio and log ratio
animal_dat[, ratio     := mean_ipsi / mean_contra]
animal_dat[, log_ratio := log(ratio)]
animal_dat[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

cat("\nAnimal-level summary:\n")
print(animal_dat[order(treatment),
                 .(mouse_id, treatment, mean_ipsi, mean_contra,
                   ratio     = round(ratio, 3),
                   log_ratio = round(log_ratio, 3),
                   n_slices)])

fwrite(animal_dat, file.path(OUT_DIR, "animal_ratios.csv"))

# -----------------------------------------------------------------------------
# 3. Primary model: log_ratio ~ treatment
# log(ipsi/contra) = 0 means no lateralisation
# mScarlet intercept captures baseline injection artefact
# -----------------------------------------------------------------------------
cat("\n── Primary model: log_ratio ~ treatment ──\n")
m_log <- lm(log_ratio ~ treatment, data = animal_dat)
print(summary(m_log))

em_log <- emmeans(m_log, ~ treatment)
ct_log <- as.data.table(
  contrast(em_log, method = "trt.vs.ctrl", ref = "mScarlet", adjust = "BH")
)
ct_log[, sig := ifelse(p.value < 0.05, "*", "")]

cat("\nContrasts vs mScarlet (log ratio, BH-adjusted):\n")
print(ct_log[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

fwrite(ct_log, file.path(OUT_DIR, "contrasts_log_ratio.csv"))

# -----------------------------------------------------------------------------
# 4. One-sample t-tests: is each group's log_ratio != 0?
# Tests whether each treatment shows a lateralised effect in its own right
# -----------------------------------------------------------------------------
cat("\nOne-sample t-tests: log_ratio != 0 per treatment?\n")
one_sample <- animal_dat[, {
  tt <- t.test(log_ratio, mu = 0)
  .(mean_log_ratio = mean(log_ratio),
    t_stat         = tt$statistic,
    p_one_sample   = tt$p.value,
    n              = .N)
}, by = treatment]
one_sample[, p_adj := p.adjust(p_one_sample, method = "BH")]
one_sample[, sig   := ifelse(p_adj < 0.05, "*", "")]
print(one_sample[order(treatment)])
fwrite(one_sample, file.path(OUT_DIR, "one_sample_tests.csv"))

# -----------------------------------------------------------------------------
# 5. Secondary model on raw ratio (for interpretability in reporting)
# -----------------------------------------------------------------------------
cat("\n── Secondary model: ratio ~ treatment ──\n")
m_raw <- lm(ratio ~ treatment, data = animal_dat)
print(summary(m_raw))

em_raw <- emmeans(m_raw, ~ treatment)
ct_raw <- as.data.table(
  contrast(em_raw, method = "trt.vs.ctrl", ref = "mScarlet", adjust = "BH")
)
ct_raw[, sig := ifelse(p.value < 0.05, "*", "")]

cat("\nContrasts vs mScarlet (raw ratio, BH-adjusted):\n")
print(ct_raw[, .(contrast, estimate, SE, df, t.ratio, p.value, sig)])

fwrite(ct_raw, file.path(OUT_DIR, "contrasts_raw_ratio.csv"))

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------

# 6a: Dotplot of individual animal ratios with group means ± SE
group_means <- animal_dat[, .(
  mean_ratio = mean(ratio),
  se_ratio   = sd(ratio) / sqrt(.N)
), by = treatment]

p1 <- ggplot(animal_dat, aes(x = treatment, y = ratio, colour = treatment)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.1, size = 3.5, alpha = 0.9) +
  geom_pointrange(data = group_means,
                  aes(x = treatment, y = mean_ratio,
                      ymin = mean_ratio - se_ratio,
                      ymax = mean_ratio + se_ratio),
                  colour = "black", size = 0.6, linewidth = 1) +
  scale_colour_manual(values = PALETTE) +
  labs(
    title    = "Ipsilateral / contralateral enwrapment ratio",
    subtitle = "Ratio < 1 = ipsilateral reduction; dashed line = no lateralisation",
    x        = NULL,
    y        = "Ipsi / contra enwrapment ratio",
    caption  = "Points = individual animals; crossbar = mean ± SE\nlm: ratio ~ treatment, BH-FDR"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "fig_ratio_dotplot.pdf"),
       p1, width = 7, height = 5)

# 6b: Log ratio version
p2 <- ggplot(animal_dat, aes(x = treatment, y = log_ratio, colour = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.1, size = 3.5, alpha = 0.9) +
  scale_colour_manual(values = PALETTE) +
  labs(
    title    = "Log ipsi/contra enwrapment ratio by treatment",
    subtitle = "Log ratio < 0 = ipsilateral reduction",
    x        = NULL,
    y        = "log(ipsi / contra enwrapment)",
    caption  = "lm: log_ratio ~ treatment, BH-FDR"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "fig_log_ratio_dotplot.pdf"),
       p2, width = 7, height = 5)

# -----------------------------------------------------------------------------
# 7. Summary printout
# -----------------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════\n")
cat("Script 19 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Files written:\n")
cat("  animal_ratios.csv\n")
cat("  contrasts_log_ratio.csv\n")
cat("  contrasts_raw_ratio.csv\n")
cat("  one_sample_tests.csv\n")
cat("  fig_ratio_dotplot.pdf\n")
cat("  fig_log_ratio_dotplot.pdf\n")
cat("══════════════════════════════════════════════════════\n")