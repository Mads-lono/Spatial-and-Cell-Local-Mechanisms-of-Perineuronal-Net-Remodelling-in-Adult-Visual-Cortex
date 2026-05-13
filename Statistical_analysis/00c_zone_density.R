# =============================================================================
# 00c_zone_density.R
# =============================================================================
# Complements scripts 00 and 00b with a structurally independent outcome:
# cell count per mm². Density captures whether the intervention reduced
# the NUMBER of detectable PNN/PV structures, independently of enwrapment.
#
# INPUT:  results/zone_density/zone_density.csv   — from compute_zone_density.ipynb
# OUTPUT: results/00c_zone_density/
#   zone_density_lmm_results.csv      — Part A: per-zone LMMs (WFA + PV)
#   zone_density_ratio_results.csv    — Part B: Core/Outside ratio LMM
#   zone_density_interaction.csv      — Part C: zone × treatment interaction
#   fig_density_profile.pdf/png
#   fig_density_ratio.pdf/png
#   fig_density_heatmap.pdf/png
#   residual_check_00c.pdf/png


suppressPackageStartupMessages({
  library(data.table)
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
DENSITY_CSV <- "/path/to/results/zone_density/zone_density.csv"
RESULTS_DIR <- "/path/to/results"
OUT_DIR     <- file.path(RESULTS_DIR, "00c_zone_density")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENTS      <- c("mScarlet","ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1","C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS

# Animal to exclude (e.g. confirmed injection failure). Set to NULL to include all.
EXCLUDE_ANIMAL <- NULL   # e.g. "GROUP_ANIMAL_4"

ZONE_ORDER  <- c("Core","Penumbra","Outside")
STAIN_ORDER <- c("WFA","PV")

PALETTE <- c(
  mScarlet       = "#888888",
  ADAMTS4        = "#4e9af1",
  ADAMTS4_MD     = "#f17c4e",
  ADAMTS15       = "#4ef196",
  C6ST1          = "#c44ef1",
  C6ST1_ADAMTS15 = "#f1c44e"
)

cat("── Script 00c: Zone-level cell density ──────────────────────────────────────\n")
cat("   Parts: (A) per-zone LMMs  (B) Core/Outside ratio\n")
cat("          (C) zone × treatment interaction\n\n")

# ── Helpers ───────────────────────────────────────────────────────────────────
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
    c("Estimate","Std. Error","df","t value","Pr(>|t|)"),
    c("beta",    "se",        "df","t",      "p"))
  res[, treatment  := gsub("^treatment", "", term)]
  res[, ci_lo      := beta - qt(0.975, df) * se]
  res[, ci_hi      := beta + qt(0.975, df) * se]
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
# LOAD & PREPARE
# =============================================================================
cat("Loading zone_density.csv...\n")
dt <- fread(DENSITY_CSV)
setnames(dt, "animal_id", "mouse_id", skip_absent = TRUE)
if (!is.null(EXCLUDE_ANIMAL))
  dt <- dt[mouse_id != EXCLUDE_ANIMAL]
dt <- dt[treatment  %in% TREATMENTS]
dt <- dt[zone       %in% ZONE_ORDER]
dt <- dt[is.finite(density_mm2) & density_mm2 > 0]

dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
dt[, zone      := factor(zone,      levels = ZONE_ORDER)]
dt[, staining  := factor(staining,  levels = STAIN_ORDER)]
dt[, log_dens  := log(density_mm2)]

cat(sprintf("  Rows: %s  |  Animals: %d  |  Areas: %d\n\n",
    format(nrow(dt), big.mark=","), uniqueN(dt$mouse_id), uniqueN(dt$brain_area)))

cat("── Mean density (cells/mm²) by treatment × zone × staining ─────────────────\n")
desc <- dt[, .(mean_dens = mean(density_mm2), sd = sd(density_mm2), n = .N),
           by = .(staining, treatment, zone)]
desc_cast <- dcast(desc, staining + treatment ~ zone, value.var = "mean_dens",
                   fun.aggregate = mean)
desc_cast <- desc_cast[order(staining, match(treatment, TREATMENT_ORDER))]
num_cols <- names(desc_cast)[sapply(desc_cast, is.numeric)]
desc_cast[, (num_cols) := lapply(.SD, round, 1), .SDcols = num_cols]
print(desc_cast)

# =============================================================================
# PART A — PER-ZONE DENSITY LMMs
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART A: Per-zone density LMMs\n")
cat("        log(density_mm2) ~ treatment + (1|mouse_id), per zone × staining\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

zone_res <- rbindlist(lapply(STAIN_ORDER, function(st) {
  rbindlist(lapply(ZONE_ORDER, function(z) {
    d <- dt[staining == st & zone == z]
    if (uniqueN(d$mouse_id) < 3 || nrow(d) < 10) {
      cat(sprintf("  SKIP [%s × %s]: too few rows\n", st, z))
      return(NULL)
    }
    cat(sprintf("── %s × %s (%d rows, %d animals) ──\n",
        st, z, nrow(d), uniqueN(d$mouse_id)))
    fit <- fit_lmm(log_dens ~ treatment + (1 | mouse_id), d,
                   paste(st, z, sep=" × "))
    res <- extract_contrasts(fit, model_label = paste(st, z, sep=" × "))
    res[, staining := st]
    res[, zone     := z]
    res
  }), fill = TRUE)
}), fill = TRUE)

zone_res[, zone     := factor(zone,     levels = ZONE_ORDER)]
zone_res[, staining := factor(staining, levels = STAIN_ORDER)]
zone_res <- add_fdr(zone_res, suffix = "_adj",
                    by_col = c("staining","zone"))

cat("\n── Part A: significant hits (FDR < 0.10) ───────────────────────────────────\n")
sig_a <- zone_res[sig_adj %in% c("***","**","*",".")]
if (nrow(sig_a) > 0) {
  print(sig_a[, .(staining, zone, treatment, beta, pct_change, p, p_adj_adj, sig_adj)])
} else {
  cat("  No hits reaching p_adj < 0.10 in any zone × staining.\n")
}

fwrite(zone_res, file.path(OUT_DIR, "zone_density_lmm_results.csv"))
cat(sprintf("  Saved: zone_density_lmm_results.csv\n\n"))

# =============================================================================
# PART B — CORE / OUTSIDE DENSITY RATIO LMM
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART B: Core/Outside density ratio LMM\n")
cat("        log(density_core/density_outside) ~ treatment + (1|mouse_id)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

ratio_res <- rbindlist(lapply(STAIN_ORDER, function(st) {
  core_out <- merge(
    dt[staining == st & zone == "Core",
       .(treatment, mouse_id, brain_area, dens_core    = density_mm2)],
    dt[staining == st & zone == "Outside",
       .(treatment, mouse_id, brain_area, dens_outside = density_mm2)],
    by = c("treatment","mouse_id","brain_area")
  )
  core_out <- core_out[dens_outside > 0]
  core_out[, log_ratio := log(dens_core / dens_outside)]
  core_out <- core_out[is.finite(log_ratio)]
  core_out[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

  cat(sprintf("── %s: %d area-level pairs, %d animals ──\n",
      st, nrow(core_out), uniqueN(core_out$mouse_id)))

  desc_r <- core_out[, .(mean_ratio = mean(exp(log_ratio)),
                          sd = sd(exp(log_ratio)), n = .N),
                     by = treatment]
  print(desc_r[order(match(treatment, TREATMENT_ORDER))])
  cat("\n")

  if (uniqueN(core_out$mouse_id) < 3 || nrow(core_out) < 10) {
    cat(sprintf("  SKIP %s ratio LMM: too few rows\n", st))
    return(NULL)
  }

  fit <- fit_lmm(log_ratio ~ treatment + (1 | mouse_id),
                 core_out, paste(st, "Core/Outside"))
  res <- extract_contrasts(fit, model_label = paste(st, "Core/Outside"))
  res[, staining := st]
  res[, zone     := "Core/Outside (ratio)"]
  res
}), fill = TRUE)

if (nrow(ratio_res) > 0) {
  ratio_res[, staining := factor(staining, levels = STAIN_ORDER)]
  ratio_res <- add_fdr(ratio_res, suffix = "_adj", by_col = "staining")

  cat("\n── Part B results ───────────────────────────────────────────────────────────\n")
  print(ratio_res[, .(staining, treatment, beta, pct_change, p, p_adj_adj, sig_adj)])
  fwrite(ratio_res, file.path(OUT_DIR, "zone_density_ratio_results.csv"))
  cat("  Saved: zone_density_ratio_results.csv\n\n")
}

# =============================================================================
# PART C — ZONE × TREATMENT INTERACTION
# =============================================================================
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("PART C: Zone × treatment interaction\n")
cat("        log(density_mm2) ~ treatment * zone + (1|mouse_id), per staining\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

ix_res <- rbindlist(lapply(STAIN_ORDER, function(st) {
  d <- dt[staining == st]
  cat(sprintf("── %s (%d rows) ──\n", st, nrow(d)))
  fit <- fit_lmm(log_dens ~ treatment * zone + (1 | mouse_id), d,
                 paste(st, "interaction"))
  anova_ix <- anova(fit, ddf = "Satterthwaite")
  ix_p <- anova_ix["treatment:zone", "Pr(>F)"]
  cat(sprintf("  treatment:zone F p = %.4f  %s\n\n", ix_p,
      ifelse(ix_p < 0.05, "✓ Zone-specific effect", "NS")))

  emm  <- emmeans(fit, ~ treatment | zone)
  con  <- contrast(emm, method = "trt.vs.ctrl", ref = "mScarlet")
  res  <- as.data.table(summary(con, infer = TRUE))
  old_names <- intersect(names(res),
    c("estimate","SE","lower.CL","upper.CL","asymp.LCL","asymp.UCL","t.ratio","z.ratio","p.value"))
  new_names <- c(
    estimate  = "beta", SE = "se",
    lower.CL  = "ci_lo", upper.CL = "ci_hi",
    asymp.LCL = "ci_lo", asymp.UCL = "ci_hi",
    t.ratio   = "t",    z.ratio   = "t",
    p.value   = "p"
  )[old_names]
  setnames(res, old_names, unname(new_names))
  res[, treatment   := gsub(" - mScarlet", "", contrast)]
  res[, pct_change  := (exp(beta) - 1) * 100]
  res[, staining    := st]
  res[, interaction_p := ix_p]
  res
}), fill = TRUE)

if (nrow(ix_res) > 0) {
  ix_res[, zone     := factor(zone,     levels = ZONE_ORDER)]
  ix_res[, staining := factor(staining, levels = STAIN_ORDER)]
  ix_res <- add_fdr(ix_res, suffix = "_adj", by_col = c("staining","zone"))

  cat("\n── Part C: significant contrasts ───────────────────────────────────────────\n")
  sig_c <- ix_res[sig_adj %in% c("***","**","*",".")]
  if (nrow(sig_c) > 0) {
    print(sig_c[, .(staining, zone, treatment, beta, pct_change, p, p_adj_adj, sig_adj)])
  } else {
    cat("  No hits reaching p_adj < 0.10.\n")
  }

  fwrite(ix_res, file.path(OUT_DIR, "zone_density_interaction.csv"))
  cat("  Saved: zone_density_interaction.csv\n\n")
}

# =============================================================================
# FIGURES
# =============================================================================

# ── Density profile line plot ──────────────────────────────────────────────────
profile_means <- dt[, .(
  mean_dens = mean(density_mm2, na.rm = TRUE),
  se_dens   = sd(density_mm2,   na.rm = TRUE) / sqrt(.N)
), by = .(staining, treatment, zone)]
profile_means[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

p1 <- ggplot(profile_means,
             aes(x = zone, y = mean_dens,
                 colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_dens - se_dens, ymax = mean_dens + se_dens),
                width = 0.15, linewidth = 0.6) +
  facet_wrap(~ staining, scales = "free_y", nrow = 1,
             labeller = labeller(staining = c(WFA = "WFA (PNN)", PV = "PV"))) +
  scale_colour_manual(values = PALETTE, name = NULL) +
  scale_y_continuous(name = "Mean cell density (cells/mm²)") +
  scale_x_discrete(name = "Injection zone") +
  labs(
    title    = "Cell density across injection zones",
    subtitle = "Mean ± SEM across animals."
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank())

# ── Core/Outside ratio forest plot (WFA — primary outcome) ────────────────────
if (nrow(ratio_res) > 0) {
  plot_ratio <- ratio_res[staining == "WFA"]
  plot_ratio[, treatment := factor(treatment, levels = rev(TREATMENT_ORDER[-1]))]

  p2 <- ggplot(plot_ratio,
               aes(x = pct_change, y = treatment, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = (exp(ci_lo)-1)*100, xmax = (exp(ci_hi)-1)*100),
                  width = 0.25, linewidth = 0.9) +
    geom_point(size = 4) +
    geom_text(aes(x = (exp(ci_hi)-1)*100 + 2, label = sig_adj),
              hjust = 0, size = 4, show.legend = FALSE) +
    scale_colour_manual(values = PALETTE[-1], guide = "none") +
    scale_x_continuous(
      name   = "% change in Core/Outside density ratio vs reference",
      expand = expansion(mult = c(0.05, 0.20))) +
    scale_y_discrete(name = NULL) +
    labs(
      title    = "Core/Outside WFA cell density ratio",
      subtitle = "Negative = fewer WFA cells at injection core relative to distal tissue.\nBH-FDR corrected per staining."
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

# ── Interaction heatmap (WFA) ─────────────────────────────────────────────────
if (nrow(ix_res) > 0) {
  ix_plot <- ix_res[staining == "WFA"]
  ix_plot[, treatment_f := factor(treatment, levels = TREATMENT_ORDER[-1])]

  p3 <- ggplot(ix_plot,
               aes(x = zone, y = treatment_f, fill = pct_change)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig_adj != "ns", sig_adj, "")),
              size = 4, fontface = "bold") +
    scale_fill_gradient2(low  = "#d73027", mid = "white", high = "#4575b4",
                         midpoint = 0, name = "% change\nvs reference") +
    scale_x_discrete(name = "Injection zone") +
    scale_y_discrete(name = NULL) +
    labs(
      title    = "Zone × treatment interaction (WFA density)",
      subtitle = sprintf(
        "treatment:zone F p = %.4f  %s\nStars = BH-FDR < 0.05 per zone.",
        unique(ix_plot$interaction_p),
        ifelse(unique(ix_plot$interaction_p) < 0.05, "✓ Zone-specific", "NS"))
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title  = element_text(face = "bold"),
          panel.grid  = element_blank())
}

# ── Residual check ────────────────────────────────────────────────────────────
fit_resid <- lmer(log_dens ~ treatment + (1 | mouse_id),
                  data = dt[staining == "WFA" & zone == "Core"],
                  control = lmerControl(optimizer = "bobyqa"))
p_resid <- ggplot(
  data.frame(fitted = fitted(fit_resid), residual = residuals(fit_resid)),
  aes(x = fitted, y = residual)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, colour = "red", linewidth = 0.7) +
  labs(title = "Residual check — WFA Core density LMM",
       x = "Fitted", y = "Residuals") +
  theme_bw(base_size = 10)

# ── Save figures ─────────────────────────────────────────────────────────────
for (ext in c("pdf","png")) {
  dpi <- ifelse(ext == "pdf", 300, 150)
  ggsave(file.path(OUT_DIR, paste0("fig_density_profile.", ext)),
         p1, width = 11, height = 5, dpi = dpi)
  if (nrow(ratio_res) > 0)
    ggsave(file.path(OUT_DIR, paste0("fig_density_ratio.", ext)),
           p2, width = 9, height = 5, dpi = dpi)
  if (nrow(ix_res) > 0)
    ggsave(file.path(OUT_DIR, paste0("fig_density_heatmap.", ext)),
           p3, width = 7, height = 5, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0("residual_check_00c.", ext)),
         p_resid, width = 6, height = 4, dpi = dpi)
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("SCRIPT 00c SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Part A — Per-zone hits (FDR, any zone × staining): %d\n",
    nrow(zone_res[sig_adj %in% c("***","**","*",".")])))
if (nrow(ratio_res) > 0)
  cat(sprintf("  Part B — Core/Outside ratio hits (FDR): %d\n",
      nrow(ratio_res[sig_adj %in% c("***","**","*",".")])))
if (nrow(ix_res) > 0) {
  for (st in STAIN_ORDER) {
    ip <- unique(ix_res[staining == st]$interaction_p)
    if (length(ip) > 0)
      cat(sprintf("  Part C — %s treatment:zone F p = %.4f  %s\n",
          st, ip, ifelse(ip < 0.05, "✓ Zone-specific", "NS")))
  }
}
cat(sprintf("\n── Outputs → %s\n", OUT_DIR))
cat("   zone_density_lmm_results.csv  zone_density_ratio_results.csv\n")
cat("   zone_density_interaction.csv  fig_density_* (pdf + png)\n")
