# =============================================================================
# 16_loo_sensitivity.R
# Leave-one-out Bayesian sensitivity analysis
#
# Reruns the two primary Bayesian models from script 13 with each animal
# removed in turn, verifying that no single animal drives the key posteriors.
#
# Models:
#   LOO-M1: frac_enwrapped ~ treatment + (1|mouse_id)
#   LOO-M2: composite ~ treatment + (1|mouse_id)
#
# Parallelisation: N_PARALLEL LOO iterations run simultaneously,
# each using MCMC_CORES chains.
#
# Summary outputs:
#   loo_stability_enwrap.csv    — p_negative across all LOO iterations (M1)
#   loo_stability_composite.csv — p_negative across all LOO iterations (M2)
#   loo_influence.csv           — max change in p_negative per treatment
#   loo_diagnostics.csv         — R-hat and divergences for all LOO models
#   fig_loo_stability.pdf       — stability plots and influence heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(brms)
  library(ggplot2)
  library(patchwork)
})

if (!requireNamespace("posterior", quietly = TRUE))
  install.packages("posterior", repos = "https://cloud.r-project.org")
library(posterior)

if (!requireNamespace("future", quietly = TRUE))
  install.packages("future", repos = "https://cloud.r-project.org")
if (!requireNamespace("future.apply", quietly = TRUE))
  install.packages("future.apply", repos = "https://cloud.r-project.org")
library(future)
library(future.apply)

# ── Paths ─────────────────────────────────────────────────────────────────────
COMPOSITE_CSV <- "/path/to/results/12_pnn_integrity_composite/composite_section_data.csv"
RESULTS_DIR   <- "/path/to/results"
OUT_DIR       <- file.path(RESULTS_DIR, "16_loo_sensitivity")
LOO_CACHE_DIR <- file.path(OUT_DIR, "loo_cache")
dir.create(OUT_DIR,       recursive = TRUE, showWarnings = FALSE)
dir.create(LOO_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENTS      <- c("mScarlet","ADAMTS4","ADAMTS4_MD",
                     "ADAMTS15","C6ST1","C6ST1_ADAMTS15")
TREATMENT_ORDER <- TREATMENTS
MCMC_CHAINS     <- 4L    # chains per model
MCMC_ITER       <- 2000L
MCMC_WARMUP     <- 1000L
MCMC_CORES      <- 4L    # cores per model (= MCMC_CHAINS)
N_PARALLEL      <- 4L    # LOO models run simultaneously
MCMC_SEED       <- 2025L
ADAPT_DELTA     <- 0.95
ROPE_ENWRAP     <- 0.05
ROPE_COMPOSITE  <- 0.10

PALETTE <- c(
  mScarlet        = "#888888",
  ADAMTS4         = "#4e9af1",
  ADAMTS4_MD      = "#f17c4e",
  ADAMTS15        = "#4ef196",
  C6ST1           = "#c44ef1",
  C6ST1_ADAMTS15  = "#f1c44e"
)

priors_weak <- c(
  prior(normal(0, 1), class = b),
  prior(normal(0, 1), class = sd),
  prior(normal(0, 1), class = sigma)
)

cat("Script 16: Leave-one-out Bayesian sensitivity analysis\n")
cat(sprintf("Chains: %d | Iter: %d | Warmup: %d\n",
    MCMC_CHAINS, MCMC_ITER, MCMC_WARMUP))
cat(sprintf("Parallel LOO: %d models simultaneously (%d total cores)\n\n",
    N_PARALLEL, N_PARALLEL * MCMC_CORES))

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading section data...\n")
dt <- fread(COMPOSITE_CSV)
dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
cat(sprintf("  %d sections, %d animals\n\n",
    nrow(dt), uniqueN(dt$mouse_id)))

all_animals <- unique(dt$mouse_id)
cat(sprintf("Animals (%d total):\n", length(all_animals)))
cat(paste(sort(all_animals), collapse = ", "), "\n\n")

# ── Helpers ───────────────────────────────────────────────────────────────────
extract_post_loo <- function(fit, model_name, left_out, rope_bound) {
  if (is.null(fit)) return(NULL)
  draws      <- as_draws_df(fit)
  treat_cols <- grep("^b_treatment", names(draws), value = TRUE)
  rbindlist(lapply(treat_cols, function(col) {
    vals  <- draws[[col]]
    treat <- str_remove(col, "^b_treatment")
    ci89  <- quantile(vals, c(0.055, 0.945))
    data.table(
      model        = model_name,
      left_out     = left_out,
      treatment    = treat,
      post_mean    = round(mean(vals),  4),
      post_sd      = round(sd(vals),    4),
      ci_89_lo     = round(ci89[[1]],   4),
      ci_89_hi     = round(ci89[[2]],   4),
      p_negative   = round(mean(vals < 0),           3),
      p_meaningful = round(mean(vals < -rope_bound), 3)
    )
  }), fill = TRUE)
}

extract_diag_loo <- function(fit, model_name, left_out) {
  if (is.null(fit)) return(NULL)
  rh   <- rhat(fit)
  ndiv <- sum(nuts_params(fit)$Value[
    nuts_params(fit)$Parameter == "divergent__"])
  data.table(
    model       = model_name,
    left_out    = left_out,
    max_rhat    = round(max(rh, na.rm = TRUE), 4),
    n_divergent = ndiv,
    converged   = max(rh, na.rm = TRUE) < 1.01 & ndiv == 0
  )
}

# ── Full-data reference models ────────────────────────────────────────────────
cat("Fitting full-data reference models...\n")

m1_full <- brm(
  frac_enwrapped ~ treatment + (1 | mouse_id),
  data = dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(LOO_CACHE_DIR, "m1_full"), silent = 2
)

m2_full <- brm(
  composite ~ treatment + (1 | mouse_id),
  data = dt, prior = priors_weak,
  chains = MCMC_CHAINS, iter = MCMC_ITER, warmup = MCMC_WARMUP,
  cores = MCMC_CORES, seed = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  file = file.path(LOO_CACHE_DIR, "m2_full"), silent = 2
)

full_m1 <- extract_post_loo(m1_full, "M1_enwrap",    "none", ROPE_ENWRAP)
full_m2 <- extract_post_loo(m2_full, "M2_composite", "none", ROPE_COMPOSITE)

cat("\nFull-data M1 (enwrapment):\n")
print(full_m1[, .(treatment, post_mean, p_negative, p_meaningful)])
cat("\nFull-data M2 (composite):\n")
print(full_m2[, .(treatment, post_mean, p_negative, p_meaningful)])

# ── Parallel LOO ──────────────────────────────────────────────────────────────
cat(sprintf("\nStarting parallel LOO (%d workers)...\n", N_PARALLEL))
future::plan(multisession, workers = N_PARALLEL)

loo_results <- future.apply::future_lapply(
  seq_along(all_animals),
  function(i_anim) {

    suppressPackageStartupMessages({
      library(data.table)
      library(stringr)
      library(brms)
      library(posterior)
    })

    anim   <- all_animals[i_anim]
    dt_loo <- dt[mouse_id != anim]
    dt_loo[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

    treat_counts <- dt_loo[, .N, by = treatment]
    if (any(treat_counts$N == 0)) {
      message(sprintf("  SKIP %s: treatment group empty", anim))
      return(NULL)
    }

    safe_name <- str_replace_all(anim, "[/_]", "-")

    m1_loo <- tryCatch(
      brm(frac_enwrapped ~ treatment + (1 | mouse_id),
          data    = dt_loo,
          prior   = priors_weak,
          chains  = MCMC_CHAINS, iter    = MCMC_ITER,
          warmup  = MCMC_WARMUP, cores   = MCMC_CORES,
          seed    = MCMC_SEED,
          control = list(adapt_delta = ADAPT_DELTA),
          file    = file.path(LOO_CACHE_DIR, sprintf("m1_loo_%s", safe_name)),
          silent  = 2),
      error = function(e) {
        message(sprintf("M1 error [%s]: %s", anim, e$message)); NULL
      }
    )

    m2_loo <- tryCatch(
      brm(composite ~ treatment + (1 | mouse_id),
          data    = dt_loo,
          prior   = priors_weak,
          chains  = MCMC_CHAINS, iter    = MCMC_ITER,
          warmup  = MCMC_WARMUP, cores   = MCMC_CORES,
          seed    = MCMC_SEED,
          control = list(adapt_delta = ADAPT_DELTA),
          file    = file.path(LOO_CACHE_DIR, sprintf("m2_loo_%s", safe_name)),
          silent  = 2),
      error = function(e) {
        message(sprintf("M2 error [%s]: %s", anim, e$message)); NULL
      }
    )

    list(
      anim    = anim,
      post_m1 = extract_post_loo(m1_loo, "M1_enwrap",    anim, ROPE_ENWRAP),
      post_m2 = extract_post_loo(m2_loo, "M2_composite", anim, ROPE_COMPOSITE),
      diag_m1 = extract_diag_loo(m1_loo, "M1_enwrap",    anim),
      diag_m2 = extract_diag_loo(m2_loo, "M2_composite", anim)
    )
  },
  future.seed = TRUE
)

future::plan(sequential)
cat("LOO iterations complete.\n\n")

# ── Unpack results ────────────────────────────────────────────────────────────
loo_results_clean <- Filter(Negate(is.null), loo_results)

all_loo_m1 <- Filter(Negate(is.null), lapply(loo_results_clean, `[[`, "post_m1"))
all_loo_m2 <- Filter(Negate(is.null), lapply(loo_results_clean, `[[`, "post_m2"))
all_diag   <- Filter(Negate(is.null),
  c(lapply(loo_results_clean, `[[`, "diag_m1"),
    lapply(loo_results_clean, `[[`, "diag_m2")))

# ── Compile and save ──────────────────────────────────────────────────────────
loo_m1_dt <- rbind(full_m1, rbindlist(all_loo_m1, fill = TRUE), fill = TRUE)
loo_m2_dt <- rbind(full_m2, rbindlist(all_loo_m2, fill = TRUE), fill = TRUE)
diag_dt   <- rbindlist(all_diag, fill = TRUE)

fwrite(loo_m1_dt, file.path(OUT_DIR, "loo_stability_enwrap.csv"))
fwrite(loo_m2_dt, file.path(OUT_DIR, "loo_stability_composite.csv"))
fwrite(diag_dt,   file.path(OUT_DIR, "loo_diagnostics.csv"))

# ── Influence analysis ────────────────────────────────────────────────────────
compute_influence <- function(loo_dt, model_name, full_dt) {
  rbindlist(lapply(TREATMENTS[-1], function(treat) {
    full_p   <- full_dt[treatment == treat, p_negative]
    loo_vals <- loo_dt[left_out != "none" & treatment == treat, p_negative]
    if (length(full_p) == 0 || length(loo_vals) == 0) return(NULL)
    most_inf_idx  <- which.max(abs(loo_vals - full_p))
    most_inf_anim <- loo_dt[
      left_out != "none" & treatment == treat][most_inf_idx, left_out]
    data.table(
      model            = model_name,
      treatment        = treat,
      full_p_neg       = round(full_p, 3),
      min_loo_p_neg    = round(min(loo_vals), 3),
      max_loo_p_neg    = round(max(loo_vals), 3),
      max_change       = round(max(abs(loo_vals - full_p)), 3),
      most_influential = most_inf_anim
    )
  }), fill = TRUE)
}

inf_m1 <- compute_influence(loo_m1_dt, "M1_enwrap",    full_m1)
inf_m2 <- compute_influence(loo_m2_dt, "M2_composite", full_m2)
influence_dt <- rbind(inf_m1, inf_m2, fill = TRUE)
fwrite(influence_dt, file.path(OUT_DIR, "loo_influence.csv"))

cat("── LOO influence summary ────────────────────────────────────────────────\n")
print(influence_dt[order(model, -max_change)])

# ── Stability report ──────────────────────────────────────────────────────────
cat("\n── Stability: ADAMTS15 and C6ST1_ADAMTS15 ──────────────────────────────\n")
for (model_name in c("M1_enwrap","M2_composite")) {
  dt_use <- if (model_name == "M1_enwrap") loo_m1_dt else loo_m2_dt
  cat(sprintf("\n%s:\n", model_name))
  for (treat in c("ADAMTS15","C6ST1_ADAMTS15")) {
    sub <- dt_use[treatment == treat]
    cat(sprintf("  %s:\n", treat))
    cat(sprintf("    Full:      p_neg=%.3f | mean=%.4f\n",
        sub[left_out == "none", p_negative],
        sub[left_out == "none", post_mean]))
    cat(sprintf("    LOO range: p_neg=[%.3f, %.3f] | mean=[%.4f, %.4f]\n",
        sub[left_out != "none", min(p_negative)],
        sub[left_out != "none", max(p_negative)],
        sub[left_out != "none", min(post_mean)],
        sub[left_out != "none", max(post_mean)]))
    below_80 <- sub[left_out != "none" & p_negative < 0.80]
    if (nrow(below_80) > 0) {
      cat(sprintf("    WARNING: p_neg < 0.80 when removing: %s\n",
          paste(below_80$left_out, collapse = ", ")))
    } else {
      cat("    STABLE: p_neg >= 0.80 in all LOO variants\n")
    }
  }
}

# ── Convergence summary ───────────────────────────────────────────────────────
cat("\n── LOO convergence ─────────────────────────────────────────────────────\n")
n_conv  <- diag_dt[converged == TRUE,  .N]
n_total <- nrow(diag_dt)
cat(sprintf("  %d / %d models converged\n", n_conv, n_total))
if (any(!diag_dt$converged)) {
  cat("  Non-converged:\n")
  print(diag_dt[converged == FALSE, .(model, left_out, max_rhat, n_divergent)])
}

# ── Figures ───────────────────────────────────────────────────────────────────
cat("\nGenerating figures...\n")

plot_stability <- function(loo_dt, title_str, y_label) {
  active <- c("ADAMTS15","C6ST1_ADAMTS15","ADAMTS4_MD","ADAMTS4")
  sub    <- loo_dt[treatment %in% active]
  sub[, treatment := factor(treatment, levels = TREATMENT_ORDER)]
  sub[, is_full   := left_out == "none"]

  ggplot(sub, aes(x = treatment, y = p_negative, colour = treatment)) +
    geom_jitter(data = sub[is_full == FALSE],
                width = 0.12, size = 2.5, alpha = 0.65) +
    geom_point(data = sub[is_full == TRUE], size = 6, shape = 18) +
    geom_hline(yintercept = 0.80, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = 0.95, linetype = "dotted",
               colour = "grey40", linewidth = 0.4) +
    annotate("text", x = 0.55, y = 0.82, label = "0.80",
             size = 3, colour = "grey50") +
    annotate("text", x = 0.55, y = 0.97, label = "0.95",
             size = 3, colour = "grey50") +
    scale_colour_manual(values = PALETTE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title    = title_str,
         subtitle = "Diamond = full data  |  dots = LOO variants  |  dashed = 0.80",
         x = NULL, y = y_label) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

p1 <- plot_stability(loo_m1_dt,
  "LOO stability — enwrapment (M1)",
  "P(treatment reduces frac_enwrapped)")

p2 <- plot_stability(loo_m2_dt,
  "LOO stability — composite (M2)",
  "P(treatment reduces composite)")

# Heatmap
heat_dt <- loo_m1_dt[
  left_out != "none" &
  treatment %in% c("ADAMTS15","C6ST1_ADAMTS15","ADAMTS4","ADAMTS4_MD")
]
heat_dt[, treatment := factor(treatment, levels = TREATMENT_ORDER)]

p3 <- ggplot(heat_dt,
             aes(x = left_out, y = treatment, fill = p_negative)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", p_negative)), size = 2.8) +
  scale_fill_gradient2(
    low = "#d73027", mid = "lightyellow", high = "#1a9850",
    midpoint = 0.80, limits = c(0, 1), name = "P(negative)") +
  labs(title = "LOO p_negative heatmap — enwrapment model",
       x = "Animal left out", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7))

combined <- (p1 / p2 / p3) + plot_layout(heights = c(1, 1, 1.3))
ggsave(file.path(OUT_DIR, "fig_loo_stability.pdf"),
       combined, width = 10, height = 14)
ggsave(file.path(OUT_DIR, "fig_loo_stability.png"),
       combined, width = 10, height = 14, dpi = 150)
cat("Saved fig_loo_stability.pdf\n")

cat("\nScript 16 complete\n")
cat(sprintf("Outputs: %s\n", OUT_DIR))
cat("Stan .rds files cached — rerunning loads from cache.\n")