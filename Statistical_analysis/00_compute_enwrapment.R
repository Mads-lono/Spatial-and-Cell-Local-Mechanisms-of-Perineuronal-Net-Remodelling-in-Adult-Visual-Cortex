# =============================================================================
# 00_compute_enwrapment.R
# =============================================================================
# SOURCE: cells_with_zones.csv (raw output from CPN detection pipeline)
# PURPOSE: Compute PV enwrapment by PNNs from raw cell coordinates.
#          A PV cell is classified as "enwrapped" if its nearest PNN cell is
#          within COLOC_THRESH pixels (Euclidean distance, hires pixel space).
#
# OUTPUTS (all in RESULTS_DIR/00_enwrapment/):
#   slice_enwrapment.csv      — per slice × hemisphere: frac, n_pv, ratio
#   animal_enwrapment.csv     — per animal: mean ratio across slices (primary)
#   zone_enwrapment.csv       — per slice × zone × hemisphere: for zone analyses
#   atlas_enwrapment.csv      — per animal × area × layer: for atlas analyses
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(RANN)       # nn2() = fast approximate nearest neighbour
  library(stringr)
})

# ── Parameters ────────────────────────────────────────────────────────────────
CELLS_CSV   <- "/path/to/cells_with_zones.csv"
RESULTS_DIR <- "/path/to/results"

COLOC_THRESH  <- 30L   # pixels — PV is enwrapped if nearest PNN ≤ this
MIN_PV_SLICE  <- 5L    # minimum PV cells per slice × hemisphere for ratio
MIN_PV_ANIMAL <- 20L   # minimum PV cells pooled per animal × area × layer

TREATMENTS <- c("ADAMTS4", "ADAMTS4_MD", "ADAMTS15",
                "C6ST1", "C6ST1_ADAMTS15", "mScarlet")

# Brain area filter — adjust regex to match your atlas naming convention
VIS_AREA_REGEX <- "isual"

# Hemisphere assumed to be injection side; "left" is the default here —
# adjust to match your experimental convention
IPSI_HEMI  <- "left"
CONTRA_HEMI <- "right"

# ── Output directory ──────────────────────────────────────────────────────────
out_dir <- file.path(RESULTS_DIR, "00_enwrapment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output → %s\n", out_dir))

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("\nLoading cells_with_zones.csv ...\n")
t0 <- proc.time()

LOAD_COLS <- c("mouse_id", "slice_id", "cell_type", "hemisphere",
               "x_hires", "y_hires", "brain_area", "zone", "treatment")

dt <- fread(CELLS_CSV, select = LOAD_COLS)
setnames(dt, "cell_type", "stain")

cat(sprintf("  Loaded %s rows in %.1fs\n",
            format(nrow(dt), big.mark = ","),
            (proc.time() - t0)[["elapsed"]]))

dt <- dt[str_detect(brain_area, regex(VIS_AREA_REGEX, ignore_case = TRUE))]
dt <- dt[treatment %in% TREATMENTS]
cat(sprintf("  After region + treatment filter: %s rows\n",
            format(nrow(dt), big.mark = ",")))

pv  <- dt[stain == "PV"]
pnn <- dt[stain == "PNN"]
cat(sprintf("  PV: %s  |  PNN: %s\n",
            format(nrow(pv),  big.mark = ","),
            format(nrow(pnn), big.mark = ",")))

# =============================================================================
# 2. CORE ENWRAPMENT FUNCTION
# =============================================================================
# For a set of PV cells and PNN cells (same slice + hemisphere),
# returns the fraction of PV cells with nearest PNN ≤ COLOC_THRESH.
# Uses RANN::nn2() with k=1 (equivalent to scipy.spatial.cKDTree.query).
compute_frac_enwrapped <- function(pv_xy, pnn_xy) {
  if (nrow(pv_xy) == 0L || nrow(pnn_xy) == 0L) return(NA_real_)
  nn <- nn2(data = pnn_xy, query = pv_xy, k = 1L)
  mean(nn$nn.dists[, 1L] <= COLOC_THRESH)
}

# =============================================================================
# 3. SLICE-LEVEL ENWRAPMENT
# =============================================================================
cat("\n── Computing slice-level enwrapment ──────────────────────────────────────\n")

slice_keys <- unique(pv[, .(treatment, mouse_id, slice_id)])
cat(sprintf("  Processing %d slice × animal combinations...\n", nrow(slice_keys)))

t1 <- proc.time()
slice_records <- vector("list", nrow(slice_keys))

for (i in seq_len(nrow(slice_keys))) {
  key    <- slice_keys[i]
  pv_sl  <- pv[mouse_id == key$mouse_id & slice_id == key$slice_id]
  pnn_sl <- pnn[mouse_id == key$mouse_id & slice_id == key$slice_id]

  rec <- list(
    treatment = key$treatment,
    mouse_id  = key$mouse_id,
    slice_id  = key$slice_id
  )

  for (hemi in c(IPSI_HEMI, CONTRA_HEMI)) {
    label  <- if (hemi == IPSI_HEMI) "ipsi" else "contra"
    pv_h   <- pv_sl[hemisphere == hemi,  .(x_hires, y_hires)]
    pnn_h  <- pnn_sl[hemisphere == hemi, .(x_hires, y_hires)]
    n_pv   <- nrow(pv_h)
    frac   <- if (n_pv >= MIN_PV_SLICE && nrow(pnn_h) > 0L)
                compute_frac_enwrapped(pv_h, pnn_h)
              else NA_real_

    rec[[paste0("n_pv_", label)]] <- n_pv
    rec[[paste0("frac_", label)]] <- frac
  }

  frac_i <- rec[["frac_ipsi"]]
  frac_c <- rec[["frac_contra"]]
  rec[["ratio"]] <- if (!is.na(frac_i) && !is.na(frac_c) && frac_c > 0)
                      frac_i / frac_c
                    else NA_real_

  slice_records[[i]] <- rec

  if (i %% 100L == 0L)
    cat(sprintf("  %d / %d slices (%.0fs elapsed)\n",
                i, nrow(slice_keys), (proc.time() - t1)[["elapsed"]]))
}

slice_df <- rbindlist(slice_records)
slice_df <- slice_df[!is.na(ratio)]
cat(sprintf("  Valid slice ratios: %d / %d\n", nrow(slice_df), nrow(slice_keys)))
fwrite(slice_df, file.path(out_dir, "slice_enwrapment.csv"))
cat("  Saved: slice_enwrapment.csv\n")

# =============================================================================
# 4. ANIMAL-LEVEL ENWRAPMENT
# =============================================================================
cat("\n── Computing animal-level means ──────────────────────────────────────────\n")

animal_df <- slice_df[, .(
  ratio     = mean(ratio,     na.rm = TRUE),
  frac_ipsi = mean(frac_ipsi, na.rm = TRUE),
  n_slices  = .N
), by = .(treatment, mouse_id)]

fwrite(animal_df, file.path(out_dir, "animal_enwrapment.csv"))
cat(sprintf("  Animal records: %d\n", nrow(animal_df)))
cat("  Saved: animal_enwrapment.csv\n")
print(animal_df[order(treatment, mouse_id)])

# =============================================================================
# 5. ZONE-LEVEL ENWRAPMENT
# =============================================================================
cat("\n── Computing zone-level enwrapment ───────────────────────────────────────\n")
# NOTE: The contralateral hemisphere has no injection zone labels —
# all its cells are "Outside" because the vector never reached it.
# Outcome here is frac_enwrapped = fraction of ipsilateral PV cells enwrapped
# within each zone. Between-animal normalisation is handled downstream by LMMs.

zone_keys <- unique(pv[!is.na(zone) & zone != "",
                        .(treatment, mouse_id, slice_id, zone)])
cat(sprintf("  Processing %d slice × zone combinations...\n", nrow(zone_keys)))

t2 <- proc.time()
zone_records <- vector("list", nrow(zone_keys))

for (i in seq_len(nrow(zone_keys))) {
  key   <- zone_keys[i]
  pv_i  <- pv[mouse_id == key$mouse_id & slice_id == key$slice_id &
               zone == key$zone & hemisphere == IPSI_HEMI, .(x_hires, y_hires)]
  pnn_i <- pnn[mouse_id == key$mouse_id & slice_id == key$slice_id &
                zone == key$zone & hemisphere == IPSI_HEMI, .(x_hires, y_hires)]

  n_pv <- nrow(pv_i)
  if (n_pv < MIN_PV_SLICE || nrow(pnn_i) == 0L) {
    zone_records[[i]] <- NULL
    next
  }

  frac <- compute_frac_enwrapped(pv_i, pnn_i)
  zone_records[[i]] <- list(
    treatment      = key$treatment,
    mouse_id       = key$mouse_id,
    slice_id       = key$slice_id,
    zone           = key$zone,
    frac_enwrapped = frac,
    n_pv_ipsi      = n_pv
  )

  if (i %% 200L == 0L)
    cat(sprintf("  %d / %d zone-slices (%.0fs)\n",
                i, nrow(zone_keys), (proc.time() - t2)[["elapsed"]]))
}

zone_df <- rbindlist(Filter(Negate(is.null), zone_records))
zone_df <- zone_df[!is.na(frac_enwrapped)]
fwrite(zone_df, file.path(out_dir, "zone_enwrapment.csv"))
cat(sprintf("  Valid zone records: %d\n", nrow(zone_df)))
cat("  Saved: zone_enwrapment.csv\n")

cat("\n  Zone × treatment coverage (n slices):\n")
print(dcast(zone_df[, .N, by = .(treatment, zone)],
            treatment ~ zone, value.var = "N", fill = 0L))

# =============================================================================
# 6. ATLAS-LEVEL ENWRAPMENT  (area × layer per animal)
# =============================================================================
cat("\n── Computing atlas-stratified enwrapment (area × layer) ──────────────────\n")
cat("  (Pools all slices per animal × area × layer — avoids per-slice sparsity)\n")

# Parse area and layer from brain_area string (e.g. "Primary visual area, layer 5")
pv[,  c("area", "layer") := {
  parts <- str_split_fixed(brain_area, ", ", 2)
  list(str_trim(parts[, 1]), str_to_lower(str_trim(parts[, 2])))
}]
pnn[, c("area", "layer") := {
  parts <- str_split_fixed(brain_area, ", ", 2)
  list(str_trim(parts[, 1]), str_to_lower(str_trim(parts[, 2])))
}]

atlas_keys <- unique(pv[, .(treatment, mouse_id, area, layer)])
cat(sprintf("  Processing %d animal × area × layer strata...\n", nrow(atlas_keys)))

t3 <- proc.time()
atlas_records <- vector("list", nrow(atlas_keys))

for (i in seq_len(nrow(atlas_keys))) {
  key   <- atlas_keys[i]
  pv_k  <- pv[mouse_id == key$mouse_id & area == key$area & layer == key$layer]
  pnn_k <- pnn[mouse_id == key$mouse_id & area == key$area & layer == key$layer]

  rec <- list(
    treatment = key$treatment,
    mouse_id  = key$mouse_id,
    area      = key$area,
    layer     = key$layer
  )

  ok <- TRUE
  for (hemi in c(IPSI_HEMI, CONTRA_HEMI)) {
    label  <- if (hemi == IPSI_HEMI) "ipsi" else "contra"
    pv_h   <- pv_k[hemisphere == hemi,  .(x_hires, y_hires)]
    pnn_h  <- pnn_k[hemisphere == hemi, .(x_hires, y_hires)]
    n_pv   <- nrow(pv_h)

    if (n_pv < MIN_PV_ANIMAL || nrow(pnn_h) == 0L) { ok <- FALSE; break }

    frac <- compute_frac_enwrapped(pv_h, pnn_h)
    rec[[paste0("n_pv_", label)]] <- n_pv
    rec[[paste0("frac_", label)]] <- frac
  }

  if (!ok) { atlas_records[[i]] <- NULL; next }

  frac_i <- rec[["frac_ipsi"]]
  frac_c <- rec[["frac_contra"]]
  rec[["ratio"]] <- if (!is.na(frac_i) && !is.na(frac_c) && frac_c > 0)
                      frac_i / frac_c
                    else NA_real_

  atlas_records[[i]] <- rec

  if (i %% 500L == 0L)
    cat(sprintf("  %d / %d strata (%.0fs)\n",
                i, nrow(atlas_keys), (proc.time() - t3)[["elapsed"]]))
}

atlas_df <- rbindlist(Filter(Negate(is.null), atlas_records))
atlas_df <- atlas_df[!is.na(ratio)]

# Abbreviate area names for plots — edit to match your atlas naming convention
AREA_SHORT <- c(
  "Primary visual area"         = "VISp",
  "Lateral visual area"         = "VISl",
  "Anterolateral visual area"   = "VISal",
  "posteromedial visual area"   = "VISpm",
  "Anteromedial visual area"    = "VISam",
  "Posterolateral visual area"  = "VISpl"
)
atlas_df[, area_short := AREA_SHORT[area]]
atlas_df[is.na(area_short), area_short := substr(area, 1, 5)]

fwrite(atlas_df, file.path(out_dir, "atlas_enwrapment.csv"))
cat(sprintf("  Valid atlas records: %d\n", nrow(atlas_df)))
cat("  Saved: atlas_enwrapment.csv\n")

# =============================================================================
# 7. SUMMARY
# =============================================================================
total_elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("\n══ All outputs written to %s\n", out_dir))
cat(sprintf("   slice_enwrapment.csv  : %d rows\n",  nrow(slice_df)))
cat(sprintf("   animal_enwrapment.csv : %d rows\n",  nrow(animal_df)))
cat(sprintf("   zone_enwrapment.csv   : %d rows\n",  nrow(zone_df)))
cat(sprintf("   atlas_enwrapment.csv  : %d rows\n",  nrow(atlas_df)))
cat(sprintf("   Total elapsed         : %.1f min\n", total_elapsed / 60))
