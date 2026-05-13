# AAV-Mediated ECM Remodelling in Mouse Primary Visual Cortex
### Analysis pipeline — master's thesis, University of Oslo

This repository contains all custom scripts used in a master's thesis examining the effects of AAV-mediated overexpression of extracellular matrix (ECM)-remodelling enzymes on perineuronal nets (PNNs) and parvalbumin-positive (PV) interneurons in the adult mouse primary visual cortex (V1).

---

## Study design

Six AAV constructs were injected unilaterally into V1 of adult male C57BL/6Rj mice (n = 4 per group, 24 animals total):

| Group | Construct | Role |
|---|---|---|
| mScarlet | Fluorescent reporter only | Control |
| ADAMTS-4 | Aggrecanase (active) | PNN degradation |
| ADAMTS-4-MD | Aggrecanase (active) + Monocular Deprivation | PNN remodelling |
| ADAMTS-15 | Aggrecanase (active) | PNN degradation |
| C6ST-1 | Chondroitin-6-sulfotransferase | PNN remodelling |
| C6ST-1 + ADAMTS-15 | Combined | Combinatorial remodelling |

Primary outcomes were **PV–PNN enwrapment fraction** (proportion of PV cells enwrapped by WFA-positive PNN) and **diffuse WFA fluorescence** across the injection zone and surrounding cortex.

---

## Repository structure

The pipeline has two sequential phases: image preprocessing (scripts `00a`–`12`) followed by statistical analysis (scripts `00`–`28`). These are separate numbering systems.

```
.
├── README.md
│
├── preprocessing/                   Image preprocessing pipeline
│   ├── 00a_vsi_to_tiff.py
│   ├── 00b_generate_quint_thumbnails.py
│   ├── 00c_split_and_convert_8bit.py
│   ├── 01_contrast_enhancement_imagej.py
│   ├── 02_flat_to_png.py
│   ├── 03_extract_c3.py
│   ├── 04_generate_injection_masks.py
│   ├── 05_tissue_masker.py
│   ├── 06_batch_pnn.sh
│   ├── 07_csvsplitter.py
│   ├── 08_diffu_spot_quantification.py
│   ├── 09_cell_fluorescence_analysis.py
│   ├── 10_merge_nested_results.py
│   ├── 11_assign_zones.py
│   └── 12_compute_zone_density.py
│
├── diagnostics/                     QC and one-time correction utilities
│   ├── Orientation_switch.py
│   ├── fix_c3_csv_coordinates.py
│   └── mask_qc.py
│
├── analysis/                        Statistical analysis pipeline
│   ├── 00_compute_enwrapment.R
│   ├── 00b_enwrapment_by_zone.R
│   ├── 01_primary_lmm.R
│   ├── 02_zone_gradient.R
│   ├── 03_atlas_stratified.R
│   ├── 08_off_target_survey.R
│   ├── 09_cell_intensity_lmm.R
│   ├── 11_density_lmm.R
│   ├── 12_pnn_integrity_composite.R
│   ├── 13_bayesian_lmm.R
│   ├── 14_scarlet_expression_covariate.py
│   ├── 15_spatial_gradient.py
│   ├── 16_loo_sensitivity.R
│   ├── 17_pv_intensity_lmm.R
│   ├── 19_ipsi_contra_ratio.R
│   ├── 20_adamts15_vs_c6st1_adamts15.R
│   ├── 21_substrate_conditioning.R
│   ├── 22_adamts4_md_effect.R
│   ├── 23_bayesian_additive.R
│   ├── 24_virus_cell_enwrapment.R
│   ├── 25_cell_level_deeper.R
│   ├── 26_rq3_cell_factorial.R
│   ├── 27_zone_density.R
│   └── 28_viral_load_heterogeneity.R
│
├── collect_all_results.py           Aggregates all result CSVs into one file
│
├── figures/                         Figure generation
   ├── thesis_figures.py
   ├── fig_zone_overlay.py
   ├── fig_enwrapment_comparison_40x.py
   └── fig_virus_proximal_40x.py

```

> Scripts 04–07 and 10 of the statistical analysis pipeline are not included; they address questions not reported in the thesis. Numbering preserved for script referencing. 

---

## Pipeline overview

### Phase 1 — Image preprocessing (`00a`–`12`)

Scripts are run once per animal in numerical order. Each script's header documents its inputs, outputs, and required software.

| Script | Description |
|---|---|
| `00a_vsi_to_tiff.py` | **Fiji/Jython.** Converts raw `.vsi` microscope files (Olympus CellSens) to TIFF with systematic `sXXX` slice indexing. Run inside Fiji. |
| `00b_generate_quint_thumbnails.py` | **Fiji/Jython.** Produces contrast-enhanced 1500 px JPEG thumbnails for atlas registration in QUINT (QuickNII / VisuAlign). Run inside Fiji. |
| `00c_split_and_convert_8bit.py` | Splits multi-channel 16-bit TIFFs into single-channel files (`-C1`, `-C2`, `-C3`) and converts to 8-bit. Source files are never modified. |
| `01_contrast_enhancement_imagej.py` | **Fiji/Jython.** Background subtraction (rolling ball) and histogram normalisation on C1 and C2 channels. Run inside Fiji. |
| `02_flat_to_png.py` | Decodes VisuAlign `.flat` binary atlas overlays to RGB PNGs using the Rainbow 2017 colour lookup table. |
| `03_extract_c3.py` | Extracts the injection-marker channel (C3) from each multi-channel TIFF into a flat folder for Cellpose input. |
| `04_generate_injection_masks.py` | Generates injection zone masks using Cellpose (CPSAM). Produces `MaskLoose` (penumbra boundary), `MaskStrict` (core), segmentation labels, and cell count CSVs per slice. |
| `05_tissue_masker.py` | Converts ilastik Simple Segmentation outputs to 1-bit binary tissue masks. |
| `06_batch_pnn.sh` | **Bash.** Batch CPN Faster R-CNN detection and rank-learning rescoring for all split 8-bit enhanced TIFFs. Run once per channel (PNN / PV). |
| `07_csvsplitter.py` | Filters CPN detections by rescore threshold and splits the localizations CSV into per-slice coordinate CSVs. |
| `08_diffu_spot_quantification.py` | Region-level diffuse fluorescence quantification (WFA and PV channels) across all brain areas and injection zones, with hemisphere assignment from VisuAlign anchoring vectors. |
| `09_cell_fluorescence_analysis.py` | Per-cell fluorescence intensity analysis: assigns each detected cell to a brain region and hemisphere, and extracts patch-level intensity metrics. |
| `10_merge_nested_results.py` | Merges per-animal diffuse fluorescence and cell fluorescence CSVs into combined hemisphere-level and zone-level datasets. |
| `11_assign_zones.py` | Assigns each cell to an injection zone (Core / Penumbra / Outside) by mask lookup at scaled coordinates. Produces `cells_with_zones.csv`. |
| `12_compute_zone_density.py` | Computes PNN and PV cell density (cells/mm²) per animal × brain area × zone. |

### Diagnostics

One-time QC and correction utilities run between preprocessing and analysis.

| Script | Description |
|---|---|
| `Orientation_switch.py` | 26-stage SSIM-based QC pipeline: matches generated thumbnails against reference thumbnails, resolves section ordering and orientation mismatches, applies horizontal flips, and produces `FINAL_TIFF_TRANSFORMATION_LOG.csv`. |
| `fix_c3_csv_coordinates.py` | Corrects `Global_X` coordinates in C3 Cellpose CSVs for slices that were horizontally flipped during orientation standardisation. |
| `mask_qc.py` | Identifies injection mask false positives from the distribution of MaskLoose white-pixel percentages per animal. Animals with predominantly large masks are likely injection failures. |

### Phase 2 — Statistical analysis (`00`–`28`)

All scripts read from `ALL_RESULTS_COMBINED.csv` or specific upstream CSVs. Scripts can be run selectively; each header documents upstream dependencies.

| Script | Description |
|---|---|
| `00_compute_enwrapment.R` | Primary enwrapment classification: nearest-neighbour distance from each PV cell to WFA-positive PNN centroids; computes per-animal and per-slice enwrapment fractions. |
| `00b_enwrapment_by_zone.R` | Zone-stratified enwrapment: repeats primary analysis split by injection zone. |
| `01_primary_lmm.R` | Primary slice-level linear mixed-effects model (log-transformed enwrapment ratio, random intercept for animal). |
| `02_zone_gradient.R` | Zone-by-treatment gradient model; tests whether enwrapment effects differ across injection zones. |
| `03_atlas_stratified.R` | Layer-stratified contrasts: fits LMMs separately for each cortical layer from the Allen Brain Atlas parcellation. |
| `08_off_target_survey.R` | Surveys enwrapment in off-target regions to assess injection specificity. |
| `09_cell_intensity_lmm.R` | LMM for WFA and PV cell-level fluorescence intensity. |
| `11_density_lmm.R` | LMM for PNN and PV cell density at hemisphere level. |
| `12_pnn_integrity_composite.R` | Composite PNN structural integrity analysis. |
| `13_bayesian_lmm.R` | Bayesian multilevel models (brms/Stan), including off-target specificity model. |
| `14_scarlet_expression_covariate.py` | Tests mScarlet expression level as a continuous covariate for enwrapment outcomes. |
| `15_spatial_gradient.py` | Distance-decay gradient of enwrapment fraction as a function of distance from the injection centre. |
| `16_loo_sensitivity.R` | Leave-one-out sensitivity analysis: refits the primary model excluding one animal at a time. |
| `17_pv_intensity_lmm.R` | LMM for PV interneuron fluorescence intensity. |
| `19_ipsi_contra_ratio.R` | Ipsilateral vs contralateral enwrapment ratio analysis. |
| `20_adamts15_vs_c6st1_adamts15.R` | Direct pairwise comparison of ADAMTS-15 and C6ST-1 + ADAMTS-15 groups. |
| `21_substrate_conditioning.R` | Tests whether prior C6ST-1 expression modifies ADAMTS-15 effects. |
| `22_adamts4_md_effect.R` | Isolates the effect of the catalytically inactive ADAMTS-4-MD construct. |
| `23_bayesian_additive.R` | Bayesian factorial model of the C6ST-1 × ADAMTS-15 interaction. |
| `24_virus_cell_enwrapment.R` | Cell-level analysis of enwrapment in virus-positive vs virus-negative PV cells. |
| `25_cell_level_deeper.R` | Deeper cell-level analyses including distance-bin enwrapment fractions. |
| `26_rq3_cell_factorial.R` | Factorial analysis for cell-level RQ3 outcomes. |
| `27_zone_density.R` | Zone-stratified cell density LMM (WFA and PV). |
| `28_viral_load_heterogeneity.R` | Models enwrapment as a function of local viral load as a continuous predictor. |

### Utility

| Script | Description |
|---|---|
| `collect_all_results.py` | Recursively finds all result CSVs, filters out large raw tables, and concatenates the remainder into `ALL_RESULTS_COMBINED.csv`. Run before figure generation. |

### Figure generation

| Script | Description |
|---|---|
| `thesis_figures.py` | All figures used in the thesis. Reads from `ALL_RESULTS_COMBINED.csv`. |
| `fig_zone_overlay.py` | Five-panel injection zone overlay figure: individual channels, RGB composite, and zone-coloured composite per slice. |
| `fig_enwrapment_comparison_40x.py` | Enwrapment comparison at 40×: full-slide overview with crop box, 40× crop with PV overlay coloured by enwrapment status, and enwrapment concept schematic. |
| `fig_virus_proximal_40x.py` | Virus-proximal / virus-distal classification figure at 40×, with concentric distance-bin rings and a two-panel proximity schematic. |

---

## Dependencies

### R packages
```r
install.packages(c("lme4", "lmerTest", "emmeans", "brms", "RANN",
                   "tidyverse", "data.table"))
```

Bayesian models (scripts `13`, `23`, `24`, `26`) require [Stan](https://mc-stan.org/), accessed via `brms`. Installation instructions at [mc-stan.org](https://mc-stan.org/users/interfaces/rstan).

### Python packages
```bash
pip install pandas numpy matplotlib scipy scikit-image tifffile pillow tqdm cellpose
```

### External tools
- **Fiji / ImageJ** — required for scripts `00a`, `00b`, `01` (Jython runtime)
- **CPN** (Cell Profiler Nucleus detector) — required for `06_batch_pnn.sh`
- **ilastik** — required upstream of `05_tissue_masker.py`
- **QuickNII / VisuAlign** — atlas registration; outputs consumed by `02`, `08`, `09`


---

## Usage

1. Set the path variables at the top of each script. All scripts use a `BASE`, `DATA_PATH`, or equivalent configuration block — no paths are embedded elsewhere.

2. Run the preprocessing pipeline in order (`00a` → `12`), once per animal.

3. Run statistical analysis scripts as needed. Each header documents which upstream outputs it requires.

4. Run `collect_all_results.py` to produce `ALL_RESULTS_COMBINED.csv`:
   ```bash
   python collect_all_results.py
   ```

5. Generate figures:
   ```bash
   python thesis_figures.py
   python fig_zone_overlay.py
   python fig_enwrapment_comparison_40x.py
   python fig_virus_proximal_40x.py
   ```

---

## Data availability

Raw imaging data and per-cell quantification tables are not included due to file size. The aggregated result tables required to reproduce all statistical outputs and figures are available on request.

---

## Citation

If you use or adapt these scripts, please cite the associated thesis:

> [Mads Lønø]. *[Spatial and Cell-Local Mechanisms of Perineuronal Net Remodelling in Adult Visual Cortex
]*. Master's thesis, University of Oslo, 2026.

