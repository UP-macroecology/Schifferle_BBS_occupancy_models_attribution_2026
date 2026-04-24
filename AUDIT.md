# Repository notes — coauthor read-through

**From:** Guillermo Fandos
**To:** Katrin
**Date:** 2026-04-15
**Branch:** `guille-feedback`

Hi Katrin — as you asked, I sat down with the repo and the README with a reviewer's hat on, focused on readability and reproducibility rather than the science. Nothing here is a blocker; these are just suggestions for tightening things up before submission.

Alongside this memo, there's a `README_DRAFT.md` in the same branch with three self-contained additions that could slot into the current README (a workflow diagram, a "Reproducing the analysis" section, and a header/citation block). It's additive, not a replacement — take any subset, or none.

Overall the repository is in good shape. The script numbering is clear, the dependency chain between scripts is linear and easy to follow, the README narrative tracks the paper's Methods closely, and the separation between data prep / fitting / evaluation / attribution is clean. Most of the points below are about smoothing friction for an external reader (a reviewer, a new collaborator, or a future user pulling the code from Zenodo) rather than anything being wrong scientifically. They're ranked roughly by how much friction they'd cause someone opening the repo cold.

---

## Section 1 — Path portability

### 1.1 Hard-coded NAS paths across the fitting and attribution scripts

Most scripts from the fitting step onwards (the 2_*, 3_*, 4_*, and 5_* series, plus some of the data-prep scripts) reference `//NAS-2-P-SN-01.ibb.uni-potsdam.de/daten$/...`. These paths only resolve inside the IBB Potsdam network, so anyone outside the institute would need to edit each script by hand before running anything.

One low-effort fix would be a single `config.R` (or `paths.R`) at the repo root that defines `dir_data`, `dir_models`, `dir_results` based on `Sys.info()` or an environment variable, with each script starting with `source("config.R")`. The existing logic stays untouched.

### 1.2 Mixed Windows-local paths (`T:/`, `C:/`) in several scripts

`1_2b` and most of the 2_* scripts (including `2_3a`, `2_3b`, `2_4c`) plus `4_0` contain hard-coded Windows paths. In the DOM-fitting and prediction scripts, most of these are already commented out as a local fallback after the HPC-oriented `set_cmdstan_path(path = NULL)` — that pattern works. The two that would actually break on a fresh machine are `2_3a` (line 23) and `2_4c` (line 19), where `cmdstanr::set_cmdstan_path("C:/Users/schifferle1/...")` is live, not commented.

Easiest fix: in `2_3a` and `2_4c`, flip to the same pattern you already use in `2_1`, `2_2`, `2_4b` — `set_cmdstan_path(path = NULL)` with the Windows path moved to a commented local fallback.

### 1.3 SLURM shell scripts (`1_2c_attrici_*.sh`) with cluster-specific paths

These point to `/mnt/ibb_share/...` and assume a specific SLURM environment. Without a header explaining "these jobs are written for cluster X, paths must be adapted before submission", someone coming to the repo fresh won't know what to do with them.

A short header in each `.sh` (plus maybe a `scripts/shell/README.md` if they get grouped) explaining the cluster context would fix it.

### 1.4 Absolute GitHub links in the README

Every script reference in the README uses `https://github.com/UP-macroecology/Schifferle_BBS_occupancy_models_2023/blob/master/scripts/...`. If the repo is forked, moved, or archived to Zenodo, all of those links break, and locally-rendered README links leave the repo.

A single find-and-replace to relative links (`[script](scripts/X_Y.R)`) solves it in ten minutes. Flagged in `README_DRAFT.md`.

---

## Section 2 — Clarity for an external reader

### 2.1 The README explains the *what* but not the *how to run it*

The current README is a faithful abridged version of the paper's Methods. It tells the reader what each script does scientifically, but doesn't answer the practical questions someone has on first contact with the repo:

- Which datasets must be downloaded before starting? (BBS, ISIMIP climate, ISIMIP land use, BirdLife ranges, BCR shapefile)
- Where do they go? What folder structure does the code expect?
- In what order are the scripts run, given the branching points (factual/counterfactual, `round <- 2` reruns)?
- What software is required? (R 4.3.1, cmdstan, ATTRICI 1.1.1, package versions, GPU, cluster)
- How long does it take? Which steps require HPC?

A short "Reproducing the analysis" section covers this — sketched in Block 3 of `README_DRAFT.md`.

### 2.2 Workflow branch points are not signposted

A handful of places require rerunning a script with a different parameter, and those branch points are only mentioned in prose:

- `1_3_dataprep_match_BBS_routes_env_data.R` is run twice (`data <- "factual"` and `data <- "counterfactual"`).
- `2_1`, `2_2`, `2_4b` are rerun with `round <- 2` for non-converged species.
- `2_3a` and `2_3b` are rerun "with adjusted file paths" (paths not specified).

A small mermaid flow diagram at the top of the Workflow section makes this obvious in seconds. Sketched as Block 2 in `README_DRAFT.md`.

### 2.3 Empty header fields in the README

Lines 5 and 11 of the current README contain literal placeholders ("Authors, institutions, funding: ..." and an empty `### Abstract`) — worth filling in before submission. `README_DRAFT.md` keeps the placeholders explicit so they are not forgotten.

### 2.4 Script headers: present but inconsistent

Most scripts open with a header comment stating the purpose, but the format varies: some include author/date, some don't; some list inputs/outputs, some don't; `2_3b.qmd` has no header at all. For a reader opening any single script in isolation, a common template applied gradually would help — something like:

```r
# Script: 2_1_fit_DOMs_full_model.R
# Purpose: Fit single-species DOMs with flocker for 192 selected species
# Inputs:  data/processed/route_year_env_data_factual.RData
#          data/processed/species_selection_final.RData
# Outputs: results/DOMs_full/<species>.rds (one per species)
# Runs on: HPC (NAS Potsdam); ~X CPU-h per species
# Author / Date: KS / 2024
```

---

## Section 3 — Technical reproducibility

### 3.1 Missing `set.seed()` in stochastic steps

Three points where the result depends on randomness without an explicit seed:

- **`2_4a_fit_DOMs_CV_fold_assignment.R`** — `blockCV::spatialBlock()` assigns folds randomly. Without a seed, the folds change on each run, so the CV metrics in the paper aren't bit-reproducible. This is the one a methodology reviewer is most likely to pick up.
- **`1_1_dataprep_BBS_route_selection.R`** — spatial thinning has an internal per-iteration seed, but the wrapper uses `furrr::future_map`; with `runs > 1` the result depends on scheduling order.
- **MCMC in `2_1`, `2_2`, `2_4b`** — brms/cmdstan accepts `seed = ...`. Not clear from the scripts whether it's being passed; if not, posterior predictive checks and exact diagnostics shift between runs (the posterior itself converges to the same target).

The fix is small: declare a seed at the top of each stochastic script and pass it explicitly to `spatialBlock(seed = ...)` and `brm(seed = ...)`. One line per script.

### 3.2 No `sessionInfo()` and no package version management

No script ends with `sessionInfo()`. The package list lives in the README as a static text block (R 4.3.1, brms 2.21.0, flocker 1.0-0, cmdstanr 0.7.1, ATTRICI 1.1.1...) — useful, but manual and easy to drift.

A few options in order of effort:

- (a) Add `sessionInfo()` at the end of each script (writing to `results/sessionInfo/<script>.txt`).
- (b) Adopt `renv` for a fully reproducible snapshot. Anyone can `renv::restore()` and obtain the exact versions.
- (c) A `DESCRIPTION` file or `dependencies.R` with `install.packages()` calls pinned to versions.

Option (b) is the strongest answer to a reviewer; (c) is a defensible minimum.

### 3.3 External data: incomplete inventory and no download instructions

The README cites the data sources (BBS, ISIMIP GSWP3-W5E5, ISIMIP3a landuse, BirdLife, BCR shapefile) but doesn't spell out:

- Which exact files to download (ISIMIP has hundreds of variables — which subset?).
- Approximate sizes.
- Which sources require registration or formal access requests (BirdLife needs an application; BBS needs terms acceptance).
- Where the files are expected to sit locally.

For a reviewer trying to rerun the pipeline, "downloaded the data from ISIMIP" isn't a reproducible specification.

A table in the README (`Dataset | Source | Files | Size | Access`) covers most of this. Sketched in Block 3 of `README_DRAFT.md`, with placeholders where I didn't want to guess the actual sizes or folder paths — those are yours to fill in.

### 3.4 ATTRICI installation and Zenodo DOI

Two loose ends worth flagging before submission:

- **ATTRICI 1.1.1** is named but not documented for install. It's a Python CLI with its own dependencies, and the SLURM `.sh` jobs assume it's already installed.
- **Zenodo DOI:** required by GCB. Currently absent from the README — which is expected at this stage, but worth flagging here as a TODO so it doesn't slip through.

A short "Setting up ATTRICI" paragraph (install command, upstream link, exact version) plus an explicit DOI placeholder as a reminder — both sketched in Block 1 and Block 3 of `README_DRAFT.md`.

---

## Summary

None of this is a blocker. Rough priorities if you want to address any of it before submission:

| Section | Effort | Priority |
|---|---|---|
| 1. Path portability | Low–medium | High |
| 2. Reader clarity | Low | Medium |
| 3. Reproducibility | Low–medium | High |

`README_DRAFT.md` in this branch covers the README-level points (1.4, 2.1, 2.2, 2.3, 3.3) and flags 3.4 as reminders. The rest (1.1–1.3, 2.4, 3.1, 3.2) need small edits in the scripts themselves and are yours to decide on.

Happy to help with any of it if useful — just let me know.

— Guille
