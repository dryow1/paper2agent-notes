# 005 — Wrapping the smallest useful slice of TISSUE (Phase 1, note 5 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Builds on:** [004](004-tissue-tutorial.md) — same repo checkout, **same `tissue-env-pinned` environment, no new venv**
**Result:** **Two functions wrap the tutorial path and reproduce its numbers exactly.** 3/3 calls succeeded, **7/7 wrong-input cases returned a clean `ValueError`**. 14 s, peak 407 MB.

## Short version

Note 004 concluded TISSUE was "closer to reusable than Scanpy, but needs a validation layer the package does not provide." This note builds that layer for the smallest slice that is still useful — **predict a held-out gene, then calibrate the uncertainty** — and it took **two functions, ~300 lines, no rewrite**. The wrappers add no maths: every scientific call is upstream TISSUE. Coverage came out at **0.7883**, identical to note 004's direct run.

## What was wrapped

`TISSUE_Agent/src/tissue_tools.py`, two functions, shaped like the Scanpy tools in `Scanpy_Agent/src/tools/clustering.py` (keyword arguments, fresh output directory per call, an `.h5ad` written so the two chain, a dict of absolute artifact paths plus a numeric summary).

### 1. `tissue_predict_spatial_gene(...)`

```python
tissue_predict_spatial_gene(
    spatial_counts_path, locations_path, scrna_counts_path, *,
    target_gene="plp1", method="spage", n_folds=3, n_pv=10, output_dir=None,
) -> dict
```

Wraps `load_paired_datasets` → `preprocess_data` → `predict_gene_expression`. Holds the target gene out of the spatial matrix and **keeps the measured values in `obs['<gene>_measured']`**, which is what lets the next step score coverage without re-loading the raw files.

### 2. `tissue_calibrate_prediction_intervals(...)`

```python
tissue_calibrate_prediction_intervals(
    data_path, *, target_gene=None, method=None, alpha_level=0.23,
    n_neighbors=15, graph_method="fixed_radius",
    grouping_method="kmeans_gene_cell", k=4, k2=2, output_dir=None,
) -> dict
```

Wraps `build_spatial_graph` → `conformalize_spatial_uncertainty` → `conformalize_prediction_interval`. Reports **coverage over the calibration genes** (the quantity TISSUE actually targets) and, because step 1 stashed the measured column, coverage for the held-out gene too, plus a figure. `target_gene` and `method` default to whatever step 1 recorded in `uns['tissue_tool']`, so the common case needs only `data_path`.

**Not wrapped, deliberately:** the other ~32 functions — multiple-imputation testing, cell filtering, weighted PCA, Tangram/gimVI/SpatialDE paths. Note 004 exercised MI testing and cell filtering; they are not part of this slice.

## What ran

Driver: `TISSUE_Agent/src/run_tissue_tools.py`, on the repo's built-in **420 KB** dataset (3,405 cells × 32 genes). **Nothing downloaded.**

| Call | Result |
|---|---|
| `predict` (plp1, spage, 3 folds, 10 PVs) | ✅ 6.9 s — 31 calibration genes, **Pearson r = 0.3097** vs measured |
| `calibrate` (alpha 0.23 → 77 % intervals) | ✅ 3.3 s — **coverage 0.7883** over 31 calibration genes (nominal 0.77); held-out gene 0.8537; mean widths 12.59 / 7.55 |
| `calibrate` (alpha 0.5 → 50 % intervals) | ✅ 3.2 s — **coverage 0.5397** (nominal 0.50); held-out gene 0.6135; mean widths 4.13 / 2.62 |

**Every number matches note 004's direct tutorial run exactly** (0.7883, 0.8537, 12.5879, 7.5545, r = 0.3097). That is the check that matters: the wrapper changed nothing.

The `alpha_level` switch does what it claims — halving the nominal coverage cut the mean interval width from 12.59 to 4.13, and the plotted band is visibly narrower.

### Error handling — the point of the exercise

All seven wrong-input cases returned a `ValueError` before any computation started:

| Case | Message |
|---|---|
| missing spatial file | `spatial_counts_path does not exist or is not a file: .../nope.txt` |
| gene not in data | `target_gene 'notagene' not found in the spatial dataset (33 genes). Gene names are compared in lowercase.` |
| unsupported method | `method 'tangram' is not available in this environment; supported: ['spage', 'knn']. tangram/gimvi need extra packages that requirements.txt leaves commented out.` |
| `n_folds` too large | `n_folds (999) cannot exceed the 31 calibration genes left after holding out 'plp1'` |
| calibrate on a non-`.h5ad` | `data_path must be one of ('.h5ad',), got: .../Spatial_count.txt` |
| alpha out of range | `alpha_level must be strictly between 0 and 1, got 1.5` |
| wrong method key | `obsm['knn_predicted_expression'] missing; run tissue_predict_spatial_gene first (available obsm keys: ['spage_predicted_expression', 'spatial'])` |

Without the wrapper each of these surfaces as a `KeyError`, `IndexError` or a numpy broadcast error from inside TISSUE, or — worse — not at all: passing `n_folds=999` to raw `predict_gene_expression` fails deep in cross-validation, and asking for a gene that isn't there fails during the intersection with no mention of the gene.

## What failed

**Nothing failed in this run** (3/3 calls, 7/7 error cases as intended). Two honest qualifications:

1. **One error case tested a different guard than intended.** "Calibrate on an un-predicted `.h5ad`" was meant to hit the missing-`obsm` check, but the file I passed was a `.txt`, so the suffix guard caught it first. The missing-`obsm` path is covered by the last case (`wrong method key`) instead — so both guards are exercised, just not by the cases I had labelled for them.
2. **The scale mismatch from note 004 is inherited, not fixed.** SpaGE predictions (mean 0.52) and measured counts (mean 2.58, max 87) are on different scales; TISSUE subtracts them anyway (`tissue/main.py:652`). The wrapper surfaces both means in its return dict so a caller can see it, but it does not correct it — that would be inventing maths the paper does not do.

## RAM and disk

- **RAM:** peak **407 MB** in-process, **289 MB** as measured by `cap10g`. 0 OOM kills, nowhere near the 10 GB cap.
- **Disk:** **no new environment** — reused `tissue-env-pinned` (693 MB) from note 004. New bytes: **25 KB of source**, plus 23 MB of run outputs in `TISSUE_Agent/out/` (two 10 MB `calibrated.h5ad` files and a 2 MB `predicted.h5ad`; the `.h5ad` carries the full spatial graph, which is why it is 5× the input).
- `TISSUE_Agent/` total is 729 MB, of which 693 MB is the venv. Delete `TISSUE_Agent/out/` to reclaim 23 MB.

## How to call it again

```bash
cd ~/Work/paper2agent && . ./env.sh
export TISSUE_ROOT="$P2A/TISSUE_Agent"

# Both tools plus the seven error cases (~14 s, ~400 MB)
PYTHONPATH="$TISSUE_ROOT/src" TISSUE_OUT="$TISSUE_ROOT/out" \
XDG_CACHE_HOME="$TISSUE_ROOT/tmp/cache" NUMBA_CACHE_DIR="$TISSUE_ROOT/tmp/cache/numba" \
MPLCONFIGDIR="$TISSUE_ROOT/tmp/mpl" MPLBACKEND=Agg \
  cap10g "$TISSUE_ROOT/tissue-env-pinned/bin/python" "$TISSUE_ROOT/src/run_tissue_tools.py"

cat "$TISSUE_ROOT/out/tools_results.json"
```

From Python:

```python
import sys; sys.path.insert(0, "<TISSUE_Agent>/src")
from tissue_tools import tissue_predict_spatial_gene, tissue_calibrate_prediction_intervals

pred = tissue_predict_spatial_gene(
    spatial_counts_path="<repo>/tests/data/Spatial_count.txt",
    locations_path="<repo>/tests/data/Locations.txt",
    scrna_counts_path="<repo>/tests/data/scRNA_count.txt",
    target_gene="plp1", method="spage", n_folds=3, output_dir="out",
)
cal = tissue_calibrate_prediction_intervals(
    data_path=pred["artifacts"][0]["path"], alpha_level=0.23, output_dir="out",
)
print(cal["coverage_calibration_genes"])   # 0.7883
```

The module finds the TISSUE checkout itself (`$TISSUE_REPO`, else the sibling `repo/TISSUE`) and puts it on `sys.path`, so no `PYTHONPATH=.` from the repo root and no `cd` into it — which is what note 004's driver needed.

## Was a rewrite needed? No

The instruction was to stop if wrapping demanded a rewrite big enough to be a new project. It did not:

- **~300 lines for two functions**, of which roughly half is validation and result formatting. Nine upstream TISSUE calls, unchanged.
- **The `.h5ad` chain works out of the box.** I checked first: `obsm['spage_predicted_expression']` survives a write/read round-trip as a `DataFrame` with column names intact, and so do `obsp['spatial_connectivities']` and `uns`. That was the one thing that could have forced a redesign (a custom serialisation format), and it didn't.
- **In-place mutation under method-derived string keys** — note 004's main complaint — turned out to be containable: the wrapper records `method` in `uns['tissue_tool']` and derives the keys itself, so callers never type `"spage_predicted_expression_lo"`.

What a fuller conversion would still need, unchanged from note 004: reference values to test against (none exist upstream — `test.py` asserts nothing numeric; the numbers here are the start of a set), the pinned Python 3.10 environment shipped alongside, and the same treatment for the downstream functions, whose return types are inconsistent (`detect_uncertain_cells` returns row indices where a boolean mask is the natural guess — the bug I hit in note 004).

## Caveats (honest)

- **One dataset, one gene, one method.** Only `spage` was run; `knn` is allowed by the validation but untested. The wrappers have never seen data other than the built-in 420 KB set.
- **No test suite.** Correctness rests on the numbers matching note 004's independent run — not on independent verification of TISSUE itself.
- **Not MCP tools.** These are plain Python functions. No `fastmcp` decoration, no server, nothing registered — unlike `Scanpy_Agent`.
- The seven error cases are the ones I thought to write; they are not an exhaustive audit of the input space.
- `n_folds=3` throughout (the repo's own `test.py` setting), not the README's 10.

## Next

Phase 1 note 6 of 10. Still open: Stages 5–6 for `Scanpy_Agent` (clean-environment install, `USAGE.md`, ZIP) from note 002. If TISSUE goes further, the next steps are a pytest suite pinning these numbers, then `fastmcp` decoration to match the Scanpy server.
