# 004 — TISSUE tutorial on the built-in test data (Phase 1, note 4 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Paper:** Sun, E.D., Ma, R., Navarro Negredo, P. et al. *TISSUE: uncertainty-calibrated prediction of single-cell spatial transcriptomics improves downstream analyses.* Nat Methods (2024). doi:10.1038/s41592-024-02184-y
**Result:** **The official example runs and passes.** `test.py` in **13 s at 420 MB**; the full tutorial path in **11 s at 277 MB**. Nothing came close to the 500 MB download or 10 GB RAM limits. No CUDA, no Hyprland/theme/system changes; everything under `~/Work/paper2agent`.

## Short version

TISSUE ships its own tiny dataset (**420 KB**, in the repo), so there was nothing to download beyond an 18.7 MB shallow clone. The whole pipeline — predict a held-out gene, calibrate uncertainty, get prediction intervals, run multiple-imputation testing, filter uncertain cells — finishes in seconds on CPU. **The central claim checks out on this data: calibrated 77 % prediction intervals achieved 78.8 % empirical coverage over the calibration genes.**

The catch is the environment, not the compute: TISSUE is pinned to a **2022-era stack** and will not run on a current one.

## Inputs and limits

| | |
|---|---|
| Repo | `https://github.com/sunericd/TISSUE`, shallow clone (`--depth 1`), commit `ffc3599` (2024-03-10, "updated tutorial"), MIT |
| Clone size | **14 MB** on disk (GitHub reports 18.7 MB) at `TISSUE_Agent/repo/TISSUE` |
| Data | **The repo's own `tests/data/`, 420 KB total** — `Spatial_count.txt` (235 KB), `scRNA_count.txt` (94 KB), `Locations.txt` (82 KB). A subset of a dataset from the paper. **Nothing downloaded**; the large spatial dumps (`tissue-figures-and-analyses`) were not touched. |
| Dimensions | Spatial 3,405 cells × 32 genes; scRNA-seq 1,000 cells × 32 genes; 32 shared genes |
| Limits | CPU only, every run under `tools/bin/cap10g`; peak **420 MB** (cap 10 GB, 0 OOM kills) |
| Downloads | 14 MB clone + Python packages. **No data file over 500 MB; no data file at all.** |

## What ran

### 1. The official self-test — passed

`python test.py` from the repo root is the smallest official example. It exercises data loading, preprocessing, SpaGE prediction, spatial graph, conformal calibration, prediction intervals, multiple-imputation t-test and cell filtering:

```
Testing TISSUE data loading...
Testing TISSUE preprocessing...
Testing TISSUE spatial gene expression prediction...
Testing TISSUE calibration...
Testing TISSUE multiple imputation t-test...
Testing TISSUE cell filtering
TISSUE tests passed!
```

**13 s, peak 420 MB.** Note `test.py` only asserts "no exception" — it prints no numbers and checks no values.

### 2. The README tutorial path, with the numbers recorded

Because `test.py` proves nothing quantitative, I re-ran the same official calls (README Tutorials 1–3) in `TISSUE_Agent/run_tutorial.py` and recorded what a user would actually read. **11 s, peak 277 MB.**

| Stage | Result |
|---|---|
| Load + preprocess | 3,405 × 32 spatial, 1,000 × 32 scRNA-seq, 32 shared genes |
| SpaGE prediction of held-out `plp1` (3 folds, 10 PVs) | Pearson **r = 0.310**, Spearman **0.185** vs the measured values |
| Conformal calibration, `alpha_level=0.23` (77 % intervals) | **Coverage over the 31 calibration genes: 0.788** (nominal 0.77); mean interval width 12.6 |
| Same intervals, held-out `plp1` only | Coverage **0.854**, mean width 7.55 |
| Multiple-imputation t-test, A vs B (10 imputations) | t = −1.611, **p = 0.117**, pooled mean diff −0.113 |
| Cell filtering (`proportion="otsu"`) | keeps **2,938 / 3,405** cells (13.7 % dropped); `filtered_PCA` keeps the same 2,938 |

The A/B split is the tutorial's own arbitrary "first half vs second half of cells", so **p = 0.117 is a mechanism check, not a biological result** — there is no real condition in this data.

The **0.788 vs 0.770** coverage is the result worth keeping: that is TISSUE's core claim reproducing on its own tiny dataset. The held-out gene over-covers (0.854), which is ordinary single-gene variation around an aggregate guarantee.

Artifacts: `TISSUE_Agent/out/tutorial_results.json` and `plp1_prediction_intervals.png`.

## What failed

### The current scientific Python stack — TISSUE will not run on it

First attempt was Python 3.12 with current packages (squidpy 1.8.3, scanpy 1.12.4, numpy 2.5.3, anndata 0.13.4). It died on the first call:

```
TypeError: AnnData.__init__() got an unexpected keyword argument 'dtype'
```

`tissue/main.py:73` calls `ad.AnnData(X=df, dtype='float64')`; `dtype=` was deprecated in anndata 0.10 and **removed** in 0.11. That 1.2 GB environment was useless and I deleted it.

### The pinned stack needs one fix the repo does not mention

`requirements.txt` pins `squidpy==1.2.3` and the README says Python 3.8. On Python 3.10 that resolves to scanpy 1.9.8 / numpy 1.22.4 / anndata 0.10.5, but the run still failed:

```
ModuleNotFoundError: No module named 'pkg_resources'
```

The era-appropriate `numba` imports `pkg_resources`, which modern virtualenvs no longer ship. **Fix: `pip install "setuptools<81"`.** That is an undocumented extra step for anyone installing TISSUE into a venv rather than the conda env the README assumes. After it, everything passed.

### A scale mismatch inside the official pipeline

Worth recording because it is not a mistake I made and the README does not mention it. SpaGE predictions and the measured counts are **not on the same scale**:

| | measured `plp1` | predicted `plp1` |
|---|---|---|
| mean | 2.577 | 0.519 |
| max | 87.0 | 1.543 |

TISSUE subtracts them anyway — `residuals = adata[:, calib_genes].X - adata.obsm[predicted][calib_genes].values` (`tissue/main.py:652`) — so the residuals absorb the offset and the prediction intervals come out wide (mean width 12.6 on data whose calibration genes are mostly small counts). Coverage still lands on target because the conformal step calibrates against those same inflated residuals. The README only ever compares predicted vs actual **visually**, after `log1p` and 95th-percentile clipping, which hides this. It does not invalidate the coverage result, but any per-cell interval is wider than the gene's own dynamic range suggests.

**I got this wrong first:** my initial run reported "coverage 0.854" as *the* result, computed only on the held-out gene, and separately miscounted kept cells as 5,422,156 by summing index values from `detect_uncertain_cells` (which returns **row indices**, not a boolean mask). Both are fixed above.

## RAM and disk

- **RAM:** peak **420 MB** (`test.py`), **277 MB** (tutorial driver). 0 OOM kills. This never needed anything close to 16 GB — it would run on a 2 GB machine.
- **Disk added:** **707 MB** total under `TISSUE_Agent/` —
  - `tissue-env-pinned/` **693 MB** (Python 3.10 + squidpy 1.2.3 and deps)
  - `repo/TISSUE` 14 MB, `out/` 80 KB, `tmp/` 568 KB
  - The failed Python 3.12 environment was another **1.2 GB**; **deleted**.
- `tools/` (shared uv cache + Python installs) grew to 2.4 GB. 204 GB still free.

## How to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh
export TISSUE_ROOT="$P2A/TISSUE_Agent"

# 1. Shallow clone (14 MB)
mkdir -p "$TISSUE_ROOT/repo"
git clone --depth 1 https://github.com/sunericd/TISSUE.git "$TISSUE_ROOT/repo/TISSUE"

# 2. Pinned environment — Python 3.10, NOT a current stack (693 MB)
uv venv --python 3.10 "$TISSUE_ROOT/tissue-env-pinned"
uv pip install --python "$TISSUE_ROOT/tissue-env-pinned/bin/python" "squidpy==1.2.3"
uv pip install --python "$TISSUE_ROOT/tissue-env-pinned/bin/python" "setuptools<81"   # numba needs pkg_resources

# 3. The official self-test (~13 s, ~420 MB)
cd "$TISSUE_ROOT/repo/TISSUE"
XDG_CACHE_HOME="$TISSUE_ROOT/tmp/cache" NUMBA_CACHE_DIR="$TISSUE_ROOT/tmp/cache/numba" \
MPLCONFIGDIR="$TISSUE_ROOT/tmp/mpl" MPLBACKEND=Agg \
  cap10g "$TISSUE_ROOT/tissue-env-pinned/bin/python" test.py

# 4. The tutorial path with recorded numbers (~11 s, ~277 MB)
PYTHONPATH=. TISSUE_OUT="$TISSUE_ROOT/out" \
XDG_CACHE_HOME="$TISSUE_ROOT/tmp/cache" NUMBA_CACHE_DIR="$TISSUE_ROOT/tmp/cache/numba" \
MPLCONFIGDIR="$TISSUE_ROOT/tmp/mpl" MPLBACKEND=Agg \
  cap10g "$TISSUE_ROOT/tissue-env-pinned/bin/python" "$TISSUE_ROOT/run_tutorial.py"

cat "$TISSUE_ROOT/out/tutorial_results.json"
```

`PYTHONPATH=.` is required in step 4 because the driver lives outside the repo root while `tissue/` must be importable from the cwd.

## Reusable like the Scanpy tools, or still a one-off script?

**In between, and closer to reusable than the Scanpy tutorial was.** The comparison that matters is note 002, where nine MCP tools had to be *extracted* from a notebook.

**Already library-shaped:**
- TISSUE is an installed package (`tissue.main`, `tissue.downstream`) with **~34 public functions**, not a notebook. No extraction step needed — the hosted `paper2agent-tissue-mcp` wraps 6 of them and note 001 confirmed that server answers.
- Functions take an AnnData plus explicit keyword arguments (`method`, `n_folds`, `alpha_level`, `grouping_method`, `k`, `proportion`), which is exactly the switch surface the Scanpy tools ended up exposing.
- Every function has a real docstring with a Parameters/Returns block — enough to generate tool schemas from.

**What still makes it one-off:**
- **Everything mutates AnnData in place under string keys** (`spage_predicted_expression`, `..._lo`, `..._hi`, `spage_A_B_pvalue`). The method name is baked into the key, so a caller must know the convention to find results. The Scanpy tools solved this by writing a fresh `.h5ad` per call and returning artifact paths; TISSUE would need the same wrapper.
- **Return types are inconsistent and undocumented at the edges.** `detect_uncertain_cells` returns row indices where a boolean mask is the natural guess — the mistake I made above. A tool wrapper must normalise these.
- **No validation.** Pass a gene that isn't there, or an unbuilt spatial graph, and you get a numpy/pandas error from deep inside, not a clear message. The Scanpy tools' best property in note 003 was that every cross-dataset failure came back as a clean `ValueError` naming the missing column. TISSUE has nothing equivalent.
- **The pinned-stack problem is the real blocker.** A Paper2Agent conversion would have to ship the Python 3.10 + squidpy 1.2.3 + `setuptools<81` environment, and could not share the Scanpy agent's environment.
- **`test.py` asserts nothing numeric**, so there are no reference values to verify a wrapper against. Note 002 got bitwise-identical references from a real notebook run; here I'd have to establish them first — the numbers in this note are a starting point.

**Verdict:** converting TISSUE would be *less* work than Scanpy (no notebook extraction, real function signatures, a hosted 6-tool server to compare against) but needs its own pinned environment and a validation layer the package does not provide. The tutorial itself is firmly a one-off script; the library underneath is not.

## Caveats (honest)

- Only the **SpaGE** prediction method was exercised. Tangram, Harmony-kNN, gimVI and SpatialDE are commented out in `requirements.txt` and were not installed or run.
- `n_folds=3` (from `test.py`), not the README's `n_folds=10` — faster, and it is the repo's own test setting.
- One dataset, one target gene (`plp1`), 32 genes. r = 0.310 is modest; the README itself calls the result "not too bad, especially considering that we used a downsampled scRNAseq dataset".
- The coverage numbers come from my driver, not from TISSUE's own output — TISSUE reports intervals, not coverage. The arithmetic is a direct comparison of measured values against `_lo`/`_hi`.
- Nothing here was converted into MCP tools. No agent was built; this was a tutorial run only.

## Next

Phase 1 note 5 of 10. Still open: Stages 5–6 for `Scanpy_Agent` (clean-environment install, `USAGE.md`, ZIP) from note 002. If TISSUE is converted next, the environment pin and a validation layer are the two things to budget for.
