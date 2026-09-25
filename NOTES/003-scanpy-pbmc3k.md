# 003 — Reusing the nine Scanpy tools on pbmc3k (Phase 1, note 3 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU used**
**Project:** `~/Work/paper2agent/Scanpy_Agent/` (existing tools from note 002; no new paper, no new clone, no new env)
**Result:** **All nine tools reused cleanly on pbmc3k.** 12 of 15 calls succeeded; the 3 failures were deliberate unadapted calls kept to record exactly what pbmc3k is missing.

## Short version

The nine MCP tools built from the bone-marrow tutorial ran end to end on Scanpy's pbmc3k with **no code change**. Three arguments had to be dropped or adapted, all for the same reason: **pbmc3k has no `sample` column and an older, smaller gene annotation**. The whole run took **37 seconds** and peaked at **1.15 GB**, against ~6 GB and hours for the bone-marrow run.

## Inputs and limits

| | |
|---|---|
| Data | `sc.datasets.pbmc3k()` → `pbmc3k_raw.h5ad`, **5.6 MB**, from `https://exampledata.scverse.org/scanpy/pbmc3k_raw.h5ad` |
| Where it landed | `~/Work/paper2agent/data/pbmc3k_raw.h5ad` — scanpy's `settings.datasetdir` defaults to `./data` under the cwd and **ignores `XDG_CACHE_HOME`** from `project.env`. Still inside the project; nothing was written to `~/.cache`. |
| Size | 2,700 cells × 32,738 genes raw; 2,700 × 13,714 after QC (no cell dropped at `min_genes=100`) |
| Limits | CPU only, every command under `tools/bin/cap10g`; peak **1,153 MB** (cap 10 GB, 0 OOM kills); download 5.6 MB (limit 500 MB) |
| Environment | Unchanged: `Scanpy_Agent/scanpy-env`, Python 3.12.14, scanpy `1.13.0a3.dev0+g0d5fd1623`. No install, no CUDA, no rapids, no second env. |
| Not touched | Hyprland, theme, system config, `.pipeline/` markers, `notebooks/`, `tests/`, `src/` |

## What ran

Driver: `Scanpy_Agent/tmp/pbmc3k/run_pbmc3k.py`. It calls the tools through an **in-process fastmcp `Client`** — the same path the verifier's tests use — so argument validation and the JSON round-trip are exercised, not just the Python functions.

| # | Tool | Call | Result |
|---|---|---|---|
| 1 | `scanpy_compute_qc_and_filter` | as tutorial | ✅ 1.7 s — kept 2,700/2,700 cells, 13,714/32,738 genes; 13 mt, 99 ribo, 4 hb genes |
| 2a | `scanpy_detect_doublets_scrublet` | `batch_key="sample"` | ⛔ **expected failure** (see below) |
| 2b | `scanpy_detect_doublets_scrublet` | `batch_key=None` | ✅ 6.3 s — **32 doublets** / 2,700 cells (1.19 %) |
| 3 | `scanpy_normalize_log1p` | as tutorial | ✅ 0.7 s |
| 4 | `scanpy_select_highly_variable_genes` | `batch_key` dropped | ✅ 1.0 s — 2,000 HVGs, flavor `seurat` |
| 5 | `scanpy_run_pca` | `color_keys=["pct_counts_mt"]` | ✅ 1.6 s — 50 PCs, variance ratio first 5 = 0.1457 / 0.0473 / 0.0321 / 0.0126 / 0.0106 |
| 6 | `scanpy_build_neighbors_umap` | `color_keys=["pct_counts_mt"]` | ✅ 11.4 s — slowest step |
| 7a | `scanpy_cluster_leiden` | res 1.0, `n_iterations=2`, QC panel | ✅ 1.3 s — **10 clusters** |
| 7b | `scanpy_cluster_leiden` | res 0.02 / 0.5 / 2.0, `n_iterations=-1` | ✅ 1.2 s — **3 / 8 / 15 clusters** |
| 8a | `scanpy_plot_marker_gene_dotplot` | tutorial marker dict | ⛔ **expected failure** (see below) |
| 8b | `scanpy_plot_marker_gene_dotplot` | 8 genes dropped | ✅ 0.5 s — 15 marker sets × 3 groups |
| 8c | `scanpy_plot_marker_gene_dotplot` | tutorial `cluster_labels` | ⛔ **expected failure** (see below) |
| 8d | `scanpy_plot_marker_gene_dotplot` | placeholder labels over real categories | ✅ 1.1 s — `cell_type_lvl1` written, 0 unlabeled cells |
| 9 | `scanpy_rank_marker_genes` | `groupby="leiden_res_0.50"`, `focus_group="7"` | ✅ 3.1 s — 8 groups, **109,712 CSV rows** |

Every tool produced its figures and they are real plots, not blanks (no PNG under 2 KB; the UMAP panel was inspected by eye and shows the expected PBMC structure). **Determinism:** step 9 was run twice with identical arguments and the two `rank_genes_groups.csv` files are **byte-identical**.

## What failed versus the bone-marrow run

All three failures are the **same root cause in two forms** — pbmc3k is one sample with an older gene annotation — and each is a clean, explicit `ValueError` from the tool's own precondition checks, not a crash deep in scanpy. That is the tools behaving well.

1. **No `sample` column.** `batch_key="sample"` →
   `batch_key not found in adata.obs: ['sample']; available: [...]`
   pbmc3k arrives with **zero obs columns**. Adapted once: `batch_key=None` everywhere (steps 2, 4). Consequence: no per-batch doublet breakdown and no batch-aware HVG selection — the bone-marrow run's `predicted_doublets_per_batch` (207 in s1d1, 6 in s1d3) has no counterpart here.
2. **Eight marker genes absent.** The tutorial's marker dict →
   `marker genes not found in var_names: ['GYPA', 'HBB', 'HBM', 'IGHD', 'IGHM', 'IGKC', 'JCHAIN', 'TRBC2']`
   Five (`IGHD IGHM IGKC JCHAIN TRBC2`) are missing from pbmc3k's 2015-era 10x hg19 annotation outright; three (`GYPA HBB HBM`) survive the annotation but are dropped by `filter_genes(min_cells=3)` — PBMCs have essentially no erythroid cells. Adapted once: drop those 8 genes. **No marker set was emptied**, so all 15 sets still plot. Saved at `tmp/pbmc3k/results/marker_genes_pbmc3k.json`.
3. **Tutorial's cluster-label map doesn't fit.** →
   `cluster_labels keys not among obs['leiden_res_0.02'] categories ['0', '1', '2']: ['3', '4']`
   The bone-marrow run had 5 clusters at res 0.02; pbmc3k has 3. Adapted once with **placeholder labels** (`cluster_0/1/2`) purely to exercise the annotation code path. **These are not cell-type calls** — I did not assign PBMC identities, because doing so would be inventing biology the run didn't establish.

Nothing else differed. No tool needed an edit, and no tool failed for a reason other than a missing input column or gene.

### Scale comparison

| | bone marrow (s1d1+s1d3) | pbmc3k |
|---|---|---|
| Cells × genes (after QC) | 17,041 × 23,427 | 2,700 × 13,714 |
| Download | 43 MB | 5.6 MB |
| Peak RAM | ~6,086 MB (verifier) / 4,991 MB (notebook) | **1,153 MB** |
| Pipeline wall time | ~60 min reference execution | **37 s** |
| Doublets | 213 (1.25 %) | 32 (1.19 %) |
| Leiden res 1.0 | 26 clusters | 10 clusters |
| Leiden 0.02 / 0.5 / 2.0 | 5 / 17 / 36 | 3 / 8 / 15 |
| Rank-genes CSV rows | 398,259 | 109,712 |

The clusters look like PBMCs, which is the sanity check that matters: top Wilcoxon genes are `LYZ/S100A9/S100A8` (monocytes, group 3), `CD74/HLA-DRA/CD79A` (B, group 1), `NKG7/GNLY/GZMB` (NK, group 4), `FCGR3A/LST1` (CD16 mono, group 5), and group 7 — 12 cells — is `GNG11/SDPR/PF4/PPBP`, the platelets pbmc3k is known for.

## RAM and disk

- **RAM:** peak 1,153 MB under `cap10g`, 0 OOM kills. Well under the 10 GB cap; the bone-marrow run needed ~5–6× more.
- **Disk added by this run:** **165 MB** total — 159 MB of chained `.h5ad` outputs, 5.6 MB download, ~76 KB of JSON records. The `.h5ad` files dominate because every tool writes a full copy; deleting `Scanpy_Agent/tmp/pbmc3k/outputs/` reclaims all 159 MB.
- Project total is now ~9.4 GB, still mostly the older `notebooks/` (2.3 GB), `tests/` (2.5 GB) and `tmp/` (4.0 GB) from note 002.

## Did the switches reuse cleanly?

**Yes.** The parameters that were meant to be switches behaved as switches:

- `batch_key` — defaulting to `None` (a note-002 design decision that differed from the tutorial) is exactly what made the single-sample dataset work with no edit. Good call in hindsight.
- `marker_genes`, `groupby`, `cluster_labels`, `annotation_key`, `focus_group` — all user-supplied, all accepted pbmc3k-appropriate values.
- `resolutions` + `n_iterations` — both tutorial patterns (single res 1.0 with `n_iterations=2`, sweep with `-1`) worked unchanged.
- `mt_prefix` stayed `"MT-"` (pbmc3k is human); `ribo_prefixes`, `hb_pattern`, `min_genes`, `min_cells`, `n_top_genes`, `n_comps`, `n_neighbors` all at tutorial defaults.

The only thing the tools **assume** and pbmc3k can't provide is a batch column, and the tools ask for it by name rather than requiring it. The failure messages list the available columns, which is what made each adaptation a one-line fix.

**One rough edge worth fixing later (not fixed here):** `scanpy_plot_marker_gene_dotplot` rejects the whole call if *any* marker gene is missing. For cross-dataset reuse an `ignore_missing_genes: bool = False` switch would turn a three-step adaptation into one argument. Same for `cluster_labels` — a partial map is already allowed for *categories*, but unknown *keys* are a hard error.

## Commands to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

# 1. Fetch pbmc3k (5.6 MB) and write the raw .h5ad the driver reads
cap10g "$PROJECT_PYTHON" -c "
import scanpy as sc, os
a = sc.datasets.pbmc3k()
a.write_h5ad(os.environ['PROJECT_ROOT']+'/tmp/pbmc3k/00_raw.h5ad', compression='gzip')
print(a.shape)"

# 2. Run all nine tools in tutorial order (~40 s, ~1.2 GB peak)
cap10g "$PROJECT_PYTHON" tmp/pbmc3k/run_pbmc3k.py 2>tmp/pbmc3k/run.err | tee tmp/pbmc3k/run.log
grep cap10g: tmp/pbmc3k/run.err          # peak RAM and OOM count

# 3. Read the per-step records
cat tmp/pbmc3k/results/summary.json
ls tmp/pbmc3k/results/                    # one JSON per call, incl. the 3 expected failures

# Reclaim the 159 MB of chained .h5ad outputs
rm -rf tmp/pbmc3k/outputs
```

## Where things are

- Driver: `Scanpy_Agent/tmp/pbmc3k/run_pbmc3k.py`
- Console log: `tmp/pbmc3k/run.log`, stderr (incl. `cap10g` peak line): `tmp/pbmc3k/run.err`
- Per-call JSON records: `tmp/pbmc3k/results/*.json`, roll-up in `summary.json`, adapted markers in `marker_genes_pbmc3k.json`
- Artifacts (figures, `.h5ad`, ranked-gene CSVs): `tmp/pbmc3k/outputs/<tool>_<timestamp>_<id>/`
- Raw input: `tmp/pbmc3k/00_raw.h5ad`; download cache: `~/Work/paper2agent/data/pbmc3k_raw.h5ad`
- Tools themselves: unchanged at `Scanpy_Agent/src/tools/clustering.py`

## Caveats (honest)

- **No reference comparison.** Unlike note 002, there is no executed-notebook reference for pbmc3k, so tool outputs were **not** compared against independent scanpy calls. What was checked: the tools ran, the numbers are internally consistent and biologically plausible, figures are non-blank, and repeated identical calls are byte-identical. That is a reuse check, not a correctness proof.
- **Step 8d's labels are placeholders**, not cell-type annotations. See failure 3.
- The pytest suite and the MCP acceptance suite from note 002 were **not** rerun — they are pinned to bone-marrow reference files and would not exercise pbmc3k anyway.
- The tools were called in-process rather than over real stdio. Note 002 already verified the stdio path with the same server.
- `sc.datasets.pbmc3k()` writes to `./data` relative to the cwd, so the cache location depends on where you run it from. Harmless here, but it is not governed by `project.env`.

## Next

Phase 1 note 4 of 10. Still open from note 002: Stages 5–6 for `Scanpy_Agent` (clean-environment install, `USAGE.md`, ZIP), and the TISSUE paper as the next conversion. This run needed no new paper and started none.
