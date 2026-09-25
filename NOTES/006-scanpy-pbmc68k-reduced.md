# 006 — The nine Scanpy tools on pbmc68k_reduced (Phase 1, note 6 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Builds on:** [002](002-scanpy-clustering.md) (the nine tools) and [003](003-scanpy-pbmc3k.md) (first reuse). **Same `scanpy-env`, no new venv, no tool code changed.**
**Result:** **All nine tools ran.** 13 of 15 calls succeeded; the 2 failures were deliberate unadapted calls. 32 s, peak 708 MB. **One silent-wrong-answer found** — see below.

## Short version

`pbmc68k_reduced` is a **processed** dataset, not a raw one: `X` is z-scaled (range −2.03 to 28.41) and it already ships `layers['counts']`, `raw`, `X_pca`, `X_umap`, a neighbour graph and louvain labels. That makes it a harder reuse test than pbmc3k, which was raw counts with empty `obs`.

The tools handled it after **one prep step** — put the integer counts back in `X` and drop the precomputed results. The important finding is what happened *without* that step: **the QC tool accepted scaled data and produced meaningless numbers without complaining.** No switch was broken, so no tool code was changed.

## Inputs and limits

| | |
|---|---|
| Data | `sc.datasets.pbmc68k_reduced()` — **bundled inside the installed scanpy package, nothing downloaded at all** |
| Size | 700 cells × 765 genes. Saved as `00_raw.h5ad`, **2.2 MB** |
| As shipped | `X` scaled (−2.03 … 28.41); `layers['counts']` int32 raw counts (0 … 259); `raw.X` log-normalised (0 … 6.49); `obsm` = `X_pca`, `X_umap`; `obsp` = connectivities/distances; `obs` = `bulk_labels`, `n_genes`, `percent_mito`, `n_counts`, `S_score`, `G2M_score`, `phase`, `louvain` |
| Limits | CPU only, `cap10g` on every command; peak **708 MB** (cap 10 GB, 0 OOM kills); **0 bytes downloaded** |
| Environment | Unchanged `Scanpy_Agent/scanpy-env`, scanpy `1.13.0a3.dev0+g0d5fd1623` |
| Tool code | **Unchanged.** No switch was broken. |

## What ran

Driver: `Scanpy_Agent/tmp/pbmc68k/run_pbmc68k.py`, calling the tools through an in-process fastmcp `Client` (same path as note 003).

### Pass A — as shipped, scaled `X`

| # | Tool | Result |
|---|---|---|
| A1 | `scanpy_compute_qc_and_filter` | ⚠️ **"ok" — but the numbers are nonsense.** See below. |
| A2 | `scanpy_normalize_log1p` | ⛔ `adata.layers['counts'] already exists; refusing to overwrite it` — **correct refusal**, the guard did its job |

### Pass B — after restoring counts (`X ← layers['counts']`, drop layers/raw/obsm/obsp/uns)

| # | Tool | Result |
|---|---|---|
| B1 | `scanpy_compute_qc_and_filter` | ✅ 10.9 s — 700/700 cells, 765/765 genes kept; 1 mt, 4 ribo, 0 hb genes |
| B2 | `scanpy_detect_doublets_scrublet` | ✅ 3.9 s — **3 doublets** / 700 cells (0.43 %) |
| B3 | `scanpy_normalize_log1p` | ✅ 0.1 s |
| B4a | `scanpy_select_highly_variable_genes` (`n_top_genes=2000`) | ⚠️ ✅ — returned **765 of 765**, i.e. every gene. Silently under-delivers; see below. |
| B4b | same, `n_top_genes=400` | ✅ 0.2 s — 400 HVGs |
| B5 | `scanpy_run_pca` | ✅ 0.6 s — 50 PCs on 400 HVGs; variance ratio first 5 = 0.244 / 0.100 / 0.047 / 0.036 / 0.030, total 0.675 |
| B6 | `scanpy_build_neighbors_umap` | ✅ 8.2 s — slowest step, as on pbmc3k |
| B7a | `scanpy_cluster_leiden` res 1.0 + QC panel | ✅ 0.5 s — **10 clusters** |
| B7b | `scanpy_cluster_leiden` res 0.02/0.5/2.0 | ✅ 0.3 s — **3 / 9 / 15 clusters** |
| B8a | `scanpy_plot_marker_gene_dotplot`, tutorial markers | ⛔ **expected failure**, 34 genes missing |
| B8b | same, adapted markers | ✅ 0.2 s — 13 sets × 3 groups |
| B9a | `scanpy_rank_marker_genes` on `leiden_res_0.50` | ✅ 1.1 s — 9 groups, 6,885 rows |
| B9b | `scanpy_rank_marker_genes` on **`bulk_labels`** | ✅ 0.8 s — 10 groups, 7,650 rows |

All 16 figures render with content (none under 2 KB); the three-panel UMAP shows clean separated clusters.

**B9b is the new capability.** pbmc3k had no cell-type column, so `groupby` could only ever take a Leiden key. Here the shipped `bulk_labels` annotation let the same tool rank markers against real immune cell types, and the results are sensible: `CD3D/CD52/LDHB` for CD4+ T Reg, `CCL5/NKG7/GZMA` for CD8+ Cytotoxic T, `MZB1/IGJ/TNFRSF17` for the plasma-cell cluster. That is a `groupby` switch doing exactly what it was designed for on a dataset the tools had never seen.

## The one real problem: QC runs happily on scaled data

`scanpy_compute_qc_and_filter` accepted the scaled matrix and reported success — `"QC computed; kept 700/700 cells and 765/765 genes"`. The metrics it wrote are meaningless:

| metric | on scaled `X` | should be |
|---|---|---|
| `total_counts` | −109.05 … 212.73, mean **−0.35** | a positive count |
| cells with **negative** `total_counts` | **389 of 700** | 0 |
| `pct_counts_mt` | **−1166.7 %** … **+444.2 %** | 0–100 % |
| `n_genes_by_counts` | 763 … 765 | varies per cell |

The `min_genes=100` filter then removed nothing, because in a scaled matrix essentially every entry is non-zero, so every cell "detects" all 765 genes.

**This is a missing precondition, not a broken switch.** The tool documents `data_path` as "Raw-count AnnData" but never checks it. `scanpy_normalize_log1p` has exactly the right guard (`uns['log1p']` / `layers['counts']`) and it fired correctly in A2 — the QC tool has no equivalent. Per the "don't change tool code unless a switch is clearly broken" rule, **I changed nothing**; a switch working wrongly is not the same as a switch being broken. Recommended fix for later, one guard at the top of the tool:

```python
if adata.X.min() < 0:
    raise ValueError("X has negative values: data appear scaled/z-scored, not raw counts")
```

This is more dangerous than any failure in note 003, because every error there was loud. This one is silent, and a user chaining tools would carry the bad QC columns all the way to a UMAP panel coloured by `pct_counts_mt`.

### Second, milder case: `n_top_genes=2000` on 765 genes

Asking for 2,000 HVGs from a 765-gene dataset returned all 765, reported as `"Selected 765 highly variable genes of 765"`. That is upstream scanpy's behaviour, and the message is honest about the count, but the *switch* silently did nothing — HVG selection that selects everything is not selection. A caller comparing across datasets would not notice. Worth a warning when `n_top_genes >= n_vars`.

## Which switches needed `ignore_missing_genes` or similar

The rough edge flagged at the end of note 003 bit again, harder:

| | pbmc3k | pbmc68k_reduced |
|---|---|---|
| Genes in dataset | 13,714 | **765** |
| Tutorial marker genes missing | 8 | **34** |
| Marker sets emptied | 0 | **2** (`Erythroblast`, `Plasmablast`) |
| Sets left with a single gene | 0 | **7** |

`scanpy_plot_marker_gene_dotplot` still rejects the whole call if *any* marker gene is absent, so the adaptation was again: compute the intersection myself, drop empty sets, re-call. **An `ignore_missing_genes: bool = False` switch would have turned this into one argument on both datasets.** With 765 genes it matters more — a third of the marker panel is gone, and two cell types can no longer be scored at all, which the caller should be *told* rather than having to discover by set arithmetic. A switch that reported `dropped_genes` and `emptied_sets` in its return dict would make the degradation visible instead of silent.

No other switch needed adapting. `batch_key` was left at its `None` default (pbmc68k_reduced has no batch column either; `bulk_labels` is a cell-type annotation, not a batch). `mt_prefix`, `ribo_prefixes`, `hb_pattern`, `min_genes`, `min_cells`, `n_comps`, `n_neighbors`, `resolutions`, `n_iterations`, `focus_group` all took tutorial defaults unchanged.

## Compared with pbmc3k

| | pbmc3k (note 003) | pbmc68k_reduced |
|---|---|---|
| Ships as | raw counts, **empty `obs`** | **scaled X** + counts layer + raw + PCA/UMAP/graph/louvain |
| Download | 5.6 MB | **0 — bundled in the scanpy package** |
| Cells × genes | 2,700 × 13,714 | 700 × 765 |
| Prep needed before tool 1 | none | **one step: `X ← layers['counts']`, drop precomputed results** |
| Wall time / peak RAM | 37 s / 1,153 MB | **32 s / 708 MB** |
| Doublets | 32 (1.19 %) | 3 (0.43 %) |
| Leiden res 1.0 | 10 clusters | 10 clusters |
| Leiden 0.02 / 0.5 / 2.0 | 3 / 8 / 15 | 3 / 9 / 15 |
| Rank-genes rows | 109,712 | 6,885 |
| Marker genes dropped | 8, 0 sets emptied | 34, 2 sets emptied |
| Failures | 3, all loud | 2 loud + **1 silent wrong answer** |
| New capability | — | `groupby="bulk_labels"` — ranking against shipped cell types |

Cluster counts landing so close (10 at res 1.0 on both; 3/9/15 vs 3/8/15) across datasets 4× apart in cell count and 18× in gene count is a reassuring sign that the Leiden switches behave consistently.

The headline difference: **pbmc3k tested whether the tools handle missing inputs; pbmc68k_reduced tested whether they handle inputs of the wrong kind.** They pass the first test cleanly and fail the second quietly.

## RAM, time, disk

- **RAM:** peak **708 MB** under `cap10g` (862 MB in-process maxrss), 0 OOM kills — well under the 10 GB cap and lighter than pbmc3k's 1,153 MB.
- **Time:** **32 s** for all 15 calls. Slowest: `scanpy_compute_qc_and_filter` at 10.9 s (first call, includes scanpy import and numba warm-up) and `scanpy_build_neighbors_umap` at 8.2 s.
- **Disk:** **23 MB** added, all under `Scanpy_Agent/tmp/pbmc68k/` — chained `.h5ad` outputs, 16 figures, JSON records. **Nothing downloaded.** `rm -rf Scanpy_Agent/tmp/pbmc68k/outputs` reclaims most of it.

## Commands to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

# 1. Write the shipped dataset out (no download; it is inside the scanpy package)
cap10g "$PROJECT_PYTHON" -c "
import scanpy as sc, os
sc.datasets.pbmc68k_reduced().write_h5ad(
    os.environ['PROJECT_ROOT']+'/tmp/pbmc68k/00_raw.h5ad', compression='gzip')"

# 2. Both passes, all nine tools (~32 s, ~700 MB peak)
cap10g "$PROJECT_PYTHON" tmp/pbmc68k/run_pbmc68k.py 2>tmp/pbmc68k/run.err | tee tmp/pbmc68k/run.log
grep cap10g: tmp/pbmc68k/run.err

cat tmp/pbmc68k/results/summary.json
rm -rf tmp/pbmc68k/outputs          # reclaim ~23 MB
```

## Caveats (honest)

- **No reference comparison**, same as note 003. Outputs were not checked against independent scanpy calls; what was checked is that the tools run, numbers are internally consistent and biologically plausible, and figures are non-blank.
- **A driver bug of mine, fixed:** the first prep attempt did `del a.layers[key]` over all keys. `adata.layers` on this dataset has a `None` key that **aliases `.X` itself**, so deleting it wiped the matrix and the run died with `'NoneType' object has no attribute 'max'`. Restricted to string keys. This is an anndata quirk, not a tool problem.
- **Pass B's prep discards the dataset's own PCA/UMAP/louvain** so the tools recompute from counts. The Leiden clusters here are therefore *not* the shipped `louvain` labels, and the two were not compared.
- Neither the note-002 pytest suite nor the MCP acceptance suite was rerun; both are pinned to bone-marrow reference files.
- The tools were called in-process, not over real stdio (note 002 verified that path).

## Next

Phase 1 note 7 of 10. Carried forward: Stages 5–6 for `Scanpy_Agent` from note 002, and now **two concrete tool fixes worth making together** — a raw-count guard on `scanpy_compute_qc_and_filter`, and `ignore_missing_genes` on `scanpy_plot_marker_gene_dotplot` (asked for in note 003, more clearly needed after this run).
