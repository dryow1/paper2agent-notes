# use-001 — Scanpy 9 tools on paul15 (mouse blood progenitors)

**Date:** 2026-09-25 · **Machine:** Dell laptop, no GPU, CPU only · **Not a Phase 1 note** — a use of the finished tools.
**Outcome:** **all nine tools ran, 12/12 steps as expected, 30 s, peak 954 MB.** No tool code changed, no new environment, 10 MB downloaded.

## Which file

| | |
|---|---|
| Dataset | `sc.datasets.paul15()` — mouse haematopoietic progenitors, Paul et al. 2015 |
| Download | `paul15.h5`, **10 MB**, to `~/Work/paper2agent/data/` (limit 500 MB) |
| Working copy | `Scanpy_Agent/tmp/paul15/00_shipped.h5ad`, 4 MB |
| Size | **2,730 cells × 3,451 genes** |
| Condition on arrival | **Raw integer counts** (min 0, max 168, all integral), per-cell totals 275–11,712 |
| Annotation it ships with | `paul15_clusters` — **19 published groups** (`1Ery`…`19Lymph`) |
| Driver | `Scanpy_Agent/tmp/paul15/run_paul15.py` (not committed) |

**No counts had to be restored.** Unlike pbmc68k_reduced (note 006), paul15 arrives as genuine raw counts, so the notes 007–008 guards accepted it unchanged. That was checked, not assumed — see the guard step below.

## Which buttons ran

| Step | Tool | Result |
|---|---|---|
| 1 | `scanpy_compute_qc_and_filter` | 1.4 s — kept **2,730/2,730 cells, 3,451/3,451 genes**; `mt_prefix="mt-"` for mouse |
| 2 | `scanpy_detect_doublets_scrublet` | 4.0 s — **1 doublet** in 2,730 cells (0.04 %) |
| 3 | `scanpy_normalize_log1p` | 0.7 s |
| **G** | `scanpy_compute_qc_and_filter` **re-fed the normalised output** | 0.2 s — ⛔ **refused, as designed** (see below) |
| 4 | `scanpy_select_highly_variable_genes` | 0.9 s — 2,000 HVGs of 3,451 |
| 5 | `scanpy_run_pca` | 2.1 s — 50 PCs; PC1 carries 11.1 % of variance, first five 11.1/2.7/1.3/1.1/0.8 % |
| 6 | `scanpy_build_neighbors_umap` | 11.0 s — slowest step |
| 7a | `scanpy_cluster_leiden` res 1.0 | 1.2 s — 10 clusters |
| 7b | `scanpy_cluster_leiden` res 0.02/0.5/2.0 | 1.3 s — **2 / 7 / 15** clusters |
| 8 | `scanpy_plot_marker_gene_dotplot` | 0.4 s — 7 marker sets × 19 published groups, one call |
| 9a | `scanpy_rank_marker_genes` on `paul15_clusters` | 2.2 s — 19 groups, 65,569 rows |
| 9b | `scanpy_rank_marker_genes` on `leiden_res_0.50` | 2.0 s — 7 groups, 24,157 rows |

Settings were tutorial defaults except `mt_prefix="mt-"` (mouse), `batch_key` left `None` (paul15 has no batch column), and `ignore_missing_genes=True` on the dotplot.

**The guard step (G) is the notable one.** I deliberately fed step 3's log-normalised output back into QC, and it was refused:

> `adata.uns['log1p'] is present: the data have already been log-transformed, so QC metrics computed from them are not counts. This object has layers['counts']: set adata.X = adata.layers['counts'] and retry.`

That is the note 007–008 guard working on a dataset it had never seen, including the hint pointing at the counts layer the pipeline itself created.

## What showed up, in plain language

**The cells form a branching continuum, not tidy islands.** The UMAP is one connected mass running from erythroid cells at one end, through shared progenitors in the middle, to granulocyte/monocyte cells at the other — the shape you expect from a differentiation process caught mid-flight. Only two small groups sit clearly apart: the dendritic cells (`11DC`, 30 cells) and the lymphoid cells (`19Lymph`, 31 cells).

**The marker genes the tools found are textbook, and nobody told them the answer.** Ranking genes against the 19 published groups produced:

- **Red-cell groups (`1Ery`–`6Ery`)** → `Car1`, `Car2`, `Hba-a2` (haemoglobin), `Blvrb`, `Ermap`, `Prdx2`. These are haemoglobin and red-cell-membrane genes.
- **Platelet precursors (`8Mk`)** → `Itga2b` (CD41), `Pf4` (platelet factor 4), `Nrgn`, `Vwf`. Exclusive to that one group in the dotplot.
- **Monocyte and neutrophil groups (`14Mo`, `15Mo`, `16Neu`)** → `Mpo`, `Elane`, `Prtn3`, `Ctsg` — the granule enzymes those cells package.
- **Eosinophils (`18Eos`, 9 cells)** → `Prg2`, `Prg3`, `Epx` (eosinophil peroxidase), `Cebpe`. Nine cells, and still a clean signature.
- **Dendritic cells (`11DC`)** → `Cd74`, `H2-Ab1`, `Cst3` — MHC class II, the antigen-presenting machinery.
- **Lymphoid (`19Lymph`)** → `Ccl5`, `Cd52`, `Ctsw`.
- **Stem/progenitor groups (`9GMP`, `10GMP`)** → `Cd34` and `Flt3` light up in the dotplot, as progenitor markers should.

**One disagreement worth flagging.** Group `17Neu` is labelled a neutrophil group, but the two genes that stand out for it in my panel are `Prss34` and `Mcpt8` — **basophil/mast-cell** markers — and its own top-ranked genes (`Srgn`, `Lmo4`, `Prss34`) do not look neutrophil-like. It has 22 cells. I am not claiming the published label is wrong; I am recording that the markers I chose put it with the basophils, and that a 22-cell group is where I would look first if this mattered.

**Our own clustering is coarser than the published labels, and sensibly so.** Leiden at resolution 0.5 gave 7 clusters against the paper's 19:

| Leiden cluster | Cells | Mostly which published group |
|---|---|---|
| 0 | 520 | `7MEP` (31 %) — the mixed progenitor middle |
| 1 | 512 | `14Mo` (57 %) |
| 2 | 841 | `2Ery` (39 %) — swallows the whole `1Ery`–`6Ery` series |
| 3 | 364 | `16Neu` (43 %) |
| 4 | 434 | `10GMP` (32 %) |
| 5 | 30 | `19Lymph` (**100 %**) |
| 6 | 29 | `11DC` (**100 %**) |

The two tiny clusters are perfectly pure; the big ones merge the paper's fine-grained stages, which is what a single resolution on a continuum does. At resolution 0.02 it collapses to just 2 clusters — essentially "red-cell side" and "white-cell side".

## RAM and time

- **30 s wall clock** for all 12 steps. Slowest: neighbours/UMAP at 11.0 s.
- **Peak 954 MB** under `cap10g` (1,024 MB in-process), 0 OOM kills — about 10 % of the 10 GB cap.
- Disk: 114 MB of run outputs under `Scanpy_Agent/tmp/paul15/`, plus the 10 MB download. `rm -rf Scanpy_Agent/tmp/paul15/outputs` reclaims most of it.

## What still failed, or degraded

1. **No mitochondrial, ribosomal or haemoglobin QC genes exist in this dataset.** `n_mt_genes: 0`, `n_ribo_genes: 0`, `n_hb_genes: 0`, so `pct_counts_mt` and `pct_counts_ribo` are **identically zero for all 2,730 cells**. paul15 ships a pre-filtered 3,451-gene set with no `mt-` genes at all. The QC violin and the QC scatter's colour axis are therefore flat and carry no information. **Nothing errored — the tool did what it was asked — but that third violin panel is meaningless here**, and a user skimming figures could mistake "no mitochondrial contamination" for a finding. This is the same *class* of problem as note 006's silent wrong answer, just milder: the guards check that X is counts, not that the requested QC genes exist.
2. **Two marker genes are absent:** `Vpreb1` and `Igll1`, so the Lymphoid set plotted with `Dntt` alone. Handled in one call by `ignore_missing_genes=True` and reported as `dropped 2 missing genes` — exactly the note 007 switch doing its job, with no hand-editing.
3. **`min_genes=100` filtered nothing** (the sparsest cell still detects 213 genes), and `min_cells=3` removed no genes. The tutorial's thresholds are simply not binding on an already-curated matrix.
4. **Scrublet found 1 doublet in 2,730 cells.** Plausible for pre-filtered data, but low enough that I would not read anything into it.

Nothing crashed; no tool needed an edit; the guards fired exactly once, on purpose.
