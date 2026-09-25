# 007 — Two guards for the Scanpy tools (Phase 1, note 7 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Fixes:** the silent-wrong-answer from [006](006-scanpy-pbmc68k-reduced.md) and the marker-panel rough edge first flagged in [003](003-scanpy-pbmc3k.md)
**Result:** **Both fixed, 44 lines across 2 of the 9 tools. 8/8 demonstration cases behaved as expected, and the existing 50-test suite still passes 50/50.** Same `scanpy-env`, no new venv, nothing downloaded.

## Short version

Two problems, two guards, no rewrite. `scanpy_compute_qc_and_filter` now refuses scaled data instead of reporting success on garbage, and `scanpy_plot_marker_gene_dotplot` gained an `ignore_missing_genes` switch so a marker panel with absent genes takes one call instead of hand-editing the dict. The other seven tools were not touched.

## Fix 1 — QC refuses scaled data

**The problem (note 006):** fed `pbmc68k_reduced`'s scaled `X`, the tool reported `"QC computed; kept 700/700 cells and 765/765 genes"` while writing negative `total_counts` for 389 of 700 cells and `pct_counts_mt` between −1167 % and +444 %.

**The fix:** a precondition helper, called immediately after the read and *before* the output directory is created, so a rejected call leaves nothing behind.

```python
def _require_nonnegative_counts(adata: ad.AnnData) -> None:
    """Reject scaled/z-scored matrices, which produce meaningless QC metrics.

    sc.pp.calculate_qc_metrics sums X per cell, so a scaled matrix yields negative
    'total_counts' and out-of-range 'pct_counts_*' while reporting success.
    """
    X = adata.X
    xmin = float(X.min()) if X.size else 0.0
    if xmin < 0:
        raise ValueError(
            f"X contains negative values (min {xmin:.4g}): the data look scaled/z-scored, "
            "not raw counts, and QC metrics computed from them are meaningless. "
            "Supply raw counts - e.g. adata.X = adata.layers['counts'] if a counts layer exists."
        )
```

The message names the fix, because on this dataset the fix really is one line — the counts layer is right there.

## Fix 2 — `ignore_missing_genes`

**The problem (notes 003 and 006):** any absent marker gene rejected the whole call, so both reuse runs had to compute the intersection by hand, drop emptied sets and call a second time. On pbmc68k_reduced that meant editing out 34 of ~60 genes.

**The fix:** a new keyword, default `False` so existing behaviour is unchanged.

```python
ignore_missing_genes: Annotated[bool, "Drop marker genes absent from var_names instead of
    failing; dropped genes and emptied sets are reported"] = False
```

Three deliberate properties:

- **The strict default now names the way out.** The error gained `"; pass ignore_missing_genes=True to drop these 34 genes and plot the rest"`.
- **The degradation is reported, not hidden.** `dropped_genes` and `emptied_marker_sets` are always in the returned dict, and the message says so inline: `"Dotplot of 13 marker sets across 3 groups of 'leiden_res_0.02' (dropped 34 missing genes, emptied 2 sets: ['Erythroblast', 'Plasmablast'])"`. That was the point of note 006's complaint — a caller should be *told*, not left to discover it by set arithmetic.
- **Total loss still fails.** If no gene survives, it raises rather than plotting an empty figure.

**Size of both changes:** `src/tools/clustering.py` went from 520 to **564 lines (+44)**. Two of nine tools touched; seven untouched.

## What ran

Driver: `Scanpy_Agent/tmp/guards/run_guards.py`, reusing the clustered `.h5ad` outputs from notes 003 and 006. **6 s, peak 423 MB.** Each case declares the outcome it expects, so a silent regression would show as `UNEXPECTED`.

| # | Case | Expected | Got |
|---|---|---|---|
| G1a | QC on pbmc68k_reduced's scaled `X` | error | ✅ `X contains negative values (min -2.032): the data look scaled/z-scored...` |
| G1b | QC on pbmc68k counts-restored | ok | ✅ kept 700/700 cells, 765/765 genes |
| G1c | QC on pbmc3k raw | ok | ✅ kept 2700/2700 cells, 13714/32738 genes |
| G2a | pbmc68k dotplot, tutorial panel, default | error | ✅ lists 34 genes **+ names the switch** |
| G2b | pbmc68k dotplot, `ignore_missing_genes=True` | ok | ✅ 13 sets × 3 groups, dropped 34, emptied `['Erythroblast', 'Plasmablast']` |
| G2a | pbmc3k dotplot, tutorial panel, default | error | ✅ lists 8 genes + names the switch |
| G2b | pbmc3k dotplot, `ignore_missing_genes=True` | ok | ✅ **15 sets** × 3 groups, dropped 8, no set emptied |
| G2c | panel where no gene exists | error | ✅ `no marker gene is present in var_names; all 1 sets are empty` |

**8/8 as expected.** The figures render correctly — the pbmc68k dotplot shows the 13 surviving sets with `Erythroblast` and `Plasmablast` absent, exactly as reported.

The two things asked for are both demonstrated: **G1a is the old silent-wrong case now erroring**, and **G2b is the missing-gene panel working in one call with no second manual pass** — the unedited tutorial marker dict goes straight in on both datasets.

## Regression check

Changing shipped tool code means the note-002 verification suite has to still pass. It does:

| Run | Result | Time | Peak RAM |
|---|---|---|---|
| The 2 directly affected test files | **15 passed** | 62 s | 2,167 MB |
| The other 8 test files | **35 passed** | 453 s | 6,129 MB |
| **Total** | **50/50 passed, 0 failed** | ~8.6 min | 6,129 MB |

That matches note 002's original 50/50 and its ~6,086 MB peak. Neither guard changed any result on the bone-marrow reference data — as expected, since raw counts are non-negative and the tutorial panel is fully present there.

**A process note:** the first attempt to run the suite in the background produced a 1-byte log. `cap10g` uses `systemd-run --user --scope`, and the transient scope was torn down when the launching shell exited, killing pytest after one test. Re-run in the foreground in two halves to fit the command timeout. Worth remembering: **`cap10g` and background execution do not mix.**

## RAM, time, disk

- **Demonstration:** 6 s, peak **423 MB** under `cap10g`, 0 OOM kills.
- **Regression suite:** ~8.6 min, peak **6,129 MB** — the heaviest thing in this project, but still under the 10 GB cap.
- **Disk:** ~17 MB added under `Scanpy_Agent/tmp/guards/` (two copied clustered `.h5ad` inputs plus small outputs). **Nothing downloaded.** Source grew by 44 lines.

## How to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

# Both guards, 8 cases (~6 s, ~420 MB). Needs the note 003 / 006 outputs present.
cap10g "$PROJECT_PYTHON" tmp/guards/run_guards.py 2>tmp/guards/run.err | tee tmp/guards/run.log
cat tmp/guards/results/summary.json

# Regression suite, in two halves (cap10g will not survive backgrounding)
cap10g "$PROJECT_PYTHON" -m pytest \
  tests/code/clustering/test_scanpy_compute_qc_and_filter.py \
  tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py -q      # 15 passed, ~1 min
cap10g "$PROJECT_PYTHON" -m pytest tests/code/clustering -q \
  --ignore=tests/code/clustering/test_scanpy_compute_qc_and_filter.py \
  --ignore=tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py  # 35 passed, ~7.5 min
```

## Caveats (honest)

- **The QC guard catches scaled data, not log-normalised data.** It tests for negative values. A log-normalised matrix (like this dataset's `raw.X`, range 0–6.49) is non-negative, so it would still pass and still produce wrong QC numbers, just less obviously wrong. Closing that needs an integrality check, which risks false positives on float-typed count matrices — I did not attempt it here. **This remains an open hole.**
- **No new tests were written.** The two guards are demonstrated by the 8-case driver, not by additions to the pytest suite. The suite was run to prove nothing broke, not to cover the new paths. Adding 3–4 cases to `test_scanpy_compute_qc_and_filter.py` and `test_scanpy_plot_marker_gene_dotplot.py` is the obvious follow-up.
- `ignore_missing_genes` changes which genes are plotted, so **a dotplot made with it is not comparable to one made from the full panel** — the returned `dropped_genes` is what makes that auditable.
- The tools were called in-process via fastmcp `Client`, not over real stdio. The note-002 MCP acceptance suite was **not** re-run, so the stdio path has not been re-verified since the change. Its cases do not exercise either new code path, but that is an assumption, not a measurement.
- Only the two tools named here were changed; the other seven are byte-identical to note 002.

## Next

Phase 1 note 8 of 10. Carried forward: Stages 5–6 for `Scanpy_Agent` from note 002 (clean-environment install, `USAGE.md`, ZIP), pytest cases for the two new guards, and the log-normalised-input hole above.
