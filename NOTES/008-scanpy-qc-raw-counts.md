# 008 — Closing the QC raw-counts hole (Phase 1, note 8 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Closes:** the open hole left by [007](007-scanpy-tool-guards.md) — the QC guard caught scaled data but waved log-normalised data through
**Result:** **Closed. 10/10 guard cases behaved as expected, 50/50 test suite still passes, both happy paths still run.** Same `scanpy-env`, no new venv, nothing downloaded.

## Short version

Note 007's guard tested one thing: negative values. That caught z-scored data and nothing else, so `pbmc68k_reduced`'s own `.raw.X` — log-normalised, range 0–6.49, perfectly non-negative — still sailed through and still produced QC numbers that are not counts.

The guard now applies **four rules, each keyed on a positive signature of a specific transformation** rather than on "not integral". That is the whole design decision: *not integral* would have been one line and would have rejected legitimately fractional count matrices. `_require_nonnegative_counts` became `_require_raw_counts`, +51 lines (564 → **615**).

## The four rules

| # | Signature | Catches | Why false positives are unlikely |
|---|---|---|---|
| 1 | `uns['log1p']` present | anything scanpy log-transformed | scanpy sets this stamp itself; nothing else does |
| 2 | any negative value | scaled / z-scored | counts are never negative |
| 3 | fractional values **and** max < 50 | log1p'd data | `log1p(100000)` is only 11.5, while a real count matrix with a max under 50 is integral — it takes **both** conditions |
| 4 | fractional values **and** per-cell totals near-constant (`std/mean < 1e-3`) | `normalize_total(target_sum=…)` | real cells never all have identical totals |

Integral, non-negative data returns immediately at rule 2's exit — the common path costs one `min()` and one modulo scan.

Rule 3 is the one that closes the note 007 hole, and rule 1 alone is not enough: a log1p'd object whose `uns` stamp was dropped (round-tripped through some other tool, or built by hand) has no stamp to find. Both cases are tested below.

The error messages now end with a hint chosen from the object itself:

```
X has fractional values and a maximum of only 6.489: the data look log-transformed,
not raw counts, so 'total_counts' would not be a count. Supply the raw count matrix instead.
```

and when the object actually carries counts:

```
X contains negative values (min -2.032): the data look scaled/z-scored, not raw counts,
and QC metrics computed from them are meaningless. This object has layers['counts']:
set adata.X = adata.layers['counts'] and retry.
```

`_counts_hint()` names `layers['counts']` only when it exists, and mentions `.raw` with a caveat that it is often log-normalised — which, on `pbmc68k_reduced`, it is.

## What ran

Driver: `Scanpy_Agent/tmp/guards/run_rawcounts.py`. **9 s, peak 521 MB.** Five cases that must be refused, five that must be accepted — the second half is the false-positive watchlist, which is where the risk actually lives.

### Must be refused

| # | Input | Result |
|---|---|---|
| R1 | pbmc3k `normalize_total` + `log1p` | ✅ caught by **rule 1** (`uns['log1p']`) |
| R2 | same, with the `uns['log1p']` stamp deleted | ✅ caught by **rule 3** — `maximum of only 5.955` |
| R3 | pbmc3k `normalize_total(target_sum=1e4)`, no log | ✅ caught by **rule 4** — `near-identical per-cell totals (~1e+04)` |
| R4 | pbmc68k_reduced as shipped (z-scaled) | ✅ caught by **rule 2**, note 007's case still fails, now with the counts-layer hint |
| R5 | **pbmc68k_reduced's own `.raw.X`** (log-normalised, max 6.49) | ✅ caught by **rule 3** — **this is the note 007 hole, now closed** |

### Must be accepted (false-positive watchlist)

| # | Input | Result |
|---|---|---|
| A1 | pbmc3k raw integer counts | ✅ kept 2700/2700 cells, 13714/32738 genes — note 003 happy path unchanged |
| A2 | pbmc68k_reduced counts-restored | ✅ kept 700/700 cells, 765/765 genes — note 006 happy path unchanged |
| A3 | integer counts stored as `float32` | ✅ accepted — the obvious false positive an integrality check on dtype would cause |
| A4 | ambient-corrected counts (SoupX/decontX style: fractional, large max, varying totals) | ✅ accepted — fractional but genuinely counts |
| A5 | real counts clipped to max 11 (integer, below rule 3's threshold) | ✅ accepted — proves rule 3 needs *both* conditions, not just a low max |

**10/10 as expected.** A3, A4 and A5 exist specifically because each is a case that a naive "reject non-integers" or "reject small maxima" rule would have broken.

## What failed

**One case failed on the first run, and it was my fixture, not the guard.** `A5` was originally a synthetic 200 × 300 random-integer matrix with invented gene names. The guard accepted it correctly, but the tool then died downstream in plotting:

```
ToolError: ... Positions outside range of features.
```

I checked the guard in isolation before touching anything — `_require_raw_counts` returned cleanly, confirming it was not a false positive — then replaced the fixture with real pbmc3k counts clipped to a maximum of 11, which exercises the same rule-3 boundary on data that can actually be plotted. The lesson is about fixture realism, not about the guard.

## Regression check

| Run | Result | Time | Peak RAM |
|---|---|---|---|
| Raw-counts matrix (this note) | **10/10 as expected** | 9 s | 521 MB |
| Note 007's guard demo, re-run | **8/8 as expected** | 6 s | 556 MB |
| 2 directly affected test files | **15 passed** | 50 s | 1,832 MB |
| The other 8 test files | **35 passed** | 364 s | 5,771 MB |
| **Suite total** | **50/50 passed, 0 failed** | ~6.9 min | 5,771 MB |

Nothing regressed. The bone-marrow reference data are raw integer counts, so rules 2–4 never fire on them, and the QC tests feed either `00_raw_concat.h5ad` or a 10x `.h5` — both counts.

## RAM, time, disk

- **Guard work:** 9 s at 521 MB; the whole verification round including the suite, ~7.5 min, peak 5,771 MB — under the 10 GB cap, 0 OOM kills.
- **Disk:** **61 MB of fixtures** under `Scanpy_Agent/tmp/guards/fixtures/` (seven `.h5ad` files; `corrected_counts.h5ad` is 23 MB because shaving the ambient fraction turns a sparse integer matrix into dense-ish float64). `Scanpy_Agent/tmp/guards/` totals 136 MB. **Nothing downloaded.** Source grew 51 lines.
- `rm -rf Scanpy_Agent/tmp/guards/fixtures` reclaims the 61 MB; the driver rebuilds them on demand.

## How to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

# The 10-case raw-counts matrix (~9 s, ~520 MB). Needs the note 003 / 006 outputs present.
cap10g "$PROJECT_PYTHON" tmp/guards/run_rawcounts.py 2>tmp/guards/raw.err | tee tmp/guards/raw.log
cat tmp/guards/results_raw/summary.json

# Note 007's cases still pass too
cap10g "$PROJECT_PYTHON" tmp/guards/run_guards.py 2>/dev/null | tail -2

# Regression suite, in two halves (cap10g does not survive backgrounding - see note 007)
cap10g "$PROJECT_PYTHON" -m pytest \
  tests/code/clustering/test_scanpy_compute_qc_and_filter.py \
  tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py -q       # 15 passed, ~50 s
cap10g "$PROJECT_PYTHON" -m pytest tests/code/clustering -q \
  --ignore=tests/code/clustering/test_scanpy_compute_qc_and_filter.py \
  --ignore=tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py # 35 passed, ~6 min
```

## Caveats (honest)

- **The guard is heuristic, and one gap is deliberate.** Fractional data with a large maximum *and* varying per-cell totals — e.g. TPM/RPKM/CPM values that were not scaled to a constant sum — is still accepted. It looks exactly like ambient-corrected counts (A4), which must be accepted, so no rule can separate them from the matrix alone. Rejecting both would trade a silent-wrong-answer for a false positive on a legitimate input. **Documented, not fixed.**
- **Rule 3's threshold of 50 is a judgement call**, not a derived constant. It sits well above any plausible `log1p` output (~11.5 for 100,000 counts) and below any fractional matrix I would expect to be counts. A count matrix with a max under 50 is fine as long as it is integral (A5 proves it).
- **Still no pytest cases for the guards.** Same gap as note 007: both guards are covered by the standalone driver, not by the suite, so a future edit to `clustering.py` would not be caught by `pytest` alone. The 10 cases here are written to be portable into `test_scanpy_compute_qc_and_filter.py` — that is the obvious next step and I did not take it.
- **The MCP acceptance suite was not re-run**, so the real-stdio path has not been re-verified since note 007. Its cases use raw counts and do not exercise the new code, but that is an assumption.
- Only `scanpy_compute_qc_and_filter` changed. The other eight tools are untouched since note 007.

## Next

Phase 1 note 9 of 10. Carried forward: Stages 5–6 for `Scanpy_Agent` from note 002 (clean-environment install, `USAGE.md`, ZIP), and **pytest cases for the three guards now in place** — the note 007 pair plus this one — which is the last thing keeping them outside the project's own safety net.
