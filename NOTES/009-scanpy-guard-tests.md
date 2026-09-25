# 009 — Pytest cases for the three guards (Phase 1, note 9 of 10)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Closes:** the gap left open by [007](007-scanpy-tool-guards.md) and [008](008-scanpy-qc-raw-counts.md) — three guards existed but only standalone drivers exercised them
**Result:** **14 new tests, all passing; suite is now 64/64.** Both guards verified by mutation: disabling each one fails exactly the tests that should fail. **No new fixture needed to enter git.** No tool code changed, no new venv, nothing downloaded.

## Short version

Notes 007 and 008 added three guards and demonstrated them with throwaway scripts under `tmp/`. That left them outside the project's own safety net: a future edit to `clustering.py` would not have been caught by `pytest`. This note moves that coverage into the suite — one new file, `tests/code/clustering/test_tool_guards.py`, 14 tests, **13 s** to run.

The tests do not just assert "guard fires". Five of the fourteen exist to prove the guards **do not** fire on inputs that must be accepted, which is where the actual risk lives.

## What was added

`Scanpy_Agent/tests/code/clustering/test_tool_guards.py` (9.9 KB). It follows the existing suite's conventions — `asyncio_mode = auto`, the `call` / `call_error` / `check_result_contract` / `out_dir` helpers, in-process fastmcp `Client`.

### Fixtures: nothing new in git

Every fixture is derived **in-process** from the `sub4_*` subsamples the suite already creates on demand, then written to pytest's `tmp_path`:

| Builder | Source | Cut to |
|---|---|---|
| `_small_counts()` | `subsample("raw")` | 200 cells × 2,000 genes |
| `_small_lognorm()` | `subsample("normalized")` (carries `uns['log1p']`) | same corner |
| `_small_clustered()` | `subsample("annotated")`, every 20th cell | 214 cells, all 5 res-0.02 groups; genes cut to the marker panel + 300 others |

Each transformed case (scaled, log1p-without-stamp, `normalize_total`, float-dtype, ambient-corrected, clipped) is built from those with two or three lines of scanpy or numpy in the test itself. **So the answer to "any tiny fixture that must live in git" is: none.** Nothing was added to `tests/data/`, and the repository gains one 9.9 KB Python file.

Keeping the slices small is what makes the file fast — 13 s for 14 tests, against ~6 minutes for the 35 tests that work on full reference checkpoints.

### The 14 tests

**QC raw-counts guard — must refuse (5):**

| Test | Rule exercised |
|---|---|
| `test_qc_refuses_scaled_x` | rule 2, negative values (note 006's silent wrong answer) |
| `test_qc_refuses_log1p_stamped` | rule 1, `uns['log1p']` |
| `test_qc_refuses_log1p_without_stamp` | rule 3, fractional + max < 50 (**note 007's hole**) |
| `test_qc_refuses_normalize_total_output` | rule 4, near-constant per-cell totals |
| `test_qc_error_hint_names_counts_layer` | the message names `layers['counts']` when one exists |

`test_qc_refuses_scaled_x` deliberately deletes `uns['log1p']` after scaling, so it tests the negative-value rule rather than being short-circuited by the stamp.

**QC raw-counts guard — must accept (4):**

| Test | Why it matters |
|---|---|
| `test_qc_accepts_raw_counts` | the happy path, plus a non-blank violin figure |
| `test_qc_accepts_float_dtype_counts` | integral counts stored as `float32` — what a dtype check would break |
| `test_qc_accepts_fractional_corrected_counts` | ambient-corrected counts: fractional but real (asserts the fixture clears the rule-3 ceiling, so the test cannot pass for the wrong reason) |
| `test_qc_accepts_integer_counts_below_rule3_ceiling` | counts with max 11 — proves rule 3 needs fractional **and** small |

**Marker `ignore_missing_genes` (5):**

| Test | What it pins |
|---|---|
| `test_dotplot_full_panel_unchanged_by_default` | happy path: 15 sets, `dropped_genes == []` |
| `test_dotplot_missing_gene_error_names_the_switch` | strict default still errors, and names `ignore_missing_genes=True` |
| `test_dotplot_ignore_missing_genes_reports_drops` | one call: plot produced, `dropped_genes` and `emptied_marker_sets` correct, loss stated in the message |
| `test_dotplot_ignore_missing_genes_matches_manual_pruning` | the switch does **exactly** what notes 003/006 did by hand |
| `test_dotplot_all_genes_missing_still_errors` | total loss fails loudly instead of plotting an empty figure |

## Proof the tests can actually fail

A test that cannot fail proves nothing, so I mutated each guard and checked which tests noticed. The source was restored from a backup afterwards and `git diff` confirmed it was byte-identical.

| Mutation | Result |
|---|---|
| Replace `_require_raw_counts(adata)` with `pass` | **5 failed, 9 passed** — exactly the five "must refuse" QC tests |
| Change `if missing and not ignore_missing_genes:` to `if missing:` (switch ignored) | **3 failed, 11 passed** — the three tests that depend on the switch doing work |

The second mutation leaves `test_dotplot_full_panel_unchanged_by_default` and `test_dotplot_missing_gene_error_names_the_switch` passing, which is correct: neither needs the switch to function, only to exist in the error message.

## Suite status

| Run | Result | Time | Peak RAM |
|---|---|---|---|
| New guard file alone | **14 passed** | 13 s | 620 MB |
| Guard file + the 2 tools it covers | **29 passed** | 62 s | 2,271 MB |
| The other 8 test files | **35 passed** | 348 s | 6,449 MB |
| **Total** | **64/64 passed, 0 failed** | ~6.9 min | 6,449 MB |

The suite was 50 tests from note 002 through note 008; it is now **64**.

## RAM, time, disk

- **RAM:** the new file peaks at **620 MB**; the full suite still peaks at 6,449 MB on the unrelated full-reference tests. Under the 10 GB cap, 0 OOM kills.
- **Disk added by this note:** ~12 MB of test outputs under `tests/results/clustering/guard_*`. Fixtures live in `tmp_path` and are removed by pytest.
- **Worth flagging:** `Scanpy_Agent/tests/results/` has grown to **11 GB** across 73 directories. Each tool call writes a fresh timestamped subdirectory, so every suite run accumulates rather than overwrites, and notes 007–009 ran the suite four times. `rm -rf Scanpy_Agent/tests/results` reclaims all of it; the suite recreates what it needs. Disk is at 41 GB used, 195 GB free, so nothing is at risk — but this grows unboundedly if the suite is run often.

## How to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

# Just the guards (~13 s, ~620 MB)
cap10g "$PROJECT_PYTHON" -m pytest tests/code/clustering/test_tool_guards.py -q

# Full suite, in two halves (cap10g does not survive backgrounding - see note 007)
cap10g "$PROJECT_PYTHON" -m pytest \
  tests/code/clustering/test_tool_guards.py \
  tests/code/clustering/test_scanpy_compute_qc_and_filter.py \
  tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py -q       # 29 passed, ~1 min
cap10g "$PROJECT_PYTHON" -m pytest tests/code/clustering -q \
  --ignore=tests/code/clustering/test_tool_guards.py \
  --ignore=tests/code/clustering/test_scanpy_compute_qc_and_filter.py \
  --ignore=tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py # 35 passed, ~6 min

rm -rf tests/results      # reclaims ~11 GB; the suite recreates what it needs
```

## Caveats (honest)

- **The committed test file is not runnable from a fresh clone of the notes repo.** It imports `clustering_verify_helpers`, needs `conftest.py`, and needs `notebooks/clustering/ref/*.h5ad` to build the `sub4_*` subsamples — none of which are in git, because they are gigabytes of `.h5ad`. The file is committed to be *read* and to live alongside the source it guards, not to run standalone. Running it requires the full working tree.
- **Mutation testing covered the two guard bodies, not every branch.** I disabled each guard wholesale; I did not mutate rule 3's threshold or rule 4's tolerance individually, so those constants are exercised but not independently pinned.
- **The deliberate gap from note 008 is still open and still untested**, because it cannot be tested: fractional data with a large maximum and varying totals (TPM/CPM) is indistinguishable from ambient-corrected counts and is accepted by design.
- The new tests use the in-process fastmcp `Client`, like the rest of the suite. The real-stdio MCP acceptance suite was **not** re-run.
- No tool code changed in this note. `clustering.py` is byte-identical to note 008 — verified with `git diff` after the mutation experiments.

## Next

Phase 1 note 10 of 10 — the last one. The only thing still carried forward from note 002 is **Stages 5–6 for `Scanpy_Agent`**: a clean-environment install from `requirements.txt`, `USAGE.md`, and the ZIP with its relocation test. That, plus the 11 GB of `tests/results`, is what remains before this project is tidy.
