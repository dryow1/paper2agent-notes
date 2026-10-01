# 012 — A numeric baseline for the TISSUE wrappers

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Builds on:** [004](004-tissue-tutorial.md), [005](005-tissue-tools.md) — same checkout, **same
`tissue-env-pinned`, no new venv, nothing installed, nothing downloaded**
**Result:** **32 assertions, 32 passed.** Every number note 005 printed reproduced exactly, to
the four decimal places it printed them. 15.7 s, peak **407 MB** — the same peak note 005 recorded.

## Short version

Note 005 ended with an admission: *"No test suite. Correctness rests on the numbers matching
note 004's independent run — not on independent verification."* Those numbers lived in prose, in
a file nothing reads. Upstream was no help — `repo/TISSUE/test.py` contains **zero** assert
statements, which I verified rather than assumed (`grep -c assert` → 0).

This note converts note 005's printed numbers into
`TISSUE_Agent/tests/test_tissue_numeric_baseline.py`: 32 named tests that **re-run the real
pipeline** and compare each reported value against the recorded one. It is not a smoke test.
Nothing asserts "it ran"; every assertion names a quantity and the value it must equal.

The headline is that the numbers held. Six days after note 005, in a fresh process, SpaGE
prediction and conformal calibration returned **bit-for-bit the same figures at 4 dp**.

## Why the numbers were reproducible

Not luck. Upstream fixes every seed on the path:

| Where | Seed |
|---|---|
| `predict_gene_expression(..., random_seed=444)` | default argument, `main.py:268` |
| `np.random.seed(random_seed)` before shuffling calibration genes | `main.py:342` |
| every `PCA(...)` and `KMeans(...)` in the calibration path | `random_state=444`, `main.py:607, 741, 747, 761, 769, 770, 810` |

The wrappers override none of them, and pass no seed of their own. So exact reproducibility is
the *expected* behaviour, and this file is what makes a regression in it visible.

## What the test checks

| Group | Tests | What it pins |
|---|---|---|
| Inputs | 2 | the three `tests/data/*.txt` are present; upstream `test.py` still has no asserts to reconcile against |
| `tissue_predict_spatial_gene` | 8 | 3,405 cells · 31 calibration genes · **Pearson r = 0.3097** · predicted mean 0.5188 · measured mean 2.5771 · recorded settings · a readable `.h5ad` |
| `tissue_calibrate_prediction_intervals`, α = 0.23 | 6 | **coverage 0.7883** (nominal 0.77) · width 12.5879 · held-out gene 0.8537 · width 7.5545 |
| same, α = 0.5 | 5 | **coverage 0.5397** (nominal 0.50) · width 4.1276 · held-out gene 0.6135 · width 2.6243 |
| Properties, not just equalities | 4 | coverage ≥ nominal at both α · halving nominal coverage must narrow intervals ≥ 3× · target-gene intervals narrower than the calibration set · step 2 inherits `target_gene`/`method` from `uns['tissue_tool']` |
| The seven refusals | 7 | each wrong input still raises `ValueError`, matched on message text |

The property tests matter more than the equality tests. An equality test fails when a number
moves; a property test fails when the *science* breaks — if conformal coverage stopped reaching
its nominal level, `coverage 0.7883` could still pass while the method had become worthless.

One caveat of note 005 is pinned deliberately as a caveat:
`test_predict_scale_mismatch_is_still_present` asserts that measured counts remain more than 4×
the predicted mean. It **asserts the known defect, it does not endorse it** — SpaGE predictions
(mean 0.52) and measured counts (mean 2.58) are on different scales and TISSUE subtracts them
anyway (`tissue/main.py:652`). If someone fixes that upstream, this test fails loudly, which is
the correct outcome: the baseline would then need rewriting, not quietly inheriting.

## The pinned environment has no pytest

`tissue-env-pinned` has no `pytest`, and installing one is a download the shop order forbids.
So the file is a **normal pytest module that also runs standalone**, using no pytest-only
features (no decorators, no fixtures — a module-level cached `results()` instead):

```bash
cd ~/Work/paper2agent
TISSUE_Agent/tissue-env-pinned/bin/python TISSUE_Agent/tests/test_tissue_numeric_baseline.py
# -> 32 passed
```

Both paths were exercised. Under real pytest 9.1.1 it collects and passes as 32 tests —
**`32 passed, 37 warnings in 13.88s`**. To get pytest onto the 3.10 interpreter without
installing anything, I symlinked the five pure-Python packages it needs (`pytest`, `_pytest`,
`pluggy`, `iniconfig`, `packaging`, plus `pygments` and the `py.py` shim) out of
`Scanpy_Agent/scanpy-env` into a scratch directory on `PYTHONPATH`. **That shim is a throwaway,
not part of the project** — it proves the file is real pytest, and the standalone path is how
anyone should actually run it.

## What this baseline does *not* establish

The important limitation, stated plainly: **this is self-consistency, not correctness.**

- It pins note 005's numbers, which pinned note 004's, which came from running the tutorial.
  **Nobody has checked any of them against the TISSUE publication.** A number can be stable and
  wrong, and if the tutorial path itself is wrong, this test now defends the error.
- **One dataset, one gene, one method.** `plp1`, `spage`, the built-in 420 KB set. The `knn`
  method passes validation and has still never been run. Checking it needs no download — it is
  simply not done.
- **4 decimal places** — the precision note 005 printed, and therefore the strongest honest
  claim. The wrappers `round(..., 4)` on the way out, so finer drift is invisible here.
- It is welded to the pinned 3.10 stack. On a modern environment TISSUE does not import at all
  (`AnnData.__init__() got an unexpected keyword argument 'dtype'`, note 004), so this file
  cannot tell you whether the numbers would survive a dependency upgrade.
- The seven refusals are matched on message *fragments*. Rewording an error message fails the
  test — intended, since the messages are the wrappers' actual product, but it means a message
  improvement shows up as a failure to be read, not a bug.

## Two corrections to earlier notes

1. **Note 005 says the dataset is "3,405 cells × 32 genes". The spatial file has 33 gene
   columns.** Traced it: `Spatial_count.txt` carries 33 genes; after `load_paired_datasets`
   prevalence-filters the scRNA side to 32 and the lowercased intersection runs, **32 genes are
   shared** — spatial `flt1` is absent from the filtered scRNA set. Holding out `plp1` leaves the
   **31** calibration genes both notes report. So "32" was the shared count mislabelled as the
   file's gene count; every downstream number is unaffected.
2. **Note 005's claim that the wrappers write 23 MB per run is a disk leak, not a fact to
   preserve.** Each baseline run writes ~22 MB of `.h5ad` to a temporary directory. Note 010
   records **11 GB** lost to exactly this pattern of accumulating per-call output directories, so
   the test now deletes its own output on exit (`atexit`), with `TISSUE_KEEP_OUTPUT=1` to keep it
   for inspection. Three stray directories from developing this note were cleaned up.

## Verification of the test itself

A test that cannot fail is worse than no test. Negative control: I changed the expected
`pearson_r_vs_measured` from `0.3097` to `0.3100` — a 3 × 10⁻⁴ perturbation — and the suite
failed with the file and quantity named:

```
AssertionError: pearson_r_vs_measured drifted: got 0.309700, note 005 recorded 0.3100
                (difference 3.00e-04)
1 failed, 31 deselected
```

The file was then restored and confirmed byte-identical. The runner also distinguishes
`AssertionError` (a drifted number) from any other exception (a broken pipeline), because those
are different problems with different responses.

## Cost

- **Time:** 15.7 s standalone, 13.88 s under pytest. Matches note 005's 14 s.
- **RAM:** peak **407 MB** — identical to note 005. Nowhere near the 10 GB cap. No GPU.
- **Disk:** **one 16 KB source file.** No new environment, no new dataset, no downloads. Net disk
  change is negative: the cleanup fix reclaimed 66 MB of stray temporary output.

## Next

Open item 4 of note 010 is closed: TISSUE has reference values of its own, and they are checked
by something that runs.

Not in git yet. `.gitignore` admits `TISSUE_Agent/run_tutorial.py` and `src/*.py` only, so the
test needs a `!TISSUE_Agent/tests/*.py` exception before it can be committed — the same situation
`MANIFEST.md` and `USAGE.md` were in.

Still open from note 010: **Stage 5** (`requirements.txt` has never been installed), the
**TPM/CPM hole** in Scanpy QC, and the Scanpy **ZIP + relocation test**. The cheapest real
extension of *this* work is the `knn` method — allowed by the validation, never executed, and it
needs nothing downloaded.
