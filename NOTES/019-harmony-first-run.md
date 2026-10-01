# 019 — First Harmony run on the two bone-marrow lanes

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 cores, **no GPU**
**Environment:** `Scanpy_Agent/scanpy-env` (Python 3.12.14) only. `tissue-env-pinned` not touched.
**Nothing installed.** No pip, no `uv pip install`, no download. `NOTES/017` not regenerated.
**Result:** **It finished.** 3.8 s, peak **492 MB** — nowhere near the 8 GiB cap. Batches mixed
about **halfway** to perfect; cell-type structure **held**. And two things turned out not to be
what note 018 assumed — see "Two corrections".

> **This does not answer the paper question.** One run, one dataset, two lanes, one method, one
> parameter set, no comparison against any alternative. It establishes that the method runs here
> and produces a plausible result. That is all.

## Exact command

The run is a script rather than a one-liner, because the point was to avoid loading the
23,427-gene matrix. The call itself:

```python
sc.pp.harmony_integrate(adata, "sample", basis="X_pca",
                        adjusted_basis="X_pca_harmony", rng=0)
```

invoked as:

```bash
cd ~/Work/paper2agent
MPLBACKEND=Agg Scanpy_Agent/scanpy-env/bin/python <script>
```

`adata` was built by reading **only** `obsm/X_pca` and `obs` out of the reference files with
`h5py` + `anndata.io.read_elem`:

```python
with h5py.File(f"{REF}/05_pca.h5ad", "r") as f:
    x_pca = np.asarray(read_elem(f["obsm/X_pca"]), dtype=np.float32)   # (17041, 50)
    obs_pca = read_elem(f["obs"])                                      # has 'sample'
with h5py.File(f"{REF}/08_annotated.h5ad", "r") as f:
    obs_ann = read_elem(f["obs"])                                      # has the labels
assert list(obs_pca.index) == list(obs_ann.index)
```

**That assertion is why joining labels across two files is legitimate here** — both are 17,041 ×
23,427 and the check confirms identical cell order, so `08_annotated`'s labels describe
`05_pca`'s cells.

### On the 8 GiB cap

Never approached. Harmony corrects the **embedding**, not the expression matrix, so the working set
is 17,041 × 50 floats — about 3.4 MB. Reading `X` would have pulled in a 23,427-gene sparse matrix
for no benefit, so it was never read. **Peak RSS 492 MB**, most of it the scanpy/numpy import
itself.

## Input files and hashes

Both verified against `Scanpy_Agent/tests/MANIFEST.md` immediately before the run, and both match:

| File | SHA-256 | Shape | Role |
|---|---|---|---|
| `notebooks/clustering/ref/05_pca.h5ad` | `467ed3b75bcfba9a5413f0b0de48a3218694c0c9abe5d2c8767ffad81811dc9c` | 17041 × 23427 | `obsm['X_pca']` (17041×50) and `obs['sample']` |
| `notebooks/clustering/ref/08_annotated.h5ad` | `e5cc44f2a698e683d53c7eb9823b83f6cfe5f708521cc993c6b460a8787df84d` | 17041 × 23427 | `obs['cell_type_lvl1']`, `obs['leiden_res_*']` |

Both descend from the two figshare lanes already on disk and hashed
(`s1d1` `322c30a7…`, `s1d3` `ca40b287…`). **Nothing was downloaded.**

The two lanes, as they appear in `obs['sample']`:

| Lane | Cells |
|---|---|
| `s1d1` | 8,713 |
| `s1d3` | 8,328 |

(Note the cell counts differ from the raw files' 8,785 / 8,340 — QC and doublet filtering in notes
002's chain removed 84 cells before PCA.)

## Did it finish?

**Yes.** `sc.pp.harmony_integrate` returned normally in **3.8 s**, writing
`obsm['X_pca_harmony']` with shape `(17041, 50)`. No error, no warning that changed the outcome,
no OOM, no GPU involved.

## The numbers — both halves

Both halves use the same measurement so they are comparable: build a k=15 nearest-neighbour graph
on the embedding, then ask what fraction of each cell's neighbours share its label, averaged over
cells. Run on `X_pca` (before) and `X_pca_harmony` (after).

### Half 1 — did batches mix?

Label = the lane. **Lower is better**; the floor is what perfect mixing would give
(∑p², here 0.5003, since the lanes are nearly equal size).

| | same-lane neighbour fraction |
|---|---|
| Before (`X_pca`) | **0.6969** |
| After (`X_pca_harmony`) | **0.6029** |
| Perfect mixing | 0.5003 |

**47.8 % of the gap to perfect mixing was closed.** Before the run, a cell's neighbours were
69.7 % same-lane against a 50.0 % baseline — a clear batch effect. After, 60.3 %. Real movement in
the right direction, and **plainly incomplete**: residual lane structure remains.

### Half 2 — did cell-type structure hold?

Label = cell type. **Higher is better, and the test is that it does not drop.** Measured at three
granularities, because `cell_type_lvl1` has only 4 categories and 57 % of cells are one of them
(`Lymphocytes`, 9,763 of 17,041) — a near-1.0 score there is easy and weak evidence.

| Label | Categories | Before | After | Δ | Chance |
|---|---|---|---|---|---|
| `cell_type_lvl1` | 4 | 0.9955 | **0.9956** | **+0.0001** | 0.3960 |
| `leiden_res_0.50` | 17 | 0.9599 | **0.9585** | **−0.0013** | 0.1224 |
| `leiden_res_2.00` | 36 | 0.8749 | **0.8598** | **−0.0151** | 0.0470 |

**Structure held.** The coarse labels are untouched; the finest set loses 1.5 percentage points,
which is the honest cost and is visible only because the finer labels were checked. A single
`cell_type_lvl1` number would have reported "+0.0001, no loss at all" and been misleading — the
note 006 lesson applied deliberately.

Cell-type composition, for the record:

| `cell_type_lvl1` | Cells |
|---|---|
| Lymphocytes | 9,763 |
| Erythroid | 3,502 |
| B Cells | 2,269 |
| Monocytes | 1,507 |

**Caveat on the leiden rows:** `leiden_res_0.50` and `leiden_res_2.00` were computed on the
**uncorrected** PCA in note 002's chain. They are therefore partly a description of the batch
structure Harmony is trying to remove, so preserving them is a *differently biased* test, not a
stricter one. Some of that −0.0151 may be Harmony correctly merging clusters that were split by
lane. This note does not try to separate those two causes.

## Two corrections

### 1. `harmonypy` was not used. Note 018's install was not needed for this.

`scanpy` in this environment is `1.13.0a3.dev0+g0d5fd1623`, which ships **its own vendored
Harmony** at `scanpy/preprocessing/_harmony/`. And `sc.external.pp.harmony_integrate` — the path
note 018 checked — is a deprecated shim:

```python
@deprecated(Deprecation("1.13.0", "Import from sc.pp instead"))
def harmony_integrate(*args, **kwargs):
    from ...preprocessing import harmony_integrate
    return harmony_integrate(*args, **kwargs)
```

`grep` finds no `import harmonypy` anywhere in scanpy's harmony module or in `external/pp` — only
a docstring remark that the defaults follow harmonypy 2.0.0. **So the 11.6 MiB install recorded in
note 018 did not enable this run, and nothing in this note depended on it.** Note 018's claim that
`harmony_integrate` is "now backed by a real harmonypy" is **wrong for this environment** and
should be read with this correction.

Why the mistake was easy: note 013 established that in `tissue-env-pinned`'s **scanpy 1.9.8**,
`harmony_integrate` genuinely did raise `ImportError` without `harmonypy`. That remains true
*there*. This environment has a four-minor-version-newer scanpy where it is no longer true, and
the version difference was not checked before installing.

### 2. The default is not reproducible. `rng=0` is.

`sc.pp.harmony_integrate` takes `rng: SeedLike | RNGLike | None = None`. With the default, **two
calls in one process on identical input gave different embeddings — max absolute difference
0.641** — and the batch-mixing figure wandered across runs (0.6001, 0.5973, 0.6029). With `rng=0`
it is exact:

```
reproducible with rng=0: True (max abs diff 0)
```

Every number in this note is from the `rng=0` run. **A future ticket that pins these figures the
way note 012 pinned TISSUE's must pass `rng` explicitly**, or it will produce a flaky test. For
the separate `harmonypy` package, `run_harmony` defaults to `random_state=0` and was verified
deterministic at both `ncores=1` and `ncores=0` — but it is not what ran here.

## What failed

Nothing failed. Three things went differently than planned and are recorded rather than tidied:

- The first two runs used the **deprecated** `sc.external.pp` path and the **default rng**, so
  their numbers were not reproducible. They are superseded by the `rng=0` run above and are quoted
  only as evidence of the non-determinism.
- `cell_type_lvl1` alone was too coarse to be worth reporting as "structure held". Caught by
  checking three granularities, not by the first number looking good.
- The premise inherited from note 018 — that `harmonypy` was the enabling dependency — was false
  in this environment, and only surfaced because the source was read when the determinism question
  came up.

## Caveats

- **One run, one dataset, two lanes, one method.** No alternative was run; nothing is compared.
- **Default parameters throughout** except `rng=0`. `theta`, `sigma`, `n_clusters`,
  `max_iter_harmony` were all left at their defaults, and no sensitivity check was done.
- **The metric is one of many.** Same-label kNN fraction at k=15 is a reasonable, transparent proxy
  for both halves, but it is not kBET, iLISI, ARI or silhouette, and the ranking of methods can
  differ by metric. k was not varied.
- **No biology is claimed.** These are two teaching lanes of bone marrow. Nothing here is evidence
  about haematopoiesis.
- **The corrected embedding was not saved into the project** — it exists only in a scratch file.
  Nothing in `notebooks/`, `tests/` or `MANIFEST.md` changed.
- **`cell_type_lvl1` is itself a product of note 002's annotation step**, not ground truth. Both
  halves of the measurement rest on labels this project produced.

## Next

**Still not answered:** whether Harmony is the right choice, how it compares to anything else,
whether 47.8 % gap closure is good, and whether the finest-grained 1.5-point structure loss
matters. Those need more than one run.

The cheapest honest follow-ups, in order: vary `rng` to see how much the figures move run to run;
vary `theta` (Harmony's own mixing-strength knob) to see the mixing/structure trade-off this single
point sits on; and only then consider whether this belongs as a tenth Scanpy tool.

Unchanged and still open: **Stage 5** (122 pins vs 178 installed packages), the **TPM/CPM hole** in
Scanpy QC, the Scanpy **ZIP + relocation test**, and note 013's **`n_neighbors` passthrough** for
TISSUE's `knn` — which this note does not touch, since TISSUE runs in the other environment.
