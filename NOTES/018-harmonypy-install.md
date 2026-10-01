# 018 — Installing harmonypy into scanpy-env

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Environment touched:** `Scanpy_Agent/scanpy-env` (Python 3.12.14) — **mutated, deliberately,
with permission**. No fresh environment, no Stage 5, no pip.
**Result:** **One package added, nothing else moved.** 177 → 178. The diff against note 017 is a
single line.

> **Corrected by [019](019-harmony-first-run.md): this install was not needed.** `scanpy` in this
> environment is `1.13.0a3.dev0`, which ships **its own vendored Harmony**
> (`scanpy/preprocessing/_harmony/`), and `sc.external.pp.harmony_integrate` is a deprecated shim
> forwarding to `sc.pp.harmony_integrate`. No `import harmonypy` exists anywhere in scanpy's
> harmony module. **The Harmony run in note 019 never touched `harmonypy`.** The *measurements*
> below — the one-line diff, the 177 → 178 count, the unchanged versions — are all still accurate;
> what is wrong is this note's conclusion that the install enabled anything. Note 013's
> `ImportError` was real, but in `tissue-env-pinned`'s **scanpy 1.9.8**, four minor versions older.
> That version difference was not checked before installing.

## Short version

Note 016 refused paper 2 because `harmonypy` was absent from both environments and adding it meant
a download into an environment whose pins have never been proven to rebuild it. That refusal was
lifted by explicit instruction. This note records what the install actually cost, measured against
the rollback map note 017 took for exactly this purpose.

It cost one line. That is the good outcome, and it was checked rather than assumed.

## The install

```
$ uv pip install harmonypy --python Scanpy_Agent/scanpy-env/bin/python
Using Python 3.12.14 environment at: Scanpy_Agent/scanpy-env
Resolved 2 packages in 341ms
Downloading harmonypy (11.6MiB)
 Downloaded harmonypy
Prepared 1 package in 2.31s
Installed 1 package in 41ms
 + harmonypy==2.0.2
```

`uv` only — `pip` is not installed in this venv and was not installed to do this (note 017).
Download: **11.6 MiB**, well inside the 500 MB cap.

## Pre-flight, so the diff means something

Before installing, the environment was confirmed byte-identical to the recorded map:

```
$ uv pip freeze --python Scanpy_Agent/scanpy-env/bin/python | diff - NOTES/017-scanpy-env-pins.txt
(no output — 177 packages, exact match)
$ Scanpy_Agent/scanpy-env/bin/python -c "import importlib.util as u; print(u.find_spec('harmonypy'))"
None
```

Without that step the diff below would be suggestive rather than conclusive: any drift could have
predated the install. It did not.

## The diff

Against `NOTES/017-scanpy-env-pins.txt` (the 177-package pre-install map):

```diff
46a47
> harmonypy==2.0.2
```

**One line added. Nothing removed. No version of any existing package changed.** 177 → 178.

Frozen both ways again, as in note 017, because `pip freeze` is unavailable here:

| Method | Before | After |
|---|---|---|
| `uv pip freeze` (uv 0.12.17) | 177 | **178** |
| stdlib `importlib.metadata.distributions()` | 177 | **178** |

The stdlib cross-check reports the same single addition. Its only other apparent difference is the
`scanpy` entry — `scanpy @ file:///…/repo/scanpy` from uv versus
`scanpy==1.13.0a3.dev0+g0d5fd1623` from `importlib.metadata` — which is the representation
difference already documented in note 017 point 1, **not** a change to the environment.

## Why nothing else moved

Not luck. `harmonypy==2.0.2` declares one runtime dependency, and it was already present:

```
requires: numpy
requires: pandas; extra == "test"
requires: pytest>=8.4.2; extra == "test"
requires: scipy; extra == "test"
```

The other three are `test` extras, which a plain install does not pull. So **there were no new
direct dependencies to account for** — the stop condition in the ticket ("anything besides
harmonypy and its new direct deps") had nothing to trigger on.

## Nothing broke

A mutated environment is only safe if what was already there still works. Checked:

| | Version | Status |
|---|---|---|
| `numpy` | 2.5.3 | **unchanged** |
| `anndata` | 0.13.4 | **unchanged** |
| `scanpy` | `1.13.0a3.dev0+g0d5fd1623` | **unchanged** — still the `0d5fd16` path install |
| `harmonypy` | 2.0.2 | imports |
| `scanpy.external.pp.harmony_integrate` | — | importable — ~~**and now backed by a real harmonypy**~~ **false, see [019](019-harmony-first-run.md)** |

That last row was the one note 013 warned about: `harmony_integrate` *imported* fine even when
`harmonypy` was missing, because in scanpy 1.9.8 it is a lazy wrapper that only raises when called.
Its importability was never evidence of anything — **and I then drew exactly the wrong conclusion
from it.** In *this* environment's scanpy 1.13.0a3 the function is not harmonypy-backed at all, so
"now the dependency behind it is actually present" was wrong twice over: the dependency was never
behind it. Note 019 read the source and settled it. The honest version of this row is: the import
succeeds, it succeeded before the install too, and it tells you nothing.

## What this changes, and what it does not

**Resolved:** note 013's blocker 2 for TISSUE's `knn` — the missing `harmonypy` — in *this*
environment. **In practice this resolved nothing**, since TISSUE does not run in this environment
(next bullet) and scanpy's Harmony does not use `harmonypy` at all ([019](019-harmony-first-run.md)).
The package is installed and, as of note 019, still unused by anything.

**Not resolved, and worth being precise about:**

- **TISSUE's `knn` is unaffected.** TISSUE runs in `tissue-env-pinned` (Python **3.10**), a
  different environment, where `harmonypy` is still absent. `available_methods()` there will
  continue to exclude `knn`, correctly. Nothing in note 013 needs changing.
- **Blocker 1 stands.** `tissue_predict_spatial_gene` still has no way to pass the `n_neighbors`
  argument `knn_impute` requires. That is our code, and it is still not fixed.
- **No Harmony run.** Nothing was integrated, no batch correction was performed, no numbers exist.
  The paper 2 question from note 016 is *unblocked*, not answered. **Superseded:**
  [019](019-harmony-first-run.md) ran it — 47.8 % of the batch gap closed, structure held — and it
  still does not answer the paper question.
- **Stage 5 got harder.** `requirements.txt` holds **122** pins; the environment now holds **178**
  packages. Note 017 measured that gap at 122 vs 177; this install widened it by one, and the new
  package is not in the pin file at all. Anyone later attempting the clean-environment install will
  not reproduce this environment, and now has one more reason why.

## Rollback

`NOTES/017-scanpy-env-pins.txt` remains the **pre-install** map — that is the point of it, and it
was deliberately not regenerated. To undo:

```bash
cd ~/Work/paper2agent && . ./env.sh
uv pip uninstall harmonypy --python Scanpy_Agent/scanpy-env/bin/python
uv pip freeze --python Scanpy_Agent/scanpy-env/bin/python | diff - NOTES/017-scanpy-env-pins.txt
# no output = back to the 177-package state note 017 recorded
```

Untested — written from the install's inverse, not run. The diff on line 3 is what would confirm it.

## Caveats

- **`harmonypy` has not been called**, only imported. Import success is weak evidence; note 013 is
  the reason to say so explicitly.
- **The post-install freeze was not written into the project.** It lives in a scratch file only, so
  there is currently no committed record of the 178-package state — just this note's diff.
- **No hashes.** Same limit as note 017: `harmonypy==2.0.2` names a version, not the 11.6 MiB of
  bytes that arrived.
- **One download was made**, by permission. It is the first package install into `scanpy-env` since
  the environment was built in note 002.

## Next

**Paper 2 via Harmony is now unblocked but unrun.** The honest next step is the smallest real one:
run `harmony_integrate` on `ref/05_pca.h5ad` — the two bone-marrow lanes already on disk and
already hashed in `tests/MANIFEST.md` — and check both halves, that batches mix *and* that
cell-type structure survives. Note 006 is the reason for the second half: a result that looks
better is not the same as a result that is right.

Still open and unchanged: **Stage 5** (now 122 pins against 178 packages), the **TPM/CPM hole** in
Scanpy QC, and the Scanpy **ZIP + relocation test**.
