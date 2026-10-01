# 013 — TISSUE `knn`: two blockers, and a refusal that was lying

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Builds on:** [005](005-tissue-tools.md), [012](012-tissue-numeric-baseline.md) — same checkout,
same `tissue-env-pinned`, **nothing installed, nothing downloaded**
**Result:** **`knn` cannot run here, and the ticket to pin its numbers was stopped.** The
refusal message that claimed otherwise is fixed. **37 assertions, 37 passed.**

## Short version

The ticket was to run `knn` on the same three on-disk files note 012 used and pin its numbers
the same way. It did not get that far, on purpose: `knn` is blocked twice over, and one blocker
needs a package that is not in the pinned environment. Per the standing order, nothing was
installed.

The more useful finding is what the attempt exposed. `tissue_tools.py` declared

```python
SUPPORTED_METHODS = ("spage", "knn")
```

and refused anything else with *"method 'tangram' **is not available in this environment**;
supported: ['spage', 'knn']"*. **That message was false.** `knn` was never runnable in this
environment, and the guard asserted it was. Note 005 recorded that text as one of its seven
clean refusals, and note 012's `test_refuses_unavailable_method` pinned the wording — so the
baseline was certifying a false claim. That is the defect this note closes.

## The two blockers

### 1. No way to pass `n_neighbors` — ours, and still unfixed

```
TypeError: knn_impute() missing 1 required positional argument: 'n_neighbors'
```

- `tissue.main.knn_impute(spatial_adata, RNAseq_adata, genes_to_predict, n_neighbors, **kwargs)`
  (`main.py:384`) — `n_neighbors` is **positional with no default**.
- `predict_gene_expression` forwards only `**kwargs` to it (`main.py:349`).
- `tissue_predict_spatial_gene` accepts exactly `target_gene, method, n_folds, n_pv,
  output_dir` — **no passthrough**. There is no argument a caller can set to reach
  `n_neighbors`.

A related wart: the wrapper hard-passes `n_pv=10`, which is a **SpaGE** parameter. For `knn` it
lands in `knn_impute`'s `**kwargs` and is silently ignored, so the one thing the wrapper does
forward is the wrong one for this method.

This blocker is in our code, not upstream's, and **it is not fixed here** — fixing it would
mean adding a parameter that still cannot be exercised, because of blocker 2. Adding untestable
code to satisfy a ticket is how the note-005 situation arose in the first place.

### 2. `harmonypy` is not installed — needs a download

Even with `n_neighbors` supplied, `knn_impute` calls Harmony unconditionally:

```python
# main.py:403-404
sc.tl.pca(joint_adata)
harmony_integrate(joint_adata, 'batch', verbose=False)
```

and scanpy's implementation raises on demand:

```python
# scanpy/external/pp/_harmony_integrate.py:78-80
import harmonypy
...
raise ImportError("\nplease install harmonypy:\n\n\tpip install harmonypy")
```

`harmonypy` is absent (`ModuleNotFoundError`); scanpy 1.9.8 is present. The trap is that
`from scanpy.external.pp import harmony_integrate` **succeeds** — it is a lazy wrapper that
only fails when called. So an availability check based on importing it would have reported
`knn` as fine, which is probably how the original claim got made.

## The fix: compute availability, do not assert it

`SUPPORTED_METHODS` was a hand-maintained tuple, and a hand-maintained claim about the
environment goes stale silently. It is replaced by a probe:

- `KNOWN_METHODS = ("spage", "knn", "gimvi", "tangram")` — what **TISSUE implements**, with an
  explicit comment that membership does *not* imply usable.
- `_method_blocker(method)` — returns **why** a method cannot run here, or `None` if it can.
  For `knn` it always reports the missing `n_neighbors` passthrough, and appends the
  `harmonypy` reason only when `importlib.util.find_spec("harmonypy") is None`.
- `available_methods()` — `[m for m in KNOWN_METHODS if _method_blocker(m) is None]`, which in
  this environment is exactly `['spage']`.

The refusals now read:

```
method 'knn' cannot be run here: tissue_predict_spatial_gene has no way to pass the
'n_neighbors' argument that tissue.main.knn_impute requires (it is positional with no
default); and harmonypy is not installed, so the harmony_integrate call inside knn_impute
raises ImportError. Runnable methods in this environment: ['spage'].

method 'tangram' cannot be run here: tangram needs extra packages that the repo's
requirements.txt leaves commented out, and they are not installed. Runnable methods in
this environment: ['spage'].

method 'notamethod' is not a TISSUE prediction method; TISSUE implements
['spage', 'knn', 'gimvi', 'tangram'] and this environment can run ['spage'].
```

Three properties worth stating. The message **names both blockers**, not just the first one hit
— a refusal that hides one of two problems sends the reader to fix the wrong thing. It
**distinguishes "not a TISSUE method" from "a real method this environment cannot run"**, which
the old single branch conflated. And if `harmonypy` is ever installed, the message **narrows by
itself** to the passthrough blocker, instead of going stale.

## Tests

`TISSUE_Agent/tests/test_tissue_numeric_baseline.py` goes from 32 to **37 tests**. The 32
numeric assertions from note 012 are untouched and still pass, which is the evidence that this
refactor changed no science. Five added:

| Test | What it pins |
|---|---|
| `test_refuses_unavailable_method` | **updated** — pinned the old false `"is not available in this environment"`; now `"cannot be run here"` |
| `test_refuses_unknown_method` | a bogus method is reported as not-a-TISSUE-method, a different failure from an unavailable one |
| `test_knn_is_not_advertised_as_runnable` | `available_methods() == ["spage"]` and `knn` is not in it — the actual bug, asserted directly |
| `test_knn_refusal_names_both_blockers` | the message contains **both** `n_neighbors` and `harmonypy` |
| `test_knn_refusal_is_a_valueerror_not_a_typeerror` | `knn` is refused by the guard, not by a `TypeError` from inside upstream |
| `test_harmonypy_really_is_absent` | the premise. If harmonypy appears, the *new* message becomes the false one and this fails first |

That last one matters: every honest claim about an environment should fail loudly when the
environment changes under it. The old code had no such test, which is why its claim rotted
unnoticed for eight days.

```
37 passed          # standalone, pinned env
37 passed, 37 warnings in 13.92s   # under pytest 9.1.1
```

## What this note does not do

- **`knn` is still unrun.** No numbers for it exist, here or anywhere. The ticket to pin them is
  not complete — it is blocked, and the block is recorded rather than worked around.
- **The `n_neighbors` passthrough is not added.** It is a two-line change, but it cannot be
  tested without `harmonypy`, so it waits for a ticket that permits the install.
- **No judgement on whether `knn` is worth running.** It routes through Harmony integration,
  so its numbers would not be comparable to SpaGE's on any axis except final coverage.
- Nothing was installed, nothing downloaded, no new environment.

## Corrections to note 012

Note 012 reports "32 assertions, 32 passed" and lists "The seven refusals". Both are now stale:
the file holds **37** tests, and one of the seven refusals it recorded was pinning a false
message. The numeric baseline it established — every value in its tables — is unaffected and
still verified.

## Next

Open item 4 of note 010 stays closed. Still open: **Stage 5** (`requirements.txt` has never been
installed), the **TPM/CPM hole** in Scanpy QC, and the Scanpy **ZIP + relocation test**.

For TISSUE specifically, the honest options are now: install `harmonypy` (one small pure-Python
package) and finish the `knn` ticket properly, including the passthrough — or leave `knn` alone
and say so in the record, which is what this note does.
