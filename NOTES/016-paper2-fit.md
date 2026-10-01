# 016 — What paper 2 could be, and what it would cost

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**, 237 GB disk
**Method:** reading what is already on disk. **No install, no download, no clone, no conversion run.**
**Result:** **There is no fourth hosted MCP to probe.** One candidate paper is the best fit on this
disk — **Harmony** — but the check that decided its cost came back `None`: `harmonypy` is in neither
environment. **Refused until the human explicitly allows one pip install into `scanpy-env`.**

## The hosted route is exhausted

Paper2Agent's own README lists exactly three connectable MCP servers:

| Server | Status in these notes |
|---|---|
| `Paper2Agent-scanpy-mcp.hf.space` | **Done.** Converted locally in 002; 9 tools, reused in 003 and 006, guarded in 007–008, tested in 009. |
| `Paper2Agent-tissue-mcp.hf.space` | **Done.** Tutorial in 004, 2 wrappers in 005, numeric baseline in 012, `knn` refused in 013. |
| `Paper2Agent-alphagenome-mcp.hf.space` | **Ruled out.** 015: live and 22 tools, but the 3 that run the model need a Google DeepMind API key, and artifacts 404. |

Its README names the same three as its only worked examples, and its demos section shows the same
three. **So "pick a fourth hosted MCP" is not an available move** — there is no fourth. Paper 2 has
to be a paper wrapped from scratch, which makes the environment question the whole decision.

## Why the environment question decides everything

This laptop holds two environments and they do not mix:

- `Scanpy_Agent/scanpy-env` — Python **3.12**, 868 MB, 123 pins. Note 002/010: **never reinstalled
  from `requirements.txt`**; Stage 5 is still open.
- `TISSUE_Agent/tissue-env-pinned` — Python **3.10**, 693 MB. Note 004: a modern stack breaks TISSUE
  outright (`AnnData.__init__() got an unexpected keyword argument 'dtype'`), and it needs
  `setuptools<81`.

Note 010 prices a fresh environment at **~900 MB of downloads** for the Scanpy side and ~700 MB for
TISSUE. Disk is not the constraint — 237 GB free, the whole factory is 11 GB. **The constraint is
the standing order: no downloads without an explicit ticket.** So a candidate paper is cheap only
if its dependencies already sit inside an environment that exists.

That test eliminates most of the obvious candidates before any science is considered: anything on
PyTorch (scvi-tools, most deep-learning methods) means a new multi-gigabyte environment and a GPU
this machine does not have. Those are not close calls.

## The candidate: Harmony

**Paper:** Korsunsky et al., *Fast, sensitive and accurate integration of single-cell data with
Harmony*, **Nature Methods** 16, 1289–1296 (2019) —
<https://www.nature.com/articles/s41592-019-0619-0>

> Confirm the DOI and the Python port's repository URL before acting on this note. I did not fetch
> either; they are written from memory and this note is not evidence for them. The rest of this
> section rests on files on this disk.

### Why it is next

Harmony is the one candidate where the hardest part of every previous note — getting data and an
environment — is **already solved on this disk**. It is a batch-integration method: given cells
from two experiments, it corrects the embedding so the same cell type from different batches lands
together. The two bone-marrow lanes the whole Scanpy half was built on, `s1d1` and `s1d3`, are
*exactly* that problem — two samples, already concatenated into `00_raw_concat.h5ad` with a
`sample` column, already hashed in `tests/MANIFEST.md`, already carried through PCA in
`ref/05_pca.h5ad`. Note 002 ran the tutorial's batch-aware steps but **never corrected the batch
effect**, so this is the natural next question on data that needs no download, with a published
method, a small pure-Python implementation, and a scanpy entry point (`sc.external.pp.
harmony_integrate`) that the nine existing tools already sit beside. It also has a reason beyond
convenience: note 013 hit `harmonypy` as a missing dependency of TISSUE's `knn` path, so knowing
whether it is installable and behaves here answers two open questions at once.

### What would run on this machine

- **CPU only, no GPU.** Harmony is iterative clustering plus linear correction on a PCA embedding —
  no neural network, no CUDA.
- **The input is already here.** `ref/05_pca.h5ad` is 137 MB and holds `X_pca` for 17,125 cells ×
  36,601 genes across the two lanes. Harmony corrects the **embedding**, not the full matrix, so the
  working set is roughly 17,125 × 50 floats — a few MB. Nothing resembling the 6.4 GB peak of the
  current pytest suite.
- **Nothing to download.** Zero bytes. The 43 MB of figshare lanes at the root of the chain are
  already on disk and hashed.
- **A real before/after exists.** Note 006's silent-wrong-answer lesson applies directly: the honest
  check is whether cells mix *and* whether cell-type structure survives, not just whether the UMAP
  looks prettier.

### What would not run here

- **Any benchmark from the paper.** The publication evaluates across several large datasets; those
  are hundreds of MB to GB and are out under the 500 MB download cap.
- **Comparisons against the alternatives** (Seurat CCA, scVI, LIGER, BBKNN). Each is its own
  environment, and the deep-learning ones want a GPU. A single-method run is all this machine
  supports.
- **The R implementation.** There is no R toolchain in this project, and adding one is a new
  environment by any definition.
- **Any claim about biology.** Same limit as every other note: these are teaching files, two lanes,
  one tissue.

### Do we need a new pip/conda env?

**No new environment — but possibly one new package, and that is still a download.**

The honest position, split so neither half is overstated:

1. **No new env is required.** Harmony's Python port is pure Python and the target interpreter is
   the existing `scanpy-env` 3.12, which already has scanpy, anndata and numpy. This is the only
   candidate I found that clears that bar.
2. **Measured: 2026-10-01 check printed `None` — `harmonypy` is absent from `scanpy-env`.**
   Together with note 013, which found it absent from `tissue-env-pinned`, it is now confirmed
   missing from **both** environments on this machine.
3. There were two possible outcomes, and the measurement landed on the one that is not free:
   - ~~**Already present** → no install, no download, no new env.~~ **Did not happen.**
   - **Absent** ← **this is the case.** One small pure-Python package. Not a new environment, but
     **a download**, and it would mutate `scanpy-env`, whose pins have never been proven
     reproducible (Stage 5 open) — so the mutation would be unrecoverable by any tested path.

**Verdict: paper 2 via Harmony is refused until the human explicitly allows one pip install into
`scanpy-env`.** Not because the paper is a bad fit — it is still the best fit on this disk — but
because the measurement came back `None`, and the only route forward is a download into an
environment that has never been rebuilt from its own pins. That is a decision for the human, not a
detail to absorb quietly.

## First command block, for later — not run now

```bash
cd ~/Work/paper2agent && . ./env.sh                                    # project-local uv/python/caches
Scanpy_Agent/scanpy-env/bin/python -c "import importlib.util as u; print(u.find_spec('harmonypy'))"
grep -n "s1d1\|05_pca" Scanpy_Agent/tests/MANIFEST.md                  # the inputs, with their hashes
(cd Scanpy_Agent && ./scanpy-env/bin/python -m pytest tests/code/test_manifest_hashes.py -q)
# line 2 printed None -> STOP, that is a download. Printed a spec -> ticket may proceed.
```

Line 2 is the decision. Line 4 proves the baseline is intact *before* anything is touched, so a
later failure can be attributed. Nothing above installs, downloads or writes.

## Caveats

- **This note chooses nothing.** It is a fit assessment; `SHOP.md` reserves picking work for the
  human, and no paper has been started.
- **The paper link and repository URL are from memory, not fetched.** Confirm both first.
- **I did not survey exhaustively.** I eliminated PyTorch-based methods on the environment rule and
  stopped at the first candidate that cleared it. There may be a better fit I did not reach.
- **The `harmonypy`-in-`scanpy-env` check has now been run** and printed `None`. It was the one
  unresolved fact in the first draft of this note; the verdict above reflects the result, not a
  guess.
- No conversion was attempted, and per note 015 none should be attempted for a paper whose tools
  cannot be executed here.

## Next

Unchanged and still ahead of any new paper: **Stage 5** (`requirements.txt` has never been
installed — and note that installing `harmonypy` into `scanpy-env` would make Stage 5 *harder*, by
moving the environment further from its untested pins), the **TPM/CPM hole** in Scanpy QC, and the
Scanpy **ZIP + relocation test**.

If paper 2 is wanted sooner, run line 2 of the block above and tell me what it printed. That one
line decides whether this is a free ticket or a refused one.
