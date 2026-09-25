# 010 — Operator sheet (Phase 1, note 10 of 10 — last)

**Date:** 2026-09-25
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**, 237 GB disk
**What this is:** not a new paper. How to run what notes 001–009 built, what is authoritative, what will fail here, and what is safe to delete.

**Verified cold today:** both toolsets import and run in their existing environments. **No new venv was created. Nothing downloaded. `Scanpy_Agent/tests/results` deleted — see [Disk](#disk).**

---

## Cold-start check (run today, 16 s total)

| Toolset | Command result | Time | Peak RAM |
|---|---|---|---|
| Scanpy | **9 tools exposed**; `scanpy_compute_qc_and_filter` on pbmc3k → `kept 2700/2700 cells and 13714/32738 genes` | 4 s | 352 MB |
| TISSUE | `tissue_predict_spatial_gene` → `Predicted 'plp1' for 3405 cells`; `tissue_calibrate_prediction_intervals` → `coverage 0.7883 over 31 calibration genes` | 12 s | 643 MB |

TISSUE's 0.7883 is identical to notes 004 and 005, so the chain still reproduces its published number from a cold process.

---

## Running the Scanpy 9 tools

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

cap10g "$PROJECT_PYTHON" -c "
import sys, asyncio; sys.path.insert(0,'src')
from tools.clustering import clustering_mcp
from fastmcp import Client
async def go():
    async with Client(clustering_mcp) as c:
        r = await c.call_tool('scanpy_compute_qc_and_filter',
            {'data_path': 'tmp/pbmc3k/00_raw.h5ad', 'output_dir': 'tmp/scratch'})
        print(r.data['message'])
asyncio.run(go())"
```

**The nine, in pipeline order.** Each reads an `.h5ad`, writes a new one plus figures, so they chain on `artifacts[0]['path']`:

`scanpy_compute_qc_and_filter` → `scanpy_detect_doublets_scrublet` → `scanpy_normalize_log1p` → `scanpy_select_highly_variable_genes` → `scanpy_run_pca` → `scanpy_build_neighbors_umap` → `scanpy_cluster_leiden` → `scanpy_plot_marker_gene_dotplot` → `scanpy_rank_marker_genes`

**As a real MCP server** (not registered; this is the command):

```bash
cd ~/Work/paper2agent && claude mcp add --scope project scanpy -- \
  "$PROJECT_PYTHON" "$PROJECT_ROOT/src/scanpy_mcp.py"
```

**Full worked examples:** `tmp/pbmc3k/run_pbmc3k.py` (note 003) and `tmp/pbmc68k/run_pbmc68k.py` (note 006). Both are drivers, not committed source.

## Running the TISSUE 2 tools

Different environment — Python 3.10, **not** the Scanpy one. They cannot share.

```bash
cd ~/Work/paper2agent && . ./env.sh && export TISSUE_ROOT="$P2A/TISSUE_Agent"

PYTHONPATH="$TISSUE_ROOT/src" XDG_CACHE_HOME="$TISSUE_ROOT/tmp/cache" \
NUMBA_CACHE_DIR="$TISSUE_ROOT/tmp/cache/numba" MPLCONFIGDIR="$TISSUE_ROOT/tmp/mpl" \
MPLBACKEND=Agg cap10g "$TISSUE_ROOT/tissue-env-pinned/bin/python" -c "
from tissue_tools import tissue_predict_spatial_gene, tissue_calibrate_prediction_intervals, _tissue_repo
d = _tissue_repo()/'tests'/'data'
p = tissue_predict_spatial_gene(
    spatial_counts_path=str(d/'Spatial_count.txt'), locations_path=str(d/'Locations.txt'),
    scrna_counts_path=str(d/'scRNA_count.txt'), n_folds=3, output_dir='/tmp/t')
print(p['message'])
print(tissue_calibrate_prediction_intervals(data_path=p['artifacts'][0]['path'], output_dir='/tmp/t')['message'])"
```

`tissue_tools.py` locates the TISSUE checkout itself (`$TISSUE_REPO`, else `TISSUE_Agent/repo/TISSUE`), so no `cd` into the repo and no `PYTHONPATH=.` — that was only needed for note 004's driver.

**The four env vars are not optional decoration:** without them numba, matplotlib and pooch write caches into `$HOME`, which breaks the "everything under ~/Work/paper2agent" rule.

---

## Source of truth

**Authoritative — in git ([github.com/dryow1/paper2agent-notes](https://github.com/dryow1/paper2agent-notes)):**

| File | What it is |
|---|---|
| `Scanpy_Agent/src/tools/clustering.py` | the 9 Scanpy tools, 615 lines, including the three guards |
| `Scanpy_Agent/src/scanpy_mcp.py` | the MCP server entry point |
| `Scanpy_Agent/tests/code/clustering/test_tool_guards.py` | 14 guard tests (note 009) |
| `TISSUE_Agent/src/tissue_tools.py` | the 2 TISSUE wrappers |
| `TISSUE_Agent/src/run_tissue_tools.py` | their runner + 7 error cases |
| `TISSUE_Agent/run_tutorial.py` | note 004's tutorial driver |
| `NOTES/001`–`010` | the record |

**Authoritative but *not* in git** (too large; lost if the working tree goes):

| Path | Size | Why it matters |
|---|---|---|
| `Scanpy_Agent/notebooks/clustering/ref/*.h5ad` | part of 2.3 GB | the reference checkpoints every numeric test compares against |
| `Scanpy_Agent/tests/code/clustering/` (other 10 files) | 316 KB | the 50 original tests + `clustering_verify_helpers.py` |
| `Scanpy_Agent/reports/` | 264 KB | scanner/executor/verification/acceptance records |
| `Scanpy_Agent/tmp/build/requirements.txt` | 123 pins | the compiled environment for a clean install |

**This is the biggest structural weakness of the factory.** The committed test file cannot run from a clone: it needs `conftest.py`, `clustering_verify_helpers.py` and the reference `.h5ad` files, none of which are in git. The repo holds the *record and the source*, not a runnable checkout.

---

## What will fail on this machine

| Thing | What happens | Evidence |
|---|---|---|
| **GPU papers** | No NVIDIA GPU, `nvidia-smi` absent. Nothing CUDA/rapids was installed, by rule. Any paper whose tutorial requires a GPU is out. | 001 |
| **Huge spatial dumps** | 500 MB download cap. TISSUE's figures-and-analyses repo was never cloned; only the 420 KB built-in set was used. | 004 |
| **Current scientific Python for TISSUE** | `TypeError: AnnData.__init__() got an unexpected keyword argument 'dtype'` — removed in anndata 0.11. TISSUE needs the pinned 3.10 stack. | 004 |
| **TISSUE in a modern venv** | `ModuleNotFoundError: No module named 'pkg_resources'` — old numba needs `setuptools<81`, undocumented upstream. | 004 |
| **Scaled or log-normalised data into QC** | Now a clear `ValueError`. Before note 007 it silently returned negative `total_counts` for 389 of 700 cells. | 006, 007, 008 |
| **TPM/CPM into QC** | **Still accepted silently.** Indistinguishable from ambient-corrected counts, which must be accepted. Known, deliberate, unfixed. | 008 |
| **Marker panels from another tissue** | Rejected by default; pass `ignore_missing_genes=True`. pbmc3k drops 8 genes, pbmc68k_reduced drops 34 and empties 2 sets. | 003, 006, 007 |
| **Datasets with no batch column** | `batch_key` defaults to `None`, so they just work. Passing `"sample"` gives a clear error naming the available columns. | 003, 006 |
| **`cap10g` + background execution** | `systemd-run --user --scope` dies with the launching shell; pytest got killed after one test, leaving a 1-byte log. Run heavy jobs in the foreground. | 007 |
| **Peak memory** | Heaviest real job is the full pytest suite at **6.4 GB**. Fits 15 GiB with room, but not alongside much else. | 009 |

---

## Disk

**Cleaned today: `Scanpy_Agent/tests/results` deleted — free space went 195 → 204 GiB, so ~9 GiB reclaimed.** (`du` had reported the directory as 11 GB; the gap is filesystem accounting.) It held 73 directories of accumulated tool outputs: every tool call writes a fresh timestamped subdirectory, so each suite run adds rather than overwrites, and notes 007–009 ran the suite four times. The suite recreates what it needs.

**Factory total now: 11 GB** (was 22 GB), of 237 GB.

| Path | Size | Safe to delete? |
|---|---|---|
| `Scanpy_Agent/tmp/outputs` | 3.9 G | **Yes** — note 002's acceptance outputs; the reports summarising them are kept |
| `Scanpy_Agent/notebooks/clustering/replay` | 1.1 G | **Yes** — the fresh-process replay, already verified bitwise identical |
| `tools/uv-cache` | 2.2 G | **Yes**, but then any reinstall re-downloads |
| `Scanpy_Agent/tests/data` | 343 M | **Yes** — `sub4_*` subsamples, rebuilt on demand by `subsample()` |
| `Scanpy_Agent/tmp/pbmc3k`, `tmp/pbmc68k`, `tmp/guards` | 325 M | **Yes** — notes 003/006/008 run outputs |
| `TISSUE_Agent/out` | 23 M | **Yes** — note 005 outputs |
| `Scanpy_Agent/notebooks/clustering` (rest) | 1.2 G | **No** — holds `ref/*.h5ad`, the reference every test compares against |
| `Scanpy_Agent/scanpy-env` | 868 M | **No** — rebuilding costs ~900 MB of downloads |
| `TISSUE_Agent/tissue-env-pinned` | 693 M | **No** — same, and the pinned stack is fiddly to rebuild |
| `Scanpy_Agent/repo`, `TISSUE_Agent/repo` | 105 M | No — pinned checkouts, cheap to keep |

One command for the safe set (**frees ~7.5 GB**):

```bash
cd ~/Work/paper2agent && rm -rf \
  Scanpy_Agent/tmp/outputs Scanpy_Agent/notebooks/clustering/replay \
  Scanpy_Agent/tests/data Scanpy_Agent/tmp/pbmc3k Scanpy_Agent/tmp/pbmc68k \
  Scanpy_Agent/tmp/guards TISSUE_Agent/out tools/uv-cache
```

That would leave the factory at ~3.5 GB: two environments, two checkouts, the reference data, the source and the notes.

---

## The clean-environment install — skipped, and why

Note 002 stopped before Stage 5 (install from `requirements.txt` in a fresh runtime, `USAGE.md`, ZIP + relocation test). It is **still not done**, and I did not do it here:

- A fresh Scanpy environment from the 123-pin `tmp/build/requirements.txt` costs **~900 MB**, and a fresh TISSUE one **~700 MB**. The brief was no new gigabyte-scale venv, so this is out of scope by instruction.
- What exists and is untested: `tmp/build/requirements.txt` (123 pins), `constraints.txt`, a built scanpy wheel under `tmp/build/wheels`, and `scanpy-src`. **None has ever been installed.** Note 002 flagged this and it remains true.
- The light check above is a weaker claim: the tools import and run **in the environment they were built in**. It says nothing about whether the pins reproduce that environment elsewhere.

Anyone picking this up should treat "does `requirements.txt` actually work?" as the first open question.

---

## State at the end of Phase 1

**Done:** Scanpy tutorial → 9 tested MCP tools (002), reused on two unseen datasets (003, 006), three input guards added (007, 008) and covered by 14 tests (009); TISSUE tutorial run (004) and its smallest useful slice wrapped as 2 validated functions (005).

**Suite:** 64/64 passing. **Heaviest job:** 6.4 GB peak. **Nothing needed a GPU.**

**Open, in the order I would take them:**

1. `requirements.txt` has never been installed — Stage 5 of note 002.
2. The reference `.h5ad` files and the original 10 test files are not in git; losing the working tree loses the ability to verify anything.
3. TPM/CPM into QC is still silently accepted (008).
4. TISSUE has no reference values of its own — `test.py` asserts nothing numeric — so note 005's numbers are the only baseline.
5. `USAGE.md` and the ZIP delivery were never written.
