# USAGE — running Scanpy_Agent from this working tree

Nine MCP tools extracted from the scanpy "Preprocessing and clustering" tutorial
([commit `0d5fd16`](https://github.com/scverse/scanpy/blob/0d5fd16234865619d2f5097d33fc4281900a2bc2/docs/tutorials/basics/clustering.ipynb)),
built and verified in note 002, hardened by three input guards in notes 007–008.

**Scope of this document.** It describes running the tools **from this working tree, in the
environment they were built in** (`Scanpy_Agent/scanpy-env`). That is the only configuration
anyone has ever run them in.

> **What has never been tested:** installing from `tmp/build/requirements.txt` into a fresh
> environment (Stage 5 of note 002, still open). The 123 pins, the `constraints.txt` and the
> built wheel under `tmp/build/` exist but **have never been installed anywhere**. Nothing is
> registered with `claude mcp add` either. So "the tools work" means *in this tree, on this
> machine* — it is not a claim that the pins reproduce this environment elsewhere. Treat
> "does `requirements.txt` actually work?" as an open question, not a documented step.

---

## 1. Run it

Everything below is run from the repo root, `~/Work/paper2agent`. The interpreter is the
project's own venv; there is nothing to activate and nothing to install.

```bash
cd ~/Work/paper2agent
Scanpy_Agent/scanpy-env/bin/python --version      # Python 3.12.14
```

### As an MCP server (stdio)

```bash
cd ~/Work/paper2agent/Scanpy_Agent
./scanpy-env/bin/python src/scanpy_mcp.py
```

Server name `scanpy`, transport stdio, fastmcp 4.0.3. It mounts the `clustering` sub-server
without a prefix, so the nine tools appear under their own names. Confirm the mount without
starting the server:

```bash
cd ~/Work/paper2agent/Scanpy_Agent && ./scanpy-env/bin/python -c "
import sys, asyncio; sys.path.insert(0, 'src'); sys.argv = ['x']
from scanpy_mcp import mcp
print(mcp.name, len(asyncio.run(mcp.list_tools())))"
# -> scanpy 9
```

Registering it with a client is **not** done and not covered here. It would be a
`claude mcp add` pointing at that interpreter and that script, but that command has never
been run against this server, so it is left for whoever needs it.

### As plain Python

The tools are ordinary functions; the `@clustering_mcp.tool()` decorator does not get in the
way. This is how notes 003, 006 and use-001 drove them.

```bash
cd ~/Work/paper2agent
Scanpy_Agent/scanpy-env/bin/python -c "
import sys; sys.path.insert(0, 'Scanpy_Agent/src')
from tools.clustering import scanpy_compute_qc_and_filter
r = scanpy_compute_qc_and_filter(
    data_path='Scanpy_Agent/notebooks/clustering/data/00_raw_concat.h5ad',
    output_dir='Scanpy_Agent/tmp/demo')
print(r['n_obs_before'], '->', r['n_obs'], 'cells')"
```

### The nine tools, in pipeline order

Each takes `data_path` (an `.h5ad`; the first also accepts a 10x `.h5`) and an optional
`output_dir`. Each returns a dict with an `artifacts` list of `{description, path}` plus
tool-specific counts, and writes a **fresh timestamped subdirectory** per call.

| # | Tool | Key arguments beyond `data_path` |
|---|---|---|
| 1 | `scanpy_compute_qc_and_filter` | `mt_prefix='MT-'` (`'mt-'` for mouse), `ribo_prefixes=['RPS','RPL']`, `hb_pattern='^HB[^(P)]'`, `min_genes=100`, `min_cells=3` |
| 2 | `scanpy_detect_doublets_scrublet` | `batch_key=None` (tutorial: `'sample'`) |
| 3 | `scanpy_normalize_log1p` | `target_sum=None` (median total count) |
| 4 | `scanpy_select_highly_variable_genes` | `n_top_genes=2000`, `batch_key=None` |
| 5 | `scanpy_run_pca` | `n_comps=None` (scanpy default 50), `color_keys=None` |
| 6 | `scanpy_build_neighbors_umap` | `n_neighbors=15`, `n_pcs=None`, `color_keys=None` |
| 7 | `scanpy_cluster_leiden` | `resolutions=[1.0]`, `n_iterations=2`, `key_added='leiden'`, `qc_color_keys=None` |
| 8 | `scanpy_plot_marker_gene_dotplot` | **`marker_genes`** and **`groupby`** required; `cluster_labels`, `annotation_key='cell_type'`, `ignore_missing_genes=False` |
| 9 | `scanpy_rank_marker_genes` | **`groupby`** required; `n_genes=5`, `focus_group=None` |

With several resolutions, tool 7 writes `obs` keys as `{key_added}_res_{res:4.2f}` — so
`resolutions=[0.02, 0.5]` gives `leiden_res_0.02` and `leiden_res_0.50`. Those are the names
tools 8 and 9 want for `groupby`.

**Where outputs land.** If you omit `output_dir`, output goes to
`$TMPDIR/scanpy_mcp_outputs/<tool>_<timestamp>_<hex>/` — **outside this repo**. Always pass
`output_dir` for work you intend to keep or find again. Every call makes a new directory
rather than overwriting, so repeated runs accumulate; see §3.

### Three guards that will stop you

`scanpy_compute_qc_and_filter` refuses input that is demonstrably not raw counts, because
`sc.pp.calculate_qc_metrics` just sums `X` and returns a plausible-looking number for
anything. Each rule keys on a positive signature of a transformation, so ambient-corrected
and otherwise fractional *counts* still pass:

1. `uns['log1p']` present — scanpy's own stamp. Definitive.
2. Negative values — scaled or z-scored.
3. Fractional values with a maximum under 50 — the signature of `log1p`.
4. Fractional values with near-constant per-cell totals — `normalize_total()` was run.

The error names the problem and, when the object has `layers['counts']` or a `.raw` slot,
points at it. The usual fix is `adata.X = adata.layers['counts']`.

**TPM/CPM input is still accepted silently.** Known, deliberate, unfixed — it is
indistinguishable from ambient-corrected counts, which must pass. Do not read a clean QC run
as proof the input was counts.

Tool 8 rejects marker genes absent from `var_names` by default; pass
`ignore_missing_genes=True` for a panel from another tissue (pbmc3k drops 8 genes,
pbmc68k_reduced drops 34 and empties 2 sets). Tools 6–9 require what the previous step
produced and say which step is missing.

### The INBOX door

For a one-file QC check with no clustering, don't drive the tools directly — use the door,
which inspects first and refuses in one line rather than guessing:

```bash
cd ~/Work/paper2agent
Scanpy_Agent/scanpy-env/bin/python Scanpy_Agent/src/inbox.py --inspect-only
Scanpy_Agent/scanpy-env/bin/python Scanpy_Agent/src/inbox.py
```

Every run writes a JSON receipt (default `RESULTS/inbox-last.json`). Full instructions are in
`INBOX/HOW-TO.md`; the standing rules it obeys are in `SHOP.md`.

---

## 2. Check the baseline

`tests/MANIFEST.md` records the SHA-256, byte size, shape and source of the 32 files the
factory cannot rebuild. It is the thing that makes their loss or silent replacement
*detectable*; it does not back them up.

```bash
cd ~/Work/paper2agent/Scanpy_Agent
./scanpy-env/bin/python -m pytest tests/code/test_manifest_hashes.py -q
# -> 36 passed in ~1s   (4 structural checks + 32 files, nothing skipped)
```

Read the result like this:

| Outcome | Meaning |
|---|---|
| **36 passed** | every recorded file is present and byte-identical to the baseline. |
| **`s` / skipped** | a recorded file is **absent**. Expected for anything listed in MANIFEST.md's "deliberately leaves out" section after a cleanup; alarming for a `ref/*.h5ad`. |
| **failed** | a file is present but its hash moved: corruption, or a different file under the same name. **Every number in `NOTES/` was measured against the old one and is now unverifiable.** Stop and find out which, rather than re-recording the hash. |

The test holds no hashes of its own — it parses MANIFEST.md, so the manifest stays the single
source of truth. Four structural tests guard the parser itself, including one that fails if
the table yields fewer than 32 rows; without it a broken regex would make every file check
pass vacuously.

### The rest of the suite

```bash
cd ~/Work/paper2agent/Scanpy_Agent
./scanpy-env/bin/python -m pytest --collect-only -q | tail -1   # 132 tests collected
./scanpy-env/bin/python -m pytest -q                            # the whole suite
```

132 tests: **50** original tool tests (note 002) + **14** guard tests (009) + **32** INBOX door
tests + **36** manifest tests. The 64 clustering tests peak at **6.4 GB RAM** and compare
against `notebooks/clustering/ref/*.h5ad`, so they need the reference data present and room to
run. Run heavy jobs in the **foreground**: `cap10g` uses `systemd-run --user --scope`, which
dies with the launching shell and leaves a truncated log.

The manifest and door tests are cheap and need no reference data; the manifest test is the one
worth running after any cleanup.

---

## 3. What not to delete

`Scanpy_Agent/tests/MANIFEST.md` lists every path below with a hash, and
`tests/code/test_manifest_hashes.py` tells you afterwards whether a cleanup took something it
should not have. **Run it after deleting anything.**

### Do not delete

| Path | Size | Why |
|---|---|---|
| `Scanpy_Agent/tmp/cache/scverse_tutorials/s1d{1,3}_filtered_feature_bc_matrix.h5` | 43 MB | **The root of everything.** The two figshare lanes every reference file is derived from. Re-fetching them is a download, which the standing order forbids without a ticket. |
| `Scanpy_Agent/notebooks/clustering/ref/*.h5ad` + `data/00_raw_concat.h5ad` | ~1.1 GB | the reference checkpoints all 64 clustering tests compare against. |
| `Scanpy_Agent/tests/code/clustering/` | 88 KB of source | the 50 tests, `clustering_verify_helpers.py` and `conftest.py` — 13 files. **Only `test_tool_guards.py` is in git**; the other twelve exist nowhere else. |
| `Scanpy_Agent/scanpy-env` | 868 MB | rebuilding costs ~900 MB of downloads, and the pins have never been proven to reproduce it. |
| `TISSUE_Agent/tissue-env-pinned` | 693 MB | same, and the pinned 3.10 stack is fiddly: TISSUE needs `setuptools<81` and breaks on anndata ≥ 0.11. |
| `Scanpy_Agent/tmp/build/` | small | the 123 pins and the wheel — the only record of the intended environment. |
| `Scanpy_Agent/reports/` | 264 KB | the scanner, executor, verification and acceptance records from note 002. |
| `data/`, `INBOX/*.h5ad`, `TISSUE_Agent/repo/TISSUE/tests/data/` | ~50 MB | inputs already here that the standing order forbids fetching again. |

### Safe to delete (the suite or a rerun recreates them)

`Scanpy_Agent/tmp/outputs` (3.9 GB) · `notebooks/clustering/replay` (1.1 GB, already verified
bitwise identical to `ref/`) · `tests/data/sub4_*.h5ad` (343 MB, rebuilt on demand by
`subsample()`) · `tmp/pbmc3k`, `tmp/pbmc68k`, `tmp/guards`, `tmp/adamson` · `TISSUE_Agent/out`
· `INBOX_OUT` · `tools/uv-cache` (2.2 GB, but any reinstall then re-downloads).

Deleting all of those frees roughly 7.5 GB and leaves ~3.5 GB: two environments, two
checkouts, the reference data, the source and the notes. Afterwards, the manifest test should
report **36 passed** — those paths are deliberately not recorded in it. Any *skip* means
something from the first table went too.

### The structural weakness, stated plainly

**This repo publishes the record and the source, not a runnable checkout.** `.gitignore`
admits the notes, the extracted tools, four test files and the manifest; everything else here
is gigabytes of `.h5ad` and two venvs. `test_tool_guards.py` is committed to be *read*
alongside the code it guards — it cannot run from a fresh clone, because it needs
`conftest.py`, `clustering_verify_helpers.py` and the reference data, none of which are in
git. Losing this working tree loses the ability to verify anything in `NOTES/`; the manifest
only guarantees you would find out.
