# use-003 — the INBOX door

**Date:** 2026-09-25 · **Machine:** Dell laptop, no GPU, CPU only · **Not a Phase 1 note.**
**Built:** one way in for a file a human names. **No new environment, no new method, nothing downloaded.**
**Required dry run (INBOX empty): refused in 12 MB of memory and fetched nothing.**

## What the door is

`Scanpy_Agent/src/inbox.py` — about 250 lines that do one job: take a file a person put in
`~/Work/paper2agent/INBOX`, look at it, and either run the existing Scanpy quality check on it or
say no in one line.

It adds **no new biology**. The only computation it performs is `scanpy_compute_qc_and_filter`,
the same tool from note 002 with the locks from notes 007–008. Everything else is looking and
deciding.

### The rules, in order

| # | Check | Refusal reads like |
|---|---|---|
| 1 | INBOX exists and holds **exactly one** thing | `INBOX is empty … I will not download a dataset to have something to do.` |
| 2 | It is a file, not a folder | `'filtered_matrix' is a folder … A 10x mtx folder needs converting to .h5ad first; the tools do not read that layout.` |
| 3 | Suffix is `.h5ad` or 10x `.h5` | `'notes.txt' has suffix '.txt'; I accept .h5ad or .h5.` |
| 4 | **≤ 80 MB** | `'big.h5ad' is 85 MB, over the 80 MB door limit — this box is a 16 GB laptop with no GPU; subset it or use a bigger machine.` |
| 5 | **Inspect** — shape, values, layers, obs/var, precomputed results | *(never refuses; always printed)* |
| 6 | **Memory look-ahead ≤ 8 GB** | `300,000 cells x 30,000 genes would peak near 135,600 MB, over the 8,000 MB budget…` |
| 7 | `X` is **raw counts** | `X is not raw counts — uns['log1p'] stamp present … supply the count matrix.` |
| 8 | Run QC, then **stop** | `door: QC only. Nothing was clustered, annotated or named.` |

**Inspection always happens before any decision that costs memory**, and its output is printed
whether the file is accepted or not. `--inspect-only` stops after step 5, so you can look at a
file without committing to anything.

### Why a 10x folder is refused

The existing tools read `.h5ad` and a 10x Genomics `.h5` file — that is what
`tools/clustering.py::_read` accepts, and nothing more. The `matrix.mtx` + `barcodes` +
`features` folder layout is **not** supported. The door could convert it silently, but the
conversion involves choices (which features file, whether to make names unique, gzipped or not)
that belong to the person, not to a door. So it refuses and names the reason.

### The memory look-ahead

The file size cap alone is not enough — a sparse file can be tiny on disk and impossible in
memory. The estimate is fitted to the four runs already recorded in these notes:

| Dataset | Shape | Observed peak | Estimate |
|---|---|---|---|
| pbmc68k_reduced | 700 × 765 | 708 MB | 608 MB |
| paul15 | 2,730 × 3,451 | 954 MB | 741 MB |
| pbmc3k | 2,700 × 13,714 | 1,153 MB | 1,155 MB |
| bone marrow | 17,041 × 23,427 | ~6,100 MB | 6,585 MB |

`peak ≈ 600 MB + 15 MB × (cells × genes ÷ 1,000,000)`. The 600 MB floor is Python and scanpy
loading; the slope is the data. The budget is **8 GB**, not the 10 GB `cap10g` ceiling, so the
laptop stays usable while a job runs.

**This is a rough guide, not a guarantee** — it is a straight line fitted to four points, all of
them small. It will be wrong for unusual data, and it errs toward over-estimating, which is the
safe direction for refusing.

## The dry run — INBOX empty

The required test. INBOX contained only its `README.md`:

```
$ cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py
door: reading /home/kee/Work/paper2agent/INBOX
door: REFUSED - INBOX is empty (/home/kee/Work/paper2agent/INBOX) - drop one .h5ad or 10x .h5
       file in it. I will not download a dataset to have something to do.
cap10g: exit=2 peak_mb=12 oom_kills=0
```

**It refused. It did not fetch anything.** Peak memory 12 MB — it never even imported scanpy,
because there was nothing to read. Exit code 2 means "refused", distinct from 3 ("could not read
the file") and 0 ("ran").

## The other doors, checked

All using files already on this disk. **Nothing was downloaded, and INBOX was left empty
afterwards.**

| Case | Outcome |
|---|---|
| Raw-count `.h5ad` (pbmc3k, 7.4 MB) | ✅ **accepted** — inspected as `raw counts (non-negative and integral)`, look-ahead 1,926 MB, QC ran: `kept 2700/2700 cells and 13714/32738 genes`, three artifacts written |
| Log-normalised `.h5ad` | ⛔ refused at rule 7 — `uns['log1p'] stamp present` |
| `notes.txt` | ⛔ refused at rule 3 |
| Two `.h5ad` files | ⛔ refused at rule 1 |
| Folder containing `matrix.mtx` | ⛔ refused at rule 2, with the conversion hint |
| 85 MB file | ⛔ refused at rule 4 |
| 300,000 × 30,000, only **8 MB on disk** | ⛔ **passed the size cap, refused by the look-ahead** at 135,600 MB estimated |

That last row is why both checks exist. A file can be small and still impossible.

## What the door deliberately does not do

- **It does not download.** An empty INBOX is a refusal. There is no fallback demo dataset.
- **It does not cluster, annotate, or name cell types.** It runs quality control and stops. No
  groups are invented.
- **It does not convert formats**, guess a mitochondrial prefix (`--mt-prefix` defaults to human
  `MT-`; pass `mt-` for mouse), or repair a file that is the wrong kind.
- **It does not swallow errors.** A refusal is a one-line reason and exit 2. An unreadable file is
  exit 3 with the exception text, not a traceback dump.

## Honest limits

- **The accept path has only been exercised on files already here** — pbmc3k and the fixtures from
  note 008. It has never seen a real study's data, which is the entire point of building it and
  also the reason its behaviour on one is unproven.
- **The memory model is four points and a straight line.** Treat it as a guard rail, not a
  prediction.
- **80 MB is a judgement call**, not a derived number. It is comfortably above every file used in
  these notes (largest: 43 MB) and comfortably below anything that would strain the box.
- **No tests were written for the door itself.** The seven cases above were run by hand and are
  recorded here; they are not in the pytest suite, so nothing would catch a regression.
- The door imports the tools in-process. It does not go through the MCP stdio path.

## How to use it

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env

# 1. Put exactly one .h5ad (or 10x .h5) in INBOX/
# 2. Look before leaping
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py --inspect-only
# 3. Run QC (add --mt-prefix mt- for mouse)
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py --report INBOX_OUT/report.json
```

Artifacts land in `INBOX_OUT/`. The door is now waiting; **it has nothing to do until someone
names a file.**
