# use-005 — the door leaves a receipt

**Date:** 2026-09-25 · **Machine:** Dell laptop, no GPU, CPU only · **Not a Phase 1 note.**
**Result: every run now writes one flat JSON receipt, accept or refuse. 32 door tests pass in 2.4 s.**
No downloads, no new method, **nothing was put in INBOX**.

## What changed

The door already accepted a `--report` flag, but it wrote a nested structure, only when asked.
Now:

- **A receipt is written on every run**, default `RESULTS/inbox-last.json`, overridable with
  `--report FILE`, suppressible with `--no-report`.
- **The keys are flat and always present.** A reader never has to check whether a field exists;
  unknown is `null`.
- **`exit_code` is in the receipt**, which meant restructuring `main()` — the old version returned
  from inside `try`, so the exit code was not known when the report was written.

### The shape of a receipt

```json
{
  "timestamp": "2026-09-25T15:40:16",
  "decision": "refuse",
  "reason": "INBOX is empty (…) - drop one .h5ad or 10x .h5 file in it. I will not download a dataset to have something to do.",
  "exit_code": 2,
  "inbox": "/home/kee/Work/paper2agent/INBOX",
  "file_name": null,
  "size_mb": null,
  "shape": null,
  "layers": null,
  "raw_counts": null,
  "raw_counts_reason": null,
  "estimated_ram_mb": null,
  "ram_budget_mb": 8000,
  "artifacts": null
}
```

`decision` is one of **accept**, **refuse**, **error**, or **inspected** (`--inspect-only`, which
is honestly neither an accept nor a refuse). `shape` is `{"cells": n, "genes": m}`.

**The nulls are the point.** A refusal on an empty tray genuinely does not know a shape, because
nothing was read. Recording `null` rather than omitting the key says "unknown", not "forgotten".

An accepted run fills everything in:

```
decision: accept        exit_code: 0           file_name: counts.h5ad
size_mb: 0.553          shape: {'cells': 700, 'genes': 765}
layers: []              raw_counts: raw counts
estimated_ram_mb: 608   ram_budget_mb: 8000    artifacts: 3 paths
reason: QC computed; kept 700/700 cells and 765/765 genes
```

## Two bugs found while testing

Both were found by the new tests, not by reading the code.

**1. `HOW-TO.md` was being treated as a dropped file.** The previous version skipped exactly one
filename, `README.md`. Adding `INBOX/HOW-TO.md` in use-004 silently broke the empty-tray
refusal — the door started saying `'HOW-TO.md' has suffix '.md'; I accept .h5ad or .h5` instead
of `INBOX is empty`. Still a refusal, so nothing dangerous happened, but the wrong reason. Now
**any `.md` file in the tray is treated as furniture**, and a test pins it:
`test_tray_instructions_are_not_mistaken_for_a_dropped_file`.

**2. An oversize refusal recorded no size.** The line was:

```python
report["size_mb"] = round(check_size(path), 4)
```

`check_size` raises when the file is too big, so the assignment never ran and the receipt for the
one refusal that is *entirely about size* had `size_mb: null`. Fixed by recording the size first
and checking it second. The size is also rounded to 4 decimal places now, because 2 recorded a
4 KB test file as `0.0`.

## Tests

**32 tests, 2.4 s, peak 276 MB** — 22 from use-004 plus **10 new ones for the receipt**:

| Test | What it pins |
|---|---|
| `test_empty_tray_still_writes_a_receipt` | empty tray → receipt exists, `decision: refuse`, `exit_code: 2`, and file fields are `null` |
| `test_receipt_always_has_every_key` | the key set is exactly the 14 documented keys |
| `test_oversize_refusal_writes_a_receipt_with_the_size` | `size_mb > 0` — the regression above |
| `test_not_counts_refusal_records_what_was_inspected` | refused at the last gate, so `shape`, `raw_counts` and `estimated_ram_mb` are all filled |
| `test_unreadable_file_receipt_says_error_not_refuse` | `decision: error`, `exit_code: 3` |
| `test_inspect_only_receipt_is_neither_accept_nor_refuse` | `decision: inspected`, no artifacts |
| `test_no_report_writes_nothing` | `--no-report` is honoured |
| `test_receipt_default_path_is_used_when_no_flag` | the default path is used and parent directories are created |
| `test_bad_receipt_path_does_not_break_the_refusal` | an unwritable receipt path does **not** turn a clean refusal into a crash |
| `test_tray_instructions_are_not_mistaken_for_a_dropped_file` | the `.md` regression above |

Tests never write into the real `RESULTS/` — the helper adds `--no-report` unless a test asks for
a receipt, and those that do point it at `tmp_path`.

**No real accept fixture was required.** Nine of the ten use only 10×5 mechanism fixtures or an
empty folder. The accept path's receipt was checked by hand (output above) using the smallest
teaching file already on disk, copied to a temp tray and deleted; the pytest accept test skips
when that file is absent.

## What was not tested

- **A receipt from a real study's file** — same gap as always: nothing real has been through the door.
- **Overwrite behaviour.** The default receipt is one file, replaced each run. There is no history,
  and no test that a second run replaces the first cleanly.
- **Concurrent runs** writing the same receipt path — undefined.
- **The receipt on the 10x `.h5` branch**, which is still entirely unexercised.
- The `artifacts` list is asserted only for its length in the manual check, not per-path.

## RAM and time

| Run | Result | Time | Peak RAM |
|---|---|---|---|
| Door tests (32) | **32 passed** | 2.4 s | 276 MB |
| Run twice, to check stability | 32 passed both times | ~2.4 s | — |
| Manual accept, end to end | `decision: accept`, 3 artifacts | ~3 s | 267 MB |

## How to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env

# Empty tray: still refuses, still leaves a receipt
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py
cat RESULTS/inbox-last.json

cd "$PROJECT_ROOT" && cap10g "$PROJECT_PYTHON" -m pytest tests/code/inbox -q   # 32 passed
```

`RESULTS/inbox-last.json` is a runtime artifact and is **not** committed — only `RESULTS/*.md` is.
