# use-004 — tests for the INBOX door

**Date:** 2026-09-25 · **Machine:** Dell laptop, no GPU, CPU only · **Not a Phase 1 note.**
**Result: 22 tests, all passing, 2.4 s, peak 275 MB.** Three rules mutation-checked. No new environment, no new method, **nothing downloaded, no INBOX contents invented**.

The door shipped in [use-003](use-003-inbox.md) with its seven cases run by hand and nothing in
the suite. This closes that gap: `Scanpy_Agent/tests/code/inbox/test_inbox_door.py`.

## What passed

**22 tests in 2.4 s.** The whole file is fast because almost none of it touches real data — the
door's job is deciding, and most decisions can be checked without a dataset.

### Empty tray (4 tests)

| Test | What it pins |
|---|---|
| `test_empty_tray_refuses` | exit **2**, message says "empty" |
| `test_empty_tray_downloads_nothing` | `urllib.request.urlopen` **and** the door's own `inspect` are replaced with functions that fail the test if called; also asserts the tray is still empty afterwards |
| `test_readme_alone_still_counts_as_empty` | the shipped `README.md` is not mistaken for a dropped file |
| `test_missing_tray_refuses` | a non-existent folder refuses rather than crashing |

The download test is the one that matters: it does not merely check the exit code, it makes any
network call an assertion failure.

### Wrong number or kind of item (4 tests)

Two files → refused, naming both. A `.txt` → refused. A folder containing `matrix.mtx` → refused
**with the conversion hint**. A corrupt `.h5ad` → exit **3** with a message, not a traceback —
exit 3 ("couldn't read it") is deliberately distinct from exit 2 ("refused on purpose").

### Size cap (3 tests)

**No 80 MB blob was created.** The rule under test is "a file whose size exceeds the limit is
refused", so the test shrinks `inbox.MAX_FILE_MB` to ~100 bytes for one case and points the door
at a tiny file. The limit's actual value is asserted separately
(`test_size_limit_constant_is_what_the_docs_claim`: it is 80), so the two halves together cover
what a real oversize file would prove, at zero disk cost. The refusal test also asserts the door
**did not modify or delete** the file it turned away.

### Memory look-ahead (4 tests)

A pure function of shape, so no file is needed at all. A fabricated 300,000 × 30,000 shape is
refused; pbmc3k's shape is allowed. One parametrised test checks the estimate stays **within 2×
of every peak actually recorded** in notes 003, 006 and use-001:

| Shape | Observed | Estimated |
|---|---|---|
| 700 × 765 | 708 MB | 608 MB |
| 2,730 × 3,451 | 954 MB | 741 MB |
| 2,700 × 13,714 | 1,153 MB | 1,155 MB |
| 17,041 × 23,427 | ~6,100 MB | 6,585 MB |

A fourth test asserts the 8 GB budget stays below the 10 GB `cap10g` ceiling.

### Raw-counts lock (2 tests)

A file carrying the `uns['log1p']` stamp is refused with "not raw counts". Separately, when the
object has a counts layer, the refusal names `layers['counts']` as the fix.

### Inspect-only and accept (5 tests)

`--inspect-only` prints the verdict and the look-ahead, does **not** run QC, and creates no
output directory. The accept path copies the **smallest real teaching file already on this disk**
(pbmc68k_reduced counts-restored, 566 KB, 700 × 765) into a temp tray, runs the door end to end,
and asserts a QC artifact was written and the "nothing was clustered" line appeared. The temp
tray dies with the test; **the real INBOX is never touched**.

That test **skips** if the teaching file is absent, because no `.h5ad` is committed to the repo.
On a fresh clone it will skip, and the suite reports 21 passed, 1 skipped.

## Mutation check

Three rules were broken on purpose to confirm the tests notice. The source was restored from a
backup and `git diff` confirmed it byte-identical.

| Mutation | Caught by |
|---|---|
| Size cap never fires | `test_oversize_refuses_without_building_a_blob` |
| Memory look-ahead ignored | `test_lookahead_refuses_a_shape_too_big_for_this_box` |
| Empty tray allowed through | `test_readme_alone_still_counts_as_empty` (plus the empty-tray tests) |

**One mutation taught me something.** With the size cap disabled, the oversize test fell through
to actually running QC on the 10×5 mechanism fixture, which died inside matplotlib with
`Positions outside range of features`. The test still failed — correctly — but via exit 3 rather
than exit 2. The lesson: **the tiny 10×5 fixture is too degenerate to survive real QC**, so it is
only ever used on paths that refuse *before* QC runs. The genuine accept path is covered solely
by the real teaching file. That is a real limit of this test file, not a hypothetical one.

## What was not tested

- **The door on a real study's data.** Still the central untested thing, and the reason the door
  exists. Every input here is either a 10×5 mechanism fixture or a teaching file.
- **The accept path on a fresh clone** — it skips, because no data file is in git.
- **The 10x `.h5` branch.** The door accepts that format, but no 10x `.h5` file was used; only
  `.h5ad`. The reference `.h5` files used in note 002 are not in the repo.
- **Concurrency** — two people using the tray at once is undefined and untested.
- **Symlinks, permissions, unreadable directories** — not exercised.
- **The `--report` JSON output** is written by the door but never asserted on.
- **Real oversize behaviour.** The 80 MB rule is proven by shrinking the limit, not by an 80 MB
  file. If `stat()` ever behaved oddly on a genuinely huge file, these tests would not see it.
- The mechanism fixtures are **not biology** — a 10×5 matrix of zeros. Nothing here says anything
  about data quality in any real sense.

## RAM and time

| Run | Result | Time | Peak RAM |
|---|---|---|---|
| Door tests alone | **22 passed** | 2.4 s | 275 MB |
| Run three times over, to check for order dependence | 22 passed each time | ~2.4 s | — |
| Door tests + the note 009 guard tests | **36 passed** | 16 s | 523 MB |

Everything under `cap10g`, 0 OOM kills. The door tests are the cheapest thing in this project by
a wide margin — the full clustering suite takes ~7 minutes and peaks at 6.4 GB.

## How to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env && cd "$PROJECT_ROOT"
cap10g "$PROJECT_PYTHON" -m pytest tests/code/inbox -q          # 22 passed, ~2.4 s
```
