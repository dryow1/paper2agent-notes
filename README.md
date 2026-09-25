# paper2agent notes

Field notes from running [Paper2Agent](https://github.com/jmiao24/Paper2Agent) on a 16 GB Dell
laptop with **no GPU** — everything CPU-only, capped at 10 GB RAM. Phase 1, ten notes.

**Start with [010 — the operator sheet](NOTES/010-operator-sheet.md)**: how to run everything from
a cold start, what is authoritative, what fails on this machine, what is safe to delete.

- [001 — Paper2Agent framework](NOTES/001-paper2agent-framework.md): what the framework is, what
  installed, what the hosted MCP servers actually return.
- [002 — Scanpy clustering → MCP tools](NOTES/002-scanpy-clustering.md): converting Scanpy's
  "Preprocessing and clustering" tutorial into 9 tested MCP tools (stages 1–4 of 6).
- [003 — Reusing those tools on pbmc3k](NOTES/003-scanpy-pbmc3k.md): the same 9 tools run
  unchanged on a different dataset in 37 s at 1.2 GB peak — what reused cleanly and what didn't.
- [004 — TISSUE tutorial](NOTES/004-tissue-tutorial.md): the spatial-uncertainty paper's own
  420 KB example, 13 s at 420 MB — calibration reproduces, but only on a pinned 2022 stack.
- [005 — Wrapping TISSUE](NOTES/005-tissue-tools.md): its smallest useful slice as 2 validated
  functions; reproduces the tutorial's numbers, 7/7 bad inputs give a clear error.
- [006 — The tools on pbmc68k_reduced](NOTES/006-scanpy-pbmc68k-reduced.md): a *processed*
  dataset this time, which found a silent wrong answer — QC accepting scaled data.
- [007 — Two guards](NOTES/007-scanpy-tool-guards.md): QC refuses scaled input;
  `ignore_missing_genes` saves hand-editing marker panels. 44 lines, suite still green.
- [008 — Closing the QC hole](NOTES/008-scanpy-qc-raw-counts.md): four rules that catch
  log-normalised and target-summed data while still accepting fractional real counts.
- [009 — Guard tests](NOTES/009-scanpy-guard-tests.md): 14 pytest cases, mutation-checked so
  they can actually fail. Suite 50 → 64.
- [010 — Operator sheet](NOTES/010-operator-sheet.md): the handover.

Notes, plus the tool source they describe. No data, environments or reference `.h5ad` files —
which means the committed tests are here to be read, not run. Results, limits and failures are
recorded as they happened, including the ones I got wrong first.
