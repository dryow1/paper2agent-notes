# paper2agent notes

Field notes from running [Paper2Agent](https://github.com/jmiao24/Paper2Agent) on a 16 GB Dell
laptop with **no GPU** — everything CPU-only, capped at 10 GB RAM.

- [001 — Paper2Agent framework](NOTES/001-paper2agent-framework.md): what the framework is, what
  installed, what the hosted MCP servers actually return.
- [002 — Scanpy clustering → MCP tools](NOTES/002-scanpy-clustering.md): converting Scanpy's
  "Preprocessing and clustering" tutorial into 9 tested MCP tools (stages 1–4 of 6).
- [003 — Reusing those tools on pbmc3k](NOTES/003-scanpy-pbmc3k.md): the same 9 tools run
  unchanged on a different dataset in 37 s at 1.2 GB peak — what reused cleanly and what didn't.
- [004 — TISSUE tutorial](NOTES/004-tissue-tutorial.md): the spatial-uncertainty paper's own
  420 KB example, 13 s at 420 MB — calibration reproduces, but only on a pinned 2022 stack.

Notes only: no code, data or environments. Results, limits and failures are recorded as they
happened.
