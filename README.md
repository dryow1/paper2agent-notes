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

Days the tools were actually used on something: [use-001 — Scanpy on paul15](RESULTS/use-001-scanpy.md)
and [use-002 — TISSUE calibration](RESULTS/use-002-tissue.md).

---

## What this factory can and cannot do

**Two cupboards, eleven tools.** The **Scanpy** cupboard holds **9 tools** that take a table of
cells and genes and walk it through the standard single-cell pipeline: check quality, spot
doublets, normalise, pick the informative genes, reduce dimensions, draw a map, group the cells,
and name the genes that mark each group. The **TISSUE** cupboard holds **2 tools** that take
spatial measurements plus a detailed reference, guess a gene that was never measured in place,
and — the real point — say how much to trust each guess. The two cupboards do not share a
workshop: they need different, incompatible versions of Python and its libraries.

**Everything here was built on tiny public teaching files.** The biggest thing ever downloaded
was 43 MB; most were under 10 MB. These are the example datasets that ship with the software so
people can learn it. **Nothing here has been tested on a real study's data**, and no claim in
these notes is evidence about biology — they are evidence that the machinery behaves as
described on small, well-understood inputs.

**There are locks on the doors, and they are the most useful thing built.** Early on, the quality
tool accepted data that had already been transformed and cheerfully returned nonsense — negative
"total counts" for more than half the cells, with no complaint. Twice. That is the worst kind of
failure, because nothing looks wrong. So the tools now **refuse** data that has been logged,
scaled, or rescaled to a fixed total when raw counts are what's required, and the refusal says
which fix to apply. They deliberately still accept fractional data that really is counts, because
rejecting that would break honest work. **No silent quality check**: if a marker gene is missing,
you are told which ones and how many, not quietly given a smaller plot. Fourteen tests exist
purely to prove the locks still hold, and were checked by breaking each lock on purpose to
confirm the tests noticed.

**Two use days so far.** The Scanpy tools ran on **paul15**, mouse bone marrow — 2,730 cells, 30
seconds — and recovered textbook groups nobody told them about: haemoglobin genes in the red-cell
cells, platelet genes in the platelet precursors, granule enzymes in the immune cells. One odd
jar was **kept rather than tidied away**: a 22-cell group labelled "neutrophil" whose genes
actually look like a different cell type entirely. It is written down as a disagreement, not
smoothed over or resolved. The TISSUE tools ran on their **411 KB** built-in file: asked for 77 %,
90 % and 50 % confidence, they were right 78.8 %, 90.8 % and 54.0 % of the time — close, and
erring toward caution. Those numbers have now come out identical on three separate runs through
two different code paths.

**This laptop is the ceiling.** 16 GB of memory, **no graphics card**, a 256 GB disk. The heaviest
job run here peaked at 6.4 GB. That is fine for teaching files and small studies; it is **not** a
machine for large papers. Anything needing a GPU is out entirely — that was a hard rule
throughout, and nothing CUDA-related was ever installed. Downloads were capped at 500 MB, so the
large public spatial datasets were never fetched. A real study's dataset would likely exceed all
three limits at once.

**What happens next needs a human to name a file.** The tools work and the locks hold, but every
run so far has been on data chosen because it was small and convenient. There is no queue of real
work. **The next real job needs someone to point at a specific dataset and say what question it
should answer** — until then, the right move is to stop rather than invent another demonstration.
