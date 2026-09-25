# 002 — Scanpy "Preprocessing and clustering" → MCP tools (Phase 1, note 2 of 10)

**Date:** 2026-09-22
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, Intel Iris Xe graphics, **no NVIDIA GPU**
**Project:** `~/Work/paper2agent/Scanpy_Agent/`
**Result:** **Partial: Stages 1–4 of 6 passed.** You stopped the run at Stage 5. Stage 5 (clean-environment install check) and Stage 6 (USAGE.md and ZIP) were not run.
I didn't change Hyprland, the theme or any system config. Everything, including the data and caches, stayed under `~/Work/paper2agent`.

## Short version

The Scanpy clustering tutorial ran end to end on this Dell. The paper2agent skill turned it into **9 MCP tools**. They match the tutorial's own results exactly and pass independent tests and a real MCP server check. It needed **about 5–6 GB of RAM** at peak, **43 MB of data**, **no GPU**, and roughly 3 hours with agents.

## Inputs and limits

| | |
|---|---|
| Source | https://github.com/scverse/scanpy, commit `0d5fd16234865619d2f5097d33fc4281900a2bc2` (2026-09-17). Shallow clone in `repo/scanpy`, 91 MB |
| Tutorial | `docs/tutorials/basics/clustering.ipynb`, title "Preprocessing and clustering" (exact title match). The legacy `clustering-2017.ipynb` is excluded because its title differs. |
| Data | Two 10x samples from figshare doi:10.6084/m9.figshare.22716739.v1: s1d1 22.8 MB and s1d3 20.4 MB, **43 MB total**, md5 verified. Cached in `Scanpy_Agent/tmp/cache`. |
| Size | 17,125 cells × 36,601 genes raw, 17,041 × 23,427 after QC |
| Limits applied | CPU only; **10 GB cap** enforced per command with `tools/bin/cap10g` (a temporary per-command memory cap, no swap); data under 500 MB |
| Environment | Python 3.12.14 via local `uv`, scanpy `1.13.0a3.dev0+g0d5fd1623` built from the pinned checkout, fastmcp 4.0.3, igraph/leidenalg, scikit-image. `uv pip check` is clean. |

## What ran (stages and agents)

| Stage | Who | Result | Peak RAM |
|---|---|---|---|
| 0 Setup | me | Project, clone, routing record: Python route, `gpu_requirement: none` | – |
| 1 Environment | env-manager agent (~7 min) | ✅ | ~1 GB smoke test |
| 1 Selection | scanner agent (~33 min) | ✅ 9 tools selected; `setup` gate passed | – |
| 2 Reference execution | executor agent (~60 min, 2 of 5 attempts) | ✅ Real notebook run with papermill, 15 figures, a reference file saved after each step, fresh-process replay **bitwise identical**; `execution` gate passed | **4,991 MB** |
| 3a Implementation | implementer agent (~15 min) | ✅ `src/tools/clustering.py`, 514 lines of thin scanpy calls with no invented math; `extraction` gate passed | 5,057 MB |
| 3b Independent verification | a separate fresh verifier agent (~61 min) | ✅ **50/50 pytest passed**, 1 repair; `verification` gate passed | **6,086 MB** |
| 4 MCP integration | me | ✅ `src/scanpy_mcp.py` over real stdio: **20/20 acceptance cases** (11 successful calls covering all 9 tools, plus 9 error cases), every returned file compared with the references | 4,700 MB |
| 5 Clean-runtime install | – | ⛔ **Skipped** (you stopped it) | – |
| 6 USAGE.md and ZIP delivery | – | ⛔ **Not started** | – |

## The 9 tools

Each tool reads an `.h5ad` file and writes a new `.h5ad` plus the tutorial's figures, so the tools chain together:

1. `scanpy_compute_qc_and_filter`: QC metrics (mt/ribo/hb), violin and scatter plots, filter_cells(min_genes=100), filter_genes(min_cells=3). Also accepts a 10x `.h5`.
2. `scanpy_detect_doublets_scrublet`: Scrublet with a user `batch_key`.
3. `scanpy_normalize_log1p`: keeps raw counts in `layers['counts']`, then normalize_total and log1p.
4. `scanpy_select_highly_variable_genes`: 2,000 HVGs, flavor `seurat`, optional batch_key.
5. `scanpy_run_pca`: PCA plus variance-ratio and PCA plots.
6. `scanpy_build_neighbors_umap`: neighbors graph and UMAP.
7. `scanpy_cluster_leiden`: Leiden (igraph) at one or several resolutions, plus UMAP panels.
8. `scanpy_plot_marker_gene_dotplot`: dotplot of a **user-supplied** marker-gene dict, with an optional cluster→label map.
9. `scanpy_rank_marker_genes`: Wilcoxon rank_genes_groups, dotplot, and a CSV of all groups.

## Numbers matched the reference

- Tutorial-scale results: 213 doublets (207 in s1d1, 6 in s1d3), 2,000 HVGs, Leiden **26** clusters at resolution 1 and **5 / 17 / 36** at resolutions 0.02 / 0.5 / 2.0. The rank-genes table has 398,259 rows.
- Every tool's output file equals the reference file from the real notebook run: identical cell and gene names, identical labels, and a max difference of **0.0** on matrices, PCA/UMAP coordinates and scores. The comparison tolerance was rtol 1e-6.
- The rank-genes CSV matches the reference exactly. Group "7" top genes match the notebook printout (CD79A 39.370392, CD37, HLA-DRA, CD74, MS4A1).
- The same comparison repeated on files produced through the real MCP server passed too: `reports/acceptance-scientific-comparison-project.json`.
- On changed inputs (different min_genes, batch_key, n_top_genes, n_comps, resolution, marker lists, groupby), each tool's result matched a direct scanpy call. These tests used an every-4th-cell subsample.

## Caveats (honest)

- **Plot pixels were only partly compared.** Tool figures are pixel-identical to the same scanpy plot calls made directly, but **not** pixel-compared with the notebook's own images. The tools skip `set_figure_params`, so colors and dpi differ for 10 or fewer clusters. The QC violin plot uses random jitter, so it's only checked to exist and not be blank.
- **Old committed outputs disagree on cluster numbering.** The notebook's saved outputs came from an older environment. Their group "7" is NK cells; at this commit NK is group "8". I treated this as renumbering and used fresh references.
- **One bug found and fixed by the verifier:** the marker dotplot's cluster→label map turned every label into NaN when cluster categories were integers.
- **Defaults differ from the tutorial in two places, by design:** `batch_key` defaults to None (pass `"sample"` to reproduce the tutorial), and Leiden's multi-resolution run needs `n_iterations=-1`.
- The multi-resolution acceptance case doesn't produce the single-resolution `leiden` column. That column is covered by the separate single-resolution cases and the pytest suite.
- **Not done, because Stage 5 was skipped:** a fresh-environment install from `requirements.txt`, `USAGE.md`, the ZIP and its relocation test.
  - A scanpy wheel (2 MB) and a compiled `requirements.txt` (123 pins) exist in `tmp/build/` but were **never installed or tested**.
  - Nothing is registered with `claude mcp add`.
- Disk: the reference files are large. `notebooks/` is about 2.3 GB; delete `notebooks/clustering/replay/` (about 1.1 GB) if you need space.

## Where things are

- Server: `Scanpy_Agent/src/scanpy_mcp.py`, tools: `src/tools/clustering.py`
- Tests: `tests/code/clustering/` (log: `tests/logs/clustering/pytest_final.log`)
- Evidence:
  - `reports/`: scanner, executor, implementation, verification and acceptance reports; `agent-runs.json`; workflow-gate outputs; `mcp-project-environment.json`
  - `notebooks/clustering/`: executed notebook, driver, references, figures
- Stage markers: `.pipeline/`, up to `mcp_integration_done`

## Commands to repeat

```bash
cd ~/Work/paper2agent && . ./env.sh && . Scanpy_Agent/project.env && cd "$PROJECT_ROOT"

# Rerun the verifier's test suite (~9 min, ~6 GB peak)
cap10g "$PROJECT_PYTHON" -m pytest tests/code/clustering -q

# Rerun the real-stdio MCP acceptance check (~5 min)
rm -rf tmp/outputs/acceptance
cap10g "$PROJECT_PYTHON" "$SKILL_ROOT/scripts/verify_mcp_server.py" --server src/scanpy_mcp.py \
  --expected reports/expected-mcp-tools.json --cases reports/mcp-acceptance-cases.json \
  --require-all-tools --cwd "$PROJECT_ROOT" --call-timeout 900 --report reports/mcp-project-environment.json
cap10g "$PROJECT_PYTHON" tests/summary/compare_acceptance_outputs.py tmp/outputs/acceptance

# Check the recorded workflow gates
"$PROJECT_PYTHON" "$SKILL_ROOT/scripts/verify_workflow.py" --project-root "$PROJECT_ROOT" --through verification

# Optional (NOT done): use the server in Claude Code, scoped to this folder only
cd ~/Work/paper2agent && claude mcp add --scope project scanpy -- \
  "$PROJECT_PYTHON" "$PROJECT_ROOT/src/scanpy_mcp.py"

# To finish Stages 5–6 later: resume /paper2agent in ~/Work/paper2agent and say
# "resume Scanpy_Agent from Stage 5". Markers and records are in .pipeline/ and reports/.
```

## Next paper

**TISSUE** (Sun et al., uncertainty-calibrated spatial transcriptomics; https://github.com/sunericd/TISSUE).
- Small Python package, CPU, and it's the README's own example. The hosted `paper2agent-tissue-mcp` (6 tools) is a reference to compare against.
- Check its tutorial data size first. Stop if it's over 500 MB.
- Prompt:
  ```
  /paper2agent Convert https://github.com/sunericd/TISSUE into tested MCP tools in TISSUE_Agent.
  Resource limits: CPU only, no GPU, stay under 10 GB RAM, do not download datasets larger than 500 MB.
  ```
- Plan for about 3 hours with agents. Say up front whether to stop after Stage 4 as this run did.
