# 001 — Paper2Agent framework (Phase 1, note 1 of 10)

**Date:** 2026-09-22
**Machine:** Dell laptop, Omarchy Linux (kernel 7.2.5), 15 GiB RAM (about 11 GiB free at idle), 30 GiB swap, 8 CPU cores, **no NVIDIA GPU** (`nvidia-smi` isn't installed), 214 GB free disk.
**Scope rule:** Everything is under `~/Work/paper2agent/`. I didn't change Hyprland, the theme, shell rc files or the global Claude config.

## What I cloned

- `https://github.com/jmiao24/Paper2Agent` → `~/Work/paper2agent/src`, using a shallow clone (depth 1)
- Commit `8c2d059165ef8cdcb70dbea76655b9c2b55b38e6` (2026-09-17), 13 MB
- The README Quick Start says to install the `skills/paper2agent` skill into the coding agent, then prompt it with a paper URL, a code repo and an output directory.

## What actually ran

| Step | Result |
|---|---|
| Clone | ✅ OK |
| Install the skill | ✅ Installed as a **project skill** at `~/Work/paper2agent/.claude/skills/paper2agent/` instead of `~/.claude/skills/`, to stay inside this folder. `claude -p` started in `~/Work/paper2agent` confirmed the skill is visible (it answered "YES"). |
| Hosted MCP: Scanpy | ✅ `https://paper2agent-scanpy-mcp.hf.space/mcp` handshake OK, **7 tools** (quality_control, normalize_data, select_features, reduce_dimensionality, build_neighborhood_graph, cluster_cells, annotate_cell_types) |
| Hosted MCP: TISSUE | ✅ `https://paper2agent-tissue-mcp.hf.space/mcp` handshake OK, **6 tools** (predict_spatial_gene_expression, calibrate_uncertainties…, multiple_imputation_hypothesis_testing, 2× cell filtering, tissue_weighted_pca) |
| Hosted MCP: AlphaGenome | ✅ `https://paper2agent-alphagenome-mcp.hf.space/mcp` handshake OK, **22 tools** |
| Real tool call: AlphaGenome `create_genomic_interval` (chr11:116837600-116837700) | ✅ `isError:false`. It returned a message and a CSV path, but the CSV stays on the **server** (`/data/tmp_outputs/...`). |
| Real tool call: Scanpy `quality_control` with no data | ⚠️ Returned the error I expected: "Path to h5ad file or 10X data directory must be provided". The server is live and checks its inputs. |
| `uv` (the skill requires it) | ✅ Installed uv 0.12.17 locally to `tools/bin`. Python and caches also go under `tools/` via `env.sh`. `tools/` is 232 MB. |
| Python 3.12 venv + `fastmcp==4.0.3` (the skill's MCP runtime) | ✅ Installs and imports (in `probe/.venv-smoke`) |
| The skill's own self-test, `paper2skill/scripts/test_paper_bundle.py` (27 tests) | ⚠️ **Incomplete, see below** |

### Skill self-test progression
All of these ran with `uv run --python 3.12 --with pytest ...`:
1. `--with pytest`: collection error, `No module named 'pymupdf'`
2. `+ pymupdf`: **18 failed / 9 passed**, `No module named 'PIL'`
3. `+ pymupdf + pillow`: **6 failed / 21 passed**. All 6 fail with `No module named 'pymupdf4llm'`.
4. `+ pymupdf4llm`: **not run.** You stopped it.

The 6 remaining failures come from one missing dependency, not from bugs. The last run should pass. The paper→skill half needs `pymupdf`, `pillow` and `pymupdf4llm`, and the repo doesn't list them in a requirements file.

## What failed and why

- **Hardware:** Nothing hit the 16 GB or GPU limit, because I ran no conversion. The limits still apply to future runs:
  - Any repo whose selected tutorial *requires* CUDA is blocked on this machine. The skill classifies each repo's GPU need as `required`, `optional`, `none` or `unknown`. `required` means stop.
  - Anything that needs more than about 11 GiB of working RAM will spill into swap or crash. Stop and record it.
- **Auth:** AlphaGenome's prediction tools (`predict_dna_sequence`, `score_variants_batch`, …) need an `api_key` argument, a Google DeepMind AlphaGenome key. I don't have one, so I didn't try them.
- **Missing data / design limit:** The hosted Scanpy and TISSUE tools take **file paths on the Hugging Face server**. You can't send a local `.h5ad` or `.txt` to them. So the hosted demos show tool shape only. Running them on your own data means running the MCP locally.
- **Tooling gaps (fixed locally):** No `uv`, and system Python 3.14 has no pip. uv-managed Python 3.12 works around both.
- **Not done:** No full `/paper2agent` conversion yet. I didn't register any MCP server with `claude mcp add`.

## Exact commands to repeat

```bash
# 0. Setup (one time)
mkdir -p ~/Work/paper2agent && cd ~/Work/paper2agent
git clone --depth 1 https://github.com/jmiao24/Paper2Agent src
mkdir -p .claude/skills/paper2agent && cp -R src/skills/paper2agent/. .claude/skills/paper2agent/
mkdir -p tools && curl -LsSf https://astral.sh/uv/install.sh | \
  env UV_INSTALL_DIR="$HOME/Work/paper2agent/tools/bin" UV_NO_MODIFY_PATH=1 sh
# env.sh already exists: sets PATH, UV_PYTHON_INSTALL_DIR, UV_CACHE_DIR, UV_TOOL_DIR under ~/Work/paper2agent

# 1. Every session
cd ~/Work/paper2agent && . ./env.sh

# 2. Probe a hosted MCP server (handshake + list tools)
U=https://paper2agent-scanpy-mcp.hf.space/mcp
H=(-H 'Content-Type: application/json' -H 'Accept: application/json, text/event-stream')
curl -s -D hdr -X POST "$U" "${H[@]}" -d '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"probe","version":"0"}}}'
SID=$(grep -i '^mcp-session-id' hdr | awk '{print $2}' | tr -d '\r')
curl -s -X POST "$U" "${H[@]}" -H "mcp-session-id: $SID" -d '{"jsonrpc":"2.0","method":"notifications/initialized"}'
curl -s -X POST "$U" "${H[@]}" -H "mcp-session-id: $SID" -d '{"jsonrpc":"2.0","id":2,"method":"tools/list"}'

# 3. Skill self-test (the final step, not run yet)
cd .claude/skills/paper2agent/paper2skill/scripts
uv run --python 3.12 --with pytest --with pymupdf --with pillow --with pymupdf4llm \
  python -m pytest -q test_paper_bundle.py

# 4. Optional: use a hosted server inside Claude Code, scoped to this folder only
cd ~/Work/paper2agent && claude mcp add --scope project --transport http scanpy https://paper2agent-scanpy-mcp.hf.space/mcp
```

Probe outputs from this session are in `~/Work/paper2agent/probe/`: `*.init`, `*.tools`, `ag_call.out` and `sc_call.out`.

## Next paper to try

**Scanpy, "Preprocessing and clustering" tutorial only** (Wolf et al., *Genome Biology* 2018; https://github.com/scverse/scanpy).
- Small, CPU-only, pip-installable, and the example dataset is a few MB of PBMC data.
- The README has the exact prompt, and the hosted Scanpy MCP gives 7 reference tools to compare against.
- Run from `~/Work/paper2agent` after `. ./env.sh`:
  ```
  /paper2agent Convert https://github.com/scverse/scanpy into tested MCP tools in Scanpy_Agent.
  Focus on the "Preprocessing and clustering" tutorial. Resource limits: CPU only, no GPU,
  stay under 10 GB RAM, do not download datasets larger than 500 MB.
  ```
- Fallback: TISSUE (https://github.com/sunericd/TISSUE). It's also CPU-based and has a hosted reference, but I haven't checked its data size.
- Avoid for now: AlphaGenome, which needs an API key and targets a huge genomic model.
