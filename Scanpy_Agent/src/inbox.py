"""The door: run the existing Scanpy QC on one human-placed file, or refuse and say why.

Contract
--------
* One file in ``~/Work/paper2agent/INBOX``: a ``.h5ad``, or a 10x Genomics ``.h5``.
* Inspect first, always. Shape, layers, obs/var, and whether raw counts are present are
  reported before anything is computed.
* Then either run ``scanpy_compute_qc_and_filter`` (which applies the notes 007-008 locks)
  or refuse with **one line** saying why.
* An empty INBOX is a refusal, never a reason to fetch a demo dataset.
* Nothing is clustered, annotated or named. This door does QC and stops.

Every refusal is a normal exit with a reason, not a traceback.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
DEFAULT_INBOX = PROJECT_ROOT.parent / "INBOX"

#: The tools read .h5ad and a 10x Genomics .h5 (tools/clustering.py::_read). A 10x *folder*
#: (matrix.mtx + barcodes + features) is NOT supported, so the door refuses it rather than
#: quietly converting - conversion is a decision for the person, not the door.
ACCEPTED_SUFFIXES = (".h5ad", ".h5")

#: Hard size ceiling. Bigger files are a conversation, not a drop.
MAX_FILE_MB = 80

#: Memory look-ahead, fitted to the runs recorded in NOTES/003, 006 and use-001:
#:   pbmc3k    2700x13714 -> 1,153 MB observed
#:   pbmc68k    700x765   ->   708 MB observed
#:   paul15    2730x3451  ->   954 MB observed
#:   bone marrow 17041x23427 -> ~6,100 MB observed
#: peak_MB ~= FLOOR + SLOPE * (n_obs * n_vars / 1e6); the floor is imports and scanpy itself.
RAM_FLOOR_MB = 600
RAM_SLOPE_MB_PER_MILLION = 15
#: Budget, not the cap. cap10g allows 10 GB; stopping at 8 GB leaves the box usable.
RAM_BUDGET_MB = 8000


class Refusal(Exception):
    """A clean 'no' with one line of reason."""


def _fmt_mb(n: float) -> str:
    return f"{n:,.0f} MB"


def find_input(inbox: Path) -> Path:
    """Pick the single acceptable file in INBOX, or refuse."""
    if not inbox.is_dir():
        raise Refusal(f"no INBOX at {inbox} - create it and put one .h5ad or 10x .h5 file in it.")
    entries = [p for p in sorted(inbox.iterdir()) if not p.name.startswith(".") and p.name != "README.md"]
    if not entries:
        raise Refusal(
            f"INBOX is empty ({inbox}) - drop one .h5ad or 10x .h5 file in it. "
            "I will not download a dataset to have something to do."
        )
    files = [p for p in entries if p.is_file() and p.suffix in ACCEPTED_SUFFIXES]
    folders = [p for p in entries if p.is_dir()]
    if len(entries) > 1:
        raise Refusal(
            f"INBOX holds {len(entries)} items ({', '.join(p.name for p in entries[:5])}"
            f"{'...' if len(entries) > 5 else ''}) - leave exactly one."
        )
    if folders:
        d = folders[0]
        hint = (
            " A 10x mtx folder needs converting to .h5ad first; the tools do not read that layout."
            if any((d / n).exists() for n in ("matrix.mtx", "matrix.mtx.gz"))
            else ""
        )
        raise Refusal(f"'{d.name}' is a folder; I accept a single .h5ad or 10x .h5 file.{hint}")
    if not files:
        bad = entries[0]
        raise Refusal(f"'{bad.name}' has suffix '{bad.suffix or '(none)'}'; I accept {' or '.join(ACCEPTED_SUFFIXES)}.")
    return files[0]


def check_size(path: Path) -> float:
    mb = path.stat().st_size / (1024 * 1024)
    if mb > MAX_FILE_MB:
        raise Refusal(
            f"'{path.name}' is {_fmt_mb(mb)}, over the {MAX_FILE_MB} MB door limit - "
            "this box is a 16 GB laptop with no GPU; subset it or use a bigger machine."
        )
    return mb


def inspect(path: Path) -> dict:
    """Read the file and describe it. No computation on the data beyond summary statistics."""
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    if path.suffix == ".h5ad":
        adata = ad.read_h5ad(path)
    else:
        import scanpy as sc

        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()

    X = adata.X
    values = X.data if sp.issparse(X) else np.asarray(X).reshape(-1)
    has_values = values.size > 0
    xmin = float(values.min()) if has_values else 0.0
    xmax = float(values.max()) if has_values else 0.0
    integral = bool(has_values and not np.any(values % 1))
    totals = np.asarray(X.sum(axis=1)).reshape(-1)
    mean_total = float(totals.mean()) if totals.size else 0.0
    constant_totals = bool(mean_total > 0 and float(totals.std()) / mean_total < 1e-3)

    # Mirror _require_raw_counts without importing the tool: say what the lock will decide.
    if "log1p" in adata.uns:
        verdict, why = "not counts", "uns['log1p'] stamp present (already log-transformed)"
    elif xmin < 0:
        verdict, why = "not counts", f"negative values (min {xmin:.4g}) - looks scaled/z-scored"
    elif integral:
        verdict, why = "raw counts", "non-negative and integral"
    elif xmax < 50:
        verdict, why = "not counts", f"fractional with max only {xmax:.4g} - looks log-transformed"
    elif constant_totals:
        verdict, why = "not counts", f"fractional with near-identical per-cell totals (~{mean_total:.4g})"
    else:
        verdict, why = "count-like", "fractional but large-valued with varying totals (e.g. corrected counts)"

    layers = [k for k in adata.layers.keys() if isinstance(k, str)]
    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "X_type": type(X).__name__,
        "X_dtype": str(getattr(X, "dtype", "?")),
        "X_min": round(xmin, 4),
        "X_max": round(xmax, 4),
        "X_integral": integral,
        "mean_total_per_cell": round(mean_total, 2),
        "layers": layers,
        "has_counts_layer": "counts" in layers,
        "has_raw": adata.raw is not None,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "obsm": list(adata.obsm.keys()),
        "counts_verdict": verdict,
        "counts_reason": why,
    }


def ram_lookahead(info: dict) -> dict:
    cells_x_genes = info["n_obs"] * info["n_vars"]
    est = RAM_FLOOR_MB + RAM_SLOPE_MB_PER_MILLION * (cells_x_genes / 1e6)
    return {
        "cells_x_genes_millions": round(cells_x_genes / 1e6, 1),
        "estimated_peak_mb": round(est),
        "budget_mb": RAM_BUDGET_MB,
        "within_budget": est <= RAM_BUDGET_MB,
    }


def decide(info: dict, ram: dict) -> None:
    """Refuse before doing any work if the file is the wrong kind or too big for this box."""
    if not ram["within_budget"]:
        raise Refusal(
            f"{info['n_obs']:,} cells x {info['n_vars']:,} genes would peak near "
            f"{_fmt_mb(ram['estimated_peak_mb'])}, over the {_fmt_mb(RAM_BUDGET_MB)} budget on this 16 GB box - "
            "subset the cells or genes first."
        )
    if info["counts_verdict"] == "not counts":
        fix = (
            " The file has layers['counts']: set adata.X = adata.layers['counts'] and re-drop it."
            if info["has_counts_layer"]
            else " QC needs raw counts; supply the count matrix."
        )
        raise Refusal(f"X is not raw counts - {info['counts_reason']}.{fix}")


async def run_qc(path: Path, out_dir: Path, mt_prefix: str) -> dict:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
    from fastmcp import Client

    from tools.clustering import clustering_mcp

    async with Client(clustering_mcp) as client:
        result = await client.call_tool(
            "scanpy_compute_qc_and_filter",
            {"data_path": str(path), "mt_prefix": mt_prefix, "output_dir": str(out_dir)},
        )
    return result.data


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inbox", default=str(DEFAULT_INBOX), help=f"folder to read (default {DEFAULT_INBOX})")
    ap.add_argument("--inspect-only", action="store_true", help="report what the file is, then stop")
    ap.add_argument("--mt-prefix", default="MT-", help="mitochondrial gene prefix ('MT-' human, 'mt-' mouse)")
    ap.add_argument("--output-dir", default=None, help="where QC artifacts go (default INBOX/../INBOX_OUT)")
    ap.add_argument("--report", default=None, help="also write the report as JSON here")
    args = ap.parse_args()

    inbox = Path(args.inbox).expanduser().resolve()
    report: dict = {"when": datetime.now().isoformat(timespec="seconds"), "inbox": str(inbox)}
    print(f"door: reading {inbox}")

    try:
        path = find_input(inbox)
        report["file"] = str(path)
        report["file_mb"] = round(check_size(path), 2)
        print(f"door: found {path.name} ({report['file_mb']} MB)")

        # --- inspect first, always ---
        info = inspect(path)
        report["inspect"] = info
        print(f"  shape         : {info['n_obs']:,} cells x {info['n_vars']:,} genes ({info['X_type']}, {info['X_dtype']})")
        print(f"  values        : min {info['X_min']}, max {info['X_max']}, integral={info['X_integral']}")
        print(f"  layers        : {info['layers'] or 'none'} | raw slot: {info['has_raw']}")
        print(f"  obs columns   : {info['obs_columns'] or 'none'}")
        print(f"  precomputed   : {info['obsm'] or 'none'}")
        print(f"  counts verdict: {info['counts_verdict']} ({info['counts_reason']})")

        ram = ram_lookahead(info)
        report["ram_lookahead"] = ram
        print(f"  RAM look-ahead: ~{_fmt_mb(ram['estimated_peak_mb'])} peak "
              f"(budget {_fmt_mb(ram['budget_mb'])}) -> {'ok' if ram['within_budget'] else 'TOO BIG'}")

        if args.inspect_only:
            report["outcome"] = "inspected"
            print("door: --inspect-only, stopping here.")
            return 0

        decide(info, ram)

        out = Path(args.output_dir) if args.output_dir else inbox.parent / "INBOX_OUT"
        out.mkdir(parents=True, exist_ok=True)
        print(f"door: accepted -> running scanpy_compute_qc_and_filter (mt_prefix={args.mt_prefix!r})")
        result = asyncio.run(run_qc(path, out, args.mt_prefix))
        report["outcome"] = "qc_ran"
        report["qc"] = {k: v for k, v in result.items() if k != "reference"}
        print(f"door: {result['message']}")
        for a in result["artifacts"]:
            print(f"  - {a['description']}: {a['path']}")
        print("door: QC only. Nothing was clustered, annotated or named.")
        return 0

    except Refusal as r:
        report["outcome"] = "refused"
        report["reason"] = str(r)
        print(f"door: REFUSED - {r}")
        return 2
    except Exception as exc:  # noqa: BLE001 - a bad file should not look like a crash in the door
        report["outcome"] = "error"
        report["reason"] = f"{type(exc).__name__}: {exc}"
        print(f"door: ERROR reading the file - {type(exc).__name__}: {exc}")
        return 3
    finally:
        if args.report:
            Path(args.report).write_text(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    raise SystemExit(main())
