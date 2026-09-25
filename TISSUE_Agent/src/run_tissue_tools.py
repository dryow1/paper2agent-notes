"""Exercise the two TISSUE wrappers on the repo's built-in 420 KB dataset.

Runs the happy path (predict -> calibrate) and a set of wrong-input cases, recording the
error message each one produces. Writes results to $TISSUE_OUT (default TISSUE_Agent/out).
"""

from __future__ import annotations

import json
import os
import resource
import time
from pathlib import Path

import warnings

warnings.filterwarnings("ignore")

from tissue_tools import (  # noqa: E402
    _tissue_repo,
    tissue_calibrate_prediction_intervals,
    tissue_predict_spatial_gene,
)

OUT = Path(os.environ.get("TISSUE_OUT", Path(__file__).resolve().parents[1] / "out")).resolve()
OUT.mkdir(parents=True, exist_ok=True)
DATA = _tissue_repo() / "tests" / "data"
SPATIAL, LOC, SCRNA = DATA / "Spatial_count.txt", DATA / "Locations.txt", DATA / "scRNA_count.txt"
REC: dict = {"data_dir": str(DATA), "calls": [], "errors": []}


def rss() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


def run(label: str, fn, **kw) -> dict | None:
    t0 = time.time()
    try:
        res = fn(**kw)
        REC["calls"].append(
            {"label": label, "seconds": round(time.time() - t0, 1), "maxrss_mb": rss(),
             "summary": {k: v for k, v in res.items() if k != "artifacts"},
             "artifacts": res["artifacts"]}
        )
        print(f"[{label}] ok in {REC['calls'][-1]['seconds']}s (maxrss {rss()} MB)\n    {res['message']}", flush=True)
        return res
    except Exception as exc:  # noqa: BLE001 - the point is to record the message
        REC["calls"].append({"label": label, "failed": f"{type(exc).__name__}: {exc}"})
        print(f"[{label}] FAILED: {type(exc).__name__}: {exc}", flush=True)
        return None


def expect_error(label: str, fn, **kw) -> None:
    """A wrong-input case: record the message the wrapper produces."""
    try:
        fn(**kw)
    except ValueError as exc:
        REC["errors"].append({"case": label, "error_type": "ValueError", "message": str(exc)})
        print(f"[error case: {label}]\n    ValueError: {exc}", flush=True)
    except Exception as exc:  # noqa: BLE001
        REC["errors"].append({"case": label, "error_type": type(exc).__name__, "message": str(exc)})
        print(f"[error case: {label}] NOT a ValueError -> {type(exc).__name__}: {exc}", flush=True)
    else:
        REC["errors"].append({"case": label, "error_type": None, "message": "NO ERROR RAISED"})
        print(f"[error case: {label}] NO ERROR RAISED", flush=True)


# --- happy path --------------------------------------------------------------------------
pred = run(
    "predict", tissue_predict_spatial_gene,
    spatial_counts_path=str(SPATIAL), locations_path=str(LOC), scrna_counts_path=str(SCRNA),
    target_gene="plp1", method="spage", n_folds=3, n_pv=10, output_dir=str(OUT / "tools"),
)
cal = None
if pred is not None:
    cal = run(
        "calibrate", tissue_calibrate_prediction_intervals,
        data_path=pred["artifacts"][0]["path"], alpha_level=0.23, output_dir=str(OUT / "tools"),
    )

# A second resolution of the same switch, to check the keyword actually does something.
if pred is not None:
    run(
        "calibrate_alpha0.5", tissue_calibrate_prediction_intervals,
        data_path=pred["artifacts"][0]["path"], alpha_level=0.5, output_dir=str(OUT / "tools"),
    )

# --- wrong-input cases -------------------------------------------------------------------
expect_error(
    "missing spatial file", tissue_predict_spatial_gene,
    spatial_counts_path=str(DATA / "nope.txt"), locations_path=str(LOC), scrna_counts_path=str(SCRNA),
)
expect_error(
    "gene not in data", tissue_predict_spatial_gene,
    spatial_counts_path=str(SPATIAL), locations_path=str(LOC), scrna_counts_path=str(SCRNA),
    target_gene="notagene",
)
expect_error(
    "unsupported method", tissue_predict_spatial_gene,
    spatial_counts_path=str(SPATIAL), locations_path=str(LOC), scrna_counts_path=str(SCRNA),
    method="tangram",
)
expect_error(
    "n_folds above calibration genes", tissue_predict_spatial_gene,
    spatial_counts_path=str(SPATIAL), locations_path=str(LOC), scrna_counts_path=str(SCRNA),
    n_folds=999,
)
expect_error(
    "calibrate on un-predicted h5ad", tissue_calibrate_prediction_intervals,
    data_path=str(SPATIAL),
)
if pred is not None:
    expect_error(
        "alpha out of range", tissue_calibrate_prediction_intervals,
        data_path=pred["artifacts"][0]["path"], alpha_level=1.5,
    )
    expect_error(
        "wrong method key", tissue_calibrate_prediction_intervals,
        data_path=pred["artifacts"][0]["path"], method="knn",
    )

REC["maxrss_mb_final"] = rss()
REC["n_ok"] = sum("summary" in c for c in REC["calls"])
REC["n_failed"] = sum("failed" in c for c in REC["calls"])
REC["n_error_cases_clean"] = sum(e["error_type"] == "ValueError" for e in REC["errors"])
REC["n_error_cases"] = len(REC["errors"])
(OUT / "tools_results.json").write_text(json.dumps(REC, indent=2, default=str))
print(
    f"\n{REC['n_ok']} ok, {REC['n_failed']} failed; "
    f"{REC['n_error_cases_clean']}/{REC['n_error_cases']} wrong-input cases gave a clean ValueError; "
    f"maxrss {REC['maxrss_mb_final']} MB",
    flush=True,
)
