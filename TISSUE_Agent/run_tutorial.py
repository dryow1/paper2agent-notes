"""Run the official TISSUE tutorial path on the repo's own tests/data and record real numbers.

Same calls as the repo's test.py (and README tutorials 1-3), but instead of only asserting
"no exception", this records the quantities a user would actually read: prediction accuracy
for the held-out gene, empirical coverage of the calibrated prediction intervals, the
multiple-imputation t-test result, and how many cells the uncertainty filter keeps.

Run from the repo root (tissue/ and tests/ must be in the cwd).
"""

from __future__ import annotations

import json
import os
import resource
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import warnings  # noqa: E402

warnings.filterwarnings("ignore")

import tissue.downstream  # noqa: E402
import tissue.main  # noqa: E402

OUT = Path(os.environ.get("TISSUE_OUT", "../../out")).resolve()
OUT.mkdir(parents=True, exist_ok=True)
TARGET = "plp1"
METHOD = "spage"
PRED_KEY = f"{METHOD}_predicted_expression"
REC: dict = {"target_gene": TARGET, "method": METHOD, "steps": []}


def step(name: str, t0: float) -> None:
    REC["steps"].append(
        {"step": name, "seconds": round(time.time() - t0, 1),
         "maxrss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024}
    )
    print(f"[{name}] {REC['steps'][-1]['seconds']}s, maxrss {REC['steps'][-1]['maxrss_mb']} MB", flush=True)


# --- Tutorial 1: load, preprocess, predict -------------------------------------------------
t = time.time()
adata, rna = tissue.main.load_paired_datasets(
    "tests/data/Spatial_count.txt", "tests/data/Locations.txt", "tests/data/scRNA_count.txt"
)
adata.var_names = [x.lower() for x in adata.var_names]
rna.var_names = [x.lower() for x in rna.var_names]
tissue.main.preprocess_data(rna, standardize=False, normalize=True)
shared = np.intersect1d(adata.var_names, rna.var_names)
adata = adata[:, shared].copy()
REC["spatial_shape"] = list(adata.shape)
REC["rnaseq_shape"] = list(rna.shape)
REC["n_shared_genes"] = int(len(shared))
# Hold out the target gene so its true values can score the prediction.
truth = np.asarray(adata[:, TARGET].X).ravel().copy()
adata = adata[:, [g for g in shared if g != TARGET]].copy()
step("load+preprocess", t)

t = time.time()
tissue.main.predict_gene_expression(adata, rna, [TARGET], method=METHOD, n_folds=3, n_pv=10)
pred = np.asarray(adata.obsm[PRED_KEY][TARGET]).ravel()
REC["prediction"] = {
    "n_folds": 3,
    "n_pv": 10,
    "pearson_r_vs_truth": round(float(np.corrcoef(pred, truth)[0, 1]), 4),
    "spearman_r_vs_truth": round(
        float(np.corrcoef(np.argsort(np.argsort(pred)), np.argsort(np.argsort(truth)))[0, 1]), 4
    ),
    "pred_mean": round(float(pred.mean()), 4),
    "truth_mean": round(float(truth.mean()), 4),
}
step("predict", t)

# --- Tutorial 2: spatial graph, calibration, prediction intervals ---------------------------
t = time.time()
tissue.main.build_spatial_graph(adata, method="fixed_radius", n_neighbors=15)
tissue.main.conformalize_spatial_uncertainty(
    adata, PRED_KEY, calib_genes=adata.var_names, grouping_method="kmeans_gene_cell", k=4, k2=2
)
ALPHA = 0.23
tissue.main.conformalize_prediction_interval(
    adata, PRED_KEY, calib_genes=adata.var_names, alpha_level=ALPHA, compute_wasserstein=True
)
lo = np.asarray(adata.obsm[f"{PRED_KEY}_lo"][TARGET]).ravel()
hi = np.asarray(adata.obsm[f"{PRED_KEY}_hi"][TARGET]).ravel()
covered = (truth >= lo) & (truth <= hi)
# TISSUE calibrates over the measured (calibration) genes, not the held-out one, and builds
# residuals as adata.X - predicted (main.py:652), i.e. against raw counts. So the aggregate
# coverage over calibration genes is the number TISSUE actually targets; the held-out gene is
# a single-gene spot check that can deviate.
calib = [g for g in adata.var_names]
truth_c = np.asarray(adata[:, calib].X)
lo_c = adata.obsm[f"{PRED_KEY}_lo"][calib].values
hi_c = adata.obsm[f"{PRED_KEY}_hi"][calib].values
cov_c = (truth_c >= lo_c) & (truth_c <= hi_c)
REC["prediction_intervals"] = {
    "alpha_level": ALPHA,
    "nominal_coverage": round(1 - ALPHA, 2),
    "empirical_coverage_target_gene": round(float(covered.mean()), 4),
    "empirical_coverage_calibration_genes": round(float(cov_c.mean()), 4),
    "n_calibration_genes": int(len(calib)),
    "mean_interval_width_target_gene": round(float((hi - lo).mean()), 4),
    "mean_interval_width_calibration_genes": round(float((hi_c - lo_c).mean()), 4),
    "n_cells": int(len(truth)),
    # Scale note: SpaGE predictions and the measured counts are not on the same scale
    # (see *_scale_check below); TISSUE subtracts them anyway, which widens the intervals.
    "scale_check": {
        "measured_target_mean": round(float(truth.mean()), 4),
        "measured_target_max": round(float(truth.max()), 4),
        "predicted_target_mean": round(float(pred.mean()), 4),
        "predicted_target_max": round(float(pred.max()), 4),
    },
}
step("calibrate+intervals", t)

# --- Tutorial 3: multiple-imputation hypothesis testing -------------------------------------
t = time.time()
adata.obs["condition"] = ["A" if i < round(adata.shape[0] / 2) else "B" for i in range(adata.shape[0])]
keys = tissue.downstream.multiple_imputation_testing(
    adata, PRED_KEY, calib_genes=adata.var_names, condition="condition",
    group1="A", group2="B", n_imputations=10, return_keys=True,
)
mi = {}
for k in (keys if isinstance(keys, (list, tuple)) else [keys]):
    for store in (adata.uns, adata.var, adata.obsm):
        if k in store:
            v = store[k]
            try:
                mi[str(k)] = round(float(np.asarray(v[TARGET]).ravel()[0]), 6)
            except Exception:  # noqa: BLE001 - key shape varies by test type
                mi[str(k)] = str(type(v))
REC["multiple_imputation_ttest"] = {"n_imputations": 10, "result_keys": [str(k) for k in np.atleast_1d(keys)], "values": mi}
step("mi_ttest", t)

# --- Cell filtering + filtered PCA ----------------------------------------------------------
t = time.time()
unc = adata.obsm[f"{PRED_KEY}_hi"].values - adata.obsm[f"{PRED_KEY}_lo"].values
keep = tissue.downstream.detect_uncertain_cells(
    unc, proportion="otsu", stratification=adata.obs["condition"].values
)
keep_pca = tissue.downstream.filtered_PCA(
    adata, METHOD, proportion="otsu", stratification=adata.obs["condition"].values, return_keep_idxs=True
)
# Both return row INDICES (not a boolean mask), so count them with len().
REC["cell_filtering"] = {
    "proportion": "otsu",
    "n_cells": int(adata.shape[0]),
    "n_kept_detect_uncertain_cells": int(len(keep)),
    "n_kept_filtered_PCA": int(len(keep_pca)),
    "pct_dropped": round(100 * (1 - len(keep) / adata.shape[0]), 1),
}
step("cell_filtering", t)

# One figure: truth vs prediction with calibrated intervals, ordered by true expression.
order = np.argsort(truth)
fig, ax = plt.subplots(figsize=(8, 4))
ax.fill_between(range(len(truth)), lo[order], hi[order], alpha=0.3, label=f"{int((1-ALPHA)*100)}% PI")
ax.plot(truth[order], lw=1, label="measured")
ax.plot(pred[order], lw=0.8, alpha=0.8, label="SpaGE predicted")
ax.set_xlabel("cells, ordered by measured expression")
ax.set_ylabel(f"{TARGET} expression")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "plp1_prediction_intervals.png", dpi=150, facecolor="white")
plt.close("all")

REC["maxrss_mb_final"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
(OUT / "tutorial_results.json").write_text(json.dumps(REC, indent=2))
print(json.dumps({k: v for k, v in REC.items() if k != "steps"}, indent=2))
print(f"\nwrote {OUT}/tutorial_results.json and plp1_prediction_intervals.png", flush=True)
