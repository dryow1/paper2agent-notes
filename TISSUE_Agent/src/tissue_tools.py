"""Two thin wrappers over the TISSUE tutorial path (see NOTES/004-tissue-tutorial.md).

Scope: the smallest useful slice of TISSUE — predict a held-out gene, then calibrate
uncertainty into prediction intervals. Everything else in tissue.main / tissue.downstream
(~34 functions) is deliberately not wrapped.

Shaped like the Scanpy tools in Scanpy_Agent/src/tools/clustering.py:

  * keyword arguments with tutorial defaults,
  * a clear ValueError naming the missing file / gene / key before any computation,
  * a fresh output directory per call,
  * an .h5ad written out so the two tools chain,
  * a dict of absolute artifact paths plus a compact numeric summary.

The wrappers add no maths of their own: every scientific call is upstream TISSUE.

Reference: https://github.com/sunericd/TISSUE (commit ffc3599)
"""

from __future__ import annotations

import os
import sys
import tempfile
import uuid
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REFERENCE = "https://github.com/sunericd/TISSUE/blob/ffc35998004ecf44eaf2aedf774c2f6a2292e7b9/README.md"

#: Prediction methods whose dependencies the pinned environment actually has.
#: tangram / gimvi are commented out in the repo's requirements.txt and are not installed.
SUPPORTED_METHODS = ("spage", "knn")


def _tissue_repo() -> Path:
    """Locate the TISSUE checkout so `import tissue` works from anywhere.

    Order: $TISSUE_REPO, then the sibling repo/TISSUE next to this file's parent.
    """
    candidates = []
    if os.environ.get("TISSUE_REPO"):
        candidates.append(Path(os.environ["TISSUE_REPO"]))
    candidates.append(Path(__file__).resolve().parents[1] / "repo" / "TISSUE")
    for c in candidates:
        if (c / "tissue" / "main.py").is_file():
            return c.resolve()
    raise ValueError(
        "TISSUE checkout not found (looked for tissue/main.py in: "
        + ", ".join(str(c) for c in candidates)
        + "). Clone it with: git clone --depth 1 https://github.com/sunericd/TISSUE.git "
        "<TISSUE_Agent>/repo/TISSUE, or set TISSUE_REPO."
    )


def _import_tissue():
    repo = _tissue_repo()
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    import tissue.main  # noqa: PLC0415

    return tissue.main


def _output_dir(output_dir: str | None, tool: str) -> Path:
    base = Path(output_dir) if output_dir else Path(tempfile.gettempdir()) / "tissue_tool_outputs"
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out = base / f"{tool}_{stamp}_{uuid.uuid4().hex[:8]}"
    out.mkdir(parents=True, exist_ok=False)
    return out.resolve()


def _require_file(path: str, what: str, suffixes: tuple[str, ...] = ()) -> Path:
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"{what} does not exist or is not a file: {path}")
    if suffixes and p.suffix not in suffixes:
        raise ValueError(f"{what} must be one of {suffixes}, got: {path}")
    return p.resolve()


def _require_gene(gene: str, names, what: str) -> None:
    """TISSUE lowercases all gene names internally, so compare lowercased."""
    names = [str(n).lower() for n in names]
    if gene not in names:
        near = [n for n in names if n.startswith(gene[:3])][:8]
        raise ValueError(
            f"target_gene '{gene}' not found in {what} ({len(names)} genes)."
            + (f" Genes starting '{gene[:3]}': {near}." if near else "")
            + " Gene names are compared in lowercase."
        )


def _artifact(description: str, path: Path) -> dict:
    return {"description": description, "path": str(Path(path).resolve())}


def tissue_predict_spatial_gene(
    spatial_counts_path: str,
    locations_path: str,
    scrna_counts_path: str,
    *,
    target_gene: str = "plp1",
    method: str = "spage",
    n_folds: int = 3,
    n_pv: int = 10,
    output_dir: str | None = None,
) -> dict:
    """Predict one held-out spatial gene from paired spatial + scRNA-seq data.

    Holds `target_gene` out of the spatial matrix, predicts it with the chosen method, and
    keeps the measured values in obs['<gene>_measured'] so the calibration step can score
    coverage. Writes an .h5ad to feed tissue_calibrate_prediction_intervals.

    Wraps: tissue.main.load_paired_datasets, preprocess_data, predict_gene_expression.
    """
    tm = _import_tissue()
    import anndata as ad  # noqa: PLC0415

    spatial = _require_file(spatial_counts_path, "spatial_counts_path")
    locations = _require_file(locations_path, "locations_path")
    scrna = _require_file(scrna_counts_path, "scrna_counts_path")
    gene = str(target_gene).lower()
    method = str(method).lower()
    if method not in SUPPORTED_METHODS:
        raise ValueError(
            f"method '{method}' is not available in this environment; supported: {list(SUPPORTED_METHODS)}. "
            "tangram/gimvi need extra packages that requirements.txt leaves commented out."
        )
    if n_folds < 2:
        raise ValueError(f"n_folds must be at least 2, got {n_folds}")

    adata, rna = tm.load_paired_datasets(str(spatial), str(locations), str(scrna))
    adata.var_names = [str(x).lower() for x in adata.var_names]
    rna.var_names = [str(x).lower() for x in rna.var_names]
    _require_gene(gene, adata.var_names, "the spatial dataset")
    _require_gene(gene, rna.var_names, "the scRNA-seq dataset")

    tm.preprocess_data(rna, standardize=False, normalize=True)
    shared = list(np.intersect1d(adata.var_names, rna.var_names))
    if gene not in shared:
        raise ValueError(f"target_gene '{gene}' is not in both datasets; {len(shared)} genes are shared")
    n_conf = len(shared) - 1
    if n_conf < 2:
        raise ValueError(
            f"need at least 3 shared genes (1 target + 2 calibration), found {len(shared)}: {shared}"
        )
    if n_folds > n_conf:
        raise ValueError(
            f"n_folds ({n_folds}) cannot exceed the {n_conf} calibration genes left after holding out '{gene}'"
        )

    out = _output_dir(output_dir, "tissue_predict")
    adata = adata[:, shared].copy()
    measured = np.asarray(adata[:, gene].X).ravel().copy()
    adata = adata[:, [g for g in shared if g != gene]].copy()

    tm.predict_gene_expression(adata, rna, [gene], method=method, n_folds=n_folds, n_pv=n_pv)

    pred_key = f"{method}_predicted_expression"
    adata.obs[f"{gene}_measured"] = measured
    adata.uns["tissue_tool"] = {"target_gene": gene, "method": method, "n_folds": n_folds, "n_pv": n_pv}
    predicted = np.asarray(adata.obsm[pred_key][gene]).ravel()

    h5ad = out / "predicted.h5ad"
    adata.write_h5ad(h5ad)
    return {
        "message": (
            f"Predicted '{gene}' for {adata.n_obs} cells with {method} "
            f"({n_conf} calibration genes, {n_folds} folds)"
        ),
        "reference": REFERENCE,
        "artifacts": [_artifact(f"AnnData with obsm['{pred_key}'] and the measured gene in obs", h5ad)],
        "n_obs": int(adata.n_obs),
        "n_calibration_genes": int(adata.n_vars),
        "target_gene": gene,
        "method": method,
        "prediction_key": pred_key,
        "pearson_r_vs_measured": round(float(np.corrcoef(predicted, measured)[0, 1]), 4),
        "predicted_mean": round(float(predicted.mean()), 4),
        "measured_mean": round(float(measured.mean()), 4),
    }


def tissue_calibrate_prediction_intervals(
    data_path: str,
    *,
    target_gene: str | None = None,
    method: str | None = None,
    alpha_level: float = 0.23,
    n_neighbors: int = 15,
    graph_method: str = "fixed_radius",
    grouping_method: str = "kmeans_gene_cell",
    k: int = 4,
    k2: int = 2,
    output_dir: str | None = None,
) -> dict:
    """Calibrate uncertainty on predictions and turn it into prediction intervals.

    Reports empirical coverage over the calibration genes — the quantity TISSUE actually
    targets — and, when the measured column is present, for the held-out gene too.

    Wraps: tissue.main.build_spatial_graph, conformalize_spatial_uncertainty,
    conformalize_prediction_interval.
    """
    tm = _import_tissue()
    import anndata as ad  # noqa: PLC0415

    path = _require_file(data_path, "data_path", suffixes=(".h5ad",))
    if not 0 < alpha_level < 1:
        raise ValueError(f"alpha_level must be strictly between 0 and 1, got {alpha_level}")

    adata = ad.read_h5ad(path)
    saved = dict(adata.uns.get("tissue_tool", {}))
    method = str(method or saved.get("method") or "spage").lower()
    gene = str(target_gene or saved.get("target_gene") or "plp1").lower()

    pred_key = f"{method}_predicted_expression"
    if pred_key not in adata.obsm:
        raise ValueError(
            f"obsm['{pred_key}'] missing; run tissue_predict_spatial_gene first "
            f"(available obsm keys: {list(adata.obsm.keys())})"
        )
    columns = [str(c).lower() for c in adata.obsm[pred_key].columns]
    if gene not in columns:
        raise ValueError(
            f"target_gene '{gene}' is not a column of obsm['{pred_key}'] "
            f"({len(columns)} predicted genes: {columns[:8]}{'...' if len(columns) > 8 else ''})"
        )

    out = _output_dir(output_dir, "tissue_calibrate")
    tm.build_spatial_graph(adata, method=graph_method, n_neighbors=n_neighbors)
    tm.conformalize_spatial_uncertainty(
        adata, pred_key, calib_genes=adata.var_names, grouping_method=grouping_method, k=k, k2=k2
    )
    tm.conformalize_prediction_interval(
        adata, pred_key, calib_genes=adata.var_names, alpha_level=alpha_level, compute_wasserstein=True
    )

    calib = list(adata.var_names)
    truth_c = np.asarray(adata[:, calib].X)
    lo_c = adata.obsm[f"{pred_key}_lo"][calib].values
    hi_c = adata.obsm[f"{pred_key}_hi"][calib].values
    result = {
        "message": "",
        "reference": REFERENCE,
        "artifacts": [],
        "n_obs": int(adata.n_obs),
        "target_gene": gene,
        "method": method,
        "alpha_level": float(alpha_level),
        "nominal_coverage": round(1 - float(alpha_level), 4),
        "n_calibration_genes": len(calib),
        "coverage_calibration_genes": round(float(((truth_c >= lo_c) & (truth_c <= hi_c)).mean()), 4),
        "mean_interval_width_calibration_genes": round(float((hi_c - lo_c).mean()), 4),
    }

    lo = np.asarray(adata.obsm[f"{pred_key}_lo"][gene]).ravel()
    hi = np.asarray(adata.obsm[f"{pred_key}_hi"][gene]).ravel()
    pred = np.asarray(adata.obsm[pred_key][gene]).ravel()
    measured_col = f"{gene}_measured"
    if measured_col in adata.obs.columns:
        # The held-out gene was not part of calibration, so this is an independent spot check.
        measured = np.asarray(adata.obs[measured_col]).ravel()
        result["coverage_target_gene"] = round(float(((measured >= lo) & (measured <= hi)).mean()), 4)
        result["mean_interval_width_target_gene"] = round(float((hi - lo).mean()), 4)
        order = np.argsort(measured)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.fill_between(range(len(measured)), lo[order], hi[order], alpha=0.3,
                        label=f"{round((1 - alpha_level) * 100)}% PI")
        ax.plot(measured[order], lw=1, label="measured")
        ax.plot(pred[order], lw=0.8, alpha=0.8, label=f"{method} predicted")
        ax.set_xlabel("cells, ordered by measured expression")
        ax.set_ylabel(f"{gene} expression")
        ax.legend(fontsize=8)
        fig.tight_layout()
        figpath = out / f"prediction_intervals_{gene}.png"
        fig.savefig(figpath, dpi=150, facecolor="white")
        plt.close("all")
        result["artifacts"].append(_artifact(f"measured vs predicted {gene} with intervals", figpath))
    else:
        result["coverage_target_gene"] = None
        result["measured_column_note"] = (
            f"obs['{measured_col}'] absent, so held-out coverage and the figure were skipped; "
            "it is written by tissue_predict_spatial_gene"
        )

    h5ad = out / "calibrated.h5ad"
    adata.write_h5ad(h5ad)
    result["artifacts"].insert(0, _artifact(f"AnnData with obsm['{pred_key}_lo'/'_hi']", h5ad))
    result["message"] = (
        f"Calibrated {round((1 - alpha_level) * 100)}% intervals: "
        f"coverage {result['coverage_calibration_genes']} over {len(calib)} calibration genes"
    )
    return result
