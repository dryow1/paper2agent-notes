"""Tools extracted from scanpy/docs/tutorials/basics/clustering.ipynb (Preprocessing and clustering).

Each tool reads an AnnData file, calls the public scanpy API with the tutorial's
arguments, writes the updated AnnData plus the tutorial's figures into a fresh
output directory, and returns absolute artifact paths with a compact summary.
"""

from __future__ import annotations

import contextlib
import sys
import tempfile
import uuid
from datetime import datetime
from pathlib import Path
from typing import Annotated, Literal

import matplotlib

matplotlib.use("Agg")
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from fastmcp import FastMCP  # noqa: E402

import scanpy as sc  # noqa: E402

REFERENCE = "https://github.com/scverse/scanpy/blob/0d5fd16234865619d2f5097d33fc4281900a2bc2/docs/tutorials/basics/clustering.ipynb"
clustering_mcp = FastMCP(name="clustering")


def _output_dir(output_dir: str | None, tool: str) -> Path:
    """Create a fresh, collision-resistant directory for one tool invocation."""
    base = Path(output_dir) if output_dir else Path(tempfile.gettempdir()) / "scanpy_mcp_outputs"
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out = base / f"{tool}_{stamp}_{uuid.uuid4().hex[:8]}"
    out.mkdir(parents=True, exist_ok=False)
    return out.resolve()


def _read(data_path: str, allow_10x_h5: bool = False) -> ad.AnnData:
    path = Path(data_path)
    if not path.is_file():
        raise ValueError(f"data_path does not exist or is not a file: {data_path}")
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    if allow_10x_h5 and path.suffix == ".h5":
        # Tutorial cell 4: sc.read_10x_h5 followed by var_names_make_unique
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        return adata
    allowed = ".h5ad or 10x Genomics .h5" if allow_10x_h5 else ".h5ad"
    raise ValueError(f"data_path must be an {allowed} file, got: {data_path}")


def _write(adata: ad.AnnData, out: Path, name: str) -> Path:
    path = out / name
    adata.write_h5ad(path, compression="gzip")
    return path


def _save_figure(path: Path) -> Path:
    """Save the figure left open by a scanpy plotting call (show=False)."""
    plt.gcf().savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    return path


def _require_obs(adata: ad.AnnData, keys: list[str], what: str) -> None:
    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        raise ValueError(f"{what} not found in adata.obs: {missing}; available: {list(adata.obs.columns)}")


def _require_color_keys(adata: ad.AnnData, keys: list[str]) -> None:
    missing = [k for k in keys if k not in adata.obs.columns and k not in adata.var_names]
    if missing:
        raise ValueError(f"color keys not found in adata.obs columns or var_names: {missing}")


def _counts_hint(adata: ad.AnnData) -> str:
    """Point at a counts layer / raw slot when this object actually has one."""
    if "counts" in adata.layers:
        return " This object has layers['counts']: set adata.X = adata.layers['counts'] and retry."
    if adata.raw is not None:
        return " This object has a .raw slot, but check it holds counts (it is often log-normalized)."
    return " Supply the raw count matrix instead."


def _require_raw_counts(adata: ad.AnnData) -> None:
    """Reject matrices that are demonstrably not raw counts.

    sc.pp.calculate_qc_metrics just sums X per cell, so it returns a plausible-looking
    number for any input: scaled data give negative 'total_counts', log-normalized data
    give totals that are not counts at all, and neither is flagged.

    Each rule below keys on a positive signature of a transformation rather than on
    "not integral", so ambient-corrected or otherwise fractional count matrices - which
    are still counts - are allowed through.
    """
    X = adata.X
    values = X.data if sp.issparse(X) else np.asarray(X).reshape(-1)
    if values.size == 0:
        return

    # 1. scanpy stamps uns['log1p'] when it log-transforms. Definitive.
    if "log1p" in adata.uns:
        raise ValueError(
            "adata.uns['log1p'] is present: the data have already been log-transformed, so "
            "QC metrics computed from them are not counts." + _counts_hint(adata)
        )

    # 2. Negative values mean scaled/z-scored.
    xmin = float(values.min())
    if xmin < 0:
        raise ValueError(
            f"X contains negative values (min {xmin:.4g}): the data look scaled/z-scored, "
            "not raw counts, and QC metrics computed from them are meaningless."
            + _counts_hint(adata)
        )

    if not np.any(values % 1):
        return  # integral and non-negative: counts

    xmax = float(values.max())
    # 3. Fractional values with a small ceiling are the signature of log1p: log1p of even
    #    100,000 counts is only 11.5, whereas a real count matrix with a max below 50 is integral.
    if xmax < 50:
        raise ValueError(
            f"X has fractional values and a maximum of only {xmax:.4g}: the data look "
            "log-transformed, not raw counts, so 'total_counts' would not be a count."
            + _counts_hint(adata)
        )

    # 4. Fractional values with near-constant per-cell totals mean normalize_total() was run.
    totals = np.asarray(X.sum(axis=1)).reshape(-1)
    mean_total = float(totals.mean())
    if mean_total > 0 and float(totals.std()) / mean_total < 1e-3:
        raise ValueError(
            f"X has fractional values and near-identical per-cell totals (~{mean_total:.4g}): "
            "the data look normalized to a fixed target sum, not raw counts."
            + _counts_hint(adata)
        )


def _require_neighbors(adata: ad.AnnData) -> None:
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        raise ValueError(
            "neighbor graph missing (uns['neighbors'] / obsp['connectivities']); run scanpy_build_neighbors_umap first"
        )


def _artifact(description: str, path: Path) -> dict:
    return {"description": description, "path": str(path.resolve())}


@clustering_mcp.tool()
def scanpy_compute_qc_and_filter(
    data_path: Annotated[str, "Raw-count AnnData .h5ad (cells x genes) or a 10x Genomics .h5 matrix"],
    mt_prefix: Annotated[str, "Prefix marking mitochondrial genes ('MT-' human, 'mt-' mouse)"] = "MT-",
    ribo_prefixes: Annotated[list[str], "Prefixes marking ribosomal genes"] = ["RPS", "RPL"],  # noqa: B006
    hb_pattern: Annotated[str, "Regular expression marking hemoglobin genes"] = "^HB[^(P)]",
    min_genes: Annotated[int, "Keep cells with at least this many detected genes"] = 100,
    min_cells: Annotated[int, "Keep genes detected in at least this many cells"] = 3,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Annotate mt/ribo/hb genes, compute QC metrics and plots, then filter low-quality cells and rare genes.
    Raw-count .h5ad/.h5 → QC-annotated filtered .h5ad plus QC violin and scatter figures.
    """
    adata = _read(data_path, allow_10x_h5=True)
    _require_raw_counts(adata)
    out = _output_dir(output_dir, "qc_filter")
    n_obs_raw, n_vars_raw = adata.shape
    with contextlib.redirect_stdout(sys.stderr):
        # mitochondrial, ribosomal and hemoglobin genes (tutorial cell 9)
        adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
        adata.var["ribo"] = adata.var_names.str.startswith(tuple(ribo_prefixes))
        adata.var["hb"] = adata.var_names.str.contains(hb_pattern)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
        # QC plots are drawn before filtering, as in the tutorial (cells 12, 14)
        sc.pl.violin(
            adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, show=False
        )
        violin = _save_figure(out / "qc_violin.png")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        scatter = _save_figure(out / "qc_scatter.png")
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
    h5ad = _write(adata, out, "qc_filtered.h5ad")
    return {
        "message": f"QC computed; kept {adata.n_obs}/{n_obs_raw} cells and {adata.n_vars}/{n_vars_raw} genes",
        "reference": REFERENCE,
        "artifacts": [
            _artifact("QC-annotated, filtered AnnData", h5ad),
            _artifact("violin plot of n_genes_by_counts, total_counts, pct_counts_mt (before filtering)", violin),
            _artifact("scatter total_counts vs n_genes_by_counts colored by pct_counts_mt (before filtering)", scatter),
        ],
        "n_obs_before": int(n_obs_raw),
        "n_vars_before": int(n_vars_raw),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_mt_genes": int(adata.var["mt"].sum()),
        "n_ribo_genes": int(adata.var["ribo"].sum()),
        "n_hb_genes": int(adata.var["hb"].sum()),
    }


@clustering_mcp.tool()
def scanpy_detect_doublets_scrublet(
    data_path: Annotated[str, "AnnData .h5ad with raw counts in .X (e.g. QC-filtered output)"],
    batch_key: Annotated[str | None, "Existing obs column with sample/batch labels; Scrublet runs per batch (tutorial: 'sample')"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Score and call doublets with Scrublet, optionally separately for each batch.
    Raw-count .h5ad → .h5ad with obs doublet_score/predicted_doublet and uns['scrublet'].
    """
    adata = _read(data_path)
    if batch_key is not None:
        _require_obs(adata, [batch_key], "batch_key")
    out = _output_dir(output_dir, "scrublet")
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cell 18; predicted doublets are annotated, not removed
        sc.pp.scrublet(adata, batch_key=batch_key)
    h5ad = _write(adata, out, "doublets.h5ad")
    result = {
        "message": f"Scrublet predicted {int(adata.obs['predicted_doublet'].sum())} doublets among {adata.n_obs} cells",
        "reference": REFERENCE,
        "artifacts": [_artifact("AnnData with doublet_score and predicted_doublet in obs", h5ad)],
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_predicted_doublets": int(adata.obs["predicted_doublet"].sum()),
        "doublet_rate": float(adata.obs["predicted_doublet"].mean()),
    }
    if batch_key is not None:
        result["predicted_doublets_per_batch"] = {
            str(k): int(v)
            for k, v in adata.obs.groupby(batch_key, observed=True)["predicted_doublet"].sum().items()
        }
    return result


@clustering_mcp.tool()
def scanpy_normalize_log1p(
    data_path: Annotated[str, "AnnData .h5ad with raw counts in .X"],
    target_sum: Annotated[float | None, "Total counts per cell after scaling; None uses the median total count (tutorial)"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Keep raw counts in layers['counts'], then normalize total counts per cell and log1p-transform.
    Raw-count .h5ad → log-normalized .h5ad with layers['counts'] and uns['log1p'].
    """
    adata = _read(data_path)
    if "log1p" in adata.uns:
        raise ValueError("adata.uns['log1p'] exists: data already appear log-transformed; supply raw counts")
    if "counts" in adata.layers:
        raise ValueError("adata.layers['counts'] already exists; refusing to overwrite it")
    out = _output_dir(output_dir, "normalize")
    with contextlib.redirect_stdout(sys.stderr):
        # Saving count data (tutorial cell 22)
        adata.layers["counts"] = adata.X.copy()
        # Normalizing to median total counts, then logarithmize (tutorial cell 23)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    h5ad = _write(adata, out, "normalized.h5ad")
    return {
        "message": f"Normalized and log1p-transformed {adata.n_obs} cells x {adata.n_vars} genes",
        "reference": REFERENCE,
        "artifacts": [_artifact("log-normalized AnnData (raw counts in layers['counts'])", h5ad)],
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "target_sum": "median" if target_sum is None else float(target_sum),
    }


@clustering_mcp.tool()
def scanpy_select_highly_variable_genes(
    data_path: Annotated[str, "Log-normalized AnnData .h5ad"],
    n_top_genes: Annotated[int, "Number of highly variable genes to select"] = 2000,
    batch_key: Annotated[str | None, "Existing obs column for batch-aware selection (tutorial: 'sample')"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Annotate highly variable genes with the Seurat dispersion method and plot dispersion versus mean.
    Log-normalized .h5ad → .h5ad with var['highly_variable'] etc. plus the HVG figure.
    """
    adata = _read(data_path)
    if batch_key is not None:
        _require_obs(adata, [batch_key], "batch_key")
    out = _output_dir(output_dir, "hvg")
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cells 25-26; flavor pinned to the ScanpyV1 default used by the tutorial
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key, flavor="seurat")
        sc.pl.highly_variable_genes(adata, show=False)
        fig = _save_figure(out / "highly_variable_genes.png")
    h5ad = _write(adata, out, "hvg.h5ad")
    n_hvg = int(adata.var["highly_variable"].sum())
    return {
        "message": f"Selected {n_hvg} highly variable genes of {adata.n_vars}",
        "reference": REFERENCE,
        "artifacts": [
            _artifact("AnnData with highly variable gene annotation in var", h5ad),
            _artifact("dispersion vs mean expression of highly variable genes", fig),
        ],
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_highly_variable": n_hvg,
        "flavor": "seurat",
    }


@clustering_mcp.tool()
def scanpy_run_pca(
    data_path: Annotated[str, "Log-normalized AnnData .h5ad with var['highly_variable']"],
    n_comps: Annotated[int | None, "Number of principal components; None uses the scanpy default (50)"] = None,
    color_keys: Annotated[list[str] | None, "obs columns or genes to color PC scatter plots (tutorial: ['sample', 'pct_counts_mt'])"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Run PCA on highly variable genes and plot the variance ratio and PC1-2 / PC3-4 scatter plots.
    HVG-annotated .h5ad → .h5ad with obsm['X_pca'], varm['PCs'], uns['pca'] plus figures.
    """
    adata = _read(data_path)
    if "highly_variable" not in adata.var.columns:
        raise ValueError("var['highly_variable'] missing; run scanpy_select_highly_variable_genes first")
    color_keys = list(color_keys or [])
    _require_color_keys(adata, color_keys)
    out = _output_dir(output_dir, "pca")
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cell 28: PCA on the highly variable genes (default mask)
        sc.tl.pca(adata, n_comps=n_comps)
        n_pcs = adata.obsm["X_pca"].shape[1]
        sc.pl.pca_variance_ratio(adata, n_pcs=min(50, n_pcs), log=True, show=False)
        variance = _save_figure(out / "pca_variance_ratio.png")
        # Tutorial cell 32: each key shown on PC1/PC2 and PC3/PC4
        if color_keys:
            sc.pl.pca(
                adata,
                color=[k for k in color_keys for _ in range(2)],
                dimensions=[(0, 1), (2, 3)] * len(color_keys),
                ncols=2,
                size=2,
                show=False,
            )
        else:
            sc.pl.pca(adata, dimensions=[(0, 1), (2, 3)], ncols=2, size=2, show=False)
        scatter = _save_figure(out / "pca.png")
    h5ad = _write(adata, out, "pca.h5ad")
    ratio = adata.uns["pca"]["variance_ratio"]
    return {
        "message": f"PCA computed with {n_pcs} components on {int(adata.var['highly_variable'].sum())} HVGs",
        "reference": REFERENCE,
        "artifacts": [
            _artifact("AnnData with PCA results", h5ad),
            _artifact("variance ratio per PC (log scale)", variance),
            _artifact("PC1/PC2 and PC3/PC4 scatter plots", scatter),
        ],
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_comps": int(n_pcs),
        "variance_ratio_first5": [round(float(v), 5) for v in ratio[:5]],
        "variance_ratio_total": float(ratio.sum()),
    }


@clustering_mcp.tool()
def scanpy_build_neighbors_umap(
    data_path: Annotated[str, "AnnData .h5ad with obsm['X_pca']"],
    n_neighbors: Annotated[int, "Number of nearest neighbors in the kNN graph"] = 15,
    n_pcs: Annotated[int | None, "Number of PCs to use; None uses all of X_pca"] = None,
    color_keys: Annotated[list[str] | None, "obs columns or genes to color the UMAP (tutorial: ['sample'])"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Build the kNN neighborhood graph from PCA and embed it in two dimensions with UMAP.
    PCA .h5ad → .h5ad with obsp distances/connectivities and obsm['X_umap'] plus a UMAP figure.
    """
    adata = _read(data_path)
    if "X_pca" not in adata.obsm:
        raise ValueError("obsm['X_pca'] missing; run scanpy_run_pca first")
    color_keys = list(color_keys or [])
    _require_color_keys(adata, color_keys)
    out = _output_dir(output_dir, "neighbors_umap")
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cells 34, 36, 38
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=color_keys or None, size=2, show=False)
        fig = _save_figure(out / "umap.png")
    h5ad = _write(adata, out, "neighbors_umap.h5ad")
    return {
        "message": f"kNN graph (n_neighbors={n_neighbors}) and UMAP computed for {adata.n_obs} cells",
        "reference": REFERENCE,
        "artifacts": [
            _artifact("AnnData with neighbor graph and UMAP embedding", h5ad),
            _artifact("UMAP embedding", fig),
        ],
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_neighbors": int(n_neighbors),
        "umap_shape": list(adata.obsm["X_umap"].shape),
    }


@clustering_mcp.tool()
def scanpy_cluster_leiden(
    data_path: Annotated[str, "AnnData .h5ad with a neighbor graph (and X_umap for plots)"],
    resolutions: Annotated[list[float], "Leiden resolution(s); tutorial uses [1.0] and then [0.02, 0.5, 2.0]"] = [1.0],  # noqa: B006
    n_iterations: Annotated[int, "Leiden iterations; tutorial uses 2 for the single run and -1 (until convergence) for multi-resolution"] = 2,
    key_added: Annotated[str, "obs key for a single resolution; with several, keys are '{key_added}_res_{res:4.2f}'"] = "leiden",
    qc_color_keys: Annotated[list[str] | None, "Extra obs columns/genes shown next to clusters on UMAP (e.g. predicted_doublet, pct_counts_mt)"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Cluster cells with the igraph Leiden algorithm at one or more resolutions and plot clusters on UMAP.
    Neighbor-graph .h5ad → .h5ad with categorical obs cluster column(s) plus UMAP figures.
    """
    if not resolutions:
        raise ValueError("resolutions must contain at least one value")
    adata = _read(data_path)
    _require_neighbors(adata)
    qc_color_keys = list(qc_color_keys or [])
    _require_color_keys(adata, qc_color_keys)
    if len(resolutions) == 1:
        keys = [key_added]
    else:
        keys = [f"{key_added}_res_{res:4.2f}" for res in resolutions]
    out = _output_dir(output_dir, "leiden")
    artifacts = []
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cells 41 and 52 (flavor pinned to 'igraph' as in the tutorial)
        for res, key in zip(resolutions, keys, strict=True):
            sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=n_iterations)
        if "X_umap" in adata.obsm:
            # Tutorial cells 42 / 54
            if len(keys) == 1:
                sc.pl.umap(adata, color=keys, show=False)
            else:
                sc.pl.umap(adata, color=keys, legend_loc="on data", show=False)
            artifacts.append(_artifact("UMAP colored by Leiden clusters", _save_figure(out / "umap_leiden.png")))
            if qc_color_keys:
                # Tutorial cells 45-46: re-assess QC metrics / doublets per cluster
                sc.pl.umap(adata, color=[keys[0], *qc_color_keys], wspace=0.5, ncols=2, show=False)
                artifacts.append(
                    _artifact(
                        "UMAP of clusters next to QC keys", _save_figure(out / "umap_qc_reassessment.png")
                    )
                )
    h5ad = _write(adata, out, "clustered.h5ad")
    artifacts.insert(0, _artifact("AnnData with Leiden cluster labels in obs", h5ad))
    n_clusters = {k: int(adata.obs[k].nunique()) for k in keys}
    return {
        "message": f"Leiden clustering done: {n_clusters}"
        + ("" if "X_umap" in adata.obsm else " (no obsm['X_umap']; plots skipped)"),
        "reference": REFERENCE,
        "artifacts": artifacts,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "cluster_keys": keys,
        "n_clusters": n_clusters,
        "cluster_sizes": {k: {str(c): int(n) for c, n in adata.obs[k].value_counts().items()} for k in keys},
    }


@clustering_mcp.tool()
def scanpy_plot_marker_gene_dotplot(
    data_path: Annotated[str, "Log-normalized AnnData .h5ad with a categorical obs grouping (e.g. Leiden clusters)"],
    marker_genes: Annotated[dict[str, list[str]], "Marker gene sets: cell-type label -> list of genes in var_names"],
    groupby: Annotated[str, "Existing categorical obs column to group cells by (e.g. 'leiden_res_0.02')"],
    cluster_labels: Annotated[dict[str, str] | None, "Optional mapping of groupby category -> cell-type label written to obs"] = None,
    annotation_key: Annotated[str, "New obs column for cluster_labels (tutorial: 'cell_type_lvl1')"] = "cell_type",
    ignore_missing_genes: Annotated[bool, "Drop marker genes absent from var_names instead of failing; dropped genes and emptied sets are reported"] = False,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Plot mean expression of marker gene sets per cluster and optionally annotate clusters with given labels.
    Clustered .h5ad + marker dict → marker dotplot figure (and annotated .h5ad if labels are given).
    """
    if not marker_genes:
        raise ValueError("marker_genes must contain at least one gene set")
    adata = _read(data_path)
    _require_obs(adata, [groupby], "groupby")
    genes = [g for gs in marker_genes.values() for g in gs]
    missing = sorted({g for g in genes if g not in adata.var_names})
    dropped_genes: list[str] = []
    emptied_sets: list[str] = []
    if missing and not ignore_missing_genes:
        raise ValueError(
            f"marker genes not found in var_names: {missing}"
            + (f"; pass ignore_missing_genes=True to drop these {len(missing)} genes and plot the rest"
               if len(missing) < len(set(genes)) else "")
        )
    if missing:
        # Keep the caller's panel intact apart from what this dataset cannot show, and say so.
        dropped_genes = missing
        kept = {k: [g for g in gs if g in adata.var_names] for k, gs in marker_genes.items()}
        emptied_sets = sorted(k for k, gs in kept.items() if not gs)
        marker_genes = {k: gs for k, gs in kept.items() if gs}
        if not marker_genes:
            raise ValueError(
                f"no marker gene is present in var_names; all {len(emptied_sets)} sets are empty "
                f"after dropping {len(dropped_genes)} missing genes"
            )
    categories = [str(c) for c in adata.obs[groupby].astype("category").cat.categories]
    if cluster_labels is not None:
        unknown = sorted(set(cluster_labels) - set(categories))
        if unknown:
            raise ValueError(f"cluster_labels keys not among obs['{groupby}'] categories {categories}: {unknown}")
        if annotation_key in adata.obs.columns:
            raise ValueError(f"obs['{annotation_key}'] already exists; choose another annotation_key")
    out = _output_dir(output_dir, "marker_dotplot")
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in groupby)
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cells 60 / 63
        sc.pl.dotplot(adata, marker_genes, groupby=groupby, standard_scale="var", show=False)
        fig = _save_figure(out / f"marker_dotplot_{safe}.png")
    artifacts = [_artifact(f"marker gene dotplot grouped by {groupby}", fig)]
    result = {
        "message": f"Dotplot of {len(marker_genes)} marker sets across {len(categories)} groups of '{groupby}'"
        + (
            f" (dropped {len(dropped_genes)} missing genes"
            + (f", emptied {len(emptied_sets)} sets: {emptied_sets}" if emptied_sets else "")
            + ")"
            if dropped_genes
            else ""
        ),
        "reference": REFERENCE,
        "artifacts": artifacts,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "groupby": groupby,
        "n_groups": len(categories),
        "n_marker_sets": len(marker_genes),
        "dropped_genes": dropped_genes,
        "emptied_marker_sets": emptied_sets,
    }
    if cluster_labels is not None:
        # Tutorial cell 62: unmapped clusters stay NaN (not filled). JSON keys are strings, so key the
        # mapping by the actual category values (identical to the tutorial map for string categories).
        by_category = {
            c: cluster_labels[str(c)]
            for c in adata.obs[groupby].astype("category").cat.categories
            if str(c) in cluster_labels
        }
        adata.obs[annotation_key] = adata.obs[groupby].map(by_category)
        h5ad = _write(adata, out, "annotated.h5ad")
        artifacts.insert(0, _artifact(f"AnnData with cell-type labels in obs['{annotation_key}']", h5ad))
        result["annotation_key"] = annotation_key
        result["label_counts"] = {
            str(k): int(v) for k, v in adata.obs[annotation_key].value_counts(dropna=True).items()
        }
        result["unmapped_groups"] = sorted(set(categories) - set(cluster_labels))
        result["n_unlabeled_cells"] = int(adata.obs[annotation_key].isna().sum())
    return result


@clustering_mcp.tool()
def scanpy_rank_marker_genes(
    data_path: Annotated[str, "Log-normalized AnnData .h5ad with a categorical obs grouping (e.g. Leiden clusters)"],
    groupby: Annotated[str, "Existing categorical obs column whose groups are compared (e.g. 'leiden_res_0.50')"],
    n_genes: Annotated[int, "Top genes per group shown in the dotplot and returned in the summary"] = 5,
    focus_group: Annotated[str | None, "Optional group whose top genes are plotted on the UMAP (tutorial: '7')"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created"] = None,
) -> dict:
    """Rank marker genes for each group versus the rest with the Wilcoxon test and plot the top genes.
    Clustered .h5ad → .h5ad with uns['rank_genes_groups'], full ranked-gene CSV and dotplot/UMAP figures.
    """
    adata = _read(data_path)
    _require_obs(adata, [groupby], "groupby")
    categories = [str(c) for c in adata.obs[groupby].astype("category").cat.categories]
    if focus_group is not None:
        if focus_group not in categories:
            raise ValueError(f"focus_group '{focus_group}' not among obs['{groupby}'] categories {categories}")
        if "X_umap" not in adata.obsm:
            raise ValueError("focus_group plotting needs obsm['X_umap']; run scanpy_build_neighbors_umap first")
    out = _output_dir(output_dir, "rank_genes")
    artifacts = []
    with contextlib.redirect_stdout(sys.stderr):
        # Tutorial cell 67; method and ScanpyV1 preset-dependent defaults pinned explicitly
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon",
            corr_method="benjamini-hochberg",
            mean_in_log_space=True,
        )
        # Tutorial cell 69 (computes uns['dendrogram_<groupby>'] implicitly)
        sc.pl.rank_genes_groups_dotplot(adata, groupby=groupby, standard_scale="var", n_genes=n_genes, show=False)
        artifacts.append(_artifact("dotplot of top ranked genes per group", _save_figure(out / "rank_genes_groups_dotplot.png")))
        df = sc.get.rank_genes_groups_df(adata, group=None)
        csv = out / "rank_genes_groups.csv"
        df.to_csv(csv, index=False)
        artifacts.insert(0, _artifact("ranked genes for all groups (names, scores, logfoldchanges, pvals, pvals_adj)", csv))
        if focus_group is not None:
            # Tutorial cells 72-73
            top = sc.get.rank_genes_groups_df(adata, group=focus_group).head(n_genes)["names"]
            sc.pl.umap(adata, color=[*top, groupby], legend_loc="on data", frameon=False, ncols=3, show=False)
            safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in focus_group)
            artifacts.append(
                _artifact(
                    f"UMAP of top {n_genes} genes of group {focus_group}",
                    _save_figure(out / f"umap_top_genes_group_{safe}.png"),
                )
            )
    h5ad = _write(adata, out, "ranked.h5ad")
    artifacts.insert(0, _artifact("AnnData with rank_genes_groups results", h5ad))
    top_genes = {
        g: sub["names"].head(n_genes).tolist() for g, sub in df.groupby("group", observed=True, sort=False)
    }
    return {
        "message": f"Wilcoxon marker ranking for {len(categories)} groups of '{groupby}'",
        "reference": REFERENCE,
        "artifacts": artifacts,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "groupby": groupby,
        "method": "wilcoxon",
        "n_groups": len(categories),
        "n_table_rows": int(len(df)),
        "top_genes_per_group": {str(k): v for k, v in top_genes.items()},
    }
