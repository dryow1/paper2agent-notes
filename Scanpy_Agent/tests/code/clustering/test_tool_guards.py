"""Verification of the input guards added in notes 007-008.

Two guards are covered:

* ``_require_raw_counts`` on ``scanpy_compute_qc_and_filter`` - refuses matrices that are
  demonstrably not raw counts (scaled, log1p'd, normalized to a target sum) while still
  accepting count matrices that happen to be float-typed or fractional.
* ``ignore_missing_genes`` on ``scanpy_plot_marker_gene_dotplot`` - lets a marker panel
  with absent genes plot in one call, reporting what it dropped.

Every fixture is derived in-process from the existing ``sub4_*`` subsamples and written to
pytest's ``tmp_path``; nothing new is stored in the repository. Fixtures are cut down to a
few hundred cells so the whole file runs in well under a minute.
"""

import json

import anndata as ad
import numpy as np
import scanpy as sc

from clustering_verify_helpers import (
    MARKERS_JSON,
    artifact,
    assert_png_nontrivial,
    call,
    call_error,
    check_result_contract,
    out_dir,
    subsample,
)

QC = "scanpy_compute_qc_and_filter"
DOTPLOT = "scanpy_plot_marker_gene_dotplot"
MARKERS = json.loads(MARKERS_JSON.read_text())

N_CELLS, N_GENES = 200, 2000


def _small_counts() -> ad.AnnData:
    """A few hundred real raw-count cells: small, but still plots like real data."""
    return ad.read_h5ad(subsample("raw"))[:N_CELLS, :N_GENES].copy()


def _small_lognorm() -> ad.AnnData:
    """The same corner of the log-normalized checkpoint, so uns['log1p'] is present."""
    return ad.read_h5ad(subsample("normalized"))[:N_CELLS, :N_GENES].copy()


def _small_clustered() -> ad.AnnData:
    """Clustered cells spread across all five res-0.02 groups, kept to the marker genes."""
    a = ad.read_h5ad(subsample("annotated"))[::20].copy()
    wanted = {g for gs in MARKERS.values() for g in gs}
    keep = [v for v in a.var_names if v in wanted] + [v for v in a.var_names if v not in wanted][:300]
    return a[:, keep].copy()


def _write(adata: ad.AnnData, tmp_path, name: str) -> str:
    path = tmp_path / f"{name}.h5ad"
    adata.write_h5ad(path)
    return str(path)


# --------------------------------------------------------------------------------------
# QC guard: inputs that must be refused
# --------------------------------------------------------------------------------------


async def test_qc_refuses_scaled_x(tmp_path):
    """Note 006's silent wrong answer: scaled X gave negative total_counts and no error."""
    a = _small_counts()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.scale(a)
    del a.uns["log1p"]  # isolate the negative-value rule from the log1p stamp
    msg = await call_error(QC, {"data_path": _write(a, tmp_path, "scaled"), "output_dir": out_dir("guard_scaled")})
    assert "negative values" in msg
    assert "scaled/z-scored" in msg


async def test_qc_refuses_log1p_stamped(tmp_path):
    """Rule 1: scanpy's own uns['log1p'] stamp."""
    msg = await call_error(
        QC, {"data_path": _write(_small_lognorm(), tmp_path, "lognorm"), "output_dir": out_dir("guard_log1p")}
    )
    assert "log1p" in msg
    assert "log-transformed" in msg


async def test_qc_refuses_log1p_without_stamp(tmp_path):
    """Rule 3: the note 007 hole - log-normalized data whose uns stamp is gone."""
    a = _small_lognorm()
    a.uns.pop("log1p", None)
    msg = await call_error(
        QC, {"data_path": _write(a, tmp_path, "lognorm_nostamp"), "output_dir": out_dir("guard_log1p_nostamp")}
    )
    assert "fractional values" in msg
    assert "log-transformed" in msg


async def test_qc_refuses_normalize_total_output(tmp_path):
    """Rule 4: normalize_total leaves every cell with the same total."""
    a = _small_counts()
    sc.pp.normalize_total(a, target_sum=1e4)
    msg = await call_error(
        QC, {"data_path": _write(a, tmp_path, "normalized_only"), "output_dir": out_dir("guard_normtotal")}
    )
    assert "near-identical per-cell totals" in msg


async def test_qc_error_hint_names_counts_layer(tmp_path):
    """The message points at layers['counts'] when the object actually has one."""
    a = _small_counts()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.scale(a)
    msg = await call_error(
        QC, {"data_path": _write(a, tmp_path, "scaled_with_counts"), "output_dir": out_dir("guard_hint")}
    )
    assert "layers['counts']" in msg


# --------------------------------------------------------------------------------------
# QC guard: inputs that must still be accepted (the false-positive watchlist)
# --------------------------------------------------------------------------------------


async def test_qc_accepts_raw_counts(tmp_path):
    """Happy path: real integer counts still pass and still filter."""
    base = out_dir("guard_raw_ok")
    res = await call(QC, {"data_path": _write(_small_counts(), tmp_path, "counts"), "output_dir": base})
    check_result_contract(res, base)
    assert res["n_obs"] > 0 and res["n_vars"] > 0
    assert res["n_obs"] <= N_CELLS and res["n_vars"] <= N_GENES
    assert_png_nontrivial(artifact(res, ".png", "qc_violin"))


async def test_qc_accepts_float_dtype_counts(tmp_path):
    """Integral values stored as float32 are counts; a dtype check would wrongly reject them."""
    a = _small_counts()
    a.X = a.X.astype(np.float32)
    base = out_dir("guard_float_ok")
    res = await call(QC, {"data_path": _write(a, tmp_path, "float_counts"), "output_dir": base})
    check_result_contract(res, base)
    assert res["n_obs"] > 0


async def test_qc_accepts_fractional_corrected_counts(tmp_path):
    """Ambient-corrected counts are fractional but still counts: large max, varying totals."""
    a = _small_counts()
    rng = np.random.default_rng(0)
    X = a.X.tocsr().astype(np.float64)
    X.data = X.data * rng.uniform(0.85, 1.0, size=X.data.shape)
    a.X = X
    assert float(X.data.max()) >= 50, "fixture must clear the rule-3 ceiling"
    base = out_dir("guard_corrected_ok")
    res = await call(QC, {"data_path": _write(a, tmp_path, "corrected"), "output_dir": base})
    check_result_contract(res, base)
    assert res["n_obs"] > 0


async def test_qc_accepts_integer_counts_below_rule3_ceiling(tmp_path):
    """Rule 3 needs fractional AND small: integral counts under the ceiling must pass."""
    a = _small_counts()
    X = a.X.tocsr().astype(np.float32)
    X.data = np.minimum(X.data, 11.0)
    a.X = X
    base = out_dir("guard_lowmax_ok")
    res = await call(QC, {"data_path": _write(a, tmp_path, "low_max_counts"), "output_dir": base})
    check_result_contract(res, base)
    assert res["n_obs"] > 0


# --------------------------------------------------------------------------------------
# Marker dotplot: ignore_missing_genes
# --------------------------------------------------------------------------------------

BOGUS = "NOT_A_GENE_XYZ"


async def test_dotplot_full_panel_unchanged_by_default(tmp_path):
    """Happy path: with every marker present, the default still plots all 15 sets."""
    base = out_dir("guard_dotplot_full")
    res = await call(
        DOTPLOT,
        {"data_path": _write(_small_clustered(), tmp_path, "clustered"), "marker_genes": MARKERS,
         "groupby": "leiden_res_0.02", "output_dir": base},
    )
    check_result_contract(res, base)
    assert res["n_marker_sets"] == 15
    assert res["dropped_genes"] == [] and res["emptied_marker_sets"] == []


async def test_dotplot_missing_gene_error_names_the_switch(tmp_path):
    """The strict default stays strict, but now tells the caller how to proceed."""
    markers = {**MARKERS, "extra": [BOGUS]}
    msg = await call_error(
        DOTPLOT,
        {"data_path": _write(_small_clustered(), tmp_path, "clustered"), "marker_genes": markers,
         "groupby": "leiden_res_0.02", "output_dir": out_dir("guard_dotplot_strict")},
    )
    assert BOGUS in msg
    assert "ignore_missing_genes=True" in msg


async def test_dotplot_ignore_missing_genes_reports_drops(tmp_path):
    """One call, unedited panel: the plot is produced and the loss is reported."""
    markers = {**MARKERS, "extra": [BOGUS]}
    base = out_dir("guard_dotplot_ignore")
    res = await call(
        DOTPLOT,
        {"data_path": _write(_small_clustered(), tmp_path, "clustered"), "marker_genes": markers,
         "groupby": "leiden_res_0.02", "ignore_missing_genes": True, "output_dir": base},
    )
    check_result_contract(res, base)
    assert res["dropped_genes"] == [BOGUS]
    assert res["emptied_marker_sets"] == ["extra"]
    assert res["n_marker_sets"] == 15  # the emptied set is gone, the other 15 remain
    assert "dropped 1 missing genes" in res["message"]
    assert_png_nontrivial(artifact(res, ".png"))


async def test_dotplot_ignore_missing_genes_matches_manual_pruning(tmp_path):
    """The switch must do exactly what notes 003/006 did by hand, not something else."""
    path = _write(_small_clustered(), tmp_path, "clustered")
    args = {"data_path": path, "groupby": "leiden_res_0.02"}
    auto = await call(
        DOTPLOT,
        {**args, "marker_genes": {**MARKERS, "extra": [BOGUS]}, "ignore_missing_genes": True,
         "output_dir": out_dir("guard_dotplot_auto")},
    )
    manual = await call(DOTPLOT, {**args, "marker_genes": MARKERS, "output_dir": out_dir("guard_dotplot_manual")})
    assert auto["n_marker_sets"] == manual["n_marker_sets"]
    assert auto["n_groups"] == manual["n_groups"]


async def test_dotplot_all_genes_missing_still_errors(tmp_path):
    """Dropping everything must fail loudly rather than plot an empty figure."""
    msg = await call_error(
        DOTPLOT,
        {"data_path": _write(_small_clustered(), tmp_path, "clustered"),
         "marker_genes": {"nothing": [BOGUS, f"{BOGUS}_2"]}, "groupby": "leiden_res_0.02",
         "ignore_missing_genes": True, "output_dir": out_dir("guard_dotplot_empty")},
    )
    assert "no marker gene is present" in msg
