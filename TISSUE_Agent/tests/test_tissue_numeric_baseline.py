"""Numeric baseline for the two TISSUE wrappers — note 012, pinning note 005's numbers.

Upstream TISSUE has no numeric tests of its own: `repo/TISSUE/test.py` contains **zero**
assert statements. Note 005's printed numbers were therefore the only baseline in existence,
and they lived in prose. This file turns them into assertions.

It does not check "it ran". It re-runs the real pipeline on the repo's own 420 KB dataset and
compares every reported number against the value note 005 recorded, to the four decimal places
note 005 printed. Any drift in SpaGE, the conformal calibration, scikit-learn's KMeans/PCA or
numpy's RNG fails a named test.

Determinism: upstream fixes `random_seed=444` in `predict_gene_expression` and
`random_state=444` in every PCA/KMeans call, and the wrappers override none of them. So these
numbers are expected to be exactly reproducible; that expectation is what the file tests.

Running it
----------
The pinned environment has no pytest (installing one is a download, which the shop order
forbids), so the file runs standalone by default and reports the same pass/fail line:

    cd ~/Work/paper2agent
    TISSUE_Agent/tissue-env-pinned/bin/python TISSUE_Agent/tests/test_tissue_numeric_baseline.py

It is a normal pytest module as well: with pytest importable it collects as 24 tests. It uses
no pytest-only features precisely so that both paths work.

Cost: ~14 s, ~400 MB peak. Writes ~22 MB to a temporary directory and deletes it on exit; set
TISSUE_KEEP_OUTPUT=1 to keep it for inspection. Nothing is downloaded.
"""

import atexit
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TISSUE_AGENT = ROOT / "TISSUE_Agent"
DATA = TISSUE_AGENT / "repo" / "TISSUE" / "tests" / "data"

sys.path.insert(0, str(TISSUE_AGENT / "src"))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(TISSUE_AGENT / "tmp" / "mpl"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(TISSUE_AGENT / "tmp" / "cache" / "numba"))

# ---------------------------------------------------------------------------
# The baseline. Every value below is from NOTES/005-tissue-tools.md, cross-checked
# against the full-precision copy in TISSUE_Agent/out/tools_results.json.
# ---------------------------------------------------------------------------

SPATIAL = DATA / "Spatial_count.txt"
LOCATIONS = DATA / "Locations.txt"
SCRNA = DATA / "scRNA_count.txt"

TARGET_GENE = "plp1"
METHOD = "spage"
N_FOLDS = 3
N_PV = 10

PREDICT_BASELINE = {
    "n_obs": 3405,
    "n_calibration_genes": 31,
    "target_gene": "plp1",
    "method": "spage",
    "prediction_key": "spage_predicted_expression",
    "pearson_r_vs_measured": 0.3097,
    "predicted_mean": 0.5188,
    "measured_mean": 2.5771,
}

CALIBRATE_BASELINE = {
    0.23: {
        "nominal_coverage": 0.77,
        "n_calibration_genes": 31,
        "coverage_calibration_genes": 0.7883,
        "mean_interval_width_calibration_genes": 12.5879,
        "coverage_target_gene": 0.8537,
        "mean_interval_width_target_gene": 7.5545,
    },
    0.5: {
        "nominal_coverage": 0.5,
        "n_calibration_genes": 31,
        "coverage_calibration_genes": 0.5397,
        "mean_interval_width_calibration_genes": 4.1276,
        "coverage_target_gene": 0.6135,
        "mean_interval_width_target_gene": 2.6243,
    },
}

# note 005 recorded four decimal places, so that is the precision this file can honestly claim
PLACES = 4


# ---------------------------------------------------------------------------
# One pipeline run, shared by every test below.
# ---------------------------------------------------------------------------

_RESULTS = {}


def results():
    """Run predict once and calibrate at both alpha levels; cache for the whole module."""
    if _RESULTS:
        return _RESULTS

    missing = [p.name for p in (SPATIAL, LOCATIONS, SCRNA) if not p.is_file()]
    if missing:
        raise AssertionError(
            "input files missing from %s: %s — these ship with the pinned TISSUE checkout and "
            "must not be re-downloaded; restore the checkout instead" % (DATA, missing)
        )

    from tissue_tools import (  # imported late so a missing checkout fails as a clear error
        tissue_calibrate_prediction_intervals,
        tissue_predict_spatial_gene,
    )

    outdir = tempfile.mkdtemp(prefix="tissue_baseline_")
    # Each run writes ~22 MB of .h5ad. Note 010 records 11 GB lost to exactly this pattern of
    # accumulating per-call output directories, so clean up unless asked to keep them.
    if not os.environ.get("TISSUE_KEEP_OUTPUT"):
        atexit.register(shutil.rmtree, outdir, ignore_errors=True)

    predict = tissue_predict_spatial_gene(
        spatial_counts_path=str(SPATIAL),
        locations_path=str(LOCATIONS),
        scrna_counts_path=str(SCRNA),
        target_gene=TARGET_GENE,
        method=METHOD,
        n_folds=N_FOLDS,
        n_pv=N_PV,
        output_dir=outdir,
    )
    predicted_h5ad = predict["artifacts"][0]["path"]

    calibrations = {}
    for alpha in sorted(CALIBRATE_BASELINE):
        calibrations[alpha] = tissue_calibrate_prediction_intervals(
            data_path=predicted_h5ad, alpha_level=alpha, output_dir=outdir
        )

    _RESULTS.update(predict=predict, calibrations=calibrations, outdir=outdir)
    return _RESULTS


def close(actual, expected, what):
    """Compare at the precision note 005 recorded, and say what moved if it did not."""
    assert isinstance(actual, (int, float)) and not isinstance(actual, bool), (
        "%s is %r, not a number" % (what, actual)
    )
    assert not math.isnan(actual), "%s is NaN" % what
    assert round(float(actual), PLACES) == round(float(expected), PLACES), (
        "%s drifted: got %.6f, note 005 recorded %.4f (difference %.2e)"
        % (what, actual, expected, abs(actual - expected))
    )


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def test_input_files_present():
    for path in (SPATIAL, LOCATIONS, SCRNA):
        assert path.is_file(), "missing input: %s" % path


def test_upstream_still_has_no_numeric_tests():
    """The reason this file exists. If upstream gains asserts, compare against them."""
    upstream = TISSUE_AGENT / "repo" / "TISSUE" / "test.py"
    if not upstream.is_file():
        return
    assert "assert" not in upstream.read_text(), (
        "upstream test.py now contains assertions — reconcile them with this baseline"
    )


# ---------------------------------------------------------------------------
# Step 1: tissue_predict_spatial_gene
# ---------------------------------------------------------------------------


def test_predict_cell_count():
    close(results()["predict"]["n_obs"], PREDICT_BASELINE["n_obs"], "predict n_obs")


def test_predict_calibration_gene_count():
    close(
        results()["predict"]["n_calibration_genes"],
        PREDICT_BASELINE["n_calibration_genes"],
        "predict n_calibration_genes",
    )


def test_predict_records_its_settings():
    summary = results()["predict"]
    for key in ("target_gene", "method", "prediction_key"):
        assert summary[key] == PREDICT_BASELINE[key], (
            "predict %s is %r, note 005 recorded %r" % (key, summary[key], PREDICT_BASELINE[key])
        )


def test_predict_pearson_r_vs_measured():
    close(
        results()["predict"]["pearson_r_vs_measured"],
        PREDICT_BASELINE["pearson_r_vs_measured"],
        "pearson_r_vs_measured",
    )


def test_predict_predicted_mean():
    close(results()["predict"]["predicted_mean"], PREDICT_BASELINE["predicted_mean"], "predicted_mean")


def test_predict_measured_mean():
    close(results()["predict"]["measured_mean"], PREDICT_BASELINE["measured_mean"], "measured_mean")


def test_predict_scale_mismatch_is_still_present():
    """Note 005's honest caveat, pinned so a silent 'fix' upstream is noticed.

    SpaGE predictions and measured counts are on different scales and TISSUE subtracts them
    anyway (tissue/main.py:652). This asserts the gap, it does not endorse it.
    """
    summary = results()["predict"]
    assert summary["measured_mean"] > 4 * summary["predicted_mean"], (
        "the note 004/005 scale mismatch has changed: predicted_mean %.4f, measured_mean %.4f"
        % (summary["predicted_mean"], summary["measured_mean"])
    )


def test_predict_writes_a_readable_h5ad():
    artifacts = results()["predict"]["artifacts"]
    assert artifacts, "predict returned no artifacts"
    path = Path(artifacts[0]["path"])
    assert path.is_file() and path.suffix == ".h5ad", "predict artifact is not an .h5ad: %s" % path


# ---------------------------------------------------------------------------
# Step 2: tissue_calibrate_prediction_intervals, at both alpha levels
# ---------------------------------------------------------------------------


def _check_calibration(alpha, key):
    summary = results()["calibrations"][alpha]
    close(summary[key], CALIBRATE_BASELINE[alpha][key], "alpha=%s %s" % (alpha, key))


def test_calibrate_a023_nominal_coverage():
    _check_calibration(0.23, "nominal_coverage")


def test_calibrate_a023_calibration_gene_count():
    _check_calibration(0.23, "n_calibration_genes")


def test_calibrate_a023_coverage_calibration_genes():
    _check_calibration(0.23, "coverage_calibration_genes")


def test_calibrate_a023_mean_width_calibration_genes():
    _check_calibration(0.23, "mean_interval_width_calibration_genes")


def test_calibrate_a023_coverage_target_gene():
    _check_calibration(0.23, "coverage_target_gene")


def test_calibrate_a023_mean_width_target_gene():
    _check_calibration(0.23, "mean_interval_width_target_gene")


def test_calibrate_a050_nominal_coverage():
    _check_calibration(0.5, "nominal_coverage")


def test_calibrate_a050_coverage_calibration_genes():
    _check_calibration(0.5, "coverage_calibration_genes")


def test_calibrate_a050_mean_width_calibration_genes():
    _check_calibration(0.5, "mean_interval_width_calibration_genes")


def test_calibrate_a050_coverage_target_gene():
    _check_calibration(0.5, "coverage_target_gene")


def test_calibrate_a050_mean_width_target_gene():
    _check_calibration(0.5, "mean_interval_width_target_gene")


# ---------------------------------------------------------------------------
# Properties the numbers must satisfy, not just equal
# ---------------------------------------------------------------------------


def test_coverage_reaches_nominal_at_both_alphas():
    """TISSUE's actual claim: conformal intervals attain at least their nominal coverage."""
    for alpha, summary in results()["calibrations"].items():
        nominal = CALIBRATE_BASELINE[alpha]["nominal_coverage"]
        assert summary["coverage_calibration_genes"] >= nominal, (
            "alpha=%s: coverage %.4f fell below nominal %.2f"
            % (alpha, summary["coverage_calibration_genes"], nominal)
        )


def test_lower_nominal_coverage_gives_narrower_intervals():
    """The alpha_level switch must do what note 005 claims: 12.59 -> 4.13."""
    wide = results()["calibrations"][0.23]["mean_interval_width_calibration_genes"]
    narrow = results()["calibrations"][0.5]["mean_interval_width_calibration_genes"]
    assert narrow < wide, (
        "50%% intervals (%.4f) are not narrower than 77%% intervals (%.4f)" % (narrow, wide)
    )
    assert narrow < wide / 2, (
        "note 005 recorded roughly a 3x narrowing (12.5879 -> 4.1276); got %.4f -> %.4f"
        % (wide, narrow)
    )


def test_target_gene_intervals_narrower_than_calibration_genes():
    for alpha, summary in results()["calibrations"].items():
        assert (
            summary["mean_interval_width_target_gene"]
            < summary["mean_interval_width_calibration_genes"]
        ), "alpha=%s: target-gene intervals are not narrower than the calibration set" % alpha


def test_calibrate_defaults_to_the_settings_predict_recorded():
    """Step 2 with only data_path must pick up target_gene and method from uns['tissue_tool']."""
    for alpha, summary in results()["calibrations"].items():
        assert summary["target_gene"] == PREDICT_BASELINE["target_gene"], (
            "alpha=%s did not inherit target_gene" % alpha
        )
        assert summary["method"] == PREDICT_BASELINE["method"], (
            "alpha=%s did not inherit method" % alpha
        )


# ---------------------------------------------------------------------------
# The seven refusals note 005 recorded, by message
# ---------------------------------------------------------------------------


def _expect_value_error(call, fragment, case):
    from tissue_tools import (
        tissue_calibrate_prediction_intervals,
        tissue_predict_spatial_gene,
    )

    namespace = {
        "predict": tissue_predict_spatial_gene,
        "calibrate": tissue_calibrate_prediction_intervals,
    }
    try:
        call(namespace)
    except ValueError as exc:
        assert fragment in str(exc), (
            "%s raised ValueError but the message changed:\n  got:      %s\n  expected: ...%s..."
            % (case, exc, fragment)
        )
        return
    except Exception as exc:  # noqa: BLE001 - the whole point is that it must be a ValueError
        raise AssertionError(
            "%s raised %s, not ValueError: %s" % (case, type(exc).__name__, exc)
        ) from exc
    raise AssertionError("%s did not raise at all" % case)


def _predict_kwargs(**overrides):
    kwargs = dict(
        spatial_counts_path=str(SPATIAL),
        locations_path=str(LOCATIONS),
        scrna_counts_path=str(SCRNA),
        target_gene=TARGET_GENE,
        method=METHOD,
        n_folds=N_FOLDS,
    )
    kwargs.update(overrides)
    return kwargs


def test_refuses_missing_spatial_file():
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(spatial_counts_path=str(DATA / "nope.txt"))),
        "spatial_counts_path does not exist or is not a file",
        "missing spatial file",
    )


def test_refuses_unknown_gene():
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(target_gene="notagene")),
        "not found in the spatial dataset",
        "gene not in data",
    )


def test_refuses_unavailable_method():
    """Note 005 pinned 'is not available in this environment'; note 013 found that message
    was false for knn, so the wording changed. See test_knn_refusal_names_both_blockers."""
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(method="tangram")),
        "cannot be run here",
        "unsupported method",
    )


def test_refuses_unknown_method():
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(method="notamethod")),
        "is not a TISSUE prediction method",
        "unknown method",
    )


def test_knn_is_not_advertised_as_runnable():
    """The note 013 bug: SUPPORTED_METHODS claimed knn was available in this environment.

    knn is unavailable for two independent reasons, so it must not appear in the runnable
    list. If harmonypy is ever installed AND the n_neighbors passthrough is added, this test
    is the one to revisit — deliberately, not by accident.
    """
    from tissue_tools import available_methods

    assert "knn" not in available_methods(), (
        "knn is listed as runnable, but it needs harmonypy and an n_neighbors passthrough "
        "(note 013); runnable=%r" % (available_methods(),)
    )
    assert available_methods() == ["spage"], (
        "runnable methods changed: %r — reconcile with note 013" % (available_methods(),)
    )


def test_knn_refusal_names_both_blockers():
    """A refusal that hides one of two blockers sends the reader to fix the wrong thing."""
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(method="knn")),
        "n_neighbors",
        "knn refusal names the missing passthrough",
    )
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(method="knn")),
        "harmonypy",
        "knn refusal names the missing package",
    )


def test_knn_refusal_is_a_valueerror_not_a_typeerror():
    """Before note 013 this path reached upstream and died with
    `TypeError: knn_impute() missing 1 required positional argument: 'n_neighbors'`."""
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(method="knn")),
        "cannot be run here",
        "knn refused before any computation",
    )


def test_harmonypy_really_is_absent():
    """The premise of the refusal. If this fails, the knn message is now the false one."""
    import importlib.util

    assert importlib.util.find_spec("harmonypy") is None, (
        "harmonypy is now installed — rerun the knn ticket and update note 013"
    )


def test_refuses_too_many_folds():
    _expect_value_error(
        lambda ns: ns["predict"](**_predict_kwargs(n_folds=999)),
        "cannot exceed the 31 calibration genes",
        "n_folds above calibration genes",
    )


def test_refuses_non_h5ad_to_calibrate():
    _expect_value_error(
        lambda ns: ns["calibrate"](data_path=str(SPATIAL)),
        "data_path must be one of ('.h5ad',)",
        "calibrate on a non-.h5ad",
    )


def test_refuses_alpha_out_of_range():
    _expect_value_error(
        lambda ns: ns["calibrate"](
            data_path=results()["predict"]["artifacts"][0]["path"], alpha_level=1.5
        ),
        "alpha_level must be strictly between 0 and 1",
        "alpha out of range",
    )


def test_refuses_method_with_no_prediction():
    _expect_value_error(
        lambda ns: ns["calibrate"](
            data_path=results()["predict"]["artifacts"][0]["path"], method="knn"
        ),
        "run tissue_predict_spatial_gene first",
        "wrong method key",
    )


# ---------------------------------------------------------------------------
# Standalone runner, for the pinned environment that has no pytest.
# ---------------------------------------------------------------------------


def _main():
    tests = [(name, obj) for name, obj in sorted(globals().items())
             if name.startswith("test_") and callable(obj)]
    failures = []
    for name, fn in tests:
        try:
            fn()
            print("  ok   %s" % name)
        except AssertionError as exc:
            failures.append((name, exc))
            print("  FAIL %s\n       %s" % (name, exc))
        except Exception as exc:  # noqa: BLE001
            failures.append((name, exc))
            print("  ERROR %s\n       %s: %s" % (name, type(exc).__name__, exc))

    print()
    if failures:
        print("%d failed, %d passed" % (len(failures), len(tests) - len(failures)))
        return 1
    print("%d passed" % len(tests))
    return 0


if __name__ == "__main__":
    sys.exit(_main())
