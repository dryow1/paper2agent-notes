"""Coverage for the INBOX door (src/inbox.py) only.

Scope is the door's decisions, not the Scanpy tools it calls - those have their own 64 tests.

Two deliberate choices:

* The oversize rule is checked by shrinking ``inbox.MAX_FILE_MB`` for one test rather than
  writing an 80 MB file. The rule under test is "st_size over the limit refuses"; the value
  of the limit is a constant, asserted separately.
* The memory look-ahead is checked against a fabricated ``info`` dict, because it is a pure
  function of shape. No giant file is created.

Fixtures are mechanism fixtures - a 10x5 matrix of zeros - never invented study data. The
"accept" path uses the smallest real teaching file already on this disk and skips if absent.
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import inbox  # noqa: E402

REFUSED, RAN, READ_ERROR = 2, 0, 3

#: Smallest real counts file on this box (pbmc68k_reduced, counts restored, 700 x 765).
#: Not in git, so tests that need it skip when it is missing.
TEACHING_FILE = PROJECT_ROOT / "tmp" / "pbmc68k" / "01_counts_restored.h5ad"


def run_door(tray: Path, *extra: str) -> int:
    """Invoke the door's CLI against one tray and return its exit code."""
    argv = ["inbox.py", "--inbox", str(tray), *extra]
    old, sys.argv = sys.argv, argv
    try:
        return inbox.main()
    finally:
        sys.argv = old


def tiny_h5ad(path: Path, *, log1p_stamp: bool = False) -> Path:
    """A 10x5 mechanism fixture. Counts unless asked to look log-transformed."""
    import anndata as ad
    import numpy as np

    a = ad.AnnData(X=np.zeros((10, 5), dtype="float32"))
    if log1p_stamp:
        a.uns["log1p"] = {"base": None}
    a.write_h5ad(path)
    return path


# ---------------------------------------------------------------------------------------
# Empty tray
# ---------------------------------------------------------------------------------------


def test_empty_tray_refuses(tmp_path, capsys):
    assert run_door(tmp_path) == REFUSED
    out = capsys.readouterr().out
    assert "REFUSED" in out
    assert "empty" in out


def test_empty_tray_downloads_nothing(tmp_path, monkeypatch, capsys):
    """The refusal must not reach the network, and must not conjure a demo dataset."""
    import urllib.request

    def explode(*a, **k):  # noqa: ANN002, ANN003
        raise AssertionError("the door attempted a network call on an empty tray")

    monkeypatch.setattr(urllib.request, "urlopen", explode)
    monkeypatch.setattr(inbox, "inspect", explode)  # nothing should be read either

    assert run_door(tmp_path) == REFUSED
    assert list(tmp_path.iterdir()) == [], "the door created files in an empty tray"


def test_readme_alone_still_counts_as_empty(tmp_path, capsys):
    """INBOX ships a README; it must not be mistaken for a dropped file."""
    (tmp_path / "README.md").write_text("instructions")
    assert run_door(tmp_path) == REFUSED
    assert "empty" in capsys.readouterr().out


def test_missing_tray_refuses(tmp_path, capsys):
    assert run_door(tmp_path / "nope") == REFUSED
    assert "no INBOX" in capsys.readouterr().out


# ---------------------------------------------------------------------------------------
# Wrong number or kind of item
# ---------------------------------------------------------------------------------------


def test_two_files_refuse(tmp_path, capsys):
    (tmp_path / "one.h5ad").touch()
    (tmp_path / "two.h5ad").touch()
    assert run_door(tmp_path) == REFUSED
    out = capsys.readouterr().out
    assert "2 items" in out and "exactly one" in out


def test_wrong_suffix_refuses(tmp_path, capsys):
    (tmp_path / "notes.txt").write_text("hello")
    assert run_door(tmp_path) == REFUSED
    assert "I accept .h5ad or .h5" in capsys.readouterr().out


def test_folder_refuses_with_conversion_hint(tmp_path, capsys):
    """A 10x mtx folder is a layout the Scanpy tools cannot read; say so, do not convert."""
    d = tmp_path / "filtered_matrix"
    d.mkdir()
    (d / "matrix.mtx").touch()
    assert run_door(tmp_path) == REFUSED
    out = capsys.readouterr().out
    assert "is a folder" in out
    assert "convert" in out.lower()


def test_unreadable_file_is_an_error_not_a_crash(tmp_path, capsys):
    """A corrupt .h5ad exits 3 with a message, not a traceback."""
    (tmp_path / "broken.h5ad").write_text("this is not HDF5")
    assert run_door(tmp_path) == READ_ERROR
    assert "ERROR reading the file" in capsys.readouterr().out


# ---------------------------------------------------------------------------------------
# Size cap
# ---------------------------------------------------------------------------------------


def test_oversize_refuses_without_building_a_blob(tmp_path, monkeypatch, capsys):
    """Shrink the limit instead of writing 80 MB; the rule is 'st_size over limit refuses'."""
    f = tiny_h5ad(tmp_path / "small.h5ad")
    monkeypatch.setattr(inbox, "MAX_FILE_MB", 0.0001)  # ~100 bytes
    assert run_door(tmp_path) == REFUSED
    out = capsys.readouterr().out
    assert "over the" in out and "door limit" in out
    assert f.exists(), "the door must not modify or remove the file it refuses"


def test_size_limit_constant_is_what_the_docs_claim():
    assert inbox.MAX_FILE_MB == 80


def test_check_size_passes_a_small_file(tmp_path):
    mb = inbox.check_size(tiny_h5ad(tmp_path / "ok.h5ad"))
    assert 0 < mb < inbox.MAX_FILE_MB


# ---------------------------------------------------------------------------------------
# Memory look-ahead (pure function of shape - no giant file needed)
# ---------------------------------------------------------------------------------------


def _info(n_obs: int, n_vars: int, **over) -> dict:
    base = {"n_obs": n_obs, "n_vars": n_vars, "counts_verdict": "raw counts",
            "counts_reason": "non-negative and integral", "has_counts_layer": False}
    base.update(over)
    return base


def test_lookahead_refuses_a_shape_too_big_for_this_box():
    """300k x 30k is only megabytes on disk when sparse, and impossible in 16 GB."""
    info = _info(300_000, 30_000)
    ram = inbox.ram_lookahead(info)
    assert not ram["within_budget"]
    with pytest.raises(inbox.Refusal) as exc:
        inbox.decide(info, ram)
    assert "over the" in str(exc.value) and "budget" in str(exc.value)


def test_lookahead_allows_a_normal_shape():
    info = _info(2_700, 13_714)  # pbmc3k
    ram = inbox.ram_lookahead(info)
    assert ram["within_budget"]
    inbox.decide(info, ram)  # must not raise


@pytest.mark.parametrize(
    ("n_obs", "n_vars", "observed_mb"),
    [(700, 765, 708), (2_730, 3_451, 954), (2_700, 13_714, 1_153), (17_041, 23_427, 6_100)],
)
def test_lookahead_tracks_the_runs_it_was_fitted_to(n_obs, n_vars, observed_mb):
    """Within 2x of every peak recorded in notes 003/006 and use-001. A guard rail, not a model."""
    est = inbox.ram_lookahead(_info(n_obs, n_vars))["estimated_peak_mb"]
    assert 0.5 * observed_mb <= est <= 2.0 * observed_mb


def test_budget_sits_below_the_cap10g_ceiling():
    assert inbox.RAM_BUDGET_MB < 10_000


# ---------------------------------------------------------------------------------------
# Raw-counts lock
# ---------------------------------------------------------------------------------------


def test_log_transformed_file_refuses(tmp_path, capsys):
    tiny_h5ad(tmp_path / "logged.h5ad", log1p_stamp=True)
    assert run_door(tmp_path) == REFUSED
    out = capsys.readouterr().out
    assert "not raw counts" in out
    assert "log1p" in out


def test_decide_points_at_the_counts_layer_when_there_is_one():
    info = _info(10, 5, counts_verdict="not counts", counts_reason="uns['log1p'] stamp present",
                 has_counts_layer=True)
    with pytest.raises(inbox.Refusal) as exc:
        inbox.decide(info, inbox.ram_lookahead(info))
    assert "layers['counts']" in str(exc.value)


# ---------------------------------------------------------------------------------------
# Inspect-only, and the accept path
# ---------------------------------------------------------------------------------------


def test_inspect_only_reports_and_stops(tmp_path, capsys):
    tiny_h5ad(tmp_path / "tiny.h5ad")
    assert run_door(tmp_path, "--inspect-only") == RAN
    out = capsys.readouterr().out
    assert "counts verdict" in out and "RAM look-ahead" in out
    assert "running scanpy_compute_qc_and_filter" not in out
    assert not (tmp_path.parent / "INBOX_OUT").exists()


@pytest.mark.skipif(not TEACHING_FILE.is_file(), reason=f"no local teaching file at {TEACHING_FILE}")
def test_accept_path_runs_qc_on_the_smallest_teaching_file(tmp_path, capsys):
    """Copied into a temp tray and removed with it; nothing is left in the real INBOX."""
    import shutil

    tray = tmp_path / "tray"
    tray.mkdir()
    shutil.copy(TEACHING_FILE, tray / "counts.h5ad")
    assert run_door(tray, "--output-dir", str(tmp_path / "out")) == RAN
    out = capsys.readouterr().out
    assert "counts verdict: raw counts" in out
    assert "QC computed" in out
    assert "Nothing was clustered, annotated or named" in out
    assert list((tmp_path / "out").glob("qc_filter_*/qc_filtered.h5ad")), "no QC artifact written"
