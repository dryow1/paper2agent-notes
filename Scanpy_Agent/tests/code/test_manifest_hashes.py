"""Re-check the SHA-256 of every file listed in tests/MANIFEST.md that is present.

MANIFEST.md is the only source of truth; this file holds no hashes of its own. Files that
are absent are reported as skips, not failures: note 010 lists several of the recorded paths
as safe to delete, and a deleted file is a known state. A file that is present but whose
hash has moved is a failure, because it means the baseline every number in NOTES/ was
measured against is no longer the file on disk.

Reads only. Nothing is downloaded and nothing is written.
"""

import hashlib
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
MANIFEST = ROOT / "Scanpy_Agent" / "tests" / "MANIFEST.md"

# | `path` *(optional marker)* | `sha256` | bytes | shape | source |
ROW = re.compile(
    r"^\|\s*`(?P<path>[^`]+)`[^|]*\|\s*`(?P<sha>[0-9a-f]{64})`\s*\|"
    r"\s*(?P<size>\d+)\s*\|\s*(?P<shape>[^|]*?)\s*\|"
)


def parse_manifest():
    """Every table row of MANIFEST.md, as (path, sha256, size, shape)."""
    rows = []
    for line in MANIFEST.read_text().splitlines():
        match = ROW.match(line)
        if match:
            rows.append(
                (
                    match.group("path"),
                    match.group("sha"),
                    int(match.group("size")),
                    match.group("shape"),
                )
            )
    return rows


ROWS = parse_manifest()


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_manifest_exists():
    assert MANIFEST.is_file(), "tests/MANIFEST.md is missing: the baseline is unrecorded"


def test_manifest_is_not_empty():
    """A manifest that parses to nothing would make every other test below vacuously pass."""
    assert len(ROWS) >= 32, "MANIFEST.md parsed to %d rows, expected at least 32" % len(ROWS)


def test_manifest_paths_are_unique():
    paths = [row[0] for row in ROWS]
    duplicates = {p for p in paths if paths.count(p) > 1}
    assert not duplicates, "listed more than once in MANIFEST.md: %s" % sorted(duplicates)


def test_manifest_paths_are_relative_to_repo_root():
    offenders = [p for p, _, _, _ in ROWS if p.startswith("/") or ".." in p]
    assert not offenders, "paths must be relative to the repo root: %s" % offenders


@pytest.mark.parametrize("rel,expected_sha,expected_size,_shape", ROWS, ids=[r[0] for r in ROWS])
def test_recorded_file_still_matches(rel, expected_sha, expected_size, _shape):
    path = ROOT / rel
    if not path.exists():
        pytest.skip("not on disk: %s" % rel)

    actual_size = path.stat().st_size
    assert actual_size == expected_size, (
        "%s is %d bytes, MANIFEST.md records %d" % (rel, actual_size, expected_size)
    )

    actual_sha = sha256_of(path)
    assert actual_sha == expected_sha, (
        "%s hashes to %s, MANIFEST.md records %s — same name, different file"
        % (rel, actual_sha, expected_sha)
    )
