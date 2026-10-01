# MANIFEST — the files the factory cannot rebuild

Recorded 2026-10-01. Closes open item 2 of note 010: *the reference `.h5ad` files and the
original test files are not in git; losing the working tree loses the ability to verify
anything.* This file does not save them. It makes their loss, or their silent replacement,
**detectable** — and it is small enough to live in git, which they are not.

Checked by `tests/code/test_manifest_hashes.py`, which re-hashes every file listed here that
is present and skips the ones that are not. A row that no longer matches is either a corrupted
file or a different file wearing the same name; both make every number in `NOTES/` unverifiable.

Shapes are `n_obs x n_vars` for `.h5ad` and 10x `.h5` (read from HDF5 metadata, never loaded),
`rows x columns` for `.txt`, line counts for `.py`. Sizes are exact bytes.

**32 files, 1.19 GiB. Nothing was downloaded to build this.**

---

## A. The reference chain

Every numeric test in `tests/code/clustering/` compares against these. The chain runs
top to bottom: two downloaded 10x lanes, their concatenation, then the nine checkpoints
the tutorial replay wrote. **None of it is in git** (about 1.1 GB of `.h5ad`), and none of
it can be rebuilt without the two lanes at the top.

| File | SHA-256 | Bytes | Shape | Source |
|---|---|---|---|---|
| `Scanpy_Agent/tmp/cache/scverse_tutorials/s1d1_filtered_feature_bc_matrix.h5` | `322c30a7a4905f7f113472442d4aa2c81a1ad736c86651f6c0b81e5b2ff94ac8` | 22848647 | 8785 x 36741 | figshare doi:10.6084/m9.figshare.22716739.v1, sample s1d1; md5 verified on download (002) |
| `Scanpy_Agent/tmp/cache/scverse_tutorials/s1d3_filtered_feature_bc_matrix.h5` | `ca40b287fac57ac2048f11f56f1b630c7c188718908fc217c6baf2e49bdb1982` | 20432298 | 8340 x 36741 | figshare doi:10.6084/m9.figshare.22716739.v1, sample s1d3; md5 verified on download (002) |
| `Scanpy_Agent/notebooks/clustering/data/00_raw_concat.h5ad` | `12468122dfc36aa836d3d11ba65d21d035add3df9de6a67dd0e762c40c8f16a1` | 66390263 | 17125 x 36601 | built by the tutorial replay: s1d1 + s1d3 concatenated (002) |
| `Scanpy_Agent/notebooks/clustering/ref/01_qc_filtered.h5ad` | `8eba9473b00ea17888cb0369c2e31fcfa00e704f0b599a3991bb5d8789add708` | 67065775 | 17041 x 23427 | tutorial replay step 1 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/02_scrublet.h5ad` | `5dd5aa9e6360fb28f95255adee21e580cc41309ad7f294506be8cdbc51262505` | 67382903 | 17041 x 23427 | tutorial replay step 2 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/03_normalized.h5ad` | `4f5920cd72b28805af9d709a1c0b143333f0bc6182220f7dbc3ee8dc53cf7772` | 132989791 | 17041 x 23427 | tutorial replay step 3 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/04_hvg.h5ad` | `0f386d7af10606b5ce4255f6d2cbcc7152ab2699a482148763cf4290022a7433` | 133460779 | 17041 x 23427 | tutorial replay step 4 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/05_pca.h5ad` | `467ed3b75bcfba9a5413f0b0de48a3218694c0c9abe5d2c8767ffad81811dc9c` | 137299735 | 17041 x 23427 | tutorial replay step 5 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/06_umap.h5ad` | `68259623d5e9a4e3502efd8d0464b1e105d4513836a1eed2a21ddea282123b05` | 141323317 | 17041 x 23427 | tutorial replay step 6 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/07_leiden.h5ad` | `ae624050bfbd53be9c748a3e625b7519ec0686a0e3774dfaf155f17b8508d45b` | 141411514 | 17041 x 23427 | tutorial replay step 7 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/08_annotated.h5ad` | `e5cc44f2a698e683d53c7eb9823b83f6cfe5f708521cc993c6b460a8787df84d` | 141422706 | 17041 x 23427 | tutorial replay step 8 (002) |
| `Scanpy_Agent/notebooks/clustering/ref/09_rank_genes.h5ad` | `12135886c40370ccdefc045125746eaeede6dee4aa0ea89c728d448c71ebf3d7` | 158558006 | 17041 x 23427 | tutorial replay step 9 (002) |
| `Scanpy_Agent/notebooks/clustering/data/cluster_labels_res0.02_tutorial.json` | `207181a9a6411f2eb899af4398dac66c1e740b79644e9856f128abd7f61fb8e9` | 102 | 5 keys | cluster labels the tutorial itself prints, kept for comparison (002) |
| `Scanpy_Agent/notebooks/clustering/data/marker_genes_tutorial.json` | `1c4780a7f1be130daf8db683d804c6665833024f7f50ab9cbec788c5478a2660` | 1194 | 15 keys | marker panel the tutorial itself lists, kept for comparison (002) |

---

## B. The test code that is not in git

`.gitignore` admits `test_tool_guards.py` and `test_inbox_door.py` only. The twelve files
below are the 50 original tests from note 002 plus the helpers and `conftest.py` they need.
Losing them loses the ability to verify the nine tools at all.

| File | SHA-256 | Bytes | Shape | Source |
|---|---|---|---|---|
| `Scanpy_Agent/tests/code/clustering/clustering_verify_helpers.py` | `60e1d38dcf8266d7be8b7eb9236f77624dc17c11bd0b568d662e881b48645101` | 11037 | 278 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/conftest.py` | `2a0749d712b1314a7d6b2769f583f64116cd7c312fc7f2d1c7298aa3e459fe95` | 285 | 8 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_pipeline_composition.py` | `e7e06b09db67ebbf114621688c4abe4129522007d19977e86e60e206ccb520d1` | 4360 | 110 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_build_neighbors_umap.py` | `eaeb2e7980b151a871419dd40a2591dc63936e5a4a374cf2495af9d1b8f9a348` | 2433 | 65 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_cluster_leiden.py` | `237a556db0348ee0699f4033de9deb72a4824a04cb8cec17ca53c15023adb5d9` | 5542 | 128 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_compute_qc_and_filter.py` | `b512d5c9b61876a601f8ecf97e984d58b2c7e7adbc4b9161825e646790109e6e` | 6201 | 149 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_detect_doublets_scrublet.py` | `9b0bafae3a15cefc36198e4a6434292f14d6a74bfb00dcecd759e549634ee1b3` | 3982 | 100 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_normalize_log1p.py` | `5496e892cf29e96bba2de1be4a3a3ad39dff5c8da69de8134358d48c8e805e19` | 2830 | 80 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_plot_marker_gene_dotplot.py` | `b09cd0ce3f193443467ff66c8753a31f83fc9283b018e6a4482d591cc3545579` | 6763 | 174 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_rank_marker_genes.py` | `f169e34a32fe99da736aada2e57d2b3cd57153b671591395f2912aa1df4730e3` | 5242 | 110 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_run_pca.py` | `5b2beaf7224e47c3e64d90b2f1e84918827bbdf18648add3141cc34aa241e415` | 3060 | 79 lines | written in 002 |
| `Scanpy_Agent/tests/code/clustering/test_scanpy_select_highly_variable_genes.py` | `532acc7654cff1d3d4cbde5a1a1e9cc34d6e6655ab822a2c67a7679228bb15df` | 2982 | 72 lines | written in 002 |

---

## C. Other inputs on disk

Datasets already here that the standing shop order forbids fetching again. Recorded so a
future run can tell the same file from a different one.

| File | SHA-256 | Bytes | Shape | Source |
|---|---|---|---|---|
| `data/paul15.h5` | `6161984f758dd464992edc23f1a8ab89b2081600130c6aa93a1ab12b3ada5bfb` | 10297693 | unknown | sc.datasets.paul15(), mouse haematopoietic progenitors, Paul et al. 2015 (use-001) |
| `data/pbmc3k_raw.h5ad` | `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1` | 5855727 | unreadable: ValueError | sc.datasets.pbmc3k() from https://exampledata.scverse.org/scanpy/pbmc3k_raw.h5ad (011) |
| `INBOX/AdamsonWeissman2016_GSM2406675_10X001.h5ad` | `119e3c1cf7dede4e13f887b86f9bcd797a9dc29213ee57d36aa80012d93f1c1c` | 34557246 | 5768 x 35635 | dropped in INBOX by hand: one CRISPR lane of Adamson 2016, K562 (003-adamson) |
| `TISSUE_Agent/repo/TISSUE/tests/data/Spatial_count.txt` | `39e75bf727f26662a3f25dc7d4cf52b35bc6ade73fcc33fafe5f21736df87958` | 240574 | 3405 x 33 | ships with the pinned TISSUE checkout, tests/data (004, 005) |
| `TISSUE_Agent/repo/TISSUE/tests/data/scRNA_count.txt` | `d77d0546fe15c26a5668525f4e3ebfedcfc35be7c83d24933dd4776e3d0916c3` | 96048 | 33 x 1001 | ships with the pinned TISSUE checkout, tests/data (004, 005) |
| `TISSUE_Agent/repo/TISSUE/tests/data/Locations.txt` | `3297b8e374656088f98cc991432591a74224ea7b3ff72b0999eb590f5f015b06` | 84373 | 3405 x 2 | ships with the pinned TISSUE checkout, tests/data (004, 005) |

---

## What this manifest deliberately leaves out

- **`Scanpy_Agent/tests/data/sub4_*.h5ad`** (12 files, 343 MB) — derived subsamples, rebuilt on
  demand by `subsample()`. Note 010 lists them as safe to delete. Pinning them would also pin
  the determinism of `subsample()`, which is worth doing, but it is a different ticket.
- **`Scanpy_Agent/notebooks/clustering/replay/*.h5ad`** (9 files, 1.1 GB) — the fresh-process
  replay, already verified bitwise identical to `ref/`. Note 010 lists it as safe to delete.
- **Run outputs** under `Scanpy_Agent/tmp/`, `TISSUE_Agent/out/` and `INBOX_OUT/` — every tool
  call writes a fresh timestamped directory, so these grow without bound and carry no baseline.
- **The two virtual environments** — 1.5 GB, and `tmp/build/requirements.txt` has still never
  been installed (open item 1). A hash of a venv would say nothing useful about reproducing it.

## Rebuilding a row

The top two rows of group A are downloads, and the standing shop order forbids fetching them.
If one is lost, say so in a note rather than re-downloading it without a ticket. Everything
else in group A is derived from them by the note 002 replay; group B is source, recoverable
only from a backup of the working tree.

