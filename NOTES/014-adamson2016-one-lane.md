# 014 — Adamson 2016, one CRISPR lane

**Date:** 2026-09-25 · **Machine:** Dell laptop, 16 GB, no GPU, CPU only · **No downloads, no CUDA.**
**File:** `INBOX/AdamsonWeissman2016_GSM2406675_10X001.h5ad`, 33 MB
**What it is:** one lane of Adamson 2016 — **K562 cells, a human leukaemia cell line grown in a lab**. Not patient blood, not an organoid. Every cell is the same cell type. Seven genes were switched off with CRISPR, one gene per cell, plus a control.

> **Numbering note:** this was first filed as `003`, which collided with
> `NOTES/003-scanpy-pbmc3k.md` from the earlier Phase 1 series. Renumbered to `014` — the next
> free number — so the two no longer share one. Note 003 keeps its original number and nothing
> that cites it had to change.

---

## The three answers

**1. Raw counts? — Yes.**
Whole numbers, 1 to 564, nothing negative, no sign of prior processing. The door accepted it on the first look.

**2. QC pass? — Yes, but one standard check could not be done.**
All **5,768 cells kept**, and 14,691 of 35,635 genes retained (the rest were detected in fewer than three cells). No cell was poor enough to drop — the weakest still had 520 genes detected, well above the threshold of 100.

**The mitochondrial check is dead on this file.** The 13 mitochondrial genes are listed but every one of them holds **zero counts** — they were stripped out before the file was deposited. The dataset's own `percent_mito` column says 0.00 for all 5,768 cells, and ours agrees. That check normally tells you which cells were dying; here it tells you nothing. Not an error, but a real gap in the QC.

**3. Do the groups match which gene was broken? — No. Not at all.**

This is the clear result of the run, and it is a negative one.

---

## What the grouping showed

Clustering the cells four different ways and comparing each to the perturbation labels:

| Grouping | Groups found | Agreement with labels (ARI) | Same, after shuffling the labels |
|---|---|---|---|
| coarse (res 0.1) | 2 | −0.0002 | +0.0000 |
| (res 0.5) | 6 | −0.0015 | −0.0004 |
| (res 1.0) | 8 | 0.0018 | +0.0002 |
| fine (res 2.0) | 24 | 0.0170 | −0.0002 |

Adjusted Rand: **1.0 = perfect match, 0.0 = no better than chance.** Every value is ~0. The final column is the same measurement after randomly shuffling the labels — the real scores sit essentially on top of the chance scores.

Looking at it directly is even plainer. Across the whole lane, 31 % of cells carry the control guide and each of the seven targets carries 8–12 %. Here is what each group at resolution 0.5 actually contains (% of cells):

| Group | control | BHLHE40 | CREB1 | DDIT3 | EP300 | SNAI1 | SPI1 | ZNF326 |
|---|---|---|---|---|---|---|---|---|
| 0 | 33.4 | 8.7 | 9.3 | 8.1 | 11.9 | 8.8 | 10.9 | 8.9 |
| 1 | 30.2 | 8.8 | 9.2 | 7.8 | 7.7 | 8.7 | 14.8 | 12.9 |
| 2 | 31.3 | 10.2 | 7.8 | 8.8 | 13.0 | 10.3 | 10.1 | 8.4 |
| 3 | 29.4 | 10.3 | 9.3 | 8.4 | 11.0 | 10.2 | 12.4 | 9.1 |
| 4 | 32.7 | 6.4 | 7.7 | 7.7 | 10.9 | 11.5 | 14.1 | 9.0 |
| **whole lane** | **30.8** | **9.6** | **8.8** | **8.3** | **11.0** | **9.8** | **12.1** | **9.7** |

Every row is a near-copy of the bottom row. **Each group is just a random handful of the lane.** (A sixth group of 7 cells is too small to read anything into.)

**What the groups track instead is how deeply each cell was sequenced.** Median counts per cell run from 10,872 in one group down to 2,907 in another, and detected genes from 2,726 down to 713. The clustering is largely sorting cells by depth and quality — the usual dominant signal — not by which gene was broken.

## The labels themselves are real

This matters, because "no match" could mean the labels are junk. They are not. For each of the seven targeted genes, its own expression in the cells labelled for it, against the control cells:

| Gene switched off | Level in targeted cells | Level in control cells | Effect |
|---|---|---|---|
| SPI1 | 0.012 | 0.209 | **~17× lower** |
| ZNF326 | 0.070 | 0.384 | ~5.5× lower |
| DDIT3 | 0.035 | 0.189 | ~5.5× lower |
| BHLHE40 | 0.041 | 0.221 | ~5.4× lower |
| CREB1 | 0.031 | 0.101 | ~3.2× lower |
| EP300 | 0.024 | 0.037 | ~1.6× lower |
| SNAI1 | **0.000** | 0.031 | absent in targeted cells |

Every single one is down where it should be. **The CRISPR worked, and the labels mean what they say.**

So the finding is not "the experiment failed". It is: **switching off one of these genes does not change the cell enough to make it cluster separately from cells with a different gene switched off.** These are subtle effects inside one cell line, and ordinary unsupervised grouping is too blunt to see them.

## Stop

That is the question that was asked, and it is answered. I am not going further — no differential expression per guide, no per-perturbation signature scoring, no attempt to find the effects by a method chosen after seeing that clustering failed. Choosing a new method because the first one gave a null result is how false positives are manufactured.

## Numbers for the record

- **Cells scored:** 5,752. Excluded: 6 labelled `*` (unassigned) and **10 with no label at all** — missing values that survive as blanks rather than text, easy to miss.
- **Labels:** 8 real groups — control `62(mod)_pBA581` (1,769 cells) plus SPI1 (696), EP300 (632), SNAI1 (562), ZNF326 (557), BHLHE40 (553), CREB1 (506), DDIT3 (477).
- **Genes selected for clustering:** 2,000 most variable, then 50 principal components.
- **Time and memory:** door + QC 7 s at 638 MB; clustering and comparison 42 s at 1,211 MB. Well under the 10 GB cap. No GPU.
- **Artifacts:** QC output in `INBOX_OUT/qc_filter_20260925-163609_4a7a36b6/`; clustering outputs and `summary.json` in `Scanpy_Agent/tmp/adamson/`.

## Honest limits

- **One lane, one cell line, 5,768 cells.** Nothing here generalises to the rest of Adamson 2016.
- **Absence of evidence.** ARI ~0 says the groups do not line up with the labels. It does not prove the perturbations have no transcriptional effect — a targeted test could still find one, and the knockdown table above shows there is certainly *something* to find.
- **The mitochondrial QC could not run**, so "these cells are healthy" is not something this run established.
- The comparison used adjusted Rand and mutual information from scikit-learn, plus a shuffle baseline. No new biology method was introduced.
- I did not verify this file against the original GEO record; I took it as given.
