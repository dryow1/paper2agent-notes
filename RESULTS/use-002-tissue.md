# use-002 — TISSUE's 2 tools on the built-in spatial example

**Date:** 2026-09-25 · **Machine:** Dell laptop, no GPU, CPU only · **Not a Phase 1 note** — a use of the finished wrappers.
**Outcome:** **both tools ran, 17 s, peak 310 MB. The headline number, 0.7883, matches notes 004 and 005 exactly.** No new environment, **nothing downloaded** — same 411 KB file as before.

## What we asked the tools

The dataset is the one that ships inside the TISSUE repository: **411 KB of text**, already on disk from note 004.

| File | Size | What it holds |
|---|---|---|
| `Spatial_count.txt` | 235 KB | gene counts for 3,405 cells whose physical positions are known |
| `Locations.txt` | 82 KB | the x/y position of each of those cells |
| `scRNA_count.txt` | 94 KB | 1,000 reference cells measured in far more detail, but with no positions |

The problem this addresses: the spatial method measures **where** cells are but reads only a handful of genes. The reference method reads **many** genes but loses the positions. TISSUE bridges the two.

**Question 1 — to `tissue_predict_spatial_gene`:** *"Here is a gene, `plp1`. Pretend you never measured it in the spatial cells. Using the reference data, guess what its level would be in each of the 3,405 cells."* We hid the real answer so it could be scored.

**Question 2 — to `tissue_calibrate_prediction_intervals`:** *"For each guess, give a range you are N% sure the true value falls inside."* We asked at three confidence levels: **77 %** (the tutorial's default), **90 %**, and **50 %**.

That second question is the whole point of TISSUE. Plenty of methods will guess a number; this one is supposed to tell you **how much to trust each guess**.

## The certainty number

The test that matters: when the tool says "77 % sure", is it right 77 % of the time?

| We asked for | It actually got right (31 calibration genes) | Error | Typical range width |
|---|---|---|---|
| **77 %** | **78.8 %** | +1.8 pts | 12.6 |
| **90 %** | **90.8 %** | +0.8 pts | 27.5 |
| **50 %** | **54.0 %** | +4.0 pts | 4.1 |

**All three land within about 4 points of what was promised, and each is slightly *more* cautious than advertised — it errs toward wider ranges, not narrower.** The widths also move the right way: demanding more confidence forces wider ranges (4.1 → 12.6 → 27.5), which is the honest trade-off. A method that gave you 90 % confidence at the same width as 50 % would be lying.

For the one gene actually held out, `plp1`, coverage was 85.4 % / 93.4 % / 61.4 % at the three levels — consistently above target. That is a single gene rather than an average over 31, so wobble is expected.

**The raw guesses themselves are mediocre.** Predicted vs measured `plp1` correlates at **r = 0.31**. The value of this run is not that the prediction is good — it is that **the tool's stated uncertainty is honest about the prediction being mediocre.**

## Does it match notes 004–005?

**Yes, exactly.** Every figure that both runs produced is identical to the last decimal:

| Number | Note 004 | Note 005 | use-002 |
|---|---|---|---|
| Coverage, 77 %, calibration genes | 0.7883 | 0.7883 | **0.7883** |
| Coverage, 77 %, held-out `plp1` | 0.8537 | 0.8537 | **0.8537** |
| Mean width, 77 % | 12.5879 | 12.5879 | **12.5879** |
| Coverage, 50 % | — | 0.5397 | **0.5397** |
| Correlation r | 0.3097 | 0.3097 | **0.3097** |

Three runs, on three days, through two different code paths — note 004 called the TISSUE library directly, notes 005 and use-002 went through the wrappers — and the numbers do not move. The 90 % level is new here and was not run before.

## Time and RAM

| Step | Time | Peak RAM |
|---|---|---|
| `tissue_predict_spatial_gene` | 6.8 s | — |
| `tissue_calibrate_prediction_intervals` × 3 levels | 3.2 s each | — |
| **Whole run** | **17 s** | **310 MB** (`cap10g`), 411 MB in-process |

0 OOM kills, about 3 % of the 10 GB cap. This would run on a phone-class machine. Temporary outputs (32 MB) were deleted afterwards; nothing was left behind.

## What a non-scientist should not over-read

1. **"78.8 % vs 77 % — so it works" is too strong a conclusion from one small file.** This is 3,405 cells, 32 genes, one prediction method, one held-out gene. It shows the machinery behaves as advertised *here*. It is not independent validation of the paper.
2. **A wide range is not a bad answer — it is the honest one.** At 90 % confidence the typical range spans 27.5 units for a gene whose measured values mostly sit near zero. That looks useless, and in a sense it is: it is the tool correctly reporting that it cannot pin this gene down. Do not read wide ranges as a malfunction.
3. **The predictions and the measurements are not on the same scale.** Predicted `plp1` averages 0.52; measured averages 2.58 and reaches 87. TISSUE subtracts one from the other anyway, which is what makes the ranges so wide. This was flagged in note 004 and is not fixed — it is inherited from the method, not introduced by our wrappers.
4. **None of this says `plp1` does anything.** No biology was tested. There is no treatment, no disease, no comparison between conditions. The A/B split used in note 004 was an arbitrary "first half vs second half of the cells", so its p-value was a plumbing check, never a result.
5. **"Calibrated" means the ranges are the right width on average, not that any single cell's guess is right.** A cell can sit well inside its range and still have a badly wrong predicted value.
6. **This only runs on a pinned 2022 software stack.** On current libraries TISSUE fails outright (note 004). Reproducing these numbers requires that exact environment — which is itself a caveat about how durable the result is.
