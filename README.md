# AmpliTrim

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/paytonyau/amplitrim/blob/main/LICENSE)
[![Tests](https://github.com/paytonyau/amplitrim/actions/workflows/tests.yml/badge.svg)](https://github.com/paytonyau/amplitrim/actions/workflows/tests.yml)

**Choose DADA2 truncation lengths and `maxEE` from your own reads.**

AmpliTrim is for **Illumina paired-end amplicon reads** (tested on 16S and ITS data from MiSeq and NovaSeq runs). In DADA2 you pick `truncLenF`, `truncLenR` and `maxEE`. Cut too short, and your longest amplicons cannot merge. Cut too long, and noisy read ends cost you reads. Most people choose by eye from quality plots. AmpliTrim tests candidate settings on your own reads, recommends the best, and shows why. It does not change your files: it only recommends settings.

```
FASTQ  →  cutadapt  →  AmpliTrim  →  DADA2
```

Run AmpliTrim after primer removal and before DADA2 (in R or in QIIME 2). It needs paired reads (R1 and R2), not long reads or single-end data.

## Why use it

- **It uses every read pair.** After primer removal, reads often differ in length. AmpliTrim does not need them equal, so nothing is thrown away.
- **It predicts the result.** It replicates DADA2's quality filter exactly and merges a sample of your pairs with vsearch, so you see how many will pass and merge.
- **It reads your data, not a table.** It measures your amplicon lengths, checks that your longest amplicons can merge, and handles ITS and other regions where lengths vary.
- **It is ready to use.** You get copy-paste settings for DADA2, in R or in QIIME 2's `dada2 denoise-paired`, and a stable `summary.json` for pipelines.

## How it compares

We compared AmpliTrim with an existing automatic tool for choosing truncation lengths, in real DADA2 runs on public datasets and at the same `maxEE`.

| | AmpliTrim | Existing tool |
|---|---|---|
| Read pairs it can use | all | 8% to 98% (it needs reads of one length) |
| Pairs merged by DADA2 with its cut | real samples (4 datasets): ahead on 2, tied on 2 (mean +0.6 points)<br>mock communities (5 datasets): ahead on 2, tied on 3 (mean +0.6 points)<br>wastewater (1 dataset): behind by 1.0 point\* | the reference |
| Error in the predicted filter pass rate | 0.8 points at most (mean 0.2) | 0.2 to 7.7 points |
| ITS reads | supported | cannot be run |

\* Wastewater is shown on its own. Its reads are noisy: it is the only real dataset in this comparison with poor read quality, an estimated 54% of its pairs pass DADA2's default filter, and the loss is to read quality, not to short reads (only 0.6% of its pairs are shorter than the cut). Its reads were also trimmed to mixed lengths, and every pair was in reverse orientation. It is the one case where the other tool was ahead: AmpliTrim chose a cut 6 nt longer, which cost about a point.

Merging is a tie: only mock community B (+2.6) is beyond the 2-point margin we set in advance, and the two datasets used in development (mock A and soil V4) both tied. The gains are in the other rows.

<details>
<summary>Per dataset: pairs merged by DADA2, AmpliTrim's cut minus the other tool's (points)</summary>

| Type | Dataset (16S) | Samples | Role | Difference |
|---|---|---|---|---|
| Real samples | V1-V2 real samples | 33 | held out | +1.1 |
| Real samples | V4 soil | 36 | development | 0.0 |
| Real samples | V3-V4 mouse gut with spike-ins | 12 | held out | +1.3 |
| Real samples | V3-V4 surface swabs | 36 | held out | 0.0 |
| Mock communities | V1-V3 mock community A | 6 | development | 0.0 |
| Mock communities | V1-V2 inter-laboratory mock | 6 | held out | 0.0 |
| Mock communities | V1-V3 mock community B | 8 | held out | +2.6 |
| Mock communities | V4-V6 inter-laboratory mock | 5 | held out | +0.3 |
| Mock communities | V3-V4 ZymoBIOMICS mock (D6300) | 10 | held out | 0.0 |
| Wastewater | V3-V4 wastewater | 9 | held out | -1.0 |

</details>

What held, what did not, and the limits are on the [validation page](https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md), in the [archive](https://github.com/paytonyau/amplitrim-validation) that also holds the scripts and the criteria written down before each run.

## Quick start

```bash
pip install .      # Python 3.9 or later; vsearch must also be on your PATH
amplitrim --input-dir reads/ --work-dir work --out-dir report
```

Example output (from simulated reads; your numbers will differ):

```
AmpliTrim 0.1.0
Analysing reads...

Recommended DADA2 settings
  truncLenF = 231   truncLenR = 217   maxEE = 2   minOverlap = 12

  99.7% of the 12,000 read pairs analysed should pass filtering and merge.
  Covers amplicons up to 436 nt (your longest: 434 nt).

Report: report/REPORT.md
```

`report/REPORT.md` has the settings with copy-paste snippets for DADA2 in R and in QIIME 2's `dada2 denoise-paired`, what to expect, and anything worth checking. To try it without your own data, make a small simulated dataset (from a clone of this repository; the generator needs numpy):

```bash
python tests/simulate.py demo --samples 2 --pairs 1500 --seed 11
amplitrim --input-dir demo/primer_trimmed --work-dir work --out-dir out
```

These reads are simulated, so use them to see what the report looks like, not to judge the tool. Expect lengths that cover the longest amplicon (434 nt), a `maxEE` of 2, and more than 70% of the pairs passing and merging.

## What you need

- **Python 3.9 or later** and **[vsearch](https://github.com/torognes/vsearch)** on your `PATH` ([install steps](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#install), including Windows).
- **Paired FASTQ files** in one folder, named like `sample_R1.fastq.gz` and `sample_R2.fastq.gz`.
- **Primers and adapters already removed** (for example with cutadapt), and no quality trimming. R1 must be the forward read ([details](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#input)).

## Options

Required: `--input-dir`, `--work-dir` and `--out-dir`. Optional: `--its` (for amplicons that vary a lot in length, such as ITS) and `--verbose`. AmpliTrim chooses everything else, including `maxEE` ([user guide](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#arguments)).

## How the settings are chosen

| Setting | How |
|---|---|
| `truncLenF`, `truncLenR` | The shortest total that can merge your longest amplicons (the amplicon plus 12), split between the two reads where quality allows. Real merges are tested with vsearch |
| `maxEE` | 2 (DADA2's default) for clean reads. For noisy reads, a limit for each read that covers three quarters of its reads, never below 2, with a note if it goes above 5 ([how](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#how-maxee-is-chosen)) |
| `minOverlap` | 12 (DADA2's default) |

The outputs are `REPORT.md` (read this one), `summary.json` (for scripts; [every key](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#summaryjson-for-scripts-and-workflow-managers)) and, if matplotlib is installed (`pip install ".[plot]"`), `summary.png` (a figure of your amplicon lengths and how errors build up).

## Limits

- **It is for DADA2.** QIIME 2's other denoiser, Deblur, uses different settings (one trim length for merged reads), which AmpliTrim does not recommend ([guide](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#what-amplitrim-recommends)).
- **The yield is an optimistic estimate.** It is measured before DADA2 denoises the reads. Use it to compare settings, and run DADA2 once to confirm.
- **It checks quality and mergeability, not accuracy.**
- **A tail of much longer amplicons** (for example non-target 18S) can force long, noisy reads. Only the longest 0.2% are set aside, so remove such reads first ([guide](https://github.com/paytonyau/amplitrim/blob/main/docs/usage.md#limitations)).

## Citation, licence and contributing

MIT (see `LICENSE`). Please cite AmpliTrim (`CITATION.cff`) and [vsearch](https://github.com/torognes/vsearch). vsearch (GPLv3 and BSD-2-Clause) runs as a separate program, and its own licence applies. numpy and matplotlib are installed separately. For bug reports, questions and contributions, see [`CONTRIBUTING.md`](https://github.com/paytonyau/amplitrim/blob/main/CONTRIBUTING.md).
