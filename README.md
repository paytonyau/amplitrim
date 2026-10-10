# AmpliTrim

**Choose DADA2 truncation lengths and maxEE from your own reads.**

`truncLenF`, `truncLenR` and `maxEE` decide how many read pairs survive filtering, how many of those merge, and whether your longest amplicons can merge at all. They are usually picked by eye from quality plots. AmpliTrim tests candidate settings on your reads and recommends them, with a short report that shows why. It works for any amplicon region whose read pairs overlap (16S, ITS, ...): amplicon lengths are measured from your reads, not assumed.

## Quick start

```bash
pip install .      # Python 3.9 or later; vsearch must also be on your PATH
amplitrim --input-dir reads/ --work-dir work --out-dir report
```

```
AmpliTrim 0.1.0
Analysing reads...

Recommended DADA2 settings
  truncLenF = 231   truncLenR = 217   maxEE = 2   minOverlap = 12

  99.7% of the 12,000 read pairs analysed should pass filtering and merge.
  Covers amplicons up to 436 nt (your longest: 434 nt).

Report: report/REPORT.md
```

Open `report/REPORT.md`: the settings with copy-paste snippets for DADA2 (R) and QIIME 2, what to expect, and anything worth checking.

## What you need

- **Python 3.9 or later** and **[vsearch](https://github.com/torognes/vsearch)** on your `PATH` ([install steps](docs/usage.md#install), including Windows).
- **Paired FASTQ files** in one folder, named like `sample_R1.fastq.gz` and `sample_R2.fastq.gz`.
- **Primers and adapters already removed** (for example with cutadapt). Do **not** quality-trim the reads, and the forward read must be R1 in every pair ([details](docs/usage.md#input)).

## Arguments

**Required**

| Argument | Meaning |
|---|---|
| `--input-dir DIR` | Folder of paired FASTQ files |
| `--work-dir DIR` | Folder for intermediate files |
| `--out-dir DIR` | Folder for the report |

**Optional**

| Argument | Meaning |
|---|---|
| `--its` | Your amplicons vary a lot in length (ITS, for example). By default AmpliTrim assumes they are about the same length (16S, 18S, most other regions) and searches for truncation lengths. A fixed truncation would discard the short ITS amplicons, and DADA2's ITS workflow does not use one, so with `--its` AmpliTrim recommends no fixed truncation length (`truncLen` 0) and evaluates exactly that. Each read is still trimmed at its first base of quality 2 or below, and pairs are kept by `minLen`, no Ns and `maxEE`. The report shows where the pairs go and what the forward reads alone would keep |
| `--verbose` | Show progress |

That is all: AmpliTrim chooses the rest, including `maxEE`.

## What it recommends

| Setting | How it is chosen |
|---|---|
| `truncLenF`, `truncLenR` | Searched on your reads: the shortest total length that can merge your longest amplicons, split between the two reads where quality allows, with the best merge yield |
| `maxEE` | Per read, the limit that covers three quarters of the reads at the recommended cut: `int` of the 75th percentile of their expected errors, plus 1; never below 2 (DADA2's default) and never above 5. Clean reads stay at 2 |
| `minOverlap` | `12` (DADA2's default) |

The `maxEE` rule is the least tested part: the percentile and the cap were set after seeing the validation data, it has not been tested on independent data, and its effect on DADA2's results has not yet been confirmed end to end ([what it does and how far to trust it](docs/usage.md#how-maxee-is-chosen)).

## Output

| File | What it is |
|---|---|
| `REPORT.md` | **Read this.** The recommended settings, what to expect, and warnings |
| `summary.json` | The same result for scripts and workflow managers (below) |
| `summary.png` | A figure: your amplicon lengths, usable pairs by total length, how errors accumulate |

Exit status: `0` a recommendation was made, `2` nothing tested could merge your longest amplicons, `3` bad input or a missing tool.

## Use in a pipeline

`summary.json` is stable (`schema_version` 1: keys may be added, never removed or renamed). The settings are in `recommended`, and `maxEE_rule` says whether `maxEE` was relaxed:

```json
{
  "schema_version": 1,
  "maxEE_rule": {"relaxed": true, "pass_rate_at_maxEE_2": 0.54, "percentile": 0.75, "cap": 5, ...},
  "recommended": {"truncLenF": 222, "truncLenR": 226, "maxEE_f": 2.0, "maxEE_r": 5.0, "minOverlap": 12, "estimated_yield": 0.95, ...},
  "warnings": []
}
```

```python
import json
rec = json.load(open("report/summary.json"))["recommended"]
print(rec["truncLenF"], rec["truncLenR"], rec["maxEE_f"], rec["maxEE_r"])    # 222 226 2.0 5.0
```

In QIIME 2's `dada2 denoise-paired` these are `--p-trunc-len-f`, `--p-trunc-len-r`, `--p-max-ee-f` and `--p-max-ee-r`. All the keys, and a Nextflow process sketch, are in the [user guide](docs/usage.md#summaryjson-for-scripts-and-workflow-managers).

## How it works

1. A shared, deterministic subsample of the reads is taken, and amplicon lengths are measured from the pairs that overlap.
2. The shortest total length (`truncLenF + truncLenR`) that can merge your longest amplicons is computed exactly: the amplicon plus the minimum overlap.
3. Totals are tested upward from there. For each total, several forward/reverse splits with a high filter pass rate get a DADA2-style filter and a vsearch merge test.
4. Among settings whose yield is within a small tolerance of the best, the one with the shortest overlap is recommended, because it relies least on the noisy read ends.

More in the [user guide](docs/usage.md#how-it-works).

## Limitations

- **The yield is optimistic: use it to compare settings, not as the retained fraction.** It is measured with vsearch before DADA2 denoises the reads, so DADA2 can keep fewer. Run DADA2 once with the recommended settings to confirm.
- **Only quality and mergeability are assessed**, not taxonomic accuracy.
- **A tail of much longer amplicons than the rest** (for example non-target 18S amplified by the same primers) can force long, noisy reads on the whole dataset: AmpliTrim sets aside only the longest 0.2% of amplicons and has no option to change that ([guide](docs/usage.md#limitations)).
- **Tested on Illumina MiSeq and NovaSeq data only.**

## Validation

Tested against real DADA2 runs on 12 public datasets (six mock communities and six real datasets: V1-V2, V1-V3, V3-V4, V4, V4-V5, V4-V6 and ITS), with the criteria written down before each run. AmpliTrim predicted the share of pairs passing DADA2's filter to within 0.8 points on all 12; its truncation merged as many pairs as FIGARO's at the same `maxEE` (equivalent, not better); the yield estimate overshoots, mostly by 7 points or less (up to 32 on a soil study). The scripts, the pre-set criteria and the write-up are in the [validation archive](https://github.com/paytonyau/amplitrim-validation).

## Citation and licence

MIT (see `LICENSE`). Please cite AmpliTrim (see `CITATION.cff`, or the repository's *Cite this repository* button) and [vsearch](https://github.com/torognes/vsearch), whose own page gives the reference to cite.

**Third-party software.** AmpliTrim runs vsearch (dual-licensed, GPLv3 and BSD-2-Clause) as a separate program. Its source code is not part of this repository, it is not distributed with it, and its own licence governs its use and distribution. The Python libraries AmpliTrim imports (numpy and matplotlib) are installed separately under their own licences.
