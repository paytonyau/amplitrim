# AmpliTrim: user guide

The short version is in the [README](../README.md); this is the full guide.

## Install

Requirements: Python 3.9 or later and [vsearch](https://github.com/torognes/vsearch) on `PATH`.

**Linux and macOS**

```bash
conda create -n amplitrim python=3.12 -y && conda activate amplitrim     # or any Python 3.9+ environment
conda install -c bioconda vsearch                                         # or: apt install vsearch / brew install vsearch
pip install .
```

**Windows** (Anaconda or Miniconda Prompt). The Python packages all have Windows builds on PyPI, so use `pip` inside a conda
environment. vsearch comes from its own Windows download.

```bat
conda create -n amplitrim python=3.12 -y
conda activate amplitrim
cd C:\path\to\amplitrim
pip install -e .
```

(or `conda env create -f environment.yml`, which does the same). Then install vsearch:

1. Download `vsearch-2.32.0-win-x86_64.zip` (or the newest `win-x86_64` zip) from
   <https://github.com/torognes/vsearch/releases> and unzip it, for example to `C:\tools`. Keep its `bin` folder intact:
   `vsearch.exe` needs `zlib1.dll` and `libbz2.dll` next to it. 
2. Put that folder on the PATH, for the current prompt: `set PATH=%PATH%;C:\tools\vsearch-2.32.0-win-x86_64\bin` (to make it
   permanent, add it under Windows' Environment Variables rather than with `setx`, which can truncate a long PATH).
3. Check it with `vsearch --version`.

Then `amplitrim --input-dir C:\data\reads --work-dir C:\data\work --out-dir C:\data\report`. The work
folder is opened by the Windows build of vsearch, which reads file names in the system's ANSI code page: accented letters from
your code page work and anything else is refused with a clear message (use a plain folder such as `C:\amplitrim_work`). Quote
paths that contain spaces. WSL2 with a Linux conda is the alternative if anything here fails.

Tested on Linux and on Windows 11 (conda, Python 3.12, vsearch 2.32.0), including paths with spaces and accented letters, the error messages and exit codes. macOS is expected to work.  The GitHub Actions workflow in `.github/workflows/` runs the tests on Linux, macOS and Windows (Python 3.9 and 3.12).

A run on 12,000-30,000 read pairs takes roughly 5-10 seconds on one core, up to twice that when `maxEE` is relaxed (the search is repeated). Memory use
is about 3 KB per pair analysed.

## Input

A folder of paired FASTQ files (gzipped or not) named `sample_R1.fastq.gz` / `sample_R2.fastq.gz`, `sample_1.fastq.gz`,
SRA-style `sample_1_fastq.gz`, or Illumina-style `sample_S1_L001_R1_001.fastq.gz`.

- **Adapters and primers must already be removed** (for example with cutadapt).
- **Removing them is a separate step.** cutadapt does it (`pip install cutadapt`; it has Windows builds). With known
  primers, `-g`/`-G` remove them and anything before them, and `-a`/`-A` with the reverse-complemented primers
  remove read-through and the adapter after it:
  `cutadapt --times 2 -g FWD -G REV -a RC_REV -A RC_FWD --discard-untrimmed -m 50 -o out_R1.fastq.gz -p out_R2.fastq.gz in_R1.fastq.gz in_R2.fastq.gz`
  (`--times 2` matters: by default cutadapt removes only one adapter per read, so a read with a primer and read-through would keep one).
  With adapters only, use `-a ADAPTER -A ADAPTER` (TruSeq: `AGATCGGAAGAGC`; Nextera: `CTGTCTCTTATACACATCT`; check
  your library prep). 
- **Do not quality-trim the reads.** AmpliTrim measures how quality falls along the reads and which
  amplicon lengths overlap, and DADA2 does its own quality control with `truncLen`, `maxEE` and `truncQ`.
  Trimming first would hide exactly what AmpliTrim needs to see.
- **The forward read must be R1 in every pair.** If your pairs were sequenced in both orientations, put them the
  same way round first (swap R1 and R2 in the reversed pairs).
  Otherwise R1 and R2 swap roles in some pairs and the truncation lengths no longer mean what they say.
- Read lengths may differ, for example after adapter read-through trimming. DADA2 discards reads shorter
  than `truncLen`, and the report says how many pairs that costs.

Run AmpliTrim after adapter and primer removal and before `filterAndTrim`.

## Arguments

Three are required (`--input-dir`, `--work-dir`, `--out-dir`) and two are optional (`--its`, `--verbose`); see the [README](../README.md#options).

`--work-dir` holds intermediate files (the subsample, trimmed and merged reads) and can be deleted afterwards. It is required, and it has no default inside the output folder, for a reason that applies to Windows. The Windows build of vsearch opens files by their names in the system's ANSI code page. Accented letters from your code page work, and any other character is refused with a clear message. Use a plain folder such as `C:\amplitrim_work`.

**`--its`.** By default AmpliTrim treats the amplicons as roughly the same length (16S, 18S and most other regions): it measures their lengths and searches for truncation lengths. Where amplicon lengths vary a lot (ITS), a fixed `truncLen` would discard every read shorter than it, so DADA2's own ITS workflow does not use one. `--its` evaluates exactly that workflow. Each read is cut at its first base of quality 2 or below (`truncQ = 2`, DADA2's permissive default, so only the very worst tails are cut). A pair is then kept if both reads are at least 50 nt long (`minLen`), contain no N and are within `maxEE`. A read with a noisy tail is not shortened: it is dropped when its expected errors exceed `maxEE`. AmpliTrim reports the expected share of pairs that pass and merge, and recommends no fixed truncation length. The ITS report also shows **where the pairs go** (too short after the cut, containing an N, over `maxEE`, passing but not merging). It also says what share of the reads would pass the same filter if you used the forward reads alone. DADA2's author calls this "totally reasonable" when the reverse reads are poor (issue #236, 2017). On the one ITS dataset we tested, the filter rows and the forward-only share match DADA2's own filter exactly, and the merging rows are an optimistic upper bound. The figures, and where DADA2 loses the rest of the pairs, are on the [validation page](https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md#its---its).

## What AmpliTrim recommends

| Setting | How it is chosen |
|---|---|
| `truncLenF`, `truncLenR` | Searched on your reads (see [How it works](#how-it-works)) |
| `maxEE` | Per read, the limit that covers three quarters of the reads at the recommended cut, never below 2 and never above 5 ([how maxEE is chosen](#how-maxee-is-chosen)) |
| `minOverlap`, `maxMismatch` | `12` and `0`, DADA2's defaults for `mergePairs` |

These are DADA2 settings, for DADA2 in R or in QIIME 2's DADA2 plugin. QIIME 2 also offers Deblur, which is a different method: it takes single or already merged reads, trims them all to one length (`--p-trim-length`), and has no `truncLenF`, `truncLenR` or `maxEE`, because its quality filtering is a separate step. AmpliTrim's recommendations therefore do not apply to Deblur, and we have not tested AmpliTrim with it. The amplicon lengths in the report may help you choose a trim length, but AmpliTrim does not recommend one. See the QIIME 2 documentation for Deblur.

## Output

| File | Contents |
|---|---|
| `REPORT.md` | The recommended settings, what to expect, warnings, what other `maxEE` values would keep, and how the settings were chosen |
| `summary.json` | The same result for programs (below) |
| `summary.png` | Your amplicon lengths, the share of usable read pairs by total trimmed length, and how expected errors accumulate along the reads, with the `maxEE` line and the recommended cuts |

Exit status: `0` recommendation made, `2` nothing tested could merge your longest amplicons (the files are still written, and `recommended` is `null`), `3` bad input or missing tool (nothing is written; the message is on stderr).

### `summary.json`, for scripts and workflow managers

`schema_version` is 1: keys may be added within a version, never removed or renamed.

| Key | Meaning |
|---|---|
| `schema_version` | `1` |
| `maxEE_rule.relaxed` | `true` if the rule gave a limit above 2 and `maxEE` was relaxed from 2 |
| `maxEE_rule.pass_rate_at_maxEE_2` | The share of pairs that passed at `maxEE` 2 with the best truncation at 2 (`null` if none did) |
| `maxEE_rule.percentile`, `.cap` | The percentile and the cap of the rule (0.75 and 5) |
| `maxEE_rule.cut_measured_at`, `.percentile_expected_errors` | The cut at which the reads were measured, and the 75th percentile of the expected errors of the forward and the reverse reads there (`null` for ITS) |
| `recommended` | `null` if nothing could merge your longest amplicons; otherwise the settings below |
| `recommended.truncLenF`, `.truncLenR` | DADA2 `truncLen` (0 and 0 with `--its`: no fixed truncation length) |
| `recommended.maxEE_f`, `.maxEE_r` | The `maxEE` of the recommendation, always one number each |
| `recommended.maxEE` | The same as one number when both reads share it, a pair otherwise (kept for older scripts; prefer `maxEE_f` and `maxEE_r`) |
| `recommended.minOverlap`, `.maxMismatch` | The `mergePairs` settings assumed (12 and 0) |
| `recommended.estimated_yield` | Expected share of pairs that pass and merge (an optimistic estimate; see Limitations) |
| `recommended.longest_mergeable_insert` | The longest amplicon these settings can merge, in nt |
| `recommended.minLen` | Only with `--its` |
| `its` | Only with `--its`: where the pairs go, as shares of the pairs analysed (`too_short`, `has_n`, `over_maxEE`, `pass`, `pass_not_merged`, `pass_and_merged`) and `forward_only_pass`, the share of forward reads that would pass the same filter |
| `longest_amplicon`, `pairs_analysed`, `data`, `warnings` | What was measured, and anything worth checking |
| `maxEE_table` | What other `maxEE` values would keep, at the recommended truncation |

In QIIME 2's `dada2 denoise-paired`, `truncLenF`, `truncLenR`, `maxEE_f` and `maxEE_r` are `--p-trunc-len-f`, `--p-trunc-len-r`, `--p-max-ee-f` and `--p-max-ee-r`. Reading the settings in Python:

```python
import json
rec = json.load(open("report/summary.json"))["recommended"]
print(rec["truncLenF"], rec["truncLenR"], rec["maxEE_f"], rec["maxEE_r"])
```

A Nextflow process, as a sketch (it has not been run in a real pipeline; adapt it). Exit status 2 is a failure to Nextflow by default; use `errorStrategy` if the pipeline should continue without a recommendation.

```groovy
process AMPLITRIM {
    input:
    path reads                      // a folder of paired FASTQ files, primers and adapters already removed

    output:
    path "report/summary.json", emit: summary
    path "report/REPORT.md",    emit: report

    script:
    """
    amplitrim --input-dir ${reads} --work-dir work --out-dir report
    """
}

workflow {
    AMPLITRIM(Channel.fromPath(params.reads, type: 'dir'))
    settings = AMPLITRIM.out.summary.map { new groovy.json.JsonSlurper().parseText(it.text).recommended }
    // settings.truncLenF, settings.truncLenR, settings.maxEE_f, settings.maxEE_r go to the DADA2 step
}
```

## How it works

1. A shared, deterministic subsample of the reads is taken.
2. Amplicon lengths are measured from the pairs that overlap at full read length.
3. The shortest total length (`truncLenF + truncLenR`) that can merge your longest amplicons is computed
   exactly: the amplicon plus the minimum overlap. Shorter totals are never tested.
4. Totals are tested upward from there. For each total, the forward/reverse split matters, so several splits
   with a high quality-filter pass rate are tested. That pass rate is computed exactly for every split at
   once. Each candidate gets a DADA2-style filter (`maxEE`, `truncQ`, `maxN`) and a vsearch merge test. The
   search stops once a longer total no longer improves the yield.
5. Among the settings whose yield is within the tolerance of the best, the one with the shortest overlap
   is recommended.

## When it matters

What matters is the overlap left after truncation, not the region. With a short overlap (for example V3-V4
amplicons on 2x250 reads, 12-20 nt) the choice decides the outcome: in simulated reads vsearch merged 39%
of pairs at a 12 nt overlap, 93% at 16 nt and 99.5% at 20 nt. With 130 nt or more of overlap (V3-V4 on
2x300 reads, V1-V2 or V4 on 2x250) essentially every pair merges at any sensible truncation, and there is
little to decide. Amplicons longer than the reads can span cannot be merged by any setting; AmpliTrim then
reports that no setting can cover them, and warns when many of your pairs look like this.

## Choosing the cut-points

Two things pull against each other. Species differ in the length of their 16S region, and the longest regions need the longest reads
to merge at all: two reads truncated to a total of T nt can merge an amplicon only if it is at most T - 12 nt (DADA2 needs a 12 nt
overlap), so the settings must be long enough for the long ones. But quality falls toward the read ends (the reverse read usually first), and
the more of that tail you keep, the more reads fail the `maxEE` filter and the more mismatches sit in the overlap.

What AmpliTrim does: it takes the longest amplicon from your own reads (setting aside the longest 0.2%), which gives the
shortest total that can cover it; it counts exactly, for every forward/reverse split of the total, how many pairs pass your `maxEE`, so each
read is cut where its own quality allows; and among settings with similar yield it prefers the smallest overlap, which relies least on the noisy ends.
The third panel of `summary.png` shows the quality side: how expected errors accumulate along each read, the `maxEE` line, and where the cuts fall.

Leaving amplicons out loses those taxa for good (in a test on a mock community no taxon longer than T - 12 nt was ever found), so the choice is between long-16S taxa and cleaner reads. AmpliTrim sets aside only the longest 0.2% of amplicons and has no option to change that. A small share of much longer amplicons (for example eukaryotic 18S amplified by the same primers) can therefore force long, noisy reads on the whole dataset. On one validation mock community, leaving out the longest 1% would have raised the merged share from 57% to 90%. If that is your case, remove the off-target reads before running AmpliTrim.

## How maxEE is chosen

AmpliTrim recommends `maxEE` as well as the truncation, in two steps. It first finds the best truncation at `maxEE` 2 (DADA2's default, which DADA2's author calls "a good general purpose value"). Then, at that cut, it measures the expected errors of every read. For each read direction it sets the limit that covers three quarters of its reads. This is `int` of the 75th percentile of their expected errors, plus 1, never below 2 and never above 5. If that is (2,2), nothing changes; otherwise the search is repeated at it, because a looser limit can let a different cut win. ITS is not adjusted. Reads are not cut at a common position, so there is no cut to measure at (each read is only trimmed where its quality first falls to 2 or below). The rule has also not been tested on reads of different lengths. It is unlikely to matter much there: on the ITS dataset, relaxing `maxEE` passed 12 points more pairs but merged only 2.7 points more. If nothing passes at `maxEE` 2, the reads are measured at the shortest cut that covers your amplicons. The report shows the share passing at other values at the recommended truncation, and what a looser limit costs in the expected errors of the kept reads.

The 75th percentile and the limit of 5 are AmpliTrim's own choices: DADA2's author says that relaxing `maxEE` somewhat is also OK, but gives no number.

**How much to trust it.** The 75th percentile and the cap were chosen after we had seen the validation data, and the rule has not been tested on independent data. In a pre-registered test it did **not** meet its adoption criteria. On the five datasets where it relaxes `maxEE`, it merged more pairs than the default on all five, and about as many as relaxing to `2,5` on average, but fewer than `2,5` on three of the five. We keep it as a trade-off that we state openly: it adapts to each read, and it keeps its limits within the range that DADA2's author calls acceptable. The report's table shows what a looser limit costs in noisier reads. AmpliTrim does not measure accuracy, so if accuracy matters for your study, check a looser `maxEE` on a mock community. The numbers, and what DADA2's author has said about `maxEE`, are on the [validation page](https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md#maxee).

## Limitations

- **The yield is optimistic: use it to compare settings, not as the retained fraction.** It is measured on the raw reads with vsearch, before DADA2 denoises them. DADA2 requires zero mismatches in the overlap and cannot confirm a sequence it sees about once, so diverse samples keep fewer reads. On a real 36-sample soil study the estimate was 98% for the recommended settings while DADA2 merged 64%. What the estimate predicted well was the cost of changing the settings. Run DADA2 once with the recommended settings to confirm them.
- **The yield counts chimeras, and pairs that merge at the wrong offset, as usable.** In simulation, merges at the wrong length were 0.00% at the 12 nt minimum overlap, except when a long perfect repeat lay across the overlap. In that case a longer overlap helped. On real soil V4 reads they were 0.00% at every cut tested, against a reference of full-length merges.
- **Only read quality and mergeability are assessed, not taxonomic accuracy.**
- **A tail of much longer amplicons** can force long reads on the whole dataset (see [Choosing the cut-points](#choosing-the-cut-points)).
- **Reads of mixed lengths** cost pairs: DADA2 discards reads shorter than `truncLen`, and the report says how many. On two validation datasets with trimmed reads, no fixed truncation merged more pairs than the recommended cut.
- **It needs overlap and clean input:** primers and adapters removed, the forward read as R1, no quality trimming (see Input).
- **Validated on Illumina MiSeq and NovaSeq data only;** other platforms are untested.

## Development

```bash
pip install -e ".[test]"
pytest tests/
```

The tests need vsearch on `PATH`, except `tests/test_maxee.py` (the `maxEE` rule), which replaces vsearch with a stand-in.
