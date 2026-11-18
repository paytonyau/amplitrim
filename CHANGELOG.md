# Changelog

## 0.4.0 (18 November 2026)

- Every setting can be set by the user: `--max-ee`, `--trunc-len` (16S: that cut is evaluated instead of searched), `--trunc-q`, `--min-len`,
  `--max-mismatch`, `--max-pairs`, and for ITS `--trim-window`, `--trim-quality` and `--trim-advice-below`. The defaults are the tested values;
  the report and `summary.json` (`settings_given`) name each value that differs, with the tested one, and the DADA2 and QIIME 2 snippets use them.
- `truncQ` is a setting of the filter simulation (16S and ITS); `summary.json` gives it under `recommended.truncQ` (a key added).
- The trim's predicted pass rate is called checked only at the tested settings (window 4, mean quality 15, truncQ 2, minLen 50).
- Run without these options, 0.4 gives exactly what 0.3 gives; on the twelve 16S datasets of the paper it recommends what 0.2.0 did.
- The validation's fourth night on ITS (register sections 43 and 44; a fourth mock community of 189 strains, Pauvert et al. 2019):
  the 0.3 rules held on four mock communities; window 4 at mean quality 15 was the best of twelve trims; `minLen` 150 lost short amplicons;
  `truncQ` 7 was safe but is not advised (DADA2 issue #1884); the trim is not advised for single-end reads; the ITS yield stays an upper
  bound. The report says these things where they apply, and when a user chooses such a value.
- The trim's predicted pass rate is called checked at every value the validation tested, one setting changed at a time.

## 0.3.0

ITS reads (`--its`). Each change was tested in the validation archive (register sections 33 to 40) on three fungal mock communities (the ATCC
mycobiome standard, the piglet isolates and mixes of Arfken et al. 2023, and the communities of Bakker 2018).

- `maxEE` is computed from the reads for ITS as well: the rule of 16S, measured on each read after DADA2's `truncQ = 2` cut (there is no common cut
  for ITS). On good reads it gives 2,2; on noisy reads that still mostly pass it relaxes maxEE, and the analysis is repeated at it.
- When fewer than half the read pairs pass DADA2's filter at maxEE 2,2, the report advises trimming the poor read ends first and predicts the effect.
  The predicted pass rate after the trim was checked against DADA2 on eleven ITS datasets (the largest difference 0.11 points) (validation archive, register sections 41 and 42).
- `amplitrim trim INPUT_DIR OUTPUT_DIR` cuts each read before the first window of 4 bases whose mean quality is below 15 (the tested values;
  `--window` and `--quality` change them, untested). It writes `amplitrim_trim.json` with what was cut.
- The ITS report has a section on what is not advised (raising `truncQ` or `maxMismatch`, choosing settings by the number of ASVs, and the
  forward-read and merge-first options for full-length ITS).
- `summary.json`: `maxEE_rule.measured_on`, `maxEE_rule.computed_max_ee` and `trim_advice` are added (no key is removed or renamed).
- 16S runs are unchanged.

## 0.2.0

- `--min-overlap` sets the minimum overlap between the forward and the reverse read, from 12 to 20 nt (default 12, DADA2's default for `mergePairs`). It sets the shortest total length (the longest amplicon plus the overlap), is used in the merge tests, and appears as `minOverlap` in the report and in `summary.json`. Use the same value in your DADA2 `mergePairs` call. The validation was run at 12.
- `amplitrim.api.recommend()` is a programmatic entry point that returns the result, the path of the report and the parsed `summary.json`, and raises `PipelineError` or `RuntimeError` for bad input or a missing vsearch. The command line now calls it, so a script, a workflow wrapper and the command line give the same answer.
- Three optional arguments (`--its`, `--min-overlap`, `--verbose`).
- A comment in `config.py` no longer claims command-line options that did not exist.
- The run time and memory of a run were measured (34 to 54 seconds and at most 292 MB on three real datasets), and the default minimum overlap of 12 nt was tested against 10, 15, 20 and 25 nt on four mock communities (validation archive, register sections 27 and 28).
- Warnings are shortened to one line each.

## 0.1.0

First release.

- Recommends DADA2 `truncLenF`, `truncLenR` and `maxEE` for paired-end amplicon reads. It measures amplicon lengths from the reads, replicates DADA2's quality filter exactly for every forward/reverse split, and tests real merges with vsearch.
- `maxEE` is recommended, not given: per read, the limit that covers three quarters of the reads at the recommended cut (never below 2, with no upper limit; the report notes a limit above 5); clean reads stay at 2.
- `--its` for amplicons that vary in length (no fixed truncation length; each read is still trimmed at its first base of quality 2 or below). Its report shows where the pairs go and what the forward reads alone would keep.
- Three required arguments and two optional ones (`--its`, `--verbose`).
- One runtime dependency, numpy. matplotlib is optional (`pip install "amplitrim[plot]"`) and only draws `summary.png`: without it the report and `summary.json` are written as usual.
- A command that makes a small simulated dataset for trying the tool (`tests/simulate.py`; see the README). The evidence for each recommendation is on the [validation page](https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md) of the validation archive.
- Report (`REPORT.md`, `summary.png`, `summary.json`) with copy-paste snippets for DADA2 and QIIME 2. `summary.json` is meant for scripts and workflow managers: it has a `schema_version`, always-numeric `maxEE_f` and `maxEE_r`, and a `maxEE_rule` that says whether `maxEE` was relaxed.
