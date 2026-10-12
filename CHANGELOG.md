# Changelog

## 0.1.0

First release.

- Recommends DADA2 `truncLenF`, `truncLenR` and `maxEE` for paired-end amplicon reads. It measures amplicon lengths from the reads, replicates DADA2's quality filter exactly for every forward/reverse split, and tests real merges with vsearch.
- `maxEE` is recommended, not given: per read, the limit that covers three quarters of the reads at the recommended cut (never below 2, with no upper limit; the report notes a limit above 5); clean reads stay at 2.
- `--its` for amplicons that vary in length (no fixed truncation length; each read is still trimmed at its first base of quality 2 or below). Its report shows where the pairs go and what the forward reads alone would keep.
- Three required arguments and two optional ones (`--its`, `--verbose`).
- One runtime dependency, numpy. matplotlib is optional (`pip install "amplitrim[plot]"`) and only draws `summary.png`: without it the report and `summary.json` are written as usual.
- A command that makes a small simulated dataset for trying the tool (`tests/simulate.py`; see the README). The evidence for each recommendation is on the [validation page](https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md) of the validation archive.
- Report (`REPORT.md`, `summary.png`, `summary.json`) with copy-paste snippets for DADA2 and QIIME 2. `summary.json` is meant for scripts and workflow managers: it has a `schema_version`, always-numeric `maxEE_f` and `maxEE_r`, and a `maxEE_rule` that says whether `maxEE` was relaxed.
