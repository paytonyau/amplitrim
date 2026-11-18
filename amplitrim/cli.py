"""Command line interface."""

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .api import recommend
from .pipeline import PipelineError
from .report import screen_summary

DESCRIPTION = """Recommend DADA2 truncLenF, truncLenR and maxEE for paired-end amplicon reads by testing candidate
settings on your own reads: how many pairs pass a DADA2-style quality filter and merge, and whether the
settings can merge your longest amplicons. The reads should have adapters and primers removed and should
not be quality-trimmed; their lengths may differ."""


def _overlap(text: str) -> int:
    """argparse type: an integer from 12 to 20."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is not a whole number")
    if not 12 <= value <= 20:
        raise argparse.ArgumentTypeError("the minimum overlap must be between 12 and 20")
    return value


def _max_ee(text: str):
    """--max-ee 2 or --max-ee 2,4 (forward, reverse)."""
    parts = [float(x) for x in text.split(",")]
    if len(parts) not in (1, 2) or min(parts) <= 0:
        raise argparse.ArgumentTypeError("one number, or forward,reverse, each above 0")
    return parts[0] if len(parts) == 1 else tuple(parts)


def _two_lengths(text: str):
    """--trunc-len 240,200 (forward, reverse)."""
    try:
        f, r = (int(x) for x in text.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError("two whole numbers, forward,reverse") from None
    if min(f, r) < 1:
        raise argparse.ArgumentTypeError("both lengths at least 1")
    return (f, r)


# 0.4: the settings a user can set, with the tested default named in the help; None = not given (AmpliTrim's own)
SETTINGS = [("--max-ee", "max_ee", _max_ee, "F[,R]", "use this maxEE instead of the one AmpliTrim computes"),
            ("--trunc-len", "trunc_len", _two_lengths, "F,R", "16S: evaluate this cut instead of searching (not with --its)"),
            ("--trunc-q", "trunc_q", int, "Q", "DADA2's truncQ (default 2; other values untested)"),
            ("--min-len", "min_len", int, "N", "DADA2's minLen, with --its (default 50; other values untested)"),
            ("--max-mismatch", "max_mismatch", int, "N", "DADA2's mergePairs maxMismatch (default 0; others raised false sequences on ITS)"),
            ("--max-pairs", "max_pairs", int, "N", "read pairs analysed, at most (default 50,000; 0 = all)"),
            ("--trim-window", "trim_window", int, "N", "with --its: the window of the trim AmpliTrim advises and predicts (default 4)"),
            ("--trim-quality", "trim_quality", float, "Q", "with --its: the trim's mean quality (default 15)"),
            ("--trim-advice-below", "trim_advice_below", float, "SHARE", "with --its: advise the trim below this share passing at maxEE 2,2 (default 1: always; 0.5 is the 0.3 rule)")]


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="amplitrim", description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="example:  amplitrim --input-dir reads --work-dir work --out-dir report")
    ap.add_argument("--version", action="version", version=f"amplitrim {__version__}")
    req = ap.add_argument_group("required arguments")
    req.add_argument("--input-dir", required=True, type=Path,
                     help="folder of paired FASTQs named sample_R1/_R2.fastq.gz or "
                          "sample_S1_L001_R1_001.fastq.gz")
    req.add_argument("--work-dir", required=True, type=Path, help="folder for intermediate files")
    req.add_argument("--out-dir", required=True, type=Path, help="folder for the report")
    opt = ap.add_argument_group("optional arguments")
    opt.add_argument("--its", action="store_true",
                     help="the amplicons vary a lot in length (ITS, for example): recommend no fixed truncation length, as DADA2's "
                          "ITS workflow does (each read is still trimmed at its first base of quality 2 or below), instead of truncation lengths")
    opt.add_argument("--min-overlap", type=_overlap, default=12, metavar="N",
                     help="minimum overlap in nt between the forward and reverse reads, 12 to 20 (default 12, DADA2's mergePairs default); "
                          "use the same value as the minOverlap of your DADA2 mergePairs call")
    opt.add_argument("--verbose", action="store_true", help="show progress and detailed messages")
    dada2 = ap.add_argument_group("settings you can set yourself (0.4; the report marks a value other than the tested one)")
    for flag, dest, kind, meta, text in SETTINGS:
        dada2.add_argument(flag, dest=dest, type=kind, default=None, metavar=meta, help=text)
    return ap


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "trim":                                            # amplitrim trim INPUT_DIR OUTPUT_DIR (0.3)
        from .trim import main as trim_main
        return trim_main(argv[1:])
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(message)s")
    print(f"AmpliTrim {__version__}\nAnalysing reads...", flush=True)
    try:
        given = {dest: getattr(args, dest) for _, dest, _, _, _ in SETTINGS if getattr(args, dest) is not None}
        result, report, _summary = recommend(args.input_dir, args.work_dir, args.out_dir, its=args.its, min_overlap=args.min_overlap,
                                             say=(lambda msg: print(msg, flush=True)) if args.verbose else None, **given)
    except (PipelineError, RuntimeError) as err:
        print(f"Error: {err}", file=sys.stderr)
        return 3
    print("\n".join(screen_summary(result, report)))
    return 0 if result.recommendation else 2


if __name__ == "__main__":
    sys.exit(main())
