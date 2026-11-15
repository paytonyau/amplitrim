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
        result, report, _summary = recommend(args.input_dir, args.work_dir, args.out_dir, its=args.its, min_overlap=args.min_overlap,
                                             say=(lambda msg: print(msg, flush=True)) if args.verbose else None)
    except (PipelineError, RuntimeError) as err:
        print(f"Error: {err}", file=sys.stderr)
        return 3
    print("\n".join(screen_summary(result, report)))
    return 0 if result.recommendation else 2


if __name__ == "__main__":
    sys.exit(main())
