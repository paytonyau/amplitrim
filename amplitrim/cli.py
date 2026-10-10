"""Command line interface."""

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .pipeline import PipelineError, Settings, run
from .report import screen_summary, write_outputs

DESCRIPTION = """Recommend DADA2 truncLenF, truncLenR and maxEE for paired-end amplicon reads by testing candidate
settings on your own reads: how many pairs pass a DADA2-style quality filter and merge, and whether the
settings can merge your longest amplicons. The reads should have adapters and primers removed and should
not be quality-trimmed; their lengths may differ."""


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
    opt.add_argument("--verbose", action="store_true", help="show progress and detailed messages")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(message)s")
    settings = Settings(input_dir=args.input_dir, work_dir=args.work_dir, marker="its" if args.its else "16s")
    print(f"AmpliTrim {__version__}\nAnalysing reads...", flush=True)
    try:
        result = run(settings, say=(lambda msg: print(msg, flush=True)) if args.verbose else (lambda msg: None))
    except (PipelineError, RuntimeError) as err:
        print(f"Error: {err}", file=sys.stderr)
        return 3
    report = write_outputs(result, args.out_dir)
    print("\n".join(screen_summary(result, report)))
    return 0 if result.recommendation else 2


if __name__ == "__main__":
    sys.exit(main())
