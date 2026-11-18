"""amplitrim trim: cut the poor ends of ITS reads with a sliding window, before DADA2.

    amplitrim trim INPUT_DIR OUTPUT_DIR [--window 4] [--quality 15]

Each read is cut before the first window of `window` bases whose mean quality is below `quality`; a read with no such window is kept whole. Every
FASTQ (.fastq, .fq, optionally .gz) in INPUT_DIR is written to OUTPUT_DIR under the same name, read for read, so the pairs stay in order and the
forward and reverse files can go to DADA2 as before (no truncation, truncQ = 2, minLen = 50; reads trimmed below minLen are dropped there).

This is the rule tested in the validation archive (register sections 33 to 44): window 4 and mean quality 15 were safe on four fungal mock
communities and kept more read pairs on every mock dataset, 5 to 57 points more (0.4.1: the report advises it for all paired-end ITS reads). Twelve trims were compared (section 44: windows 2 to 20 at mean
quality 15, and qualities 10 to 25 at window 4): this one kept the most pairs while staying safe; windows 2, 3 and 6 and qualities of 18 and above were not
safe. It is not advised for single-end (forward-only) reads, where it let false sequences through. The cut is the same code as the validation's (`sliding_cut`).
"""
import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Iterator, Tuple

WINDOW, QUALITY = 4, 15.0
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def sliding_cut(qual: str, window: int = WINDOW, q: float = QUALITY, offset: int = 33) -> int:
    """The length to keep: everything before the first window of `window` bases whose mean quality is below q (the whole read when there is none)."""
    s = [ord(c) - offset for c in qual]
    for i in range(0, max(0, len(s) - window + 1)):
        if sum(s[i:i + window]) / window < q:
            return i
    return len(s)


def fastq_records(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, p, qq = fh.readline().rstrip("\n"), fh.readline().rstrip("\n"), fh.readline().rstrip("\n")
            yield h.rstrip("\n"), s, p, qq


def trim_file(src: Path, dst: Path, window: int = WINDOW, q: float = QUALITY) -> dict:
    """Trims one FASTQ file into dst (gzip when dst ends in .gz); returns the reads, the reads cut, and the bases before and after."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")
    reads = cut = before = after = 0
    with gzip.open(part, "wt", newline="\n") if dst.name.endswith(".gz") else open(part, "w", newline="\n") as o:
        for h, s, p, qq in fastq_records(src):
            k = sliding_cut(qq, window, q)
            reads += 1
            cut += k < len(s)
            before += len(s)
            after += k
            o.write(f"{h}\n{s[:k]}\n{p}\n{qq[:k]}\n")
    part.replace(dst)
    return {"file": src.name, "reads": reads, "reads_cut": cut, "bases_before": before, "bases_after": after}


def is_fastq(p: Path) -> bool:
    return p.is_file() and p.name.lower().endswith(FASTQ_SUFFIXES)


def trim_folder(src: Path, dst: Path, window: int = WINDOW, q: float = QUALITY) -> list:
    """Every FASTQ in src trimmed into dst under the same name. Returns one summary per file."""
    files = sorted(p for p in Path(src).iterdir() if is_fastq(p))
    return [trim_file(f, Path(dst) / f.name, window, q) for f in files]


def describe(rows: list) -> str:
    """One line per file: the reads cut and the mean length before and after."""
    out = []
    for r in rows:
        n = max(1, r["reads"])
        out.append(f"  {r['file']}: {r['reads']:,} reads, {r['reads_cut'] / n:.0%} cut, mean length {r['bases_before'] / n:.0f} -> {r['bases_after'] / n:.0f} nt")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="amplitrim trim", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", type=Path)
    ap.add_argument("output_dir", type=Path)
    ap.add_argument("--window", type=int, default=WINDOW, help="bases per window (default 4, the tested value)")
    ap.add_argument("--quality", type=float, default=QUALITY, help="the mean quality below which a window cuts the read (default 15, the tested value)")
    a = ap.parse_args(argv)
    if a.window < 1 or not 0 <= a.quality <= 60:
        ap.error("--window must be at least 1 and --quality between 0 and 60")
    if not a.input_dir.is_dir():
        ap.error(f"not a folder: {a.input_dir}")
    if a.input_dir.resolve() == a.output_dir.resolve():
        ap.error("the output folder must differ from the input folder")
    rows = trim_folder(a.input_dir, a.output_dir, a.window, a.quality)
    if not rows:
        print(f"No FASTQ files in {a.input_dir}")
        return 2
    (a.output_dir / "amplitrim_trim.json").write_text(json.dumps({"window": a.window, "quality": a.quality, "files": rows}, indent=1), encoding="utf-8")
    untested = "" if (a.window, a.quality) == (WINDOW, QUALITY) else " (not the tested values 4 and 15)"
    print(f"Trimmed {len(rows)} files with window {a.window}, mean quality {a.quality:g}{untested} into {a.output_dir}:")
    print(describe(rows))
    print("Next: DADA2 on these files as before (no truncation, truncQ = 2, minLen = 50, maxEE 2,2).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
