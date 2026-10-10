"""A stand-in for vsearch and a generator of reads with a bad reverse tail, shared by the tests that must run without vsearch."""
import gzip
from types import SimpleNamespace

import numpy as np

INSERTS = sorted([423] * 500 + [428] * 3000 + [430] * 1000 + [434] * 500)


def _first_len(path):
    with open(path) as fh:
        fh.readline()
        return len(fh.readline().strip())


def install(monkeypatch, pipeline):
    """Merge a pair when the amplicon plus the minimum overlap fits in the two (possibly truncated) reads; amplicon lengths come from INSERTS."""
    def merge(r1, r2, out_dir, tag, min_overlap, max_diffs, kind, keep_merged=False):
        total = _first_len(r1) + _first_len(r2)
        ok = [i for i in INSERTS if i + min_overlap <= total]
        return SimpleNamespace(merge_rate=len(ok) / len(INSERTS) * (0.97 if kind == "strict" else 1.0),
                               mean_fragment_len=sum(ok) / len(ok) if ok else float("nan"), mean_overlap_diffs=0.2)
    monkeypatch.setattr(pipeline, "run_vsearch", merge)
    monkeypatch.setattr(pipeline, "estimate_insert_sizes", lambda r1, r2, out_dir, min_overlap: list(INSERTS)[:3000])
    monkeypatch.setattr(pipeline.shutil, "which", lambda name: "/usr/bin/" + name)


def reads_with_a_bad_reverse_tail(folder, share_good, bad_q=15, lengths=(231, 227), n_in_forward=0.0, crash_in_forward=0.0):
    """Forward reads are clean. The reverse read is clean for 100 bases, then Q25 (about 0.3 expected errors) in `share_good` of the pairs and Q`bad_q` in the rest.
    `lengths` are the read lengths; `n_in_forward` is the share of forward reads with an N (of high quality, so that it is not cut) at base 100;
    `crash_in_forward` is the share of forward reads whose quality falls to Q2 at base 30, so that the cut leaves them too short."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(5)
    for name in "AB":
        for mate in (1, 2):
            with gzip.open(folder / f"sample{name}_R{mate}.fastq.gz", "wt") as fh:
                for i in range(2000):
                    tail = 25 if rng.random() < share_good else bad_q
                    ln = lengths[mate - 1]
                    quals = [37] * ln if mate == 1 else [37] * 100 + [tail] * (ln - 100)
                    if mate == 1 and rng.random() < crash_in_forward:
                        quals = [37] * 30 + [2] * (ln - 30)
                    seq = list(rng.choice(list("ACGT"), len(quals)))
                    if mate == 1 and rng.random() < n_in_forward:
                        seq[99] = "N"
                    seq = "".join(seq)
                    fh.write(f"@s{name}_{i}\n{seq}\n+\n{''.join(chr(33 + q) for q in quals)}\n")
    return folder
