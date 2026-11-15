"""amplitrim trim: the sliding-window cut of the validation (window 4, mean quality 15), the folder trim and the command."""
import gzip
import json

import pytest

from amplitrim import cli
from amplitrim.trim import sliding_cut, trim_folder

q = lambda xs: "".join(chr(33 + x) for x in xs)


def test_the_cut_is_before_the_first_window_whose_mean_is_below_15():
    assert sliding_cut(q([40] * 10 + [2] * 4 + [40] * 4)) == 9                     # windows from 7 and 8 average 30.5 and 21; from 9: 11.5
    assert sliding_cut(q([30, 30, 30, 13, 30, 30, 30, 30])) == 8                    # one low base is averaged out: kept whole
    assert sliding_cut(q([30] * 5 + [15] * 4)) == 9                                 # a mean of exactly 15 is not below 15
    assert sliding_cut(q([38, 38, 37, 36, 35, 12, 34, 33, 30, 28, 14, 12, 10, 11, 9, 8])) == 10    # the worked example of the user guide
    assert sliding_cut(q([10] * 3)) == 3 and sliding_cut("") == 0                   # shorter than a window: nothing to judge


def fq(path, quals):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("".join(f"@r{i}\n{'ACGT' * (len(x) // 4)}{'A' * (len(x) % 4)}\n+\n{q(x)}\n" for i, x in enumerate(quals)))


def test_the_folder_trim_keeps_the_pairs_in_order_and_counts_what_it_cut(tmp_path):
    fq(tmp_path / "in" / "s_R1.fastq.gz", [[37] * 100, [37] * 100])
    fq(tmp_path / "in" / "s_R2.fastq.gz", [[37] * 50 + [10] * 50, [37] * 100])
    rows = trim_folder(tmp_path / "in", tmp_path / "out")
    assert [(r["file"], r["reads"], r["reads_cut"]) for r in rows] == [("s_R1.fastq.gz", 2, 0), ("s_R2.fastq.gz", 2, 1)]
    with gzip.open(tmp_path / "out" / "s_R2.fastq.gz", "rt") as fh:
        lines = fh.read().splitlines()
    assert [lines[0], len(lines[1]), lines[4], len(lines[5])] == ["@r0", 50, "@r1", 100]


def test_the_command(tmp_path, capsys):
    fq(tmp_path / "in" / "s_R1.fastq.gz", [[37] * 100])
    fq(tmp_path / "in" / "s_R2.fastq.gz", [[37] * 60 + [8] * 40])
    assert cli.main(["trim", str(tmp_path / "in"), str(tmp_path / "out")]) == 0
    out = capsys.readouterr().out
    assert "Trimmed 2 files with window 4, mean quality 15 into" in out and "s_R2.fastq.gz: 1 reads, 100% cut, mean length 100 -> 60 nt" in out and "not the tested" not in out
    meta = json.loads((tmp_path / "out" / "amplitrim_trim.json").read_text())
    assert meta["window"] == 4 and meta["quality"] == 15 and meta["files"][1]["bases_after"] == 60
    assert cli.main(["trim", str(tmp_path / "in"), str(tmp_path / "out2"), "--quality", "20"]) == 0 and "not the tested values" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        cli.main(["trim", str(tmp_path / "in"), str(tmp_path / "in")])
    (tmp_path / "empty").mkdir()
    assert cli.main(["trim", str(tmp_path / "empty"), str(tmp_path / "out3")]) == 2
