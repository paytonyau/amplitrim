"""0.4: every DADA2 setting can be set by the user; the default is the tested value, and a value that differs is marked in the report and summary."""
import gzip
import json

import numpy as np
import pytest

from amplitrim import cli, pipeline
from amplitrim.pipeline import PipelineError
from amplitrim.report import screen_summary, write_outputs
from amplitrim.trimmer import ReadPairs
from fake_vsearch import install, reads_with_a_bad_reverse_tail


@pytest.fixture(autouse=True)
def stand_in_for_vsearch(monkeypatch):
    install(monkeypatch, pipeline)


def run(folder, tmp_path, **kw):
    res = pipeline.run(pipeline.Settings(input_dir=folder, work_dir=tmp_path / "w", **kw))
    write_outputs(res, tmp_path / "o")
    return res, (tmp_path / "o" / "REPORT.md").read_text(), json.loads((tmp_path / "o" / "summary.json").read_text())


def test_truncq_cuts_each_read_at_its_first_base_at_or_below_the_value(tmp_path):
    q = lambda xs: "".join(chr(33 + x) for x in xs)
    for m in (1, 2):
        with gzip.open(tmp_path / f"s_R{m}.fastq.gz", "wt") as fh:
            fh.write(f"@a\n{'A' * 100}\n+\n{q([37] * 40 + [8] + [37] * 59)}\n")
    two, ten = ReadPairs.load(tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz"), ReadPairs.load(tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz", 10)
    assert two.low1[0] == 100 and ten.low1[0] == 40 and ten.trunc_q == 10


def test_defaults_given_explicitly_change_nothing(reads, tmp_path):
    a, _, sa = run(reads, tmp_path / "a")
    b, _, sb = run(reads, tmp_path / "b", trunc_q=2, max_mismatch=0, min_len=50, trim_window=4, trim_quality=15.0, trim_advice_below=1.0, min_overlap=12)
    assert b.given == {} and "settings_given" not in sb and sa["recommended"] == sb["recommended"] and sa["maxEE_rule"] == sb["maxEE_rule"]


def test_a_16s_cut_given_is_evaluated_instead_of_searched(reads, tmp_path):
    base, _, _ = run(reads, tmp_path / "a")
    tf, tr = base.recommendation.trim_f - 10, base.recommendation.trim_r
    res, report, summary = run(reads, tmp_path / "b", trunc_len=(tf, tr))
    assert (res.recommendation.trim_f, res.recommendation.trim_r) == (tf, tr) and len(res.candidates) == 1 and res.given["truncLen"] == ((tf, tr), None)
    assert "## Settings you chose" in report and f"truncLen = {tf},{tr}: given by you instead of computed by AmpliTrim." in report
    assert summary["settings_given"]["truncLen"] == {"value": [tf, tr], "tested": None} and summary["recommended"]["truncLenF"] == tf
    short, _, _ = run(reads, tmp_path / "c", trunc_len=(60, 60))
    assert any("is shorter than the" in w for w in short.warnings)
    with pytest.raises(PipelineError, match="longer than the reads"):
        run(reads, tmp_path / "d", trunc_len=(900, 100))
    with pytest.raises(PipelineError, match="not with --its"):
        run(reads, tmp_path / "e", trunc_len=(100, 100), marker="its")


def test_maxee_and_maxmismatch_given(reads, tmp_path):
    res, report, summary = run(reads, tmp_path, max_ee=(2.0, 3.0), max_mismatch=1)
    assert res.cfg.max_ee == (2.0, 3.0) and res.rule is None and res.cfg.dada2_max_mismatch == 1
    assert summary["recommended"]["maxMismatch"] == 1 and "maxMismatch = 1" in report and "maxMismatch = 1: not the default (0). In the validation: it let false sequences through on ITS." in report
    assert "maxEE = 2,3: given by you instead of computed by AmpliTrim." in report


def test_truncq_and_minlen_for_its_reach_the_snippets_and_are_marked(reads, tmp_path):
    res, report, summary = run(reads, tmp_path, marker="its", trunc_q=5, min_len=100)
    assert summary["recommended"]["truncQ"] == 5 and summary["recommended"]["minLen"] == 100
    assert "truncQ = 5, minLen = 100)" in report and "--p-trunc-q 5" in report and "truncQ = 5: not the default (2). In the validation: it changed nothing on the data tested (3 to 6)." in report
    assert summary["maxEE_rule"]["measured_on"] == "reads after the truncQ = 5 cut" and "You set: truncQ = 5; minLen = 100" in " ".join(screen_summary(res, tmp_path / "R"))


def test_the_trim_settings_change_the_advice_and_its_prediction_is_not_called_checked(tmp_path):
    poor = reads_with_a_bad_reverse_tail(tmp_path / "poor", share_good=0.4, bad_q=12)
    res, report, summary = run(poor, tmp_path / "a", marker="its", trim_window=8, trim_quality=20.0)
    t = res.trim_advice
    assert (t["window"], t["quality"], t["checked"]) == (8, 20.0, False) and "amplitrim trim INPUT_DIR TRIMMED_DIR --window 8 --quality 20" in report
    assert "this prediction is not checked against DADA2" in report and summary["trim_advice"]["checked"] is False
    none, _, _ = run(poor, tmp_path / "b", marker="its", trim_advice_below=0.3)      # 40% pass: above 0.3, so no advice
    assert none.trim_advice is None and none.given["the share below which the trim is advised"] == (0.3, 1.0)


def test_max_pairs_caps_the_pairs_analysed(reads, tmp_path):
    res, _, summary = run(reads, tmp_path, max_pairs=500)
    assert res.sub.n_pairs <= 500 and summary["settings_given"]["the read pairs analysed"]["value"] == 500


def test_values_outside_their_range_are_refused(reads, tmp_path):
    for kw in ({"trunc_q": 50}, {"max_mismatch": -1}, {"trim_window": 0}, {"trim_advice_below": 2.0}):
        with pytest.raises(PipelineError):
            run(reads, tmp_path / str(len(kw)), **kw)


def test_the_command_line(reads, tmp_path, capsys):
    rc = cli.main(["--input-dir", str(reads), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o"), "--its", "--trunc-q", "3", "--max-ee", "2,4"])
    out = capsys.readouterr().out
    summary = json.loads((tmp_path / "o" / "summary.json").read_text())
    assert rc in (0, 2) and "You set: truncQ = 3; maxEE = 2,4 (untested values are marked in the report)." in out
    assert summary["settings_given"]["truncQ"] == {"value": 3, "tested": 2} and summary["recommended"]["maxEE_r"] == 4.0
    with pytest.raises(SystemExit):
        cli.main(["--input-dir", str(reads), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o"), "--max-ee", "0"])


def test_truncq_changes_what_the_its_filter_keeps(tmp_path):
    """Half the forward reads hold one base of quality 8 at position 30: truncQ 10 cuts them to 30 bases, below minLen, so those pairs are lost;
    truncQ 2 keeps them."""
    rng = np.random.default_rng(4)
    folder = tmp_path / "q8"
    folder.mkdir()
    for m in (1, 2):
        with gzip.open(folder / f"s_R{m}.fastq.gz", "wt") as fh:
            for i in range(1000):
                q = [37] * 200
                if m == 1 and i % 2 == 0:
                    q[30] = 8
                fh.write(f"@r{i}\n{''.join(rng.choice(list('ACGT'), 200))}\n+\n{''.join(chr(33 + x) for x in q)}\n")
    two = pipeline.run(pipeline.Settings(input_dir=folder, work_dir=tmp_path / "a", marker="its"))
    ten = pipeline.run(pipeline.Settings(input_dir=folder, work_dir=tmp_path / "b", marker="its", trunc_q=10))
    assert two.its_loss["too_short"] == 0 and ten.its_loss["too_short"] == pytest.approx(0.5, abs=0.01)


def test_the_validations_findings_are_given_for_values_a_user_chooses(reads, tmp_path):
    res, report, _ = run(reads, tmp_path / "a", marker="its", trunc_q=7, min_len=150, max_mismatch=1)
    assert "truncQ = 7: not the default (2). In the validation: safe on the mock data tested, but not advised" in report and "DADA2 issue #1884" in report
    assert "minLen = 150: not the default (50). In the validation: it lost short amplicons" in report
    assert "maxMismatch = 1: not the default (0). In the validation: it let false sequences through on ITS." in report
    res, report, _ = run(reads, tmp_path / "b", marker="its", trunc_q=11)
    assert "truncQ = 11: not the default (2). In the validation: it lost true members (8 and above)." in report
    res, report, _ = run(reads, tmp_path / "c", marker="its", min_len=75)
    assert "minLen = 75: not the tested value (50); results at this value are untested." in report


def test_the_trim_prediction_is_checked_within_the_tested_settings(tmp_path):
    poor = reads_with_a_bad_reverse_tail(tmp_path / "poor", share_good=0.4, bad_q=12)
    for kw, checked in (({"trim_window": 8}, True), ({"trim_quality": 12.0}, True), ({"min_len": 100}, True), ({"trim_window": 7}, False),
                        ({"trim_window": 8, "min_len": 100}, False), ({"trunc_q": 3}, False)):
        res, _, _ = run(poor, tmp_path / str(len(str(kw)) + len(kw)) / str(abs(hash(str(kw))) % 10000), marker="its", **kw)
        assert res.trim_advice["checked"] is checked, kw
