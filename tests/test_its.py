"""ITS mode: where the pairs go, and what the forward reads alone would keep. vsearch is replaced by a stand-in (see fake_vsearch.py), so these run without it."""
import json

import numpy as np
import pytest

from amplitrim import pipeline
from amplitrim.report import write_outputs
from amplitrim.trimmer import ReadPairs
from fake_vsearch import install, reads_with_a_bad_reverse_tail


@pytest.fixture(autouse=True)
def stand_in_for_vsearch(monkeypatch):
    install(monkeypatch, pipeline)


def run_its(folder, tmp_path, name="w"):
    res = pipeline.run(pipeline.Settings(input_dir=folder, work_dir=tmp_path / name, marker="its"))
    write_outputs(res, tmp_path / (name + "_out"))
    return res, (tmp_path / (name + "_out") / "REPORT.md").read_text(), json.loads((tmp_path / (name + "_out") / "summary.json").read_text())


def test_where_the_pairs_go_adds_up_and_agrees_with_the_filter(reads, tmp_path):
    res, report, summary = run_its(reads, tmp_path)
    x = res.its_loss
    assert x["too_short"] + x["has_n"] + x["over_maxEE"] + x["pass"] == pytest.approx(1.0)                    # every pair is in exactly one row
    assert x["pass"] == pytest.approx(res.recommendation.ee_pass_rate) and x["pass_and_merged"] == pytest.approx(res.recommendation.yield_est)
    assert x["pass_not_merged"] == pytest.approx(x["pass"] - x["pass_and_merged"]) and x["pass_not_merged"] >= 0
    assert 0 <= x["overlap_at_full_length"] <= 1 and x["pairs"] == res.sub.n_pairs
    assert summary["its"] == x and summary["schema_version"] == 1
    for need in ("## Where the pairs go", "too short after the cut at the first base of quality 2 or below", "contain an N", "pass the filter but do not merge", "do not merge even at full read length"):
        assert need in report, need


def test_the_forward_reads_alone_are_suggested_when_the_reverse_reads_are_the_problem(tmp_path):
    poor = reads_with_a_bad_reverse_tail(tmp_path / "poor", share_good=0.4)                                      # clean forward reads, a noisy reverse tail in 60% of the pairs
    res, report, summary = run_its(poor, tmp_path)
    x = res.its_loss
    assert x["pass"] == pytest.approx(0.4, abs=0.05) and x["forward_only_pass"] > 0.99
    assert x["over_maxEE"] == pytest.approx(0.6, abs=0.05) and x["too_short"] == 0 and x["has_n"] == 0
    assert "**Using the forward reads alone,**" in report and "totally reasonable" in report and "issue #236" in report and "does not measure what this does to your ASVs" in report


def test_nothing_is_suggested_when_the_forward_reads_keep_no_more_than_the_pairs(tmp_path):
    clean = reads_with_a_bad_reverse_tail(tmp_path / "clean", share_good=1.0)                                     # both reads clean: every pair passes and merges
    res, report, _ = run_its(clean, tmp_path)
    assert res.its_loss["pass_and_merged"] == pytest.approx(1.0) and res.its_loss["forward_only_pass"] == pytest.approx(1.0)
    assert "## Where the pairs go" in report and "Using the forward reads alone" not in report


def test_when_nothing_passes_the_report_still_says_where_the_pairs_went(tmp_path):
    dead = reads_with_a_bad_reverse_tail(tmp_path / "dead", share_good=0.0, bad_q=8)                              # every reverse read carries far more than 2 expected errors
    res, report, summary = run_its(dead, tmp_path)
    assert res.recommendation is None and res.its_loss["pass"] == 0 and res.its_loss["over_maxEE"] == pytest.approx(1.0)
    assert "## No recommendation" in report and "## Where the pairs go" in report and "100.0%" in report
    assert "**Using the forward reads alone,**" in report                                                         # the most useful line when nothing passes
    assert summary["recommended"] is None and summary["its"]["forward_only_pass"] > 0.99


def test_the_forward_only_share_is_never_below_the_pair_share(reads, tmp_path):
    """A pair needs both reads to pass, so the forward reads alone always keep at least as many."""
    for k, share in enumerate((0.0, 0.3, 0.7, 1.0)):
        folder = reads_with_a_bad_reverse_tail(tmp_path / f"d{k}", share_good=share)
        res, _, _ = run_its(folder, tmp_path, name=f"w{k}")
        assert res.its_loss["forward_only_pass"] >= res.its_loss["pass"] - 1e-12
    res, _, _ = run_its(reads, tmp_path, name="sim")
    assert res.its_loss["forward_only_pass"] >= res.its_loss["pass"] - 1e-12


def test_16s_reports_have_no_its_block(reads, tmp_path):
    res = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "w"))
    write_outputs(res, tmp_path / "o")
    assert "its" not in json.loads((tmp_path / "o" / "summary.json").read_text()) and "Where the pairs go" not in (tmp_path / "o" / "REPORT.md").read_text()


def test_untruncated_forward_matches_a_hand_built_reference():
    def rec(i, seq, quals): return (f"r{i}", seq, "+", "".join(chr(33 + q) for q in quals))
    fwd = [rec(0, "A" * 100, [40] * 100),                                       # kept
           rec(1, "A" * 100, [40] * 30 + [2] * 70),                             # cut at base 31: 30 nt, under minLen
           rec(2, "A" * 69 + "N" + "A" * 30, [40] * 100),                       # an N inside the read
           rec(3, "A" * 100, [15] * 100),                                       # 100 * 0.0316 = 3.2 expected errors
           rec(4, "A" * 40, [40] * 40),                                         # shorter than minLen
           rec(5, "A" * 100, [40] * 60 + [2] * 40),                             # cut at base 61: 60 nt, 0.01 expected errors: kept, although the whole read would fail
           rec(6, "A" * 79 + "N" + "A" * 20, [40] * 60 + [2] * 40)]             # the N lies after the cut at base 61, so it does not count: kept
    rev = [rec(i, "C" * 100, [40] * 100) for i in range(7)]
    pairs = ReadPairs(fwd, rev)
    assert pairs.untruncated_forward(2.0, 50).tolist() == [True, False, False, False, False, True, True]
    assert pairs.untruncated_forward(4.0, 50).tolist() == [True, False, False, True, False, True, True]                # a looser limit admits the noisy read
    assert pairs.untruncated_forward((2.0, 0.1), 50).tolist() == [True, False, False, False, False, True, True]        # the reverse limit does not matter for the forward reads
    kept_pairs = ~(lambda t: t[2] | t[3] | t[4])(pairs.untruncated(2.0, 50))
    assert not (kept_pairs & ~pairs.untruncated_forward(2.0, 50)).any()                                               # a kept pair always has a kept forward read


def test_pairs_with_an_n_are_counted_in_their_own_row(tmp_path):
    with_n = reads_with_a_bad_reverse_tail(tmp_path / "n", share_good=1.0, n_in_forward=0.2)                       # 20% of the forward reads have an N before the end of the read
    res, report, summary = run_its(with_n, tmp_path)
    x = res.its_loss
    assert x["has_n"] == pytest.approx(0.2, abs=0.03) and x["pass"] == pytest.approx(1.0 - x["has_n"], abs=0.01) and x["too_short"] == 0
    assert x["too_short"] + x["has_n"] + x["over_maxEE"] + x["pass"] == pytest.approx(1.0)
    assert f"| contain an N | {x['has_n']:.1%} |" in report and x["forward_only_pass"] == pytest.approx(x["pass"], abs=1e-9)    # the N is in the forward read, so the forward reads alone do not escape it


def test_pairs_that_pass_the_filter_but_do_not_merge_are_counted(tmp_path):
    short_reads = reads_with_a_bad_reverse_tail(tmp_path / "short", share_good=1.0, lengths=(220, 220))             # 440 nt for amplicons of 423 to 434: the longest ones cannot overlap by 12
    res, report, _ = run_its(short_reads, tmp_path)
    x = res.its_loss
    assert x["pass"] == pytest.approx(1.0, abs=0.01) and x["pass_and_merged"] == pytest.approx(0.7, abs=0.01)         # the stand-in merges the 70% of amplicons that fit
    assert x["pass_not_merged"] == pytest.approx(0.3, abs=0.01) and f"| pass the filter but do not merge | {x['pass_not_merged']:.1%} |" in report
    assert x["pass"] == pytest.approx(x["pass_and_merged"] + x["pass_not_merged"])


def test_pairs_cut_too_short_by_the_first_base_of_quality_2_are_counted(tmp_path):
    crashed = reads_with_a_bad_reverse_tail(tmp_path / "crash", share_good=1.0, crash_in_forward=0.15)             # 15% of the forward reads fall to Q2 at base 30: 30 nt, under minLen 50
    res, report, _ = run_its(crashed, tmp_path)
    x = res.its_loss
    assert x["too_short"] == pytest.approx(0.15, abs=0.03) and x["pass"] == pytest.approx(1.0 - x["too_short"], abs=0.01) and x["has_n"] == 0
    assert x["too_short"] + x["has_n"] + x["over_maxEE"] + x["pass"] == pytest.approx(1.0)
    assert f"(under {res.min_len} nt) | {x['too_short']:.1%} |" in report

