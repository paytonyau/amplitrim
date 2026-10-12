"""AmpliTrim recommends maxEE: per read, the limit that covers three quarters of the reads at the recommended cut (int of the 75th percentile of their expected errors
plus 1), never below 2, with no upper limit; where that is not (2,2) the search is repeated at it.

vsearch is replaced by a stand-in that merges a pair when the amplicon plus the minimum overlap fits in the two truncated reads, so these tests run whether or not
vsearch is installed, and what they check is the rule and nothing else."""
import json

import numpy as np
import pytest

from amplitrim import pipeline
from amplitrim.report import screen_summary, write_outputs
from fake_vsearch import install, reads_with_a_bad_reverse_tail as _reads_with_a_bad_reverse_tail


@pytest.fixture(autouse=True)
def stand_in_for_vsearch(monkeypatch):
    install(monkeypatch, pipeline)


def test_clean_reads_keep_dadas_default_and_the_report_says_why(reads, tmp_path):
    res = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "w"))
    assert res.cfg.max_ee == 2.0 and res.relaxed is False and tuple(res.rule["max_ee"]) == (2, 2)
    assert max(res.rule["percentile"]) < 2 and res.pass_rate_default == res.recommendation.ee_pass_rate
    write_outputs(res, tmp_path / "o")
    report = (tmp_path / "o" / "REPORT.md").read_text()
    assert "maxEE is 2, DADA2's default" in report and "already covers three quarters of the reads" in report
    assert "relaxed" not in " ".join(screen_summary(res, tmp_path / "R.md"))


def test_clean_reads_give_the_same_result_as_asking_for_maxEE_2(reads, tmp_path):
    auto = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "a"))
    given = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "b", max_ee=2.0))
    assert (auto.recommendation.trim_f, auto.recommendation.trim_r) == (given.recommendation.trim_f, given.recommendation.trim_r)
    assert [(c.trim_f, c.trim_r, c.yield_est) for c in auto.candidates] == [(c.trim_f, c.trim_r, c.yield_est) for c in given.candidates]


def test_noisy_reverse_reads_get_a_looser_reverse_limit_and_a_forward_limit_that_stays_2(tmp_path):
    poor = _reads_with_a_bad_reverse_tail(tmp_path / "poor", share_good=0.4)
    res = pipeline.run(pipeline.Settings(input_dir=poor, work_dir=tmp_path / "w"))
    assert res.relaxed and res.cfg.max_ee == (2.0, 4.0)                       # forward: clean, so 2; reverse: three quarters of the reads carry 3.8 expected errors, so 4
    assert res.rule["percentile"][0] < 1 and 3.5 < res.rule["percentile"][1] < 4 and res.rule["raw"] == [1, 4] and res.rule["above_warning"] is False
    assert res.pass_rate_default == pytest.approx(0.4, abs=0.05) and res.recommendation.ee_pass_rate > 0.95
    note = " ".join(screen_summary(res, tmp_path / "REPORT.md"))
    assert "maxEE = 2,4" in note and "not 2" in note and "40%" in note
    write_outputs(res, tmp_path / "o")
    report = (tmp_path / "o" / "REPORT.md").read_text()
    assert "not DADA2's default of 2" in report and "75th percentile" in report and "never below 2" in report and "never above" not in report and "has not been tested on independent data" in report and "In a DADA2 test" not in report
    assert "never below 2). Only 40% of the read pairs pass at maxEE 2" in report                              # a sentence starts with a capital
    summary = json.loads((tmp_path / "o" / "summary.json").read_text())
    rule = summary["maxEE_rule"]
    assert rule["relaxed"] is True and rule["percentile"] == 0.75 and rule["cap"] is None and rule["warn_above"] == 5 and rule["pass_rate_at_maxEE_2"] == pytest.approx(res.pass_rate_default)
    assert rule["cut_measured_at"] == res.rule["cut"] and rule["percentile_expected_errors"] == res.rule["percentile"]
    assert (summary["recommended"]["maxEE_f"], summary["recommended"]["maxEE_r"]) == (2.0, 4.0)


def test_reads_where_nothing_passes_maxEE_2_are_measured_at_the_shortest_cut_that_covers_the_amplicons(tmp_path):
    dead = _reads_with_a_bad_reverse_tail(tmp_path / "dead", share_good=0.0)
    res = pipeline.run(pipeline.Settings(input_dir=dead, work_dir=tmp_path / "w"))
    assert res.relaxed and res.cfg.max_ee == (2.0, 4.0) and res.recommendation is not None and res.pass_rate_default is None
    assert sum(res.rule["cut"]) == res.floor                                  # no recommendation at maxEE 2, so the rule looked at the cut of the coverage floor
    assert any("none of the read pairs pass" in ln for ln in screen_summary(res, tmp_path / "R.md"))         # inside a sentence, after a colon
    write_outputs(res, tmp_path / "o")
    assert "never below 2). None of the read pairs pass at maxEE 2" in (tmp_path / "o" / "REPORT.md").read_text()


def test_there_is_no_upper_limit_and_the_report_notes_a_limit_above_5(tmp_path):
    noisy = _reads_with_a_bad_reverse_tail(tmp_path / "noisy", share_good=0.4, bad_q=12)        # three quarters of the reverse reads carry about 7 expected errors
    res = pipeline.run(pipeline.Settings(input_dir=noisy, work_dir=tmp_path / "w"))
    raw = res.rule["raw"][1]
    assert res.relaxed and raw > 5 and res.rule["max_ee"][1] == raw and res.cfg.max_ee == (2.0, float(raw))        # the limit follows the reads: nothing holds it at 5
    assert res.rule["above_warning"] is True and any("above 5: these reads are very noisy" in w for w in res.warnings)
    write_outputs(res, tmp_path / "o")
    text = (tmp_path / "o" / "REPORT.md").read_text()
    assert "above 5: these reads are very noisy" in text and "capped" not in text
    rule = json.loads((tmp_path / "o" / "summary.json").read_text())["maxEE_rule"]
    assert rule["cap"] is None and rule["warn_above"] == 5                                            # the key stays (the schema never loses one); it says there is no cap
    assert any("above 5: these reads are very noisy" in w for w in json.loads((tmp_path / "o" / "summary.json").read_text())["warnings"])


def test_a_limit_at_or_below_5_gives_no_note(tmp_path):
    res = pipeline.run(pipeline.Settings(input_dir=_reads_with_a_bad_reverse_tail(tmp_path / "mid", share_good=0.4, bad_q=22), work_dir=tmp_path / "w"))
    if res.rule is not None and not res.rule["above_warning"]:
        assert not any("very noisy" in w for w in res.warnings)


def test_reads_so_noisy_that_nothing_passes_maxEE_2_get_the_limit_they_need_and_a_note(tmp_path):
    hopeless = _reads_with_a_bad_reverse_tail(tmp_path / "hopeless", share_good=0.0, bad_q=12)   # every reverse read carries about 7 expected errors
    res = pipeline.run(pipeline.Settings(input_dir=hopeless, work_dir=tmp_path / "w"))
    assert res.relaxed and res.recommendation is not None and res.pass_rate_default is None and res.rule["max_ee"][1] > 5        # with no cap, the limit is what three quarters of the reads need
    assert any("above 5: these reads are very noisy" in w for w in res.warnings)
    write_outputs(res, tmp_path / "o")
    assert (tmp_path / "o" / "summary.json").is_file()


def test_if_the_relaxed_search_finds_nothing_the_result_stays_at_maxEE_2(tmp_path, monkeypatch):
    """The guard that remains: if the search at the rule's limit finds no recommendation (as real vsearch could do), the report keeps maxEE 2."""
    real, calls = pipeline._run_once, []

    def second_finds_nothing(s, say=lambda msg: None):
        res = real(s, say); calls.append(s.max_ee)
        if len(calls) == 2: res.recommendation = None
        return res
    monkeypatch.setattr(pipeline, "_run_once", second_finds_nothing)
    res = pipeline.run(pipeline.Settings(input_dir=_reads_with_a_bad_reverse_tail(tmp_path / "noisy", share_good=0.4), work_dir=tmp_path / "w"))
    assert len(calls) == 2 and not res.relaxed and res.cfg.max_ee == 2.0 and res.rule is not None


def test_a_maxEE_given_in_the_settings_is_used_exactly(tmp_path):
    poor = _reads_with_a_bad_reverse_tail(tmp_path / "poor", share_good=0.4)
    res = pipeline.run(pipeline.Settings(input_dir=poor, work_dir=tmp_path / "w", max_ee=2.0))
    assert not res.relaxed and res.cfg.max_ee == 2.0 and res.rule is None and res.recommendation.ee_pass_rate < 0.5


def test_reads_that_cannot_span_the_amplicons_are_not_measured_and_not_run_twice(reads, tmp_path, monkeypatch):
    """When no truncation can merge the longest amplicon, maxEE cannot help, so the rule is not computed and the search is not repeated."""
    monkeypatch.setattr(pipeline, "estimate_insert_sizes", lambda *a, **k: [700] * 1000)
    runs = []
    original = pipeline._run_once
    monkeypatch.setattr(pipeline, "_run_once", lambda s, say=None: runs.append(s.max_ee) or original(s, say or (lambda m: None)))
    res = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "w"))
    assert res.recommendation is None and not res.relaxed and res.rule is None and runs == [2.0]


def test_its_reads_are_left_at_the_default_because_there_is_no_cut(reads, tmp_path):
    res = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "w", marker="its"))
    assert res.cfg.max_ee == 2.0 and not res.relaxed and res.rule is None
    write_outputs(res, tmp_path / "o")
    text = (tmp_path / "o" / "REPORT.md").read_text()
    assert "does not adjust it for ITS" in text and "no cut to measure them at" in text and "only trimmed where its quality first falls to 2 or below" in text


def test_the_rule_arithmetic_on_known_expected_errors():
    """int(75th percentile) + 1, never below 2, with no upper limit, per read."""
    class Pairs:
        def __init__(self, ef, er): self._ef, self._er = np.array(ef), np.array(er)
        def cut_errors(self, tf, tr): return self._ef, self._er
    clean = pipeline._rule_max_ee(Pairs([0.2] * 100, [0.9] * 100), 100, 100)
    assert clean["raw"] == [1, 1] and clean["max_ee"] == [2, 2] and clean["above_warning"] is False
    mid = pipeline._rule_max_ee(Pairs([2.4] * 100, [3.9] * 100), 100, 100)
    assert mid["raw"] == [3, 4] and mid["max_ee"] == [3, 4] and mid["above_warning"] is False
    skew = pipeline._rule_max_ee(Pairs([0.1] * 80 + [9.0] * 20, [0.1] * 70 + [9.0] * 30), 100, 100)      # the 75th percentile: forward clean, reverse not
    assert skew["max_ee"] == [2, 10] and skew["raw"][1] == 10 and skew["above_warning"] is True               # no cap: the limit is what the reads need
    assert pipeline._rule_max_ee(Pairs([], []), 100, 100) is None
    assert pipeline._collapse(3, 3) == 3.0 and pipeline._collapse(2, 4) == (2.0, 4.0)


def test_cut_errors_counts_expected_errors_over_the_cut_and_only_for_reads_long_enough():
    from amplitrim.trimmer import ReadPairs
    def rec(i, q, ln): return (f"r{i}", "ACGT" * (ln // 4) + "A" * (ln % 4), "+", chr(33 + q) * ln)
    pairs = ReadPairs([rec(i, 20, 100 if i < 40 else 250) for i in range(100)], [rec(i, 30, 250) for i in range(100)])         # 40 forward reads are only 100 nt
    ef, er = pairs.cut_errors(200, 150)
    assert len(ef) == 60 and len(er) == 100                                    # the short forward reads do not reach the cut
    assert np.allclose(ef, 200 * 0.01) and np.allclose(er, 150 * 0.001)         # Q20 is 0.01 expected errors per base, Q30 is 0.001


def test_without_matplotlib_the_recommender_still_recommends_and_says_how_to_get_the_figure(reads, tmp_path, monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "matplotlib", None)                                              # a machine where `import matplotlib` fails
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    res = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "w"))
    report = write_outputs(res, tmp_path / "o")
    assert {p.name for p in (tmp_path / "o").iterdir()} == {"REPORT.md", "summary.json"}                # no figure, everything else as usual
    text = report.read_text()
    assert "![summary]" not in text and "## Recommended settings" in text and "## Limits" in text
    rec = json.loads((tmp_path / "o" / "summary.json").read_text())["recommended"]
    assert rec["truncLenF"] == res.recommendation.trim_f and rec["truncLenR"] == res.recommendation.trim_r
    assert any('pip install "amplitrim[plot]"' in line for line in screen_summary(res, report))
