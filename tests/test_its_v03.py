"""ITS in 0.3: maxEE computed from the reads after DADA2's truncQ = 2 cut, and (0.4.1: for all paired-end reads; 0.3: only when fewer than half the pairs pass at
maxEE 2,2 (validation archive, register sections 33 to 40)."""
import gzip
import json

import numpy as np
import pytest

from amplitrim import pipeline
from amplitrim.report import screen_summary, write_outputs
from amplitrim.trimmer import ReadPairs
from fake_vsearch import install, reads_with_a_bad_reverse_tail


@pytest.fixture(autouse=True)
def stand_in_for_vsearch(monkeypatch):
    install(monkeypatch, pipeline)


def run_its(folder, tmp_path, **kw):
    res = pipeline.run(pipeline.Settings(input_dir=folder, work_dir=tmp_path / "w", marker="its", **kw))
    write_outputs(res, tmp_path / "o")
    return res, (tmp_path / "o" / "REPORT.md").read_text(), json.loads((tmp_path / "o" / "summary.json").read_text())


def test_clean_its_reads_get_the_rule_measured_and_it_gives_2(reads, tmp_path):
    res, report, summary = run_its(reads, tmp_path, trim_advice_below=0.5)                 # the 0.3 rule: no trim advice where half or more pass
    assert res.rule is not None and res.rule["max_ee"] == [2, 2] and res.rule["cut"] is None and not res.relaxed and res.trim_advice is None
    assert res.cfg.max_ee == 2.0 and "after DADA2's truncQ = 2 cut" in report and "both under 2" in report and "does not adjust it for ITS" not in report
    assert "## Not advised for ITS" in report and "Raising `truncQ` above 2" in report and "DADA2 issue #1884" in report and "Trimming single-end" in report and "## Trim the poor read ends first" not in report
    assert summary["maxEE_rule"]["measured_on"] == "reads after the truncQ = 2 cut" and summary["maxEE_rule"]["computed_max_ee"] == [2, 2] and "trim_advice" not in summary


def test_0_4_1_advises_the_trim_for_reads_that_mostly_pass_too(reads, tmp_path):
    res, report, summary = run_its(reads, tmp_path)                                         # about 86% pass at maxEE 2,2: the 0.3 rule gave no advice here
    t = res.trim_advice
    assert t is not None and t["below"] == 1.0 and t["pass_rate_default"] > 0.5 and t["rule_after"] == [2, 2] and res.cfg.max_ee == 2.0
    assert "## Trim the poor read ends first" in report and "AmpliTrim advises cutting the poor read ends of paired-end ITS reads first" in report
    assert "Only **" not in report and "optional here" not in report and "after the trim the rule gives 2" in report          # 86% is below the 90% wording line
    assert "kept more pairs on every mock dataset (5 to 57 points" in report and summary["trim_advice"]["rule_after"] == [2, 2]
    assert "trim the poor read ends first, then use these settings" in "\n".join(screen_summary(res, tmp_path / "REPORT.md"))


def test_reads_that_already_pass_well_get_the_note_and_the_trim_is_called_optional(tmp_path):
    good = reads_with_a_bad_reverse_tail(tmp_path / "good", share_good=0.95, bad_q=12)     # about 95% pass
    res, report, summary = run_its(good, tmp_path)
    p = res.trim_advice["pass_rate_default"]
    assert p >= 0.9 and "Your reads already pass well" in report and "the trim is optional here" in report
    assert "5.1 to 6.2 points more pairs" in report and "between -1.2 and +2.2 points" in report and f"at most the {1 - p:.0%} of pairs that fail" in report
    assert "(optional: your reads already pass well)" in "\n".join(screen_summary(res, tmp_path / "REPORT.md"))


def test_the_rule_after_the_trim_is_reported_when_it_is_not_2_2(reads, tmp_path):
    res, _, _ = run_its(reads, tmp_path)
    res.trim_advice["rule_after"] = [2, 4]
    write_outputs(res, tmp_path / "o2")
    report = (tmp_path / "o2" / "REPORT.md").read_text()
    assert "on the trimmed reads the rule would give maxEE 2,4, but the validation tested the trim with maxEE 2,2 only" in report and "after the trim the rule gives 2" not in report


def test_noisy_its_reads_that_mostly_pass_get_maxee_computed_from_the_reads(tmp_path):
    noisy = reads_with_a_bad_reverse_tail(tmp_path / "noisy", share_good=0.6, bad_q=15)       # 40% of reverse reads carry about 4 expected errors
    res, report, summary = run_its(noisy, tmp_path, trim_advice_below=0.5)                    # the 0.3 rule; by default 0.4.1 advises the trim here
    assert res.pass_rate_default == pytest.approx(0.6, abs=0.05) and res.trim_advice is None
    assert res.relaxed and res.cfg.max_ee == (2.0, 5.0) and res.rule["raw"] == [1, 5] and res.recommendation.ee_pass_rate > 0.95
    assert "not DADA2's default of 2: after DADA2's truncQ = 2 cut" in report and "safe on four fungal mock communities" in report
    assert (summary["recommended"]["maxEE_f"], summary["recommended"]["maxEE_r"]) == (2.0, 5.0) and summary["maxEE_rule"]["relaxed"] is True
    assert "maxEE = 2,5" in " ".join(screen_summary(res, tmp_path / "REPORT.md"))


def test_its_reads_that_mostly_fail_get_the_trim_advice_and_its_prediction(tmp_path):
    poor = reads_with_a_bad_reverse_tail(tmp_path / "poor", share_good=0.4, bad_q=12)        # 60% of reverse reads fall to Q12 after base 100
    res, report, summary = run_its(poor, tmp_path)
    t = res.trim_advice
    assert t is not None and res.pass_rate_default == pytest.approx(0.4, abs=0.05) and t["pass_rate_default"] == res.pass_rate_default
    assert t["pass_rate"] > 0.95 and t["window"] == 4 and t["quality"] == 15 and t["checked"] is True
    assert res.cfg.max_ee == 2.0 and not res.relaxed and res.rule["max_ee"][1] > 2                        # maxEE stays 2 for the trimmed reads; the rule is the second option
    assert "## Trim the poor read ends first" in report and "amplitrim trim INPUT_DIR TRIMMED_DIR" in report and "checked against DADA2 on eleven ITS datasets" in report
    assert "second option is maxEE 2," in report and "for the reads after trimming their poor ends" in report
    screen = "\n".join(screen_summary(res, tmp_path / "REPORT.md"))
    assert "trim the poor read ends first" in screen and "amplitrim trim INPUT_DIR TRIMMED_DIR" in screen
    assert summary["trim_advice"]["checked"] is True and summary["trim_advice"]["pass_rate"] > 0.95 and summary["recommended"]["maxEE"] == 2.0


def test_a_16s_run_is_unchanged(reads, tmp_path):
    res = pipeline.run(pipeline.Settings(input_dir=reads, work_dir=tmp_path / "w"))
    assert res.trim_advice is None and res.marker == "16s" and (res.rule is None or res.rule["cut"] is not None)


def fq(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("".join(f"@r{i}\n{'A' * len(q)}\n+\n{''.join(chr(33 + x) for x in q)}\n" for i, q in enumerate(records)))


def test_untruncated_errors_count_up_to_the_first_q2_base_and_only_reads_long_enough(tmp_path):
    fwd = [[30] * 60 + [2] + [30] * 39, [20] * 80, [30] * 40]                      # cut at 60; uncut 80; too short (40 < 50)
    rev = [[10] * 70, [30] * 100, [2] + [30] * 99]                                  # uncut 70; uncut 100; cut at 0 (too short)
    fq(tmp_path / "s_R1.fastq.gz", fwd)
    fq(tmp_path / "s_R2.fastq.gz", rev)
    ef, er = ReadPairs.load(tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz").untruncated_errors(50)
    assert np.allclose(sorted(ef), sorted([60 * 1e-3, 80 * 1e-2])) and np.allclose(sorted(er), sorted([70 * 0.1, 100 * 1e-3]))
