"""The command-line tool: what it accepts, what it prints, what it writes, and how it fails."""
import json
import re

import pytest

from amplitrim import cli
from conftest import need_vsearch

REQUIRED = ["--input-dir", "in", "--work-dir", "work", "--out-dir", "out"]


@pytest.mark.parametrize("option", ["--region", "--primer-f-len", "--splits", "--a-values", "--max-ee", "--marker", "--min-overlap",
                                    "--figaro-script", "--explore-max-ee", "--explore-coverage", "--exhaustive", "--longest-insert"])
def test_removed_options_are_rejected(option):
    """AmpliTrim recommends; it does not take the settings it is meant to recommend, and it has no tuning options."""
    with pytest.raises(SystemExit):
        cli.main(REQUIRED + [option, "x"])


def test_options_are_three_required_arguments_then_two_optional_ones():
    ap = cli.build_parser()
    groups = {g.title: [a.option_strings[0] for a in g._group_actions] for g in ap._action_groups}
    assert groups["required arguments"] == ["--input-dir", "--work-dir", "--out-dir"]
    assert groups["optional arguments"] == ["--its", "--verbose"]
    assert [a.option_strings[0] for a in ap._actions if a.required] == ["--input-dir", "--work-dir", "--out-dir"]


def test_its_is_a_flag_and_the_default_is_the_16s_search():
    args = cli.build_parser().parse_args(REQUIRED)
    assert args.its is False and args.verbose is False
    assert cli.build_parser().parse_args(REQUIRED + ["--its"]).its is True


def test_missing_vsearch_is_reported_plainly(monkeypatch, capsys):
    monkeypatch.setattr("amplitrim.pipeline.shutil.which", lambda name: None)
    assert cli.main(REQUIRED) == 3
    assert "vsearch was not found" in capsys.readouterr().err


@need_vsearch
def test_writes_exactly_three_outputs(cli_run):
    assert cli_run.rc == 0
    assert {p.name for p in cli_run.out.iterdir()} == {"REPORT.md", "summary.png", "summary.json"}


@need_vsearch
def test_screen_shows_the_settings_in_a_few_plain_lines(cli_run):
    lines = [ln for ln in cli_run.stdout.splitlines() if ln.strip()]
    assert any(ln.startswith("Recommended DADA2 settings") for ln in lines)
    assert any(re.fullmatch(r"\s+truncLenF = \d+\s+truncLenR = \d+\s+maxEE = 2\s+minOverlap = 12", ln) for ln in lines)
    assert any("read pairs analysed should pass filtering and merge" in ln for ln in lines)
    assert any(ln.startswith("Report: ") for ln in lines)
    assert len(lines) <= 9
    assert not re.search(r"\b(INFO|DEBUG)\b|20\d\d-\d\d-\d\d|FIGARO|nt: ", cli_run.stdout)   # no log noise or progress


@need_vsearch
def test_screen_and_summary_json_agree_and_cover_the_longest_amplicon(cli_run):
    summary = json.loads((cli_run.out / "summary.json").read_text())
    rec = summary["recommended"]
    assert f"truncLenF = {rec['truncLenF']}" in cli_run.stdout
    assert rec["longest_mergeable_insert"] >= 434 and rec["estimated_yield"] > 0.7   # simulated max insert: 434
    assert summary["maxEE_rule"]["relaxed"] is False and summary["maxEE_rule"]["pass_rate_at_maxEE_2"] > 0.7


@need_vsearch
def test_its_flag_recommends_no_truncation_and_writes_the_same_three_files(reads, tmp_path, capsys):
    args = ["--input-dir", str(reads), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o"), "--its"]
    assert cli.main(args) == 0
    out = capsys.readouterr().out
    assert "no fixed truncation (truncLen = 0)   maxEE = 2   minLen = 50   minOverlap = 12" in out
    assert {p.name for p in (tmp_path / "o").iterdir()} == {"REPORT.md", "summary.png", "summary.json"}
    rec = json.loads((tmp_path / "o" / "summary.json").read_text())["recommended"]
    assert (rec["truncLenF"], rec["truncLenR"], rec["minLen"]) == (0, 0, 50)
