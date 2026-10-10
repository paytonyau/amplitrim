"""The commands that the README shows for trying the tool on a small simulated dataset work, and give the recommendation that the README describes."""
import gzip
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from conftest import need_vsearch                    # noqa: E402
from amplitrim import cli                            # noqa: E402


def readme_commands():
    """The generator command and the AmpliTrim command that the README shows together."""
    block = re.search(r"```bash\n(python tests/simulate\.py[^\n]*\namplitrim[^\n]*)\n```", (ROOT / "README.md").read_text(encoding="utf-8"))
    assert block, "the README must show the simulate and amplitrim commands together in one bash block"
    gen, run = block.group(1).split("\n")
    return shlex.split(gen)[2:], shlex.split(run)[1:]


def test_the_readme_command_makes_the_reads_that_the_cli_tests_run_on(tmp_path, sim):
    """Run the README's generator command: it must write the same reads (once decompressed) as the fixture that the CLI tests use."""
    args, _ = readme_commands()
    assert args[0] == "demo"
    done = subprocess.run([sys.executable, str(ROOT / "tests" / "simulate.py"), str(tmp_path / "demo"), *args[1:]], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    made = sorted(p.name for p in (tmp_path / "demo" / "primer_trimmed").glob("*.fastq.gz"))
    assert made == ["sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz", "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"]
    for name in made:
        with gzip.open(tmp_path / "demo" / "primer_trimmed" / name, "rt") as a, gzip.open(sim / "primer_trimmed" / name, "rt") as b:
            text = a.read()
            assert text == b.read(), name
        assert text.count("\n") == 4 * 1500                                                              # 1,500 reads of four lines
    assert (tmp_path / "demo" / "truth.tsv").read_text() == (sim / "truth.tsv").read_text()


def test_the_readme_amplitrim_command_is_valid_and_reads_the_folder_that_the_generator_writes():
    gen, run = readme_commands()
    try:
        parsed = cli.build_parser().parse_args(run)
    except SystemExit as err:                                                                            # argparse exits on an option that does not exist
        raise AssertionError(f"the README's amplitrim command is not valid for the command line: {run}") from err
    assert parsed.input_dir == Path(gen[0]) / "primer_trimmed" and parsed.work_dir == Path("work") and parsed.out_dir == Path("out")


@need_vsearch
def test_the_readme_run_gives_the_recommendation_it_describes(tmp_path, sim):
    """The data the README makes are the data of the fixture (first test), so run on the fixture: the README's claims about the recommendation."""
    rc = cli.main(["--input-dir", str(sim / "primer_trimmed"), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o")])
    assert rc == 0 and {p.name for p in (tmp_path / "o").iterdir()} == {"REPORT.md", "summary.png", "summary.json"}
    summary = json.loads((tmp_path / "o" / "summary.json").read_text())
    assert summary["recommended"]["longest_mergeable_insert"] >= 434                                      # covers the longest simulated amplicon
    assert summary["maxEE_rule"]["relaxed"] is False                                                      # maxEE stays at 2
    assert summary["recommended"]["estimated_yield"] > 0.70                                               # more than 70% of the pairs
