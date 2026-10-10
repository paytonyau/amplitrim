"""Shared fixtures. Tests that need vsearch are skipped without it."""
import contextlib
import io
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from simulate import write_dataset  # noqa: E402

need_vsearch = pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")


@pytest.fixture(scope="session")
def sim(tmp_path_factory):
    """Simulated mock community with ground truth (true insert lengths 423-434)."""
    d = tmp_path_factory.mktemp("sim")
    write_dataset(d, n_samples=2, pairs_per_sample=1500, seed=11)
    return d


@pytest.fixture(scope="session")
def reads(sim):
    """The simulated reads as the tool expects them: adapters and primers already removed."""
    return sim / "primer_trimmed"


def _run_cli(reads, root, extra):
    from amplitrim.cli import main
    inp = root / "input"
    inp.mkdir()
    for i, r1 in enumerate(sorted(reads.glob("*_R1.fastq.gz")), 1):          # Illumina-style names
        base = r1.name.replace("_R1.fastq.gz", "")
        shutil.copy(r1, inp / f"{base}_S{i}_L001_R1_001.fastq.gz")
        shutil.copy(r1.with_name(f"{base}_R2.fastq.gz"), inp / f"{base}_S{i}_L001_R2_001.fastq.gz")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = main(["--input-dir", str(inp), "--work-dir", str(root / "work"), "--out-dir", str(root / "out")] + extra)
    return SimpleNamespace(rc=rc, out=root / "out", work=root / "work", stdout=buf.getvalue())


def _run_api(reads, root, **settings):
    """The same run through the Python interface, for the options the command line does not offer."""
    from amplitrim.pipeline import Settings, run
    from amplitrim.report import screen_summary, write_outputs
    res = run(Settings(input_dir=reads, work_dir=root / "work", **settings))
    report = write_outputs(res, root / "out")
    return SimpleNamespace(rc=0 if res.recommendation else 2, out=root / "out", work=root / "work",
                           stdout="\n".join(screen_summary(res, report)))


@pytest.fixture(scope="session")
def cli_run(reads, tmp_path_factory):
    """One real run of the command-line tool."""
    return _run_cli(reads, tmp_path_factory.mktemp("run"), [])


@pytest.fixture(scope="session")
def cli_run_ee25(reads, tmp_path_factory):
    """One real run at a given maxEE of 2,5 (the Python interface takes a maxEE; the command line chooses it)."""
    return _run_api(reads, tmp_path_factory.mktemp("runee"), max_ee=(2.0, 5.0))


@pytest.fixture(scope="session")
def cli_run_explore(reads, tmp_path_factory):
    """One real run that also finds the best truncation for several maxEE values (Python interface only)."""
    return _run_api(reads, tmp_path_factory.mktemp("runsug"), explore_max_ee=True)


@pytest.fixture(scope="session")
def cli_run_coverage(reads, tmp_path_factory):
    """One real run that also finds the best setting for several coverage levels (Python interface only)."""
    return _run_api(reads, tmp_path_factory.mktemp("runcov"), explore_coverage=True)
