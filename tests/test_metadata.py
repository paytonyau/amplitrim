"""Release metadata: the README works on PyPI, the version is written once and agrees everywhere, the built metadata is complete and claims only what CI tests, and the source archive can run the tests. No network."""
import atexit
import functools
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import amplitrim                                                                                 # noqa: E402

BLOB = "https://github.com/paytonyau/amplitrim/blob/main/"
read = lambda name: (ROOT / name).read_text(encoding="utf-8")


def slug(h):
    return re.sub(r"\s+", "-", re.sub(r"[^\w\s-]", "", h.strip().lower()))


@functools.lru_cache(maxsize=1)
def built():
    """Build the source archive and the wheel in a temporary copy of the repository; returns (sdist path, wheel path)."""
    try:
        import setuptools
    except ImportError:
        pytest.skip("setuptools is not installed")
    if int(setuptools.__version__.split(".")[0]) < 68:
        pytest.skip("setuptools older than 68 cannot read the project table")
    work = Path(tempfile.mkdtemp(prefix="amplitrim_build_"))
    atexit.register(shutil.rmtree, str(work), True)
    src, out = work / "src", work / "out"
    shutil.copytree(ROOT, src, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc", "*.egg-info", "build", "dist", ".pytest_cache"))
    out.mkdir()
    code = "import setuptools.build_meta as b; print('SDIST=' + b.build_sdist(r'%s')); print('WHEEL=' + b.build_wheel(r'%s'))" % (out, out)
    done = subprocess.run([sys.executable, "-c", code], cwd=str(src), capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-600:]
    sdist, wheel = re.search(r"^SDIST=(.+)$", done.stdout, flags=re.M), re.search(r"^WHEEL=(.+)$", done.stdout, flags=re.M)      # setuptools prints its log to the same stream
    assert sdist and wheel, done.stdout[-400:]
    return out / sdist.group(1).strip(), out / wheel.group(1).strip()


def wheel_metadata():
    z = zipfile.ZipFile(built()[1])
    return z, z.read([n for n in z.namelist() if n.endswith("METADATA")][0]).decode("utf-8")


def fields(meta, key):
    return [l.split(": ", 1)[1] for l in meta.split("\n\n", 1)[0].splitlines() if l.startswith(key + ": ")]


def test_the_readme_links_work_on_pypi():
    """PyPI shows the README without the repository, so a relative link breaks there. Links into the repository must also point at something that exists."""
    body = re.sub(r"```.*?```", "", read("README.md"), flags=re.S)
    targets = re.findall(r"\]\(([^)\s]+)\)", body)
    assert not [t for t in targets if not t.startswith(("https://", "#", "mailto:"))], "relative link in the README"
    for t in targets:
        if not t.startswith(BLOB):
            continue
        path, _, anchor = t[len(BLOB):].partition("#")
        assert (ROOT / path).exists(), t
        if anchor and path.endswith(".md"):
            heads = {slug(re.sub(r"^#+\s*", "", l)) for l in read(path).splitlines() if l.startswith("#")}
            assert anchor in heads, t


def test_the_version_is_written_once_and_agrees_everywhere():
    v = amplitrim.__version__
    pyproject = read("pyproject.toml")
    assert 'dynamic = ["version"]' in pyproject and not re.search(r'^version\s*=\s*"', pyproject, flags=re.M)      # not a second copy in pyproject
    assert re.search(r'^version:\s*"?([^"\n]+)"?', read("CITATION.cff"), flags=re.M).group(1) == v
    assert re.search(r"^## %s\b" % re.escape(v), read("CHANGELOG.md"), flags=re.M)
    assert "AmpliTrim %s\n" % v in read("README.md")                                                  # the example output in the README
    assert read("docs/index.html").count("AmpliTrim %s" % v) >= 2                                      # the example output and the footer of the page


def test_the_built_metadata_is_complete_and_claims_only_what_ci_tests():
    z, meta = wheel_metadata()
    v = amplitrim.__version__
    cff = read("CITATION.cff")
    author = "%s %s" % (re.search(r'given-names:\s*"([^"]+)"', cff).group(1), re.search(r'family-names:\s*"([^"]+)"', cff).group(1))
    assert fields(meta, "Name") == ["amplitrim"] and fields(meta, "Version") == [v] and fields(meta, "Requires-Python") == [">=3.9"]
    assert fields(meta, "Author") == [author]                                                          # the name in the citation file
    assert fields(meta, "Description-Content-Type") == ["text/markdown"] and fields(meta, "License-File") == ["LICENSE"]
    runtime = {r.replace(" ", "") for r in fields(meta, "Requires-Dist") if "extra ==" not in r}
    assert runtime == {"numpy>=1.24", "matplotlib>=3.7"}, runtime                                      # the tests' pandas and pytest stay in the extra
    classifiers = set(fields(meta, "Classifier"))
    for c in ("Environment :: Console", "Operating System :: OS Independent", "Topic :: Scientific/Engineering :: Bio-Informatics", "License :: OSI Approved :: MIT License"):
        assert c in classifiers, c
    ci = read(".github/workflows/tests.yml"); matrix = ci[ci.index("matrix:"):ci.index("runs-on:")]
    tested = set(re.findall(r'"(3\.\d+)"', matrix))
    claimed = {c.rsplit(" :: ", 1)[1] for c in classifiers if re.search(r"Python :: 3\.\d+$", c)}
    assert claimed == tested, (claimed, tested)                                                        # a Python version is claimed only if CI runs it
    urls = {u.split(", ", 1)[0]: u.split(", ", 1)[1] for u in fields(meta, "Project-URL")}
    assert {"Homepage", "Documentation", "Source", "Issues", "Changelog", "Citation"} <= set(urls), set(urls)
    assert all(u.startswith("https://github.com/paytonyau/") for u in urls.values())
    for u in urls.values():
        if u.startswith(BLOB):
            assert (ROOT / u[len(BLOB):].split("#")[0]).exists(), u
    names = z.namelist()
    assert any(n.endswith("LICENSE") for n in names) and not any(n.startswith(("tests", "docs")) for n in names)
    assert "amplitrim = amplitrim.cli:main" in z.read([n for n in names if n.endswith("entry_points.txt")][0]).decode()


def test_the_source_archive_can_run_the_tests():
    """The archive carries the changelog, the citation file, the guide and the whole tests folder (helpers and fixtures too, not only test_*.py)."""
    t = tarfile.open(str(built()[0]))
    names = {m.name.split("/", 1)[1] for m in t.getmembers() if m.isfile() and "/" in m.name}
    need = {"LICENSE", "README.md", "CHANGELOG.md", "CITATION.cff", "CONTRIBUTING.md", "pyproject.toml", "docs/usage.md", "docs/index.html",
            "tests/conftest.py", "tests/simulate.py", "tests/fake_vsearch.py"}
    assert not need - names, need - names
    assert {"tests/" + p.name for p in (ROOT / "tests").glob("test_*.py")} <= names
    assert any(n.startswith("tests/fixtures/") for n in names)
    assert not [n for n in names if "__pycache__" in n or n.endswith(".pyc") or n.startswith((".git/", ".github/"))]
