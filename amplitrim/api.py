"""Programmatic entry point, shared by the command line and by wrappers (Galaxy, QIIME 2, scripts)."""
import json
from pathlib import Path
from typing import Callable, Optional, Tuple

from .pipeline import PipelineError, Result, Settings, run
from .report import write_outputs

__all__ = ["recommend", "PipelineError"]


def recommend(input_dir, work_dir, out_dir, its: bool = False, min_overlap: int = 12, say: Optional[Callable[[str], None]] = None, **settings) -> Tuple[Result, Path, dict]:
    """Analyse the paired reads in input_dir, write REPORT.md, summary.json (and summary.png when matplotlib is installed)
    to out_dir, and return (result, path of the report, the parsed summary.json).
    Raises PipelineError or RuntimeError for bad input or a missing tool (vsearch)."""
    settings = Settings(input_dir=Path(input_dir), work_dir=Path(work_dir), marker="its" if its else "16s", min_overlap=min_overlap, **settings)     # 0.4: any setting of Settings
    result = run(settings, say=say or (lambda msg: None))
    report = write_outputs(result, Path(out_dir))
    summary = json.loads((Path(out_dir) / "summary.json").read_text(encoding="utf-8"))
    return result, report, summary
