from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from ..benchmark_results import BenchmarkResults


@dataclass
class ReportMetadata:
    title: str
    benchmark_id: str
    description: str = ""
    model_description: str = ""
    code_name: str = ""
    code_version: str = ""
    notes: str = ""


@dataclass
class ResultSource:
    name: str
    kind: str
    results: BenchmarkResults
    tally_names: Optional[Sequence[str]] = None


@dataclass
class ReportConfig:
    output_dir: Path
    include_yaml: bool = True
    include_pdf: bool = False
    plot_tallies: Optional[Sequence[str]] = None
    verbosity: int = 2


@dataclass
class PlotStyle:
    x_label: str = "index"
    y_label: str = "value"
    ce_y_label: str = "C/E"
    x_scale: str = "linear"
    y_scale: str = "auto"
    ce_y_scale: str = "linear"
    title: Optional[str] = None
    ce_title: Optional[str] = None
    auto_scale_threshold: float = 100.0


@dataclass
class PlotSpec:
    tally_name: str
    experiment: BenchmarkResults
    calculation: BenchmarkResults
    style: PlotStyle = field(default_factory=PlotStyle)

    @property
    def reference(self) -> BenchmarkResults:
        return self.experiment

    @property
    def candidate(self) -> BenchmarkResults:
        return self.calculation


@dataclass
class Report:
    metadata: ReportMetadata
    sources: list[ResultSource]
    plots: list[PlotSpec] = field(default_factory=list)
    data: dict = field(default_factory=dict)
