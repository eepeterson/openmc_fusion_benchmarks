import importlib.util
import sys
import sys
import types
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml
import h5py


def _module_available(name: str) -> bool:
    if name in sys.modules:
        return True
    try:
        return importlib.util.find_spec(name) is not None
    except ValueError:
        return False


# Stub heavy dependencies so importing the package works in minimal test envs.
if not _module_available("openmc"):
    openmc_stub = types.ModuleType("openmc")
    openmc_stub.__path__ = []
    for cls_name in (
        "StatePoint",
        "Tally",
        "Filter",
        "CellFilter",
        "SurfaceFilter",
        "MaterialFilter",
        "EnergyFilter",
        "ParticleFilter",
        "Materials",
        "Material",
        "Geometry",
        "Settings",
        "Tallies",
        "Model",
        "DAGMCUniverse",
    ):
        setattr(openmc_stub, cls_name, type(cls_name, (), {}))
    sys.modules.setdefault("openmc", openmc_stub)

    openmc_data_stub = types.ModuleType("openmc.data")
    openmc_data_stub.zam = lambda _name: (1, 1, 0)
    sys.modules.setdefault("openmc.data", openmc_data_stub)

try:
    import pydagmc as _pydagmc  # noqa: F401
except Exception:
    pydagmc_stub = types.ModuleType("pydagmc")
    pydagmc_stub.Model = type("Model", (), {})
    sys.modules["pydagmc"] = pydagmc_stub

try:
    import cad_to_dagmc as _cad_to_dagmc  # noqa: F401
except Exception:
    cad_stub = types.ModuleType("cad_to_dagmc")
    cad_stub.CadToDagmc = type("CadToDagmc", (), {})
    sys.modules["cad_to_dagmc"] = cad_stub

if not _module_available("sandy"):
    sys.modules.setdefault("sandy", types.ModuleType("sandy"))

from openmc_fusion_benchmarks.benchmark_results import BenchmarkResults
from openmc_fusion_benchmarks.report import (
    PlotStyle,
    ReportConfig,
    ReportMetadata,
    ResultSource,
    Report,
    build_report,
)
from openmc_fusion_benchmarks.report import renderers
from openmc_fusion_benchmarks.report.renderers import (
    _collect_observable_metrics,
    _observable_metrics_for_verbosity,
    _quality_metric_description,
    _quality_metric_equation,
    _quality_metrics_for_verbosity,
    render_plots_for_report,
    render_yaml,
)
from openmc_fusion_benchmarks.report import plots as report_plots
from openmc_fusion_benchmarks.tallies import BaseTally


def _write_structured_group(filepath: Path, group: str, tally_name: str) -> None:
    ds = xr.Dataset(
        {
            "mean": xr.DataArray(
                np.arange(4.0).reshape(2, 1, 2),
                dims=("cell", "nuclide", "score"),
                coords={
                    "cell": [0, 1],
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux", "heating"], dtype="U"),
                },
            ),
            "mc_std": xr.DataArray(
                np.full((2, 1, 2), 0.1),
                dims=("cell", "nuclide", "score"),
                coords={
                    "cell": [0, 1],
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux", "heating"], dtype="U"),
                },
            ),
        }
    )
    ds["mean"].attrs["tally_id"] = 1
    ds["mean"].attrs["tally_name"] = tally_name
    ds["mc_std"].attrs["tally_id"] = 1
    ds["mc_std"].attrs["tally_name"] = tally_name
    ds.to_netcdf(filepath, mode="a" if filepath.exists() else "w", engine="h5netcdf", group=group)


def _write_spec_and_metadata(filepath: Path) -> None:
    spec = {"metadata": {"title": "Spec Title", "description": "Spec description"}}
    spec_bytes = yaml.safe_dump(spec).encode("utf-8")

    with h5py.File(filepath, "a") as handle:
        if "specifications" in handle:
            del handle["specifications"]
        spec_group = handle.create_group("specifications")
        spec_group.attrs["format"] = "yaml"
        spec_group.create_dataset("yaml", data=np.bytes_(spec_bytes))

        if "run_metadata" in handle:
            del handle["run_metadata"]
        meta_group = handle.create_group("run_metadata")
        meta_group.attrs["code_name"] = "openmc"
        meta_group.attrs["code_version"] = "0.15.2"


@pytest.fixture
def results_pair(tmp_path):
    exp_path = tmp_path / "reference_results.h5"
    calc_path = tmp_path / "benchmark_results.h5"
    _write_structured_group(exp_path, group="tally_1", tally_name="tally_1")
    _write_structured_group(calc_path, group="tally_1", tally_name="tally_1")
    return BenchmarkResults.from_file(exp_path), BenchmarkResults.from_file(calc_path)


def test_build_report(results_pair, tmp_path):
    exp, calc = results_pair
    metadata = ReportMetadata(title="Test", benchmark_id="bench")
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)

    report = build_report(metadata, sources, config)

    assert report.metadata.title == "Test"
    assert report.data["benchmark_id"] == "bench"
    assert report.data["sources"][0]["name"] == "experiment"
    assert report.plots[0].tally_name == "tally_1"
    assert isinstance(report.plots[0].style, PlotStyle)


def test_build_report_without_metadata(results_pair, tmp_path):
    exp, calc = results_pair
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)

    report = build_report(sources, config)

    assert report.metadata.title
    assert report.data["sources"][0]["name"] == "experiment"
    assert report.data["verbosity"] == 2


def test_build_report_includes_specs_and_run_metadata(tmp_path):
    exp_path = tmp_path / "reference_results.h5"
    calc_path = tmp_path / "benchmark_results.h5"
    _write_structured_group(exp_path, group="tally_1", tally_name="tally_1")
    _write_structured_group(calc_path, group="tally_1", tally_name="tally_1")
    _write_spec_and_metadata(exp_path)
    _write_spec_and_metadata(calc_path)

    exp = BenchmarkResults.from_file(exp_path)
    calc = BenchmarkResults.from_file(calc_path)
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)

    report = build_report(sources, config)

    assert report.data["specifications"]["metadata"]["title"] == "Spec Title"
    assert report.data["run_metadata"]["code_name"] == "openmc"


def test_render_yaml(results_pair, tmp_path):
    exp, calc = results_pair
    metadata = ReportMetadata(title="Test", benchmark_id="bench")
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)

    report = build_report(metadata, sources, config)
    output_path = render_yaml(report, tmp_path / "report.yaml")

    assert output_path.exists()
    payload = yaml.safe_load(output_path.read_text())
    assert payload["benchmark_id"] == "bench"
    assert payload["sources"][1]["name"] == "calculation"


def test_render_plots_for_report_records_plot_entries(results_pair, tmp_path, monkeypatch):
    exp, calc = results_pair
    metadata = ReportMetadata(title="Test", benchmark_id="bench")
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)
    report = build_report(metadata, sources, config)

    class DummyArtifacts:
        def __init__(self, tally):
            self.tally_name = tally
            self.absolute_plot = tmp_path / f"{tally}_absolute.png"
            self.ce_plot = tmp_path / f"{tally}_ce.png"

    def _fake_build_plot_artifacts(tally_name, experiment, calculation, output_dir):
        return DummyArtifacts(tally_name)

    def _fake_render_plots(artifacts, experiment, calculation, style=None):
        artifacts.absolute_plot.write_text("ok")
        artifacts.ce_plot.write_text("ok")

    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.build_plot_artifacts",
        _fake_build_plot_artifacts,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.render_plots",
        _fake_render_plots,
    )

    entries = render_plots_for_report(report, tmp_path / "plots")

    assert entries
    assert report.data["plots"][0]["tally"] == "tally_1"
    assert Path(entries[0]["absolute_plot"]).exists()


def test_render_plots_for_report_quality_entries(results_pair, tmp_path, monkeypatch):
    exp, calc = results_pair
    metadata = ReportMetadata(title="Test", benchmark_id="bench")
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path, verbosity=2)
    report = build_report(metadata, sources, config)

    class DummyArtifacts:
        def __init__(self, tally):
            self.tally_name = tally
            self.absolute_plot = tmp_path / f"{tally}_absolute.png"
            self.ce_plot = tmp_path / f"{tally}_ce.png"

    class DummyQualityArtifacts:
        def __init__(self, tally, metrics):
            self.tally_name = tally
            self.metric_plots = {m: tmp_path / f"{tally}_{m}.png" for m in metrics}

    def _fake_build_plot_artifacts(tally_name, experiment, calculation, output_dir):
        return DummyArtifacts(tally_name)

    def _fake_render_plots(artifacts, experiment, calculation, style=None):
        artifacts.absolute_plot.write_text("ok")
        artifacts.ce_plot.write_text("ok")

    def _fake_build_quality_plot_artifacts(tally_name, output_dir, metrics):
        return DummyQualityArtifacts(tally_name, metrics)

    def _fake_render_quality_plots(artifacts, experiment, calculation, metrics, style=None):
        for path in artifacts.metric_plots.values():
            path.write_text("ok")

    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.build_plot_artifacts",
        _fake_build_plot_artifacts,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.render_plots",
        _fake_render_plots,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.build_quality_plot_artifacts",
        _fake_build_quality_plot_artifacts,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.render_quality_plots",
        _fake_render_quality_plots,
    )

    render_plots_for_report(report, tmp_path / "plots")

    quality = report.data.get("quality_plots", [])
    assert quality
    assert quality[0]["tally"] == "tally_1"
    assert "ce" in quality[0]["metrics"]


def test_plot_spec_reference_candidate_aliases(results_pair, tmp_path):
    exp, calc = results_pair
    metadata = ReportMetadata(title="Test", benchmark_id="bench")
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)
    report = build_report(metadata, sources, config)

    plot = report.plots[0]
    assert plot.reference is plot.experiment
    assert plot.candidate is plot.calculation


def test_collect_observable_metrics(results_pair, tmp_path):
    exp, calc = results_pair
    metadata = ReportMetadata(title="Test", benchmark_id="bench")
    sources = [
        ResultSource(name="experiment", kind="experiment", results=exp),
        ResultSource(name="calculation", kind="calculation", results=calc),
    ]
    config = ReportConfig(output_dir=tmp_path)
    report = build_report(metadata, sources, config)

    entries = _collect_observable_metrics(report)
    assert entries
    assert entries[0]["tally"] == "tally_1"
    assert "rms_relative_deviation" in entries[0]


def test_quality_metrics_for_verbosity():
    assert _quality_metrics_for_verbosity(0) == ["ce"]
    assert _quality_metrics_for_verbosity(1) == ["ce", "chi2_contribution"]
    assert "relative_deviation" in _quality_metrics_for_verbosity(2)
    assert "normalized_residual" in _quality_metrics_for_verbosity(3)


def test_observable_metrics_for_verbosity():
    assert _observable_metrics_for_verbosity(0) == ["rms_relative_deviation"]
    assert "reduced_chi2" in _observable_metrics_for_verbosity(1)
    assert "mean_abs_normalized_residual" in _observable_metrics_for_verbosity(3)


def test_quality_metric_equations_and_descriptions():
    assert _quality_metric_equation("ce")
    assert _quality_metric_description("ce")
    assert _quality_metric_equation("unknown") == ""
    assert _quality_metric_description("unknown") == ""


def test_plot_utilities_auto_scale_and_default_x():
    values = np.array([1.0, 1000.0])
    assert report_plots._auto_scale(values, threshold=100.0) == "log"
    assert report_plots._auto_scale(np.array([-1.0, 2.0]), threshold=100.0) == "linear"

    da = xr.DataArray(
        np.array([1.0, 2.0, 3.0]),
        dims=("energy",),
        coords={"energy": np.array([0.1, 1.0, 10.0])},
    )
    tally = BaseTally(da)
    x_vals = report_plots._default_x(tally)
    assert np.allclose(x_vals, np.array([0.1, 1.0, 10.0]))


def test_plot_utilities_flatten_tally():
    da = xr.DataArray(
        np.arange(6.0).reshape(2, 3),
        dims=("cell", "energy"),
    )
    tally = BaseTally(da)
    flat = report_plots._flatten_tally(tally)
    assert flat.shape == (6,)


def _install_fake_matplotlib(monkeypatch):
    class DummyLine:
        def __init__(self, color="#000"):
            self._color = color

        def get_color(self):
            return self._color

    class DummyAx:
        def imshow(self, *_args, **_kwargs):
            return None

        def axis(self, *_args, **_kwargs):
            return None

        def set_title(self, *_args, **_kwargs):
            return None

        def text(self, *_args, **_kwargs):
            return None

        def bar(self, *_args, **_kwargs):
            return None

        def grid(self, *_args, **_kwargs):
            return None

        def set_xticks(self, *_args, **_kwargs):
            return None

        def set_xticklabels(self, *_args, **_kwargs):
            return None

    class DummyFig:
        def text(self, *_args, **_kwargs):
            return None

        def add_axes(self, *_args, **_kwargs):
            return DummyAx()

        def tight_layout(self, *_args, **_kwargs):
            return None

    class DummyPdfPages:
        def __init__(self, path):
            self.path = Path(path)
            self.path.write_text("pdf")

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return None

        def savefig(self, *_args, **_kwargs):
            return None

    class DummyPlt:
        def figure(self, *_args, **_kwargs):
            return DummyFig()

        def imread(self, *_args, **_kwargs):
            return np.zeros((2, 2, 3))

        def close(self, *_args, **_kwargs):
            return None

        def subplots(self, rows, cols, **_kwargs):
            fig = DummyFig()
            axes = np.array([[DummyAx() for _ in range(cols)] for _ in range(rows)])
            return fig, axes

        def plot(self, *_args, **_kwargs):
            return [DummyLine()]

        def fill_between(self, *_args, **_kwargs):
            return None

        def axhline(self, *_args, **_kwargs):
            return DummyLine()

        def xlabel(self, *_args, **_kwargs):
            return None

        def ylabel(self, *_args, **_kwargs):
            return None

        def yscale(self, *_args, **_kwargs):
            return None

        def xscale(self, *_args, **_kwargs):
            return None

        def title(self, *_args, **_kwargs):
            return None

        def legend(self, *_args, **_kwargs):
            return None

        def tight_layout(self, *_args, **_kwargs):
            return None

        def savefig(self, path, **_kwargs):
            Path(path).write_text("png")

    import types

    fake_pdf = types.SimpleNamespace(PdfPages=DummyPdfPages)
    fake_plt = DummyPlt()
    monkeypatch.setitem(sys.modules, "matplotlib.backends.backend_pdf", fake_pdf)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", fake_plt)


def test_render_pdf_with_quality_and_observable(tmp_path, monkeypatch):
    _install_fake_matplotlib(monkeypatch)

    report = Report(
        metadata=ReportMetadata(title="Title", benchmark_id="bench"),
        sources=[],
        plots=[],
        data={
            "verbosity": 2,
            "specifications": {"metadata": {"description": "desc"}, "materials": []},
            "sources": [{"kind": "experiment", "run_metadata": {"code_name": "openmc"}}],
        },
    )

    quality_plot = tmp_path / "ce.png"
    quality_plot.write_text("img")
    plot_entries = [
        {"tally": "t1", "absolute_plot": str(quality_plot), "ce_plot": str(quality_plot)}
    ]

    def _fake_render_plots_for_report(_report, _plots_dir):
        _report.data["quality_plots"] = [
            {"tally": "t1", "metrics": {"ce": str(quality_plot)}}
        ]
        return plot_entries

    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers.render_plots_for_report",
        _fake_render_plots_for_report,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.report.renderers._collect_observable_metrics",
        lambda _report: [{"tally": "t1", "rms_relative_deviation": 0.1}],
    )

    output_path = renderers.render_pdf(report, tmp_path / "report.pdf", tmp_path / "plots")
    assert output_path.exists()


def test_render_plots_and_quality_plots(tmp_path, monkeypatch):
    _install_fake_matplotlib(monkeypatch)

    mean = xr.DataArray(
        np.array([1.0, 2.0, 3.0]),
        dims=("cell",),
        coords={"cell": [0, 1, 2]},
    )
    std = xr.DataArray(np.array([0.1, 0.1, 0.1]), dims=("cell",), coords=mean.coords)
    exp = BaseTally(mean, std)
    calc = BaseTally(mean * 1.1, std)

    artifacts = report_plots.build_plot_artifacts("tally", exp, calc, tmp_path)
    report_plots.render_plots(artifacts, exp, calc, style=PlotStyle())
    assert artifacts.absolute_plot.exists()
    assert artifacts.ce_plot.exists()

    quality = report_plots.build_quality_plot_artifacts("tally", tmp_path, ["ce"])
    report_plots.render_quality_plots(quality, exp, calc, ["ce"], style=PlotStyle())
    assert quality.metric_plots["ce"].exists()


def test_build_spec_sections_and_format_helpers():
    spec = {
        "metadata": {"title": "Demo", "id": "demo", "geometry": "skip", "settings": "skip"},
        "materials": [{"id": 1, "name": "mat", "density": 1.0, "composition": {"H": 1}}],
        "tallies": [{"name": "t1", "particle": "neutron", "scores": ["flux"], "filters": []}],
        "geometry": {"cad_file": "demo.step", "meshing": {}},
        "settings": {"run_mode": "fixed_source", "batches": 1, "particles_per_batch": 10, "photon_transport": False},
    }

    sections = renderers._build_spec_sections(spec, verbosity=2)
    titles = [title for title, _body in sections]
    assert "Metadata" in titles
    assert "Materials" in titles
    assert "Tallies" in titles
    assert "Geometry" in titles
    assert "Settings" in titles

    inline = renderers._format_inline({"a": 1})
    assert "\"a\":1" in inline

    formatted = renderers._format_key_fields([{"a": 1, "b": 2}], keys=["a"], max_items=1)
    assert "\"a\":1" in formatted


def test_format_run_metadata_and_reference_lines():
    run_meta = {
        "code_name": "openmc",
        "code_version": "0.15.2",
        "nuclear_data_name": "endfb",
        "nuclear_data_version": "8.0",
    }
    line = renderers._format_run_metadata(run_meta)
    assert "openmc 0.15.2" in line
    assert "endfb 8.0" in line

    report = Report(
        metadata=ReportMetadata(title="Title", benchmark_id="bench"),
        sources=[],
        plots=[],
        data={"sources": [{"kind": "experiment", "run_metadata": run_meta}]},
    )

    assert renderers._format_reference(report)
    assert renderers._format_validation(report) == ""


def test_render_helper_lines_and_description_defaults():
    report = Report(
        metadata=ReportMetadata(title="Title", benchmark_id="bench", description="desc"),
        sources=[],
        plots=[],
        data={},
    )
    assert renderers._wrap_lines("", 10) == []
    assert renderers._resolve_benchmark_description(report) == "desc"
    assert renderers._format_reference(report) == ""
    assert renderers._format_validation(report) == ""


def test_render_specifications_invalid_spec_noop(tmp_path, monkeypatch):
    _install_fake_matplotlib(monkeypatch)

    report = Report(
        metadata=ReportMetadata(title="Title", benchmark_id="bench"),
        sources=[],
        plots=[],
        data={"specifications": "not-a-dict"},
    )

    class DummyPdf:
        def savefig(self, *_args, **_kwargs):
            return None

    renderers._render_specifications(DummyPdf(), report, verbosity=2)


def test_render_specifications_with_fake_matplotlib(tmp_path, monkeypatch):
    _install_fake_matplotlib(monkeypatch)

    report = Report(
        metadata=ReportMetadata(title="Title", benchmark_id="bench"),
        sources=[],
        plots=[],
        data={
            "specifications": {
                "metadata": {"title": "Demo"},
                "materials": [],
                "tallies": [],
            }
        },
    )

    class DummyPdf:
        def savefig(self, *_args, **_kwargs):
            return None

    renderers._render_specifications(DummyPdf(), report, verbosity=1)


def test_render_quality_section_handles_missing_and_bad_images(monkeypatch):
    _install_fake_matplotlib(monkeypatch)

    def _bad_imread(_path):
        raise ValueError("bad image")

    monkeypatch.setattr("matplotlib.pyplot.imread", _bad_imread)

    class DummyPdf:
        def savefig(self, *_args, **_kwargs):
            return None

    quality_entries = [
        {"tally": "t1", "metrics": {}},
        {"tally": "t2", "metrics": {"ce": "bad.png"}},
    ]

    renderers._render_quality_section(DummyPdf(), quality_entries, verbosity=1)


def test_render_observable_summary_with_fake_matplotlib(monkeypatch):
    _install_fake_matplotlib(monkeypatch)

    class DummyPdf:
        def savefig(self, *_args, **_kwargs):
            return None

    observable_entries = [
        {"tally": "t1", "rms_relative_deviation": 0.1, "reduced_chi2": 1.2},
        {"tally": "t2", "rms_relative_deviation": 0.2, "reduced_chi2": 0.9},
    ]

    renderers._render_observable_summary(DummyPdf(), observable_entries, verbosity=1)


def test_build_spec_sections_low_verbosity_and_format_helpers():
    spec = {
        "metadata": {"title": "Demo", "settings": "skip"},
        "materials": [{"name": "mat1"}, {"name": "mat2"}],
        "tallies": [{"name": "t1"}, {"name": "t2"}],
        "geometry": {"cad_file": "demo.step"},
        "settings": {"run_mode": "fixed_source"},
    }

    sections = renderers._build_spec_sections(spec, verbosity=1)
    titles = [title for title, _body in sections]
    assert "Materials" in titles
    assert "Tallies" in titles

    formatted = renderers._format_key_fields(["raw", {"a": 1, "b": 2}], keys=["a"], max_items=1)
    assert "\"raw\"" in formatted
    assert "..." in formatted


def test_render_pdf_import_error(monkeypatch, tmp_path):
    def _raise_import(name, *args, **kwargs):
        if name.startswith("matplotlib"):
            raise ImportError("no mpl")
        return __import__(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", _raise_import)

    report = Report(
        metadata=ReportMetadata(title="Title", benchmark_id="bench"),
        sources=[],
        plots=[],
        data={},
    )

    with pytest.raises(RuntimeError, match="matplotlib is required"):
        renderers.render_pdf(report, tmp_path / "report.pdf", tmp_path / "plots")
