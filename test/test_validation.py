import importlib.util
import sys
import types

import numpy as np
import pytest
import xarray as xr
import h5py


def _module_available(name: str) -> bool:
    if name in sys.modules:
        return True
    try:
        return importlib.util.find_spec(name) is not None
    except ValueError:
        return False


# Keep imports stable when OpenMC stack is unavailable locally.
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
from openmc_fusion_benchmarks.tallies import BaseTally
from openmc_fusion_benchmarks.validation.adapters import (
    compare_benchmark_results,
    compare_tallies,
    datapoints_from_tally,
)
from openmc_fusion_benchmarks.validation import comparison as comparison_module
from openmc_fusion_benchmarks.validation.comparison import compare_point_set
from openmc_fusion_benchmarks.validation.metrics import compute_point_metrics
from openmc_fusion_benchmarks.validation.model import (
    BenchmarkStatus,
    DataPoint,
    ObservableComparison,
    PointComparison,
    PointStatus,
)


def _write_results_file(filepath, group="tally_1", tally_name="tally_1"):
    ds = xr.Dataset(
        {
            "mean": xr.DataArray(
                np.array([[1.0], [2.0]]).reshape(2, 1, 1),
                dims=("cell", "nuclide", "score"),
                coords={
                    "cell": [1, 2],
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux"], dtype="U"),
                },
            ),
            "mc_std": xr.DataArray(
                np.full((2, 1, 1), 0.1),
                dims=("cell", "nuclide", "score"),
                coords={
                    "cell": [1, 2],
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux"], dtype="U"),
                },
            ),
        }
    )
    ds["mean"].attrs["tally_id"] = 1
    ds["mean"].attrs["tally_name"] = tally_name
    ds["mc_std"].attrs["tally_id"] = 1
    ds["mc_std"].attrs["tally_name"] = tally_name
    ds.to_netcdf(filepath, mode="w", engine="h5netcdf", group=group)


def _add_run_metadata(filepath, code_name="openmc", code_version="0.14.0"):
    with h5py.File(filepath, "a") as handle:
        if "run_metadata" in handle:
            del handle["run_metadata"]
        group = handle.create_group("run_metadata")
        group.attrs["code_name"] = code_name
        group.attrs["code_version"] = code_version


def test_compute_point_metrics_grading_toggle():
    point = PointComparison(
        id="p0",
        observable_type="tally",
        experiment=DataPoint(value=1.0, uncertainty=0.1),
        calculation=DataPoint(value=1.05, uncertainty=0.1),
    )

    metrics = compute_point_metrics(point, include_grading=False)
    assert metrics.status is None
    assert metrics.c_over_e == pytest.approx(1.05)

    metrics = compute_point_metrics(point, include_grading=True)
    assert metrics.status in {PointStatus.OK, PointStatus.WARNING, PointStatus.OUTLIER}


def test_compare_point_set_point_ids_length_mismatch():
    with pytest.raises(ValueError, match="point_ids length must match"):
        compare_point_set(
            observable_name="tally",
            observable_type="tally",
            experiment_points=[{"value": 1.0, "uncertainty": 0.1}],
            calculation_points=[{"value": 1.0, "uncertainty": 0.1}],
            point_ids=["a", "b"],
        )


def test_compare_point_set_invalid_item_type():
    with pytest.raises(TypeError, match="Unsupported point type"):
        compare_point_set(
            observable_name="tally",
            observable_type="tally",
            experiment_points=["bad"],
            calculation_points=[{"value": 1.0, "uncertainty": 0.1}],
        )


def test_score_component_invalid_thresholds():
    with pytest.raises(ValueError, match="bad > good"):
        comparison_module._score_component(0.1, good=1.0, bad=0.5)


def test_aggregate_benchmark_empty_observables():
    bench = comparison_module.aggregate_benchmark(
        benchmark_id="fng",
        code_name="openmc",
        code_version="0.14.0",
        reference_source="experiment",
        observables=[],
        include_grading=False,
    )
    assert bench.benchmark_status is None


def test_aggregate_benchmark_with_grading():
    point = PointComparison(
        id="p0",
        observable_type="tally",
        experiment=DataPoint(value=1.0, uncertainty=0.1),
        calculation=DataPoint(value=2.0, uncertainty=0.1),
    )
    point.metrics = compute_point_metrics(point, include_grading=True)

    obs = comparison_module._aggregate_observable(
        "tally",
        "tally",
        [point],
        include_grading=True,
    )

    bench = comparison_module.aggregate_benchmark(
        benchmark_id="fng",
        code_name="openmc",
        code_version="0.14.0",
        reference_source="experiment",
        observables=[obs],
        include_grading=True,
    )
    assert bench.benchmark_status is not None
    assert bench.dashboard_score is not None


def test_datapoints_from_tally_flatten_dims_and_ids():
    mean = xr.DataArray(
        np.arange(6.0).reshape(2, 3),
        dims=("cell", "energy"),
        coords={"cell": [1, 2], "energy": [0.0, 1.0, 2.0]},
    )
    std = xr.DataArray(np.full((2, 3), 0.2), dims=mean.dims, coords=mean.coords)
    tally = BaseTally(mean, std)

    exp_points, calc_points, point_ids = datapoints_from_tally(
        reference=tally,
        candidate=tally,
        flatten_dims=["cell", "energy"],
    )

    assert len(exp_points) == 6
    assert len(calc_points) == 6
    assert point_ids[0].startswith("cell=1")
    assert "energy=0.0" in point_ids[0]


def test_datapoints_from_tally_missing_dim_raises():
    mean = xr.DataArray(np.arange(4.0).reshape(2, 2), dims=("cell", "energy"))
    tally = BaseTally(mean)

    with pytest.raises(ValueError, match="Flatten dims not found"):
        datapoints_from_tally(reference=tally, candidate=tally, flatten_dims=["missing"])


def test_compare_tallies_grading_toggle():
    mean = xr.DataArray(
        np.array([[1.0], [10.0]]).reshape(2, 1, 1),
        dims=("cell", "nuclide", "score"),
        coords={"cell": [1, 2], "nuclide": ["total"], "score": ["flux"]},
    )
    std = xr.DataArray(np.full((2, 1, 1), 0.01), dims=mean.dims, coords=mean.coords)
    exp = BaseTally(mean, std)
    calc = BaseTally(mean * 2.0, std)

    obs = compare_tallies(
        observable_name="tally",
        observable_type="tally",
        reference=exp,
        candidate=calc,
        include_grading=False,
    )
    assert obs.outlier_fraction is None
    assert all(p.metrics is not None and p.metrics.status is None for p in obs.points)

    obs = compare_tallies(
        observable_name="tally",
        observable_type="tally",
        reference=exp,
        candidate=calc,
        include_grading=True,
    )
    assert obs.outlier_fraction is not None
    assert any(p.metrics is not None and p.metrics.status == PointStatus.OUTLIER for p in obs.points)


def test_compare_benchmark_results_defaults_and_grading(tmp_path):
    ref_path = tmp_path / "ref.h5"
    cand_path = tmp_path / "cand.h5"
    _write_results_file(ref_path)
    _write_results_file(cand_path)
    _add_run_metadata(cand_path, code_name="openmc", code_version="0.14.0")

    reference = BenchmarkResults.from_file(ref_path)
    candidate = BenchmarkResults.from_file(cand_path)

    bench = compare_benchmark_results(
        benchmark_id="fng",
        reference_source="experiment",
        reference=reference,
        candidate=candidate,
    )
    assert bench.code_name == "openmc"
    assert bench.code_version == "0.14.0"
    assert bench.benchmark_status is None
    assert bench.dashboard_score is None
    assert bench.outlier_fraction is None

    bench = compare_benchmark_results(
        benchmark_id="fng",
        reference_source="experiment",
        reference=reference,
        candidate=candidate,
        include_grading=True,
    )
    assert bench.benchmark_status in {
        BenchmarkStatus.ACCEPTABLE,
        BenchmarkStatus.BORDERLINE,
        BenchmarkStatus.PROBLEMATIC,
    }
    assert bench.dashboard_score is not None
    assert bench.outlier_fraction is not None


def test_compare_benchmark_results_missing_reference_source(tmp_path):
    ref_path = tmp_path / "ref.h5"
    cand_path = tmp_path / "cand.h5"
    _write_results_file(ref_path)
    _write_results_file(cand_path)
    reference = BenchmarkResults.from_file(ref_path)
    candidate = BenchmarkResults.from_file(cand_path)

    with pytest.raises(ValueError, match="reference_source is required"):
        compare_benchmark_results(
            benchmark_id="fng",
            reference=reference,
            candidate=candidate,
        )


def test_compare_benchmark_results_alias_mismatch_raises(tmp_path):
    ref_path = tmp_path / "ref.h5"
    cand_path = tmp_path / "cand.h5"
    other_path = tmp_path / "other.h5"
    _write_results_file(ref_path)
    _write_results_file(cand_path)
    _write_results_file(other_path)

    reference = BenchmarkResults.from_file(ref_path)
    candidate = BenchmarkResults.from_file(cand_path)
    other = BenchmarkResults.from_file(other_path)

    with pytest.raises(ValueError, match="reference and experiment were provided but differ"):
        compare_benchmark_results(
            benchmark_id="fng",
            reference_source="experiment",
            reference=reference,
            candidate=candidate,
            experiment=other,
        )

    with pytest.raises(ValueError, match="candidate and calculation were provided but differ"):
        compare_benchmark_results(
            benchmark_id="fng",
            reference_source="experiment",
            reference=reference,
            candidate=candidate,
            calculation=other,
        )


def test_compare_benchmark_results_no_tallies_raises(tmp_path):
    ref_path = tmp_path / "ref.h5"
    cand_path = tmp_path / "cand.h5"
    _write_results_file(ref_path)
    _write_results_file(cand_path)

    reference = BenchmarkResults.from_file(ref_path)
    candidate = BenchmarkResults.from_file(cand_path)

    with pytest.raises(ValueError, match="No tallies found"):
        compare_benchmark_results(
            benchmark_id="fng",
            reference_source="experiment",
            reference=reference,
            candidate=candidate,
            tally_names=[],
        )


def test_model_repr_and_relative_uncertainty():
    point = DataPoint(value=0.0, uncertainty=0.1)
    assert point.relative_uncertainty == np.inf

    comparison = PointComparison(
        id="p1",
        observable_type="tally",
        experiment=DataPoint(value=1.0, uncertainty=0.1),
        calculation=DataPoint(value=1.1, uncertainty=0.1),
    )
    assert "uncomputed" in repr(comparison)

    comparison.metrics = compute_point_metrics(comparison, include_grading=False)
    assert "ungraded" in repr(comparison)

    obs = ObservableComparison(name="t1", observable_type="tally", points=[comparison])
    assert "ObservableComparison" in repr(obs)


def test_compare_benchmark_results_alias_conflicts(tmp_path):
    ref_path = tmp_path / "ref.h5"
    cand_path = tmp_path / "cand.h5"
    other_path = tmp_path / "other.h5"
    _write_results_file(ref_path)
    _write_results_file(cand_path)
    _write_results_file(other_path)

    reference = BenchmarkResults.from_file(ref_path)
    candidate = BenchmarkResults.from_file(cand_path)
    other = BenchmarkResults.from_file(other_path)

    with pytest.raises(ValueError, match="reference and experiment were provided but differ"):
        compare_benchmark_results(
            benchmark_id="fng",
            reference_source="experiment",
            reference=reference,
            experiment=other,
            candidate=candidate,
        )

    with pytest.raises(ValueError, match="candidate and calculation were provided but differ"):
        compare_benchmark_results(
            benchmark_id="fng",
            reference_source="experiment",
            reference=reference,
            candidate=candidate,
            calculation=other,
        )
