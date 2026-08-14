import importlib.util
import json
import sys
import types

import h5py

import numpy as np
import pytest
import xarray as xr


def _module_available(name: str) -> bool:
    if name in sys.modules:
        return True
    try:
        return importlib.util.find_spec(name) is not None
    except ValueError:
        return False


# The package __init__ imports benchmark.py, which imports these dependencies.
# Provide minimal stubs in environments where OpenMC stack is not available.
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


from openmc_fusion_benchmarks.benchmark import Benchmark
from openmc_fusion_benchmarks.benchmark_results import BenchmarkResults, OFBResults, Results
from openmc_fusion_benchmarks.tallies import Tally


def _write_structured_group(filepath, group="test_tally", tally_name="mytally"):
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
    ds.attrs["observed_tally"] = json.dumps({"name": tally_name})
    ds.to_netcdf(filepath, mode="w", engine="h5netcdf", group=group)


@pytest.fixture
def temp_results_file(tmp_path):
    filepath = tmp_path / "test_results.h5"
    _write_structured_group(filepath, group="test_tally", tally_name="mytally")
    return filepath


def test_benchmark_results_from_file(temp_results_file):
    results = BenchmarkResults.from_file(temp_results_file)
    assert results.filepath == temp_results_file
    assert results.filepath.exists()


def test_benchmark_results_file_not_found():
    with pytest.raises(FileNotFoundError):
        BenchmarkResults.from_file("/nonexistent/path/results.h5")


def test_benchmark_results_from_run_dir(tmp_path):
    filepath = tmp_path / "results.h5"
    _write_structured_group(filepath, group="tally", tally_name="tally")

    results = BenchmarkResults.from_run_dir(tmp_path, "results.h5")
    assert results.filepath.exists()
    assert results.filepath.name == "results.h5"


def test_results_aliases(temp_results_file):
    base = Results.from_file(temp_results_file)
    alias = OFBResults.from_file(temp_results_file)
    assert isinstance(base, Results)
    assert isinstance(alias, Results)


def test_benchmark_results_tallies_property(temp_results_file):
    results = BenchmarkResults.from_file(temp_results_file)
    assert results.tallies == ["test_tally"]


def test_get_tally_by_group_name(temp_results_file):
    results = BenchmarkResults.from_file(temp_results_file)
    tally = results.get_tally("test_tally")
    assert isinstance(tally, Tally)
    assert tally.name == "mytally"
    assert tally.id == 1


def test_get_tally_by_logical_name(temp_results_file):
    results = BenchmarkResults.from_file(temp_results_file)
    tally = results.get_tally("mytally")
    assert isinstance(tally, Tally)
    assert tally.name == "mytally"


def test_get_tally_missing_mean_raises(tmp_path):
    filepath = tmp_path / "bad_results.h5"
    ds = xr.Dataset({"only_var": xr.DataArray(np.arange(3.0), dims=("row",))})
    ds.to_netcdf(filepath, mode="w", engine="h5netcdf", group="bad")

    results = BenchmarkResults.from_file(filepath)
    with pytest.raises(ValueError, match="does not contain a 'mean' dataset"):
        results.get_tally("bad")


def test_get_tally_unknown_name_raises(temp_results_file):
    results = BenchmarkResults.from_file(temp_results_file)
    with pytest.raises(ValueError, match="No tally with name or group"):
        results.get_tally("missing")


def test_get_spec_consistency_report_and_only_mismatches(tmp_path):
    filepath = tmp_path / "consistency_results.h5"

    ds_ok = xr.Dataset(
        {
            "mean": xr.DataArray(
                np.ones((1, 1), dtype=float),
                dims=("nuclide", "score"),
                coords={
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux"], dtype="U"),
                },
            )
        }
    )
    ds_ok.attrs["spec_consistent"] = 1
    ds_ok.attrs["spec_consistency_issues"] = json.dumps([])
    ds_ok.attrs["observed_tally"] = json.dumps({"name": "ok_name"})
    ds_ok.to_netcdf(filepath, mode="w", engine="h5netcdf", group="ok")

    ds_bad = xr.Dataset(
        {
            "mean": xr.DataArray(
                np.ones((1, 1), dtype=float),
                dims=("nuclide", "score"),
                coords={
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux"], dtype="U"),
                },
            )
        }
    )
    ds_bad.attrs["spec_consistent"] = 0
    ds_bad.attrs["spec_consistency_issues"] = json.dumps(["scores mismatch"])
    ds_bad.attrs["observed_tally"] = json.dumps({"name": "fallback_name"})
    ds_bad.to_netcdf(filepath, mode="a", engine="h5netcdf", group="bad")

    results = BenchmarkResults.from_file(filepath)
    report = results.get_spec_consistency_report()
    mismatches = results.get_spec_consistency_report(only_mismatches=True)

    assert len(report) == 2
    assert len(mismatches) == 1
    assert mismatches[0]["group"] == "bad"
    assert mismatches[0]["tally_name"] == "fallback_name"
    assert mismatches[0]["issues"] == ["scores mismatch"]


def test_results_tallies_skips_non_tally_groups(tmp_path):
    filepath = tmp_path / "mixed_groups.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    with h5py.File(filepath, "a") as handle:
        handle.create_group("run_metadata")
        handle.create_group("specifications")
        handle.create_group("empty_group")

    results = BenchmarkResults.from_file(filepath)
    assert results.tallies == ["tally_1"]


def test_results_tallies_handles_invalid_group(tmp_path):
    filepath = tmp_path / "invalid_group.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    with h5py.File(filepath, "a") as handle:
        handle.create_group("raw_group")

    results = BenchmarkResults.from_file(filepath)
    assert results.tallies == ["tally_1"]


def test_get_spec_consistency_report_handles_empty_and_invalid_attrs(tmp_path):
    filepath = tmp_path / "consistency_empty.h5"

    ds = xr.Dataset(
        {
            "mean": xr.DataArray(
                np.ones((1, 1), dtype=float),
                dims=("nuclide", "score"),
                coords={
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux"], dtype="U"),
                },
            )
        }
    )
    ds.attrs["spec_consistent"] = np.array([], dtype=int)
    ds.attrs["spec_consistency_issues"] = b"not json"
    ds.attrs["observed_tally"] = json.dumps({"name": "fallback"}).encode("utf-8")
    ds.attrs["tally_name"] = ""
    ds.to_netcdf(filepath, mode="w", engine="h5netcdf", group="tally_1")

    results = BenchmarkResults.from_file(filepath)
    report = results.get_spec_consistency_report()

    assert report[0]["spec_consistent"] is None
    assert report[0]["issues"] == []
    assert report[0]["tally_name"] == "fallback"


def test_get_spec_snapshot_and_run_metadata(tmp_path):
    filepath = tmp_path / "metadata_results.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    fake_self = types.SimpleNamespace(
        name="dummy",
        _benchmark_spec={"metadata": {"title": "dummy"}, "tallies": []},
    )

    Benchmark._write_spec_snapshot(fake_self, filename=str(filepath))
    Benchmark._write_run_metadata(
        fake_self,
        code_name="openmc",
        code_version="0.15.2",
        nuclear_data_name=None,
        nuclear_data_version=None,
        geometry="cad",
        filename=str(filepath),
    )

    results = BenchmarkResults.from_file(filepath)
    spec = results.get_spec_snapshot()
    run_meta = results.get_run_metadata()

    assert spec["metadata"]["title"] == "dummy"
    assert run_meta["code_name"] == "openmc"
    assert run_meta["code_version"] == "0.15.2"
    assert run_meta["geometry"] == "cad"


def test_get_spec_snapshot_missing_group_raises(tmp_path):
    filepath = tmp_path / "no_spec.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    results = BenchmarkResults.from_file(filepath)
    with pytest.raises(ValueError, match="No specifications snapshot"):
        results.get_spec_snapshot()


def test_get_spec_snapshot_missing_yaml_dataset_raises(tmp_path):
    filepath = tmp_path / "bad_spec.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    with h5py.File(filepath, "a") as handle:
        if "specifications" in handle:
            del handle["specifications"]
        group = handle.create_group("specifications")
        group.attrs["format"] = "yaml"

    results = BenchmarkResults.from_file(filepath)
    with pytest.raises(ValueError, match="missing 'yaml' dataset"):
        results.get_spec_snapshot()


def test_get_spec_snapshot_invalid_yaml_raises(tmp_path):
    filepath = tmp_path / "invalid_spec.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    with h5py.File(filepath, "a") as handle:
        if "specifications" in handle:
            del handle["specifications"]
        group = handle.create_group("specifications")
        group.create_dataset("yaml", data=np.bytes_(b": not yaml"))

    results = BenchmarkResults.from_file(filepath)
    with pytest.raises(ValueError, match="Failed to parse specifications snapshot"):
        results.get_spec_snapshot()


def test_get_spec_snapshot_empty_dataset_raises(tmp_path):
    filepath = tmp_path / "empty_spec.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    with h5py.File(filepath, "a") as handle:
        if "specifications" in handle:
            del handle["specifications"]
        group = handle.create_group("specifications")
        group.create_dataset("yaml", data=np.array([], dtype="S"))

    results = BenchmarkResults.from_file(filepath)
    with pytest.raises(ValueError, match="dataset is empty"):
        results.get_spec_snapshot()


def test_get_run_metadata_missing_group_raises(tmp_path):
    filepath = tmp_path / "no_meta.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    results = BenchmarkResults.from_file(filepath)
    with pytest.raises(ValueError, match="No run metadata found"):
        results.get_run_metadata()


def test_get_run_metadata_decodes_bytes(tmp_path):
    filepath = tmp_path / "run_meta_bytes.h5"
    _write_structured_group(filepath, group="tally_1", tally_name="tally_1")

    with h5py.File(filepath, "a") as handle:
        group = handle.create_group("run_metadata")
        group.attrs["code_name"] = b"openmc"
        group.attrs["code_version"] = b"0.15.2"

    results = BenchmarkResults.from_file(filepath)
    run_meta = results.get_run_metadata()
    assert run_meta["code_name"] == "openmc"
    assert run_meta["code_version"] == "0.15.2"


def test_get_spec_consistency_report_handles_bytes_and_numpy(tmp_path):
    filepath = tmp_path / "consistency_bytes.h5"

    ds = xr.Dataset(
        {
            "mean": xr.DataArray(
                np.ones((1, 1), dtype=float),
                dims=("nuclide", "score"),
                coords={
                    "nuclide": np.array(["total"], dtype="U"),
                    "score": np.array(["flux"], dtype="U"),
                },
            )
        }
    )
    ds.to_netcdf(filepath, mode="w", engine="h5netcdf", group="tally_1")

    with h5py.File(filepath, "a") as handle:
        attrs = handle["tally_1"].attrs
        attrs["spec_consistent"] = np.array([1], dtype=int)
        attrs["spec_consistency_issues"] = json.dumps(["ok"]).encode("utf-8")
        attrs["observed_tally"] = json.dumps({"name": "from_bytes"}).encode("utf-8")

    results = BenchmarkResults.from_file(filepath)
    report = results.get_spec_consistency_report()
    assert report[0]["spec_consistent"] == 1
    assert report[0]["issues"] == ["ok"]
    assert report[0]["tally_name"] == "from_bytes"
