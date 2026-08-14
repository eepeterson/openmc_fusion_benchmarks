import importlib.util
import json
import sys
import types

import numpy as np
import xarray as xr


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


from openmc_fusion_benchmarks.tallies import BaseTally, Tally


def test_base_tally_properties_and_slice_from_coords():
    mean = xr.DataArray(
        np.arange(8.0).reshape(2, 2, 1, 2),
        dims=("realization", "cell", "nuclide", "score"),
        coords={
            "realization": ["r0", "r1"],
            "cell": [10, 20],
            "nuclide": np.array(["total"], dtype="U"),
            "score": np.array(["flux", "heating"], dtype="U"),
        },
    )
    mean.attrs["tally_id"] = 11
    mean.attrs["tally_name"] = "my_tally"

    std = xr.DataArray(np.full(mean.shape, 0.2), dims=mean.dims, coords=mean.coords)
    tally = BaseTally(mean, std)

    assert tally.id == 11
    assert tally.name == "my_tally"
    assert tally.shape == (2, 2, 1, 2)
    assert tally.dims == ("realization", "cell", "nuclide", "score")
    assert tally.scores == ["flux", "heating"]
    assert tally.nuclides == ["total"]

    sliced = tally.get_slice(scores=["heating"], nuclides=["total"], cell=20)
    assert sliced.shape == (2, 1, 1)


def test_base_tally_uses_parent_dataset_attrs_and_std_default():
    mean = xr.DataArray(np.ones((1, 1)), dims=("i", "j"))
    mean.attrs["scores"] = json.dumps(["score_attr"])
    mean.attrs["nuclides"] = json.dumps(["n_attr"])
    mean.attrs["filter_axes"] = json.dumps([{"name": "CellFilter", "num_bins": 1}])

    parent = xr.Dataset()
    parent.attrs["scores"] = json.dumps(["score_parent"])
    parent.attrs["nuclides"] = json.dumps(["n_parent"])
    parent.attrs["filter_axes"] = json.dumps([{"name": "SurfaceFilter", "num_bins": 2}])

    tally = BaseTally(mean, parent_ds=parent)
    assert tally.scores == ["score_parent"]
    assert tally.nuclides == ["n_parent"]
    assert tally.filters == [{"name": "SurfaceFilter", "num_bins": 2}]
    np.testing.assert_allclose(tally.std_dev, np.zeros_like(mean.values))


def test_tally_repr():
    mean = xr.DataArray(np.ones((1,)), dims=("i",))
    mean.attrs["tally_id"] = 3
    mean.attrs["tally_name"] = "repr_tally"
    tally = Tally(mean)
    assert "repr_tally" in repr(tally)


def test_dimension_report_with_filter_axes_metadata():
    mean = xr.DataArray(
        np.ones((2, 3, 4, 1, 2), dtype=float),
        dims=("realization", "cell", "energy", "nuclide", "score"),
        coords={
            "realization": ["r0", "r1"],
            "cell": [1, 2, 3],
            "energy": [0, 1, 2, 3],
            "nuclide": np.array(["total"], dtype="U"),
            "score": np.array(["flux", "heating"], dtype="U"),
        },
    )
    mean.attrs["tally_id"] = 21
    mean.attrs["tally_name"] = "energy_current"

    parent = xr.Dataset()
    parent.attrs["filter_axes"] = json.dumps(
        [
            {"name": "CellFilter", "axis": "cell", "num_bins": 3, "bins": [1, 2, 3]},
            {"name": "EnergyFilter", "axis": "energy", "num_bins": 4, "bins": [0.0, 1.0, 2.0, 3.0, 4.0]},
        ]
    )

    tally = BaseTally(mean, parent_ds=parent)
    report = tally.dimension_report()

    assert report["tally_id"] == 21
    assert report["tally_name"] == "energy_current"
    assert report["filter_dims"] == ["cell", "energy"]
    assert report["tmc_dims"] == ["realization"]
    assert report["filter_dim_sizes"] == {"cell": 3, "energy": 4}
    assert report["openmc_equivalent_raw_shape"] == (12, 1, 2)


def test_dimension_report_fallback_without_filter_axes():
    mean = xr.DataArray(
        np.ones((2, 5, 1, 1), dtype=float),
        dims=("sample", "mesh", "nuclide", "score"),
        coords={
            "sample": [0, 1],
            "mesh": [0, 1, 2, 3, 4],
            "nuclide": np.array(["total"], dtype="U"),
            "score": np.array(["flux"], dtype="U"),
        },
    )
    mean.attrs["tally_name"] = "fallback"

    tally = BaseTally(mean)
    report = tally.dimension_report()

    # "sample" is treated as TMC; "mesh" should be inferred as filter axis.
    assert report["tmc_dims"] == ["sample"]
    assert report["filter_dims"] == ["mesh"]
    assert report["openmc_equivalent_raw_shape"] == (5, 1, 1)


def test_format_dimension_report_includes_sections():
    mean = xr.DataArray(
        np.ones((1, 1), dtype=float),
        dims=("nuclide", "score"),
        coords={
            "nuclide": np.array(["total"], dtype="U"),
            "score": np.array(["flux"], dtype="U"),
        },
    )
    mean.attrs["tally_id"] = 99
    mean.attrs["tally_name"] = "mini"

    text = BaseTally(mean).format_dimension_report()
    assert "Tally mini (id=99)" in text
    assert "OFB dims" in text
    assert "TMC axes" in text
    assert "Filter axes" in text
    assert "OpenMC raw equivalent" in text
