from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

openmc = pytest.importorskip("openmc")

from openmc_fusion_benchmarks.backends.openmc import tallies as backend


def _make_filter(filter_cls, bins, num_bins=None):
    """Create filter objects for both real OpenMC and lightweight stubs."""
    try:
        flt = filter_cls(bins)
        # Real OpenMC objects expose read-only num_bins; do not overwrite.
        return flt
    except TypeError:
        flt = filter_cls()
        flt.bins = bins
        flt.num_bins = int(num_bins if num_bins is not None else len(np.asarray(bins).reshape(-1)))
        return flt


class DummyStatePoint:
    def __init__(self, tallies):
        self._selected = list(tallies)
        self.tallies = {i + 1: t for i, t in enumerate(self._selected)}

    def get_tally(self, name):
        for tally in self._selected:
            if tally.name == name:
                return tally
        raise ValueError(name)


class DummyTally:
    def __init__(self, tally_id, name, filters, nuclides, scores, mean_nd, std_nd):
        self.id = tally_id
        self.name = name
        self.filters = filters
        self.nuclides = nuclides
        self.scores = scores
        self.mean = np.asarray(mean_nd).reshape(-1)
        self.std_dev = np.asarray(std_nd).reshape(-1)


def test_unique_filter_dims_and_coord_helper():
    filters = [
        _make_filter(openmc.CellFilter, [1], num_bins=1),
        _make_filter(openmc.CellFilter, [2], num_bins=1),
        _make_filter(openmc.EnergyFilter, [0.0, 1.0], num_bins=1),
    ]
    assert backend._unique_filter_dims(filters) == ["cell", "cell_1", "energy"]

    scalar = backend._to_1d_coord(7)
    assert scalar.shape == (1,)
    assert scalar[0] == 7


def test_make_default_openmc_normalizer_cell_and_surface():
    mesh = SimpleNamespace(
        volumes_by_id={1: SimpleNamespace(volume=2.0), 2: SimpleNamespace(volume=4.0)},
        surfaces_by_id={10: SimpleNamespace(area=5.0)},
    )
    normalizer = backend.make_default_openmc_normalizer(mesh)

    tally = DummyTally(
        tally_id=1,
        name="norm",
        filters=[
            _make_filter(openmc.CellFilter, [1, 2], num_bins=2),
            _make_filter(openmc.SurfaceFilter, [10], num_bins=1),
        ],
        nuclides=["total"],
        scores=["flux"],
        mean_nd=np.ones((2, 1, 1)),
        std_nd=np.ones((2, 1, 1)),
    )

    mean, std = normalizer(tally, np.ones((2, 1, 1)), np.ones((2, 1, 1)))
    np.testing.assert_allclose(mean[:, 0, 0], np.array([1.0 / 10.0, 1.0 / 20.0]))
    np.testing.assert_allclose(std[:, 0, 0], np.array([1.0 / 10.0, 1.0 / 20.0]))


def test_make_default_openmc_normalizer_material_filter_raises():
    mesh = SimpleNamespace(volumes_by_id={}, surfaces_by_id={})
    normalizer = backend.make_default_openmc_normalizer(mesh)
    tally = DummyTally(
        tally_id=1,
        name="mat",
        filters=[_make_filter(openmc.MaterialFilter, [1], num_bins=1)],
        nuclides=["total"],
        scores=["flux"],
        mean_nd=np.ones((1, 1, 1)),
        std_nd=np.ones((1, 1, 1)),
    )

    with pytest.raises(NotImplementedError, match="Material filter normalization"):
        normalizer(tally, np.ones((1, 1, 1)), np.ones((1, 1, 1)))


def test_make_default_openmc_normalizer_zero_factor_raises():
    mesh = SimpleNamespace(
        volumes_by_id={1: SimpleNamespace(volume=0.0)},
        surfaces_by_id={},
    )
    normalizer = backend.make_default_openmc_normalizer(mesh)
    tally = DummyTally(
        tally_id=1,
        name="zero",
        filters=[_make_filter(openmc.CellFilter, [1], num_bins=1)],
        nuclides=["total"],
        scores=["flux"],
        mean_nd=np.ones((1, 1, 1)),
        std_nd=np.ones((1, 1, 1)),
    )

    with pytest.raises(ValueError, match="Normalization factor contains zero"):
        normalizer(tally, np.ones((1, 1, 1)), np.ones((1, 1, 1)))


def test_openmc_tally_to_dataset_and_save_append(tmp_path):
    filters = [
        _make_filter(openmc.ParticleFilter, ["neutron"], num_bins=1),
        _make_filter(openmc.CellFilter, [1, 2], num_bins=2),
    ]
    mean_nd = np.arange(4.0).reshape(1, 2, 1, 2)
    std_nd = np.full_like(mean_nd, 0.2)
    tally = DummyTally(
        tally_id=7,
        name="cell_flux",
        filters=filters,
        nuclides=["total"],
        scores=["flux", "heating"],
        mean_nd=mean_nd,
        std_nd=std_nd,
    )

    ds = backend.openmc_tally_to_dataset(
        tally=tally,
        tmc_coords={"realization": ["r0"]},
        normalizer=lambda _t, m, s: (m * 2.0, s * 2.0),
    )
    assert tuple(ds["mean"].dims) == ("realization", "particle", "cell", "nuclide", "score")
    assert ds["mean"].attrs["tally_name"] == "cell_flux"

    sp = DummyStatePoint([tally])
    fpath = tmp_path / "results.h5"
    backend.save_openmc_statepoint_tallies(
        statepoint=sp,
        filename=fpath,
        tally_names=["cell_flux"],
        tmc_coords={"realization": ["r0"]},
        append_dim="realization",
    )
    backend.save_openmc_statepoint_tallies(
        statepoint=sp,
        filename=fpath,
        tally_names=["cell_flux"],
        tmc_coords={"realization": ["r1"]},
        append_dim="realization",
    )

    loaded = xr.open_dataset(fpath, group="cell_flux", engine="h5netcdf")
    assert int(loaded.sizes["realization"]) == 2
    loaded.close()


def test_save_openmc_statepoint_tallies_group_by_id_and_validation(tmp_path):
    tally = DummyTally(
        tally_id=42,
        name="",
        filters=[_make_filter(openmc.CellFilter, [1], num_bins=1)],
        nuclides=["total"],
        scores=["flux"],
        mean_nd=np.ones((1, 1, 1)),
        std_nd=np.ones((1, 1, 1)),
    )
    sp = DummyStatePoint([tally])
    fpath = tmp_path / "id_results.h5"

    out = backend.save_openmc_statepoint_tallies(statepoint=sp, filename=fpath, group_by="id")
    assert Path(out).exists()

    loaded = xr.open_dataset(fpath, group="tally_42", engine="h5netcdf")
    assert loaded.attrs["group"] == "tally_42"
    loaded.close()

    with pytest.raises(ValueError, match="group_by"):
        backend.save_openmc_statepoint_tallies(statepoint=sp, filename=fpath, group_by="bad")


def test_build_spec_lookup_variants():
    assert backend._build_spec_lookup(None) == {}
    assert backend._build_spec_lookup({"a": {"name": "a"}}) == {"a": {"name": "a"}}

    lookup = backend._build_spec_lookup([
        {"name": "first", "scores": ["flux"]},
        "ignored",
        {"name": "second", "scores": ["heating"]},
    ])
    assert set(lookup.keys()) == {"first", "second"}


def test_openmc_tally_to_dataset_with_spec_consistency_attrs():
    filters = [
        _make_filter(openmc.ParticleFilter, ["neutron"], num_bins=1),
        _make_filter(openmc.EnergyFilter, [0.0, 1.0e6, 14.0e6], num_bins=2),
    ]
    tally = DummyTally(
        tally_id=9,
        name="energy_flux",
        filters=filters,
        nuclides=["total"],
        scores=["flux"],
        mean_nd=np.ones((1, 2, 1, 1)),
        std_nd=np.full((1, 2, 1, 1), 0.1),
    )

    spec_tally = {
        "name": "energy_flux",
        "particle": "neutron",
        "filters": [{"type": "energy", "values": [0.0, 1.0e6, 14.0e6], "units": "eV"}],
        "scores": ["flux"],
        "nuclides": ["total"],
    }

    ds = backend.openmc_tally_to_dataset(tally=tally, spec_tally=spec_tally)

    assert ds.attrs["spec_consistent"] == 1
    assert "spec_tally" in ds.attrs
    assert "spec_consistency_issues" in ds.attrs

    filter_axes = __import__("json").loads(ds.attrs["filter_axes"])
    energy_axes = [a for a in filter_axes if a["name"] == "EnergyFilter"]
    assert len(energy_axes) == 1
    assert energy_axes[0]["bins"] == [0.0, 1.0e6, 14.0e6]
    assert energy_axes[0]["units"] == "eV"
    assert energy_axes[0]["kind"] == "edges"


def test_openmc_tally_to_dataset_with_spec_mismatch_sets_issues():
    filters = [
        _make_filter(openmc.ParticleFilter, ["neutron"], num_bins=1),
        _make_filter(openmc.CellFilter, [1], num_bins=1),
    ]
    tally = DummyTally(
        tally_id=10,
        name="cell_flux",
        filters=filters,
        nuclides=["total"],
        scores=["flux"],
        mean_nd=np.ones((1, 1, 1, 1)),
        std_nd=np.full((1, 1, 1, 1), 0.1),
    )

    # Intentionally mismatch score name.
    spec_tally = {
        "name": "cell_flux",
        "particle": "neutron",
        "filters": [{"type": "cell", "values": [1]}],
        "scores": ["heating"],
        "nuclides": ["total"],
    }

    ds = backend.openmc_tally_to_dataset(tally=tally, spec_tally=spec_tally)
    assert ds.attrs["spec_consistent"] == 0
    issues = __import__("json").loads(ds.attrs["spec_consistency_issues"])
    assert any("scores mismatch" in issue for issue in issues)
