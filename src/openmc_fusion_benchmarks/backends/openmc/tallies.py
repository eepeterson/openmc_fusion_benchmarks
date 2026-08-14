from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable
import json

import h5py
import numpy as np
import openmc
import xarray as xr

from ...validate_results import normalize_filter_type, validate_tally_consistency


def make_default_openmc_normalizer(mesh: str | Path | object):
    """
    Build a standard OpenMC tally normalizer based on geometry measures.

    The returned callback normalizes each tally over its filter axes using:
    - cell filter bins -> cell volumes from a DAGMC mesh model
    - surface filter bins -> surface areas from a DAGMC mesh model

    Parameters
    ----------
    mesh:
        Either a mesh path accepted by ``pydagmc.Model`` or an existing
        object exposing ``volumes_by_id`` and ``surfaces_by_id`` mappings.
    """
    if hasattr(mesh, "volumes_by_id") and hasattr(mesh, "surfaces_by_id"):
        msh = mesh
    else:
        # Local import keeps this backend module importable even if pydagmc
        # is not installed in non-OpenMC environments.
        import pydagmc

        msh = pydagmc.Model(str(mesh))

    def normalizer(tally: openmc.Tally, mean_nd: np.ndarray, std_nd: np.ndarray):
        filters = list(tally.filters)

        # Filter dimensions are the leading axes in mean_nd/std_nd.
        for axis, flt in enumerate(filters):
            if isinstance(flt, openmc.CellFilter):
                ids = np.asarray(flt.bins, dtype=int).reshape(-1)
                factors = np.asarray([msh.volumes_by_id[i].volume for i in ids], dtype=float)
            elif isinstance(flt, openmc.SurfaceFilter):
                ids = np.asarray(flt.bins, dtype=int).reshape(-1)
                factors = np.asarray([msh.surfaces_by_id[i].area for i in ids], dtype=float)
            elif isinstance(flt, openmc.MaterialFilter):
                raise NotImplementedError(
                    "Material filter normalization is not implemented yet."
                )
            else:
                continue

            if np.any(factors == 0.0):
                raise ValueError("Normalization factor contains zero values.")

            reshape = (1,) * axis + (factors.shape[0],) + (1,) * (mean_nd.ndim - axis - 1)
            norm = factors.reshape(reshape)
            mean_nd = mean_nd / norm
            std_nd = std_nd / norm

        return mean_nd, std_nd

    return normalizer

def _unique_filter_dims(filters: list[openmc.Filter]) -> list[str]:
    """Build stable, unique filter dimension names."""
    counts: dict[str, int] = {}
    dims: list[str] = []
    for flt in filters:
        base = type(flt).__name__.replace("Filter", "").lower() or "filter"
        idx = counts.get(base, 0)
        counts[base] = idx + 1
        dims.append(base if idx == 0 else f"{base}_{idx}")
    return dims


def _serialize_filter_bins(flt: openmc.Filter):
    """Serialize OpenMC filter bins into JSON-compatible Python objects."""
    bins = flt.bins
    arr = np.asarray(bins)

        # Store energy filters as edge lists (E0..En), consistent with specifications.yaml.
    if isinstance(flt, openmc.EnergyFilter):
        if arr.ndim == 1:
            return arr.tolist()
        if arr.ndim == 2 and arr.shape[1] == 2:
            low = arr[:, 0]
            high = arr[:, 1]
            if low.size == 0:
                return []
            return np.concatenate([low[:1], high]).tolist()
        return arr.reshape(-1).tolist()

    # Scalar bins -> simple list of values.
    if arr.ndim <= 1:
        return arr.tolist()

    # Structured bins (for example mesh-like tuples) -> nested lists.
    return [list(np.asarray(b).tolist()) for b in bins]


def _build_filter_axis_metadata(
    filters: list[openmc.Filter],
    filter_dims: list[str],
    spec_tally: dict | None = None,
):
    """Create rich per-filter metadata with axis names and bin definitions."""
    spec_filters = []
    if isinstance(spec_tally, dict):
        spec_filters = list(spec_tally.get("filters", []))

    spec_idx = 0
    axes = []
    for flt, dim in zip(filters, filter_dims):
        ftype = normalize_filter_type(type(flt).__name__)
        spec_filter = None
        if ftype != "particle" and spec_idx < len(spec_filters):
            spec_filter = spec_filters[spec_idx]
            spec_idx += 1

        axis_meta = {
            "name": type(flt).__name__,
            "axis": dim,
            "num_bins": int(flt.num_bins),
            "bins": _serialize_filter_bins(flt),
        }

        if isinstance(spec_filter, dict) and "units" in spec_filter:
            axis_meta["units"] = spec_filter.get("units")

        if ftype == "energy":
            # Energy filters are represented as edge lists in OFB results.
            axis_meta["kind"] = "edges"
            axis_meta["units"] = (
                spec_filter.get("units", "eV")
                if isinstance(spec_filter, dict)
                else "eV"
            )

        axes.append(axis_meta)
    return axes


def _build_spec_lookup(spec_tallies) -> dict[str, dict]:
    """Build a lookup map from tally name to tally specification entry."""
    if spec_tallies is None:
        return {}
    if isinstance(spec_tallies, dict):
        return spec_tallies

    lookup = {}
    for entry in spec_tallies:
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        if name:
            lookup[str(name)] = entry
    return lookup

def _to_1d_coord(values) -> np.ndarray:
    """Normalize scalar/list-like coordinate input into a 1D numpy array."""
    arr = np.asarray(values)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    return arr


def openmc_tally_to_dataset(
    tally: openmc.Tally,
    tmc_coords: dict[str, Iterable] | None = None,
    normalizer: Callable[[openmc.Tally, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]] | None = None,
    spec_tally: dict | None = None,
) -> xr.Dataset:
    """
    Convert one OpenMC tally into an xarray Dataset with variables `mean` and `mc_std`.

    The output schema mirrors the one used in `uq/tmc_manager.py`:
    dimensions are `(tmc dims..., filter dims..., nuclide, score)` and metadata
    is stored in dataset attrs (`filter_axes`, `nuclides`, `scores`).
    """
    tmc_coords = tmc_coords or {}

    filters = list(tally.filters)
    filter_bins = [f.num_bins for f in filters]
    n_nuclides = max(len(tally.nuclides), 1)
    n_scores = len(tally.scores)
    nd_shape = tuple(filter_bins) + (n_nuclides, n_scores)

    mean_nd = tally.mean.reshape(nd_shape)
    std_nd = tally.std_dev.reshape(nd_shape)

    if normalizer is not None:
        mean_nd, std_nd = normalizer(tally, mean_nd, std_nd)

    tmc_dims = tuple(tmc_coords.keys())
    tmc_coord_arrays = {k: _to_1d_coord(v) for k, v in tmc_coords.items()}
    tmc_shape = tuple(len(v) for v in tmc_coord_arrays.values())

    full_shape = tmc_shape + nd_shape
    mean_full = mean_nd.reshape((1,) * len(tmc_shape) + nd_shape)
    std_full = std_nd.reshape((1,) * len(tmc_shape) + nd_shape)

    # Broadcast to provided TMC coordinate lengths (usually all 1 for a single write).
    mean_full = np.broadcast_to(mean_full, full_shape)
    std_full = np.broadcast_to(std_full, full_shape)

    filter_dims = _unique_filter_dims(filters)
    dims = tmc_dims + tuple(filter_dims) + ("nuclide", "score")

    coords: dict[str, tuple[str, np.ndarray] | np.ndarray] = {}
    for d, c in tmc_coord_arrays.items():
        coords[d] = (d, c)
    for f, d in zip(filters, filter_dims):
        coords[d] = (d, np.arange(f.num_bins))

    nuclides = [str(n) for n in tally.nuclides] if tally.nuclides else ["total"]
    scores = [str(s) for s in tally.scores]
    coords["nuclide"] = ("nuclide", np.asarray(nuclides, dtype="U"))
    coords["score"] = ("score", np.asarray(scores, dtype="U"))

    ds = xr.Dataset(
        {
            "mean": xr.DataArray(mean_full, dims=dims, coords=coords),
            "mc_std": xr.DataArray(std_full, dims=dims, coords=coords),
        }
    )

    tally_name = tally.name or f"tally_{tally.id}"
    ds["mean"].attrs["tally_id"] = int(tally.id)
    ds["mean"].attrs["tally_name"] = tally_name
    ds["mc_std"].attrs["tally_id"] = int(tally.id)
    ds["mc_std"].attrs["tally_name"] = tally_name

    filter_axes = _build_filter_axis_metadata(filters, filter_dims, spec_tally=spec_tally)
    
    ds.attrs["filter_axes"] = json.dumps(filter_axes)
    ds.attrs["nuclides"] = json.dumps(nuclides)
    ds.attrs["scores"] = json.dumps(scores)

    observed_tally = {
        "name": tally_name,
        "id": int(tally.id),
        "filters": filter_axes,
        "scores": scores,
        "nuclides": nuclides,
    }
    ds.attrs["observed_tally"] = json.dumps(observed_tally)

    if spec_tally is not None:
        ds.attrs["spec_tally"] = json.dumps(spec_tally)
        consistent, issues = validate_tally_consistency(spec_tally, observed_tally)
        # h5netcdf/netCDF attrs do not support boolean dtype reliably.
        ds.attrs["spec_consistent"] = int(bool(consistent))
        ds.attrs["spec_consistency_issues"] = json.dumps(issues)

    return ds

def save_openmc_statepoint_tallies(
    statepoint: openmc.StatePoint,
    filename: str | Path,
    tally_names: Iterable[str] | None = None,
    spec_tallies: Iterable[dict] | dict[str, dict] | None = None,
    tmc_coords: dict[str, Iterable] | None = None,
    append_dim: str | None = None,
    normalizer: Callable[[openmc.Tally, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]] | None = None,
    group_by: str = "name",
    engine: str = "h5netcdf",
) -> Path:
    """
    Write OpenMC statepoint tallies to grouped HDF5/NetCDF datasets.


    Parameters
    ----------
    group_by:
        Group naming strategy, either:
        - ``"name"``: use tally names (recommended for benchmark workflows)
        - ``"id"``: use ``tally_<id>`` groups
    """
    if group_by not in {"name", "id"}:
        raise ValueError("group_by must be either 'name' or 'id'")
    
    filename = Path(filename)
    spec_lookup = _build_spec_lookup(spec_tallies)
    if tally_names is None:
        selected = list(statepoint.tallies.values())
    else:
        selected = [statepoint.get_tally(name=name) for name in tally_names]

    for tally in selected:
        ds_new = openmc_tally_to_dataset(
            tally=tally,
            tmc_coords=tmc_coords,
            normalizer=normalizer,
            spec_tally=spec_lookup.get(str(tally.name)),
        )
        if group_by == "name" and tally.name:
            group = str(tally.name)
        else:
            group = f"tally_{int(tally.id)}"

        # Keep explicit group metadata alongside tally_id / tally_name metadata.
        ds_new.attrs["group"] = group
        ds_new["mean"].attrs["tally_group"] = group
        ds_new["mc_std"].attrs["tally_group"] = group

        if not filename.exists():
            ds_new.to_netcdf(filename, mode="w", group=group, engine=engine)
            continue

        try:
            ds_old = xr.open_dataset(filename, group=group, engine=engine)
            if append_dim is not None and append_dim in ds_old.dims and append_dim in ds_new.dims:
                ds_combined = xr.concat([ds_old, ds_new], dim=append_dim)
            else:
                ds_combined = ds_new
            ds_old.close()
        except (OSError, ValueError, KeyError):
            ds_combined = ds_new

        with h5py.File(filename, "a") as h5f:
            if group in h5f:
                del h5f[group]
        ds_combined.to_netcdf(filename, mode="a", group=group, engine=engine)

    return filename.resolve()


# # Backward-compatible aliases while call sites migrate.
# tally_to_dataset = openmc_tally_to_dataset
# save_statepoint_tallies = save_openmc_statepoint_tallies