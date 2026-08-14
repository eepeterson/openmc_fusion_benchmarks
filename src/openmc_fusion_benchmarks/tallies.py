from __future__ import annotations

import json

import numpy as np
import xarray as xr


class BaseTally:
	"""Common wrapper for structured tally datasets."""

	def __init__(self, mean_da: xr.DataArray, mc_std_da: xr.DataArray | None = None, parent_ds: xr.Dataset | None = None):
		self._da = mean_da
		self._da_mc_std = mc_std_da
		self._parent_ds = parent_ds

	@property
	def id(self):
		"""Tally ID."""
		return self._da.attrs.get("tally_id")

	@property
	def name(self):
		"""Tally name."""
		return self._da.attrs.get("tally_name")

	@property
	def data(self):
		"""Underlying mean data array values."""
		return self._da.values

	@property
	def scores(self):
		"""List of score names."""
		if "score" in self._da.coords:
			return [str(s) for s in self._da.coords["score"].values]

		if self._parent_ds is not None:
			scores_json = self._parent_ds.attrs.get("scores")
			if scores_json:
				return json.loads(scores_json)

		scores_json = self._da.attrs.get("scores")
		if scores_json:
			return json.loads(scores_json)
		return []

	@property
	def nuclides(self):
		"""List of nuclide names."""
		if "nuclide" in self._da.coords:
			return [str(n) for n in self._da.coords["nuclide"].values]

		if self._parent_ds is not None:
			nuclides_json = self._parent_ds.attrs.get("nuclides")
			if nuclides_json:
				return json.loads(nuclides_json)

		nuclides_json = self._da.attrs.get("nuclides")
		if nuclides_json:
			return json.loads(nuclides_json)
		return []

	@property
	def filters(self):
		"""List of filter metadata entries."""
		if self._parent_ds is not None:
			filt_json = self._parent_ds.attrs.get("filter_axes")
			if filt_json:
				return json.loads(filt_json)

		filt_json = self._da.attrs.get("filter_axes")
		if filt_json:
			return json.loads(filt_json)
		return []

	@property
	def shape(self):
		"""Shape of the underlying mean data array."""
		return self._da.shape

	@property
	def dims(self):
		"""Dimension names of the underlying mean data array."""
		return self._da.dims

	@property
	def mean(self):
		"""Mean array for non-TMC tallies (already the final mean)."""
		return self._da.values

	@property
	def std_dev(self):
		"""MC standard deviation array for non-TMC tallies."""
		if self._da_mc_std is not None:
			return self._da_mc_std.values
		return np.zeros_like(self._da.values)
	
	def dimension_report(self):
		"""Return a compact mapping report for OFB and OpenMC-equivalent shapes."""
		dims = list(self._da.dims)
		ofb_shape = tuple(int(self._da.sizes[d]) for d in dims)

		filter_axes = self.filters if isinstance(self.filters, list) else []
		filter_dims = [
			str(a.get("axis"))
			for a in filter_axes
			if isinstance(a, dict) and a.get("axis") in self._da.dims
		]

		if not filter_dims:
			known_tmc_dims = {"realization", "sample", "replica", "batch", "iteration", "case"}
			filter_dims = [d for d in dims if d not in known_tmc_dims and d not in {"nuclide", "score"}]

		filter_dim_sizes = {d: int(self._da.sizes[d]) for d in filter_dims}
		flat_filter_bins = int(np.prod(list(filter_dim_sizes.values()))) if filter_dim_sizes else 1

		nuclide_size = int(self._da.sizes["nuclide"]) if "nuclide" in self._da.sizes else max(len(self.nuclides), 1)
		score_size = int(self._da.sizes["score"]) if "score" in self._da.sizes else max(len(self.scores), 1)

		tmc_dims = [d for d in dims if d not in set(filter_dims) and d not in {"nuclide", "score"}]
		tmc_sizes = {d: int(self._da.sizes[d]) for d in tmc_dims}

		return {
			"tally_name": self.name,
			"tally_id": self.id,
			"ofb_dims": dims,
			"ofb_shape": ofb_shape,
			"tmc_dims": tmc_dims,
			"tmc_sizes": tmc_sizes,
			"filter_dims": filter_dims,
			"filter_dim_sizes": filter_dim_sizes,
			"nuclide_size": nuclide_size,
			"score_size": score_size,
			"openmc_equivalent_raw_shape": (flat_filter_bins, nuclide_size, score_size),
		}

	def format_dimension_report(self):
		"""Return a readable text report from ``dimension_report``."""
		report = self.dimension_report()

		def _fmt_axis_sizes(axis_names, axis_sizes):
			if not axis_names:
				return "none"
			return ", ".join(f"{name}={axis_sizes.get(name, '?')}" for name in axis_names)

		tally_name = report.get("tally_name") or "<unnamed>"
		tally_id = report.get("tally_id")
		ofb_dims = report.get("ofb_dims", [])
		ofb_shape = report.get("ofb_shape", ())
		tmc_dims = report.get("tmc_dims", [])
		tmc_sizes = report.get("tmc_sizes", {})
		filter_dims = report.get("filter_dims", [])
		filter_sizes = report.get("filter_dim_sizes", {})
		raw_shape = report.get("openmc_equivalent_raw_shape", ())

		label_w = 20
		lines = [
			f"Tally {tally_name} (id={tally_id})",
			f"{'OFB dims':<{label_w}}: {ofb_dims}",
			f"{'OFB shape':<{label_w}}: {ofb_shape}",
			f"{'TMC axes':<{label_w}}: {_fmt_axis_sizes(tmc_dims, tmc_sizes)}",
			f"{'Filter axes':<{label_w}}: {_fmt_axis_sizes(filter_dims, filter_sizes)}",
			f"{'Nuclide/score sizes':<{label_w}}: nuclide={report.get('nuclide_size')}, score={report.get('score_size')}",
			f"{'OpenMC raw equivalent':<{label_w}}: {raw_shape}  (flat_filters, nuclide, score)",
		]
		return "\n".join(lines)


	def get_slice(self, scores=None, nuclides=None, **filter_kwargs):
		"""Get a filtered xarray view of the mean tally data."""
		da = self._da

		if filter_kwargs:
			da = da.sel(**filter_kwargs)

		if scores is not None and "score" in da.dims:
			all_scores = self.scores
			score_indices = [all_scores.index(s) for s in scores]
			da = da.isel(score=score_indices)

		if nuclides is not None and "nuclide" in da.dims:
			all_nuclides = self.nuclides
			nuclide_indices = [all_nuclides.index(n) for n in nuclides]
			da = da.isel(nuclide=nuclide_indices)

		return da

class Tally(BaseTally):
	"""Generic OFB tally wrapper for non-TMC result groups."""

	def __repr__(self):
		return f"<Tally {self.id}: '{self.name}', shape={self.shape}>"

def tally_to_dataset(*args, **kwargs):
	"""Backward-compatible alias to the OpenMC backend serializer."""
	from .backends.openmc.tallies import openmc_tally_to_dataset

	return openmc_tally_to_dataset(*args, **kwargs)


def save_statepoint_tallies(*args, **kwargs):
	"""Backward-compatible alias to the OpenMC backend statepoint writer."""
	from .backends.openmc.tallies import save_openmc_statepoint_tallies

	return save_openmc_statepoint_tallies(*args, **kwargs)
