from __future__ import annotations

import json
import yaml
from pathlib import Path
from typing import Union
import h5py
import numpy as np
import xarray as xr

from .database import _resolve_database_path

from .tallies import Tally


class Results:
    """Class to handle OFB results stored in HDF5 files.

    This is the generic entry point and can be used for benchmark and non-benchmark
    workflows as long as they follow the OFB tally group schema.
        - Load from arbitrary file path: `Results.from_file(filepath)`
        - Load from run directory: `Results.from_run_dir(run_dir, filename)`
        - Load from package database: `Results.from_database(benchmark, filename)`
    """

    def __init__(self, filepath: Union[str, Path]):
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(f"Results file not found: {self.filepath}")

    @classmethod
    def from_file(cls, filepath: Union[str, Path]) -> "Results":
        """Load results from an arbitrary file path."""
        return cls(filepath)

    @classmethod
    def from_run_dir(
        cls,
        run_dir: Union[str, Path] = ".",
        filename: str = "benchmark_results.h5",
    ) -> "Results":
        """Load results from a run directory (default: current directory)."""
        path = Path(run_dir) / filename
        return cls(path)

    @classmethod
    def from_database(cls, benchmark: str, filename: str = "reference_results.h5") -> "Results":
        """Load reference results from the package database."""
        db_path = _resolve_database_path(benchmark, filename)
        return cls(db_path)

    @property
    def tallies(self):
        tallies: list[str] = []
        with h5py.File(self.filepath, "r") as f:
            for group in f.keys():
                try:
                    with xr.open_dataset(self.filepath, group=group, engine="h5netcdf") as ds:
                        if "mean" in ds:
                            tallies.append(group)
                except (OSError, ValueError, KeyError):
                    continue
        return tallies


    def get_tally(self, name: str) -> Tally:
        """
        Get a tally wrapper.

        `name` can be either:
        - an HDF5 group name (for example `tally_1`), or
        - a tally logical name stored in attrs (for example `nuclear_heating`).
        """
        group = None

        with h5py.File(self.filepath, "r") as f:
            if name in f:
                group = name
            else:
                # Fallback lookup by tally_name attribute.
                for candidate in f.keys():
                    try:
                        with xr.open_dataset(self.filepath, group=candidate, engine="h5netcdf") as ds:
                            if "mean" not in ds:
                                continue
                            if ds["mean"].attrs.get("tally_name") == name:
                                group = candidate
                                break
                    except (OSError, ValueError, KeyError):
                        continue

        if group is None:
            raise ValueError(f"No tally with name or group '{name}' found")

        ds = xr.open_dataset(self.filepath, group=group, engine="h5netcdf")
        if "mean" not in ds:
            ds.close()
            raise ValueError(f"Group '{group}' does not contain a 'mean' dataset")

        da_mean = ds["mean"]
        da_mc_std = ds["mc_std"] if "mc_std" in ds else None
        return Tally(da_mean, da_mc_std, parent_ds=ds)
    
    def get_spec_consistency_report(self, only_mismatches: bool = False):
        """
        Return spec-vs-observed consistency metadata for each tally group.

        Parameters
        ----------
        only_mismatches:
            If True, return only entries with `spec_consistent == 0`.

        Returns
        -------
        list[dict]
            One entry per tally group with fields:
            - `group`
            - `tally_name`
            - `spec_consistent` (0/1/None)
            - `issues` (list[str])
            - `spec_tally` (dict)
            - `observed_tally` (dict)
        """

        def _loads_attr(attrs, key, default):
            value = attrs.get(key)
            if value is None:
                return default
            if isinstance(value, bytes):
                value = value.decode("utf-8")
            try:
                return json.loads(value)
            except Exception:
                return default

        report = []
        with h5py.File(self.filepath, "r") as f:
            for group in f.keys():
                attrs = f[group].attrs
                spec_consistent = attrs.get("spec_consistent")
                if spec_consistent is not None:
                    try:
                        # h5 attrs may come back as numpy scalars/arrays depending on backend.
                        if isinstance(spec_consistent, np.ndarray):
                            if spec_consistent.size == 0:
                                spec_consistent = None
                            else:
                                spec_consistent = int(np.ravel(spec_consistent)[0])
                        else:
                            spec_consistent = int(spec_consistent)
                    except Exception:
                        spec_consistent = None

                entry = {
                    "group": group,
                    "tally_name": attrs.get("tally_name"),
                    "spec_consistent": spec_consistent,
                    "issues": _loads_attr(attrs, "spec_consistency_issues", []),
                    "spec_tally": _loads_attr(attrs, "spec_tally", {}),
                    "observed_tally": _loads_attr(attrs, "observed_tally", {}),
                }

                # Fallback tally_name from nested metadata when variable-level attrs were used.
                if entry["tally_name"] in (None, ""):
                    entry["tally_name"] = entry["observed_tally"].get("name")

                if only_mismatches and entry["spec_consistent"] != 0:
                    continue

                report.append(entry)

        return report

    def get_spec_snapshot(self) -> dict:
        """Return the embedded specifications.yaml snapshot if present."""
        with h5py.File(self.filepath, "r") as f:
            if "specifications" not in f:
                raise ValueError("No specifications snapshot found in results file")
            group = f["specifications"]
            if "yaml" not in group:
                raise ValueError("Specifications snapshot is missing 'yaml' dataset")
            raw = group["yaml"][()]
            if isinstance(raw, np.ndarray):
                if raw.size == 0:
                    raise ValueError("Specifications snapshot dataset is empty")
                raw = raw[()]
            if isinstance(raw, bytes):
                yaml_text = raw.decode("utf-8")
            else:
                yaml_text = str(raw)
        try:
            return json.loads(json.dumps(yaml.safe_load(yaml_text)))
        except Exception as exc:
            raise ValueError("Failed to parse specifications snapshot YAML") from exc

    def get_run_metadata(self) -> dict:
        """Return run metadata embedded in the results file, if present."""
        with h5py.File(self.filepath, "r") as f:
            if "run_metadata" not in f:
                raise ValueError("No run metadata found in results file")
            attrs = f["run_metadata"].attrs
            data: dict[str, str] = {}
            for key, value in attrs.items():
                if isinstance(value, bytes):
                    data[key] = value.decode("utf-8")
                else:
                    data[key] = str(value)
            return data


class BenchmarkResults(Results):
    """Backward-compatible alias class for benchmark-centric naming."""


# Explicit alias for OFB naming style.
OFBResults = Results