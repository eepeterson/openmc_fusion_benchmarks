from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import xarray as xr

from ..benchmark_results import BenchmarkResults
from ..tallies import BaseTally
from .comparison import aggregate_benchmark, compare_point_set
from .model import BenchmarkComparison, DataPoint, ObservableComparison


def _normalize_flatten_dims(tally: BaseTally, flatten_dims: Sequence[str] | None) -> list[str]:
    if flatten_dims is None:
        return list(tally.dims)

    dims = list(tally.dims)
    missing = [d for d in flatten_dims if d not in dims]
    if missing:
        raise ValueError(f"Flatten dims not found in tally dims: {missing}")
    return list(flatten_dims)


def _format_point_ids(dims: Sequence[str], index_values: Iterable) -> list[str]:
    ids: list[str] = []
    for entry in index_values:
        values = entry if isinstance(entry, tuple) else (entry,)
        parts = [f"{dim}={value}" for dim, value in zip(dims, values)]
        ids.append("|".join(parts))
    return ids


def _stack_dataarray(da: xr.DataArray, dims: Sequence[str]) -> xr.DataArray:
    if list(da.dims) == list(dims):
        return da.stack(point=dims)
    return da.transpose(*dims).stack(point=dims)


def datapoints_from_tally(
    reference: BaseTally | None = None,
    candidate: BaseTally | None = None,
    *,
    experiment: BaseTally | None = None,
    calculation: BaseTally | None = None,
    flatten_dims: Sequence[str] | None = None,
) -> tuple[list[DataPoint], list[DataPoint], list[str]]:
    """
    Convert two BaseTally objects into DataPoint lists and point IDs.

    Parameters
    ----------
    reference:
        Reference results tally.
    candidate:
        Calculated results tally.
    flatten_dims:
        Dimensions to flatten into points. When omitted, all tally dims are used.
    """
    if reference is None:
        reference = experiment
    if candidate is None:
        candidate = calculation

    if reference is None or candidate is None:
        raise ValueError("Both reference and candidate tallies are required.")
    if experiment is not None and reference is not experiment:
        raise ValueError("Both reference and experiment were provided but differ.")
    if calculation is not None and candidate is not calculation:
        raise ValueError("Both candidate and calculation were provided but differ.")

    if reference.dims != candidate.dims or reference.shape != candidate.shape:
        raise ValueError("Reference and candidate tallies must share dims and shape.")

    dims = _normalize_flatten_dims(reference, flatten_dims)

    exp_da = reference._da
    calc_da = candidate._da

    exp_std = reference._da_mc_std
    if exp_std is None:
        exp_std = xr.zeros_like(exp_da)

    calc_std = candidate._da_mc_std
    if calc_std is None:
        calc_std = xr.zeros_like(calc_da)

    exp_vals = _stack_dataarray(exp_da, dims)
    calc_vals = _stack_dataarray(calc_da, dims)
    exp_unc = _stack_dataarray(exp_std, dims)
    calc_unc = _stack_dataarray(calc_std, dims)

    if exp_vals.shape != calc_vals.shape:
        raise ValueError("Reference and candidate tallies produce mismatched point counts.")

    point_ids = _format_point_ids(dims, exp_vals["point"].values)

    exp_points = [
        DataPoint(value=float(val), uncertainty=float(unc))
        for val, unc in zip(exp_vals.values, exp_unc.values)
    ]
    calc_points = [
        DataPoint(value=float(val), uncertainty=float(unc))
        for val, unc in zip(calc_vals.values, calc_unc.values)
    ]

    return exp_points, calc_points, point_ids


def compare_tallies(
    observable_name: str,
    observable_type: str,
    reference: BaseTally | None = None,
    candidate: BaseTally | None = None,
    *,
    experiment: BaseTally | None = None,
    calculation: BaseTally | None = None,
    flatten_dims: Sequence[str] | None = None,
    include_grading: bool = False,
) -> ObservableComparison:
    """Compare two tallies by converting them into point metrics.

    The flatten_dims parameter controls which tally dims are treated as
    comparison points. When omitted, all dims are flattened.
    """
    exp_points, calc_points, point_ids = datapoints_from_tally(
        reference=reference,
        candidate=candidate,
        experiment=experiment,
        calculation=calculation,
        flatten_dims=flatten_dims,
    )

    return compare_point_set(
        observable_name=observable_name,
        observable_type=observable_type,
        experiment_points=exp_points,
        calculation_points=calc_points,
        point_ids=point_ids,
        include_grading=include_grading,
    )


def compare_benchmark_results(
    benchmark_id: str,
    reference: BenchmarkResults | None = None,
    candidate: BenchmarkResults | None = None,
    *,
    experiment: BenchmarkResults | None = None,
    calculation: BenchmarkResults | None = None,
    code_name: str | None = None,
    code_version: str | None = None,
    reference_source: str | None = None,
    tally_names: Sequence[str] | None = None,
    observable_type_map: dict[str, str] | None = None,
    flatten_dims_map: dict[str, Sequence[str]] | None = None,
    include_grading: bool = False,
) -> BenchmarkComparison:
    """
    Compare two BenchmarkResults objects and aggregate to benchmark-level metrics.

    Parameters
    ----------
    benchmark_id, code_name, code_version, reference_source:
        Metadata for the benchmark comparison.
        When code_name or code_version is omitted, they are pulled from the
        candidate run metadata (fallback: reference run metadata).
    reference, candidate:
        Benchmark results with matching tally groups.
    tally_names:
        Optional list of tally names/groups to compare. Defaults to reference.tallies.
    observable_type_map:
        Optional mapping from tally name to observable type.
    flatten_dims_map:
        Optional per-tally flatten dims. Tallies missing from the map use all dims.
    include_grading:
        When False, omit status/grade fields from the comparison output.
    """
    if reference is None:
        reference = experiment
    if candidate is None:
        candidate = calculation

    if reference is None or candidate is None:
        raise ValueError("Both reference and candidate results are required.")
    if experiment is not None and reference is not experiment:
        raise ValueError("Both reference and experiment were provided but differ.")
    if calculation is not None and candidate is not calculation:
        raise ValueError("Both candidate and calculation were provided but differ.")

    if code_name is None or code_version is None:
        run_metadata = None
        try:
            run_metadata = candidate.get_run_metadata()
        except ValueError:
            try:
                run_metadata = reference.get_run_metadata()
            except ValueError:
                run_metadata = None

        if run_metadata is not None:
            if code_name is None:
                code_name = run_metadata.get("code_name")
            if code_version is None:
                code_version = run_metadata.get("code_version")

    if code_name is None:
        code_name = ""
    if code_version is None:
        code_version = ""

    if reference_source is None:
        raise ValueError("reference_source is required")

    names = list(tally_names) if tally_names is not None else list(reference.tallies)
    if not names:
        raise ValueError("No tallies found to compare.")

    observable_type_map = observable_type_map or {}
    flatten_dims_map = flatten_dims_map or {}
    observables: list[ObservableComparison] = []

    for name in names:
        exp_tally = reference.get_tally(name)
        calc_tally = candidate.get_tally(name)
        observable_type = observable_type_map.get(name, "tally")

        tally_flatten_dims = flatten_dims_map.get(name)

        obs = compare_tallies(
            observable_name=name,
            observable_type=observable_type,
            reference=exp_tally,
            candidate=calc_tally,
            flatten_dims=tally_flatten_dims,
            include_grading=include_grading,
        )
        observables.append(obs)

    return aggregate_benchmark(
        benchmark_id=benchmark_id,
        code_name=code_name,
        code_version=code_version,
        reference_source=reference_source,
        observables=observables,
        include_grading=include_grading,
    )
