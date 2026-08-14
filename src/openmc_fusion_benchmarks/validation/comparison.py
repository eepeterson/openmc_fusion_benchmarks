from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Sequence

import numpy as np

from .model import (
    BenchmarkComparison,
    BenchmarkStatus,
    DataPoint,
    ObservableComparison,
    PointComparison,
)
from .metrics import compute_point_metrics

_DEFAULT_GRADING_THRESHOLDS = {
    "rms_relative_deviation": {"good": 0.05, "bad": 0.15},
    "reduced_chi2": {"good": 1.5, "bad": 3.0},
    "outlier_fraction": {"good": 0.01, "bad": 0.05},
}

_DEFAULT_SCORE_WEIGHTS = {
    "rms_relative_deviation": 0.5,
    "reduced_chi2": 0.3,
    "outlier_fraction": 0.2,
}


def _as_datapoint(item) -> DataPoint:
    """
    Normalize an input item to DataPoint.

    Expected input can be:
    - DataPoint
    - dict with keys: value, uncertainty
    - xarray / object wrapper if you adapt it later
    """
    if isinstance(item, DataPoint):
        return item
    if isinstance(item, dict):
        return DataPoint(
            value=float(item["value"]),
            uncertainty=float(item["uncertainty"]),
        )
    raise TypeError(f"Unsupported point type: {type(item)!r}")


def compare_point_set(
    observable_name: str,
    observable_type: str,
    experiment_points: Sequence,
    calculation_points: Sequence,
    point_ids: Sequence[str] | None = None,
    include_grading: bool = False,
) -> ObservableComparison:
    """
    Compare one observable point-by-point.

    Parameters
    ----------
    observable_name:
        Name of the tally / observable.
    observable_type:
        Category such as reaction_rate, spectrum, leakage.
    experiment_points:
        Reference values.
    calculation_points:
        Calculated values from the transport code.
    point_ids:
        Optional point labels. If omitted, indices are used.
    """
    if len(experiment_points) != len(calculation_points):
        raise ValueError("Experiment and calculation point counts must match.")

    n = len(experiment_points)
    ids = point_ids if point_ids is not None else [str(i) for i in range(n)]

    if len(ids) != n:
        raise ValueError("point_ids length must match number of points.")

    points: list[PointComparison] = []

    for pid, exp_raw, calc_raw in zip(ids, experiment_points, calculation_points):
        exp = _as_datapoint(exp_raw)
        calc = _as_datapoint(calc_raw)

        point = PointComparison(
            id=str(pid),
            observable_type=observable_type,
            experiment=exp,
            calculation=calc,
        )
        point.metrics = compute_point_metrics(point, include_grading=include_grading)
        points.append(point)

    return _aggregate_observable(
        observable_name,
        observable_type,
        points,
        include_grading=include_grading,
    )


def _aggregate_observable(
    name: str,
    observable_type: str,
    points: Sequence[PointComparison],
    include_grading: bool = False,
) -> ObservableComparison:
    """Compute observable-level summary metrics."""
    if not points:
        return ObservableComparison(name=name, observable_type=observable_type, points=[])

    metrics = [p.metrics for p in points if p.metrics is not None]
    if len(metrics) != len(points):
        raise ValueError("All points must have computed metrics.")

    rel_dev = np.array([m.relative_deviation for m in metrics], dtype=float)
    norm_res = np.array([m.normalized_residual for m in metrics], dtype=float)
    chi2 = np.array([m.chi2_contribution for m in metrics], dtype=float)

    total = len(points)
    within_1 = None
    within_2 = None
    beyond_3 = None
    within_3 = None
    outlier_fraction = None
    pass_count = None
    warning_count = None
    outlier_count = None

    if include_grading:
        within_1 = int(np.sum(np.abs(norm_res) <= 1.0))
        within_2 = int(np.sum(np.abs(norm_res) <= 2.0))
        beyond_3 = int(np.sum(np.abs(norm_res) > 3.0))
        within_3 = total - beyond_3
        outlier_fraction = float(beyond_3 / total)
        pass_count = within_1
        warning_count = within_2 - within_1
        outlier_count = beyond_3

    return ObservableComparison(
        name=name,
        observable_type=observable_type,
        points=list(points),
        mean_bias=float(np.mean(rel_dev)),
        mean_abs_relative_deviation=float(np.mean(np.abs(rel_dev))),
        rms_relative_deviation=float(np.sqrt(np.mean(rel_dev**2))),
        mean_abs_normalized_residual=float(np.mean(np.abs(norm_res))),
        reduced_chi2=float(np.sum(chi2) / max(total - 1, 1)),
        fraction_within_1sigma=None if within_1 is None else float(within_1 / total),
        fraction_within_2sigma=None if within_2 is None else float(within_2 / total),
        fraction_within_3sigma=None if within_3 is None else float(within_3 / total),
        outlier_fraction=outlier_fraction,
        pass_count=pass_count,
        warning_count=warning_count,
        outlier_count=outlier_count,
    )


def _score_component(value: float, good: float, bad: float) -> float:
    if bad <= good:
        raise ValueError("Score thresholds must satisfy bad > good.")
    if value <= good:
        return 1.0
    if value >= bad:
        return 0.0
    return 1.0 - (value - good) / (bad - good)


def _grade_benchmark(
    rms_relative_deviation: float,
    reduced_chi2: float,
    outlier_fraction: float,
) -> tuple[BenchmarkStatus, float]:
    thresholds = _DEFAULT_GRADING_THRESHOLDS
    weights = _DEFAULT_SCORE_WEIGHTS

    metrics = {
        "rms_relative_deviation": rms_relative_deviation,
        "reduced_chi2": reduced_chi2,
        "outlier_fraction": outlier_fraction,
    }

    status = BenchmarkStatus.ACCEPTABLE
    for name, value in metrics.items():
        good = thresholds[name]["good"]
        bad = thresholds[name]["bad"]
        if value >= bad:
            status = BenchmarkStatus.PROBLEMATIC
            break
        if value > good:
            status = BenchmarkStatus.BORDERLINE

    total_weight = sum(weights.values())
    score = 0.0
    for name, value in metrics.items():
        good = thresholds[name]["good"]
        bad = thresholds[name]["bad"]
        component = _score_component(value, good, bad)
        score += weights[name] * component
    score = 100.0 * score / total_weight if total_weight > 0 else 0.0

    return status, score


def aggregate_benchmark(
    benchmark_id: str,
    code_name: str,
    code_version: str,
    reference_source: str,
    observables: Sequence[ObservableComparison],
    include_grading: bool = False,
) -> BenchmarkComparison:
    """Aggregate observable-level comparisons into a benchmark-level result."""
    if not observables:
        return BenchmarkComparison(
            benchmark_id=benchmark_id,
            code_name=code_name,
            code_version=code_version,
            reference_source=reference_source,
            observables=[],
            benchmark_status=None if not include_grading else BenchmarkStatus.ACCEPTABLE,
        )

    total_points = sum(len(obs.points) for obs in observables)
    if total_points == 0:
        return BenchmarkComparison(
            benchmark_id=benchmark_id,
            code_name=code_name,
            code_version=code_version,
            reference_source=reference_source,
            observables=list(observables),
            benchmark_status=None if not include_grading else BenchmarkStatus.ACCEPTABLE,
        )

    weights = np.array([len(obs.points) for obs in observables], dtype=float)
    mean_biases = np.array([obs.mean_bias for obs in observables], dtype=float)
    rms_devs = np.array([obs.rms_relative_deviation for obs in observables], dtype=float)
    chi2_vals = np.array([obs.reduced_chi2 for obs in observables], dtype=float)

    weighted_mean_bias = float(np.average(mean_biases, weights=weights))
    weighted_rms_relative_deviation = float(np.average(rms_devs, weights=weights))
    global_reduced_chi2 = float(np.average(chi2_vals, weights=weights))

    outlier_fraction = None
    benchmark_status = None
    dashboard_score = None

    if include_grading:
        total_outliers = int(
            sum(obs.outlier_count for obs in observables if obs.outlier_count is not None)
        )
        outlier_fraction = float(total_outliers / total_points)
        benchmark_status, dashboard_score = _grade_benchmark(
            rms_relative_deviation=weighted_rms_relative_deviation,
            reduced_chi2=global_reduced_chi2,
            outlier_fraction=outlier_fraction,
        )

    return BenchmarkComparison(
        benchmark_id=benchmark_id,
        code_name=code_name,
        code_version=code_version,
        reference_source=reference_source,
        observables=list(observables),
        weighted_mean_bias=weighted_mean_bias,
        weighted_rms_relative_deviation=weighted_rms_relative_deviation,
        global_reduced_chi2=global_reduced_chi2,
        total_point_count=total_points,
        outlier_fraction=outlier_fraction,
        benchmark_status=benchmark_status,
        dashboard_score=dashboard_score,
    )