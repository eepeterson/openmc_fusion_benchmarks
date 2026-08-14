from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from enum import Enum
import numpy as np


class PointStatus(str, Enum):
    OK = "ok"
    WARNING = "warning"
    OUTLIER = "outlier"
    ERROR = "error"

class BenchmarkStatus(str, Enum):
    ACCEPTABLE = "acceptable"
    BORDERLINE = "borderline"
    PROBLEMATIC = "problematic"

@dataclass
class DataPoint:
    """A single experimental or calculated measurement."""
    value: float
    uncertainty: float
    
    @property
    def relative_uncertainty(self) -> float:
        """Relative uncertainty as a fraction."""
        return abs(self.uncertainty / self.value) if self.value != 0 else np.inf
    

@dataclass
class PointMetrics:
    """Metrics for a single comparison point."""
    c_over_e: float
    relative_deviation: float
    absolute_deviation: float
    combined_uncertainty: float
    normalized_residual: float
    chi2_contribution: float
    status: Optional[PointStatus] = None


@dataclass
class PointComparison:
    """One comparison point (detector, foil, energy bin, etc.)."""
    id: str
    observable_type: str  # "reaction_rate", "spectrum", "leakage", etc.
    experiment: DataPoint
    calculation: DataPoint
    metrics: Optional[PointMetrics] = None
    covariance_index: Optional[int] = None  # Index in covariance matrix if available

    def __repr__(self) -> str:
        metrics = self.metrics
        if metrics is None:
            status = "uncomputed"
        elif metrics.status is None:
            status = "ungraded"
        else:
            status = metrics.status.value
        return (
            f"<PointComparison id={self.id!r}, type={self.observable_type!r}, "
            f"calc={self.calculation.value:.6g}, exp={self.experiment.value:.6g}, "
            f"status={status}>"
        )


@dataclass
class ObservableComparison:
    """Comparison of all points for one observable (one tally, one spectrum, one foil traverse)."""
    name: str
    observable_type: str
    points: list[PointComparison]
    
    # Observable-level aggregates
    mean_bias: float = field(default=0.0)
    mean_abs_relative_deviation: float = field(default=0.0)
    rms_relative_deviation: float = field(default=0.0)
    mean_abs_normalized_residual: float = field(default=0.0)
    reduced_chi2: float = field(default=0.0)
    fraction_within_1sigma: Optional[float] = None
    fraction_within_2sigma: Optional[float] = None
    fraction_within_3sigma: Optional[float] = None
    outlier_fraction: Optional[float] = None
    pass_count: Optional[int] = None
    warning_count: Optional[int] = None
    outlier_count: Optional[int] = None

    def __repr__(self) -> str:
        return (
            f"<ObservableComparison name={self.name!r}, type={self.observable_type!r}, "
            f"points={len(self.points)}, mean_bias={self.mean_bias:.3g}, "
            f"rms_rel_dev={self.rms_relative_deviation:.3g}, "
            f"reduced_chi2={self.reduced_chi2:.3g}>"
        )


@dataclass
class BenchmarkComparison:
    """Full validation comparison for one benchmark run against reference data."""
    benchmark_id: str
    code_name: str
    code_version: str
    reference_source: str

    observables: list[ObservableComparison] = field(default_factory=list)

    weighted_mean_bias: float = field(default=0.0)
    weighted_rms_relative_deviation: float = field(default=0.0)
    global_reduced_chi2: float = field(default=0.0)
    total_point_count: int = field(default=0)
    outlier_fraction: Optional[float] = None

    benchmark_status: Optional[BenchmarkStatus] = None
    dashboard_score: Optional[float] = None

    quality_flags: Dict[str, bool] = field(default_factory=dict)
    assumptions: Dict[str, Any] = field(default_factory=dict)