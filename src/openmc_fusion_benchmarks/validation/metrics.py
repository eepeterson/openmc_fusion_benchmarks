import numpy as np
from .model import PointComparison, PointMetrics, PointStatus


def compute_point_metrics(point: PointComparison, include_grading: bool = False) -> PointMetrics:
    """Compute all metrics for a single comparison point."""
    exp = point.experiment
    calc = point.calculation
    
    c = calc.value
    e = exp.value
    u_c = calc.uncertainty
    u_e = exp.uncertainty
    
    c_over_e = c / e if e != 0 else np.inf
    rel_dev = (c - e) / e if e != 0 else np.inf
    abs_dev = abs(c - e)
    
    combined_u = np.sqrt(u_e**2 + u_c**2)
    norm_residual = (c - e) / combined_u if combined_u > 0 else np.inf
    
    chi2_contrib = (c - e)**2 / combined_u**2 if combined_u > 0 else 0
    
    status = None
    if include_grading:
        # Classify this point
        if abs(norm_residual) > 3:
            status = PointStatus.OUTLIER
        elif abs(norm_residual) > 2:
            status = PointStatus.WARNING
        else:
            status = PointStatus.OK
    
    return PointMetrics(
        c_over_e=c_over_e,
        relative_deviation=rel_dev,
        absolute_deviation=abs_dev,
        combined_uncertainty=combined_u,
        normalized_residual=norm_residual,
        chi2_contribution=chi2_contrib,
        status=status
    )