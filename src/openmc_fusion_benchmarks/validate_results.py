from __future__ import annotations

import numpy as np


def normalize_filter_type(type_name: str) -> str:
    """Normalize filter type naming across spec and backend representations."""
    name = str(type_name).strip().lower()
    if name.endswith("filter"):
        name = name[:-6]
    return name


def _filter_bins_match(spec_filter: dict, observed_axis: dict) -> bool:
    """Compare filter bin definitions from spec and observed metadata."""
    expected = spec_filter.get("values")
    observed = observed_axis.get("bins")

    if expected is None or observed is None:
        return True

    ftype = normalize_filter_type(spec_filter.get("type", ""))
    if ftype == "energy":
        try:
            return np.allclose(np.asarray(expected, dtype=float), np.asarray(observed, dtype=float))
        except Exception:
            return False

    try:
        return list(expected) == list(observed)
    except Exception:
        return False


def validate_tally_consistency(spec_tally: dict, observed_tally: dict) -> tuple[bool, list[str]]:
    """Validate observed tally metadata against repository specification."""
    issues: list[str] = []

    observed_scores = [str(s) for s in observed_tally.get("scores", [])]
    observed_nuclides = [str(n) for n in observed_tally.get("nuclides", [])]
    observed_filters = list(observed_tally.get("filters", []))

    spec_scores = [str(s) for s in spec_tally.get("scores", [])]
    if spec_scores and spec_scores != observed_scores:
        issues.append(f"scores mismatch: expected {spec_scores}, observed {observed_scores}")

    spec_nuclides = [str(n) for n in spec_tally.get("nuclides", [])]
    if spec_nuclides and spec_nuclides != observed_nuclides:
        issues.append(f"nuclides mismatch: expected {spec_nuclides}, observed {observed_nuclides}")

    expected_particle = spec_tally.get("particle")
    if expected_particle is not None:
        particle_filters = [a for a in observed_filters if normalize_filter_type(a.get("name", "")) == "particle"]
        if not particle_filters:
            issues.append("missing ParticleFilter in observed tally")
        else:
            bins = particle_filters[0].get("bins", [])
            observed_particle = str(bins[0]) if bins else None
            if str(expected_particle) != str(observed_particle):
                issues.append(
                    f"particle mismatch: expected {expected_particle}, observed {observed_particle}"
                )

    spec_filters = spec_tally.get("filters", [])
    observed_non_particle = [
        a for a in observed_filters if normalize_filter_type(a.get("name", "")) != "particle"
    ]

    expected_types = [normalize_filter_type(f.get("type", "")) for f in spec_filters]
    observed_types = [normalize_filter_type(a.get("name", "")) for a in observed_non_particle]
    if expected_types != observed_types:
        issues.append(f"filter type/order mismatch: expected {expected_types}, observed {observed_types}")


    return len(issues) == 0, issues