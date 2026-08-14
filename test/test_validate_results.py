import importlib.util
import sys
import types


def _module_available(name: str) -> bool:
    if name in sys.modules:
        return True
    try:
        return importlib.util.find_spec(name) is not None
    except ValueError:
        return False


# Allow importing package modules even when OpenMC stack is unavailable.
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


from openmc_fusion_benchmarks.validate_results import (
    _filter_bins_match,
    normalize_filter_type,
    validate_tally_consistency,
)


def test_normalize_filter_type():
    assert normalize_filter_type("CellFilter") == "cell"
    assert normalize_filter_type(" energy ") == "energy"
    assert normalize_filter_type("ParticleFilter") == "particle"


def test_filter_bins_match_energy_and_non_energy():
    assert _filter_bins_match(
        {"type": "energy", "values": [0.0, 1.0, 2.0]},
        {"bins": [0.0, 1.0, 2.0]},
    )

    assert _filter_bins_match(
        {"type": "cell", "values": [1, 2]},
        {"bins": [1, 2]},
    )

    assert not _filter_bins_match(
        {"type": "cell", "values": [1, 2]},
        {"bins": [1, 3]},
    )

    assert _filter_bins_match(
        {"type": "cell", "values": None},
        {"bins": [1, 2]},
    )
    assert _filter_bins_match(
        {"type": "cell", "values": [1, 2]},
        {"bins": None},
    )

    assert not _filter_bins_match(
        {"type": "energy", "values": ["bad"]},
        {"bins": [0.0, 1.0]},
    )


def test_validate_tally_consistency_success_and_mismatch():
    spec = {
        "particle": "neutron",
        "scores": ["flux"],
        "nuclides": ["total"],
        "filters": [{"type": "cell", "values": [1, 2]}],
    }

    observed_ok = {
        "scores": ["flux"],
        "nuclides": ["total"],
        "filters": [
            {"name": "ParticleFilter", "bins": ["neutron"]},
            {"name": "CellFilter", "bins": [1, 2]},
        ],
    }

    consistent, issues = validate_tally_consistency(spec, observed_ok)
    assert consistent
    assert issues == []

    observed_bad = {
        "scores": ["heating"],
        "nuclides": ["total"],
        "filters": [
            {"name": "ParticleFilter", "bins": ["photon"]},
            {"name": "SurfaceFilter", "bins": [1, 2]},
        ],
    }

    consistent, issues = validate_tally_consistency(spec, observed_bad)
    assert not consistent
    assert any("scores mismatch" in issue for issue in issues)
    assert any("particle mismatch" in issue for issue in issues)
    assert any("filter type/order mismatch" in issue for issue in issues)


def test_validate_tally_consistency_missing_particle_filter():
    spec = {
        "particle": "neutron",
        "scores": ["flux"],
        "nuclides": ["total"],
        "filters": [{"type": "cell", "values": [1]}],
    }

    observed = {
        "scores": ["flux"],
        "nuclides": ["total"],
        "filters": [{"name": "CellFilter", "bins": [1]}],
    }

    consistent, issues = validate_tally_consistency(spec, observed)
    assert not consistent
    assert any("missing ParticleFilter" in issue for issue in issues)
