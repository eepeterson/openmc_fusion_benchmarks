import importlib.util
import pytest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch


OPENMC_AVAILABLE = importlib.util.find_spec("openmc") is not None


from openmc_fusion_benchmarks.benchmark import OpenmcBenchmark


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_postprocess_with_open_statepoint_object():
    fake_self = SimpleNamespace(_benchmark_spec={"tallies": [{"name": "t1"}, {"name": "t2"}]})
    fake_sp = object()

    with patch("openmc_fusion_benchmarks.benchmark.make_default_openmc_normalizer", return_value="norm") as mk_norm:
        with patch("openmc_fusion_benchmarks.benchmark.save_openmc_statepoint_tallies") as save:
            with patch("openmc.StatePoint", new=object):
                OpenmcBenchmark._postprocess(fake_self, statepoint=fake_sp, mesh="mesh.h5m")

    mk_norm.assert_called_once_with("mesh.h5m")
    save.assert_called_once()
    kwargs = save.call_args.kwargs
    assert kwargs["statepoint"] is fake_sp
    assert kwargs["filename"] == "benchmark_results.h5"
    assert kwargs["tally_names"] == ["t1", "t2"]
    assert kwargs["tmc_coords"] == {"realization": ["baseline"]}


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_postprocess_with_statepoint_path_uses_context_manager():
    fake_self = SimpleNamespace(_benchmark_spec={"tallies": [{"name": "t1"}]})

    class DummyStatePoint:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return SimpleNamespace(get_tally=lambda *a, **k: None)

        def __exit__(self, *_exc):
            return False

    with patch("openmc_fusion_benchmarks.benchmark.make_default_openmc_normalizer", return_value="norm"):
        with patch("openmc_fusion_benchmarks.benchmark.save_openmc_statepoint_tallies") as save:
            with patch("openmc.StatePoint", new=DummyStatePoint):
                OpenmcBenchmark._postprocess(fake_self, statepoint=Path("statepoint.10.h5"))

    save.assert_called_once()
