import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from openmc_fusion_benchmarks.uq.tmc_manager import TMCManager, TMCStatePoint

try:
    OPENMC_AVAILABLE = importlib.util.find_spec("openmc") is not None
except ValueError:
    OPENMC_AVAILABLE = False


class FakeTally:
    def __init__(self, tid, name, filters, nuclides, scores, mean, std_dev):
        self.id = tid
        self.name = name
        self.filters = filters
        self.nuclides = nuclides
        self.scores = scores
        self.mean = mean
        self.std_dev = std_dev


class FakeStatePoint:
    def __init__(self, *args, **kwargs):
        self.tallies = kwargs.pop("tallies")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
@pytest.mark.parametrize(
    "mode, record",
    [
        ("sequential", {"perturbation": 0, "realization": 0, "statepoint": "statepoint.1.h5"}),
        ("matrix", {"mode": "matrix", "indices": [0], "statepoint": "statepoint.1.h5"}),
        ("diagonal", {"mode": "diagonal", "indices": [0], "statepoint": "statepoint.1.h5"}),
    ],
)
def test_tmc_statepoint_mode_is_saved(tmp_path, monkeypatch, mode, record):
    import openmc

    filters = [openmc.CellFilter([1])]
    mean = np.array([[[1.0]]])
    std_dev = np.array([[[0.1]]])
    tally = FakeTally(
        tid=1,
        name="tally_1",
        filters=filters,
        nuclides=["U235"],
        scores=["flux"],
        mean=mean,
        std_dev=std_dev,
    )

    def _fake_statepoint_ctor(*args, **kwargs):
        return FakeStatePoint(tallies={1: tally})

    monkeypatch.setattr(openmc, "StatePoint", _fake_statepoint_ctor)

    sp_path = tmp_path / "statepoint.1.h5"
    sp_path.write_text("stub")

    manifest = tmp_path / "tmc_manifest.jsonl"
    manifest.write_text(json.dumps(record) + "\n")

    manager = TMCManager(base_model=None, perturbations=[], realizations=1)
    manager._process_tmc(manifest_path=manifest)

    tmc_statepoint = tmp_path / "tmc_statepoint.1.h5"
    assert tmc_statepoint.exists()

    sp = TMCStatePoint(tmc_statepoint)
    assert sp.tmc_mode == mode
