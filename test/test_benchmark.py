import pytest
import yaml
from pathlib import Path
from unittest.mock import patch, mock_open

from openmc_fusion_benchmarks import Benchmark


# Minimal subclass for testing abstract class
class DummyBenchmark(Benchmark):
    def build_materials(self): return "materials"
    def build_geometry(self): return "geometry"
    def build_source(self): return "source"
    def build_settings(self): return "settings"
    def build_tallies(self): return "tallies"


@pytest.fixture
def valid_yaml():
    return yaml.dump({
        "metadata": {"title": "Dummy"},
        "materials": [],
        "geometry": {},
        "sources": [],
        "tallies": []
    })


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_benchmark_success(mock_validate, valid_yaml):
    # Full path is resolved inside the class
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    assert bench.name == "dummy"
    assert bench._benchmark_spec["metadata"]["title"] == "Dummy"
    assert bench.build_materials() == "materials"


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_benchmark_file_not_found(mock_validate):
    # Simulate missing file
    with patch.object(Path, "open", side_effect=FileNotFoundError):
        with pytest.raises(FileNotFoundError):
            DummyBenchmark("nonexistent")
