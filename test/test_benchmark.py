import importlib.util
import pytest
import yaml
import numpy as np
from pathlib import Path
from unittest.mock import patch, mock_open, Mock, MagicMock

from openmc_fusion_benchmarks import Benchmark
from openmc_fusion_benchmarks.benchmark import OpenmcBenchmark

OPENMC_AVAILABLE = importlib.util.find_spec("openmc") is not None


# Minimal subclass for testing abstract class
class DummyBenchmark(Benchmark):
    def _build_materials(self):
        return "materials"

    def _build_geometry(self):
        return "geometry"

    def _build_source(self):
        return "source"

    def _build_settings(self):
        return "settings"

    def _build_tallies(self):
        return "tallies"

    def _build_model(self):
        return "model"

    def _postprocess(self):
        return "postprocess"

    def _uncertainty_quantification(self):
        return "uncertainty_quantification"

    def run(self):
        return "run"


@pytest.fixture
def valid_yaml():
    return yaml.dump({
        "metadata": {"title": "Dummy"},
        "materials": [],
        "geometry": {},
        "sources": [],
        "tallies": [],
    })


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_benchmark_success(mock_validate, valid_yaml):
    # Full path is resolved inside the class
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    assert bench.name == "dummy"
    assert bench._benchmark_spec["metadata"]["title"] == "Dummy"
    assert bench._build_materials() == "materials"


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_benchmark_file_not_found(mock_validate):
    # Simulate missing file
    with patch.object(Path, "open", side_effect=FileNotFoundError):
        with pytest.raises(FileNotFoundError):
            DummyBenchmark("nonexistent")


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_benchmark_metadata_parsing(mock_validate, valid_yaml):
    """Test that metadata is parsed and formatted correctly."""
    metadata_yaml = yaml.dump({
        "metadata": {
            "title": "Test Benchmark",
            "type": "experimental",
            "category": "fusion",
            "version": "1.0",
            "description": "Test description",
            "date": "2025-01-01",
            "location": {
                "facility": "Test Facility",
                "city": "Test City",
                "country": "Test Country"
            },
            "references": [
                {"title": "Test Paper", "doi": "10.1234/test"}
            ],
            "authors": [
                {"name": "Test Author", "affiliation": "Test Org", "email": "test@test.com"}
            ]
        },
        "materials": [],
        "geometry": {},
        "sources": [],
        "tallies": [],
    })
    
    with patch.object(Path, "open", mock_open(read_data=metadata_yaml)):
        bench = DummyBenchmark("test")
    
    assert "Test Benchmark" in bench.metadata
    assert "experimental" in bench.metadata
    assert "fusion" in bench.metadata
    assert "Test Facility" in bench.metadata
    assert "Test Author" in bench.metadata


def test_read_metadata_with_full_spec():
    """Test _read_metadata with complete metadata including authors."""
    metadata_yaml = yaml.dump({
        "metadata": {
            "title": "Full Test",
            "type": "experimental",
            "category": "fusion",
            "version": "2.0",
            "description": "Complete test",
            "date": "2025-01-01",
            "location": {
                "facility": "Test Lab",
                "city": "Boston",
                "country": "USA"
            },
            "references": [
                {"title": "Paper 1", "doi": "10.1234/test", "url": "http://example.com"}
            ],
            "authors": [
                {"name": "Alice", "affiliation": "MIT", "email": "alice@mit.edu"},
                {"name": "Bob", "affiliation": "Harvard", "email": "bob@harvard.edu"}
            ]
        },
        "materials": [],
        "geometry": {},
        "sources": [],
        "tallies": [],
    })
    
    with patch("openmc_fusion_benchmarks.benchmark.validate_benchmark"):
        with patch.object(Path, "open", mock_open(read_data=metadata_yaml)):
            bench = DummyBenchmark("test")
    
    assert "Full Test" in bench.metadata
    assert "experimental" in bench.metadata
    assert "Alice" in bench.metadata
    assert "Bob" in bench.metadata
    assert "MIT" in bench.metadata


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_generate_report_uses_default_config(mock_validate, valid_yaml, tmp_path, monkeypatch):
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    results_path = tmp_path / "benchmark_results.h5"
    import h5py
    with h5py.File(results_path, "w"):
        pass
    monkeypatch.chdir(tmp_path)

    captured = {}

    def _fake_from_file(path):
        captured["from_file"] = path
        return Mock(filepath=Path(path))

    def _fake_from_database(name, filename="reference_results.h5"):
        captured["from_database"] = (name, filename)
        return Mock(filepath=Path(filename))

    def _fake_build_report(sources, config):
        captured["config"] = config
        captured["sources"] = sources
        return Mock()

    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.BenchmarkResults.from_file",
        _fake_from_file,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.BenchmarkResults.from_database",
        _fake_from_database,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.build_report",
        _fake_build_report,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.render_yaml",
        lambda report, path: path,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.render_pdf",
        lambda report, path, plots_dir: path,
    )

    bench._generate_report(report_config=None)

    assert captured["from_database"] == ("dummy", "experiment.h5")
    assert captured["config"].include_yaml is True
    assert captured["config"].include_pdf is True
    assert captured["config"].output_dir == Path("report")


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_generate_report_missing_results_warns(mock_validate, valid_yaml, tmp_path):
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    with pytest.warns(UserWarning, match="benchmark_results.h5 not found"):
        bench._generate_report(report_config=None)


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_write_spec_snapshot_writes_group(mock_validate, valid_yaml, tmp_path):
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    results_path = tmp_path / "benchmark_results.h5"
    import h5py
    with h5py.File(results_path, "w"):
        pass

    bench._write_spec_snapshot(str(results_path))

    import h5py

    with h5py.File(results_path, "r") as handle:
        assert "specifications" in handle
        assert "yaml" in handle["specifications"]


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_write_run_metadata_writes_attrs(mock_validate, valid_yaml, tmp_path):
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    results_path = tmp_path / "benchmark_results.h5"
    import h5py
    with h5py.File(results_path, "w"):
        pass

    bench._write_run_metadata(
        code_name="openmc",
        code_version="0.15.2",
        nuclear_data_name="endfb",
        nuclear_data_version="8.0",
        geometry="cad",
        filename=str(results_path),
    )

    import h5py

    with h5py.File(results_path, "r") as handle:
        attrs = handle["run_metadata"].attrs
        assert attrs["code_name"] == "openmc"
        assert attrs["code_version"] == "0.15.2"
        assert attrs["nuclear_data_name"] == "endfb"
        assert attrs["nuclear_data_version"] == "8.0"
        assert attrs["geometry"] == "cad"


@patch("openmc_fusion_benchmarks.benchmark.validate_benchmark")
def test_generate_report_without_reference(mock_validate, valid_yaml, tmp_path, monkeypatch):
    with patch.object(Path, "open", mock_open(read_data=valid_yaml)):
        bench = DummyBenchmark("dummy")

    results_path = tmp_path / "benchmark_results.h5"
    results_path.write_text("data")
    monkeypatch.chdir(tmp_path)

    def _fake_from_file(path):
        return Mock(filepath=Path(path))

    def _fake_from_database(_name, filename="reference_results.h5"):
        raise FileNotFoundError("missing")

    captured = {}

    def _fake_build_report(sources, config):
        captured["sources"] = sources
        return Mock()

    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.BenchmarkResults.from_file",
        _fake_from_file,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.BenchmarkResults.from_database",
        _fake_from_database,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.build_report",
        _fake_build_report,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.render_yaml",
        lambda report, path: path,
    )
    monkeypatch.setattr(
        "openmc_fusion_benchmarks.benchmark.render_pdf",
        lambda report, path, plots_dir: path,
    )

    bench._generate_report(report_config=None)

    assert captured["sources"][0].kind == "calculation"
    assert len(captured["sources"]) == 1


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_openmc_benchmark_build_materials():
    """Test OpenmcBenchmark._build_materials with atomic fractions."""
    import openmc
    
    spec = {
        "metadata": {"title": "Test"},
        "materials": [
            {
                "id": 1,
                "name": "aluminum",
                "density": {"value": 2.7, "units": "g/cm3"},
                "composition": {
                    "composition_type": "nuclide",
                    "fraction_type": "atomic",
                    "data": {"Al27": 1.0}
                }
            }
        ],
        "geometry": {"cad_file": "test.step", "meshing": {"volumes": [], "global_mesh_size_min": 0.1, "global_mesh_size_max": 10}},
        "settings": {"run_mode": "fixed_source", "batches": 10, "particles_per_batch": 100, "photon_transport": False},
        "sources": [],
        "tallies": []
    }
    
    with patch("openmc_fusion_benchmarks.benchmark.validate_benchmark"):
        with patch.object(Path, "open", mock_open(read_data=yaml.dump(spec))):
            with patch.object(OpenmcBenchmark, '_build_geometry', return_value=openmc.Geometry()):
                with patch.object(OpenmcBenchmark, '_build_settings', return_value=openmc.Settings()):
                    with patch.object(OpenmcBenchmark, '_build_tallies', return_value=openmc.Tallies()):
                        bench = OpenmcBenchmark("test")
                        materials = bench._build_materials()
    
    assert isinstance(materials, openmc.Materials)
    assert len(materials) == 1
    assert materials[0].name == "aluminum"


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_openmc_benchmark_build_materials_with_elements():
    """Test OpenmcBenchmark._build_materials with element composition."""
    import openmc
    
    spec = {
        "metadata": {"title": "Test"},
        "materials": [
            {
                "id": 2,
                "name": "steel",
                "density": {"value": 7.8, "units": "g/cm3"},
                "composition": {
                    "composition_type": "element",
                    "fraction_type": "weight",
                    "data": {"Fe": 0.98, "C": 0.02}
                }
            }
        ],
        "geometry": {"cad_file": "test.step", "meshing": {"volumes": [], "global_mesh_size_min": 0.1, "global_mesh_size_max": 10}},
        "settings": {"run_mode": "fixed_source", "batches": 10, "particles_per_batch": 100, "photon_transport": False},
        "sources": [],
        "tallies": []
    }
    
    with patch("openmc_fusion_benchmarks.benchmark.validate_benchmark"):
        with patch.object(Path, "open", mock_open(read_data=yaml.dump(spec))):
            with patch.object(OpenmcBenchmark, '_build_geometry', return_value=openmc.Geometry()):
                with patch.object(OpenmcBenchmark, '_build_settings', return_value=openmc.Settings()):
                    with patch.object(OpenmcBenchmark, '_build_tallies', return_value=openmc.Tallies()):
                        bench = OpenmcBenchmark("test")
                        materials = bench._build_materials()
    
    assert len(materials) == 1
    assert materials[0].name == "steel"


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_openmc_benchmark_invalid_fraction_type():
    """Test that invalid fraction type raises ValueError."""
    import openmc
    
    spec = {
        "metadata": {"title": "Test"},
        "materials": [
            {
                "id": 1,
                "name": "bad_material",
                "density": {"value": 1.0, "units": "g/cm3"},
                "composition": {
                    "composition_type": "nuclide",
                    "fraction_type": "invalid_type",
                    "data": {"H1": 1.0}
                }
            }
        ],
        "geometry": {"cad_file": "test.step", "meshing": {"volumes": [], "global_mesh_size_min": 0.1, "global_mesh_size_max": 10}},
        "settings": {"run_mode": "fixed_source", "batches": 10, "particles_per_batch": 100, "photon_transport": False},
        "sources": [],
        "tallies": []
    }
    
    with patch("openmc_fusion_benchmarks.benchmark.validate_benchmark"):
        with patch.object(Path, "open", mock_open(read_data=yaml.dump(spec))):
            with patch.object(OpenmcBenchmark, '_build_geometry', return_value=openmc.Geometry()):
                with patch.object(OpenmcBenchmark, '_build_settings', return_value=openmc.Settings()):
                    with patch.object(OpenmcBenchmark, '_build_tallies', return_value=openmc.Tallies()):
                        with pytest.raises(ValueError, match="Invalid fraction type"):
                            OpenmcBenchmark("test")

