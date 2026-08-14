import pytest
import yaml
from pathlib import Path
from unittest.mock import patch, mock_open, MagicMock
from openmc_fusion_benchmarks.validate_spec import validate_benchmark


@pytest.fixture
def mock_paths(tmp_path):
    """Creates mock schema and benchmark YAML files."""
    base = tmp_path / "src" / "openmc_fusion_benchmarks" / "benchmarks"
    base.mkdir(parents=True)
    schema_path = base / "benchmark_schema.yaml"
    benchmark_dir = base / "dummy_benchmark"
    benchmark_dir.mkdir()
    benchmark_path = benchmark_dir / "specifications.yaml"
    return schema_path, benchmark_path, benchmark_dir.name


def test_file_not_found(mock_paths):
    schema_path, benchmark_path, benchmark_name = mock_paths
    # Don't create the benchmark YAML file
    with pytest.raises(FileNotFoundError):
        validate_benchmark("dummy_benchmark")


@patch("openmc_fusion_benchmarks.validate_spec.yaml.safe_load")
@patch("openmc_fusion_benchmarks.validate_spec.open")
@patch("openmc_fusion_benchmarks.validate_spec.Path.is_file")
def test_valid_yaml_schema_validation_passes(mock_is_file, mock_open_file, mock_safe_load, capsys):
    # Arrange
    mock_is_file.return_value = True

    schema = {
        "$id": "https://openmc-fusion/schemas/benchmark_schema",
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "metadata": {"type": "object"},
            "materials": {"type": "array"}
        },
        "required": ["metadata", "materials"]
    }

    benchmark = {
        "metadata": {},
        "materials": []
    }

    # Simulate loading schema, then benchmark data
    mock_safe_load.side_effect = [schema, benchmark]
    mock_open_file.return_value.__enter__.return_value = MagicMock()

    # Act
    validate_benchmark("dummy_benchmark")
    out = capsys.readouterr().out

    # Assert
    assert "✅ dummy_benchmark is valid!" in out


@patch("openmc_fusion_benchmarks.validate_spec.yaml.safe_load")
@patch("openmc_fusion_benchmarks.validate_spec.open")
@patch("openmc_fusion_benchmarks.validate_spec.Path.is_file")
def test_invalid_yaml_schema_validation_fails(mock_is_file, mock_open_file, mock_safe_load, capsys):
    # Arrange
    mock_is_file.return_value = True

    schema = {
        "$id": "https://openmc-fusion/schemas/benchmark_schema",
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "metadata": {"type": "object"},
            "materials": {"type": "array"}
        },
        "required": ["metadata", "materials"]
    }

    benchmark = {
        "materials": []  # missing metadata!
    }

    mock_safe_load.side_effect = [schema, benchmark]
    mock_open_file.return_value.__enter__.return_value = MagicMock()

    # Act
    validate_benchmark("dummy_benchmark")
    out = capsys.readouterr().out

    # Assert
    assert "❌" in out
    assert "metadata" in out
