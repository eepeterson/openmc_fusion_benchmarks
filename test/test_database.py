import pytest
from pathlib import Path
from openmc_fusion_benchmarks.database import (
    _resolve_database_path,
    list_database_benchmarks,
    list_database_files
)


def test_list_database_benchmarks_not_empty():
    """Test that the database has at least one benchmark."""
    benchmarks = list_database_benchmarks()
    assert isinstance(benchmarks, list)
    assert len(benchmarks) > 0  # Should find oktavian_al


def test_list_database_benchmarks_contains_strings():
    """Test that all benchmark names are strings."""
    benchmarks = list_database_benchmarks()
    assert all(isinstance(name, str) for name in benchmarks)


def test_list_database_files_existing_benchmark():
    """Test listing files for an existing benchmark."""
    benchmarks = list_database_benchmarks()
    if benchmarks:
        # Test with the first available benchmark
        files = list_database_files(benchmarks[0])
        assert isinstance(files, list)
        # All files should end with .h5
        assert all(f.endswith('.h5') for f in files)


def test_list_database_files_nonexistent_benchmark():
    """Test listing files for a nonexistent benchmark returns empty list."""
    files = list_database_files("nonexistent_benchmark_xyz")
    assert isinstance(files, list)
    assert len(files) == 0


def test_resolve_database_path_nonexistent_file():
    """Test that resolving a nonexistent file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        _resolve_database_path("nonexistent_benchmark", "nonexistent_file.h5")


def test_resolve_database_path_existing_file():
    """Test resolving an existing database file path."""
    benchmarks = list_database_benchmarks()
    if benchmarks:
        files = list_database_files(benchmarks[0])
        if files:
            # Test with first available file
            path = _resolve_database_path(benchmarks[0], files[0])
            assert isinstance(path, Path)
            assert path.exists()
            assert path.suffix == '.h5'


def test_resolve_database_path_provides_helpful_error():
    """Test that FileNotFoundError contains helpful benchmark list."""
    try:
        _resolve_database_path("fake_benchmark", "fake_file.h5")
        assert False, "Should have raised FileNotFoundError"
    except FileNotFoundError as e:
        # Error message should contain available benchmarks
        assert "Available benchmarks:" in str(e)


def test_list_database_benchmarks_fallback(monkeypatch):
    def _raise_files(_pkg):
        raise ModuleNotFoundError("no resources")

    monkeypatch.setattr("openmc_fusion_benchmarks.database.files", _raise_files)
    benchmarks = list_database_benchmarks()
    assert "oktavian_al" in benchmarks


def test_resolve_database_path_fallback(monkeypatch):
    def _raise_files(_pkg):
        raise ModuleNotFoundError("no resources")

    monkeypatch.setattr("openmc_fusion_benchmarks.database.files", _raise_files)
    path = _resolve_database_path("oktavian_al", "experiment.h5")
    assert path.exists()


def test_list_database_files_fallback(monkeypatch):
    def _raise_files(_pkg):
        raise ModuleNotFoundError("no resources")

    monkeypatch.setattr("openmc_fusion_benchmarks.database.files", _raise_files)
    files = list_database_files("oktavian_al")
    assert "experiment.h5" in files


def test_list_database_benchmarks_as_file_failure(monkeypatch):
    def _raise_files(_pkg):
        class Dummy:
            pass
        return Dummy()

    def _raise_as_file(_path):
        raise IsADirectoryError("dir")

    monkeypatch.setattr("openmc_fusion_benchmarks.database.files", _raise_files)
    monkeypatch.setattr("openmc_fusion_benchmarks.database.as_file", _raise_as_file)

    benchmarks = list_database_benchmarks()
    assert "oktavian_al" in benchmarks


def test_list_database_files_as_file_failure(monkeypatch):
    def _raise_files(_pkg):
        class Dummy:
            pass
        return Dummy()

    def _raise_as_file(_path):
        raise IsADirectoryError("dir")

    monkeypatch.setattr("openmc_fusion_benchmarks.database.files", _raise_files)
    monkeypatch.setattr("openmc_fusion_benchmarks.database.as_file", _raise_as_file)

    files = list_database_files("oktavian_al")
    assert "experiment.h5" in files


def test_resolve_database_path_missing_dev_database(monkeypatch, tmp_path):
    def _raise_files(_pkg):
        raise ModuleNotFoundError("no resources")

    monkeypatch.setattr("openmc_fusion_benchmarks.database.files", _raise_files)

    monkeypatch.setattr(Path, "exists", lambda _self: False)

    with pytest.raises(FileNotFoundError, match="results_database"):
        _resolve_database_path("missing", "experiment.h5")
