import importlib.util
import sys
import types
from pathlib import Path

import h5py
import pytest
import yaml


def _module_available(name: str) -> bool:
    if name in sys.modules:
        return True
    try:
        return importlib.util.find_spec(name) is not None
    except ValueError:
        return False


# Provide minimal stubs so package imports work in minimal test envs.
if not _module_available("openmc"):
    openmc_stub = types.ModuleType("openmc")
    openmc_stub.__path__ = []
    sys.modules.setdefault("openmc", openmc_stub)

if not _module_available("sandy"):
    sys.modules.setdefault("sandy", types.ModuleType("sandy"))


from scripts.backfill_results_metadata import main as backfill_main


def _write_empty_results(filepath: Path) -> None:
    with h5py.File(filepath, "w") as handle:
        handle.create_group("tally_1")


def test_backfill_experiment(tmp_path, monkeypatch):
    results_path = tmp_path / "experiment.h5"
    spec_path = tmp_path / "specifications.yaml"

    _write_empty_results(results_path)
    spec_path.write_text(yaml.safe_dump({"metadata": {"title": "Test"}}), encoding="utf-8")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "backfill_results_metadata.py",
            "--file",
            str(results_path),
            "--spec-file",
            str(spec_path),
            "--kind",
            "experiment",
            "--benchmark-name",
            "test",
        ],
    )

    backfill_main()

    with h5py.File(results_path, "r") as handle:
        assert "specifications" in handle
        assert "run_metadata" in handle
        spec_group = handle["specifications"]
        assert spec_group.attrs["benchmark_name"] == "test"
        run_meta = handle["run_metadata"].attrs
        assert run_meta["kind"] == "experiment"


def test_backfill_calculation(tmp_path, monkeypatch):
    results_path = tmp_path / "calc.h5"
    spec_path = tmp_path / "specifications.yaml"

    _write_empty_results(results_path)
    spec_path.write_text(yaml.safe_dump({"metadata": {"title": "Test"}}), encoding="utf-8")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "backfill_results_metadata.py",
            "--file",
            str(results_path),
            "--spec-file",
            str(spec_path),
            "--kind",
            "calculation",
            "--code-name",
            "tripoli",
            "--code-version",
            "4.11",
            "--nuclear-data-name",
            "FENDL",
            "--nuclear-data-version",
            "3.2b",
            "--geometry",
            "cad",
        ],
    )

    backfill_main()

    with h5py.File(results_path, "r") as handle:
        assert "specifications" in handle
        assert "run_metadata" in handle
        meta = handle["run_metadata"].attrs
        assert meta["code_name"] == "tripoli"
        assert meta["code_version"] == "4.11"
        assert meta["nuclear_data_name"] == "FENDL"
        assert meta["nuclear_data_version"] == "3.2b"
        assert meta["geometry"] == "cad"
