#!/usr/bin/env python3
"""Backfill specifications and run metadata into existing results files."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import yaml


def _write_spec_snapshot(handle: h5py.File, spec_path: Path, benchmark_name: str | None) -> None:
    spec_yaml = spec_path.read_text(encoding="utf-8")
    spec_bytes = spec_yaml.encode("utf-8")

    if "specifications" in handle:
        del handle["specifications"]
    group = handle.create_group("specifications")
    group.attrs["format"] = "yaml"
    if benchmark_name:
        group.attrs["benchmark_name"] = benchmark_name
    group.create_dataset("yaml", data=np.bytes_(spec_bytes))


def _write_run_metadata(
    handle: h5py.File,
    code_name: str | None,
    code_version: str | None,
    nuclear_data_name: str | None,
    nuclear_data_version: str | None,
    geometry: str | None,
) -> None:
    if "run_metadata" in handle:
        del handle["run_metadata"]
    group = handle.create_group("run_metadata")

    if code_name:
        group.attrs["code_name"] = code_name
    if code_version:
        group.attrs["code_version"] = code_version
    if nuclear_data_name:
        group.attrs["nuclear_data_name"] = nuclear_data_name
    if nuclear_data_version:
        group.attrs["nuclear_data_version"] = nuclear_data_version
    if geometry:
        group.attrs["geometry"] = geometry


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Backfill specifications and run metadata into an OFB results file."
    )
    parser.add_argument("--file", type=Path, required=True, help="Path to results .h5 file")
    parser.add_argument("--spec-file", type=Path, required=True, help="Path to specifications.yaml")
    parser.add_argument(
        "--kind",
        choices=["experiment", "calculation"],
        required=True,
        help="Result kind (controls which metadata are written)",
    )
    parser.add_argument("--benchmark-name", type=str, default=None, help="Benchmark name label")
    parser.add_argument("--code-name", type=str, default=None, help="Code name (calculation only)")
    parser.add_argument("--code-version", type=str, default=None, help="Code version (calculation only)")
    parser.add_argument(
        "--nuclear-data-name",
        type=str,
        default=None,
        help="Nuclear data library name (calculation only)",
    )
    parser.add_argument(
        "--nuclear-data-version",
        type=str,
        default=None,
        help="Nuclear data library version (calculation only)",
    )
    parser.add_argument(
        "--geometry",
        type=str,
        default=None,
        choices=["cad", "csg"],
        help="Geometry type (calculation only)",
    )
    args = parser.parse_args()

    if not args.file.exists():
        raise FileNotFoundError(f"Results file not found: {args.file}")
    if not args.spec_file.exists():
        raise FileNotFoundError(f"Spec file not found: {args.spec_file}")

    # Validate the YAML early to fail fast if malformed.
    with args.spec_file.open("r", encoding="utf-8") as handle:
        yaml.safe_load(handle)

    with h5py.File(args.file, "a") as handle:
        _write_spec_snapshot(handle, args.spec_file, args.benchmark_name)
        if args.kind == "calculation":
            _write_run_metadata(
                handle,
                code_name=args.code_name,
                code_version=args.code_version,
                nuclear_data_name=args.nuclear_data_name,
                nuclear_data_version=args.nuclear_data_version,
                geometry=args.geometry,
            )
        else:
            if "run_metadata" in handle:
                del handle["run_metadata"]
            group = handle.create_group("run_metadata")
            group.attrs["kind"] = "experiment"

    print(f"Updated: {args.file}")


if __name__ == "__main__":
    main()
