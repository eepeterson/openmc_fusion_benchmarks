from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import xarray as xr

from .models import PlotSpec, Report, ReportConfig, ReportMetadata, ResultSource


def _resolve_plot_tallies(
    sources: Sequence[ResultSource],
    plot_tallies: Sequence[str] | None,
) -> list[str]:
    if plot_tallies is not None:
        return list(plot_tallies)

    for source in sources:
        if source.tally_names:
            return list(source.tally_names)

    if sources:
        return _list_report_tallies(sources[0].results)
    return []


def _list_report_tallies(results) -> list[str]:
    tallies: list[str] = []
    for group in results.tallies:
        try:
            with xr.open_dataset(results.filepath, group=group, engine="h5netcdf") as ds:
                if "mean" in ds:
                    tallies.append(group)
        except Exception:
            continue
    return tallies


def _derive_metadata(
    spec_snapshot: dict | None,
    run_metadata: dict | None,
) -> ReportMetadata:
    spec_meta = (spec_snapshot or {}).get("metadata", {})
    title = spec_meta.get("title") or "Benchmark report"
    benchmark_id = spec_meta.get("id") or spec_meta.get("title") or "unknown"

    description = spec_meta.get("description", "")
    model_description = spec_meta.get("model_description", "")

    code_name = (run_metadata or {}).get("code_name", "")
    code_version = (run_metadata or {}).get("code_version", "")

    return ReportMetadata(
        title=title,
        benchmark_id=benchmark_id,
        description=description,
        model_description=model_description,
        code_name=code_name,
        code_version=code_version,
    )


def build_report(
    metadata: ReportMetadata | Sequence[ResultSource],
    sources: Sequence[ResultSource] | ReportConfig,
    config: ReportConfig | None = None,
) -> Report:
    if isinstance(metadata, ReportMetadata):
        report_metadata = metadata
        report_sources = sources
        report_config = config
        if report_config is None:
            raise ValueError("ReportConfig is required when metadata is provided")
    else:
        report_metadata = None
        report_sources = metadata
        report_config = sources

    if report_config is None:
        raise ValueError("ReportConfig is required")

    source_entries = []
    spec_snapshot = None
    run_metadata = None

    for source in report_sources:
        tally_list = list(source.tally_names) if source.tally_names else _list_report_tallies(source.results)
        entry = {
            "name": source.name,
            "kind": source.kind,
            "file": str(source.results.filepath),
            "tallies": tally_list,
        }
        try:
            entry["specifications"] = source.results.get_spec_snapshot()
            if spec_snapshot is None:
                spec_snapshot = entry["specifications"]
        except Exception:
            pass
        try:
            entry["run_metadata"] = source.results.get_run_metadata()
            if run_metadata is None:
                run_metadata = entry["run_metadata"]
        except Exception:
            pass

        source_entries.append(entry)

    if report_metadata is None:
        report_metadata = _derive_metadata(spec_snapshot, run_metadata)

    data = {
        "title": report_metadata.title,
        "benchmark_id": report_metadata.benchmark_id,
        "description": report_metadata.description,
        "model_description": report_metadata.model_description,
        "code_name": report_metadata.code_name,
        "code_version": report_metadata.code_version,
        "notes": report_metadata.notes,
        "verbosity": report_config.verbosity,
        "sources": source_entries,
    }

    if spec_snapshot is not None:
        data["specifications"] = spec_snapshot
    if run_metadata is not None:
        data["run_metadata"] = run_metadata

    plot_tallies = _resolve_plot_tallies(report_sources, report_config.plot_tallies)
    plots: list[PlotSpec] = []

    reference_source = next((s for s in report_sources if s.kind == "experiment"), None)
    candidate_source = next((s for s in report_sources if s.kind == "calculation"), None)

    if reference_source is not None and candidate_source is not None:
        for tally_name in plot_tallies:
            plots.append(
                PlotSpec(
                    tally_name=tally_name,
                    experiment=reference_source.results,
                    calculation=candidate_source.results,
                )
            )

    return Report(metadata=report_metadata, sources=list(report_sources), plots=plots, data=data)
