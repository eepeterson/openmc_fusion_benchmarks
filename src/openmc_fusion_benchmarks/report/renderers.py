from __future__ import annotations

from pathlib import Path

import yaml
import numpy as np

from .models import Report
from .plots import (
    build_plot_artifacts,
    build_quality_plot_artifacts,
    quality_metric_title,
    render_plots,
    render_quality_plots,
)
from ..validation.adapters import compare_tallies


def render_yaml(report: Report, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(report.data, handle, sort_keys=False)
    return output_path


def render_plots_for_report(report: Report, output_dir: Path) -> list[dict]:
    output_dir.mkdir(parents=True, exist_ok=True)
    entries: list[dict] = []
    quality_entries: list[dict] = []
    verbosity = int(report.data.get("verbosity", 0) or 0)
    quality_metrics = _quality_metrics_for_verbosity(verbosity)
    for plot in report.plots:
        artifacts = build_plot_artifacts(
            tally_name=plot.tally_name,
            experiment=plot.reference.get_tally(plot.tally_name),
            calculation=plot.candidate.get_tally(plot.tally_name),
            output_dir=output_dir,
        )
        render_plots(
            artifacts,
            plot.reference.get_tally(plot.tally_name),
            plot.candidate.get_tally(plot.tally_name),
            style=plot.style,
        )
        entries.append(
            {
                "tally": plot.tally_name,
                "absolute_plot": str(artifacts.absolute_plot),
                "ce_plot": str(artifacts.ce_plot),
            }
        )

        if quality_metrics:
            quality_artifacts = build_quality_plot_artifacts(
                tally_name=plot.tally_name,
                output_dir=output_dir,
                metrics=quality_metrics,
            )
            render_quality_plots(
                quality_artifacts,
                plot.reference.get_tally(plot.tally_name),
                plot.candidate.get_tally(plot.tally_name),
                metrics=quality_metrics,
                style=plot.style,
            )
            quality_entries.append(
                {
                    "tally": plot.tally_name,
                    "metrics": {k: str(v) for k, v in quality_artifacts.metric_plots.items()},
                }
            )
    report.data["plots"] = entries
    report.data["quality_plots"] = quality_entries
    return entries


def render_pdf(report: Report, output_path: Path, plots_dir: Path) -> Path:
    try:
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for PDF rendering. Install it to enable PDF output.") from exc

    plot_entries = render_plots_for_report(report, plots_dir)
    observable_entries = _collect_observable_metrics(report)
    if observable_entries:
        report.data["observable_summary"] = observable_entries
    quality_entries = report.data.get("quality_plots", [])

    with PdfPages(output_path) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.95, report.metadata.title, ha="center", fontsize=14)
        fig.text(0.1, 0.9, f"Benchmark: {report.metadata.benchmark_id}", fontsize=10)

        description = _resolve_benchmark_description(report)
        y_cursor = _wrap_text(fig, 0.1, 0.86, description, fontsize=9, width_chars=80)
        y_cursor = _wrap_text(fig, 0.1, y_cursor - 0.02, report.metadata.model_description, fontsize=9, width_chars=80)

        reference_line = _format_reference(report)
        if reference_line:
            fig.text(0.1, y_cursor - 0.02, f"Reference: {reference_line}", fontsize=9)
            y_cursor -= 0.06

        validation_line = _format_validation(report)
        if validation_line:
            fig.text(0.1, y_cursor - 0.02, f"Validation case: {validation_line}", fontsize=9)
            y_cursor -= 0.06

        _wrap_text(fig, 0.1, y_cursor - 0.02, report.metadata.notes, fontsize=9, width_chars=80)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        verbosity = int(report.data.get("verbosity", 0) or 0)
        if verbosity > 0:
            _render_specifications(pdf, report, verbosity)

        for entry in plot_entries:
            fig = plt.figure(figsize=(8.5, 11))
            fig.text(0.5, 0.95, f"Tally: {entry['tally']}", ha="center", fontsize=12)
            try:
                abs_img = plt.imread(entry["absolute_plot"])
                ce_img = plt.imread(entry["ce_plot"])
                ax1 = fig.add_axes([0.1, 0.55, 0.8, 0.35])
                ax1.imshow(abs_img)
                ax1.axis("off")
                ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.35])
                ax2.imshow(ce_img)
                ax2.axis("off")
            except Exception:
                fig.text(0.1, 0.6, "Plot images could not be loaded.", fontsize=9)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        if observable_entries:
            _render_observable_summary(pdf, observable_entries, verbosity)

        if quality_entries:
            _render_quality_section(pdf, quality_entries, verbosity)

    return output_path


def _wrap_text(fig, x: float, y: float, text: str, fontsize: int, width_chars: int) -> float:
    if not text:
        return y

    lines = _wrap_lines(text, width_chars)
    line_height = 0.02
    for idx, line in enumerate(lines):
        fig.text(x, y - idx * line_height, line, fontsize=fontsize)
    return y - len(lines) * line_height


def _wrap_lines(text: str, width_chars: int) -> list[str]:
    if not text:
        return []

    import textwrap

    return textwrap.wrap(text, width=width_chars, break_long_words=True, break_on_hyphens=False)


def _resolve_benchmark_description(report: Report) -> str:
    spec = report.data.get("specifications", {})
    spec_meta = spec.get("metadata", {}) if isinstance(spec, dict) else {}
    return spec_meta.get("description") or report.metadata.description


def _format_reference(report: Report) -> str:
    for entry in report.data.get("sources", []):
        if entry.get("kind") != "experiment":
            continue
        run_meta = entry.get("run_metadata", {})
        if run_meta.get("kind") == "experiment":
            return "experiment"
        return _format_run_metadata(run_meta)
    return ""


def _format_validation(report: Report) -> str:
    for entry in report.data.get("sources", []):
        if entry.get("kind") != "calculation":
            continue
        run_meta = entry.get("run_metadata", {})
        return _format_run_metadata(run_meta)
    return ""


def _format_run_metadata(run_meta: dict) -> str:
    if not run_meta:
        return ""

    code = " ".join([run_meta.get("code_name", ""), run_meta.get("code_version", "")]).strip()
    library = " ".join(
        [run_meta.get("nuclear_data_name", ""), run_meta.get("nuclear_data_version", "")]
    ).strip()

    parts = [p for p in (code, library) if p]
    return " | ".join(parts)


def _render_specifications(pdf, report: Report, verbosity: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    spec = report.data.get("specifications", {})
    if not isinstance(spec, dict) or not spec:
        return

    sections = _build_spec_sections(spec, verbosity)
    if not sections:
        return

    fig = plt.figure(figsize=(8.5, 11))
    y_cursor = 0.95
    fig.text(0.5, y_cursor, "Specification Summary", ha="center", fontsize=13)
    y_cursor -= 0.04

    for title, body in sections:
        body_lines = _wrap_lines(body, width_chars=100)
        needed = 0.03 + (len(body_lines) * 0.02) + 0.03
        if y_cursor - needed < 0.12:
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            fig = plt.figure(figsize=(8.5, 11))
            y_cursor = 0.95

        fig.text(0.1, y_cursor, title, fontsize=10, weight="bold")
        y_cursor -= 0.03
        y_cursor = _wrap_text(fig, 0.1, y_cursor, body, fontsize=8, width_chars=100)
        y_cursor -= 0.03

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _quality_metrics_for_verbosity(verbosity: int) -> list[str]:
    verbosity = max(0, min(int(verbosity), 3))
    if verbosity == 0:
        return ["ce"]
    if verbosity == 1:
        return ["ce", "chi2_contribution"]
    if verbosity == 2:
        return ["ce", "relative_deviation", "combined_uncertainty", "chi2_contribution"]
    return [
        "ce",
        "relative_deviation",
        "absolute_deviation",
        "combined_uncertainty",
        "normalized_residual",
        "chi2_contribution",
    ]


def _observable_metrics_for_verbosity(verbosity: int) -> list[str]:
    verbosity = max(0, min(int(verbosity), 3))
    if verbosity == 0:
        return ["rms_relative_deviation"]
    if verbosity == 1:
        return ["rms_relative_deviation", "reduced_chi2"]
    if verbosity == 2:
        return [
            "mean_bias",
            "mean_abs_relative_deviation",
            "rms_relative_deviation",
            "reduced_chi2",
        ]
    return [
        "mean_bias",
        "mean_abs_relative_deviation",
        "rms_relative_deviation",
        "mean_abs_normalized_residual",
        "reduced_chi2",
    ]


def _collect_observable_metrics(report: Report) -> list[dict]:
    entries: list[dict] = []
    for plot in report.plots:
        try:
            obs = compare_tallies(
                observable_name=plot.tally_name,
                observable_type="tally",
                reference=plot.reference.get_tally(plot.tally_name),
                candidate=plot.candidate.get_tally(plot.tally_name),
                include_grading=False,
            )
        except Exception:
            continue

        entries.append(
            {
                "tally": plot.tally_name,
                "mean_bias": obs.mean_bias,
                "mean_abs_relative_deviation": obs.mean_abs_relative_deviation,
                "rms_relative_deviation": obs.rms_relative_deviation,
                "mean_abs_normalized_residual": obs.mean_abs_normalized_residual,
                "reduced_chi2": obs.reduced_chi2,
            }
        )

    return entries


def _quality_metric_equation(metric: str) -> str:
    equations = {
        "ce": r"$\frac{C}{E}$",
        "relative_deviation": r"$\frac{C - E}{E}$",
        "absolute_deviation": r"$|C - E|$",
        "combined_uncertainty": r"$\sqrt{u_E^2 + u_C^2}$",
        "normalized_residual": r"$\frac{C - E}{\sqrt{u_E^2 + u_C^2}}$",
        "chi2_contribution": r"$\frac{(C - E)^2}{u_E^2 + u_C^2}$",
    }
    return equations.get(metric, "")


def _quality_metric_description(metric: str) -> str:
    descriptions = {
        "ce": "Ratio of calculation to experiment (ideal = 1).",
        "relative_deviation": "Signed fractional bias relative to experiment.",
        "absolute_deviation": "Absolute difference between calculation and experiment.",
        "combined_uncertainty": "Combined statistical uncertainty of C and E.",
        "normalized_residual": "Deviation in units of combined uncertainty.",
        "chi2_contribution": "Pointwise contribution to chi2 misfit.",
    }
    return descriptions.get(metric, "")


def _render_quality_section(pdf, quality_entries: list[dict], verbosity: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    metrics = _quality_metrics_for_verbosity(verbosity)
    if not metrics:
        return

    for entry in quality_entries:
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.95, f"Quality evaluation - Tally: {entry['tally']}", ha="center", fontsize=12)

        metric_paths = entry.get("metrics", {})
        cols = 2 if len(metrics) > 1 else 1
        rows = int(np.ceil(len(metrics) / cols))
        top = 0.88
        bottom = 0.08
        left = 0.08
        right = 0.92
        v_gap = 0.04
        h_gap = 0.06

        cell_h = (top - bottom - (rows - 1) * v_gap) / rows
        cell_w = (right - left - (cols - 1) * h_gap) / cols

        for idx, metric in enumerate(metrics):
            row = idx // cols
            col = idx % cols
            x0 = left + col * (cell_w + h_gap)
            y0 = top - (row + 1) * cell_h - row * v_gap
            ax = fig.add_axes([x0, y0, cell_w, cell_h])
            path = metric_paths.get(metric)
            if path:
                try:
                    img = plt.imread(path)
                    ax.imshow(img)
                except Exception:
                    ax.text(0.5, 0.5, "Plot image could not be loaded.", ha="center", va="center", fontsize=8)
            else:
                ax.text(0.5, 0.5, "Plot missing.", ha="center", va="center", fontsize=8)
            title = quality_metric_title(metric)
            equation = _quality_metric_equation(metric)
            if equation:
                ax.text(0.5, 1.05, equation, ha="center", va="bottom", fontsize=8, transform=ax.transAxes)
            description = _quality_metric_description(metric)
            if description:
                ax.text(0.5, 1.0, description, ha="center", va="bottom", fontsize=7, transform=ax.transAxes)
            ax.set_title(title, fontsize=9, pad=24)
            ax.axis("off")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def _render_observable_summary(pdf, observable_entries: list[dict], verbosity: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    metrics = _observable_metrics_for_verbosity(verbosity)
    if not metrics or not observable_entries:
        return

    tallies = [entry["tally"] for entry in observable_entries]
    cols = 2 if len(metrics) > 1 else 1
    rows = int(np.ceil(len(metrics) / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(8.5, 11))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        values = [entry.get(metric, np.nan) for entry in observable_entries]
        ax.bar(range(len(tallies)), values, color="#4C72B0")
        ax.set_title(metric.replace("_", " "))
        ax.set_xticks(range(len(tallies)))
        ax.set_xticklabels(tallies, rotation=45, ha="right", fontsize=8)
        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.6)

    for idx in range(len(metrics), len(axes)):
        axes[idx].axis("off")

    fig.suptitle("Observable summary", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)


def _build_spec_sections(spec: dict, verbosity: int) -> list[tuple[str, str]]:
    sections: list[tuple[str, str]] = []
    meta = spec.get("metadata", {})

    if meta:
        fields = {k: v for k, v in meta.items() if k not in {"geometry", "settings"}}
        sections.append(("Metadata", _format_inline(fields)))

    materials = spec.get("materials", [])
    if materials:
        if verbosity >= 3:
            sections.append(("Materials", _format_inline(materials)))
        elif verbosity == 2:
            sections.append((
                "Materials",
                _format_key_fields(
                    materials,
                    keys=["id", "name", "density", "composition"],
                    max_items=10,
                ),
            ))
        else:
            names = [m.get("name", "") for m in materials if isinstance(m, dict)]
            sections.append(("Materials", f"count={len(materials)}; names={', '.join(n for n in names if n)}"))

    tallies = spec.get("tallies", [])
    if tallies:
        if verbosity >= 3:
            sections.append(("Tallies", _format_inline(tallies)))
        elif verbosity == 2:
            sections.append((
                "Tallies",
                _format_key_fields(
                    tallies,
                    keys=["name", "particle", "scores", "filters"],
                    max_items=12,
                ),
            ))
        else:
            names = [t.get("name", "") for t in tallies if isinstance(t, dict)]
            sections.append(("Tallies", ", ".join(n for n in names if n)))

    geometry = spec.get("geometry")
    if geometry is not None:
        if verbosity >= 3:
            sections.append(("Geometry", _format_inline(geometry)))
        elif verbosity == 2:
            sections.append(("Geometry", _format_key_fields([geometry], keys=["cad_file", "meshing"], max_items=1)))

    settings = spec.get("settings")
    if settings is not None:
        if verbosity >= 3:
            sections.append(("Settings", _format_inline(settings)))
        elif verbosity == 2:
            sections.append((
                "Settings",
                _format_key_fields(
                    [settings],
                    keys=["run_mode", "batches", "particles_per_batch", "photon_transport"],
                    max_items=1,
                ),
            ))

    return sections


def _format_inline(value: object) -> str:
    import json

    try:
        return json.dumps(value, ensure_ascii=True, separators=(",", ":"))
    except TypeError:
        return json.dumps(str(value), ensure_ascii=True)


def _format_key_fields(items: list, keys: list[str], max_items: int) -> str:
    rows = []
    for item in items[:max_items]:
        if not isinstance(item, dict):
            rows.append(_format_inline(item))
            continue
        row = {k: item.get(k) for k in keys if k in item}
        rows.append(_format_inline(row))
    if len(items) > max_items:
        rows.append(f"... ({len(items) - max_items} more)")
    return " | ".join(rows)
