#!/usr/bin/env python3
"""
Build preprocessing/search runtime line charts from SSM-GED summary results.

The script reads the same summary.tsv inputs as make_average_runtime_tables.py.
For current summaries that only contain *_run_ms and *_output columns, it falls
back to parsing timing lines from each *_output file.

By default it writes one HTML file:

  - phase_runtime_by_t.html

The report has two chart sections:

  - preprocessing time, default metric: preprocessing_ms
  - search time, default metric: search_ms
"""

from __future__ import annotations

import argparse
import html
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from make_average_runtime_tables import (
    chart_segments,
    dataset_key_for_row,
    dataset_sort_key,
    dataset_title,
    parse_columns,
)
from make_result_tables import (
    CHART_COLORS,
    DEFAULT_CHART_THRESHOLDS,
    DEFAULT_TIMEOUT_SECONDS,
    KNOWN_ALGORITHM_SUFFIXES,
    KNOWN_META_PREFIXES,
    choose_output_dir,
    find_input_files,
    format_axis_runtime,
    format_compact_number,
    natural_key,
    parse_algorithm_list,
    parse_chart_thresholds,
    read_rows,
    resolve_algorithm_list,
    runtime_value_for_chart,
    threshold_matches,
    to_float,
)


DEFAULT_OUTPUT_NAME = "phase_runtime_by_t.html"
DEFAULT_PREPROCESS_METRIC = "preprocessing_ms"
DEFAULT_SEARCH_METRIC = "search_ms"
SUMMARY_PREFIX = "SSM_GED_SUMMARY "
PHASE_PREFIX = "SSM_GED_PHASES "
MISSING_OUTPUT_VALUES = {"", "NA", "N/A", "NULL", "NONE"}
TIME_VALUE_RE = r"([-+]?\d+(?:\.\d+)?)"
OUTPUT_TIME_PATTERNS = (
    ("init_ms", re.compile(rf"^\s*Init Time:\s*{TIME_VALUE_RE}\s*ms\b")),
    ("search_ms", re.compile(rf"^\s*Search Time:\s*{TIME_VALUE_RE}\s*ms\b")),
    ("load_ms", re.compile(rf"^\s*Load Graphs Time:\s*{TIME_VALUE_RE}\s*ms\b")),
    ("run_ms", re.compile(rf"^\s*Run Time:\s*{TIME_VALUE_RE}\s*ms\b")),
)


@dataclass(frozen=True)
class PhaseSpec:
    key: str
    label: str
    metric: str
    count_timeout_as_value: bool


@dataclass
class SourceCounts:
    summary: int = 0
    output: int = 0
    timeout: int = 0
    missing: int = 0


AverageBuckets = Dict[str, Dict[Tuple[str, ...], Dict[str, Dict[str, Optional[float]]]]]
CountBuckets = Dict[str, Dict[Tuple[str, ...], Dict[str, Dict[str, int]]]]
SourceStats = Dict[str, SourceCounts]
OutputMetricCache = Dict[str, Optional[Dict[str, float]]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create preprocessing/search phase runtime line charts from "
            "SSM-GED summary.tsv files."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="summary.tsv files, CSV/TSV files, or directories containing summary.tsv",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default=None,
        help="output directory; default: the summary.tsv directory or ./result_report",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_NAME,
        help=f"output HTML file name or path; default: {DEFAULT_OUTPUT_NAME}",
    )
    parser.add_argument(
        "--algorithms",
        default=None,
        help=(
            "comma-separated algorithms to include; default: all algorithms with "
            "phase metrics or output columns"
        ),
    )
    parser.add_argument(
        "--preprocess-metric",
        default=DEFAULT_PREPROCESS_METRIC,
        help=f"preprocessing-stage metric; default: {DEFAULT_PREPROCESS_METRIC}",
    )
    parser.add_argument(
        "--search-metric",
        default=DEFAULT_SEARCH_METRIC,
        help=f"search-stage metric; default: {DEFAULT_SEARCH_METRIC}",
    )
    parser.add_argument(
        "--thresholds",
        default=",".join(DEFAULT_CHART_THRESHOLDS),
        help="comma-separated t values to chart; default: 0,1,2,3,4,5,6",
    )
    parser.add_argument(
        "--dataset-columns",
        default="dataset_group,dataset",
        help=(
            "comma-separated columns defining a dataset row; "
            "default: dataset_group,dataset"
        ),
    )
    parser.add_argument(
        "--exact-dataset-names",
        action="store_true",
        help=(
            "do not combine dataset names with suffixes _15, _30, _45, and _60 "
            "by prefix"
        ),
    )
    parser.add_argument(
        "--no-output-fallback",
        action="store_true",
        help="do not parse *_output files when a phase metric is missing from summary.tsv",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=DEFAULT_TIMEOUT_SECONDS,
        help=(
            "timeout assigned to timed-out search runs; "
            f"default: {format(DEFAULT_TIMEOUT_SECONDS, 'g')} seconds"
        ),
    )
    return parser.parse_args()


def output_path_for_args(args: argparse.Namespace, input_files: Sequence[Path]) -> Path:
    output = Path(args.output)
    if output.is_absolute() or output.parent != Path("."):
        return output
    return choose_output_dir(args, input_files) / output


def phase_specs(args: argparse.Namespace) -> List[PhaseSpec]:
    return [
        PhaseSpec(
            key="preprocess",
            label="Preprocessing",
            metric=args.preprocess_metric,
            count_timeout_as_value=False,
        ),
        PhaseSpec(
            key="search",
            label="Search",
            metric=args.search_metric,
            count_timeout_as_value=True,
        ),
    ]


def detect_phase_algorithms(
    fieldnames: Iterable[str],
    metrics: Sequence[str],
) -> List[str]:
    metric_suffixes = tuple(f"_{metric}" for metric in metrics) + ("_output",)
    candidates = set()

    for field in fieldnames:
        for suffix in metric_suffixes:
            if not field.endswith(suffix):
                continue
            prefix = field[: -len(suffix)]
            if not prefix or prefix in KNOWN_META_PREFIXES or prefix == "expected":
                continue
            candidates.add(prefix)

    # Keep compatibility with future summaries that add any known algorithm metric.
    for field in fieldnames:
        for suffix in KNOWN_ALGORITHM_SUFFIXES:
            if not field.endswith(suffix):
                continue
            prefix = field[: -len(suffix)]
            if prefix in candidates or not prefix or prefix in KNOWN_META_PREFIXES:
                continue
            if any(f"{prefix}_{metric}" in fieldnames for metric in metrics):
                candidates.add(prefix)

    return sorted(candidates, key=natural_key)


def parse_summary_metrics(line: str) -> Dict[str, float]:
    if not line.startswith((SUMMARY_PREFIX, PHASE_PREFIX)):
        return {}

    metrics: Dict[str, float] = {}
    for token in line.strip().split():
        if "=" not in token:
            continue
        key, raw_value = token.split("=", 1)
        value = to_float(raw_value)
        if value is not None:
            metrics[key] = value
    return metrics


def parse_output_time_metrics(line: str) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    for metric, pattern in OUTPUT_TIME_PATTERNS:
        match = pattern.match(line)
        if not match:
            continue
        value = to_float(match.group(1))
        if value is not None:
            metrics[metric] = value
    return metrics


def rewrite_absolute_path_to_workspace(path: Path) -> Optional[Path]:
    cwd = Path.cwd().resolve()
    parts = path.parts
    for index, part in enumerate(parts):
        if part != cwd.name:
            continue
        candidate = cwd.joinpath(*parts[index + 1 :])
        if candidate.is_file():
            return candidate
    return None


def resolve_output_path(row: Dict[str, str], algorithm: str) -> Optional[Path]:
    raw_value = str(row.get(f"{algorithm}_output", "")).strip()
    if raw_value.upper() in MISSING_OUTPUT_VALUES:
        return None

    raw_path = Path(raw_value)
    candidates: List[Path] = []
    if raw_path.is_absolute():
        candidates.append(raw_path)
        rewritten = rewrite_absolute_path_to_workspace(raw_path)
        if rewritten is not None:
            candidates.append(rewritten)
    else:
        source_file = str(row.get("source_file", "")).strip()
        if source_file:
            candidates.append(Path(source_file).parent / raw_path)
        candidates.append(Path.cwd() / raw_path)
        candidates.append(raw_path)

    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            return candidate
    return None


def parse_output_metrics(path: Path) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parsed = parse_output_time_metrics(line)
                if parsed:
                    metrics.update(parsed)
                    continue

                parsed = parse_summary_metrics(line)
                if parsed:
                    metrics.update(parsed)
    except OSError:
        return {}
    return metrics


def output_metrics_for_row(
    row: Dict[str, str],
    algorithm: str,
    cache: OutputMetricCache,
) -> Optional[Dict[str, float]]:
    raw_value = str(row.get(f"{algorithm}_output", "")).strip()
    if raw_value.upper() in MISSING_OUTPUT_VALUES:
        return None

    if raw_value in cache:
        return cache[raw_value]

    path = resolve_output_path(row, algorithm)
    if path is None:
        cache[raw_value] = None
        return None

    metrics = parse_output_metrics(path)
    cache[raw_value] = metrics
    return metrics


def direct_metric_value(row: Dict[str, str], algorithm: str, metric: str) -> Optional[float]:
    value = to_float(row.get(f"{algorithm}_{metric}", ""))
    if value is None or value < 0:
        return None
    return value


def phase_value_for_row(
    row: Dict[str, str],
    algorithm: str,
    phase: PhaseSpec,
    timeout_ms: float,
    use_output_fallback: bool,
    output_cache: OutputMetricCache,
) -> Tuple[Optional[float], str]:
    direct_value = direct_metric_value(row, algorithm, phase.metric)
    if direct_value is not None:
        if phase.count_timeout_as_value:
            direct_value = min(direct_value, timeout_ms)
        return direct_value, "summary"

    if use_output_fallback:
        output_metrics = output_metrics_for_row(row, algorithm, output_cache)
        if output_metrics is not None:
            output_value = output_metrics.get(phase.metric)
            if output_value is not None and output_value >= 0:
                if phase.count_timeout_as_value:
                    output_value = min(output_value, timeout_ms)
                return output_value, "output"

    if phase.count_timeout_as_value:
        timeout_value = runtime_value_for_chart(
            row,
            algorithm,
            phase.metric,
            timeout_ms,
        )
        if timeout_value is not None:
            return timeout_value, "timeout"

    return None, "missing"


def matched_threshold_for_row(
    row: Dict[str, str],
    thresholds: Sequence[str],
) -> Optional[str]:
    for threshold in thresholds:
        if threshold_matches(row.get("threshold", ""), threshold):
            return threshold
    return None


def aggregate_phase_runtimes(
    rows: Sequence[Dict[str, str]],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    dataset_columns: Sequence[str],
    combine_label_count_suffixes: bool,
    phases: Sequence[PhaseSpec],
    timeout_ms: float,
    use_output_fallback: bool,
) -> Tuple[List[Tuple[str, ...]], AverageBuckets, CountBuckets, SourceStats]:
    value_buckets: Dict[Tuple[str, Tuple[str, ...], str, str], List[float]] = defaultdict(list)
    dataset_keys = set()
    output_cache: OutputMetricCache = {}
    source_stats: SourceStats = {phase.key: SourceCounts() for phase in phases}

    for row in rows:
        threshold = matched_threshold_for_row(row, thresholds)
        if threshold is None:
            continue

        dataset_key = dataset_key_for_row(
            row,
            dataset_columns,
            combine_label_count_suffixes,
        )
        dataset_keys.add(dataset_key)

        for phase in phases:
            stats = source_stats[phase.key]
            for algorithm in algorithms:
                value, source = phase_value_for_row(
                    row,
                    algorithm,
                    phase,
                    timeout_ms,
                    use_output_fallback,
                    output_cache,
                )
                if source == "summary":
                    stats.summary += 1
                elif source == "output":
                    stats.output += 1
                elif source == "timeout":
                    stats.timeout += 1
                else:
                    stats.missing += 1

                if value is not None:
                    value_buckets[(phase.key, dataset_key, threshold, algorithm)].append(value)

    ordered_keys = sorted(dataset_keys, key=dataset_sort_key)
    averages: AverageBuckets = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    counts: CountBuckets = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    for phase in phases:
        for dataset_key in ordered_keys:
            for threshold in thresholds:
                for algorithm in algorithms:
                    values = value_buckets.get(
                        (phase.key, dataset_key, threshold, algorithm),
                        [],
                    )
                    if not values:
                        averages[phase.key][dataset_key][threshold][algorithm] = None
                        counts[phase.key][dataset_key][threshold][algorithm] = 0
                        continue
                    averages[phase.key][dataset_key][threshold][algorithm] = (
                        sum(values) / len(values)
                    )
                    counts[phase.key][dataset_key][threshold][algorithm] = len(values)

    return ordered_keys, averages, counts, source_stats


def build_phase_line_chart_svg(
    phase: PhaseSpec,
    dataset_key: Tuple[str, ...],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    averages: AverageBuckets,
    counts: CountBuckets,
) -> str:
    phase_averages = averages.get(phase.key, {})
    values = [
        value
        for threshold in thresholds
        for value in phase_averages.get(dataset_key, {}).get(threshold, {}).values()
        if value is not None and value >= 0
    ]
    title = dataset_title(dataset_key)
    if not values:
        return (
            '<div class="chart-card">'
            f'<h3>{html.escape(title)}</h3>'
            f'<div class="chart-empty">No {html.escape(phase.label.lower())} data '
            "for this dataset.</div>"
            "</div>"
        )

    positive_values = [value for value in values if value > 0]
    if not positive_values:
        return (
            '<div class="chart-card">'
            f'<h3>{html.escape(title)}</h3>'
            f'<div class="chart-empty">No positive {html.escape(phase.label.lower())} '
            "data for this dataset.</div>"
            "</div>"
        )

    min_value = min(positive_values)
    max_value = max(positive_values)
    min_exp = math.floor(math.log10(min_value))
    max_exp = math.ceil(math.log10(max_value))
    if min_exp == max_exp:
        min_exp -= 1
    while max_exp - min_exp > 8:
        min_exp += 1
    y_min = 10.0 ** min_exp
    y_max = 10.0 ** max_exp
    ticks = [10.0 ** exponent for exponent in range(min_exp, max_exp + 1)]

    left = 68
    right = 22
    top = 22
    bottom = 46
    plot_width = 620
    plot_height = 230
    width = left + plot_width + right
    height = top + plot_height + bottom
    axis_bottom = top + plot_height

    def x_position(index: int) -> float:
        if len(thresholds) == 1:
            return left + plot_width / 2.0
        return left + plot_width * index / (len(thresholds) - 1)

    def y_position(value: float) -> float:
        clipped = min(max(value, y_min), y_max)
        portion = (math.log10(clipped) - min_exp) / (max_exp - min_exp)
        return top + plot_height * (1.0 - portion)

    parts = [
        '<div class="chart-card">',
        f'<h3>{html.escape(title)}</h3>',
        '<div class="runtime-chart-legend">',
    ]
    for index, algorithm in enumerate(algorithms):
        color = CHART_COLORS[index % len(CHART_COLORS)]
        parts.append(
            '<span class="runtime-legend-item">'
            f'<i style="background:{color}"></i>{html.escape(algorithm)}</span>'
        )
    parts.extend(
        [
            "</div>",
            '<div class="chart-scroll">',
            f'<svg class="line-chart" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" '
            f'aria-label="Average {html.escape(phase.label.lower(), quote=True)} '
            f'time by t for {html.escape(title, quote=True)}">',
        ]
    )

    for tick in ticks:
        y = y_position(tick)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" '
            'stroke="#d1d5db" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" '
            f'class="axis-label">{html.escape(format_axis_runtime(tick))}</text>'
        )

    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_bottom}" '
        'stroke="#374151" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{left}" y1="{axis_bottom}" x2="{left + plot_width}" '
        f'y2="{axis_bottom}" stroke="#374151" stroke-width="1.2" />'
    )

    for index, threshold in enumerate(thresholds):
        x = x_position(index)
        parts.append(
            f'<line x1="{x:.2f}" y1="{axis_bottom}" x2="{x:.2f}" '
            f'y2="{axis_bottom + 4}" stroke="#374151" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{axis_bottom + 19}" text-anchor="middle" '
            f'class="axis-label">t={html.escape(threshold)}</text>'
        )

    for algorithm_index, algorithm in enumerate(algorithms):
        color = CHART_COLORS[algorithm_index % len(CHART_COLORS)]
        raw_points: List[Optional[Tuple[float, float]]] = []
        for threshold_index, threshold in enumerate(thresholds):
            value = phase_averages.get(dataset_key, {}).get(threshold, {}).get(algorithm)
            if value is None:
                raw_points.append(None)
            else:
                raw_points.append((x_position(threshold_index), y_position(value)))

        for segment in chart_segments(raw_points):
            if len(segment) == 1:
                continue
            point_text = " ".join(f"{x:.2f},{y:.2f}" for x, y in segment)
            parts.append(
                f'<polyline points="{point_text}" fill="none" stroke="{color}" '
                'stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" />'
            )

        for threshold_index, threshold in enumerate(thresholds):
            value = phase_averages.get(dataset_key, {}).get(threshold, {}).get(algorithm)
            if value is None:
                continue
            count = counts.get(phase.key, {}).get(dataset_key, {}).get(threshold, {}).get(
                algorithm,
                0,
            )
            x = x_position(threshold_index)
            y = y_position(value)
            tooltip = (
                f"{algorithm}; t={threshold}; {phase.label.lower()} "
                f"average={value:.4f} ms; n={count}"
            )
            parts.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}" '
                'stroke="#111827" stroke-width=".6">'
                f"<title>{html.escape(tooltip)}</title></circle>"
            )

    y_title_x = 16
    y_title_y = top + plot_height / 2.0
    parts.append(
        f'<text x="{y_title_x}" y="{y_title_y:.2f}" text-anchor="middle" '
        f'transform="rotate(-90 {y_title_x} {y_title_y:.2f})" '
        f'class="axis-title">Average {html.escape(phase.label.lower())} '
        "time (ms, log scale)</text>"
    )
    parts.append(
        f'<text x="{left + plot_width / 2.0:.2f}" y="{height - 10}" '
        'text-anchor="middle" class="axis-title">t</text>'
    )

    parts.extend(["</svg>", "</div>", "</div>"])
    return "".join(parts)


def source_summary_html(stats: SourceCounts) -> str:
    return (
        f"<span>summary cells: {stats.summary}</span>"
        f"<span>output-file cells: {stats.output}</span>"
        f"<span>timeout cells: {stats.timeout}</span>"
        f"<span>missing cells: {stats.missing}</span>"
    )


def write_html_report(
    path: Path,
    dataset_keys: Sequence[Tuple[str, ...]],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    phases: Sequence[PhaseSpec],
    averages: AverageBuckets,
    counts: CountBuckets,
    source_stats: SourceStats,
    timeout_seconds: float,
    combine_label_count_suffixes: bool,
    use_output_fallback: bool,
) -> None:
    style = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 24px; color: #111827; }
h1 { font-size: 24px; margin: 0 0 8px; }
h2 { font-size: 18px; margin: 28px 0 8px; }
h3 { font-size: 14px; margin: 0 0 7px; }
p { color: #4b5563; margin: 0 0 18px; }
.summary { display: flex; gap: 10px; flex-wrap: wrap; margin: 12px 0 18px; font-size: 12px; color: #374151; }
.summary span { background: #f9fafb; border: 1px solid #e5e7eb; padding: 4px 8px; }
.chart-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 620px), 1fr)); gap: 12px; }
.chart-card { border: 1px solid #d1d5db; border-radius: 6px; background: #fff; padding: 8px 9px 9px; min-width: 0; }
.chart-scroll { overflow-x: auto; }
.line-chart { display: block; max-width: none; }
.runtime-chart-legend { display: flex; flex-wrap: wrap; gap: 5px 10px; margin: 0 0 5px; font-size: 11px; }
.runtime-legend-item { display: inline-flex; align-items: center; gap: 5px; }
.runtime-legend-item i { display: inline-block; width: 11px; height: 11px; border: 1px solid rgba(17,24,39,.35); }
.axis-label { fill: #374151; font-size: 10px; }
.axis-title { fill: #111827; font-size: 11px; font-weight: 600; }
.chart-empty { color: #6b7280; font-size: 12px; padding: 12px 4px; }
"""
    with path.open("w", encoding="utf-8") as handle:
        handle.write("<!doctype html>\n<html><head><meta charset=\"utf-8\">\n")
        handle.write("<title>SSM-GED Phase Runtime by t</title>\n")
        handle.write(f"<style>{style}</style>\n</head><body>\n")
        handle.write("<h1>SSM-GED Phase Runtime by t</h1>\n")
        suffix_note = (
            " Dataset suffixes _15, _30, _45, and _60 are combined by prefix."
            if combine_label_count_suffixes
            else ""
        )
        fallback_note = (
            " Missing phase metrics are parsed from *_output files when possible."
            if use_output_fallback
            else " Output-file fallback is disabled."
        )
        metric_note = "; ".join(f"{phase.label}: {phase.metric}" for phase in phases)
        handle.write(
            f"<p>Metrics: {html.escape(metric_note)}. Each point is the arithmetic "
            "mean over all rows for the same dataset, t, phase, and algorithm. "
            f"Timed-out search runs and NA search runtime cells are counted as "
            f"{html.escape(format_compact_number(timeout_seconds))} seconds."
            f"{html.escape(fallback_note)}{html.escape(suffix_note)}</p>\n"
        )
        handle.write(
            '<div class="summary">'
            f"<span>datasets: {len(dataset_keys)}</span>"
            f"<span>algorithms: {html.escape(', '.join(algorithms))}</span>"
            f"<span>t values: {html.escape(', '.join(thresholds))}</span>"
            "</div>\n"
        )

        for phase in phases:
            handle.write(f"<h2>{html.escape(phase.label)} Time Lines</h2>\n")
            stats = source_stats.get(phase.key, SourceCounts())
            handle.write(f'<div class="summary">{source_summary_html(stats)}</div>\n')
            handle.write('<div class="chart-grid">\n')
            for key in dataset_keys:
                handle.write(
                    build_phase_line_chart_svg(
                        phase,
                        key,
                        algorithms,
                        thresholds,
                        averages,
                        counts,
                    )
                )
                handle.write("\n")
            handle.write("</div>\n")

        handle.write("</body></html>\n")


def validate_metric_name(value: str, option_name: str) -> None:
    if not value or value.strip() != value or " " in value:
        raise ValueError(f"{option_name} must be a simple column suffix such as run_ms")


def main() -> int:
    args = parse_args()
    try:
        validate_metric_name(args.preprocess_metric, "--preprocess-metric")
        validate_metric_name(args.search_metric, "--search-metric")

        input_files = find_input_files(args.inputs)
        rows, fieldnames = read_rows(input_files)
        if not rows:
            raise ValueError("No result rows found.")

        thresholds = parse_chart_thresholds(args.thresholds)
        if args.timeout_seconds <= 0:
            raise ValueError("--timeout-seconds must be greater than zero")
        timeout_ms = args.timeout_seconds * 1000.0

        dataset_columns = parse_columns(args.dataset_columns, "--dataset-columns")
        missing_dataset_columns = [
            column for column in dataset_columns if column not in fieldnames
        ]
        if missing_dataset_columns:
            missing = ", ".join(missing_dataset_columns)
            raise ValueError(f"Missing --dataset-columns column(s): {missing}")

        phases = phase_specs(args)
        detected_algorithms = detect_phase_algorithms(
            fieldnames,
            [phase.metric for phase in phases],
        )
        if not detected_algorithms:
            raise ValueError(
                "No algorithms detected from phase metric or *_output columns."
            )

        selected_algorithm_names = parse_algorithm_list(args.algorithms)
        if selected_algorithm_names is None:
            algorithms = detected_algorithms
        else:
            algorithms = resolve_algorithm_list(
                selected_algorithm_names,
                detected_algorithms,
            )

        combine_label_count_suffixes = not args.exact_dataset_names
        use_output_fallback = not args.no_output_fallback
        dataset_keys, averages, counts, source_stats = aggregate_phase_runtimes(
            rows,
            algorithms,
            thresholds,
            dataset_columns,
            combine_label_count_suffixes,
            phases,
            timeout_ms,
            use_output_fallback,
        )
        if not dataset_keys:
            raise ValueError(
                "No rows matched the selected t values: " + ", ".join(thresholds)
            )

        output_path = output_path_for_args(args, input_files)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        write_html_report(
            output_path,
            dataset_keys,
            algorithms,
            thresholds,
            phases,
            averages,
            counts,
            source_stats,
            args.timeout_seconds,
            combine_label_count_suffixes,
            use_output_fallback,
        )

        print(f"Input files: {len(input_files)}")
        print(f"Rows: {len(rows)}")
        print(f"Datasets: {len(dataset_keys)}")
        print(f"Algorithms: {', '.join(algorithms)}")
        for phase in phases:
            stats = source_stats[phase.key]
            print(
                f"{phase.label} metric: {phase.metric} "
                f"(summary={stats.summary}, output={stats.output}, "
                f"timeout={stats.timeout}, missing={stats.missing})"
            )
        print(f"t values: {', '.join(thresholds)}")
        print(f"Search timeout: {args.timeout_seconds:g} seconds")
        print(f"Output HTML: {output_path}")
        return 0
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
