#!/usr/bin/env python3
"""
Build an average-runtime HTML report from SSM-GED summary results.

The report is intentionally narrower than make_result_tables.py:

  - one output HTML file, named avg_runtime_by_t.html by default
  - one line chart per dataset, with t on the x axis and average runtime on
    a logarithmic y axis
  - no per-t summary tables are written to the HTML output

Runtime values are arithmetic means over all rows belonging to the same
dataset, threshold, and algorithm. Timeout-like values are handled the same way
as the runtime charts in make_result_tables.py.
"""

from __future__ import annotations

import argparse
import html
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from make_result_tables import (
    CHART_COLORS,
    DEFAULT_CHART_RUNTIME_METRIC,
    DEFAULT_CHART_THRESHOLDS,
    DEFAULT_TIMEOUT_SECONDS,
    chart_dataset_prefix,
    choose_output_dir,
    detect_algorithms,
    find_input_files,
    format_axis_runtime,
    format_compact_number,
    format_number,
    natural_key,
    parse_algorithm_list,
    parse_chart_thresholds,
    read_rows,
    resolve_algorithm_list,
    runtime_value_for_chart,
    threshold_matches,
)


DEFAULT_OUTPUT_NAME = "avg_runtime_by_t.html"
RESERVED_OUTPUT_NAMES = {"speedup_chart.html"}


AverageBuckets = Dict[Tuple[str, ...], Dict[str, Dict[str, Optional[float]]]]
CountBuckets = Dict[Tuple[str, ...], Dict[str, Dict[str, int]]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create dataset-level average-runtime line charts from "
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
            "the selected runtime metric"
        ),
    )
    parser.add_argument(
        "--runtime-metric",
        default=DEFAULT_CHART_RUNTIME_METRIC,
        help=f"runtime metric to average; default: {DEFAULT_CHART_RUNTIME_METRIC}",
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
        "--timeout-seconds",
        type=float,
        default=DEFAULT_TIMEOUT_SECONDS,
        help=(
            "timeout assigned to timed-out runs; "
            f"default: {format(DEFAULT_TIMEOUT_SECONDS, 'g')} seconds"
        ),
    )
    return parser.parse_args()


def parse_columns(value: str, option_name: str) -> List[str]:
    columns = [part.strip() for part in value.split(",") if part.strip()]
    if not columns:
        raise ValueError(f"{option_name} must contain at least one column")
    return columns


def output_path_for_args(args: argparse.Namespace, input_files: Sequence[Path]) -> Path:
    output = Path(args.output)
    if output.name in RESERVED_OUTPUT_NAMES:
        reserved = ", ".join(sorted(RESERVED_OUTPUT_NAMES))
        raise ValueError(f"--output must not reuse original chart output name(s): {reserved}")
    if output.is_absolute() or output.parent != Path("."):
        return output
    return choose_output_dir(args, input_files) / output


def dataset_key_for_row(
    row: Dict[str, str],
    dataset_columns: Sequence[str],
    combine_label_count_suffixes: bool,
) -> Tuple[str, ...]:
    key = []
    for column in dataset_columns:
        value = str(row.get(column, "")).strip()
        if combine_label_count_suffixes and column == "dataset":
            value = chart_dataset_prefix(value)
        key.append(value)
    return tuple(key)


def dataset_title(key: Tuple[str, ...]) -> str:
    return " / ".join(part or "<empty>" for part in key)


def dataset_sort_key(key: Tuple[str, ...]) -> Tuple[Tuple[object, ...], ...]:
    return tuple(natural_key(part) for part in key)


def aggregate_average_runtimes(
    rows: Sequence[Dict[str, str]],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    dataset_columns: Sequence[str],
    combine_label_count_suffixes: bool,
    runtime_metric: str,
    timeout_ms: float,
) -> Tuple[List[Tuple[str, ...]], AverageBuckets, CountBuckets]:
    threshold_set = set(thresholds)
    value_buckets: Dict[Tuple[Tuple[str, ...], str, str], List[float]] = defaultdict(list)
    dataset_keys = set()

    for row in rows:
        matched_threshold = None
        for threshold in thresholds:
            if threshold_matches(row.get("threshold", ""), threshold):
                matched_threshold = threshold
                break
        if matched_threshold is None:
            continue

        key = dataset_key_for_row(row, dataset_columns, combine_label_count_suffixes)
        dataset_keys.add(key)
        for algorithm in algorithms:
            runtime = runtime_value_for_chart(row, algorithm, runtime_metric, timeout_ms)
            if runtime is not None:
                value_buckets[(key, matched_threshold, algorithm)].append(runtime)

    ordered_keys = sorted(dataset_keys, key=dataset_sort_key)
    averages: AverageBuckets = defaultdict(lambda: defaultdict(dict))
    counts: CountBuckets = defaultdict(lambda: defaultdict(dict))

    for key in ordered_keys:
        for threshold in threshold_set:
            for algorithm in algorithms:
                values = value_buckets.get((key, threshold, algorithm), [])
                if not values:
                    averages[key][threshold][algorithm] = None
                    counts[key][threshold][algorithm] = 0
                    continue
                averages[key][threshold][algorithm] = sum(values) / len(values)
                counts[key][threshold][algorithm] = len(values)

    return ordered_keys, averages, counts


def table_columns(dataset_columns: Sequence[str], algorithms: Sequence[str]) -> List[str]:
    return list(dataset_columns) + list(algorithms)


def table_rows_for_threshold(
    dataset_keys: Sequence[Tuple[str, ...]],
    algorithms: Sequence[str],
    threshold: str,
    averages: AverageBuckets,
    counts: CountBuckets,
    dataset_columns: Sequence[str],
) -> List[Dict[str, str]]:
    rows = []
    for key in dataset_keys:
        output = {column: value for column, value in zip(dataset_columns, key)}
        for algorithm in algorithms:
            value = averages.get(key, {}).get(threshold, {}).get(algorithm)
            output[algorithm] = format_number(value) if value is not None else ""
            output[f"__title__{algorithm}"] = (
                f"n={counts.get(key, {}).get(threshold, {}).get(algorithm, 0)}"
            )
        rows.append(output)
    return rows


def is_text_column(column: str, dataset_columns: Sequence[str]) -> bool:
    return column in dataset_columns


def write_table_block(
    handle,
    columns: Sequence[str],
    rows: Sequence[Dict[str, str]],
    dataset_columns: Sequence[str],
) -> None:
    handle.write('<div class="table-block"><div class="table-wrap"><table>\n')
    handle.write("<thead><tr>")
    for column in columns:
        css_class = "text" if is_text_column(column, dataset_columns) else ""
        handle.write(f'<th class="{css_class}">{html.escape(column)}</th>')
    handle.write("</tr></thead>\n<tbody>\n")
    for row in rows:
        handle.write("<tr>")
        for column in columns:
            css_class = "text" if is_text_column(column, dataset_columns) else ""
            value = row.get(column, "")
            title = row.get(f"__title__{column}", "")
            title_attr = f' title="{html.escape(title)}"' if title else ""
            handle.write(
                f'<td class="{css_class}"{title_attr}>{html.escape(value)}</td>'
            )
        handle.write("</tr>\n")
    handle.write("</tbody></table></div></div>\n")


def chart_segments(points: Sequence[Optional[Tuple[float, float]]]) -> List[List[Tuple[float, float]]]:
    segments: List[List[Tuple[float, float]]] = []
    current: List[Tuple[float, float]] = []
    for point in points:
        if point is None:
            if current:
                segments.append(current)
                current = []
            continue
        current.append(point)
    if current:
        segments.append(current)
    return segments


def build_line_chart_svg(
    dataset_key: Tuple[str, ...],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    averages: AverageBuckets,
) -> str:
    values = [
        value
        for threshold in thresholds
        for value in averages.get(dataset_key, {}).get(threshold, {}).values()
        if value is not None and value >= 0
    ]
    title = dataset_title(dataset_key)
    if not values:
        return (
            '<div class="chart-card">'
            f'<h3>{html.escape(title)}</h3>'
            '<div class="chart-empty">No runtime data for this dataset.</div>'
            '</div>'
        )

    positive_values = [value for value in values if value > 0]
    if not positive_values:
        return (
            '<div class="chart-card">'
            f'<h3>{html.escape(title)}</h3>'
            '<div class="chart-empty">No positive runtime data for this dataset.</div>'
            '</div>'
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
            f'aria-label="Average runtime by t for {html.escape(title, quote=True)}">',
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
            value = averages.get(dataset_key, {}).get(threshold, {}).get(algorithm)
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
            value = averages.get(dataset_key, {}).get(threshold, {}).get(algorithm)
            if value is None:
                continue
            x = x_position(threshold_index)
            y = y_position(value)
            tooltip = f"{algorithm}; t={threshold}; average={value:.4f} ms"
            parts.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}" '
                'stroke="#111827" stroke-width=".6">'
                f'<title>{html.escape(tooltip)}</title></circle>'
            )

    y_title_x = 16
    y_title_y = top + plot_height / 2.0
    parts.append(
        f'<text x="{y_title_x}" y="{y_title_y:.2f}" text-anchor="middle" '
        f'transform="rotate(-90 {y_title_x} {y_title_y:.2f})" '
        'class="axis-title">Average runtime (ms, log scale)</text>'
    )
    parts.append(
        f'<text x="{left + plot_width / 2.0:.2f}" y="{height - 10}" '
        'text-anchor="middle" class="axis-title">t</text>'
    )

    parts.extend(["</svg>", "</div>", "</div>"])
    return "".join(parts)


def write_html_report(
    path: Path,
    dataset_keys: Sequence[Tuple[str, ...]],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    dataset_columns: Sequence[str],
    averages: AverageBuckets,
    counts: CountBuckets,
    runtime_metric: str,
    timeout_seconds: float,
    combine_label_count_suffixes: bool,
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
        handle.write("<title>SSM-GED Average Runtime by t</title>\n")
        handle.write(f"<style>{style}</style>\n</head><body>\n")
        handle.write("<h1>SSM-GED Average Runtime by t</h1>\n")
        suffix_note = (
            " Dataset suffixes _15, _30, _45, and _60 are combined by prefix."
            if combine_label_count_suffixes
            else ""
        )
        handle.write(
            f"<p>Metric: {html.escape(runtime_metric)}. Each cell is the "
            "arithmetic mean over all rows for the same dataset, t, and algorithm. "
            f"Timed-out runs and NA runtime cells are counted as "
            f"{html.escape(format_compact_number(timeout_seconds))} seconds."
            f"{html.escape(suffix_note)}</p>\n"
        )
        handle.write(
            '<div class="summary">'
            f"<span>datasets: {len(dataset_keys)}</span>"
            f"<span>algorithms: {html.escape(', '.join(algorithms))}</span>"
            f"<span>t values: {html.escape(', '.join(thresholds))}</span>"
            "</div>\n"
        )

        handle.write("<h2>Average runtime lines</h2>\n")
        handle.write('<div class="chart-grid">\n')
        for key in dataset_keys:
            handle.write(build_line_chart_svg(key, algorithms, thresholds, averages))
            handle.write("\n")
        handle.write("</div>\n")
        handle.write("</body></html>\n")


def main() -> int:
    args = parse_args()
    try:
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

        detected_algorithms = detect_algorithms(fieldnames, args.runtime_metric)
        if not detected_algorithms:
            raise ValueError(f"No algorithms detected for metric: {args.runtime_metric}")

        selected_algorithm_names = parse_algorithm_list(args.algorithms)
        if selected_algorithm_names is None:
            algorithms = detected_algorithms
        else:
            algorithms = resolve_algorithm_list(
                selected_algorithm_names,
                detected_algorithms,
            )

        runtime_columns = [f"{algorithm}_{args.runtime_metric}" for algorithm in algorithms]
        if not any(column in fieldnames for column in runtime_columns):
            raise ValueError(
                "No selected algorithm has the runtime metric column: "
                f"*_{args.runtime_metric}"
            )

        combine_label_count_suffixes = not args.exact_dataset_names
        dataset_keys, averages, counts = aggregate_average_runtimes(
            rows,
            algorithms,
            thresholds,
            dataset_columns,
            combine_label_count_suffixes,
            args.runtime_metric,
            timeout_ms,
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
            dataset_columns,
            averages,
            counts,
            args.runtime_metric,
            args.timeout_seconds,
            combine_label_count_suffixes,
        )

        print(f"Input files: {len(input_files)}")
        print(f"Rows: {len(rows)}")
        print(f"Datasets: {len(dataset_keys)}")
        print(f"Algorithms: {', '.join(algorithms)}")
        print(f"Runtime metric: {args.runtime_metric}")
        print(f"t values: {', '.join(thresholds)}")
        print(f"Chart timeout: {args.timeout_seconds:g} seconds")
        print(f"Output HTML: {output_path}")
        return 0
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
