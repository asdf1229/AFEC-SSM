#!/usr/bin/env python3
"""
Build filtered-candidate-count line charts from SSM-GED summary results.

The script reads the same summary.tsv inputs as make_average_runtime_tables.py.
It averages the total number of candidates left after filtering for each
dataset, threshold, and algorithm, then writes one HTML report by default:

  - filter_candidates_by_t.html

Current summary.tsv files usually do not contain candidate-count columns, so
the script falls back to parsing "Filter Candidates:" lines from each
algorithm's *_output file.
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
    choose_output_dir,
    find_input_files,
    format_axis_runtime,
    natural_key,
    parse_algorithm_list,
    parse_chart_thresholds,
    read_rows,
    resolve_algorithm_list,
    threshold_matches,
    to_float,
)


DEFAULT_OUTPUT_NAME = "filter_candidates_by_t.html"
FILTER_CANDIDATE_RE = re.compile(r"^\s*(?:-\s*)?Filter Candidates:\s*([0-9][0-9,]*)\b")
MISSING_OUTPUT_VALUES = {"", "NA", "N/A", "NULL", "NONE"}
DIRECT_CANDIDATE_SUFFIXES = (
    "filter_candidates",
    "filter_candidate_count",
    "filtered_candidates",
    "filtered_candidate_count",
)


@dataclass
class SourceCounts:
    summary: int = 0
    output: int = 0
    missing: int = 0


AverageBuckets = Dict[Tuple[str, ...], Dict[str, Dict[str, Optional[float]]]]
CountBuckets = Dict[Tuple[str, ...], Dict[str, Dict[str, int]]]
OutputCandidateCache = Dict[str, Optional[float]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create dataset-level line charts comparing total filtered "
            "candidate counts from SSM-GED summary.tsv files."
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
            "candidate-count columns or output columns"
        ),
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
        help="do not parse *_output files when candidate-count columns are missing",
    )
    parser.add_argument(
        "--linear-scale",
        action="store_true",
        help="use a linear y axis instead of the default logarithmic y axis",
    )
    return parser.parse_args()


def output_path_for_args(args: argparse.Namespace, input_files: Sequence[Path]) -> Path:
    output = Path(args.output)
    if output.is_absolute() or output.parent != Path("."):
        return output
    return choose_output_dir(args, input_files) / output


def detect_candidate_algorithms(fieldnames: Iterable[str]) -> List[str]:
    candidates = set()
    for field in fieldnames:
        if field.endswith("_output"):
            prefix = field[: -len("_output")]
            if prefix:
                candidates.add(prefix)
            continue
        for suffix in DIRECT_CANDIDATE_SUFFIXES:
            full_suffix = f"_{suffix}"
            if field.endswith(full_suffix):
                prefix = field[: -len(full_suffix)]
                if prefix:
                    candidates.add(prefix)
    return sorted(candidates, key=natural_key)


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


def parse_filter_candidates(path: Path) -> Optional[float]:
    value: Optional[float] = None
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                match = FILTER_CANDIDATE_RE.match(line)
                if not match:
                    continue
                value = float(match.group(1).replace(",", ""))
    except OSError:
        return None
    return value


def output_candidate_for_row(
    row: Dict[str, str],
    algorithm: str,
    cache: OutputCandidateCache,
) -> Optional[float]:
    raw_value = str(row.get(f"{algorithm}_output", "")).strip()
    if raw_value.upper() in MISSING_OUTPUT_VALUES:
        return None
    if raw_value in cache:
        return cache[raw_value]

    path = resolve_output_path(row, algorithm)
    if path is None:
        cache[raw_value] = None
        return None

    value = parse_filter_candidates(path)
    cache[raw_value] = value
    return value


def direct_candidate_value(row: Dict[str, str], algorithm: str) -> Optional[float]:
    for suffix in DIRECT_CANDIDATE_SUFFIXES:
        value = to_float(row.get(f"{algorithm}_{suffix}", ""))
        if value is not None and value >= 0:
            return value
    return None


def candidate_value_for_row(
    row: Dict[str, str],
    algorithm: str,
    use_output_fallback: bool,
    output_cache: OutputCandidateCache,
) -> Tuple[Optional[float], str]:
    direct_value = direct_candidate_value(row, algorithm)
    if direct_value is not None:
        return direct_value, "summary"

    if use_output_fallback:
        output_value = output_candidate_for_row(row, algorithm, output_cache)
        if output_value is not None and output_value >= 0:
            return output_value, "output"

    return None, "missing"


def matched_threshold_for_row(
    row: Dict[str, str],
    thresholds: Sequence[str],
) -> Optional[str]:
    for threshold in thresholds:
        if threshold_matches(row.get("threshold", ""), threshold):
            return threshold
    return None


def aggregate_candidate_counts(
    rows: Sequence[Dict[str, str]],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    dataset_columns: Sequence[str],
    combine_label_count_suffixes: bool,
    use_output_fallback: bool,
) -> Tuple[List[Tuple[str, ...]], AverageBuckets, CountBuckets, SourceCounts]:
    value_buckets: Dict[Tuple[Tuple[str, ...], str, str], List[float]] = defaultdict(list)
    dataset_keys = set()
    output_cache: OutputCandidateCache = {}
    source_counts = SourceCounts()

    for row in rows:
        threshold = matched_threshold_for_row(row, thresholds)
        if threshold is None:
            continue

        key = dataset_key_for_row(
            row,
            dataset_columns,
            combine_label_count_suffixes,
        )
        dataset_keys.add(key)

        for algorithm in algorithms:
            value, source = candidate_value_for_row(
                row,
                algorithm,
                use_output_fallback,
                output_cache,
            )
            if source == "summary":
                source_counts.summary += 1
            elif source == "output":
                source_counts.output += 1
            else:
                source_counts.missing += 1

            if value is not None:
                value_buckets[(key, threshold, algorithm)].append(value)

    ordered_keys = sorted(dataset_keys, key=dataset_sort_key)
    averages: AverageBuckets = defaultdict(lambda: defaultdict(dict))
    counts: CountBuckets = defaultdict(lambda: defaultdict(dict))

    for key in ordered_keys:
        for threshold in thresholds:
            for algorithm in algorithms:
                values = value_buckets.get((key, threshold, algorithm), [])
                if not values:
                    averages[key][threshold][algorithm] = None
                    counts[key][threshold][algorithm] = 0
                    continue
                averages[key][threshold][algorithm] = sum(values) / len(values)
                counts[key][threshold][algorithm] = len(values)

    return ordered_keys, averages, counts, source_counts


def format_count(value: float) -> str:
    if value >= 1000000:
        return f"{value / 1000000:g}M"
    if value >= 1000:
        return f"{value / 1000:g}k"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.2f}"


def log_scale_bounds(values: Sequence[float]) -> Tuple[float, float, List[float]]:
    positive_values = [max(value, 1.0) for value in values if value > 0]
    if not positive_values:
        return 1.0, 10.0, [1.0, 10.0]

    min_value = min(positive_values)
    max_value = max(positive_values)
    min_exp = min(0, math.floor(math.log10(min_value)))
    max_exp = max(1, math.ceil(math.log10(max_value)))
    if min_exp == max_exp:
        max_exp += 1
    while max_exp - min_exp > 8:
        min_exp += 1
    return 10.0 ** min_exp, 10.0 ** max_exp, [
        10.0 ** exponent for exponent in range(min_exp, max_exp + 1)
    ]


def linear_scale_bounds(values: Sequence[float]) -> Tuple[float, float, List[float]]:
    max_value = max(values) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    y_max = max_value * 1.08
    if y_max <= 1:
        y_max = 1.0
    step = y_max / 5.0
    return 0.0, y_max, [step * index for index in range(6)]


def build_line_chart_svg(
    dataset_key: Tuple[str, ...],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    averages: AverageBuckets,
    counts: CountBuckets,
    use_linear_scale: bool,
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
            f"<h3>{html.escape(title)}</h3>"
            '<div class="chart-empty">No filter-candidate data for this dataset.</div>'
            "</div>"
        )

    left = 68
    right = 22
    top = 22
    bottom = 46
    plot_width = 620
    plot_height = 230
    width = left + plot_width + right
    height = top + plot_height + bottom
    axis_bottom = top + plot_height

    if use_linear_scale:
        y_min, y_max, ticks = linear_scale_bounds(values)

        def y_position(value: float) -> float:
            portion = (min(max(value, y_min), y_max) - y_min) / (y_max - y_min)
            return top + plot_height * (1.0 - portion)

        axis_title = "Average filter candidates"
    else:
        y_min, y_max, ticks = log_scale_bounds(values)
        log_min = math.log10(y_min)
        log_max = math.log10(y_max)

        def y_position(value: float) -> float:
            clipped = min(max(value, y_min), y_max)
            if value == 0:
                clipped = y_min
            portion = (math.log10(clipped) - log_min) / (log_max - log_min)
            return top + plot_height * (1.0 - portion)

        axis_title = "Average filter candidates (log scale)"

    def x_position(index: int) -> float:
        if len(thresholds) == 1:
            return left + plot_width / 2.0
        return left + plot_width * index / (len(thresholds) - 1)

    parts = [
        '<div class="chart-card">',
        f"<h3>{html.escape(title)}</h3>",
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
            f'aria-label="Average filter candidates by t for '
            f'{html.escape(title, quote=True)}">',
        ]
    )

    for tick in ticks:
        y = y_position(tick)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" '
            'stroke="#d1d5db" stroke-width="1" />'
        )
        label = format_axis_runtime(tick) if not use_linear_scale else format_count(tick)
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" '
            f'class="axis-label">{html.escape(label)}</text>'
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
            count = counts.get(dataset_key, {}).get(threshold, {}).get(algorithm, 0)
            x = x_position(threshold_index)
            y = y_position(value)
            tooltip = (
                f"{algorithm}; t={threshold}; average={value:.4f}; n={count}"
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
        f'class="axis-title">{html.escape(axis_title)}</text>'
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
    averages: AverageBuckets,
    counts: CountBuckets,
    source_counts: SourceCounts,
    combine_label_count_suffixes: bool,
    use_output_fallback: bool,
    use_linear_scale: bool,
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
        handle.write("<title>SSM-GED Filter Candidates by t</title>\n")
        handle.write(f"<style>{style}</style>\n</head><body>\n")
        handle.write("<h1>SSM-GED Filter Candidates by t</h1>\n")
        suffix_note = (
            " Dataset suffixes _15, _30, _45, and _60 are combined by prefix."
            if combine_label_count_suffixes
            else ""
        )
        fallback_note = (
            " Missing candidate-count columns are parsed from *_output files when possible."
            if use_output_fallback
            else " Output-file fallback is disabled."
        )
        scale_note = "linear" if use_linear_scale else "logarithmic"
        handle.write(
            "<p>Each point is the arithmetic mean of total candidates left after "
            "filtering for the same dataset, t, and algorithm. "
            f"The y axis is {html.escape(scale_note)}."
            f"{html.escape(fallback_note)}{html.escape(suffix_note)}</p>\n"
        )
        handle.write(
            '<div class="summary">'
            f"<span>datasets: {len(dataset_keys)}</span>"
            f"<span>algorithms: {html.escape(', '.join(algorithms))}</span>"
            f"<span>t values: {html.escape(', '.join(thresholds))}</span>"
            f"<span>summary cells: {source_counts.summary}</span>"
            f"<span>output-file cells: {source_counts.output}</span>"
            f"<span>missing cells: {source_counts.missing}</span>"
            "</div>\n"
        )

        handle.write("<h2>Average Filter Candidates</h2>\n")
        handle.write('<div class="chart-grid">\n')
        for key in dataset_keys:
            handle.write(
                build_line_chart_svg(
                    key,
                    algorithms,
                    thresholds,
                    averages,
                    counts,
                    use_linear_scale,
                )
            )
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
        dataset_columns = parse_columns(args.dataset_columns, "--dataset-columns")
        missing_dataset_columns = [
            column for column in dataset_columns if column not in fieldnames
        ]
        if missing_dataset_columns:
            missing = ", ".join(missing_dataset_columns)
            raise ValueError(f"Missing --dataset-columns column(s): {missing}")

        detected_algorithms = detect_candidate_algorithms(fieldnames)
        if not detected_algorithms:
            raise ValueError(
                "No algorithms detected from candidate-count or *_output columns."
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
        dataset_keys, averages, counts, source_counts = aggregate_candidate_counts(
            rows,
            algorithms,
            thresholds,
            dataset_columns,
            combine_label_count_suffixes,
            use_output_fallback,
        )
        if not dataset_keys:
            raise ValueError(
                "No rows matched the selected t values: " + ", ".join(thresholds)
            )
        if source_counts.summary + source_counts.output == 0:
            raise ValueError(
                "No filter-candidate values were found in summary columns or output files."
            )

        output_path = output_path_for_args(args, input_files)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        write_html_report(
            output_path,
            dataset_keys,
            algorithms,
            thresholds,
            averages,
            counts,
            source_counts,
            combine_label_count_suffixes,
            use_output_fallback,
            bool(args.linear_scale),
        )

        print(f"Input files: {len(input_files)}")
        print(f"Rows: {len(rows)}")
        print(f"Datasets: {len(dataset_keys)}")
        print(f"Algorithms: {', '.join(algorithms)}")
        print(f"t values: {', '.join(thresholds)}")
        print(f"Scale: {'linear' if args.linear_scale else 'logarithmic'}")
        print(
            "Filter candidate cells: "
            f"summary={source_counts.summary}, "
            f"output={source_counts.output}, "
            f"missing={source_counts.missing}"
        )
        print(f"Output HTML: {output_path}")
        return 0
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
