#!/usr/bin/env python3
"""
Build an HTML comparison report from SSM-GED summary results.

The repository's benchmark summaries are TSV files shaped like:

  dataset_group dataset query_group query threshold status expected_count
  cde_edge_ie_count cde_edge_ie_run_ms ... treespan_count treespan_run_ms ...

This script detects algorithm columns, compares the selected algorithm against
the best non-selected algorithm for a metric, and writes one report.html by
default:

  - report.html with colored speedup cells

CSV output is optional via --write-csv.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


KNOWN_META_PREFIXES = {
    "dataset",
    "dataset_group",
    "query",
    "query_group",
    "threshold",
    "status",
    "expected",
    "source",
    "source_file",
}

KNOWN_ALGORITHM_SUFFIXES = (
    "_run_ms",
    "_total_ms",
    "_load_ms",
    "_recursion_calls",
    "_count",
    "_output",
)

DEFAULT_DISPLAY_METRICS = ("count", "run_ms", "recursion_calls")
RUNTIME_LOW_COLOR_THRESHOLD_MS = 200.0
MISSING_RUNTIME_RATIO = 10.0
DEFAULT_OURS_ALGORITHM = "cde_edge_ie"
SPEEDUP_BUCKET_EDGES = tuple(
    [i / 10.0 for i in range(1, 11)] + [float(i) for i in range(2, 11)]
)
TOO_SHORT_RUNTIME_REASON = "compared_runtimes_too_short"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a colored SSM-GED HTML report with selected-algorithm speedups."
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
        "--ours",
        default=None,
        help=(
            "primary algorithm name or prefix; default: first --algorithms item, "
            f"or {DEFAULT_OURS_ALGORITHM} when --algorithms is not set"
        ),
    )
    parser.add_argument(
        "--algorithms",
        default=None,
        help=(
            "comma-separated algorithms to include in the comparison; "
            "if --ours is omitted, the first item is the primary algorithm"
        ),
    )
    parser.add_argument(
        "--metric",
        default="run_ms",
        help="metric used for best-other comparison; default: run_ms",
    )
    parser.add_argument(
        "--group-by",
        default="dataset_group,dataset",
        help="comma-separated columns defining one comparison object; default: dataset_group,dataset",
    )
    parser.add_argument(
        "--higher-is-better",
        action="store_true",
        help="use this for metrics where larger is better; default assumes smaller is better",
    )
    parser.add_argument(
        "--cap-ratio",
        type=float,
        default=4.0,
        help="ratio where color reaches full intensity; default: 4",
    )
    parser.add_argument(
        "--write-csv",
        action="store_true",
        help="also write per-object CSV files, index.csv, and all_tables.csv",
    )
    return parser.parse_args()


def find_input_files(inputs: Sequence[str]) -> List[Path]:
    files: List[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_file():
            files.append(path)
            continue
        if not path.is_dir():
            raise FileNotFoundError(f"Input does not exist: {item}")

        direct_summary = path / "summary.tsv"
        if direct_summary.is_file():
            files.append(direct_summary)
            continue

        summaries = sorted(path.rglob("summary.tsv"))
        if not summaries:
            raise FileNotFoundError(f"No summary.tsv found under: {item}")
        files.extend(summaries)

    deduped: List[Path] = []
    seen = set()
    for file_path in files:
        resolved = file_path.resolve()
        if resolved not in seen:
            deduped.append(file_path)
            seen.add(resolved)
    return deduped


def sniff_delimiter(file_path: Path) -> str:
    with file_path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
    if "\t" in sample.splitlines()[0]:
        return "\t"
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except csv.Error:
        return ","


def read_rows(files: Sequence[Path]) -> Tuple[List[Dict[str, str]], List[str]]:
    rows: List[Dict[str, str]] = []
    field_order: List[str] = []
    seen_fields = set()

    for file_path in files:
        delimiter = sniff_delimiter(file_path)
        with file_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if reader.fieldnames is None:
                continue
            for field in reader.fieldnames:
                if field not in seen_fields:
                    field_order.append(field)
                    seen_fields.add(field)

            source = file_path.parent.name
            for row in reader:
                clean_row = {key: (value or "").strip() for key, value in row.items() if key}
                clean_row["source"] = source
                clean_row["source_file"] = str(file_path)
                rows.append(clean_row)

    for extra in ("source", "source_file"):
        if extra not in seen_fields:
            field_order.append(extra)
            seen_fields.add(extra)
    return rows, field_order


def detect_algorithms(fieldnames: Iterable[str], metric: str) -> List[str]:
    candidates = set()
    for field in fieldnames:
        for suffix in KNOWN_ALGORITHM_SUFFIXES:
            if not field.endswith(suffix):
                continue
            prefix = field[: -len(suffix)]
            if not prefix or prefix in KNOWN_META_PREFIXES:
                continue
            if prefix == "expected":
                continue
            candidates.add(prefix)

    metric_suffix = f"_{metric}"
    usable = [algo for algo in candidates if f"{algo}{metric_suffix}" in fieldnames]
    return sorted(usable, key=natural_key)


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def resolve_algorithm(name_or_prefix: str, algorithms: Sequence[str]) -> str:
    if name_or_prefix in algorithms:
        return name_or_prefix

    wanted = normalize_name(name_or_prefix)
    exact = [algo for algo in algorithms if normalize_name(algo) == wanted]
    if len(exact) == 1:
        return exact[0]

    prefix_matches = [algo for algo in algorithms if normalize_name(algo).startswith(wanted)]
    if len(prefix_matches) == 1:
        return prefix_matches[0]

    contains_matches = [algo for algo in algorithms if wanted in normalize_name(algo)]
    if len(contains_matches) == 1:
        return contains_matches[0]

    available = ", ".join(algorithms) if algorithms else "<none>"
    raise ValueError(f"Cannot resolve algorithm '{name_or_prefix}'. Available: {available}")


def parse_algorithm_list(value: Optional[str]) -> Optional[List[str]]:
    if value is None:
        return None
    algorithms = [part.strip() for part in value.split(",") if part.strip()]
    if not algorithms:
        raise ValueError("--algorithms must contain at least one algorithm")
    return algorithms


def resolve_algorithm_list(
    names_or_prefixes: Sequence[str],
    available_algorithms: Sequence[str],
) -> List[str]:
    resolved = []
    seen = set()
    for name_or_prefix in names_or_prefixes:
        algorithm = resolve_algorithm(name_or_prefix, available_algorithms)
        if algorithm in seen:
            continue
        resolved.append(algorithm)
        seen.add(algorithm)
    return resolved


def parse_group_by(value: str) -> List[str]:
    columns = [part.strip() for part in value.split(",") if part.strip()]
    if not columns:
        raise ValueError("--group-by must contain at least one column")
    return columns


def to_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if math.isnan(number):
        return None
    return number


def better_value(
    values: Sequence[Tuple[str, float]], higher_is_better: bool
) -> Optional[Tuple[str, float]]:
    if not values:
        return None
    key = (lambda item: item[1])
    return max(values, key=key) if higher_is_better else min(values, key=key)


def compute_speedup(
    ours: Optional[float],
    best_other: Optional[float],
    higher_is_better: bool,
    missing_ratio: Optional[float] = None,
) -> Optional[float]:
    if (ours is not None and ours < 0) or (best_other is not None and best_other < 0):
        return None
    if ours is None or best_other is None:
        if missing_ratio is None or (ours is None and best_other is None):
            return None
        if ours is not None:
            return missing_ratio
        if best_other is not None:
            return 1.0 / missing_ratio
        return None

    if higher_is_better:
        numerator = ours
        denominator = best_other
    else:
        numerator = best_other
        denominator = ours

    if denominator == 0:
        if numerator == 0:
            return 1.0
        return math.inf
    return numerator / denominator


def blend_channel(start: int, end: int, amount: float) -> int:
    return int(round(start + (end - start) * amount))


def is_runtime_metric(metric: str) -> bool:
    return metric.endswith("_ms")


def suppress_runtime_color(
    metric: str,
    ours_value: Optional[float],
    best_value: Optional[float],
) -> bool:
    return (
        is_runtime_metric(metric)
        and ours_value is not None
        and best_value is not None
        and 0 <= ours_value < RUNTIME_LOW_COLOR_THRESHOLD_MS
        and 0 <= best_value < RUNTIME_LOW_COLOR_THRESHOLD_MS
    )


def ratio_color(
    ratio: Optional[float],
    cap_ratio: float,
    suppress_color: bool = False,
) -> Tuple[str, str]:
    if ratio is None:
        return "", ""
    if suppress_color:
        return "", ""
    if ratio == 1:
        return "#ffffff", "tie"

    cap_ratio = max(cap_ratio, 1.000001)
    if math.isinf(ratio):
        amount = 1.0
        target = (46, 125, 50)
        direction = "green"
    elif ratio > 1:
        amount = min(1.0, math.log(ratio, cap_ratio))
        target = (46, 125, 50)
        direction = "green"
    elif ratio > 0:
        amount = min(1.0, math.log(1.0 / ratio, cap_ratio))
        target = (198, 40, 40)
        direction = "red"
    else:
        amount = 1.0
        target = (198, 40, 40)
        direction = "red"

    white = (255, 255, 255)
    rgb = tuple(blend_channel(white[i], target[i], amount) for i in range(3))
    return "#{:02x}{:02x}{:02x}".format(*rgb), direction


def format_number(value: Optional[float], digits: int = 4) -> str:
    if value is None:
        return ""
    if math.isinf(value):
        return "inf"
    return f"{value:.{digits}f}"


def format_compact_number(value: float) -> str:
    if math.isinf(value):
        return "inf"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:g}"


def natural_key(value: object) -> Tuple[object, ...]:
    parts = re.split(r"(\d+)", str(value))
    key: List[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return tuple(key)


def row_sort_key(row: Dict[str, str]) -> Tuple[object, ...]:
    parts = []
    for column in ("query_group", "query", "threshold", "status"):
        value = row.get(column, "")
        if column == "threshold":
            number = to_float(value)
            parts.append(number if number is not None else value)
        else:
            parts.append(natural_key(value))
    return tuple(parts)


def safe_filename(value: str, fallback: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    if not safe:
        safe = fallback
    return safe[:140]


def choose_output_dir(args: argparse.Namespace, input_files: Sequence[Path]) -> Path:
    if args.out_dir:
        return Path(args.out_dir)
    if len(input_files) == 1:
        return input_files[0].parent
    return Path("result_report")


def build_output_columns(
    fieldnames: Sequence[str],
    algorithms: Sequence[str],
    ours_algorithm: str,
    metric: str,
) -> List[str]:
    preferred_common = (
        "source",
        "dataset_group",
        "dataset",
        "query_group",
        "query",
        "threshold",
        "status",
        "expected_count",
    )
    columns = [column for column in preferred_common if column in fieldnames]

    ordered_algorithms = [ours_algorithm] + [algo for algo in algorithms if algo != ours_algorithm]
    metric_names = list(DEFAULT_DISPLAY_METRICS)
    if metric not in metric_names:
        metric_names.append(metric)

    for algo in ordered_algorithms:
        for metric_name in metric_names:
            column = f"{algo}_{metric_name}"
            if column in fieldnames and column not in columns:
                columns.append(column)

    best_metric_col = f"best_other_{metric}"
    ours_metric_col = f"{ours_algorithm}_{metric}"
    for column in (
        "best_other_algorithm",
        best_metric_col,
        "ours_algorithm",
        ours_metric_col,
        "speedup_vs_best_other",
        "speedup_direction",
        "speedup_color",
        "speedup_comparable",
        "speedup_skip_reason",
    ):
        if column not in columns:
            columns.append(column)
    return columns


def enrich_row(
    row: Dict[str, str],
    algorithms: Sequence[str],
    ours_algorithm: str,
    metric: str,
    higher_is_better: bool,
    cap_ratio: float,
) -> Dict[str, str]:
    output = dict(row)
    ours_value = to_float(row.get(f"{ours_algorithm}_{metric}"))
    other_values = []
    for algo in algorithms:
        if algo == ours_algorithm:
            continue
        value = to_float(row.get(f"{algo}_{metric}"))
        if value is not None:
            other_values.append((algo, value))

    best = better_value(other_values, higher_is_better)
    if best is None:
        best_algo = ""
        best_value = None
    else:
        best_algo, best_value = best

    missing_ratio = None
    if is_runtime_metric(metric) and any(algo != ours_algorithm for algo in algorithms):
        missing_ratio = MISSING_RUNTIME_RATIO

    speedup = compute_speedup(
        ours_value,
        best_value,
        higher_is_better,
        missing_ratio,
    )
    too_short_runtime = suppress_runtime_color(metric, ours_value, best_value)
    displayed_speedup = None if too_short_runtime else speedup
    if too_short_runtime:
        color = ""
        direction = TOO_SHORT_RUNTIME_REASON
    else:
        color, direction = ratio_color(
            displayed_speedup,
            cap_ratio,
        )

    output["best_other_algorithm"] = best_algo
    output[f"best_other_{metric}"] = format_number(best_value)
    output["ours_algorithm"] = ours_algorithm
    output[f"{ours_algorithm}_{metric}"] = format_number(ours_value)
    output["speedup_vs_best_other"] = format_number(displayed_speedup)
    output["speedup_direction"] = direction
    output["speedup_color"] = color
    output["speedup_comparable"] = (
        "1" if displayed_speedup is not None and displayed_speedup > 0 else ""
    )
    if too_short_runtime:
        output["speedup_skip_reason"] = TOO_SHORT_RUNTIME_REASON
    elif displayed_speedup is None:
        output["speedup_skip_reason"] = "missing_or_invalid_speedup"
    else:
        output["speedup_skip_reason"] = ""
    return output


def group_rows(
    rows: Sequence[Dict[str, str]], group_by: Sequence[str]
) -> Dict[Tuple[str, ...], List[Dict[str, str]]]:
    groups: Dict[Tuple[str, ...], List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = tuple(row.get(column, "") for column in group_by)
        groups[key].append(row)
    return groups


def write_csv(path: Path, rows: Sequence[Dict[str, str]], columns: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_group(rows: Sequence[Dict[str, str]]) -> Dict[str, str]:
    ratios = []
    faster = slower = equal = comparable = 0
    too_short_runtime = 0
    for row in rows:
        if row.get("speedup_skip_reason") == TOO_SHORT_RUNTIME_REASON:
            too_short_runtime += 1
        if row.get("speedup_comparable") != "1":
            continue

        ratio = to_float(row.get("speedup_vs_best_other"))
        if ratio is None or ratio <= 0:
            continue
        comparable += 1
        ratios.append(ratio)
        if ratio > 1:
            faster += 1
        elif ratio < 1:
            slower += 1
        else:
            equal += 1

    geomean = ""
    if ratios:
        geomean_value = math.exp(sum(math.log(ratio) for ratio in ratios) / len(ratios))
        geomean = format_number(geomean_value)

    return {
        "rows": str(len(rows)),
        "comparable_rows": str(comparable),
        "not_compared_rows": str(len(rows) - comparable),
        "too_short_runtime_rows": str(too_short_runtime),
        "ours_faster_rows": str(faster),
        "ours_slower_rows": str(slower),
        "equal_rows": str(equal),
        "geomean_speedup": geomean,
    }


def aggregate_speedups(ratios: Sequence[float]) -> Optional[float]:
    valid = [ratio for ratio in ratios if ratio is not None and ratio > 0]
    if not valid:
        return None
    if any(math.isinf(ratio) for ratio in valid):
        return math.inf
    return math.exp(sum(math.log(ratio) for ratio in valid) / len(valid))


def aggregate_metric_values(values: Sequence[float]) -> Optional[float]:
    valid = [value for value in values if value is not None and value >= 0]
    if not valid:
        return None
    if len(valid) == 1:
        return valid[0]
    if any(value == 0 for value in valid):
        return sum(valid) / len(valid)
    return math.exp(sum(math.log(value) for value in valid) / len(valid))


def format_threshold_cell(
    ratio: Optional[float],
    ours_value: Optional[float],
    best_value: Optional[float],
    best_algorithm: str,
    metric: str,
    ours_algorithm: str,
) -> str:
    if ratio is None:
        return ""

    unit = " ms" if metric.endswith("_ms") else ""
    lines = [format_number(ratio)]
    if ours_value is not None:
        lines.append(f"{ours_algorithm}: {format_number(ours_value)}{unit}")
    if best_value is not None:
        best_label = f"best({best_algorithm})" if best_algorithm else "best"
        lines.append(f"{best_label}: {format_number(best_value)}{unit}")
    return "\n".join(lines)


def build_threshold_speedup_table(
    rows: Sequence[Dict[str, str]],
    cap_ratio: float,
    metric: str,
    ours_algorithm: str,
) -> Tuple[List[str], List[Dict[str, str]]]:
    label_columns = []
    sources = {row.get("source", "") for row in rows if row.get("source", "")}
    if len(sources) > 1:
        label_columns.append("source")
    for column in ("query_group", "query"):
        if any(row.get(column, "") for row in rows):
            label_columns.append(column)

    if not label_columns:
        label_columns = ["query"]

    thresholds = sorted(
        {row.get("threshold", "") for row in rows if row.get("threshold", "")},
        key=natural_key,
    )
    threshold_columns = [f"t={threshold}" for threshold in thresholds]

    buckets: Dict[Tuple[Tuple[str, ...], str], Dict[str, object]] = defaultdict(
        lambda: {"ratios": [], "ours_values": [], "best_values": [], "best_algorithms": set()}
    )
    row_keys = set()
    for row in rows:
        key = tuple(row.get(column, "") for column in label_columns)
        threshold = row.get("threshold", "")
        ratio = to_float(row.get("speedup_vs_best_other"))
        ours_value = to_float(row.get(f"{ours_algorithm}_{metric}"))
        best_value = to_float(row.get(f"best_other_{metric}"))
        best_algorithm = row.get("best_other_algorithm", "")
        if not threshold:
            continue
        row_keys.add(key)
        bucket = buckets[(key, threshold)]
        if ratio is not None and row.get("speedup_comparable") == "1":
            bucket["ratios"].append(ratio)
        if ours_value is not None:
            bucket["ours_values"].append(ours_value)
        if best_value is not None:
            bucket["best_values"].append(best_value)
        if best_algorithm:
            bucket["best_algorithms"].add(best_algorithm)

    pivot_rows = []
    for key in sorted(row_keys, key=lambda item: tuple(natural_key(part) for part in item)):
        output = {column: value for column, value in zip(label_columns, key)}
        for threshold, threshold_column in zip(thresholds, threshold_columns):
            bucket = buckets.get((key, threshold), {})
            ratio = aggregate_speedups(bucket.get("ratios", []))
            ours_value = aggregate_metric_values(bucket.get("ours_values", []))
            best_value = aggregate_metric_values(bucket.get("best_values", []))
            best_algorithms = sorted(bucket.get("best_algorithms", []), key=natural_key)
            if len(best_algorithms) == 1:
                best_algorithm = best_algorithms[0]
            elif best_algorithms:
                best_algorithm = "mixed"
            else:
                best_algorithm = ""
            color, direction = ratio_color(
                ratio,
                cap_ratio,
                suppress_runtime_color(metric, ours_value, best_value),
            )
            output[threshold_column] = format_threshold_cell(
                ratio, ours_value, best_value, best_algorithm, metric, ours_algorithm
            )
            output[f"__color__{threshold_column}"] = color
            output[f"__direction__{threshold_column}"] = direction
        pivot_rows.append(output)

    return label_columns + threshold_columns, pivot_rows


def text_color_for_background(hex_color: str) -> str:
    if not hex_color or hex_color == "#ffffff":
        return "#111827"
    try:
        red = int(hex_color[1:3], 16)
        green = int(hex_color[3:5], 16)
        blue = int(hex_color[5:7], 16)
    except (ValueError, IndexError):
        return "#111827"
    luma = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#ffffff" if luma < 140 else "#111827"


def write_table_block(
    handle,
    columns: Sequence[str],
    rows: Sequence[Dict[str, str]],
    speedup_column: str = "speedup_vs_best_other",
) -> None:
    handle.write(
        "<div class=\"table-block\">"
        "<div class=\"scroll-sync\"><div class=\"scroll-sync-inner\"></div></div>"
        "<div class=\"table-wrap\"><table>\n<thead><tr>"
    )
    for column in columns:
        css_class = "text" if is_text_column(column) else ""
        handle.write(f"<th class=\"{css_class}\">{html.escape(column)}</th>")
    handle.write("</tr></thead>\n<tbody>\n")

    for row in rows:
        handle.write("<tr>")
        for column in columns:
            value = row.get(column, "")
            css_class = "text" if is_text_column(column) else ""
            background = row.get(f"__color__{column}", "")
            if column == speedup_column:
                background = row.get("speedup_color", background)

            style_attr = ""
            if background:
                text_color = text_color_for_background(background)
                style_attr = (
                    f" style=\"background:{html.escape(background)};"
                    f" color:{html.escape(text_color)}; font-weight:600;"
                    " white-space:pre-line; line-height:1.35; min-width:130px\""
                )
            handle.write(f"<td class=\"{css_class}\"{style_attr}>{html.escape(value)}</td>")
        handle.write("</tr>\n")
    handle.write("</tbody></table></div></div>\n")


def write_summary_chips(
    handle,
    summary: Dict[str, str],
    ours_algorithm: str,
) -> None:
    ours_label = html.escape(ours_algorithm)
    handle.write(
        "<div class=\"summary\">"
        f"<span>rows: {html.escape(summary.get('rows', ''))}</span>"
        f"<span>comparable: {html.escape(summary.get('comparable_rows', ''))}</span>"
        f"<span>not compared: {html.escape(summary.get('not_compared_rows', ''))}</span>"
        f"<span>too-short runtime: {html.escape(summary.get('too_short_runtime_rows', ''))}</span>"
        f"<span>{ours_label} faster: {html.escape(summary.get('ours_faster_rows', ''))}</span>"
        f"<span>{ours_label} slower: {html.escape(summary.get('ours_slower_rows', ''))}</span>"
        f"<span>tie: {html.escape(summary.get('equal_rows', ''))}</span>"
        f"<span>geomean speedup: {html.escape(summary.get('geomean_speedup', ''))}</span>"
        "</div>\n"
    )


def build_speedup_distribution_table(
    rows: Sequence[Dict[str, str]],
) -> Tuple[List[str], List[Dict[str, str]]]:
    labels = []
    lower = None
    for edge in SPEEDUP_BUCKET_EDGES:
        upper_label = format_compact_number(edge)
        if lower is None:
            labels.append(f"<= {upper_label}")
        else:
            labels.append(f"({format_compact_number(lower)}, {upper_label}]")
        lower = edge
    labels.append(f"> {format_compact_number(SPEEDUP_BUCKET_EDGES[-1])}")

    counts = [0 for _ in labels]
    total = 0
    for row in rows:
        if row.get("speedup_comparable") != "1":
            continue
        ratio = to_float(row.get("speedup_vs_best_other"))
        if ratio is None or ratio <= 0:
            continue

        total += 1
        bucket_index = len(SPEEDUP_BUCKET_EDGES)
        if not math.isinf(ratio):
            for index, edge in enumerate(SPEEDUP_BUCKET_EDGES):
                if ratio <= edge:
                    bucket_index = index
                    break
        counts[bucket_index] += 1

    table_rows = []
    for label, count in zip(labels, counts):
        share = f"{count / total * 100:.2f}%" if total else ""
        table_rows.append(
            {
                "speedup_range": label,
                "count": str(count),
                "share": share,
            }
        )

    table_rows.append(
        {
            "speedup_range": "total",
            "count": str(total),
            "share": "100.00%" if total else "",
        }
    )
    return ["speedup_range", "count", "share"], table_rows


def write_overall_section(
    handle,
    grouped: Sequence[Tuple[str, Dict[str, str], Sequence[Dict[str, str]]]],
    ours_algorithm: str,
) -> None:
    all_rows = [row for _, _, rows in grouped for row in rows]
    summary = summarize_group(all_rows)
    bucket_columns, bucket_rows = build_speedup_distribution_table(all_rows)

    handle.write("<h2>Overall</h2>\n")
    write_summary_chips(handle, summary, ours_algorithm)
    handle.write("<h3>Speedup distribution</h3>\n")
    write_table_block(handle, bucket_columns, bucket_rows, speedup_column="")


def write_html_report(
    path: Path,
    grouped: Sequence[Tuple[str, Dict[str, str], Sequence[Dict[str, str]]]],
    columns: Sequence[str],
    metric: str,
    ours_algorithm: str,
    higher_is_better: bool,
    cap_ratio: float,
) -> None:
    style = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 24px; color: #111827; }
h1 { font-size: 24px; margin: 0 0 8px; }
h2 { font-size: 18px; margin: 28px 0 8px; }
h3 { font-size: 14px; margin: 16px 0 8px; }
p { color: #4b5563; margin: 0 0 18px; }
table { border-collapse: collapse; width: 100%; margin-bottom: 24px; font-size: 13px; }
th, td { border: 1px solid #d1d5db; padding: 6px 8px; text-align: right; white-space: nowrap; }
th { background: #f3f4f6; position: sticky; top: 0; z-index: 1; }
td.text, th.text { text-align: left; }
.table-block { margin-bottom: 24px; }
.table-wrap { overflow-x: auto; }
.scroll-sync { overflow-x: auto; overflow-y: hidden; height: 16px; border: 1px solid #d1d5db; border-bottom: 0; background: #f9fafb; }
.scroll-sync-inner { height: 1px; }
.meta { font-size: 12px; color: #6b7280; }
.legend { display: flex; gap: 12px; flex-wrap: wrap; margin: 14px 0 22px; font-size: 12px; color: #374151; }
.legend span { border: 1px solid #d1d5db; padding: 4px 8px; }
.summary { display: flex; gap: 12px; flex-wrap: wrap; margin: 8px 0 10px; font-size: 12px; color: #374151; }
.summary span { background: #f9fafb; border: 1px solid #e5e7eb; padding: 4px 8px; }
"""
    with path.open("w", encoding="utf-8") as handle:
        ours_label = html.escape(ours_algorithm)
        handle.write("<!doctype html>\n<html><head><meta charset=\"utf-8\">\n")
        handle.write("<title>SSM-GED Result Tables</title>\n")
        handle.write(f"<style>{style}</style>\n</head><body>\n")
        handle.write("<h1>SSM-GED Result Tables</h1>\n")
        metric_direction = "higher-is-better" if higher_is_better else "lower-is-better"
        handle.write(
            f"<p>Metric: {html.escape(metric)} ({metric_direction}). "
            "Speedup values above 1 favor "
            f"{ours_label}; values below 1 favor the best other algorithm.</p>\n"
        )
        handle.write(
            f"<div class=\"legend\">"
            f"<span style=\"background:#eaf2ea\">green: {ours_label} faster</span>"
            f"<span style=\"background:#faeded\">red: {ours_label} slower</span>"
            f"<span style=\"background:#2e7d32;color:#fff\">deeper color: more extreme ratio</span>"
            f"<span>not compared: both compared runtimes &lt; "
            f"{html.escape(format_compact_number(RUNTIME_LOW_COLOR_THRESHOLD_MS))} ms</span>"
            f"<span>missing runtime: treated as a 10x gap</span>"
            "</div>\n"
        )

        hidden_columns = {
            "speedup_direction",
            "speedup_color",
            "speedup_comparable",
            "speedup_skip_reason",
        }
        display_columns = [column for column in columns if column not in hidden_columns]

        write_overall_section(handle, grouped, ours_algorithm)

        for title, summary, rows in grouped:
            handle.write(f"<h2>{html.escape(title)}</h2>\n")
            write_summary_chips(handle, summary, ours_algorithm)
            handle.write("<h3>Detailed results</h3>\n")
            write_table_block(handle, display_columns, rows)

        handle.write(
            "<script>\n"
            "(function () {\n"
            "  function initBlock(block) {\n"
            "    const top = block.querySelector('.scroll-sync');\n"
            "    const inner = block.querySelector('.scroll-sync-inner');\n"
            "    const wrap = block.querySelector('.table-wrap');\n"
            "    const table = block.querySelector('table');\n"
            "    if (!top || !inner || !wrap || !table) return;\n"
            "    const resize = () => { inner.style.width = table.scrollWidth + 'px'; };\n"
            "    let syncing = false;\n"
            "    top.addEventListener('scroll', () => {\n"
            "      if (syncing) return;\n"
            "      syncing = true;\n"
            "      wrap.scrollLeft = top.scrollLeft;\n"
            "      syncing = false;\n"
            "    });\n"
            "    wrap.addEventListener('scroll', () => {\n"
            "      if (syncing) return;\n"
            "      syncing = true;\n"
            "      top.scrollLeft = wrap.scrollLeft;\n"
            "      syncing = false;\n"
            "    });\n"
            "    resize();\n"
            "    window.addEventListener('resize', resize);\n"
            "  }\n"
            "  document.querySelectorAll('.table-block').forEach(initBlock);\n"
            "}());\n"
            "</script>\n"
            "</body></html>\n"
        )


def write_speedup_by_t_report(
    path: Path,
    grouped: Sequence[Tuple[str, Dict[str, str], Sequence[Dict[str, str]]]],
    metric: str,
    ours_algorithm: str,
    higher_is_better: bool,
    cap_ratio: float,
) -> None:
    style = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 24px; color: #111827; }
h1 { font-size: 24px; margin: 0 0 8px; }
h2 { font-size: 18px; margin: 28px 0 8px; }
h3 { font-size: 14px; margin: 16px 0 8px; }
p { color: #4b5563; margin: 0 0 18px; }
table { border-collapse: collapse; width: 100%; margin-bottom: 24px; font-size: 13px; }
th, td { border: 1px solid #d1d5db; padding: 6px 8px; text-align: right; white-space: nowrap; }
th { background: #f3f4f6; position: sticky; top: 0; z-index: 1; }
td.text, th.text { text-align: left; }
.table-block { margin-bottom: 24px; }
.table-wrap { overflow-x: auto; }
.scroll-sync { overflow-x: auto; overflow-y: hidden; height: 16px; border: 1px solid #d1d5db; border-bottom: 0; background: #f9fafb; }
.scroll-sync-inner { height: 1px; }
.legend { display: flex; gap: 12px; flex-wrap: wrap; margin: 14px 0 22px; font-size: 12px; color: #374151; }
.legend span { border: 1px solid #d1d5db; padding: 4px 8px; }
.summary { display: flex; gap: 12px; flex-wrap: wrap; margin: 8px 0 10px; font-size: 12px; color: #374151; }
.summary span { background: #f9fafb; border: 1px solid #e5e7eb; padding: 4px 8px; }
"""
    with path.open("w", encoding="utf-8") as handle:
        ours_label = html.escape(ours_algorithm)
        handle.write("<!doctype html>\n<html><head><meta charset=\"utf-8\">\n")
        handle.write("<title>SSM-GED Speedup by t</title>\n")
        handle.write(f"<style>{style}</style>\n</head><body>\n")
        handle.write("<h1>SSM-GED Speedup by t</h1>\n")
        metric_direction = "higher-is-better" if higher_is_better else "lower-is-better"
        handle.write(
            f"<p>Metric: {html.escape(metric)} ({metric_direction}). "
            "Each cell is the running-time speedup of "
            f"{ours_label} against the best other algorithm at that t.</p>\n"
        )
        handle.write(
            f"<div class=\"legend\">"
            f"<span style=\"background:#eaf2ea\">green: {ours_label} faster</span>"
            f"<span style=\"background:#faeded\">red: {ours_label} slower</span>"
            f"<span style=\"background:#2e7d32;color:#fff\">deeper color: more extreme ratio</span>"
            f"<span>not compared: both compared runtimes &lt; "
            f"{html.escape(format_compact_number(RUNTIME_LOW_COLOR_THRESHOLD_MS))} ms</span>"
            f"<span>missing runtime: treated as a 10x gap</span>"
            "</div>\n"
        )

        write_overall_section(handle, grouped, ours_algorithm)

        for title, summary, rows in grouped:
            handle.write(f"<h2>{html.escape(title)}</h2>\n")
            write_summary_chips(handle, summary, ours_algorithm)
            threshold_columns, threshold_rows = build_threshold_speedup_table(
                rows, cap_ratio, metric, ours_algorithm
            )
            write_table_block(handle, threshold_columns, threshold_rows)

        handle.write(
            "<script>\n"
            "(function () {\n"
            "  function initBlock(block) {\n"
            "    const top = block.querySelector('.scroll-sync');\n"
            "    const inner = block.querySelector('.scroll-sync-inner');\n"
            "    const wrap = block.querySelector('.table-wrap');\n"
            "    const table = block.querySelector('table');\n"
            "    if (!top || !inner || !wrap || !table) return;\n"
            "    const resize = () => { inner.style.width = table.scrollWidth + 'px'; };\n"
            "    let syncing = false;\n"
            "    top.addEventListener('scroll', () => {\n"
            "      if (syncing) return;\n"
            "      syncing = true;\n"
            "      wrap.scrollLeft = top.scrollLeft;\n"
            "      syncing = false;\n"
            "    });\n"
            "    wrap.addEventListener('scroll', () => {\n"
            "      if (syncing) return;\n"
            "      syncing = true;\n"
            "      top.scrollLeft = wrap.scrollLeft;\n"
            "      syncing = false;\n"
            "    });\n"
            "    resize();\n"
            "    window.addEventListener('resize', resize);\n"
            "  }\n"
            "  document.querySelectorAll('.table-block').forEach(initBlock);\n"
            "}());\n"
            "</script>\n"
            "</body></html>\n"
        )


def is_text_column(column: str) -> bool:
    if column in {
        "source",
        "dataset_group",
        "dataset",
        "query_group",
        "query",
        "status",
        "best_other_algorithm",
        "ours_algorithm",
        "speedup_direction",
        "speedup_color",
        "speedup_skip_reason",
        "speedup_range",
    }:
        return True
    return column.endswith("_output") or column.endswith("_file")


def main() -> int:
    args = parse_args()

    try:
        input_files = find_input_files(args.inputs)
        rows, fieldnames = read_rows(input_files)
        if not rows:
            raise ValueError("No result rows found.")

        group_by = parse_group_by(args.group_by)
        missing_group_columns = [column for column in group_by if column not in fieldnames]
        if missing_group_columns:
            missing = ", ".join(missing_group_columns)
            raise ValueError(f"Missing --group-by column(s): {missing}")

        detected_algorithms = detect_algorithms(fieldnames, args.metric)
        selected_algorithm_names = parse_algorithm_list(args.algorithms)
        if selected_algorithm_names is None:
            algorithms = detected_algorithms
            ours_algorithm = resolve_algorithm(
                args.ours or DEFAULT_OURS_ALGORITHM,
                algorithms,
            )
        else:
            algorithms = resolve_algorithm_list(
                selected_algorithm_names,
                detected_algorithms,
            )
            if args.ours:
                ours_algorithm = resolve_algorithm(args.ours, algorithms)
            else:
                ours_algorithm = algorithms[0]

        output_columns = build_output_columns(
            fieldnames, algorithms, ours_algorithm, args.metric
        )
        higher_is_better = bool(args.higher_is_better)

        enriched_rows = [
            enrich_row(
                row,
                algorithms,
                ours_algorithm,
                args.metric,
                higher_is_better,
                args.cap_ratio,
            )
            for row in rows
        ]

        out_dir = choose_output_dir(args, input_files)
        out_dir.mkdir(parents=True, exist_ok=True)

        grouped_rows = group_rows(enriched_rows, group_by)
        html_groups = []
        index_rows = []
        all_rows = []

        for index, key in enumerate(sorted(grouped_rows, key=lambda item: tuple(natural_key(x) for x in item)), 1):
            group = sorted(grouped_rows[key], key=row_sort_key)
            title = " / ".join(value or "<empty>" for value in key)

            summary = summarize_group(group)
            html_groups.append((title, summary, group))

            if args.write_csv:
                file_name = safe_filename(title, f"object_{index:03d}")
                csv_path = out_dir / f"{index:03d}_{file_name}.csv"
                write_csv(csv_path, group, output_columns)

                index_row = {column: value for column, value in zip(group_by, key)}
                index_row.update(summary)
                index_row["csv_file"] = csv_path.name
                index_rows.append(index_row)

                for row in group:
                    combined = dict(row)
                    combined["comparison_object"] = title
                    all_rows.append(combined)

        if args.write_csv:
            index_columns = list(group_by) + [
                "rows",
                "comparable_rows",
                "not_compared_rows",
                "too_short_runtime_rows",
                "ours_faster_rows",
                "ours_slower_rows",
                "equal_rows",
                "geomean_speedup",
                "csv_file",
            ]
            write_csv(out_dir / "index.csv", index_rows, index_columns)
            write_csv(out_dir / "all_tables.csv", all_rows, ["comparison_object"] + output_columns)

        report_path = out_dir / "report.html"
        speedup_by_t_path = out_dir / "speedup_by_t.html"
        write_html_report(
            report_path,
            html_groups,
            output_columns,
            args.metric,
            ours_algorithm,
            higher_is_better,
            args.cap_ratio,
        )
        write_speedup_by_t_report(
            speedup_by_t_path,
            html_groups,
            args.metric,
            ours_algorithm,
            higher_is_better,
            args.cap_ratio,
        )

        print(f"Input files: {len(input_files)}")
        print(f"Rows: {len(rows)}")
        print(f"Algorithms: {', '.join(algorithms)}")
        print(f"Our algorithm: {ours_algorithm}")
        print(f"Metric: {args.metric}")
        print(f"Output directory: {out_dir}")
        print(f"HTML report: {report_path}")
        print(f"Speedup by t report: {speedup_by_t_path}")
        if args.write_csv:
            print(f"Per-object CSV files: {len(grouped_rows)}")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
