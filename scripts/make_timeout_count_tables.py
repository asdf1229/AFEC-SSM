#!/usr/bin/env python3
"""Build an HTML table that counts runs without a captured result.

For each dataset, threshold, and algorithm, a run is counted when its selected
result field (``*_count`` by default) is missing or is not numeric.  This matches
the current ``compare.sh`` output, which writes ``NA`` for timed-out runs.
"""

from __future__ import annotations

import argparse
import html
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from make_average_runtime_tables import (
    dataset_key_for_row,
    dataset_sort_key,
    parse_columns,
)
from make_result_tables import (
    DEFAULT_CHART_THRESHOLDS,
    choose_output_dir,
    contains_timeout_marker,
    detect_algorithms,
    find_input_files,
    parse_algorithm_list,
    read_rows,
    resolve_algorithm_list,
    threshold_matches,
    to_float,
)


DEFAULT_OUTPUT_NAME = "timeout_counts_by_t.html"

DatasetKey = Tuple[str, ...]
CountBuckets = Dict[DatasetKey, Dict[str, Dict[str, int]]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count SSM-GED runs that did not capture a result and write an "
            "HTML table."
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
            "the selected result metric"
        ),
    )
    parser.add_argument(
        "--result-metric",
        default="count",
        help=(
            "metric whose missing/non-numeric value means no result was captured; "
            "default: count"
        ),
    )
    parser.add_argument(
        "--thresholds",
        default=",".join(DEFAULT_CHART_THRESHOLDS),
        help="comma-separated t values to include; default: 0,1,2,3,4,5,6",
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
    return parser.parse_args()


def parse_thresholds(value: str) -> List[str]:
    thresholds = [part.strip() for part in value.split(",") if part.strip()]
    if not thresholds:
        raise ValueError("--thresholds must contain at least one t value")
    return thresholds


def output_path_for_args(args: argparse.Namespace, input_files: Sequence[Path]) -> Path:
    output = Path(args.output)
    if output.is_absolute() or output.parent != Path("."):
        return output
    return choose_output_dir(args, input_files) / output


def result_was_not_captured(
    row: Dict[str, str],
    algorithm: str,
    result_metric: str,
) -> bool:
    """Return true when an algorithm has no usable result for this row."""
    algorithm_status = row.get(f"{algorithm}_status", "")
    if contains_timeout_marker(algorithm_status):
        return True

    raw_result = row.get(f"{algorithm}_{result_metric}", "")
    return to_float(raw_result) is None


def aggregate_timeout_counts(
    rows: Sequence[Dict[str, str]],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    dataset_columns: Sequence[str],
    combine_label_count_suffixes: bool,
    result_metric: str,
) -> Tuple[List[DatasetKey], CountBuckets, CountBuckets]:
    timeout_counts: CountBuckets = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    case_counts: CountBuckets = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    dataset_keys = set()

    for row in rows:
        matched_threshold = next(
            (
                threshold
                for threshold in thresholds
                if threshold_matches(row.get("threshold", ""), threshold)
            ),
            None,
        )
        if matched_threshold is None:
            continue

        key = dataset_key_for_row(
            row,
            dataset_columns,
            combine_label_count_suffixes,
        )
        dataset_keys.add(key)
        for algorithm in algorithms:
            case_counts[key][algorithm][matched_threshold] += 1
            if result_was_not_captured(row, algorithm, result_metric):
                timeout_counts[key][algorithm][matched_threshold] += 1

    return (
        sorted(dataset_keys, key=dataset_sort_key),
        timeout_counts,
        case_counts,
    )


def cell_html(timeout_count: int, case_count: int) -> str:
    if case_count == 0:
        return '<td class="empty" title="No cases">&mdash;</td>'

    css_class = "timeout" if timeout_count else "zero"
    title = f"{timeout_count} no-result cases out of {case_count} cases"
    return (
        f'<td class="number {css_class}" title="{html.escape(title)}">'
        f"{timeout_count}</td>"
    )


def totals_for_algorithm(
    dataset_keys: Sequence[DatasetKey],
    algorithm: str,
    threshold: str,
    counts: CountBuckets,
) -> int:
    return sum(counts[key][algorithm].get(threshold, 0) for key in dataset_keys)


def write_html_report(
    path: Path,
    dataset_keys: Sequence[DatasetKey],
    algorithms: Sequence[str],
    thresholds: Sequence[str],
    dataset_columns: Sequence[str],
    timeout_counts: CountBuckets,
    case_counts: CountBuckets,
    result_metric: str,
    combine_label_count_suffixes: bool,
) -> None:
    total_timeouts = sum(
        timeout_counts[key][algorithm].get(threshold, 0)
        for key in dataset_keys
        for algorithm in algorithms
        for threshold in thresholds
    )
    total_runs = sum(
        case_counts[key][algorithm].get(threshold, 0)
        for key in dataset_keys
        for algorithm in algorithms
        for threshold in thresholds
    )
    style = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 24px; color: #111827; }
h1 { font-size: 24px; margin: 0 0 8px; }
p { color: #4b5563; margin: 0 0 16px; }
.summary { display: flex; gap: 10px; flex-wrap: wrap; margin: 12px 0 18px; font-size: 12px; color: #374151; }
.summary span { background: #f9fafb; border: 1px solid #e5e7eb; padding: 4px 8px; }
.table-wrap { overflow-x: auto; border: 1px solid #d1d5db; border-radius: 6px; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
th, td { border-right: 1px solid #e5e7eb; border-bottom: 1px solid #e5e7eb; padding: 7px 9px; }
th:last-child, td:last-child { border-right: 0; }
thead th { background: #f3f4f6; text-align: center; white-space: nowrap; }
tbody th { background: #fafafa; text-align: left; font-weight: 600; white-space: nowrap; }
tbody tr.dataset-start > * { border-top: 2px solid #9ca3af; }
tbody tr:first-child > * { border-top: 0; }
td.number { text-align: right; font-variant-numeric: tabular-nums; }
td.timeout { background: #fee2e2; color: #991b1b; font-weight: 700; }
td.zero { color: #6b7280; }
td.empty { color: #9ca3af; text-align: center; }
tfoot th, tfoot td { background: #eef2ff; border-top: 2px solid #818cf8; font-weight: 700; }
"""

    with path.open("w", encoding="utf-8") as handle:
        handle.write('<!doctype html>\n<html><head><meta charset="utf-8">\n')
        handle.write("<title>SSM-GED Timeout Counts by t</title>\n")
        handle.write(f"<style>{style}</style>\n</head><body>\n")
        handle.write("<h1>SSM-GED Timeout Counts by t</h1>\n")
        suffix_note = (
            " Dataset suffixes _15, _30, _45, and _60 are combined by prefix."
            if combine_label_count_suffixes
            else ""
        )
        handle.write(
            f"<p>A timeout/no-result case is a run whose "
            f"<code>*_{html.escape(result_metric)}</code> value is missing or "
            f"non-numeric. Each cell is the number of such runs."
            f"{html.escape(suffix_note)}</p>\n"
        )
        handle.write(
            '<div class="summary">'
            f"<span>datasets: {len(dataset_keys)}</span>"
            f"<span>algorithms: {html.escape(', '.join(algorithms))}</span>"
            f"<span>selected runs: {total_runs}</span>"
            f"<span>timeout/no-result runs: {total_timeouts}</span>"
            "</div>\n"
        )

        handle.write('<div class="table-wrap"><table>\n<thead><tr>')
        for column in dataset_columns:
            handle.write(f'<th scope="col">{html.escape(column)}</th>')
        handle.write('<th scope="col">algorithm</th>')
        for threshold in thresholds:
            handle.write(f'<th scope="col">t={html.escape(threshold)}</th>')
        handle.write('<th scope="col">total</th></tr></thead>\n<tbody>\n')

        for key in dataset_keys:
            for algorithm_index, algorithm in enumerate(algorithms):
                row_class = ' class="dataset-start"' if algorithm_index == 0 else ""
                handle.write(f"<tr{row_class}>")
                if algorithm_index == 0:
                    for value in key:
                        handle.write(
                            f'<th scope="rowgroup" rowspan="{len(algorithms)}">'
                            f"{html.escape(value or '<empty>')}</th>"
                        )
                handle.write(f'<th scope="row">{html.escape(algorithm)}</th>')
                algorithm_timeouts = 0
                algorithm_cases = 0
                for threshold in thresholds:
                    timeout_count = timeout_counts[key][algorithm].get(threshold, 0)
                    case_count = case_counts[key][algorithm].get(threshold, 0)
                    algorithm_timeouts += timeout_count
                    algorithm_cases += case_count
                    handle.write(cell_html(timeout_count, case_count))
                handle.write(cell_html(algorithm_timeouts, algorithm_cases))
                handle.write("</tr>\n")

        handle.write("</tbody>\n<tfoot>\n")
        for algorithm_index, algorithm in enumerate(algorithms):
            handle.write("<tr>")
            if algorithm_index == 0:
                handle.write(
                    f'<th colspan="{len(dataset_columns)}" rowspan="{len(algorithms)}">'
                    "All datasets</th>"
                )
            handle.write(f'<th scope="row">{html.escape(algorithm)}</th>')
            algorithm_timeouts = 0
            algorithm_cases = 0
            for threshold in thresholds:
                timeout_count = totals_for_algorithm(
                    dataset_keys,
                    algorithm,
                    threshold,
                    timeout_counts,
                )
                case_count = totals_for_algorithm(
                    dataset_keys,
                    algorithm,
                    threshold,
                    case_counts,
                )
                algorithm_timeouts += timeout_count
                algorithm_cases += case_count
                handle.write(cell_html(timeout_count, case_count))
            handle.write(cell_html(algorithm_timeouts, algorithm_cases))
            handle.write("</tr>\n")
        handle.write("</tfoot>\n</table></div>\n</body></html>\n")


def main() -> int:
    args = parse_args()
    try:
        input_files = find_input_files(args.inputs)
        rows, fieldnames = read_rows(input_files)
        if not rows:
            raise ValueError("No result rows found.")

        thresholds = parse_thresholds(args.thresholds)
        dataset_columns = parse_columns(args.dataset_columns, "--dataset-columns")
        missing_dataset_columns = [
            column for column in dataset_columns if column not in fieldnames
        ]
        if missing_dataset_columns:
            missing = ", ".join(missing_dataset_columns)
            raise ValueError(f"Missing --dataset-columns column(s): {missing}")

        detected_algorithms = detect_algorithms(fieldnames, args.result_metric)
        if not detected_algorithms:
            raise ValueError(
                f"No algorithms detected for result metric: {args.result_metric}"
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
        dataset_keys, timeout_counts, case_counts = aggregate_timeout_counts(
            rows,
            algorithms,
            thresholds,
            dataset_columns,
            combine_label_count_suffixes,
            args.result_metric,
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
            timeout_counts,
            case_counts,
            args.result_metric,
            combine_label_count_suffixes,
        )

        total_timeouts = sum(
            timeout_counts[key][algorithm].get(threshold, 0)
            for key in dataset_keys
            for algorithm in algorithms
            for threshold in thresholds
        )
        print(f"Input files: {len(input_files)}")
        print(f"Rows: {len(rows)}")
        print(f"Datasets: {len(dataset_keys)}")
        print(f"Algorithms: {', '.join(algorithms)}")
        print(f"Result metric: {args.result_metric}")
        print(f"t values: {', '.join(thresholds)}")
        print(f"Timeout/no-result runs: {total_timeouts}")
        print(f"Output HTML: {output_path}")
        return 0
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
