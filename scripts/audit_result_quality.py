#!/usr/bin/env python3
"""Summarize zero-result and output-capped query rates from benchmark results."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import re
import sys
from typing import Iterable, Mapping


COUNT_COLUMN = re.compile(r"^(.+)_count$")


def _write_tsv(path: Path, columns: list[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_nonnegative_integer(value: str) -> int | None:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 0 else None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="compare.sh summary.tsv")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: beside summary.tsv).",
    )
    parser.add_argument(
        "--output-limit",
        type=int,
        default=1_000_000_000,
        help="Count treated as output-capped (default: 1e9).",
    )
    parser.add_argument(
        "--zero-warning-ratio",
        type=float,
        default=0.80,
        help="Warn above this zero-result ratio (default: 0.80).",
    )
    parser.add_argument(
        "--cap-warning-ratio",
        type=float,
        default=0.20,
        help="Warn above this capped-result ratio (default: 0.20).",
    )
    parser.add_argument("--strict", action="store_true", help="Exit nonzero when warnings fire.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output_limit <= 0:
        raise SystemExit("--output-limit must be positive")
    if not 0.0 <= args.zero_warning_ratio <= 1.0:
        raise SystemExit("--zero-warning-ratio must be in [0, 1]")
    if not 0.0 <= args.cap_warning_ratio <= 1.0:
        raise SystemExit("--cap-warning-ratio must be in [0, 1]")

    summary_path = args.summary.resolve()
    try:
        with summary_path.open("r", encoding="utf-8", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError("summary has no header")
            algorithms = sorted(
                match.group(1)
                for column in reader.fieldnames
                if (match := COUNT_COLUMN.match(column)) and match.group(1) != "expected"
            )
            if not algorithms and "expected_count" not in reader.fieldnames:
                raise ValueError("no final *_count or expected_count column found")
            input_rows = list(reader)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    # The reference series avoids multiplying queries by the number of
    # algorithms; per-algorithm rows remain available to expose missing runs.
    series = (["reference"] if "expected_count" in (reader.fieldnames or []) else []) + algorithms
    groups: dict[tuple[str, str, str, str, str], list[int | None]] = defaultdict(list)
    for row in input_rows:
        dataset_group = row.get("dataset_group", "")
        dataset = row.get("dataset", "")
        query_group = row.get("query_group", "")
        theta = row.get("threshold", row.get("theta", ""))
        for name in series:
            column = "expected_count" if name == "reference" else f"{name}_count"
            value = _parse_nonnegative_integer(row.get(column, ""))
            groups[(dataset_group, dataset, query_group, theta, name)].append(value)
            groups[(dataset_group, dataset, "__ALL__", theta, name)].append(value)

    output_rows: list[dict[str, object]] = []
    issue_rows: list[dict[str, object]] = []
    for (dataset_group, dataset, query_group, theta, series_name), values in sorted(groups.items()):
        available = [value for value in values if value is not None]
        missing = len(values) - len(available)
        zeros = sum(value == 0 for value in available)
        unique = sum(value == 1 for value in available)
        moderate = sum(2 <= value < 1_000_000 for value in available)
        large = sum(1_000_000 <= value < args.output_limit for value in available)
        capped = sum(value >= args.output_limit for value in available)
        denominator = len(available)
        zero_ratio = zeros / denominator if denominator else 0.0
        capped_ratio = capped / denominator if denominator else 0.0
        warnings: list[str] = []
        if denominator and zero_ratio > args.zero_warning_ratio:
            warnings.append("ZERO_RESULT_RATIO_HIGH")
        if denominator and capped_ratio > args.cap_warning_ratio:
            warnings.append("OUTPUT_CAP_RATIO_HIGH")
        if missing:
            warnings.append("MISSING_RESULT_COUNTS")
        for code in warnings:
            issue_rows.append(
                {
                    "severity": "WARNING",
                    "code": code,
                    "dataset_group": dataset_group,
                    "dataset": dataset,
                    "query_group": query_group,
                    "theta": theta,
                    "series": series_name,
                    "details": (
                        f"total={len(values)}, available={denominator}, missing={missing}, "
                        f"zero_ratio={zero_ratio:.6g}, capped_ratio={capped_ratio:.6g}"
                    ),
                }
            )
        output_rows.append(
            {
                "dataset_group": dataset_group,
                "dataset": dataset,
                "query_group": query_group,
                "theta": theta,
                "series": series_name,
                "total_queries": len(values),
                "available_counts": denominator,
                "missing_counts": missing,
                "zero_results": zeros,
                "zero_result_ratio": f"{zero_ratio:.12g}",
                "unique_results": unique,
                "moderate_results_2_to_999999": moderate,
                "large_results_below_cap": large,
                "output_capped": capped,
                "output_capped_ratio": f"{capped_ratio:.12g}",
                "warning_codes": ",".join(warnings),
            }
        )

    output_dir = (args.output_dir or summary_path.parent).resolve()
    _write_tsv(
        output_dir / "result_quality_summary.tsv",
        [
            "dataset_group",
            "dataset",
            "query_group",
            "theta",
            "series",
            "total_queries",
            "available_counts",
            "missing_counts",
            "zero_results",
            "zero_result_ratio",
            "unique_results",
            "moderate_results_2_to_999999",
            "large_results_below_cap",
            "output_capped",
            "output_capped_ratio",
            "warning_codes",
        ],
        output_rows,
    )
    _write_tsv(
        output_dir / "result_quality_issues.tsv",
        [
            "severity",
            "code",
            "dataset_group",
            "dataset",
            "query_group",
            "theta",
            "series",
            "details",
        ],
        issue_rows,
    )
    print(
        f"groups={len(output_rows)} warnings={len(issue_rows)} output={output_dir}",
        file=sys.stderr,
    )
    if args.strict and issue_rows:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
