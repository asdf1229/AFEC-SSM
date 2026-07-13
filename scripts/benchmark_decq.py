#!/usr/bin/env python3
"""Repeated, count-checked efficiency comparison including the DecQ baseline."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


SUMMARY_RE = re.compile(
    r"SSM_GED_SUMMARY\s+algorithm=(\S+)\s+count=(\d+)\s+"
    r"load_ms=([0-9.]+)\s+run_ms=([0-9.]+)\s+total_ms=([0-9.]+)"
)
INTEGER_STAT_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9 +/_-]*):\s+(\d+)\s*$", re.MULTILINE)
FLOAT_MS_STAT_RE = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9 +/()_-]*):\s+([0-9]+(?:\.[0-9]+)?)\s+ms(?:\s|$)",
    re.MULTILINE,
)
LIMIT_RE = re.compile(
    r"(?:Output|Intermediate(?: Match)?|Pattern) Limit:.*(?:\(reached\)|\breached\b)",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class Observation:
    algorithm: str
    count: int
    load_ms: float
    run_ms: float
    total_ms: float
    integer_stats: dict[str, int]
    float_ms_stats: dict[str, float]


def normalize_stat_name(label: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def run_once(
    executable: Path,
    data_graph: Path,
    query_graph: Path,
    threshold: int,
    timeout: float,
) -> Observation:
    try:
        completed = subprocess.run(
            [
                str(executable),
                "-d",
                str(data_graph),
                "-q",
                str(query_graph),
                "-t",
                str(threshold),
            ],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"{executable} timed out after {timeout:g}s") from error

    if completed.returncode != 0:
        raise RuntimeError(
            f"{executable} exited with {completed.returncode}:\n{completed.stdout}"
        )
    if LIMIT_RE.search(completed.stdout):
        raise RuntimeError(f"{executable} reached a configured resource limit")

    match = SUMMARY_RE.search(completed.stdout)
    if not match:
        raise RuntimeError(f"missing SSM_GED_SUMMARY from {executable}:\n{completed.stdout}")

    integer_stats = {
        normalize_stat_name(label): int(value)
        for label, value in INTEGER_STAT_RE.findall(completed.stdout)
    }
    # The first "Total Time" belongs to the algorithm's stats block. The
    # common runner prints another total later that also includes file parsing
    # and stats formatting, so preserve the first occurrence of every label.
    float_ms_stats: dict[str, float] = {}
    for label, value in FLOAT_MS_STAT_RE.findall(completed.stdout):
        float_ms_stats.setdefault(normalize_stat_name(label), float(value))
    return Observation(
        algorithm=match.group(1),
        count=int(match.group(2)),
        load_ms=float(match.group(3)),
        run_ms=float(match.group(4)),
        total_ms=float(match.group(5)),
        integer_stats=integer_stats,
        float_ms_stats=float_ms_stats,
    )


def summarize(values: list[float]) -> tuple[float, float, float, float]:
    return statistics.median(values), statistics.mean(values), min(values), max(values)


def stable_stat(observations: list[Observation], *names: str) -> str:
    values = {
        observation.integer_stats[name]
        for observation in observations
        for name in names
        if name in observation.integer_stats
    }
    if not values:
        return ""
    if len(values) != 1:
        raise RuntimeError(f"non-deterministic structural counter {names}: {sorted(values)}")
    return str(values.pop())


def median_float_stat(observations: list[Observation], name: str) -> str:
    values = [
        observation.float_ms_stats[name]
        for observation in observations
        if name in observation.float_ms_stats
    ]
    return f"{statistics.median(values):.4f}" if values else ""


def median_float_difference(
    observations: list[Observation], total_name: str, subtract_name: str
) -> str:
    values = [
        observation.float_ms_stats[total_name]
        - observation.float_ms_stats[subtract_name]
        for observation in observations
        if total_name in observation.float_ms_stats
        and subtract_name in observation.float_ms_stats
    ]
    return f"{statistics.median(values):.4f}" if values else ""


def stable_stat_sum(observations: list[Observation], *names: str) -> str:
    values = {
        sum(observation.integer_stats[name] for name in names)
        for observation in observations
        if all(name in observation.integer_stats for name in names)
    }
    if not values:
        return ""
    if len(values) != 1:
        raise RuntimeError(f"non-deterministic structural counter sum {names}: {sorted(values)}")
    return str(values.pop())


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run repeated cold-process DecQ comparisons, reject limited runs, "
            "and require identical unique-result counts."
        )
    )
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--thresholds", type=int, nargs="+", default=[0, 1, 2, 3])
    parser.add_argument(
        "--algorithms",
        default="decq,treespan",
        help="comma-separated executable suffixes",
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.repetitions <= 0 or args.warmups < 0:
        parser.error("--repetitions must be positive and --warmups non-negative")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if any(threshold < 0 for threshold in args.thresholds):
        parser.error("--thresholds must be non-negative")
    if not args.data.is_file() or not args.query.is_file():
        parser.error("--data and --query must name existing files")

    algorithm_keys = [key.strip() for key in args.algorithms.split(",") if key.strip()]
    if not algorithm_keys:
        parser.error("--algorithms must contain at least one key")

    executables: dict[str, Path] = {}
    for key in algorithm_keys:
        executable = (args.build_dir / f"ssm_ged_{key}").resolve()
        if not executable.is_file():
            parser.error(f"missing executable: {executable}")
        executables[key] = executable

    rows: list[dict[str, object]] = []
    for threshold in args.thresholds:
        for _ in range(args.warmups):
            for key in algorithm_keys:
                run_once(executables[key], args.data, args.query, threshold, args.timeout)

        observations: dict[str, list[Observation]] = {key: [] for key in algorithm_keys}
        for repetition in range(args.repetitions):
            offset = repetition % len(algorithm_keys)
            ordered_keys = algorithm_keys[offset:] + algorithm_keys[:offset]
            for key in ordered_keys:
                observations[key].append(
                    run_once(
                        executables[key], args.data, args.query, threshold, args.timeout
                    )
                )

        per_algorithm_counts = {
            key: sorted({observation.count for observation in values})
            for key, values in observations.items()
        }
        all_counts = {count for counts in per_algorithm_counts.values() for count in counts}
        if len(all_counts) != 1:
            raise RuntimeError(
                f"unique-result count mismatch at threshold {threshold}: "
                f"{per_algorithm_counts}"
            )
        expected_count = all_counts.pop()

        for key in algorithm_keys:
            values = observations[key]
            run_median, run_mean, run_min, run_max = summarize(
                [value.run_ms for value in values]
            )
            total_median, total_mean, _, _ = summarize(
                [value.total_ms for value in values]
            )
            algorithm_times = [
                value.float_ms_stats.get("total_time", value.run_ms) for value in values
            ]
            algorithm_median, algorithm_mean, algorithm_min, algorithm_max = summarize(
                algorithm_times
            )
            rows.append(
                {
                    "threshold": threshold,
                    "algorithm": key,
                    "unique_count": expected_count,
                    "repetitions": args.repetitions,
                    "algorithm_ms_median": f"{algorithm_median:.4f}",
                    "algorithm_ms_mean": f"{algorithm_mean:.4f}",
                    "algorithm_ms_min": f"{algorithm_min:.4f}",
                    "algorithm_ms_max": f"{algorithm_max:.4f}",
                    "runner_run_ms_median": f"{run_median:.4f}",
                    "runner_run_ms_mean": f"{run_mean:.4f}",
                    "runner_run_ms_min": f"{run_min:.4f}",
                    "runner_run_ms_max": f"{run_max:.4f}",
                    "process_total_ms_median": f"{total_median:.4f}",
                    "process_total_ms_mean": f"{total_mean:.4f}",
                    "decq_index_ms_median": median_float_stat(values, "data_label_nlf_index"),
                    "decq_decomposition_ms_median": median_float_stat(
                        values, "query_decomposition"
                    ),
                    "decq_local_pattern_ms_median": median_float_stat(
                        values, "local_pattern_generation"
                    ),
                    "decq_local_matching_ms_median": median_float_stat(
                        values, "local_exact_matching"
                    ),
                    "decq_lattice_ms_median": median_float_stat(
                        values, "global_lattice_generation"
                    ),
                    "decq_minimal_merge_ms_median": median_float_stat(
                        values, "minimal_sharing_merge"
                    ),
                    "decq_edge_validation_ms_median": median_float_stat(
                        values, "non_minimal_edge_validation"
                    ),
                    "decq_output_adapter_ms_median": median_float_stat(
                        values, "canonical_output_adapter"
                    ),
                    "decq_core_ms_median": median_float_difference(
                        values, "total_time", "canonical_output_adapter"
                    ),
                    "local_search_partial_rows": stable_stat(
                        values, "local_search_partial_rows"
                    ),
                    "materialized_match_rows": stable_stat(
                        values, "materialized_match_rows"
                    ),
                    "produced_intermediate_rows": stable_stat_sum(
                        values, "local_search_partial_rows", "materialized_match_rows"
                    ),
                    "global_pattern_rows": stable_stat(
                        values, "global_pattern_rows", "query_pattern_rows"
                    ),
                    "noncanonical_pattern_rows": stable_stat(
                        values, "non_canonical_pattern_rows"
                    ),
                    "duplicate_canonical_mappings": stable_stat(
                        values, "duplicate_canonical_mappings"
                    ),
                    "peak_live_match_rows": stable_stat(values, "peak_live_match_rows"),
                    "fragments": stable_stat(values, "fragments", "query_fragments"),
                    "local_patterns": stable_stat(values, "local_patterns"),
                    "global_patterns": stable_stat(values, "global_patterns", "query_subgraphs"),
                    "minimal_patterns": stable_stat(values, "minimal_patterns"),
                    "hash_joins": stable_stat(values, "hash_joins"),
                    "exact_match_executions": stable_stat(
                        values, "exact_match_executions", "exact_matching_executions"
                    ),
                }
            )

    fieldnames = list(rows[0]) if rows else []
    writer = csv.DictWriter(
        sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            output_writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            output_writer.writeheader()
            output_writer.writerows(rows)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
