#!/usr/bin/env python3
"""Repeated, count-checked SASUM efficiency comparison for one query."""

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
FALLBACK_RE = re.compile(r"Fallback Exact Runs:\s+(\d+)")


@dataclass(frozen=True)
class Observation:
    algorithm: str
    count: int
    load_ms: float
    run_ms: float
    total_ms: float


def run_once(
    executable: Path,
    data_graph: Path,
    query_graph: Path,
    threshold: int,
    timeout: float,
) -> Observation:
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
    if completed.returncode != 0:
        raise RuntimeError(
            f"{executable} exited with {completed.returncode}:\n{completed.stdout}"
        )
    if "(reached)" in completed.stdout:
        raise RuntimeError(f"{executable} reached a configured output/intermediate limit")
    fallback_match = FALLBACK_RE.search(completed.stdout)
    if fallback_match and int(fallback_match.group(1)) != 0:
        raise RuntimeError(
            f"{executable} used {fallback_match.group(1)} unexpected fallback exact run(s)"
        )
    match = SUMMARY_RE.search(completed.stdout)
    if not match:
        raise RuntimeError(f"missing SSM_GED_SUMMARY from {executable}:\n{completed.stdout}")
    return Observation(
        algorithm=match.group(1),
        count=int(match.group(2)),
        load_ms=float(match.group(3)),
        run_ms=float(match.group(4)),
        total_ms=float(match.group(5)),
    )


def summarize(values: list[float]) -> tuple[float, float, float, float]:
    return (
        statistics.median(values),
        statistics.mean(values),
        min(values),
        max(values),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run repeated cold-process SASUM comparisons and verify result counts."
    )
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--thresholds", type=int, nargs="+", default=[0, 1, 2, 3])
    parser.add_argument(
        "--algorithms",
        default="sasum,treespan",
        help="comma-separated algorithm keys",
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.repetitions <= 0 or args.warmups < 0:
        parser.error("--repetitions must be positive and --warmups non-negative")
    if any(threshold < 0 for threshold in args.thresholds):
        parser.error("--thresholds must be non-negative")
    if not args.data.is_file() or not args.query.is_file():
        parser.error("--data and --query must name existing files")

    algorithm_keys = [key.strip() for key in args.algorithms.split(",") if key.strip()]
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
            # Rotate execution order to reduce systematic first/last-run bias.
            offset = repetition % len(algorithm_keys)
            ordered_keys = algorithm_keys[offset:] + algorithm_keys[:offset]
            for key in ordered_keys:
                observations[key].append(
                    run_once(executables[key], args.data, args.query, threshold, args.timeout)
                )

        counts = {obs.count for values in observations.values() for obs in values}
        if len(counts) != 1:
            detail = {key: sorted({obs.count for obs in values}) for key, values in observations.items()}
            raise RuntimeError(f"result-count mismatch at threshold {threshold}: {detail}")
        expected_count = counts.pop()

        for key in algorithm_keys:
            values = observations[key]
            run_median, run_mean, run_min, run_max = summarize(
                [value.run_ms for value in values]
            )
            total_median, total_mean, _, _ = summarize(
                [value.total_ms for value in values]
            )
            rows.append(
                {
                    "threshold": threshold,
                    "algorithm": key,
                    "count": expected_count,
                    "repetitions": args.repetitions,
                    "run_ms_median": f"{run_median:.4f}",
                    "run_ms_mean": f"{run_mean:.4f}",
                    "run_ms_min": f"{run_min:.4f}",
                    "run_ms_max": f"{run_max:.4f}",
                    "total_ms_median": f"{total_median:.4f}",
                    "total_ms_mean": f"{total_mean:.4f}",
                }
            )

    fieldnames = list(rows[0]) if rows else []
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
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
    except (RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
