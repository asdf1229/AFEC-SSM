#!/usr/bin/env python3
"""Check SSM-GED executables against tiny brute-force SASUM fixtures."""

from __future__ import annotations

import argparse
import csv
import itertools
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


SUMMARY_RE = re.compile(r"SSM_GED_SUMMARY\s+.*?\bcount=(\d+)\b")
DIGEST_RE = re.compile(r"Result Digest:\s+([0-9a-f]{16}:[0-9a-f]{16})")
UINT64_MASK = (1 << 64) - 1
STAT_PATTERNS = {
    "query_subgraphs": re.compile(r"Query Subgraphs:\s+(\d+)"),
    "terminal_graphs": re.compile(r"Terminal Graphs:\s+(\d+)"),
    "base_graphs": re.compile(r"Base Graphs:\s+(\d+)"),
    "seed_graphs": re.compile(r"Seed Graphs:\s+(\d+)"),
}


@dataclass(frozen=True)
class Graph:
    labels: tuple[str, ...]
    edges: frozenset[tuple[int, int]]


def load_graph(path: Path) -> Graph:
    labels: dict[int, str] = {}
    edges: set[tuple[int, int]] = set()
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.split()
            if not fields:
                continue
            if fields[0] == "v":
                labels[int(fields[1])] = fields[2]
            elif fields[0] == "e":
                u, v = int(fields[1]), int(fields[2])
                edges.add((min(u, v), max(u, v)))
    expected_ids = list(range(len(labels)))
    if sorted(labels) != expected_ids:
        raise ValueError(f"vertex IDs in {path} must be contiguous from zero")
    return Graph(tuple(labels[i] for i in expected_ids), frozenset(edges))


def connected(vertex_count: int, edges: list[tuple[int, int]]) -> bool:
    if vertex_count == 0:
        return False
    adjacency = [[] for _ in range(vertex_count)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == vertex_count


def brute_force_mappings(
    query: Graph, data: Graph, threshold: int
) -> set[tuple[int, ...]]:
    candidates = [
        tuple(v for v, label in enumerate(data.labels) if label == query_label)
        for query_label in query.labels
    ]
    results: set[tuple[int, ...]] = set()
    for mapping in itertools.product(*candidates):
        if len(set(mapping)) != len(mapping):
            continue
        preserved: list[tuple[int, int]] = []
        missing = 0
        for u, v in query.edges:
            mapped_edge = (min(mapping[u], mapping[v]), max(mapping[u], mapping[v]))
            if mapped_edge in data.edges:
                preserved.append((u, v))
            else:
                missing += 1
        if missing <= threshold and connected(len(query.labels), preserved):
            results.add(mapping)
    return results


def mapping_hash(mapping: tuple[int, ...]) -> int:
    value = 14695981039346656037
    for vertex in mapping:
        for byte in int(vertex).to_bytes(4, byteorder="little", signed=False):
            value ^= byte
            value = (value * 1099511628211) & UINT64_MASK
    return value


def mapping_digest(mappings: set[tuple[int, ...]]) -> str:
    digest_sum = 0
    digest_xor = 0
    for mapping in mappings:
        value = mapping_hash(mapping)
        digest_sum = (digest_sum + value) & UINT64_MASK
        digest_xor ^= value
    return f"{digest_sum:016x}:{digest_xor:016x}"


def executable_result(
    executable: Path, data: Path, query: Path, threshold: int, timeout: float
) -> tuple[int, str | None, dict[str, int]]:
    try:
        completed = subprocess.run(
            [str(executable), "-d", str(data), "-q", str(query), "-t", str(threshold)],
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
    match = SUMMARY_RE.search(completed.stdout)
    if not match:
        raise RuntimeError(f"missing SSM_GED_SUMMARY from {executable}:\n{completed.stdout}")
    digest_match = DIGEST_RE.search(completed.stdout)
    digest = digest_match.group(1) if digest_match else None
    stats = {
        name: int(stat_match.group(1))
        for name, pattern in STAT_PATTERNS.items()
        if (stat_match := pattern.search(completed.stdout)) is not None
    }
    return int(match.group(1)), digest, stats


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executables", nargs="+", type=Path)
    parser.add_argument("--timeout", type=float, default=10.0)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")

    root = Path(__file__).resolve().parent
    with (root / "expected.tsv").open(encoding="utf-8", newline="") as handle:
        cases = list(csv.DictReader(handle, delimiter="\t"))
    with (root / "expected_stats.tsv").open(encoding="utf-8", newline="") as handle:
        expected_stats = {
            (row["dataset"], row["query"], int(row["threshold"])): {
                name: int(row[name]) for name in STAT_PATTERNS
            }
            for row in csv.DictReader(handle, delimiter="\t")
        }

    failures: list[str] = []
    for case in cases:
        dataset = case["dataset"]
        query_name = case["query"]
        threshold = int(case["threshold"])
        expected = int(case["expected_count"])
        data_path = root / dataset / "graph_g.txt"
        query_path = root / dataset / "query_graph" / f"{query_name}.txt"

        oracle_mappings = brute_force_mappings(
            load_graph(query_path), load_graph(data_path), threshold
        )
        oracle = len(oracle_mappings)
        oracle_digest = mapping_digest(oracle_mappings)
        case_name = f"{dataset}/{query_name}@{threshold}"
        if oracle != expected:
            failures.append(f"{case_name}: fixture expected {expected}, brute force found {oracle}")

        row = [case_name, f"oracle={oracle}"]
        for executable in args.executables:
            if not executable.is_file():
                failures.append(f"missing executable: {executable}")
                row.append(f"{executable.name}=MISSING")
                continue
            try:
                observed, observed_digest, observed_stats = executable_result(
                    executable.resolve(), data_path, query_path, threshold, args.timeout
                )
            except RuntimeError as error:
                failures.append(f"{case_name}: {error}")
                row.append(f"{executable.name}=ERROR")
                continue
            row.append(f"{executable.name}={observed}")
            if observed != oracle:
                failures.append(
                    f"{case_name}: {executable} returned {observed}, expected {oracle}"
                )
            is_sasum = executable.stem == "ssm_ged_sasum"
            if is_sasum and observed_digest is None:
                failures.append(f"{case_name}: {executable} did not report Result Digest")
            elif observed_digest is not None and observed_digest != oracle_digest:
                failures.append(
                    f"{case_name}: {executable} digest {observed_digest}, "
                    f"expected {oracle_digest}"
                )
            structural_oracle = expected_stats.get((dataset, query_name, threshold))
            if structural_oracle is not None and is_sasum:
                for name, expected_value in structural_oracle.items():
                    if observed_stats.get(name) != expected_value:
                        failures.append(
                            f"{case_name}: {executable} {name}="
                            f"{observed_stats.get(name)}, expected {expected_value}"
                        )
        print("\t".join(row))

    if failures:
        print("\nFAIL", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(f"\nPASS: {len(cases)} cases across {len(args.executables)} executable(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
