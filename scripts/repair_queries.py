#!/usr/bin/env python3
"""Regenerate only SSM-GED queries that fail configured acceptance gates.

Replacement queries are connected induced subgraphs obtained by a reproducible
random walk with restarts on already selected vertices.  Existing qualifying
query files are copied byte-for-byte into staging and are never rewritten by
the apply phase.
"""

from __future__ import annotations

import argparse
from array import array
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import sys
from typing import Iterable, Mapping, Sequence

from graph_audit_common import (
    NO_EDGE_LABEL,
    QueryGraph,
    are_label_isomorphic,
    connected_components,
    isomorphism_bucket,
    parse_query_graph,
    sha256_file,
)


QUERY_GROUP_PATTERN = re.compile(r"^(?:(sparse|dense)_)?([1-9][0-9]*)$")


MANIFEST_COLUMNS = [
    "dataset_group",
    "dataset",
    "query_group",
    "query",
    "relative_query_path",
    "action",
    "repair_reason",
    "old_sha256",
    "new_sha256",
    "bytes_unchanged",
    "vertex_count",
    "old_edge_count",
    "new_edge_count",
    "old_slack",
    "new_slack",
    "query_seed",
    "attempts",
    "isomorphic_rejections",
    "source_set_rejections",
    "data_vertices",
    "induced_verified",
    "applied",
]


@dataclass
class DataGraphCSR:
    path: Path
    graph_id: str
    label_names: list[str]
    labels: array
    offsets: array
    neighbors: array

    @property
    def vertex_count(self) -> int:
        return len(self.labels)

    @property
    def edge_count(self) -> int:
        return len(self.neighbors) // 2

    def degree(self, vertex: int) -> int:
        return self.offsets[vertex + 1] - self.offsets[vertex]

    def random_neighbor(self, vertex: int, rng: random.Random) -> int | None:
        start = self.offsets[vertex]
        stop = self.offsets[vertex + 1]
        if start == stop:
            return None
        return self.neighbors[rng.randrange(start, stop)]

    def has_edge(self, left: int, right: int) -> bool:
        start = self.offsets[left]
        stop = self.offsets[left + 1]
        index = bisect_left(self.neighbors, right, start, stop)
        return index < stop and self.neighbors[index] == right

    def vertex_label(self, vertex: int) -> str:
        return self.label_names[self.labels[vertex]]

    @classmethod
    def load(cls, path: Path) -> "DataGraphCSR":
        if array("I").itemsize != 4 or array("Q").itemsize != 8:
            raise RuntimeError("repair tool requires 32-bit I and 64-bit Q arrays")

        graph_id = ""
        label_to_id: dict[str, int] = {}
        label_names: list[str] = []
        labels = array("I")
        degrees = array("I")
        edge_count = 0
        previous_edge: tuple[int, int] | None = None
        in_graph = False

        print(f"[index 1/2] {path}", file=sys.stderr, flush=True)
        with path.open("rb") as source:
            for line_number, raw_line in enumerate(source, 1):
                if not raw_line.strip():
                    continue
                record_type = raw_line[0]
                if record_type == ord("t"):
                    if in_graph:
                        break
                    fields = raw_line.split()
                    if len(fields) < 3:
                        raise ValueError(f"{path}:{line_number}: invalid graph header")
                    graph_id = fields[2].decode("utf-8")
                    in_graph = True
                elif not in_graph:
                    continue
                elif record_type == ord("v"):
                    fields = raw_line.split()
                    if len(fields) < 3:
                        raise ValueError(f"{path}:{line_number}: invalid vertex")
                    vertex = int(fields[1])
                    if vertex != len(labels):
                        raise ValueError(
                            f"{path}:{line_number}: expected vertex {len(labels)}, got {vertex}"
                        )
                    label = fields[2].decode("utf-8")
                    label_id = label_to_id.get(label)
                    if label_id is None:
                        label_id = len(label_names)
                        label_to_id[label] = label_id
                        label_names.append(label)
                    labels.append(label_id)
                    degrees.append(0)
                elif record_type == ord("e"):
                    fields = raw_line.split()
                    if len(fields) < 3:
                        raise ValueError(f"{path}:{line_number}: invalid edge")
                    left, right = int(fields[1]), int(fields[2])
                    if left == right:
                        continue
                    edge = (left, right) if left < right else (right, left)
                    if previous_edge is not None and edge < previous_edge:
                        raise ValueError(
                            f"{path}:{line_number}: edges are not sorted by normalized endpoints"
                        )
                    if edge == previous_edge:
                        continue
                    previous_edge = edge
                    if edge[1] >= len(labels):
                        raise ValueError(f"{path}:{line_number}: unknown edge endpoint")
                    degrees[edge[0]] += 1
                    degrees[edge[1]] += 1
                    edge_count += 1
                else:
                    raise ValueError(f"{path}:{line_number}: unknown record")

        offsets = array("Q", [0])
        running = 0
        for degree in degrees:
            running += degree
            offsets.append(running)
        neighbors = array("I", [0]) * running
        positions = array("Q", offsets)
        positions.pop()

        print(
            f"[index 2/2] vertices={len(labels)} edges={edge_count}",
            file=sys.stderr,
            flush=True,
        )
        previous_edge = None
        in_graph = False
        with path.open("rb") as source:
            for raw_line in source:
                if not raw_line.strip():
                    continue
                record_type = raw_line[0]
                if record_type == ord("t"):
                    if in_graph:
                        break
                    in_graph = True
                elif not in_graph or record_type != ord("e"):
                    continue
                else:
                    fields = raw_line.split()
                    left, right = int(fields[1]), int(fields[2])
                    if left == right:
                        continue
                    edge = (left, right) if left < right else (right, left)
                    if edge == previous_edge:
                        continue
                    previous_edge = edge
                    u, v = edge
                    neighbors[positions[u]] = v
                    positions[u] += 1
                    neighbors[positions[v]] = u
                    positions[v] += 1

        for vertex in range(len(labels)):
            if positions[vertex] != offsets[vertex + 1]:
                raise RuntimeError(f"CSR fill failed at vertex {vertex}")
        del positions
        del degrees

        return cls(path, graph_id, label_names, labels, offsets, neighbors)


@dataclass(frozen=True)
class Candidate:
    data_vertices: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]
    graph: QueryGraph
    content: bytes

    @property
    def slack(self) -> int:
        return len(self.edges) - len(self.data_vertices) + 1


class LabeledGraphDeduplicator:
    def __init__(self) -> None:
        self.buckets: dict[str, list[QueryGraph]] = defaultdict(list)

    def add(self, graph: QueryGraph) -> None:
        self.buckets[isomorphism_bucket(graph)].append(graph)

    def contains_isomorphic(self, graph: QueryGraph) -> bool:
        return any(
            are_label_isomorphic(graph, existing)
            for existing in self.buckets.get(isomorphism_bucket(graph), ())
        )


def _natural_query_key(path: Path) -> tuple[str, int, str]:
    try:
        number = int(path.stem)
    except ValueError:
        number = 2**31 - 1
    return str(path.parent), number, path.stem


def _query_group_design(query_group: str) -> tuple[str, int] | None:
    match = QUERY_GROUP_PATTERN.fullmatch(Path(query_group).name)
    if not match:
        return None
    return match.group(1) or "any", int(match.group(2))


def _query_group_violation(query_group: str, graph: QueryGraph) -> str | None:
    design = _query_group_design(query_group)
    if design is None:
        return None
    density_class, declared_vertex_count = design
    if graph.vertex_count != declared_vertex_count:
        return "query_group_size_mismatch"
    doubled_edges = 2 * graph.edge_count
    tripled_vertices = 3 * graph.vertex_count
    if density_class == "sparse" and doubled_edges > tripled_vertices:
        return "query_group_density_mismatch"
    if density_class == "dense" and doubled_edges <= tripled_vertices:
        return "query_group_density_mismatch"
    return None


def _repair_target(
    query_group: str,
    graph: QueryGraph,
    repair_reason: str,
    min_slack: int,
) -> tuple[int, int]:
    design = _query_group_design(query_group)
    vertex_count = design[1] if design is not None else graph.vertex_count
    old_slack = graph.edge_count - graph.vertex_count + 1
    if repair_reason == "isomorphic_duplicate":
        return vertex_count, old_slack
    if repair_reason == "slack_below_target":
        density_class = design[0] if design is not None else "any"
        if density_class == "dense":
            target_edges = (3 * vertex_count) // 2 + 1
            return vertex_count, max(min_slack, target_edges - vertex_count + 1)
        if density_class == "sparse":
            max_slack = (3 * vertex_count) // 2 - vertex_count + 1
            if min_slack > max_slack:
                raise RuntimeError(
                    f"{query_group}: sparse density and slack>={min_slack} are incompatible"
                )
        return vertex_count, min_slack
    if repair_reason in {"query_group_size_mismatch", "query_group_density_mismatch"}:
        density_class = design[0] if design is not None else "any"
        if density_class == "dense":
            # d > 3 means 2E > 3V, so floor(3V/2)+1 is the smallest edge count.
            target_edges = (3 * vertex_count) // 2 + 1
            return vertex_count, max(min_slack, target_edges - vertex_count + 1)
        if density_class == "sparse":
            # Move minimally to the inclusive d <= 3 boundary.
            target_edges = (3 * vertex_count) // 2
            target_slack = target_edges - vertex_count + 1
            if target_slack < min_slack:
                raise RuntimeError(
                    f"{query_group}: sparse density and slack>={min_slack} are incompatible"
                )
            return vertex_count, target_slack
        return vertex_count, min_slack
    raise RuntimeError(f"unknown repair reason: {repair_reason}")


def _discover_datasets(data_root: Path) -> list[tuple[str, str, Path]]:
    roots: list[tuple[str, Path]] = []
    if (data_root / "synthetic").is_dir():
        roots.append(("synthetic", data_root / "synthetic"))
    if (data_root / "real_graphs").is_dir():
        roots.append(("real_graphs", data_root / "real_graphs"))
    if not roots:
        roots.append(("custom", data_root))
    result = []
    for group, root in roots:
        for graph_path in sorted(root.rglob("graph_g.txt")):
            if "query_graph" in graph_path.parts:
                continue
            result.append((group, graph_path.parent.name, graph_path.parent))
    return result


def _stable_seed(global_seed: int, dataset: str, query_group: str, query: str) -> int:
    value = f"{global_seed}\0{dataset}\0{query_group}\0{query}".encode()
    return int.from_bytes(hashlib.sha256(value).digest()[:8], "little")


def _weighted_choice(
    rng: random.Random, values: Sequence[tuple[int, int]], weights: Sequence[int]
) -> tuple[int, int]:
    total = sum(weights)
    choice = rng.randrange(total)
    for value, weight in zip(values, weights):
        if choice < weight:
            return value
        choice -= weight
    return values[-1]


def _sample_vertices(
    data: DataGraphCSR,
    vertex_count: int,
    target_slack: int,
    rng: random.Random,
    proposal_count: int,
) -> tuple[int, ...] | None:
    start = None
    for _ in range(128):
        candidate = rng.randrange(data.vertex_count)
        if data.degree(candidate):
            start = candidate
            break
    if start is None:
        return None

    selected = [start]
    selected_set = {start}
    current = start
    edge_count = 0
    slack = 0

    while len(selected) < vertex_count:
        proposals: dict[int, int] = {}
        samples = proposal_count + 4 * len(selected)
        for _ in range(samples):
            base = current if rng.random() < 0.45 else selected[rng.randrange(len(selected))]
            neighbor = data.random_neighbor(base, rng)
            if neighbor is None or neighbor in selected_set:
                continue
            connections = sum(data.has_edge(neighbor, other) for other in selected)
            new_slack = slack + connections - 1
            if new_slack <= target_slack:
                proposals[neighbor] = connections
        if not proposals:
            return None

        values = list(proposals.items())
        needed = target_slack - slack
        cycle_values = [value for value in values if value[1] > 1]
        if needed and cycle_values and rng.random() < 0.88:
            values = cycle_values
            weights = [(connections - 1) ** 2 + 1 for _, connections in values]
        else:
            # Keeping tree-like extensions available preserves the sparse-query
            # character and prevents every cycle from accumulating immediately.
            weights = [2 if connections == 1 else 1 for _, connections in values]
        vertex, connections = _weighted_choice(rng, values, weights)
        selected.append(vertex)
        selected_set.add(vertex)
        current = vertex
        edge_count += connections
        slack = edge_count - len(selected) + 1

    return tuple(selected) if slack == target_slack else None


def _render_query(graph_id: str, labels: Sequence[str], edges: Sequence[tuple[int, int]]) -> bytes:
    lines = [f"t # {graph_id}"]
    lines.extend(f"v {vertex} {label}" for vertex, label in enumerate(labels))
    lines.extend(f"e {left} {right}" for left, right in edges)
    return ("\n".join(lines) + "\n").encode("utf-8")


def _build_candidate(
    data: DataGraphCSR,
    graph_id: str,
    path: Path,
    data_vertices: tuple[int, ...],
) -> Candidate:
    labels = tuple(data.vertex_label(vertex) for vertex in data_vertices)
    edges = tuple(
        (left, right)
        for left in range(len(data_vertices))
        for right in range(left + 1, len(data_vertices))
        if data.has_edge(data_vertices[left], data_vertices[right])
    )
    adjacency = [set() for _ in data_vertices]
    edge_labels = {}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
        edge_labels[(left, right)] = frozenset({NO_EDGE_LABEL})
    content = _render_query(graph_id, labels, edges)
    graph = QueryGraph(
        path=path,
        graph_id=graph_id,
        labels=labels,
        adjacency=tuple(frozenset(neighbors) for neighbors in adjacency),
        edge_labels=edge_labels,
        sha256=hashlib.sha256(content).hexdigest(),
        raw_edge_records=len(edges),
        duplicate_edge_records=0,
        self_loops=0,
        edge_records_sorted=True,
    )
    return Candidate(data_vertices, edges, graph, content)


def _write_tsv(path: Path, columns: list[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _copy_or_link_graph(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    destination.symlink_to(source.resolve())


def _apply_repairs(
    data_root: Path,
    staging_root: Path,
    backup_root: Path,
    rows: list[dict[str, object]],
) -> None:
    if backup_root.exists() and any(backup_root.iterdir()):
        raise RuntimeError(f"backup directory is not empty: {backup_root}")
    repair_rows = [row for row in rows if row["action"] == "regenerated"]

    for row in repair_rows:
        relative = Path(str(row["relative_query_path"]))
        target = data_root / relative
        staged = staging_root / relative
        if sha256_file(target) != row["old_sha256"]:
            raise RuntimeError(f"source changed before apply: {target}")
        if sha256_file(staged) != row["new_sha256"]:
            raise RuntimeError(f"staged hash mismatch: {staged}")
        backup = backup_root / relative
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, backup)

    prepared: list[tuple[Path, Path, str]] = []
    for row in repair_rows:
        relative = Path(str(row["relative_query_path"]))
        target = data_root / relative
        staged = staging_root / relative
        temporary = target.with_name(target.name + ".query_repair_tmp")
        shutil.copy2(staged, temporary)
        if sha256_file(temporary) != row["new_sha256"]:
            raise RuntimeError(f"temporary hash mismatch: {temporary}")
        prepared.append((temporary, target, str(row["new_sha256"])))

    for temporary, target, _ in prepared:
        os.replace(temporary, target)
    for _, target, expected_hash in prepared:
        if sha256_file(target) != expected_hash:
            raise RuntimeError(f"post-apply hash mismatch: {target}")


def _verify_staging(
    data_root: Path,
    staging_root: Path,
    min_slack: int,
) -> int:
    manifest_path = staging_root / "query_repair_manifest.tsv"
    if not manifest_path.is_file():
        raise RuntimeError(f"repair manifest not found: {manifest_path}")
    with manifest_path.open("r", encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if not rows:
        raise RuntimeError("repair manifest is empty")

    verification_rows: list[dict[str, object]] = []
    errors = 0
    rows_by_dataset: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_dataset[(row["dataset_group"], row["dataset"])].append(row)

    for (dataset_group, dataset), dataset_rows in sorted(rows_by_dataset.items()):
        data_path = data_root / dataset_group / dataset / "graph_g.txt"
        has_regenerated = any(row["action"] == "regenerated" for row in dataset_rows)
        data = DataGraphCSR.load(data_path) if has_regenerated else None
        regenerated_source_sets: set[frozenset[int]] = set()
        for index, row in enumerate(dataset_rows, 1):
            relative = Path(row["relative_query_path"])
            source_path = data_root / relative
            staged_path = staging_root / relative
            reasons: list[str] = []
            staged_graph: QueryGraph | None = None
            if not staged_path.is_file():
                reasons.append("staged_file_missing")
            else:
                staged_sha = sha256_file(staged_path)
                if staged_sha != row["new_sha256"]:
                    reasons.append("staged_sha_mismatch")
                staged_graph = parse_query_graph(staged_path)
                if staged_graph.vertex_count != int(row["vertex_count"]):
                    reasons.append("manifest_vertex_count_mismatch")
                staged_slack = staged_graph.edge_count - staged_graph.vertex_count + 1
                if staged_slack != int(row["new_slack"]):
                    reasons.append("manifest_slack_mismatch")
                if staged_slack < min_slack:
                    reasons.append("slack_below_target")
                if len(connected_components(staged_graph)) != 1:
                    reasons.append("disconnected")
                if _query_group_violation(row["query_group"], staged_graph) is not None:
                    reasons.append("query_group_design_violation")
            if row["action"] == "kept":
                if not source_path.is_file() or sha256_file(source_path) != row["old_sha256"]:
                    reasons.append("source_sha_mismatch")
                if staged_path.is_file() and source_path.is_file():
                    if staged_path.read_bytes() != source_path.read_bytes():
                        reasons.append("kept_bytes_changed")
            elif row["action"] == "regenerated" and staged_graph is not None:
                assert data is not None
                graph = staged_graph
                data_vertices = tuple(json.loads(row["data_vertices"]))
                if len(data_vertices) != graph.vertex_count:
                    reasons.append("mapping_size_mismatch")
                elif len(set(data_vertices)) != len(data_vertices):
                    reasons.append("mapping_not_injective")
                else:
                    source_set = frozenset(data_vertices)
                    if source_set in regenerated_source_sets:
                        reasons.append("duplicate_source_vertex_set")
                    regenerated_source_sets.add(source_set)
                    for query_vertex, data_vertex in enumerate(data_vertices):
                        if not 0 <= data_vertex < data.vertex_count:
                            reasons.append("mapping_vertex_out_of_range")
                            break
                        if graph.labels[query_vertex] != data.vertex_label(data_vertex):
                            reasons.append("mapped_label_mismatch")
                            break
                    for left in range(graph.vertex_count):
                        if reasons and reasons[-1] in {
                            "mapping_vertex_out_of_range",
                            "mapped_label_mismatch",
                        }:
                            break
                        for right in range(left + 1, graph.vertex_count):
                            query_has_edge = right in graph.adjacency[left]
                            data_has_edge = data.has_edge(
                                data_vertices[left], data_vertices[right]
                            )
                            if query_has_edge != data_has_edge:
                                reasons.append("not_an_induced_subgraph")
                                break
                        if reasons and reasons[-1] == "not_an_induced_subgraph":
                            break
            else:
                reasons.append("unknown_action")

            status = "OK" if not reasons else "FAIL"
            errors += bool(reasons)
            verification_rows.append(
                {
                    "dataset_group": dataset_group,
                    "dataset": dataset,
                    "query_group": row["query_group"],
                    "query": row["query"],
                    "action": row["action"],
                    "status": status,
                    "reasons": ",".join(reasons),
                }
            )
            if index == 1 or index % 200 == 0 or index == len(dataset_rows):
                print(
                    f"[verify] {dataset} {index}/{len(dataset_rows)} errors={errors}",
                    file=sys.stderr,
                    flush=True,
                )
        if data is not None:
            del data

    _write_tsv(
        staging_root / "query_repair_verification.tsv",
        ["dataset_group", "dataset", "query_group", "query", "action", "status", "reasons"],
        verification_rows,
    )
    print(
        f"verification complete: queries={len(verification_rows)} errors={errors}",
        file=sys.stderr,
    )
    return 1 if errors else 0


def _apply_verified_staging(
    data_root: Path,
    staging_root: Path,
    backup_root: Path,
) -> None:
    manifest_path = staging_root / "query_repair_manifest.tsv"
    verification_path = staging_root / "query_repair_verification.tsv"
    if not manifest_path.is_file() or not verification_path.is_file():
        raise RuntimeError("staging manifest or verification report is missing")
    with verification_path.open("r", encoding="utf-8", newline="") as source:
        verification_rows = list(csv.DictReader(source, delimiter="\t"))
    if not verification_rows or any(row["status"] != "OK" for row in verification_rows):
        raise RuntimeError("staging verification is missing or contains failures")
    with manifest_path.open("r", encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if len(rows) != len(verification_rows):
        raise RuntimeError("manifest and verification row counts differ")

    _apply_repairs(data_root, staging_root, backup_root, rows)
    for row in rows:
        if row["action"] == "regenerated":
            row["applied"] = 1
    _write_tsv(manifest_path, MANIFEST_COLUMNS, rows)
    for name in (
        "query_repair_manifest.tsv",
        "query_repair_summary.tsv",
        "query_repair_verification.tsv",
    ):
        shutil.copy2(staging_root / name, backup_root / name)
    print(f"[applied verified staging] backup={backup_root}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("test/datasets"))
    parser.add_argument(
        "--staging-dir", type=Path, default=Path("tmp/query_repair_staging")
    )
    parser.add_argument(
        "--backup-dir", type=Path, default=Path("tmp/query_repair_backup")
    )
    parser.add_argument("--datasets", default="", help="Comma-separated dataset names.")
    parser.add_argument("--min-slack", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--max-attempts", type=int, default=50000)
    parser.add_argument("--proposal-count", type=int, default=96)
    parser.add_argument(
        "--overwrite-staging",
        action="store_true",
        help="Remove an existing staging directory before generation.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="After successful staging validation, back up and replace only bad queries.",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Verify an existing staging manifest, mappings, and induced-subgraph identity.",
    )
    parser.add_argument(
        "--apply-staged",
        action="store_true",
        help="Apply an existing successful --verify-only staging area without regenerating.",
    )
    parser.add_argument(
        "--repair-isomorphic-duplicates",
        action="store_true",
        help="Keep the first exact labeled-isomorphism class member and regenerate the rest.",
    )
    parser.add_argument(
        "--repair-query-group-violations",
        action="store_true",
        help="Repair declared size/density violations (sparse d<=3, dense d>3).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_slack < 0 or args.max_attempts <= 0 or args.proposal_count <= 0:
        raise SystemExit("slack must be nonnegative; attempts and proposal count must be positive")

    data_root = args.data_root.resolve()
    staging_root = args.staging_dir.resolve()
    backup_root = args.backup_dir.resolve()
    selected_datasets = {item.strip() for item in args.datasets.split(",") if item.strip()}
    if args.verify_only and args.apply_staged:
        raise SystemExit("--verify-only and --apply-staged are mutually exclusive")
    if args.verify_only:
        try:
            return _verify_staging(data_root, staging_root, args.min_slack)
        except (OSError, RuntimeError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
    if args.apply_staged:
        try:
            _apply_verified_staging(data_root, staging_root, backup_root)
            return 0
        except (OSError, RuntimeError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1

    datasets = _discover_datasets(data_root)
    if selected_datasets:
        datasets = [entry for entry in datasets if entry[1] in selected_datasets]
    if not datasets:
        print("error: no datasets selected", file=sys.stderr)
        return 2

    if staging_root.exists():
        if not args.overwrite_staging:
            print(f"error: staging exists; pass --overwrite-staging: {staging_root}", file=sys.stderr)
            return 2
        shutil.rmtree(staging_root)
    staging_root.mkdir(parents=True)

    manifest_rows: list[dict[str, object]] = []
    summary_counter: Counter[tuple[str, str, str]] = Counter()
    try:
        for dataset_group, dataset, dataset_dir in datasets:
            query_root = dataset_dir / "query_graph"
            query_paths = sorted(query_root.rglob("*.txt"), key=_natural_query_key)
            parsed = [(path, parse_query_graph(path)) for path in query_paths]
            keep: list[tuple[Path, QueryGraph]] = []
            repair: list[tuple[Path, QueryGraph, str]] = []
            classification_deduplicator = LabeledGraphDeduplicator()
            for path, graph in parsed:
                slack = graph.edge_count - graph.vertex_count + 1
                query_group = str(path.parent.relative_to(query_root))
                group_violation = _query_group_violation(query_group, graph)
                if slack < args.min_slack:
                    repair.append((path, graph, "slack_below_target"))
                elif args.repair_query_group_violations and group_violation is not None:
                    repair.append((path, graph, group_violation))
                elif (
                    args.repair_isomorphic_duplicates
                    and classification_deduplicator.contains_isomorphic(graph)
                ):
                    repair.append((path, graph, "isomorphic_duplicate"))
                else:
                    keep.append((path, graph))
                    classification_deduplicator.add(graph)
            print(
                f"[dataset] {dataset}: total={len(parsed)} keep={len(keep)} regenerate={len(repair)}",
                file=sys.stderr,
                flush=True,
            )

            staged_dataset = staging_root / dataset_group / dataset
            _copy_or_link_graph(dataset_dir / "graph_g.txt", staged_dataset / "graph_g.txt")
            deduplicator = LabeledGraphDeduplicator()
            accepted_hashes: set[str] = set()
            for path, graph in keep:
                deduplicator.add(graph)
                accepted_hashes.add(graph.sha256)

            repair_paths = {path for path, _, _ in repair}
            for path, graph in parsed:
                relative = path.relative_to(data_root)
                staged_path = staging_root / relative
                staged_path.parent.mkdir(parents=True, exist_ok=True)
                if path not in repair_paths:
                    shutil.copy2(path, staged_path)
                    slack = graph.edge_count - graph.vertex_count + 1
                    manifest_rows.append(
                        {
                            "dataset_group": dataset_group,
                            "dataset": dataset,
                            "query_group": str(path.parent.relative_to(query_root)),
                            "query": path.stem,
                            "relative_query_path": str(relative),
                            "action": "kept",
                            "repair_reason": "",
                            "old_sha256": graph.sha256,
                            "new_sha256": graph.sha256,
                            "bytes_unchanged": 1,
                            "vertex_count": graph.vertex_count,
                            "old_edge_count": graph.edge_count,
                            "new_edge_count": graph.edge_count,
                            "old_slack": slack,
                            "new_slack": slack,
                            "query_seed": "",
                            "attempts": 0,
                            "isomorphic_rejections": 0,
                            "source_set_rejections": 0,
                            "data_vertices": "",
                            "induced_verified": 1,
                            "applied": 0,
                        }
                    )
                    summary_counter[(dataset, str(path.parent.relative_to(query_root)), "kept")] += 1

            if not repair:
                continue
            data = DataGraphCSR.load(dataset_dir / "graph_g.txt")
            accepted_source_sets: set[frozenset[int]] = set()
            for repair_index, (path, old_graph, repair_reason) in enumerate(repair, 1):
                query_group = str(path.parent.relative_to(query_root))
                query_seed = _stable_seed(args.seed, dataset, query_group, path.stem)
                rng = random.Random(query_seed)
                old_slack = old_graph.edge_count - old_graph.vertex_count + 1
                target_vertex_count, target_slack = _repair_target(
                    query_group,
                    old_graph,
                    repair_reason,
                    args.min_slack,
                )
                chosen: Candidate | None = None
                iso_rejections = 0
                source_rejections = 0
                attempts = 0
                for attempts in range(1, args.max_attempts + 1):
                    vertices = _sample_vertices(
                        data,
                        target_vertex_count,
                        target_slack,
                        rng,
                        args.proposal_count,
                    )
                    if vertices is None:
                        continue
                    source_set = frozenset(vertices)
                    if source_set in accepted_source_sets:
                        source_rejections += 1
                        continue
                    candidate = _build_candidate(data, old_graph.graph_id, path, vertices)
                    if candidate.slack != target_slack:
                        continue
                    if candidate.graph.sha256 in accepted_hashes:
                        iso_rejections += 1
                        continue
                    if deduplicator.contains_isomorphic(candidate.graph):
                        iso_rejections += 1
                        continue
                    chosen = candidate
                    break
                if chosen is None:
                    raise RuntimeError(
                        f"{dataset}/{query_group}/{path.name}: no unique candidate after "
                        f"{args.max_attempts} attempts"
                    )

                relative = path.relative_to(data_root)
                staged_path = staging_root / relative
                staged_path.parent.mkdir(parents=True, exist_ok=True)
                staged_path.write_bytes(chosen.content)
                parsed_candidate = parse_query_graph(staged_path)
                if parsed_candidate.sha256 != chosen.graph.sha256:
                    raise RuntimeError(f"staged parse/hash mismatch: {staged_path}")
                if len(connected_components(parsed_candidate)) != 1:
                    raise RuntimeError(f"generated query is disconnected: {staged_path}")
                if _query_group_violation(query_group, parsed_candidate) is not None:
                    raise RuntimeError(f"generated query violates group design: {staged_path}")

                accepted_source_sets.add(frozenset(chosen.data_vertices))
                accepted_hashes.add(chosen.graph.sha256)
                deduplicator.add(chosen.graph)
                manifest_rows.append(
                    {
                        "dataset_group": dataset_group,
                        "dataset": dataset,
                        "query_group": query_group,
                        "query": path.stem,
                        "relative_query_path": str(relative),
                        "action": "regenerated",
                        "repair_reason": repair_reason,
                        "old_sha256": old_graph.sha256,
                        "new_sha256": chosen.graph.sha256,
                        "bytes_unchanged": 0,
                        "vertex_count": target_vertex_count,
                        "old_edge_count": old_graph.edge_count,
                        "new_edge_count": len(chosen.edges),
                        "old_slack": old_slack,
                        "new_slack": chosen.slack,
                        "query_seed": query_seed,
                        "attempts": attempts,
                        "isomorphic_rejections": iso_rejections,
                        "source_set_rejections": source_rejections,
                        "data_vertices": json.dumps(chosen.data_vertices, separators=(",", ":")),
                        "induced_verified": 1,
                        "applied": int(args.apply),
                    }
                )
                summary_counter[(dataset, query_group, "regenerated")] += 1
                summary_counter[(dataset, query_group, repair_reason)] += 1
                summary_counter[(dataset, query_group, "attempts")] += attempts
                summary_counter[(dataset, query_group, "iso_rejections")] += iso_rejections
                if repair_index == 1 or repair_index % 25 == 0 or repair_index == len(repair):
                    print(
                        f"[repair] {dataset} {repair_index}/{len(repair)} "
                        f"{query_group}/{path.name} reason={repair_reason} "
                        f"target_vertices={target_vertex_count} target_slack={target_slack} "
                        f"attempts={attempts} "
                        f"iso_reject={iso_rejections}",
                        file=sys.stderr,
                        flush=True,
                    )

            # Explicitly release the large CSR before indexing the next graph.
            del data

        manifest_rows.sort(
            key=lambda row: (
                str(row["dataset_group"]),
                str(row["dataset"]),
                str(row["query_group"]),
                int(row["query"]) if str(row["query"]).isdigit() else 2**31 - 1,
                str(row["query"]),
            )
        )
        _write_tsv(
            staging_root / "query_repair_manifest.tsv", MANIFEST_COLUMNS, manifest_rows
        )

        summary_rows = []
        keys = sorted({(dataset, group) for dataset, group, _ in summary_counter})
        for dataset, group in keys:
            regenerated = summary_counter[(dataset, group, "regenerated")]
            attempts = summary_counter[(dataset, group, "attempts")]
            summary_rows.append(
                {
                    "dataset": dataset,
                    "query_group": group,
                    "kept": summary_counter[(dataset, group, "kept")],
                    "regenerated": regenerated,
                    "slack_repairs": summary_counter[
                        (dataset, group, "slack_below_target")
                    ],
                    "isomorphic_repairs": summary_counter[
                        (dataset, group, "isomorphic_duplicate")
                    ],
                    "query_group_size_repairs": summary_counter[
                        (dataset, group, "query_group_size_mismatch")
                    ],
                    "query_group_density_repairs": summary_counter[
                        (dataset, group, "query_group_density_mismatch")
                    ],
                    "total_attempts": attempts,
                    "mean_attempts_per_regenerated": (
                        f"{attempts / regenerated:.6g}" if regenerated else "0"
                    ),
                    "isomorphic_rejections": summary_counter[
                        (dataset, group, "iso_rejections")
                    ],
                }
            )
        _write_tsv(
            staging_root / "query_repair_summary.tsv",
            [
                "dataset",
                "query_group",
                "kept",
                "regenerated",
                "slack_repairs",
                "isomorphic_repairs",
                "query_group_size_repairs",
                "query_group_density_repairs",
                "total_attempts",
                "mean_attempts_per_regenerated",
                "isomorphic_rejections",
            ],
            summary_rows,
        )

        if args.apply:
            _apply_repairs(data_root, staging_root, backup_root, manifest_rows)
            shutil.copy2(
                staging_root / "query_repair_manifest.tsv",
                backup_root / "query_repair_manifest.tsv",
            )
            shutil.copy2(
                staging_root / "query_repair_summary.tsv",
                backup_root / "query_repair_summary.tsv",
            )
            print(f"[applied] backup={backup_root}", file=sys.stderr)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    regenerated_count = sum(row["action"] == "regenerated" for row in manifest_rows)
    print(
        f"complete: queries={len(manifest_rows)} regenerated={regenerated_count} "
        f"staging={staging_root}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
