#!/usr/bin/env python3
"""Shared graph parsing and statistics for the SSM-GED data audit tools.

The project matcher uses vertex labels and an undirected simple-graph view.  The
helpers below deliberately keep format diagnostics (self-loops, duplicate edge
records, and edge labels) separate from that semantic view.
"""

from __future__ import annotations

from array import array
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import struct
from typing import Iterable, Iterator, Mapping, Sequence


NO_EDGE_LABEL = "__NO_EDGE_LABEL__"


class GraphAuditError(ValueError):
    """Raised when an audit input cannot be interpreted as one graph."""


@dataclass(frozen=True)
class QueryGraph:
    path: Path
    graph_id: str
    labels: tuple[str, ...]
    adjacency: tuple[frozenset[int], ...]
    edge_labels: Mapping[tuple[int, int], frozenset[str]]
    sha256: str
    raw_edge_records: int
    duplicate_edge_records: int
    self_loops: int
    edge_records_sorted: bool

    @property
    def vertex_count(self) -> int:
        return len(self.labels)

    @property
    def edge_count(self) -> int:
        return sum(len(neighbors) for neighbors in self.adjacency) // 2


@dataclass(frozen=True)
class DataGraphStats:
    path: Path
    graph_id: str
    sha256: str
    topology_sha256: str
    vertex_count: int
    edge_count: int
    raw_edge_records: int
    duplicate_edge_records: int
    duplicate_detection: str
    self_loops: int
    edge_records_sorted: bool
    component_count: int
    largest_component_vertices: int
    isolated_vertices: int
    mean_degree: float
    max_degree: int
    degree_p50: int
    degree_p90: int
    degree_p95: int
    degree_p99: int
    density: float
    label_frequencies: Mapping[str, int]
    label_entropy_bits: float
    label_pair_edges: frozenset[tuple[str, str]]
    labeled_edge_combinations: frozenset[tuple[str, str, str]]


def json_compact(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _meaningful_lines(path: Path) -> Iterator[tuple[int, str]]:
    with path.open("r", encoding="utf-8") as source:
        for line_number, raw_line in enumerate(source, 1):
            line = raw_line.strip()
            if line:
                yield line_number, line


def parse_query_graph(path: Path) -> QueryGraph:
    """Parse the first graph in *path* using the matcher's simple-graph view."""
    graph_id: str | None = None
    labels_by_id: dict[int, str] = {}
    edge_label_sets: dict[tuple[int, int], set[str]] = {}
    raw_edge_records = 0
    duplicate_edge_records = 0
    self_loops = 0
    sorted_edges = True
    previous_edge: tuple[int, int] | None = None

    for line_number, line in _meaningful_lines(path):
        fields = line.split()
        record_type = fields[0]
        if record_type == "t":
            if graph_id is not None:
                break
            if len(fields) < 3 or fields[1] != "#":
                raise GraphAuditError(f"{path}:{line_number}: invalid graph header")
            graph_id = fields[2]
        elif graph_id is None:
            # Match the C++ loader: ignore preamble until the first t record.
            continue
        elif record_type == "v":
            if len(fields) < 3:
                raise GraphAuditError(f"{path}:{line_number}: invalid vertex record")
            vertex_id = int(fields[1])
            label = fields[2]
            if vertex_id in labels_by_id and labels_by_id[vertex_id] != label:
                raise GraphAuditError(
                    f"{path}:{line_number}: vertex {vertex_id} has conflicting labels"
                )
            labels_by_id[vertex_id] = label
        elif record_type == "e":
            if len(fields) < 3:
                raise GraphAuditError(f"{path}:{line_number}: invalid edge record")
            raw_edge_records += 1
            u, v = int(fields[1]), int(fields[2])
            edge_label = fields[3] if len(fields) >= 4 else NO_EDGE_LABEL
            if u == v:
                self_loops += 1
                continue
            edge = (u, v) if u < v else (v, u)
            if previous_edge is not None and edge < previous_edge:
                sorted_edges = False
            previous_edge = edge
            labels_for_edge = edge_label_sets.setdefault(edge, set())
            if edge_label in labels_for_edge:
                duplicate_edge_records += 1
            elif labels_for_edge:
                # Multiple labels on one simple edge are parallel records for the
                # matcher, whose Graph representation does not retain edge labels.
                duplicate_edge_records += 1
            labels_for_edge.add(edge_label)
        else:
            raise GraphAuditError(f"{path}:{line_number}: unknown record type {record_type!r}")

    if graph_id is None:
        raise GraphAuditError(f"{path}: graph header not found")
    if labels_by_id:
        expected_ids = set(range(max(labels_by_id) + 1))
        if set(labels_by_id) != expected_ids:
            missing = sorted(expected_ids - set(labels_by_id))[:10]
            raise GraphAuditError(f"{path}: vertex ids are not dense; missing {missing}")
    labels = tuple(labels_by_id[index] for index in range(len(labels_by_id)))
    adjacency = [set() for _ in labels]
    for u, v in edge_label_sets:
        if u >= len(labels) or v >= len(labels):
            raise GraphAuditError(f"{path}: edge ({u}, {v}) references an unknown vertex")
        adjacency[u].add(v)
        adjacency[v].add(u)

    return QueryGraph(
        path=path,
        graph_id=graph_id,
        labels=labels,
        adjacency=tuple(frozenset(neighbors) for neighbors in adjacency),
        edge_labels={edge: frozenset(values) for edge, values in edge_label_sets.items()},
        sha256=sha256_file(path),
        raw_edge_records=raw_edge_records,
        duplicate_edge_records=duplicate_edge_records,
        self_loops=self_loops,
        edge_records_sorted=sorted_edges,
    )


class _DisjointSet:
    def __init__(self) -> None:
        self.parent = array("i")
        self.size = array("I")
        self.components = 0
        self.maximum_size = 0

    def add(self) -> None:
        index = len(self.parent)
        self.parent.append(index)
        self.size.append(1)
        self.components += 1
        self.maximum_size = max(self.maximum_size, 1)

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]
        self.components -= 1
        self.maximum_size = max(self.maximum_size, self.size[left_root])


def _entropy_bits(counts: Iterable[int]) -> float:
    values = list(counts)
    total = sum(values)
    if total == 0:
        return 0.0
    return -sum((count / total) * math.log2(count / total) for count in values if count)


def _nearest_rank_percentiles(histogram: Counter[int], percentiles: Sequence[float]) -> list[int]:
    total = sum(histogram.values())
    if total == 0:
        return [0 for _ in percentiles]
    targets = [max(1, math.ceil(total * percentile)) for percentile in percentiles]
    answers: list[int] = []
    cumulative = 0
    target_index = 0
    for value in sorted(histogram):
        cumulative += histogram[value]
        while target_index < len(targets) and cumulative >= targets[target_index]:
            answers.append(value)
            target_index += 1
    return answers


def scan_data_graph(path: Path) -> DataGraphStats:
    """Stream a potentially large data graph and retain only bounded metadata.

    Duplicate detection is exact when edge endpoints are sorted, which is true
    for the repository datasets.  If order is not monotone, adjacent duplicates
    are still counted but the output explicitly marks the result as partial.
    """
    digest = hashlib.sha256()
    topology_digest = hashlib.sha256(b"SSM-GED undirected simple topology\0")
    graph_id: str | None = None
    label_ids: dict[str, int] = {}
    label_names: list[str] = []
    labels = array("I")
    label_frequencies: Counter[str] = Counter()
    degrees = array("Q")
    disjoint_set = _DisjointSet()
    label_pair_edges: set[tuple[str, str]] = set()
    labeled_edge_combinations: set[tuple[str, str, str]] = set()
    raw_edge_records = 0
    unique_edges = 0
    duplicate_edge_records = 0
    self_loops = 0
    sorted_edges = True
    previous_edge: tuple[int, int] | None = None
    in_graph = False

    with path.open("rb") as source:
        for line_number, raw_line in enumerate(source, 1):
            digest.update(raw_line)
            try:
                line = raw_line.decode("utf-8").strip()
            except UnicodeDecodeError as exc:
                raise GraphAuditError(f"{path}:{line_number}: not UTF-8") from exc
            if not line:
                continue
            fields = line.split()
            record_type = fields[0]
            if record_type == "t":
                if in_graph:
                    break
                if len(fields) < 3 or fields[1] != "#":
                    raise GraphAuditError(f"{path}:{line_number}: invalid graph header")
                graph_id = fields[2]
                in_graph = True
                continue
            if not in_graph:
                continue
            if record_type == "v":
                if len(fields) < 3:
                    raise GraphAuditError(f"{path}:{line_number}: invalid vertex record")
                vertex_id = int(fields[1])
                if vertex_id != len(labels):
                    raise GraphAuditError(
                        f"{path}:{line_number}: expected vertex id {len(labels)}, got {vertex_id}"
                    )
                label = fields[2]
                if label not in label_ids:
                    label_ids[label] = len(label_names)
                    label_names.append(label)
                labels.append(label_ids[label])
                label_frequencies[label] += 1
                degrees.append(0)
                disjoint_set.add()
            elif record_type == "e":
                if len(fields) < 3:
                    raise GraphAuditError(f"{path}:{line_number}: invalid edge record")
                raw_edge_records += 1
                u, v = int(fields[1]), int(fields[2])
                edge_label = fields[3] if len(fields) >= 4 else NO_EDGE_LABEL
                if u >= len(labels) or v >= len(labels):
                    raise GraphAuditError(
                        f"{path}:{line_number}: edge ({u}, {v}) references an unknown vertex"
                    )
                if u == v:
                    self_loops += 1
                    continue
                edge = (u, v) if u < v else (v, u)
                if previous_edge is not None and edge < previous_edge:
                    sorted_edges = False
                if edge == previous_edge:
                    duplicate_edge_records += 1
                    continue
                previous_edge = edge
                unique_edges += 1
                topology_digest.update(b"E")
                topology_digest.update(struct.pack("<QQ", edge[0], edge[1]))
                degrees[edge[0]] += 1
                degrees[edge[1]] += 1
                disjoint_set.union(*edge)
                left_label = label_names[labels[edge[0]]]
                right_label = label_names[labels[edge[1]]]
                label_pair = tuple(sorted((left_label, right_label)))
                label_pair_edges.add(label_pair)
                labeled_edge_combinations.add((label_pair[0], label_pair[1], edge_label))
            else:
                raise GraphAuditError(
                    f"{path}:{line_number}: unknown record type {record_type!r}"
                )

    if graph_id is None:
        raise GraphAuditError(f"{path}: graph header not found")

    degree_histogram = Counter(degrees)
    p50, p90, p95, p99 = _nearest_rank_percentiles(
        degree_histogram, (0.50, 0.90, 0.95, 0.99)
    )
    vertex_count = len(labels)
    topology_digest.update(b"V")
    topology_digest.update(struct.pack("<Q", vertex_count))
    mean_degree = (2.0 * unique_edges / vertex_count) if vertex_count else 0.0
    density = (
        2.0 * unique_edges / (vertex_count * (vertex_count - 1))
        if vertex_count > 1
        else 0.0
    )
    return DataGraphStats(
        path=path,
        graph_id=graph_id,
        sha256=digest.hexdigest(),
        topology_sha256=topology_digest.hexdigest(),
        vertex_count=vertex_count,
        edge_count=unique_edges,
        raw_edge_records=raw_edge_records,
        duplicate_edge_records=duplicate_edge_records,
        duplicate_detection="exact_sorted" if sorted_edges else "adjacent_only_unsorted",
        self_loops=self_loops,
        edge_records_sorted=sorted_edges,
        component_count=disjoint_set.components,
        largest_component_vertices=disjoint_set.maximum_size,
        isolated_vertices=degree_histogram.get(0, 0),
        mean_degree=mean_degree,
        max_degree=max(degree_histogram, default=0),
        degree_p50=p50,
        degree_p90=p90,
        degree_p95=p95,
        degree_p99=p99,
        density=density,
        label_frequencies=dict(label_frequencies),
        label_entropy_bits=_entropy_bits(label_frequencies.values()),
        label_pair_edges=frozenset(label_pair_edges),
        labeled_edge_combinations=frozenset(labeled_edge_combinations),
    )


def connected_components(graph: QueryGraph) -> list[set[int]]:
    unseen = set(range(graph.vertex_count))
    components: list[set[int]] = []
    while unseen:
        start = next(iter(unseen))
        component = {start}
        frontier = [start]
        unseen.remove(start)
        while frontier:
            vertex = frontier.pop()
            for neighbor in graph.adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    frontier.append(neighbor)
        components.append(component)
    return components


def triangle_count(graph: QueryGraph) -> int:
    total = 0
    for u in range(graph.vertex_count):
        for v in graph.adjacency[u]:
            if v > u:
                total += sum(1 for w in graph.adjacency[u] & graph.adjacency[v] if w > v)
    return total


def query_edge_label_pairs(graph: QueryGraph) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for u, v in graph.edge_labels:
        pairs.add(tuple(sorted((graph.labels[u], graph.labels[v]))))
    return pairs


def query_labeled_edge_combinations(graph: QueryGraph) -> set[tuple[str, str, str]]:
    combinations: set[tuple[str, str, str]] = set()
    for (u, v), edge_labels in graph.edge_labels.items():
        left, right = sorted((graph.labels[u], graph.labels[v]))
        combinations.update((left, right, edge_label) for edge_label in edge_labels)
    return combinations


def _wl_colors(graph: QueryGraph) -> tuple[str, ...]:
    colors = tuple(hashlib.sha256(("label\0" + label).encode()).hexdigest() for label in graph.labels)
    # A fixed number of rounds makes independently computed colors comparable.
    for _ in range(graph.vertex_count):
        refined: list[str] = []
        for vertex, own_color in enumerate(colors):
            digest = hashlib.sha256()
            digest.update(own_color.encode())
            digest.update(b"\0")
            for neighbor_color in sorted(colors[n] for n in graph.adjacency[vertex]):
                digest.update(neighbor_color.encode())
                digest.update(b"\0")
            refined.append(digest.hexdigest())
        colors = tuple(refined)
    return colors


def isomorphism_bucket(graph: QueryGraph) -> str:
    """Return an isomorphism-invariant bucket key (not an exact certificate)."""
    colors = _wl_colors(graph)
    description = {
        "n": graph.vertex_count,
        "m": graph.edge_count,
        "labels": sorted(Counter(graph.labels).items()),
        "vertices": sorted(
            (graph.labels[u], len(graph.adjacency[u]), colors[u])
            for u in range(graph.vertex_count)
        ),
    }
    return hashlib.sha256(json_compact(description).encode()).hexdigest()


def are_label_isomorphic(left: QueryGraph, right: QueryGraph) -> bool:
    """Exact vertex-label-preserving undirected graph isomorphism test."""
    if left.vertex_count != right.vertex_count or left.edge_count != right.edge_count:
        return False
    if Counter(left.labels) != Counter(right.labels):
        return False
    left_colors = _wl_colors(left)
    right_colors = _wl_colors(right)
    if Counter(left_colors) != Counter(right_colors):
        return False

    def key(graph: QueryGraph, colors: tuple[str, ...], vertex: int) -> tuple[str, int, str]:
        return graph.labels[vertex], len(graph.adjacency[vertex]), colors[vertex]

    right_by_key: dict[tuple[str, int, str], list[int]] = {}
    for vertex in range(right.vertex_count):
        right_by_key.setdefault(key(right, right_colors, vertex), []).append(vertex)
    if Counter(key(left, left_colors, u) for u in range(left.vertex_count)) != Counter(
        key(right, right_colors, v) for v in range(right.vertex_count)
    ):
        return False

    mapping: dict[int, int] = {}
    used_right: set[int] = set()

    def compatible(u: int, v: int) -> bool:
        for mapped_left, mapped_right in mapping.items():
            if (mapped_left in left.adjacency[u]) != (mapped_right in right.adjacency[v]):
                return False
        left_frontier = Counter(
            key(left, left_colors, neighbor)
            for neighbor in left.adjacency[u]
            if neighbor not in mapping
        )
        right_frontier = Counter(
            key(right, right_colors, neighbor)
            for neighbor in right.adjacency[v]
            if neighbor not in used_right
        )
        return left_frontier == right_frontier

    def search() -> bool:
        if len(mapping) == left.vertex_count:
            return True

        best_u = -1
        best_candidates: list[int] | None = None
        best_mapped_neighbors = -1
        for u in range(left.vertex_count):
            if u in mapping:
                continue
            candidates = [
                v
                for v in right_by_key[key(left, left_colors, u)]
                if v not in used_right and compatible(u, v)
            ]
            if not candidates:
                return False
            mapped_neighbors = sum(1 for neighbor in left.adjacency[u] if neighbor in mapping)
            if (
                best_candidates is None
                or len(candidates) < len(best_candidates)
                or (
                    len(candidates) == len(best_candidates)
                    and mapped_neighbors > best_mapped_neighbors
                )
            ):
                best_u = u
                best_candidates = candidates
                best_mapped_neighbors = mapped_neighbors

        assert best_candidates is not None
        for v in best_candidates:
            mapping[best_u] = v
            used_right.add(v)
            if search():
                return True
            used_right.remove(v)
            del mapping[best_u]
        return False

    return search()


def assign_isomorphism_classes(
    graphs: Sequence[QueryGraph], dataset_name: str
) -> tuple[dict[Path, str], Counter[str]]:
    """Assign exact isomorphism classes, using WL only as a safe prefilter."""
    buckets: dict[str, list[tuple[str, QueryGraph]]] = {}
    assignments: dict[Path, str] = {}
    class_sizes: Counter[str] = Counter()
    next_class = 1
    for graph in graphs:
        bucket = isomorphism_bucket(graph)
        matched_class: str | None = None
        for class_id, representative in buckets.setdefault(bucket, []):
            if are_label_isomorphic(graph, representative):
                matched_class = class_id
                break
        if matched_class is None:
            matched_class = f"{dataset_name}:iso:{next_class:04d}"
            next_class += 1
            buckets[bucket].append((matched_class, graph))
        assignments[graph.path] = matched_class
        class_sizes[matched_class] += 1
    return assignments, class_sizes
