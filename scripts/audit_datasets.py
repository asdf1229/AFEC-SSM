#!/usr/bin/env python3
"""Audit SSM-GED data graphs and query graphs before an experiment run."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import math
from pathlib import Path
import re
import sys
from typing import Iterable, Mapping

from graph_audit_common import (
    DataGraphStats,
    GraphAuditError,
    QueryGraph,
    assign_isomorphism_classes,
    connected_components,
    json_compact,
    parse_query_graph,
    query_edge_label_pairs,
    query_labeled_edge_combinations,
    scan_data_graph,
    triangle_count,
)


QUERY_GROUP_PATTERN = re.compile(r"^(?:(sparse|dense)_)?([1-9][0-9]*)$")


DATASET_COLUMNS = [
    "dataset_group",
    "dataset",
    "data_graph_path",
    "graph_id",
    "sha256",
    "content_duplicate_group",
    "content_duplicate_group_size",
    "topology_sha256",
    "topology_duplicate_group",
    "topology_duplicate_group_size",
    "vertex_count",
    "edge_count",
    "component_count",
    "largest_component_vertices",
    "largest_component_ratio",
    "isolated_vertices",
    "mean_degree",
    "max_degree",
    "degree_p50",
    "degree_p90",
    "degree_p95",
    "degree_p99",
    "density",
    "label_count",
    "label_frequencies",
    "label_entropy_bits",
    "raw_edge_records",
    "duplicate_edge_records",
    "duplicate_detection",
    "self_loops",
    "edge_records_sorted",
]


BASE_QUERY_COLUMNS = [
    "dataset_group",
    "dataset",
    "query_group",
    "query",
    "query_path",
    "data_graph_path",
    "graph_id",
    "sha256",
    "content_duplicate_group",
    "content_duplicate_group_size",
    "isomorphism_status",
    "isomorphism_class",
    "isomorphism_class_size",
    "vertex_count",
    "edge_count",
    "component_count",
    "largest_component_vertices",
    "isolated_vertices",
    "connected",
    "empty_graph",
    "single_vertex_graph",
    "mean_degree",
    "declared_density_class",
    "declared_vertex_count",
    "query_group_recognized",
    "query_group_size_ok",
    "query_group_density_ok",
    "query_group_design_ok",
    "max_degree",
    "density",
    "slack",
    "slack_bucket",
    "meets_slack_2",
    "meets_slack_4",
    "meets_slack_6",
    "label_count",
    "label_frequencies",
    "label_entropy_bits",
    "triangle_count",
    "raw_edge_records",
    "duplicate_edge_records",
    "self_loops",
    "edge_records_sorted",
    "missing_data_labels",
    "over_capacity_labels",
    "label_capacity_shortfall_vertices",
    "missing_data_label_pairs",
    "mandatory_missing_edges_by_label_pair",
    "missing_data_labeled_edge_combinations",
    "vertex_label_infeasible",
    "zero_result_static_risk",
    "planted_mapping_status",
    "perturbation_metadata_status",
    "source_overlap_status",
    "hard_error_reasons",
]


BASE_SUMMARY_COLUMNS = [
    "dataset_group",
    "dataset",
    "experiment_roles",
    "query_count",
    "expected_query_count",
    "query_count_ok",
    "query_group_count",
    "connected_queries",
    "disconnected_queries",
    "empty_queries",
    "single_vertex_queries",
    "hard_error_queries",
    "unique_file_contents",
    "content_duplicate_queries",
    "unique_isomorphism_classes",
    "isomorphic_duplicate_queries",
    "isomorphism_status",
    "slack_invalid",
    "slack_0",
    "slack_1",
    "slack_2",
    "slack_3",
    "slack_4",
    "slack_5",
    "slack_6_plus",
    "slack_ge_2",
    "slack_ge_4",
    "slack_ge_6",
    "query_group_design_queries",
    "unrecognized_query_group_queries",
    "query_group_size_violations",
    "query_group_density_violations",
    "query_group_design_violations",
    "static_zero_result_risk_queries",
    "gate_count_800",
    "gate_unique_file_content",
    "gate_unique_labeled_graph",
    "gate_connected",
    "gate_slack_2",
    "gate_slack_4",
    "gate_slack_6",
    "gate_query_group_design",
    "gate_static_label_feasible",
    "main_ready",
    "threshold_ready",
]


QUERY_GROUP_SUMMARY_COLUMNS = [
    "dataset_group",
    "dataset",
    "query_group",
    "query_count",
    "declared_density_class",
    "declared_vertex_count",
    "min_vertex_count",
    "max_vertex_count",
    "min_edge_count",
    "max_edge_count",
    "min_mean_degree",
    "max_mean_degree",
    "min_slack",
    "max_slack",
    "slack_ge_4",
    "size_violations",
    "density_violations",
    "slack_lt_4",
    "design_violations",
    "gate_size_density",
    "gate_slack_4",
    "group_ready",
]


def _parse_csv_set(value: str) -> set[str]:
    return {item.strip() for item in value.split(",") if item.strip()}


def _discover_datasets(data_root: Path) -> list[tuple[str, str, Path]]:
    roots: list[tuple[str, Path]] = []
    if (data_root / "synthetic").is_dir():
        roots.append(("synthetic", data_root / "synthetic"))
    if (data_root / "real_graphs").is_dir():
        roots.append(("real", data_root / "real_graphs"))
    if not roots:
        roots.append(("custom", data_root))

    datasets: list[tuple[str, str, Path]] = []
    for dataset_group, root in roots:
        for graph_path in sorted(root.rglob("graph_g.txt")):
            if "query_graph" in graph_path.parts:
                continue
            dataset_dir = graph_path.parent
            datasets.append((dataset_group, dataset_dir.name, dataset_dir))
    return sorted(datasets, key=lambda item: (item[0], item[1], str(item[2])))


def _discover_queries(dataset_dir: Path) -> list[tuple[str, str, Path]]:
    query_root = dataset_dir / "query_graph"
    if not query_root.is_dir():
        return []
    queries: list[tuple[str, str, Path]] = []
    for query_path in sorted(query_root.rglob("*.txt"), key=lambda path: str(path)):
        relative = query_path.relative_to(query_root)
        query_group = str(relative.parent) if relative.parent != Path(".") else "legacy"
        queries.append((query_group, query_path.stem, query_path))
    return queries


def _entropy_bits(counts: Iterable[int]) -> float:
    values = list(counts)
    total = sum(values)
    if total == 0:
        return 0.0
    return -sum((value / total) * math.log2(value / total) for value in values if value)


def _format_float(value: float) -> str:
    return f"{value:.12g}"


def _slack_bucket(slack: int | None) -> str:
    if slack is None:
        return "invalid"
    return str(slack) if slack < 6 else "6+"


def _query_record(
    dataset_group: str,
    dataset: str,
    query_group: str,
    query_name: str,
    graph: QueryGraph,
    data_stats: DataGraphStats,
    theta_max: int,
) -> dict[str, object]:
    components = connected_components(graph)
    connected = graph.vertex_count > 0 and len(components) == 1
    slack = graph.edge_count - graph.vertex_count + 1 if connected else None
    mean_degree = 2.0 * graph.edge_count / graph.vertex_count if graph.vertex_count else 0.0
    query_group_name = Path(query_group).name
    group_match = QUERY_GROUP_PATTERN.fullmatch(query_group_name)
    if group_match:
        declared_density_class = group_match.group(1) or "any"
        declared_vertex_count = int(group_match.group(2))
        group_size_ok: int | str = int(graph.vertex_count == declared_vertex_count)
        if declared_density_class == "sparse":
            # Integer comparison makes the inclusive d <= 3 boundary exact.
            group_density_ok: int | str = int(
                2 * graph.edge_count <= 3 * graph.vertex_count
            )
        elif declared_density_class == "dense":
            group_density_ok = int(2 * graph.edge_count > 3 * graph.vertex_count)
        else:
            group_density_ok = 1
        group_design_ok: int | str = int(bool(group_size_ok) and bool(group_density_ok))
    else:
        declared_density_class = "unrecognized"
        declared_vertex_count = "NA"
        group_size_ok = "NA"
        group_density_ok = "NA"
        group_design_ok = "NA"
    label_frequencies = Counter(graph.labels)
    missing_labels = sorted(set(label_frequencies) - set(data_stats.label_frequencies))
    over_capacity = {
        label: demand - data_stats.label_frequencies.get(label, 0)
        for label, demand in label_frequencies.items()
        if demand > data_stats.label_frequencies.get(label, 0)
    }
    missing_label_pairs = sorted(query_edge_label_pairs(graph) - data_stats.label_pair_edges)
    mandatory_missing_edges = sum(
        tuple(sorted((graph.labels[u], graph.labels[v]))) not in data_stats.label_pair_edges
        for u, v in graph.edge_labels
    )
    missing_labeled_edges = sorted(
        query_labeled_edge_combinations(graph) - data_stats.labeled_edge_combinations
    )
    # Vertex-label capacity is threshold-independent and therefore proves a
    # zero-result query.  A missing endpoint-label pair only proves that the
    # corresponding query edge must consume theta; it is not by itself a zero
    # result under this project's missing-query-edge semantics.
    zero_result_risk = bool(over_capacity)

    hard_errors: list[str] = []
    if graph.vertex_count == 0:
        hard_errors.append("empty_graph")
    elif graph.vertex_count == 1:
        hard_errors.append("single_vertex_graph")
    if not connected:
        hard_errors.append("disconnected")
    if graph.self_loops:
        hard_errors.append("self_loops")
    if graph.duplicate_edge_records:
        hard_errors.append("duplicate_edges")
    if over_capacity:
        hard_errors.append("label_capacity")

    record: dict[str, object] = {
        "dataset_group": dataset_group,
        "dataset": dataset,
        "query_group": query_group,
        "query": query_name,
        "query_path": str(graph.path.resolve()),
        "data_graph_path": str(data_stats.path.resolve()),
        "graph_id": graph.graph_id,
        "sha256": graph.sha256,
        "content_duplicate_group": "",
        "content_duplicate_group_size": 1,
        "isomorphism_status": "pending",
        "isomorphism_class": "",
        "isomorphism_class_size": "",
        "vertex_count": graph.vertex_count,
        "edge_count": graph.edge_count,
        "component_count": len(components),
        "largest_component_vertices": max((len(component) for component in components), default=0),
        "isolated_vertices": sum(1 for neighbors in graph.adjacency if not neighbors),
        "connected": int(connected),
        "empty_graph": int(graph.vertex_count == 0),
        "single_vertex_graph": int(graph.vertex_count == 1),
        "mean_degree": _format_float(mean_degree),
        "declared_density_class": declared_density_class,
        "declared_vertex_count": declared_vertex_count,
        "query_group_recognized": int(group_match is not None),
        "query_group_size_ok": group_size_ok,
        "query_group_density_ok": group_density_ok,
        "query_group_design_ok": group_design_ok,
        "max_degree": max((len(neighbors) for neighbors in graph.adjacency), default=0),
        "density": _format_float(
            2.0 * graph.edge_count / (graph.vertex_count * (graph.vertex_count - 1))
            if graph.vertex_count > 1
            else 0.0
        ),
        "slack": slack if slack is not None else "NA",
        "slack_bucket": _slack_bucket(slack),
        "meets_slack_2": int(slack is not None and slack >= 2),
        "meets_slack_4": int(slack is not None and slack >= 4),
        "meets_slack_6": int(slack is not None and slack >= 6),
        "label_count": len(label_frequencies),
        "label_frequencies": json_compact(dict(label_frequencies)),
        "label_entropy_bits": _format_float(_entropy_bits(label_frequencies.values())),
        "triangle_count": triangle_count(graph),
        "raw_edge_records": graph.raw_edge_records,
        "duplicate_edge_records": graph.duplicate_edge_records,
        "self_loops": graph.self_loops,
        "edge_records_sorted": int(graph.edge_records_sorted),
        "missing_data_labels": json_compact(missing_labels),
        "over_capacity_labels": json_compact(over_capacity),
        "label_capacity_shortfall_vertices": sum(over_capacity.values()),
        "missing_data_label_pairs": json_compact(missing_label_pairs),
        "mandatory_missing_edges_by_label_pair": mandatory_missing_edges,
        "missing_data_labeled_edge_combinations": json_compact(missing_labeled_edges),
        "vertex_label_infeasible": int(bool(over_capacity)),
        "zero_result_static_risk": int(zero_result_risk),
        "planted_mapping_status": "not_provided",
        "perturbation_metadata_status": "not_provided",
        "source_overlap_status": "not_checkable_without_source_mapping",
        "hard_error_reasons": json_compact(hard_errors),
        "_hard_error": bool(hard_errors),
        "_slack_value": slack,
        "_graph": graph,
    }
    for theta in range(1, theta_max + 1):
        effective = min(theta, slack) if slack is not None else None
        record[f"effective_theta_{theta}"] = effective if effective is not None else "NA"
        record[f"theta_{theta}_truncated"] = int(slack is not None and theta > slack)
        record[f"theta_{theta}_static_infeasible"] = int(
            bool(over_capacity)
            or effective is None
            or mandatory_missing_edges > effective
        )
    return record


def _dataset_record(
    dataset_group: str,
    dataset: str,
    stats: DataGraphStats,
    duplicate_group: str,
    duplicate_size: int,
    topology_duplicate_group: str,
    topology_duplicate_size: int,
) -> dict[str, object]:
    return {
        "dataset_group": dataset_group,
        "dataset": dataset,
        "data_graph_path": str(stats.path.resolve()),
        "graph_id": stats.graph_id,
        "sha256": stats.sha256,
        "content_duplicate_group": duplicate_group,
        "content_duplicate_group_size": duplicate_size,
        "topology_sha256": stats.topology_sha256,
        "topology_duplicate_group": topology_duplicate_group,
        "topology_duplicate_group_size": topology_duplicate_size,
        "vertex_count": stats.vertex_count,
        "edge_count": stats.edge_count,
        "component_count": stats.component_count,
        "largest_component_vertices": stats.largest_component_vertices,
        "largest_component_ratio": _format_float(
            stats.largest_component_vertices / stats.vertex_count if stats.vertex_count else 0.0
        ),
        "isolated_vertices": stats.isolated_vertices,
        "mean_degree": _format_float(stats.mean_degree),
        "max_degree": stats.max_degree,
        "degree_p50": stats.degree_p50,
        "degree_p90": stats.degree_p90,
        "degree_p95": stats.degree_p95,
        "degree_p99": stats.degree_p99,
        "density": _format_float(stats.density),
        "label_count": len(stats.label_frequencies),
        "label_frequencies": json_compact(stats.label_frequencies),
        "label_entropy_bits": _format_float(stats.label_entropy_bits),
        "raw_edge_records": stats.raw_edge_records,
        "duplicate_edge_records": stats.duplicate_edge_records,
        "duplicate_detection": stats.duplicate_detection,
        "self_loops": stats.self_loops,
        "edge_records_sorted": int(stats.edge_records_sorted),
    }


def _summarize_dataset(
    dataset_group: str,
    dataset: str,
    records: list[dict[str, object]],
    expected_queries: int,
    theta_max: int,
    roles: list[str],
) -> dict[str, object]:
    query_count = len(records)
    hashes = Counter(str(record["sha256"]) for record in records)
    iso_status = "exact" if all(record["isomorphism_status"] == "exact" for record in records) else "skipped"
    iso_classes = Counter(
        str(record["isomorphism_class"])
        for record in records
        if record["isomorphism_class"]
    )
    slack_counts = Counter(record["slack_bucket"] for record in records)
    connected_count = sum(int(record["connected"]) for record in records)
    static_risk = sum(int(record["zero_result_static_risk"]) for record in records)
    design_records = [record for record in records if int(record["query_group_recognized"])]
    unrecognized_groups = query_count - len(design_records)
    size_violations = sum(not int(record["query_group_size_ok"]) for record in design_records)
    density_violations = sum(
        not int(record["query_group_density_ok"]) for record in design_records
    )
    design_violations = sum(
        not int(record["query_group_design_ok"]) for record in design_records
    )
    unique_iso = len(iso_classes) if iso_status == "exact" else "NA"
    iso_duplicates = query_count - len(iso_classes) if iso_status == "exact" else "NA"

    count_gate = query_count == expected_queries
    content_gate = len(hashes) == query_count
    iso_gate = iso_status == "exact" and len(iso_classes) == query_count
    connected_gate = connected_count == query_count
    slack_2_gate = sum(int(record["meets_slack_2"]) for record in records) == query_count
    slack_4_gate = sum(int(record["meets_slack_4"]) for record in records) == query_count
    slack_6_gate = sum(int(record["meets_slack_6"]) for record in records) == query_count
    group_design_gate = unrecognized_groups == 0 and design_violations == 0
    label_gate = static_risk == 0

    summary: dict[str, object] = {
        "dataset_group": dataset_group,
        "dataset": dataset,
        "experiment_roles": ",".join(roles) if roles else "unassigned",
        "query_count": query_count,
        "expected_query_count": expected_queries,
        "query_count_ok": int(count_gate),
        "query_group_count": len({str(record["query_group"]) for record in records}),
        "connected_queries": connected_count,
        "disconnected_queries": query_count - connected_count,
        "empty_queries": sum(int(record["empty_graph"]) for record in records),
        "single_vertex_queries": sum(int(record["single_vertex_graph"]) for record in records),
        "hard_error_queries": sum(bool(record["_hard_error"]) for record in records),
        "unique_file_contents": len(hashes),
        "content_duplicate_queries": query_count - len(hashes),
        "unique_isomorphism_classes": unique_iso,
        "isomorphic_duplicate_queries": iso_duplicates,
        "isomorphism_status": iso_status,
        "slack_invalid": slack_counts["invalid"],
        "slack_0": slack_counts["0"],
        "slack_1": slack_counts["1"],
        "slack_2": slack_counts["2"],
        "slack_3": slack_counts["3"],
        "slack_4": slack_counts["4"],
        "slack_5": slack_counts["5"],
        "slack_6_plus": slack_counts["6+"],
        "slack_ge_2": sum(int(record["meets_slack_2"]) for record in records),
        "slack_ge_4": sum(int(record["meets_slack_4"]) for record in records),
        "slack_ge_6": sum(int(record["meets_slack_6"]) for record in records),
        "query_group_design_queries": len(design_records),
        "unrecognized_query_group_queries": unrecognized_groups,
        "query_group_size_violations": size_violations,
        "query_group_density_violations": density_violations,
        "query_group_design_violations": design_violations,
        "static_zero_result_risk_queries": static_risk,
        "gate_count_800": int(count_gate),
        "gate_unique_file_content": int(content_gate),
        "gate_unique_labeled_graph": int(iso_gate),
        "gate_connected": int(connected_gate),
        "gate_slack_2": int(slack_2_gate),
        "gate_slack_4": int(slack_4_gate),
        "gate_slack_6": int(slack_6_gate),
        "gate_query_group_design": int(group_design_gate),
        "gate_static_label_feasible": int(label_gate),
        "main_ready": int(
            count_gate
            and content_gate
            and iso_gate
            and connected_gate
            and slack_4_gate
            and group_design_gate
            and label_gate
        ),
        "threshold_ready": int(
            count_gate
            and content_gate
            and iso_gate
            and connected_gate
            and slack_6_gate
            and group_design_gate
            and label_gate
        ),
    }
    for theta in range(1, theta_max + 1):
        truncated = sum(int(record[f"theta_{theta}_truncated"]) for record in records)
        distribution = Counter(str(record[f"effective_theta_{theta}"]) for record in records)
        static_infeasible = sum(
            int(record[f"theta_{theta}_static_infeasible"]) for record in records
        )
        summary[f"theta_{theta}_truncated_queries"] = truncated
        summary[f"theta_{theta}_truncated_ratio"] = _format_float(
            truncated / query_count if query_count else 0.0
        )
        summary[f"effective_theta_{theta}_distribution"] = json_compact(dict(distribution))
        summary[f"theta_{theta}_static_infeasible_queries"] = static_infeasible
    return summary


def _summarize_query_group(
    dataset_group: str,
    dataset: str,
    query_group: str,
    records: list[dict[str, object]],
) -> dict[str, object]:
    vertex_counts = [int(record["vertex_count"]) for record in records]
    edge_counts = [int(record["edge_count"]) for record in records]
    mean_degrees = [float(record["mean_degree"]) for record in records]
    slack_values = [
        int(record["_slack_value"])
        for record in records
        if record["_slack_value"] is not None
    ]
    size_violations = sum(
        not int(record["query_group_size_ok"])
        for record in records
        if int(record["query_group_recognized"])
    )
    density_violations = sum(
        not int(record["query_group_density_ok"])
        for record in records
        if int(record["query_group_recognized"])
    )
    design_violations = sum(
        not int(record["query_group_design_ok"])
        for record in records
        if int(record["query_group_recognized"])
    )
    recognized = all(int(record["query_group_recognized"]) for record in records)
    slack_ge_4 = sum(int(record["meets_slack_4"]) for record in records)
    group_design_gate = recognized and design_violations == 0
    slack_gate = slack_ge_4 == len(records)
    return {
        "dataset_group": dataset_group,
        "dataset": dataset,
        "query_group": query_group,
        "query_count": len(records),
        "declared_density_class": records[0]["declared_density_class"],
        "declared_vertex_count": records[0]["declared_vertex_count"],
        "min_vertex_count": min(vertex_counts),
        "max_vertex_count": max(vertex_counts),
        "min_edge_count": min(edge_counts),
        "max_edge_count": max(edge_counts),
        "min_mean_degree": _format_float(min(mean_degrees)),
        "max_mean_degree": _format_float(max(mean_degrees)),
        "min_slack": min(slack_values) if slack_values else "NA",
        "max_slack": max(slack_values) if slack_values else "NA",
        "slack_ge_4": slack_ge_4,
        "size_violations": size_violations,
        "density_violations": density_violations,
        "slack_lt_4": len(records) - slack_ge_4,
        "design_violations": design_violations,
        "gate_size_density": int(group_design_gate),
        "gate_slack_4": int(slack_gate),
        "group_ready": int(group_design_gate and slack_gate),
    }


def _write_tsv(path: Path, columns: list[str], records: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("test/datasets"),
        help="Dataset root (default: test/datasets).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data_audit"),
        help="Directory for TSV output (default: data_audit).",
    )
    parser.add_argument(
        "--datasets",
        default="",
        help="Optional comma-separated dataset names to scan.",
    )
    parser.add_argument(
        "--expected-queries",
        type=int,
        default=800,
        help="Expected queries per dataset (default: 800).",
    )
    parser.add_argument(
        "--theta-max",
        type=int,
        default=6,
        help="Largest requested theta to inventory (default: 6).",
    )
    parser.add_argument(
        "--main-datasets",
        default="",
        help="Comma-separated datasets designated for the theta=2 main/ablation runs.",
    )
    parser.add_argument(
        "--threshold-datasets",
        default="",
        help="Comma-separated datasets designated for the theta=1...6 run.",
    )
    parser.add_argument(
        "--skip-isomorphism",
        action="store_true",
        help="Skip exact labeled-isomorphism deduplication (faster, gate remains unverified).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit nonzero when a designated experiment role fails its readiness gates.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.expected_queries <= 0:
        raise SystemExit("--expected-queries must be positive")
    if args.theta_max <= 0:
        raise SystemExit("--theta-max must be positive")

    selected = _parse_csv_set(args.datasets)
    main_datasets = _parse_csv_set(args.main_datasets)
    threshold_datasets = _parse_csv_set(args.threshold_datasets)
    discovered = _discover_datasets(args.data_root.resolve())
    if selected:
        discovered = [entry for entry in discovered if entry[1] in selected]
    if not discovered:
        raise SystemExit(f"no datasets found under {args.data_root}")
    discovered_names = {dataset for _, dataset, _ in discovered}
    missing_roles = (main_datasets | threshold_datasets) - discovered_names
    if missing_roles:
        raise SystemExit(f"designated datasets not found: {', '.join(sorted(missing_roles))}")

    data_stats_entries: list[tuple[str, str, DataGraphStats]] = []
    all_query_records: list[dict[str, object]] = []
    records_by_dataset: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    sha_records: list[dict[str, object]] = []

    for dataset_group, dataset, dataset_dir in discovered:
        graph_path = dataset_dir / "graph_g.txt"
        print(f"[data] {dataset_group}/{dataset}", file=sys.stderr, flush=True)
        try:
            data_stats = scan_data_graph(graph_path)
        except (OSError, GraphAuditError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        data_stats_entries.append((dataset_group, dataset, data_stats))
        sha_records.append(
            {
                "artifact_type": "data_graph",
                "dataset_group": dataset_group,
                "dataset": dataset,
                "query_group": "",
                "query": "",
                "path": str(graph_path.resolve()),
                "sha256": data_stats.sha256,
            }
        )

        query_specs = _discover_queries(dataset_dir)
        graphs: list[QueryGraph] = []
        for query_index, (query_group, query_name, query_path) in enumerate(query_specs, 1):
            if query_index == 1 or query_index % 100 == 0:
                print(
                    f"[query] {dataset_group}/{dataset} {query_index}/{len(query_specs)}",
                    file=sys.stderr,
                    flush=True,
                )
            try:
                graph = parse_query_graph(query_path)
            except (OSError, GraphAuditError) as exc:
                print(f"error: {exc}", file=sys.stderr)
                return 2
            graphs.append(graph)
            record = _query_record(
                dataset_group,
                dataset,
                query_group,
                query_name,
                graph,
                data_stats,
                args.theta_max,
            )
            records_by_dataset[(dataset_group, dataset)].append(record)
            all_query_records.append(record)
            sha_records.append(
                {
                    "artifact_type": "query_graph",
                    "dataset_group": dataset_group,
                    "dataset": dataset,
                    "query_group": query_group,
                    "query": query_name,
                    "path": str(query_path.resolve()),
                    "sha256": graph.sha256,
                }
            )

        dataset_records = records_by_dataset[(dataset_group, dataset)]
        hashes = Counter(str(record["sha256"]) for record in dataset_records)
        for record in dataset_records:
            duplicate_size = hashes[str(record["sha256"])]
            record["content_duplicate_group_size"] = duplicate_size
            if duplicate_size > 1:
                record["content_duplicate_group"] = (
                    f"{dataset}:sha256:{str(record['sha256'])[:16]}"
                )

        if args.skip_isomorphism:
            for record in dataset_records:
                record["isomorphism_status"] = "skipped"
                record["isomorphism_class"] = ""
                record["isomorphism_class_size"] = ""
        else:
            print(f"[isomorphism] {dataset_group}/{dataset}", file=sys.stderr, flush=True)
            assignments, class_sizes = assign_isomorphism_classes(graphs, dataset)
            for record in dataset_records:
                graph = record["_graph"]
                assert isinstance(graph, QueryGraph)
                class_id = assignments[graph.path]
                record["isomorphism_status"] = "exact"
                record["isomorphism_class"] = class_id
                record["isomorphism_class_size"] = class_sizes[class_id]

    data_hashes = Counter(stats.sha256 for _, _, stats in data_stats_entries)
    topology_hashes = Counter(stats.topology_sha256 for _, _, stats in data_stats_entries)
    dataset_rows = []
    for dataset_group, dataset, stats in data_stats_entries:
        duplicate_size = data_hashes[stats.sha256]
        duplicate_group = (
            f"data:sha256:{stats.sha256[:16]}" if duplicate_size > 1 else ""
        )
        topology_duplicate_size = topology_hashes[stats.topology_sha256]
        topology_duplicate_group = (
            f"topology:sha256:{stats.topology_sha256[:16]}"
            if topology_duplicate_size > 1
            else ""
        )
        dataset_rows.append(
            _dataset_record(
                dataset_group,
                dataset,
                stats,
                duplicate_group,
                duplicate_size,
                topology_duplicate_group,
                topology_duplicate_size,
            )
        )

    summary_rows: list[dict[str, object]] = []
    for dataset_group, dataset, _ in data_stats_entries:
        roles = []
        if dataset in main_datasets:
            roles.extend(("main", "ablation"))
        if dataset in threshold_datasets:
            roles.append("threshold")
        summary_rows.append(
            _summarize_dataset(
                dataset_group,
                dataset,
                records_by_dataset[(dataset_group, dataset)],
                args.expected_queries,
                args.theta_max,
                roles,
            )
        )

    records_by_query_group: dict[
        tuple[str, str, str], list[dict[str, object]]
    ] = defaultdict(list)
    for record in all_query_records:
        key = (
            str(record["dataset_group"]),
            str(record["dataset"]),
            str(record["query_group"]),
        )
        records_by_query_group[key].append(record)
    query_group_rows = [
        _summarize_query_group(dataset_group, dataset, query_group, records)
        for (dataset_group, dataset, query_group), records in sorted(
            records_by_query_group.items()
        )
    ]

    theta_columns: list[str] = []
    summary_theta_columns: list[str] = []
    for theta in range(1, args.theta_max + 1):
        theta_columns.extend(
            (
                f"effective_theta_{theta}",
                f"theta_{theta}_truncated",
                f"theta_{theta}_static_infeasible",
            )
        )
        summary_theta_columns.extend(
            (
                f"theta_{theta}_truncated_queries",
                f"theta_{theta}_truncated_ratio",
                f"effective_theta_{theta}_distribution",
                f"theta_{theta}_static_infeasible_queries",
            )
        )

    output_dir = args.output_dir.resolve()
    _write_tsv(output_dir / "dataset_inventory.tsv", DATASET_COLUMNS, dataset_rows)
    _write_tsv(
        output_dir / "query_inventory.tsv",
        BASE_QUERY_COLUMNS + theta_columns,
        all_query_records,
    )
    _write_tsv(
        output_dir / "data_audit_summary.tsv",
        BASE_SUMMARY_COLUMNS + summary_theta_columns,
        summary_rows,
    )
    _write_tsv(
        output_dir / "query_group_audit_summary.tsv",
        QUERY_GROUP_SUMMARY_COLUMNS,
        query_group_rows,
    )
    _write_tsv(
        output_dir / "sha256_manifest.tsv",
        ["artifact_type", "dataset_group", "dataset", "query_group", "query", "path", "sha256"],
        sha_records,
    )

    print(f"wrote audit files to {output_dir}", file=sys.stderr)
    for summary in summary_rows:
        print(
            f"{summary['dataset']}: queries={summary['query_count']} "
            f"unique_files={summary['unique_file_contents']} "
            f"unique_iso={summary['unique_isomorphism_classes']} "
            f"connected={summary['connected_queries']} "
            f"slack>=4={summary['slack_ge_4']} slack>=6={summary['slack_ge_6']} "
            f"group_design_errors={summary['query_group_design_violations']}",
            file=sys.stderr,
        )

    if args.strict:
        failed = False
        for summary in summary_rows:
            dataset = str(summary["dataset"])
            if dataset in main_datasets and not int(summary["main_ready"]):
                failed = True
            if dataset in threshold_datasets and not int(summary["threshold_ready"]):
                failed = True
        if failed:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
