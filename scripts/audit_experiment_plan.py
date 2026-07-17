#!/usr/bin/env python3
"""Verify that algorithms, thresholds, and ablations use identical graph files."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable, Mapping

from graph_audit_common import sha256_file


EXPECTED_DATASETS = {"main": 10, "threshold": 2, "ablation": 10}
EXPECTED_THETAS = {
    "main": frozenset({2}),
    "threshold": frozenset(range(1, 7)),
    "ablation": frozenset({2}),
}


@dataclass(frozen=True)
class PlanRow:
    experiment: str
    variant: str
    dataset: str
    query_id: str
    data_graph: Path
    query_graph: Path
    requested_theta: int
    preprocessing_artifact: Path | None
    source: str


def _write_tsv(path: Path, columns: list[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_input_spec(spec: str) -> tuple[str | None, Path]:
    if "=" not in spec:
        return None, Path(spec)
    experiment, path_text = spec.split("=", 1)
    if experiment not in {"main", "threshold", "ablation", "generic"}:
        raise ValueError(
            f"invalid experiment prefix {experiment!r}; expected main, threshold, ablation, or generic"
        )
    return experiment, Path(path_text)


def _canonical_value(row: Mapping[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name, "").strip()
        if value:
            return value
    return ""


def _read_canonical(path: Path, override: str | None) -> list[PlanRow]:
    rows: list[PlanRow] = []
    with path.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            return rows
        for line_number, row in enumerate(reader, 2):
            experiment = override or _canonical_value(row, "experiment", "experiment_type")
            variant = _canonical_value(row, "variant", "algorithm", "case")
            dataset = _canonical_value(row, "dataset")
            query_id = _canonical_value(row, "query_id", "query")
            data_graph = _canonical_value(row, "data_graph", "data_graph_path", "gfile")
            query_graph = _canonical_value(row, "query_graph", "query_graph_path", "qfile")
            theta_text = _canonical_value(row, "requested_theta", "theta", "threshold")
            preprocessing_text = _canonical_value(
                row, "preprocessing_artifact", "preprocessing_path"
            )
            missing = [
                name
                for name, value in (
                    ("experiment", experiment),
                    ("variant", variant),
                    ("dataset", dataset),
                    ("query_id", query_id),
                    ("data_graph", data_graph),
                    ("query_graph", query_graph),
                    ("requested_theta", theta_text),
                )
                if not value
            ]
            if missing:
                raise ValueError(f"{path}:{line_number}: missing {', '.join(missing)}")
            rows.append(
                PlanRow(
                    experiment=experiment,
                    variant=variant,
                    dataset=dataset,
                    query_id=query_id,
                    data_graph=Path(data_graph).expanduser().resolve(),
                    query_graph=Path(query_graph).expanduser().resolve(),
                    requested_theta=int(theta_text),
                    preprocessing_artifact=(
                        Path(preprocessing_text).expanduser().resolve()
                        if preprocessing_text
                        else None
                    ),
                    source=f"{path}:{line_number}",
                )
            )
    return rows


def _read_compare_tasks(path: Path, experiment: str | None) -> list[PlanRow]:
    rows: list[PlanRow] = []
    with path.open("r", encoding="utf-8", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if not fields:
                continue
            if len(fields) != 9:
                raise ValueError(
                    f"{path}:{line_number}: compare tasks row must have 9 columns, got {len(fields)}"
                )
            _, _, _, algorithm, data_text, query_text, theta_text, _, _ = fields
            data_graph = Path(data_text).expanduser().resolve()
            query_graph = Path(query_text).expanduser().resolve()
            query_root = data_graph.parent / "query_graph"
            try:
                query_id = str(query_graph.relative_to(query_root).with_suffix(""))
            except ValueError:
                query_id = query_graph.stem
            rows.append(
                PlanRow(
                    experiment=experiment or "generic",
                    variant=algorithm,
                    dataset=data_graph.parent.name,
                    query_id=query_id,
                    data_graph=data_graph,
                    query_graph=query_graph,
                    requested_theta=int(theta_text),
                    preprocessing_artifact=None,
                    source=f"{path}:{line_number}",
                )
            )
    return rows


def _read_input(spec: str) -> list[PlanRow]:
    override, path = _parse_input_spec(spec)
    path = path.expanduser().resolve()
    with path.open("r", encoding="utf-8") as source:
        first_line = source.readline().rstrip("\n")
    first_fields = first_line.split("\t") if first_line else []
    canonical_markers = {"experiment", "variant", "data_graph", "data_graph_path"}
    if canonical_markers & set(first_fields):
        return _read_canonical(path, override)
    return _read_compare_tasks(path, override)


def _issue(
    issues: list[dict[str, object]],
    severity: str,
    code: str,
    experiment: str = "",
    variant: str = "",
    dataset: str = "",
    query_id: str = "",
    details: str = "",
) -> None:
    issues.append(
        {
            "severity": severity,
            "code": code,
            "experiment": experiment,
            "variant": variant,
            "dataset": dataset,
            "query_id": query_id,
            "details": details,
        }
    )


def _ordered_unique(values: Iterable[tuple[str, str]]) -> list[tuple[str, str]]:
    seen: set[tuple[str, str]] = set()
    result: list[tuple[str, str]] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            result.append(value)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="[EXPERIMENT=]PLAN.tsv",
        help=(
            "Canonical headered plan or compare.sh tasks.tsv. Prefix headerless tasks with "
            "main=, threshold=, ablation=, or generic=."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("experiment_plan_audit"),
        help="Output directory (default: experiment_plan_audit).",
    )
    parser.add_argument(
        "--expected-queries",
        type=int,
        default=800,
        help="Expected unique queries per dataset and variant (default: 800).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit nonzero when any ERROR is found.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.expected_queries <= 0:
        raise SystemExit("--expected-queries must be positive")
    try:
        rows = [row for input_spec in args.inputs for row in _read_input(input_spec)]
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if not rows:
        print("error: plan contains no rows", file=sys.stderr)
        return 2

    issues: list[dict[str, object]] = []
    hash_cache: dict[Path, str] = {}
    preprocessing_hashes: dict[str, str] = {}
    resolved_rows: list[dict[str, object]] = []
    usable_rows: list[tuple[PlanRow, str, str]] = []
    for row in rows:
        required_paths = [row.data_graph, row.query_graph]
        if row.preprocessing_artifact is not None:
            required_paths.append(row.preprocessing_artifact)
        missing = [str(path) for path in required_paths if not path.is_file()]
        if missing:
            _issue(
                issues,
                "ERROR",
                "FILE_MISSING",
                row.experiment,
                row.variant,
                row.dataset,
                row.query_id,
                ", ".join(missing),
            )
            continue
        for path in required_paths:
            if path not in hash_cache:
                hash_cache[path] = sha256_file(path)
        data_sha = hash_cache[row.data_graph]
        query_sha = hash_cache[row.query_graph]
        preprocessing_sha = (
            hash_cache[row.preprocessing_artifact]
            if row.preprocessing_artifact is not None
            else ""
        )
        preprocessing_hashes[row.source] = preprocessing_sha
        usable_rows.append((row, data_sha, query_sha))
        resolved_rows.append(
            {
                "experiment": row.experiment,
                "variant": row.variant,
                "dataset": row.dataset,
                "query_id": row.query_id,
                "requested_theta": row.requested_theta,
                "data_graph": str(row.data_graph),
                "data_sha256": data_sha,
                "query_graph": str(row.query_graph),
                "query_sha256": query_sha,
                "preprocessing_artifact": (
                    str(row.preprocessing_artifact) if row.preprocessing_artifact else ""
                ),
                "preprocessing_sha256": preprocessing_sha,
                "source": row.source,
            }
        )

    exact_rows = Counter(
        (
            row.experiment,
            row.variant,
            row.dataset,
            row.query_id,
            row.requested_theta,
            data_sha,
            query_sha,
        )
        for row, data_sha, query_sha in usable_rows
    )
    for key, count in exact_rows.items():
        if count > 1:
            experiment, variant, dataset, query_id, theta, _, _ = key
            _issue(
                issues,
                "ERROR",
                "DUPLICATE_PLAN_ROW",
                experiment,
                variant,
                dataset,
                query_id,
                f"theta={theta}, copies={count}",
            )

    # Every variant participating in one logical experiment case must read the
    # same data/query bytes, even when paths or filenames differ.
    logical_cases: dict[tuple[str, str, str, int], list[tuple[PlanRow, str, str]]] = defaultdict(list)
    for item in usable_rows:
        row = item[0]
        logical_cases[(row.experiment, row.dataset, row.query_id, row.requested_theta)].append(item)
    for (experiment, dataset, query_id, theta), items in logical_cases.items():
        data_hashes = {data_sha for _, data_sha, _ in items}
        query_hashes = {query_sha for _, _, query_sha in items}
        if len(data_hashes) != 1:
            _issue(
                issues,
                "ERROR",
                "CASE_DATA_GRAPH_MISMATCH",
                experiment,
                dataset=dataset,
                query_id=query_id,
                details=f"theta={theta}, hashes={sorted(data_hashes)}",
            )
        if len(query_hashes) != 1:
            _issue(
                issues,
                "ERROR",
                "CASE_QUERY_GRAPH_MISMATCH",
                experiment,
                dataset=dataset,
                query_id=query_id,
                details=f"theta={theta}, hashes={sorted(query_hashes)}",
            )
        preprocessing_values = [preprocessing_hashes[row.source] for row, _, _ in items]
        provided_preprocessing = {value for value in preprocessing_values if value}
        if provided_preprocessing and len(provided_preprocessing) != 1:
            _issue(
                issues,
                "ERROR",
                "CASE_PREPROCESSING_MISMATCH",
                experiment,
                dataset=dataset,
                query_id=query_id,
                details=f"theta={theta}, hashes={sorted(provided_preprocessing)}",
            )
        if provided_preprocessing and any(not value for value in preprocessing_values):
            _issue(
                issues,
                "ERROR",
                "CASE_PREPROCESSING_METADATA_INCOMPLETE",
                experiment,
                dataset=dataset,
                query_id=query_id,
                details=f"theta={theta}",
            )

    query_versions: dict[tuple[str, str, str, str], list[tuple[str, str]]] = defaultdict(list)
    for row, data_sha, query_sha in usable_rows:
        query_versions[(row.experiment, row.variant, row.dataset, row.query_id)].append(
            (data_sha, query_sha)
        )
    for (experiment, variant, dataset, query_id), versions in query_versions.items():
        query_hashes = {query_sha for _, query_sha in versions}
        if len(query_hashes) > 1:
            _issue(
                issues,
                "ERROR",
                "QUERY_VERSION_CHANGES_WITH_THETA",
                experiment,
                variant,
                dataset,
                query_id,
                f"hashes={sorted(query_hashes)}",
            )

    by_experiment_variant_dataset: dict[
        tuple[str, str, str], list[tuple[PlanRow, str, str]]
    ] = defaultdict(list)
    for item in usable_rows:
        row = item[0]
        by_experiment_variant_dataset[(row.experiment, row.variant, row.dataset)].append(item)

    summary_rows: list[dict[str, object]] = []
    for (experiment, variant, dataset), items in sorted(by_experiment_variant_dataset.items()):
        query_thetas: dict[str, set[int]] = defaultdict(set)
        for row, _, _ in items:
            query_thetas[row.query_id].add(row.requested_theta)
        expected_thetas = EXPECTED_THETAS.get(experiment)
        bad_theta_queries = 0
        if expected_thetas is not None:
            for query_id, actual_thetas in query_thetas.items():
                if actual_thetas != expected_thetas:
                    bad_theta_queries += 1
                    _issue(
                        issues,
                        "ERROR",
                        "QUERY_THETA_SET_MISMATCH",
                        experiment,
                        variant,
                        dataset,
                        query_id,
                        f"expected={sorted(expected_thetas)}, actual={sorted(actual_thetas)}",
                    )
        if len(query_thetas) != args.expected_queries:
            _issue(
                issues,
                "ERROR",
                "QUERY_COUNT_MISMATCH",
                experiment,
                variant,
                dataset,
                details=f"expected={args.expected_queries}, actual={len(query_thetas)}",
            )
        summary_rows.append(
            {
                "experiment": experiment,
                "variant": variant,
                "dataset": dataset,
                "row_count": len(items),
                "unique_queries": len(query_thetas),
                "expected_queries": args.expected_queries,
                "query_count_ok": int(len(query_thetas) == args.expected_queries),
                "expected_theta_set": (
                    ",".join(map(str, sorted(expected_thetas))) if expected_thetas else "policy_not_set"
                ),
                "bad_theta_queries": bad_theta_queries,
                "theta_gate": int(expected_thetas is None or bad_theta_queries == 0),
                "unique_data_hashes": len({data_sha for _, data_sha, _ in items}),
                "unique_query_hashes": len({query_sha for _, _, query_sha in items}),
            }
        )

    # All variants of one experiment must have exactly the same dataset/query set.
    experiment_variants: dict[str, dict[str, set[tuple[str, str]]]] = defaultdict(
        lambda: defaultdict(set)
    )
    experiment_dataset_sets: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for row, _, _ in usable_rows:
        experiment_variants[row.experiment][row.variant].add((row.dataset, row.query_id))
        experiment_dataset_sets[row.experiment][row.variant].add(row.dataset)
    for experiment, variants in experiment_variants.items():
        variant_names = sorted(variants)
        if variant_names:
            reference = variants[variant_names[0]]
            for variant in variant_names[1:]:
                if variants[variant] != reference:
                    missing = len(reference - variants[variant])
                    extra = len(variants[variant] - reference)
                    _issue(
                        issues,
                        "ERROR",
                        "VARIANT_QUERY_SET_MISMATCH",
                        experiment,
                        variant,
                        details=(
                            f"reference={variant_names[0]}, missing_query_keys={missing}, "
                            f"extra_query_keys={extra}"
                        ),
                    )
            reference_order = _ordered_unique(
                (row.dataset, row.query_id)
                for row, _, _ in usable_rows
                if row.experiment == experiment and row.variant == variant_names[0]
            )
            for variant in variant_names[1:]:
                variant_order = _ordered_unique(
                    (row.dataset, row.query_id)
                    for row, _, _ in usable_rows
                    if row.experiment == experiment and row.variant == variant
                )
                if variant_order != reference_order:
                    _issue(
                        issues,
                        "ERROR",
                        "VARIANT_QUERY_ORDER_MISMATCH",
                        experiment,
                        variant,
                        details=f"reference={variant_names[0]}",
                    )
        expected_dataset_count = EXPECTED_DATASETS.get(experiment)
        union_datasets = set().union(*experiment_dataset_sets[experiment].values())
        if expected_dataset_count is not None and len(union_datasets) != expected_dataset_count:
            _issue(
                issues,
                "ERROR",
                "EXPERIMENT_DATASET_COUNT_MISMATCH",
                experiment,
                details=f"expected={expected_dataset_count}, actual={len(union_datasets)}",
            )

    # A dataset name must always resolve to the same data graph bytes.
    dataset_hashes: dict[str, set[str]] = defaultdict(set)
    for row, data_sha, _ in usable_rows:
        dataset_hashes[row.dataset].add(data_sha)
    for dataset, hashes in dataset_hashes.items():
        if len(hashes) > 1:
            _issue(
                issues,
                "ERROR",
                "DATASET_VERSION_MISMATCH",
                dataset=dataset,
                details=f"hashes={sorted(hashes)}",
            )

    # Main and ablation must share not merely names, but identical query bytes.
    experiments_present = {row.experiment for row, _, _ in usable_rows}
    if {"main", "ablation"} <= experiments_present:
        query_hash_by_experiment: dict[str, dict[tuple[str, str], set[str]]] = defaultdict(
            lambda: defaultdict(set)
        )
        for row, _, query_sha in usable_rows:
            if row.experiment in {"main", "ablation"}:
                query_hash_by_experiment[row.experiment][(row.dataset, row.query_id)].add(query_sha)
        main_keys = set(query_hash_by_experiment["main"])
        ablation_keys = set(query_hash_by_experiment["ablation"])
        if main_keys != ablation_keys:
            _issue(
                issues,
                "ERROR",
                "MAIN_ABLATION_QUERY_SET_MISMATCH",
                details=(
                    f"missing_from_ablation={len(main_keys - ablation_keys)}, "
                    f"extra_in_ablation={len(ablation_keys - main_keys)}"
                ),
            )
        for dataset_query in main_keys & ablation_keys:
            main_hashes = query_hash_by_experiment["main"][dataset_query]
            ablation_hashes = query_hash_by_experiment["ablation"][dataset_query]
            if main_hashes != ablation_hashes or len(main_hashes) != 1:
                dataset, query_id = dataset_query
                _issue(
                    issues,
                    "ERROR",
                    "MAIN_ABLATION_QUERY_VERSION_MISMATCH",
                    dataset=dataset,
                    query_id=query_id,
                    details=f"main={sorted(main_hashes)}, ablation={sorted(ablation_hashes)}",
                )
        main_order = _ordered_unique(
            (row.dataset, row.query_id)
            for row, _, _ in usable_rows
            if row.experiment == "main"
        )
        ablation_order = _ordered_unique(
            (row.dataset, row.query_id)
            for row, _, _ in usable_rows
            if row.experiment == "ablation"
        )
        if main_order != ablation_order:
            _issue(
                issues,
                "ERROR",
                "MAIN_ABLATION_QUERY_ORDER_MISMATCH",
                details="ordered dataset/query sequence differs",
            )

    output_dir = args.output_dir.resolve()
    _write_tsv(
        output_dir / "experiment_plan_resolved.tsv",
        [
            "experiment",
            "variant",
            "dataset",
            "query_id",
            "requested_theta",
            "data_graph",
            "data_sha256",
            "query_graph",
            "query_sha256",
            "preprocessing_artifact",
            "preprocessing_sha256",
            "source",
        ],
        resolved_rows,
    )
    _write_tsv(
        output_dir / "experiment_plan_summary.tsv",
        [
            "experiment",
            "variant",
            "dataset",
            "row_count",
            "unique_queries",
            "expected_queries",
            "query_count_ok",
            "expected_theta_set",
            "bad_theta_queries",
            "theta_gate",
            "unique_data_hashes",
            "unique_query_hashes",
        ],
        summary_rows,
    )
    _write_tsv(
        output_dir / "experiment_plan_issues.tsv",
        ["severity", "code", "experiment", "variant", "dataset", "query_id", "details"],
        issues,
    )

    error_count = sum(issue["severity"] == "ERROR" for issue in issues)
    warning_count = sum(issue["severity"] == "WARNING" for issue in issues)
    print(
        f"resolved_rows={len(resolved_rows)} errors={error_count} warnings={warning_count} "
        f"output={output_dir}",
        file=sys.stderr,
    )
    if args.strict and error_count:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
