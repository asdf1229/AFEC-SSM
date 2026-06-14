#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import dataclasses
import datetime as _datetime
import os
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


@dataclasses.dataclass(frozen=True)
class Variant:
    name: str
    macros: Tuple[str, ...] = ()
    env: Tuple[Tuple[str, str], ...] = ()
    description: str = ""


@dataclasses.dataclass(frozen=True)
class Case:
    case_id: str
    dataset_group: str
    dataset: str
    query_group: str
    query: str
    threshold: int
    data_graph: Path
    query_graph: Path


@dataclasses.dataclass
class RunResult:
    case_id: str
    dataset_group: str
    dataset: str
    query_group: str
    query: str
    threshold: int
    variant: str
    repeat: int
    status: str
    return_code: str
    count: str = "NA"
    load_ms: str = "NA"
    run_ms: str = "NA"
    total_ms: str = "NA"
    recursion_calls: str = "NA"
    pruning_calls: str = "NA"
    terminal_tail_calls: str = "NA"
    terminal_prunes: str = "NA"
    terminal_delayed: str = "NA"
    terminal_cand_checks: str = "NA"
    elapsed_wall_ms: str = "NA"
    output: str = "NA"


DEFAULT_VARIANTS: Tuple[Variant, ...] = (
    Variant("base", description="current default CDE-Edge-IE"),
    Variant("no_spoke", macros=("DISABLE_SPOKE_FILTERING",),
        description="disable spoke candidate filtering"),
    Variant("onehop", macros=("ENABLE_ONEHOP_FILTERING",),
        description="enable one-hop filtering on top of spoke records"),
    Variant("no_terminal_buckets", env=(("CDE_EDGE_IE_TERMINAL_BUCKETS", "0"),),
        description="disable terminal-tail bucket ordering at runtime"),
)

DEFAULT_VARIANT_CONFIG = "scripts/cde_ablation_variants.tsv"


PRESETS = {
    "smoke": {
        "thresholds": "0,1",
        "query_limit": 2,
        "max_cases": 4,
        "timeout": 20.0,
    },
    "pilot": {
        "thresholds": "0,1,2",
        "query_limit": 5,
        "max_cases": 60,
        "timeout": 60.0,
    },
    "full": {
        "thresholds": "0,1,2,3,4,5,6",
        "query_limit": None,
        "max_cases": None,
        "timeout": 300.0,
    },
}


RUN_COLUMNS = (
    "case_id",
    "dataset_group",
    "dataset",
    "query_group",
    "query",
    "threshold",
    "variant",
    "repeat",
    "status",
    "return_code",
    "count",
    "load_ms",
    "run_ms",
    "total_ms",
    "recursion_calls",
    "pruning_calls",
    "terminal_tail_calls",
    "terminal_prunes",
    "terminal_delayed",
    "terminal_cand_checks",
    "elapsed_wall_ms",
    "output",
)


STAT_LABELS = {
    "recursion_calls": "Recursion Calls",
    "pruning_calls": "Pruning Calls",
    "terminal_tail_calls": "Terminal Tail Calls",
    "terminal_prunes": "Terminal Prunes",
    "terminal_delayed": "Terminal Delayed",
    "terminal_cand_checks": "Terminal Cand Checks",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_thresholds(value: str) -> List[int]:
    thresholds: List[int] = []
    for raw in value.split(","):
        raw = raw.strip()
        if not raw:
            continue
        thresholds.append(int(raw))
    if not thresholds:
        raise ValueError("threshold list is empty")
    return thresholds


def tsv(value: object) -> str:
    text = str(value)
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(tsv(v) for v in row) + "\n")


def query_sort_key(path: Path) -> Tuple[int, int, str]:
    name = path.stem
    if re.fullmatch(r"[0-9]+", name):
        return (2_147_483_647, int(name), name)
    match = re.fullmatch(r"query_([0-9]+)_([0-9]+)", name)
    if match:
        return (int(match.group(1)), int(match.group(2)), name)
    return (2_147_483_647, 2_147_483_647, name)


def query_group_sort_key(entry: Tuple[str, Path]) -> Tuple[int, int, str]:
    name = entry[0]
    rank = 100
    value = 2_147_483_647
    if name in (".", "baseline", "legacy"):
        rank = 0
        value = 0
    else:
        for prefix, next_rank in (
            ("vertices_num_", 10),
            ("avg_degree_", 20),
            ("missing_edge_threshold_", 30),
        ):
            if name.startswith(prefix):
                suffix = name[len(prefix):]
                if suffix.isdigit():
                    rank = next_rank
                    value = int(suffix)
                break
    return (rank, value, name)


def discover_dataset_dirs(data_dir: Path) -> List[Tuple[str, Path]]:
    roots: List[Tuple[str, Path]] = []
    if (data_dir / "synthetic").is_dir():
        roots.append(("synthetic", data_dir / "synthetic"))
    if (data_dir / "real_graphs").is_dir():
        roots.append(("real", data_dir / "real_graphs"))
    if not roots:
        roots.append(("custom", data_dir))

    datasets: List[Tuple[str, Path]] = []
    for group, root in roots:
        for graph in root.rglob("graph_g.txt"):
            if "query_graph" in graph.parts:
                continue
            datasets.append((group, graph.parent))
    return sorted(datasets, key=lambda item: (item[0], str(item[1])))


def discover_query_groups(query_root: Path) -> List[Tuple[str, Path]]:
    groups: List[Tuple[str, Path]] = []
    if any(query_root.glob("*.txt")):
        groups.append(("legacy", query_root))
    for child in sorted(query_root.iterdir() if query_root.is_dir() else []):
        if child.is_dir() and any(child.glob("*.txt")):
            groups.append((child.name, child))
    return sorted(groups, key=query_group_sort_key)


def thresholds_for_group(group: str, default_thresholds: Sequence[int], explicit: bool) -> List[int]:
    if not explicit:
        match = re.fullmatch(r"missing_edge_threshold_([0-9]+)", group)
        if match:
            return [int(match.group(1))]
    return list(default_thresholds)


def compile_regex(value: Optional[str]) -> Optional[re.Pattern[str]]:
    return re.compile(value) if value else None


def discover_dataset_cases(
    data_dir: Path,
    thresholds: Sequence[int],
    thresholds_explicit: bool,
    query_limit: Optional[int],
    max_cases: Optional[int],
    dataset_regex: Optional[re.Pattern[str]],
    query_group_regex: Optional[re.Pattern[str]],
) -> List[Case]:
    cases: List[Case] = []
    case_index = 0

    for dataset_group, dataset_dir in discover_dataset_dirs(data_dir):
        dataset = dataset_dir.name
        dataset_key = f"{dataset_group}/{dataset}"
        if dataset_regex and not dataset_regex.search(dataset_key) and not dataset_regex.search(str(dataset_dir)):
            continue

        graph = dataset_dir / "graph_g.txt"
        query_root = dataset_dir / "query_graph"
        if not graph.is_file() or not query_root.is_dir():
            continue

        for query_group, query_group_dir in discover_query_groups(query_root):
            if query_group_regex and not query_group_regex.search(query_group):
                continue
            query_files = sorted(query_group_dir.glob("*.txt"), key=query_sort_key)
            if query_limit is not None:
                query_files = query_files[:query_limit]

            for query_file in query_files:
                for threshold in thresholds_for_group(query_group, thresholds, thresholds_explicit):
                    case_index += 1
                    cases.append(Case(
                        case_id=f"case_{case_index:06d}",
                        dataset_group=dataset_group,
                        dataset=dataset,
                        query_group=query_group,
                        query=query_file.stem,
                        threshold=int(threshold),
                        data_graph=graph,
                        query_graph=query_file,
                    ))
                    if max_cases is not None and len(cases) >= max_cases:
                        return cases
    return cases


def find_case_file(case_dir: Path, stem: str) -> Optional[Path]:
    for candidate in (case_dir / stem, case_dir / f"{stem}.txt"):
        if candidate.is_file():
            return candidate
    return None


def extract_threshold(source: Path) -> Optional[int]:
    text = source.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"(^|[^A-Za-z0-9_])t\s*=\s*([0-9]+)", text)
    return int(match.group(2)) if match else None


def discover_tmp_cases(cases_dir: Path, max_cases: Optional[int]) -> List[Case]:
    cases: List[Case] = []
    for child in sorted(cases_dir.iterdir() if cases_dir.is_dir() else []):
        if not child.is_dir():
            continue
        graph = find_case_file(child, "graph_g")
        query = find_case_file(child, "q")
        source = find_case_file(child, "source")
        if graph is None or query is None or source is None:
            continue
        threshold = extract_threshold(source)
        if threshold is None:
            continue
        case_id = f"case_{len(cases) + 1:06d}"
        cases.append(Case(
            case_id=case_id,
            dataset_group="custom",
            dataset=child.name,
            query_group="tmp",
            query=child.name,
            threshold=threshold,
            data_graph=graph,
            query_graph=query,
        ))
        if max_cases is not None and len(cases) >= max_cases:
            break
    return cases


def parse_list_field(value: str) -> Tuple[str, ...]:
    value = value.strip()
    if not value or value == "-":
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def parse_env_field(value: str, source: str, line_number: int) -> Tuple[Tuple[str, str], ...]:
    entries: List[Tuple[str, str]] = []
    for item in parse_list_field(value):
        if "=" not in item:
            raise SystemExit(f"{source}:{line_number}: env entry must be KEY=VALUE: {item}")
        key, env_value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise SystemExit(f"{source}:{line_number}: env key is empty")
        entries.append((key, env_value.strip()))
    return tuple(entries)


def load_variant_config(path: Path) -> List[Variant]:
    variants: List[Variant] = []
    seen = set()

    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(
            (line for line in fh if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        if reader.fieldnames is None:
            raise SystemExit(f"variant config is empty: {path}")

        required = {"name", "macros", "env", "description"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise SystemExit(f"{path}: missing columns: {', '.join(sorted(missing))}")

        for line_number, row in enumerate(reader, start=2):
            name = (row.get("name") or "").strip()
            if not name:
                raise SystemExit(f"{path}:{line_number}: variant name is empty")
            if not re.fullmatch(r"[A-Za-z0-9_.=-]+", name):
                raise SystemExit(
                    f"{path}:{line_number}: variant name must use letters, numbers, _, ., =, or -: {name}"
                )
            if name in seen:
                raise SystemExit(f"{path}:{line_number}: duplicate variant name: {name}")
            seen.add(name)

            variants.append(Variant(
                name=name,
                macros=parse_list_field(row.get("macros") or ""),
                env=parse_env_field(row.get("env") or "", str(path), line_number),
                description=(row.get("description") or "").strip(),
            ))

    if not variants:
        raise SystemExit(f"variant config has no variants: {path}")
    return variants


def resolve_variant_config(repo: Path, configured_path: Optional[str]) -> Tuple[List[Variant], str]:
    raw_path = configured_path or DEFAULT_VARIANT_CONFIG
    path = Path(raw_path)
    if not path.is_absolute():
        path = repo / path

    if path.is_file():
        return load_variant_config(path), str(path)
    if configured_path:
        raise SystemExit(f"variant config not found: {path}")
    return list(DEFAULT_VARIANTS), "built-in defaults"


def resolve_variants(all_variants: Sequence[Variant], names: Optional[str]) -> List[Variant]:
    by_name = {variant.name: variant for variant in all_variants}
    if not names:
        return list(all_variants)
    selected: List[Variant] = []
    for raw in names.split(","):
        name = raw.strip()
        if not name:
            continue
        if name not in by_name:
            raise SystemExit(f"unknown variant '{name}'. Use --list-variants.")
        selected.append(by_name[name])
    if not selected:
        raise SystemExit("no variants selected")
    return selected


def run_checked(command: Sequence[str], cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=str(cwd), check=True)


def build_variant(repo: Path, build_root: Path, variant: Variant, build_parallel: int, skip_build: bool) -> Path:
    build_dir = build_root / variant.name
    exe = build_dir / "ssm_ged_cde_edge_ie"
    if skip_build:
        if not exe.is_file():
            raise SystemExit(f"missing executable for --skip-build: {exe}")
        return exe

    build_dir.mkdir(parents=True, exist_ok=True)
    definitions = ";".join(variant.macros)
    run_checked([
        "cmake",
        "-S",
        str(repo),
        "-B",
        str(build_dir),
        f"-DSSM_GED_COMPILE_DEFINITIONS={definitions}",
    ], cwd=repo)
    run_checked([
        "cmake",
        "--build",
        str(build_dir),
        "--target",
        "ssm_ged_cde_edge_ie",
        "--parallel",
        str(build_parallel),
    ], cwd=repo)
    if not exe.is_file():
        raise SystemExit(f"build did not produce expected executable: {exe}")
    return exe


def parse_output(path: Path) -> Dict[str, str]:
    metrics: Dict[str, str] = {}
    text = path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""

    summary_line = ""
    for line in text.splitlines():
        if line.startswith("SSM_GED_SUMMARY "):
            summary_line = line
    for token in summary_line.split():
        if "=" in token:
            key, value = token.split("=", 1)
            if key in {"count", "load_ms", "run_ms", "total_ms"}:
                metrics[key] = value

    for field, label in STAT_LABELS.items():
        pattern = re.compile(rf"^\s*{re.escape(label)}:\s*([0-9]+)", re.MULTILINE)
        match = pattern.search(text)
        if match:
            metrics[field] = match.group(1)

    return metrics


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.=-]+", "_", value)


def run_one(
    repo: Path,
    result_dir: Path,
    case: Case,
    variant: Variant,
    executable: Path,
    repeat: int,
    timeout: float,
) -> RunResult:
    output_dir = (
        result_dir / "outputs" / safe_name(case.dataset_group) / safe_name(case.dataset) /
        safe_name(case.query_group) / safe_name(variant.name)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{safe_name(case.query)}_t={case.threshold}_r={repeat}.txt"

    command = [
        str(executable),
        "-d",
        str(case.data_graph),
        "-q",
        str(case.query_graph),
        "-t",
        str(case.threshold),
    ]
    env = os.environ.copy()
    env.update(dict(variant.env))
    start = time.monotonic()

    try:
        with output.open("w", encoding="utf-8") as fh:
            completed = subprocess.run(
                command,
                cwd=str(repo),
                env=env,
                stdout=fh,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
        return_code = str(completed.returncode)
        status = "OK" if completed.returncode == 0 else "RunError"
    except subprocess.TimeoutExpired:
        with output.open("a", encoding="utf-8") as fh:
            fh.write(f"\nSSM_GED_ABLATION_TIMEOUT seconds={timeout}\n")
        return_code = "Timeout"
        status = "Timeout"

    elapsed_ms = (time.monotonic() - start) * 1000.0
    parsed = parse_output(output) if status == "OK" else {}
    if status == "OK" and "run_ms" not in parsed:
        status = "ParseError"

    return RunResult(
        case_id=case.case_id,
        dataset_group=case.dataset_group,
        dataset=case.dataset,
        query_group=case.query_group,
        query=case.query,
        threshold=case.threshold,
        variant=variant.name,
        repeat=repeat,
        status=status,
        return_code=return_code,
        count=parsed.get("count", "NA"),
        load_ms=parsed.get("load_ms", "NA"),
        run_ms=parsed.get("run_ms", "NA"),
        total_ms=parsed.get("total_ms", "NA"),
        recursion_calls=parsed.get("recursion_calls", "NA"),
        pruning_calls=parsed.get("pruning_calls", "NA"),
        terminal_tail_calls=parsed.get("terminal_tail_calls", "NA"),
        terminal_prunes=parsed.get("terminal_prunes", "NA"),
        terminal_delayed=parsed.get("terminal_delayed", "NA"),
        terminal_cand_checks=parsed.get("terminal_cand_checks", "NA"),
        elapsed_wall_ms=f"{elapsed_ms:.3f}",
        output=str(output.resolve()),
    )


def median_string(values: Sequence[str]) -> str:
    numeric: List[float] = []
    for value in values:
        try:
            numeric.append(float(value))
        except ValueError:
            pass
    if not numeric:
        return "NA"
    med = statistics.median(numeric)
    if med.is_integer():
        return str(int(med))
    return f"{med:.4f}"


def write_cases(path: Path, cases: Sequence[Case]) -> None:
    rows = (
        (
            case.case_id,
            case.dataset_group,
            case.dataset,
            case.query_group,
            case.query,
            case.threshold,
            case.data_graph.resolve(),
            case.query_graph.resolve(),
        )
        for case in cases
    )
    write_tsv(path, (
        "case_id",
        "dataset_group",
        "dataset",
        "query_group",
        "query",
        "threshold",
        "data_graph",
        "query_graph",
    ), rows)


def write_variants(path: Path, variants: Sequence[Variant], executables: Dict[str, Path]) -> None:
    rows = (
        (
            variant.name,
            ",".join(variant.macros),
            ",".join(f"{k}={v}" for k, v in variant.env),
            variant.description,
            executables[variant.name].resolve(),
        )
        for variant in variants
    )
    write_tsv(path, ("variant", "macros", "env", "description", "executable"), rows)


def write_runs(path: Path, results: Sequence[RunResult]) -> None:
    rows = (
        tuple(getattr(result, column) for column in RUN_COLUMNS)
        for result in results
    )
    write_tsv(path, RUN_COLUMNS, rows)


def write_summary(
    path: Path,
    failures_path: Path,
    cases: Sequence[Case],
    variants: Sequence[Variant],
    results: Sequence[RunResult],
    repeat_count: int,
) -> int:
    by_case_variant: Dict[Tuple[str, str], List[RunResult]] = {}
    for result in results:
        by_case_variant.setdefault((result.case_id, result.variant), []).append(result)

    columns = ["dataset_group", "dataset", "query_group", "query", "threshold", "status", "expected_count"]
    for variant in variants:
        columns.extend([
            f"{variant.name}_count",
            f"{variant.name}_run_ms",
            f"{variant.name}_recursion_calls",
            f"{variant.name}_output",
        ])

    rows: List[List[object]] = []
    failures: List[List[object]] = []

    for case in cases:
        row: List[object] = [
            case.dataset_group,
            case.dataset,
            case.query_group,
            case.query,
            case.threshold,
        ]
        case_status = "OK"
        expected_count = "NA"
        counts_by_variant: Dict[str, str] = {}
        variant_cells: Dict[str, Tuple[str, str, str, str]] = {}

        for variant in variants:
            run_group = sorted(
                by_case_variant.get((case.case_id, variant.name), []),
                key=lambda item: item.repeat,
            )
            ok_group = [item for item in run_group if item.status == "OK"]
            if len(run_group) != repeat_count:
                case_status = "FAIL"
                failures.append([case.case_id, variant.name, "MissingRepeat", "NA"])
            for item in run_group:
                if item.status != "OK":
                    case_status = "FAIL"
                    failures.append([case.case_id, variant.name, item.status, item.output])

            count_values = {item.count for item in ok_group if item.count != "NA"}
            if len(count_values) > 1:
                case_status = "COUNT_MISMATCH"
                failures.append([case.case_id, variant.name, "RepeatCountMismatch", ",".join(sorted(count_values))])

            count = next(iter(count_values), "NA")
            if count != "NA":
                counts_by_variant[variant.name] = count
            run_ms = median_string([item.run_ms for item in ok_group])
            recursion = median_string([item.recursion_calls for item in ok_group])
            output = ok_group[0].output if ok_group else (run_group[0].output if run_group else "NA")
            variant_cells[variant.name] = (count, run_ms, recursion, output)

        distinct_counts = {count for count in counts_by_variant.values() if count != "NA"}
        if len(distinct_counts) > 1:
            case_status = "COUNT_MISMATCH"
            failures.append([case.case_id, "ALL", "VariantCountMismatch", ",".join(sorted(distinct_counts))])

        first_variant = variants[0].name if variants else ""
        expected_count = counts_by_variant.get(first_variant, next(iter(distinct_counts), "NA"))

        row.extend([case_status, expected_count])
        for variant in variants:
            row.extend(variant_cells[variant.name])
        rows.append(row)

    write_tsv(path, columns, rows)
    write_tsv(failures_path, ("case_id", "variant", "status", "detail"), failures)
    return len(failures)


def list_variants(variants: Sequence[Variant], source: str) -> None:
    print(f"Variant config: {source}")
    for variant in variants:
        macros = ",".join(variant.macros) or "-"
        env = ",".join(f"{k}={v}" for k, v in variant.env) or "-"
        print(f"{variant.name}\tmacros={macros}\tenv={env}\t{variant.description}")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and run CDE-Edge-IE ablation variants.")
    parser.add_argument("--preset", choices=sorted(PRESETS), default="smoke",
        help="case selection preset; default: smoke")
    parser.add_argument("-d", "--data-dir", default="test/datasets",
        help="dataset root with synthetic/ and real_graphs/; default: test/datasets")
    parser.add_argument("--cases-dir",
        help="run tmp-style case folders containing graph_g, q, and source with t=x")
    parser.add_argument("--result-dir",
        help="output directory; default: result/cde_ablation_<timestamp>")
    parser.add_argument("--build-root", default="build/cde_ablations",
        help="root for per-variant build directories; default: build/cde_ablations")
    parser.add_argument("--skip-build", action="store_true",
        help="reuse existing per-variant executables under --build-root")
    parser.add_argument("--build-parallel", type=int, default=max(1, os.cpu_count() or 1),
        help="parallel jobs passed to cmake --build")
    parser.add_argument("-p", "--parallel", type=int, default=1,
        help="parallel algorithm runs; default: 1")
    parser.add_argument("-t", "--timeout", type=float,
        help="per-run timeout in seconds; default comes from --preset")
    parser.add_argument("--repeat", type=int, default=1,
        help="repeat each case/variant and summarize medians; default: 1")
    parser.add_argument("--variants",
        help="comma-separated variant names; default: all variants from --variant-config")
    parser.add_argument("--variant-config",
        help=f"TSV variant config; default: {DEFAULT_VARIANT_CONFIG} when present")
    parser.add_argument("--list-variants", action="store_true",
        help="print variants from the selected config and exit")
    parser.add_argument("--dataset-regex",
        help="only include datasets whose group/name or path matches this regex")
    parser.add_argument("--query-group-regex",
        help="only include query groups matching this regex")
    parser.add_argument("--query-limit", type=int,
        help="maximum queries per query group; default comes from --preset")
    parser.add_argument("--max-cases", type=int,
        help="maximum expanded cases before variant/repeat multiplication; default comes from --preset")
    parser.add_argument("--thresholds",
        help="comma-separated thresholds; default comes from --preset")
    parser.add_argument("--allow-failures", action="store_true",
        help="exit 0 even when a run fails or counts mismatch")
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    repo = repo_root()
    all_variants, variant_source = resolve_variant_config(repo, args.variant_config)
    if args.list_variants:
        list_variants(all_variants, variant_source)
        return 0

    if args.parallel < 1:
        raise SystemExit("--parallel must be positive")
    if args.repeat < 1:
        raise SystemExit("--repeat must be positive")

    preset = PRESETS[args.preset]
    thresholds_raw = args.thresholds if args.thresholds is not None else preset["thresholds"]
    thresholds = parse_thresholds(str(thresholds_raw))
    thresholds_explicit = args.thresholds is not None
    query_limit = args.query_limit if args.query_limit is not None else preset["query_limit"]
    max_cases = args.max_cases if args.max_cases is not None else preset["max_cases"]
    timeout = args.timeout if args.timeout is not None else float(preset["timeout"])

    timestamp = _datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = Path(args.result_dir) if args.result_dir else repo / "result" / f"cde_ablation_{timestamp}"
    if not result_dir.is_absolute():
        result_dir = repo / result_dir
    build_root = Path(args.build_root)
    if not build_root.is_absolute():
        build_root = repo / build_root

    variants = resolve_variants(all_variants, args.variants)
    dataset_regex = compile_regex(args.dataset_regex)
    query_group_regex = compile_regex(args.query_group_regex)

    if args.cases_dir:
        cases_dir = Path(args.cases_dir)
        if not cases_dir.is_absolute():
            cases_dir = repo / cases_dir
        cases = discover_tmp_cases(cases_dir, max_cases)
    else:
        data_dir = Path(args.data_dir)
        if not data_dir.is_absolute():
            data_dir = repo / data_dir
        cases = discover_dataset_cases(
            data_dir=data_dir,
            thresholds=thresholds,
            thresholds_explicit=thresholds_explicit,
            query_limit=query_limit,
            max_cases=max_cases,
            dataset_regex=dataset_regex,
            query_group_regex=query_group_regex,
        )

    if not cases:
        raise SystemExit("no cases discovered; adjust --data-dir/--cases-dir or filters")

    result_dir.mkdir(parents=True, exist_ok=True)
    print(f"Results: {result_dir}")
    print(f"Variant config: {variant_source}")
    print(f"Preset: {args.preset}; cases={len(cases)} variants={len(variants)} repeat={args.repeat}")
    print(f"Per-run timeout: {timeout}s; run parallel={args.parallel}")

    executables: Dict[str, Path] = {}
    for variant in variants:
        print(f"\n=== Build variant: {variant.name} ===", flush=True)
        executables[variant.name] = build_variant(repo, build_root, variant, args.build_parallel, args.skip_build)

    write_cases(result_dir / "cases.tsv", cases)
    write_variants(result_dir / "variants.tsv", variants, executables)

    jobs = []
    for case in cases:
        for variant in variants:
            for repeat in range(1, args.repeat + 1):
                jobs.append((case, variant, repeat))

    print(f"\nQueued runs: {len(jobs)}")
    results: List[RunResult] = []
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as pool:
        future_to_job = {
            pool.submit(
                run_one,
                repo,
                result_dir,
                case,
                variant,
                executables[variant.name],
                repeat,
                timeout,
            ): (case, variant, repeat)
            for case, variant, repeat in jobs
        }
        for future in concurrent.futures.as_completed(future_to_job):
            case, variant, repeat = future_to_job[future]
            completed += 1
            try:
                result = future.result()
            except Exception as exc:
                output = result_dir / "outputs" / "internal_errors.tsv"
                result = RunResult(
                    case_id=case.case_id,
                    dataset_group=case.dataset_group,
                    dataset=case.dataset,
                    query_group=case.query_group,
                    query=case.query,
                    threshold=case.threshold,
                    variant=variant.name,
                    repeat=repeat,
                    status="InternalError",
                    return_code="InternalError",
                    output=str(output.resolve()),
                )
                with output.open("a", encoding="utf-8") as fh:
                    fh.write(f"{case.case_id}\t{variant.name}\tr{repeat}\t{exc!r}\n")
            results.append(result)
            print(
                f"[{completed}/{len(jobs)}] {result.status} "
                f"{case.dataset_group}/{case.dataset} {case.query_group}/{case.query} "
                f"t={case.threshold} {variant.name} r{repeat} "
                f"run_ms={result.run_ms} count={result.count}",
                flush=True,
            )

    results.sort(key=lambda item: (
        item.case_id,
        next(i for i, variant in enumerate(variants) if variant.name == item.variant),
        item.repeat,
    ))

    runs_path = result_dir / "runs.tsv"
    summary_path = result_dir / "summary.tsv"
    failures_path = result_dir / "failures.tsv"
    write_runs(runs_path, results)
    failures = write_summary(summary_path, failures_path, cases, variants, results, args.repeat)

    print("\n===== CDE ABLATION DONE =====")
    print(f"Result dir:  {result_dir.resolve()}")
    print(f"Runs TSV:    {runs_path.resolve()}")
    print(f"Summary TSV: {summary_path.resolve()}")
    print(f"Failures:    {failures}")
    if failures and not args.allow_failures:
        print(f"Failure TSV: {failures_path.resolve()}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
