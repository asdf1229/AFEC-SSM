from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import audit_datasets
import audit_experiment_plan
import audit_result_quality
import repair_queries
from graph_audit_common import (
    are_label_isomorphic,
    connected_components,
    parse_query_graph,
    scan_data_graph,
)


def write_graph(path: Path, graph_id: str, labels: list[str], edges: list[tuple[int, int]]) -> None:
    lines = [f"t # {graph_id}"]
    lines.extend(f"v {vertex} {label}" for vertex, label in enumerate(labels))
    lines.extend(f"e {u} {v}" for u, v in edges)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class GraphAuditCommonTests(unittest.TestCase):
    def test_query_connectivity_slack_inputs_and_isomorphism(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.txt"
            permuted = root / "permuted.txt"
            path = root / "path.txt"
            write_graph(first, "a", ["A", "A", "B"], [(0, 1), (0, 2), (1, 2)])
            write_graph(permuted, "b", ["B", "A", "A"], [(0, 1), (0, 2), (1, 2)])
            write_graph(path, "c", ["A", "A", "B"], [(0, 1), (1, 2)])

            first_graph = parse_query_graph(first)
            permuted_graph = parse_query_graph(permuted)
            path_graph = parse_query_graph(path)
            self.assertEqual(len(connected_components(first_graph)), 1)
            self.assertEqual(first_graph.edge_count - first_graph.vertex_count + 1, 1)
            self.assertTrue(are_label_isomorphic(first_graph, permuted_graph))
            self.assertFalse(are_label_isomorphic(first_graph, path_graph))

    def test_streamed_data_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            graph_path = Path(directory) / "graph_g.txt"
            graph_path.write_text(
                "\n".join(
                    [
                        "t # data",
                        "v 0 A",
                        "v 1 B",
                        "v 2 B",
                        "v 3 C",
                        "e 0 1",
                        "e 0 1",
                        "e 1 2",
                        "e 3 3",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            stats = scan_data_graph(graph_path)
            self.assertEqual(stats.vertex_count, 4)
            self.assertEqual(stats.edge_count, 2)
            self.assertEqual(stats.duplicate_edge_records, 1)
            self.assertEqual(stats.self_loops, 1)
            self.assertEqual(stats.component_count, 2)
            self.assertEqual(stats.largest_component_vertices, 3)
            self.assertEqual(stats.isolated_vertices, 1)

    def test_repair_sampler_builds_an_induced_query_with_target_slack(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            graph_path = Path(directory) / "graph_g.txt"
            edges = [
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 7),
                (0, 2),
                (2, 4),
                (4, 6),
                (0, 7),
            ]
            write_graph(graph_path, "data", ["A"] * 8, sorted(edges))
            data = repair_queries.DataGraphCSR.load(graph_path)
            rng = __import__("random").Random(7)
            vertices = None
            for _ in range(100):
                vertices = repair_queries._sample_vertices(data, 8, 4, rng, 64)
                if vertices is not None:
                    break
            self.assertIsNotNone(vertices)
            candidate = repair_queries._build_candidate(
                data, "q", Path(directory) / "q.txt", vertices
            )
            self.assertEqual(candidate.slack, 4)
            self.assertEqual(len(candidate.edges), 11)
            self.assertIsNone(repair_queries._query_group_violation("sparse_8", candidate.graph))
            self.assertEqual(
                repair_queries._query_group_violation("dense_8", candidate.graph),
                "query_group_density_mismatch",
            )
            self.assertEqual(
                repair_queries._repair_target(
                    "dense_8", candidate.graph, "query_group_density_mismatch", 4
                ),
                (8, 6),
            )


class AuditToolTests(unittest.TestCase):
    def test_static_audit_writes_fixed_inventories(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = root / "datasets" / "real_graphs" / "Tiny"
            write_graph(
                dataset / "graph_g.txt",
                "data",
                ["A", "A", "A", "A"],
                [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            )
            write_graph(
                dataset / "query_graph" / "dense_4" / "1.txt",
                "q1",
                ["A", "A", "A", "A"],
                [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            )
            write_graph(
                dataset / "query_graph" / "sparse_4" / "2.txt",
                "q2",
                ["A", "A", "A", "A"],
                [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            )
            output = root / "audit"
            result = audit_datasets.main(
                [
                    "--data-root",
                    str(root / "datasets"),
                    "--output-dir",
                    str(output),
                    "--expected-queries",
                    "2",
                ]
            )
            self.assertEqual(result, 0)
            for name in (
                "dataset_inventory.tsv",
                "query_inventory.tsv",
                "data_audit_summary.tsv",
                "query_group_audit_summary.tsv",
                "sha256_manifest.tsv",
            ):
                self.assertTrue((output / name).is_file())
            with (output / "data_audit_summary.tsv").open(newline="") as source:
                row = next(csv.DictReader(source, delimiter="\t"))
            self.assertEqual(row["query_count"], "2")
            self.assertEqual(row["unique_file_contents"], "2")
            self.assertEqual(row["unique_isomorphism_classes"], "1")
            self.assertEqual(row["slack_ge_2"], "2")
            self.assertEqual(row["slack_ge_4"], "0")
            self.assertEqual(row["query_group_design_queries"], "2")
            self.assertEqual(row["query_group_size_violations"], "0")
            self.assertEqual(row["query_group_density_violations"], "1")
            self.assertEqual(row["gate_query_group_design"], "0")
            with (output / "query_group_audit_summary.tsv").open(newline="") as source:
                group_rows = list(csv.DictReader(source, delimiter="\t"))
            self.assertEqual({row["query_group"] for row in group_rows}, {"dense_4", "sparse_4"})
            violations = {row["query_group"]: row["density_violations"] for row in group_rows}
            self.assertEqual(violations, {"dense_4": "1", "sparse_4": "0"})

    def test_plan_audit_detects_different_query_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "Tiny" / "graph_g.txt"
            first = root / "Tiny" / "query_graph" / "1.txt"
            second = root / "Tiny" / "query_graph" / "2.txt"
            write_graph(data, "d", ["A", "A"], [(0, 1)])
            write_graph(first, "q1", ["A", "A"], [(0, 1)])
            write_graph(second, "q2", ["A", "A"], [])
            plan = root / "plan.tsv"
            plan.write_text(
                "experiment\tvariant\tdataset\tquery_id\tdata_graph\tquery_graph\trequested_theta\n"
                f"main\ta\tTiny\tq\t{data}\t{first}\t2\n"
                f"main\tb\tTiny\tq\t{data}\t{second}\t2\n",
                encoding="utf-8",
            )
            output = root / "plan_audit"
            result = audit_experiment_plan.main(
                [str(plan), "--expected-queries", "1", "--output-dir", str(output)]
            )
            self.assertEqual(result, 0)
            with (output / "experiment_plan_issues.tsv").open(newline="") as source:
                codes = {row["code"] for row in csv.DictReader(source, delimiter="\t")}
            self.assertIn("CASE_QUERY_GRAPH_MISMATCH", codes)

    def test_plan_audit_detects_order_and_theta_version_changes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "Tiny" / "graph_g.txt"
            first = root / "Tiny" / "query_graph" / "1.txt"
            second = root / "Tiny" / "query_graph" / "2.txt"
            write_graph(data, "d", ["A", "A"], [(0, 1)])
            write_graph(first, "q1", ["A", "A"], [(0, 1)])
            write_graph(second, "q2", ["A", "A"], [])
            plan = root / "plan.tsv"
            plan.write_text(
                "experiment\tvariant\tdataset\tquery_id\tdata_graph\tquery_graph\trequested_theta\n"
                f"main\ta\tTiny\tq1\t{data}\t{first}\t2\n"
                f"main\ta\tTiny\tq2\t{data}\t{second}\t2\n"
                f"main\tb\tTiny\tq2\t{data}\t{second}\t2\n"
                f"main\tb\tTiny\tq1\t{data}\t{first}\t2\n"
                f"threshold\tt\tTiny\tq\t{data}\t{first}\t1\n"
                f"threshold\tt\tTiny\tq\t{data}\t{second}\t2\n",
                encoding="utf-8",
            )
            output = root / "plan_audit"
            result = audit_experiment_plan.main(
                [str(plan), "--expected-queries", "2", "--output-dir", str(output)]
            )
            self.assertEqual(result, 0)
            with (output / "experiment_plan_issues.tsv").open(newline="") as source:
                codes = {row["code"] for row in csv.DictReader(source, delimiter="\t")}
            self.assertIn("VARIANT_QUERY_ORDER_MISMATCH", codes)
            self.assertIn("QUERY_VERSION_CHANGES_WITH_THETA", codes)

    def test_result_quality_warnings(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = root / "summary.tsv"
            summary.write_text(
                "dataset_group\tdataset\tquery_group\tquery\tthreshold\texpected_count\tafee_count\n"
                "real\tTiny\tdense\t1\t2\t0\t0\n"
                "real\tTiny\tdense\t2\t2\t0\t0\n",
                encoding="utf-8",
            )
            result = audit_result_quality.main(
                [str(summary), "--output-dir", str(root / "quality")]
            )
            self.assertEqual(result, 0)
            with (root / "quality" / "result_quality_issues.tsv").open(newline="") as source:
                codes = {row["code"] for row in csv.DictReader(source, delimiter="\t")}
            self.assertIn("ZERO_RESULT_RATIO_HIGH", codes)


if __name__ == "__main__":
    unittest.main()
