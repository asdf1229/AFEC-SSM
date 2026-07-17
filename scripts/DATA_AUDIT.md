# Data audit tools

The audit is intentionally split into static data checks, experiment-plan
checks, and post-run result-quality checks.  Static hard errors should be fixed
before benchmarking; result-distribution findings are warnings.

## 1. Static dataset and query audit

Run from the repository root:

```bash
python3 -B scripts/audit_datasets.py \
  --data-root test/datasets \
  --output-dir data_audit
```

To enforce the intended roles, name the selected datasets explicitly:

```bash
python3 -B scripts/audit_datasets.py \
  --data-root test/datasets \
  --output-dir data_audit \
  --main-datasets DATASET_1,DATASET_2,DATASET_3,DATASET_4,DATASET_5,DATASET_6,DATASET_7,DATASET_8,DATASET_9,DATASET_10 \
  --threshold-datasets DATASET_1,DATASET_2 \
  --strict
```

The static audit writes:

- `dataset_inventory.tsv`: data-graph size, connectivity, degree, density,
  labels, self-loops, duplicate edges, raw-file SHA-256, and a label-agnostic
  topology SHA-256 that reveals renamed or relabeled versions of one graph.
- `query_inventory.tsv`: one query per row, including connectivity, size,
  density, labels, triangles, `slack`, every requested/effective theta from 1
  through 6, lower bounds from absent endpoint-label pairs, static zero-result
  risks, raw hash, exact labeled-isomorphism class, and validation of query-group
  declarations such as `sparse_8` (`d <= 3`) and `dense_8` (`d > 3`). A
  size-only group such as `8` accepts any density. Edge-label combinations are
  inventoried but are informational because the current matcher discards edge
  labels.
- `data_audit_summary.tsv`: per-dataset gates and `slack=0,1,...,5,6+`
  distributions, truncation ratios, and effective-theta distributions.
- `query_group_audit_summary.tsv`: per-dataset/per-query-group counts, actual
  vertex/edge/mean-degree ranges, size and density violations, and the
  `slack>=4` gate.
- `sha256_manifest.tsv`: fixed data/query file checksum inventory.

The default isomorphism check uses only the standard library.  It uses a
Weisfeiler-Lehman invariant as a prefilter and then runs an exact
vertex-label-preserving isomorphism check, so the prefilter is never treated as
proof.  `--skip-isomorphism` is available for a quick first pass, but leaves the
isomorphism uniqueness gate unverified.

`main_ready` requires exactly 800 queries, unique file contents, unique labeled
graphs, connectivity, a recognized and satisfied query-group size/density
declaration, `slack>=4`, and static label feasibility. `threshold_ready` uses
the same gates with `slack>=6`.

## 2. Experiment-plan audit

The plan audit accepts either a headerless `tasks.tsv` produced by
`compare.sh`, or a canonical headered TSV.  Prefix a `compare.sh` task file with
its experiment role:

```bash
python3 -B scripts/audit_experiment_plan.py \
  main=result/main/tasks.tsv \
  threshold=result/threshold/tasks.tsv \
  ablation=result/ablation/tasks.tsv \
  --output-dir experiment_plan_audit \
  --strict
```

The canonical format is useful when plans from multiple runners must be
combined.  Required columns are:

```text
experiment	variant	dataset	query_id	data_graph	query_graph	requested_theta
```

`experiment` is `main`, `threshold`, or `ablation`.  `variant` is an algorithm
name for the main run and an ablation case name for the ablation run.  Paths may
be different only when their SHA-256 contents are identical.

If preprocessing produces a reusable artifact, add an optional
`preprocessing_artifact` column. Its bytes are then hashed and required to match
across variants of the same case. Headerless `compare.sh` tasks do not provide
this metadata, so preprocessing identity cannot be claimed from those files
alone.

The tool verifies:

- 10/2/10 dataset counts for main/threshold/ablation;
- 800 unique query keys per dataset and variant;
- theta 2 for every main and ablation query;
- the full theta set 1 through 6 for every threshold query;
- identical data/query hashes across algorithms and variants;
- identical query sets across variants;
- identical query order across variants;
- identical main and ablation query bytes.

It writes a resolved plan with hashes, a per-variant summary, and a detailed
issues TSV.

## 3. Result-quality audit

After a benchmark run:

```bash
python3 -B scripts/audit_result_quality.py \
  result/RUN/summary.tsv \
  --output-dir result/RUN
```

The tool writes per-query-group rows plus `query_group=__ALL__` dataset totals
to `result_quality_summary.tsv`, and writes
`result_quality_issues.tsv`.  It warns when zero results exceed 80%, output-cap
hits exceed 20%, or result counts are missing.  These thresholds can be changed
with `--zero-warning-ratio`, `--cap-warning-ratio`, and `--output-limit`.

## Metadata-dependent checks

The current query files contain renumbered query vertices but no planted source
mapping or perturbation log.  Therefore the static inventory deliberately marks
the following as unavailable instead of claiming they passed:

- validation of a saved planted embedding;
- original data-vertex overlap between queries;
- number of queries generated from each planted source;
- actual modified edges and minimum required perturbation theta.

Add those fields to generator metadata before relying on the corresponding
experimental claims.

## 4. Repair queries that fail acceptance gates

`repair_queries.py` leaves qualifying files byte-for-byte unchanged and
regenerates only queries selected by the requested repair gates. Replacements are
reproducible random-walk-with-restart samples, written as the complete induced
subgraph on the sampled data vertices.

Stage repairs without touching the source dataset:

```bash
python3 -B scripts/repair_queries.py \
  --data-root test/datasets \
  --staging-dir tmp/query_repair_staging \
  --backup-dir tmp/query_repair_backup \
  --min-slack 4 \
  --seed 20260717
```

Independently verify every retained hash and every regenerated planted mapping.
For regenerated queries this checks injectivity, labels, connectivity, slack,
the declared query-group size/density rule, and both directions of the
induced-edge condition:

```bash
python3 -B scripts/repair_queries.py \
  --data-root test/datasets \
  --staging-dir tmp/query_repair_staging \
  --min-slack 4 \
  --verify-only
```

After `query_repair_verification.tsv` contains only `OK`, back up and atomically
replace only the regenerated files:

```bash
python3 -B scripts/repair_queries.py \
  --data-root test/datasets \
  --staging-dir tmp/query_repair_staging \
  --backup-dir tmp/query_repair_backup \
  --apply-staged
```

The manifest records old/new hashes, graph sizes, slack values, per-query seeds,
attempt counts, duplicate rejections, and the planted data-vertex mapping.

To keep the first member of each exact labeled-isomorphism class and regenerate
all redundant members, add:

```bash
--repair-isomorphic-duplicates
```

To repair named group violations (`sparse_N` means `d<=3`, `dense_N` means
`d>3`, and `N` means any density), add:

```bash
--repair-query-group-violations
```

Isomorphism repairs retain the original vertex count and target the original
slack. Low-slack repairs target `--min-slack` subject to their group density
constraint. Density repairs use the nearest valid density boundary.
