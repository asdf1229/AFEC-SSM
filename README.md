# SSM-GED

This repository provides AFEE and four comparable baselines under one exact-label,
missing-query-edge, preserved-query-connectivity result semantics.

## AFEE

The primary algorithm builds as `ssm_afee`. Its ablation executables follow
the experiment names directly:

- `ssm_afee_no_split` (`AFEE_NO_SPLIT`)
- `ssm_afee_mat` (`AFEE_MAT`)
- `ssm_afee_batch` (`AFEE_BATCH`)
- `ssm_afee_random_order` (`AFEE_RANDOM_ORDER`)
- `ssm_afee_fixed_order` (`AFEE_FIXED_ORDER`)
- `ssm_afee_only_nlf` (`AFEE_ONLY_NLF`)

## TreeSpan baseline

`ssm_treespan` is a deterministic clean-room TreeSpan implementation with
initial/replacement spanning trees, exclusion sets, QISequences, prefix-shared
search, and unique-result validation.

See [docs/treespan_reproduction.md](docs/treespan_reproduction.md) for the
algorithm mapping, deterministic implementation choices, and fidelity limits.

## S3AND-SSM baseline

`ssm_s3and_ssm` adapts S3AND to SSM-GED's single exact label, total missing
query-edge distance, and preserved-query-edge connectivity. It must be named
**S3AND-SSM**, not an unqualified S3AND reproduction.

See [docs/s3and_adaptation.md](docs/s3and_adaptation.md) for the semantic
differences, fixed parameters, timing scope, and fidelity limits.

## DecQ baseline

`ssm_decq` is a clean-room reproduction of Zhu et al.'s DASFAA 2012
algorithm from “Efficient Subgraph Similarity All-Matching.” It preserves the
paper's per-pattern work counters while adapting its output to this repository's
one-maximal-result-per-vertex-mapping convention.

See [docs/decq_reproduction.md](docs/decq_reproduction.md) for the problem
comparison, fidelity limits, and repeated benchmark commands.

## SASUM baseline

This repository includes a clean-room reproduction of the 2013 SASUM
sharing-based approximate subgraph matcher.  It builds as
`build/ssm_sasum` and uses the same graph format and `-d/-q/-t` interface as
the other algorithms.

See [docs/sasum_reproduction.md](docs/sasum_reproduction.md) for the problem
comparison, fidelity limits, and repeated benchmark commands.

## flat_hash_map Dependency

This project uses `absl::flat_hash_map`, which is provided by the Abseil C++ library. Before building SSM-GED, install Abseil on your system.

Linux download/install commands:

```bash
sudo apt-get update
sudo apt-get install -y libabsl-dev
```

## Reproducible build

The project requires C++17. For a local Release build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build \
  --target ssm_afee ssm_treespan ssm_s3and_ssm ssm_sasum ssm_decq \
  --parallel
```

## Benchmark planning and resume

Use `--dry-run` to build and validate the discovered case/task plan without
executing algorithms or creating the target result directory:

```bash
./compare.sh --dry-run --skip-build \
  --build-dir build \
  --data-dir test/experiment_main \
  --algorithms afee,treespan,decq,sasum,s3and_ssm \
  --thresholds 2 \
  --result-dir result/main_comparison
```

Use `--result-dir` for a stable experiment location. A non-empty directory is
never overwritten implicitly:

```bash
./compare.sh --skip-build \
  --build-dir build \
  --data-dir test/experiment_main \
  --algorithms afee,treespan,decq,sasum,s3and_ssm \
  --thresholds 2 \
  --result-dir result/main_comparison
```

At the start of a real run, the selected executables are copied into
`result/main_comparison/bin/`. The task manifest runs only these snapshots, so
later rebuilds cannot change an experiment that is already in progress.

Use `--thresholds 2` for the main comparison and AFEE ablations. Use
`--thresholds 1,2,3,4,5,6` for the two-dataset threshold experiment. An
explicit list overrides all query-group threshold defaults and is recorded in
`run_config.tsv` for resume validation.

While tasks are running, `progress.tsv` in the result directory is atomically
replaced every 30 seconds. It contains the current `completed`, `succeeded`,
`failed`, and `pending` task counts plus a session-rate ETA. This file is meant
for background-run monitoring and is not periodically printed to the terminal
or `out.log`.

After an interruption, pass the same experiment arguments with `--resume`.
The script validates the saved manifests, timeouts, parallelism, runner hash,
and hashes of the executable snapshots under `result-dir/bin/`; the external
build directory is not executed during resume. It skips every task with a
valid terminal status and reruns only missing, incomplete, corrupt, or
output-less tasks. Terminal
failures such as `Timeout` are retained rather than retried.

```bash
./compare.sh \
  --build-dir build \
  --data-dir test/experiment_main \
  --algorithms afee,treespan,decq,sasum,s3and_ssm \
  --thresholds 2 \
  --result-dir result/main_comparison \
  --resume
```

See [docs/experimental_protocol.md](docs/experimental_protocol.md) for timing
definitions, count-checked repeated runs, provenance fields, single-process
execution, and release tagging.
