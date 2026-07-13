# SSM-GED

This repository provides four comparable baselines under one exact-label,
missing-query-edge, preserved-query-connectivity result semantics.

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
  --target ssm_treespan ssm_s3and_ssm ssm_sasum ssm_decq \
  --parallel
```

See [docs/experimental_protocol.md](docs/experimental_protocol.md) for timing
definitions, count-checked repeated runs, provenance fields, single-process
execution, and release tagging.
