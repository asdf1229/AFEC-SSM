# SSM-GED

## DecQ baseline

`ssm_ged_decq` is a clean-room reproduction of Zhu et al.'s DASFAA 2012
algorithm from “Efficient Subgraph Similarity All-Matching.” It preserves the
paper's per-pattern work counters while adapting its output to this repository's
one-maximal-result-per-vertex-mapping convention.

See [docs/decq_reproduction.md](docs/decq_reproduction.md) for the problem
comparison, fidelity limits, correctness checks, and repeated benchmark
commands.

## SASUM baseline

This repository includes a clean-room reproduction of the 2013 SASUM
sharing-based approximate subgraph matcher.  It builds as
`build/ssm_ged_sasum` and uses the same graph format and `-d/-q/-t` interface as
the other algorithms.

See [docs/sasum_reproduction.md](docs/sasum_reproduction.md) for the problem
comparison, fidelity limits, correctness oracle, and repeated benchmark
commands.

## flat_hash_map Dependency

This project uses `absl::flat_hash_map`, which is provided by the Abseil C++ library. Before building SSM-GED, install Abseil on your system.

Linux download/install commands:

```bash
sudo apt-get update
sudo apt-get install -y libabsl-dev
```
