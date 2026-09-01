# AFEC

This repository provides the implementation of **AFEC:
Anchor-Frontier Enumeration with Cost-Split Compression for Subgraph
Similarity Matching**. AFEC enumerates exact-label vertex mappings under a
missing-query-edge budget while requiring the preserved query graph to remain
connected. The repository also includes the paper's ablation variants and the
baseline implementations used in its experimental evaluation.

## Implemented algorithms

| Category | Display name | Executable | Benchmark key | Role |
| --- | --- | --- | --- | --- |
| Proposed | **AFEC** | `build/ssm_afec` | `afec` | Complete paper algorithm |
| Paper baseline | TreeSpan | `build/ssm_treespan` | `treespan` | Clean-room TreeSpan implementation |
| Paper baseline | DecQ | `build/ssm_decq` | `decq` | Clean-room DecQ implementation |
| Paper baseline | SASUM | `build/ssm_sasum` | `sasum` | Clean-room SASUM implementation |
| Additional implementation | S3AND-SSM | `build/ssm_s3and_ssm` | `s3and_ssm` | S3AND adapted to this repository's SSM semantics; not part of the paper's main comparison |

All executables use the same graph format and command-line interface:

```bash
build/ssm_afec -d path/to/graph_g.txt -q path/to/query.txt -t 2
```

## Paper ablation variants

The default build creates exactly the seven configurations used by the paper.
Their headers live in `configuration/afec_variants/`.

| Display name | Executable / key | Anchor-frontier | Cost-Split | Black--white policy | Anchor-edge order |
| --- | --- | --- | --- | --- | --- |
| **AFEC** | `ssm_afec` / `afec` | On | On | Dynamic | Dynamic |
| **AFE** | `ssm_afe` / `afe` | On | Off | Candidate-wise black expansion | Dynamic |
| **AFE-NoAF** | `ssm_afe_no_af` / `afe_no_af` | Off | Off | Candidate-wise expansion | Fixed query-vertex order |
| **AFEC-White** | `ssm_afec_white` / `afec_white` | On | On | Always white | Dynamic |
| **AFEC-Black** | `ssm_afec_black` / `afec_black` | On | On | Always black | Dynamic |
| **AFEC-Fixed** | `ssm_afec_fixed` / `afec_fixed` | On | On | Dynamic | TreeSpan-style static priorities |
| **AFEC-Random** | `ssm_afec_random` / `afec_random` | On | On | Dynamic | Reproducible random order |

These configurations form the paper's four component studies:

1. **Anchor-frontier expansion:** AFE vs. AFE-NoAF.
2. **Cost-Split compression:** AFEC vs. AFE.
3. **Dynamic black--white selection:** AFEC vs. AFEC-White vs. AFEC-Black.
4. **Dynamic anchor-edge ordering:** AFEC vs. AFEC-Fixed vs. AFEC-Random.

AFE and AFE-NoAF both disable Cost-Split. AFE immediately materializes every
candidate in an inclusion branch. AFE-NoAF additionally removes anchor-based
candidate restriction, follows a branch-independent query-vertex order, and
checks preserved-query connectivity only after a complete mapping is formed.

AFEC-Black keeps the AFEC configuration surface but always chooses concrete
black expansion. Consequently it explores the same candidate-wise successor
space as AFE; it remains a separate control for the black--white policy study.

All seven configurations enable TreeSpan's NLF filtering.

## Dependencies and build

The project requires CMake 3.16+, C++17, threads, and the Abseil C++ library
for `absl::flat_hash_map`. On Debian or Ubuntu, install Abseil with:

```bash
sudo apt-get update
sudo apt-get install -y libabsl-dev
```

Configure and build every algorithm and paper variant:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

To build only one AFEC-family configuration, select its header and build the
corresponding target:

```bash
cmake -S . -B build-afe \
  -DCMAKE_BUILD_TYPE=Release \
  -DSSM_CONFIG_HEADER=configuration/afec_variants/afe.h
cmake --build build-afe --target ssm_afe --parallel
```
