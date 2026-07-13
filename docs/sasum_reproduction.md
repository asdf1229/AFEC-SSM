# SASUM clean-room reproduction

This repository contains a clean-room implementation of Kim et al.,
"SASUM: A Sharing-Based Approach to Fast Approximate Subgraph Matching for
Large Graphs" (IEICE Transactions on Information and Systems, 2013).

The paper does not publish SASUM source code, pseudocode, its exact-matching
backend, compiler flags, or set-cover tie-breaking rules.  The implementation
therefore reproduces the algorithm described in Sections 3-4, not the authors'
unavailable executable bit for bit.

The official [J-STAGE record](https://www.jstage.jst.go.jp/article/transinf/E96.D/3/E96.D_624/_article/-char/en)
lists no supplementary material, the [KAIST repository record](https://koasas.kaist.ac.kr/handle/10203/174317)
provides the article rather than a software archive, and the [first author's
publication page](https://songhyonkim.github.io/) does not link SASUM code.  No
public software license could therefore be identified.

## Is SASUM the same problem as SSM-GED?

Yes, at the problem level used by this repository.  Both systems enumerate
injective, label-preserving mappings from every query vertex into one large
data graph.  They allow only missing query edges, ignore extra data edges, and
require the query edges preserved by a mapping to form a connected graph.

They are not implementations of general graph edit distance: neither permits
vertex insertion/deletion, edge-insertion cost, label substitution, weighted
edits, or top-k/best-only answers.

There are two interface-level differences:

| Dimension | SASUM paper | This repository |
|---|---|---|
| Result organization | Maintains `M(q, G)` for every relaxed query subgraph `q` | Reports the union, deduplicated by the full query-to-data vertex mapping |
| Threshold | Assumes `theta > 0` and a cyclic query | Also supports `theta = 0`; the common runner caps `theta` at `|E(Q)| - |V(Q)| + 1` because a connected answer cannot delete more edges |

The deduplicated union is semantics-preserving for this repository's maximal
answer definition.  Internal per-pattern match sets are still retained long
enough to execute SASUM's sharing strategy.

## Implemented paper algorithm

The `ssm_ged_sasum` executable implements the following phases:

1. Generate the connected query-subgraph lattice for zero through `theta`
   deleted edges.  Patterns retain original query vertex and edge identities;
   label-isomorphic patterns are not merged.
2. Treat the connected patterns at the effective, cycle-rank-capped `theta` as
   terminal graphs.
3. Generate base graphs by pruning one terminal edge.  If pruning a leaf edge
   isolates a query vertex, remove that vertex as Definition 6 requires.
4. Deduplicate bases by active-query-vertex and query-edge masks and record the
   terminal graphs covered by each base.
5. Apply the paper's greedy set-cover approximation, with deterministic
   generation-order tie-breaking, to select seed graphs.
6. Enumerate each seed with a non-induced exact subgraph matcher using label,
   degree, and neighborhood-label-frequency filtering.
7. Derive terminal matches either by checking the restored edge or by extending
   a pruned leaf through neighbors of its mapped anchor.  If several selected
   seeds cover one terminal, use the seed with the smallest materialized match
   set; the paper permits any covering seed but does not specify this tie-break.
8. Traverse the lattice from fewer to more edges.  For each pattern, select the
   predecessor with the fewest matches and filter on the one restored edge.
9. Deduplicate the union by complete vertex mapping for compatibility with the
   other SSM-GED algorithms.

The executable reports lattice, base generation, seed selection, exact
matching, and derivation/dedup times, as well as pattern/base/seed counts,
seed-match rows, `sum_q |M(q,G)|` as query-pattern rows, unique result mappings,
and the number of exact-matcher executions.  These counters are useful because
the paper proves only that SASUM reduces the number of exact matching calls; it
does not prove a wall-clock bound.  A smaller seed may be much less selective
and can create a large intermediate match set.

## Build and run

```bash
./build.sh

./build/ssm_ged_sasum \
  -d test/datasets/real_graphs/Lastfm/graph_g.txt \
  -q test/datasets/real_graphs/Lastfm/query_graph/query_10_2.txt \
  -t 3
```

The standard dataset runner discovers the new executable automatically.  A
two-algorithm comparison can be started with:

```bash
./compare.sh --algorithms sasum,treespan
```

`compare.sh` now marks a case as failed when successful algorithms return
different result counts; timing rows with mismatched answers must not be used.

For repeated measurements of a specific query, use the count-checked benchmark
driver.  `run_ms` excludes common graph-file parsing but includes every
algorithm-specific index and preprocessing step.

```bash
python3 scripts/benchmark_sasum.py \
  --build-dir build \
  --data test/datasets/real_graphs/Lastfm/graph_g.txt \
  --query test/datasets/real_graphs/Lastfm/query_graph/query_10_2.txt \
  --thresholds 1 2 3 \
  --algorithms sasum,treespan \
  --repetitions 7 \
  --output result/sasum_lastfm.tsv
```

The script runs algorithms as separate cold processes, rotates their execution
order, rejects output-limit truncation, verifies identical counts, and reports
median/mean/min/max algorithm time.

## Correctness checks

The fixed micro-oracle covers exact matching, one missing cycle edge,
connectivity after a bridge is missing, repeated-label automorphisms,
injectivity, and leaf-base extension:

```bash
python3 tests/sasum_oracle/run_oracle.py \
  build/ssm_ged_sasum \
  build/ssm_ged_treespan
```

The test computes an independent brute-force mapping set for every fixture.  It
checks counts for every executable and, when an executable reports the SASUM
set digest, also checks an order-independent digest of the complete mapping
set.  The paper-example fixture additionally asserts the published
`14/8/8/4` query-subgraph/terminal/base/seed structure.

## Reproduction limits

- The original SASUM exact matcher is unspecified and unavailable.  The local
  backend is a deterministic backtracking enumerator, so absolute 2013 timing
  values are not directly reproducible.
- The paper minimizes only the number of seed graphs.  It does not use a
  selectivity- or cost-aware set-cover objective.
- The published synthetic-query generation, label assignment, random seeds,
  run repetitions, and raw measurements are absent.  Reproducing chart trends
  is realistic; claiming identical chart values is not.
- The paper's own defaults were `|V(G)|=5000`, average data degree 8, 250 labels,
  `|V(Q)|=20`, average query degree 3, and `theta=1`.  Its reported machine was
  a single-threaded C++ run on a 2.33 GHz Xeon E5345 with 8 GB RAM and Fedora 12.
- `MATCH_OUTPUT_LIMIT` stops unique-result enumeration and release builds do
  not retain the final `MatchResults` vectors.  SASUM must nevertheless
  materialize seed and current lattice-layer matches to perform sharing, so a
  non-selective seed can be much larger than the final output.  The independent
  `SASUM_INTERMEDIATE_MATCH_LIMIT` (one million simultaneously live mapping
  rows by default) aborts and marks such a run invalid for benchmarking.  Set
  either limit to `0` only when deliberately requesting an unlimited run.

To change the SASUM cap for a dedicated build:

```bash
cmake -S . -B build_sasum_large \
  -DSSM_GED_COMPILE_DEFINITIONS=SASUM_INTERMEDIATE_MATCH_LIMIT=5000000
cmake --build build_sasum_large --target ssm_ged_sasum --parallel
```

For fair experiments, report result-count agreement, exact-matcher calls,
phase times, seed-match rows, per-query-pattern match rows, unique-result count,
output-limit status, and wall-clock time.  Keep compiler flags and hardware
identical across algorithms.
