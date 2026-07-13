# DecQ clean-room reproduction

This repository contains a clean-room reproduction of the DecQ algorithm from:

> Gaoping Zhu, Ke Zhu, Wenjie Zhang, Xuemin Lin, and Chuan Xiao. “Efficient
> Subgraph Similarity All-Matching.” DASFAA 2012, LNCS 7238, pp. 455-469.
> DOI: 10.1007/978-3-642-29038-1_33.

The paper PDF supplied for this reproduction is the 2012 DecQ paper. It is not
the 2013 SASUM paper, and the `ssm_ged_sasum` executable is therefore a
different baseline.

The authors' publication pages provide the
[paper](https://cgi.cse.unsw.edu.au/~zhangw/files/2012DASFAA_allmatching.pdf)
and slides, but no source-code or binary link for DecQ. No public DecQ software
license was found. The implementation here is consequently based on the paper's
definitions, equations, two algorithms, and prose rather than on original
source code.

## Is it the same problem as this repository?

The feasibility predicate is the same. Both works use an undirected,
vertex-labeled simple query graph and one large data graph. A result is a full,
injective, label-preserving query-to-data vertex mapping. Extra data edges are
ignored, at most `threshold` query edges may be absent, and the preserved query
edges must connect all mapped query vertices.

The output identity is different:

| Dimension | DecQ paper | This repository |
|---|---|---|
| Internal result | One exact-match set `M_p` for every connected query pattern `p` | One maximal preserved-edge set induced by every complete mapping |
| Repeated mapping | May occur in several `M_p` sets | Reported once |
| Public count | Naturally counts `(pattern, mapping)` rows if the sets are flattened | Counts distinct complete mappings |

Algorithm 2 in the paper returns the collection `{M_p | p in FP(q, delta)}`.
The local executable retains those pattern-specific rows for DecQ processing
and reports their total as `Global-Pattern Rows`. For compatibility with the
other SSM-GED executables, it additionally applies a canonical maximal-output
adapter and reports `Unique Results`. A mapping is emitted only for the pattern
whose missing-edge key equals the edges actually absent under that mapping.
This changes result organization, not which complete mappings are feasible.

This distinction matters in an efficiency study: compare algorithms only after
their `Unique Results` counts agree, and also report DecQ's
`Global-Pattern Rows` because they are work performed by the paper algorithm.

## Reproduced algorithm

The `ssm_ged_decq` executable implements the paper's three phases:

1. Build the vertex and edge mapping statistics used by Equations (1)-(3).
   `GenOrder` greedily chooses the next query vertex that minimizes the
   estimated match count of the induced prefix.
2. Construct an edge-disjoint binary decomposition tree by recursive connected
   bisection. A split must be approximately balanced, leave both fragments
   connected, have fragments of at least `delta + 1` edges, and have positive
   estimated gain under Equation (6). Interpreting the paper's unspecified
   “fragment size” as edge count is a documented clean-room choice.
3. Enumerate every local pattern of every leaf fragment after zero through
   `delta` edge deletions. A local deletion key that disconnects the full query
   is discarded. Isolated vertices are removed, and disconnected local patterns
   retain a separate exact-match relation for each component.
4. Enumerate the connected global-pattern lattice through level `delta` and
   visit it from level `delta` back to zero, in lexicographic missing-edge order.
5. Compute minimal-pattern matches with cached hash equi-joins along the
   decomposition tree. Disconnected relations are not Cartesian-producted
   until a later fragment bridges them.
6. Compute each non-minimal pattern by selecting its child with the smallest
   materialized match set and validating the one restored query edge.
7. Apply the maximal-output adapter described above and expose phase times,
   pattern/relation counts, intermediate rows, paper-level rows, and unique
   mappings.

## Build and run

Configure the project after adding the new source file so CMake discovers the
target:

```bash
cmake -S . -B build_decq -DCMAKE_BUILD_TYPE=Release
cmake --build build_decq \
  --target ssm_ged_decq ssm_ged_treespan \
  --parallel
```

Run it through the common graph format and command-line interface:

```bash
./build_decq/ssm_ged_decq \
  -d test/datasets/real_graphs/Lastfm/graph_g.txt \
  -q test/datasets/real_graphs/Lastfm/query_graph/query_10_2.txt \
  -t 2
```

The normal comparison runner can include the new executable:

```bash
./compare.sh --build-dir build_decq --skip-build \
  --algorithms decq,treespan
```

## Correctness and repeated timing

The existing small exhaustive oracle is semantics-generic despite its historic
directory name. It independently enumerates every injective label-preserving
mapping and checks the missing-edge and connectivity predicates:

```bash
python3 tests/sasum_oracle/run_oracle.py \
  build_decq/ssm_ged_decq \
  build_decq/ssm_ged_treespan

python3 tests/decq_structure_check.py build_decq/ssm_ged_decq
```

The second check fixes the five-edge diamond at `delta = 2`: 14 connected
global patterns, 8 minimal patterns, 18 paper-level `(pattern, mapping)` rows,
and 2 unique maximal mappings.

For a count-checked repeated comparison, use:

```bash
python3 scripts/benchmark_decq.py \
  --build-dir build_decq \
  --data test/datasets/real_graphs/Lastfm/graph_g.txt \
  --query test/datasets/real_graphs/Lastfm/query_graph/query_10_2.txt \
  --thresholds 0 1 2 3 \
  --algorithms decq,treespan \
  --repetitions 7 \
  --output result/decq_lastfm.tsv
```

The driver starts a fresh process for every observation, rotates execution
order, rejects resource-limited runs, requires identical unique-result counts,
and reports median/mean/min/max time. Use `algorithm_ms_*` for efficiency
comparisons: it is parsed from each algorithm's first internal `Total Time` and
includes algorithm-specific indexing and preprocessing, but excludes common
graph parsing and statistics formatting. `runner_run_ms_*` is retained as a
diagnostic and includes `printStats()` output overhead; it should not be the
primary comparison metric. The TSV also records DecQ phase-time medians and
structural counters. `decq_core_ms_median` subtracts only the repository's
canonical output adapter and is the closer DecQ-paper core time;
`algorithm_ms_*` includes that adapter and is the fairer metric when every
algorithm must return this repository's unique maximal mappings.

`produced_intermediate_rows` is `Local Search Partial Rows + Materialized Match
Rows`, the closest local counter to the paper's intermediate-result metric.
`Global-Pattern Rows` is reported separately and must not be added again because
those final pattern relations are already materialized rows.

## Reproduction limits

Several implementation choices required deterministic clean-room substitutes:

- The paper says that local matching extends QuickSI and applies SAPPER's
  candidate index, but does not publish that modified matcher. This version uses
  a deterministic depth-first exact matcher with label, degree, and
  neighborhood-label-frequency filtering plus the published `GenOrder`.
- The paper requires approximately balanced recursive bisection with positive
  gain, but does not specify how candidate bisections or ties are generated.
  This version enumerates deterministic connected cuts of line-graph spanning
  trees, considers only the globally best-balanced candidates, then takes the
  largest positive estimated gain. It stops rather than falling back to a less
  balanced positive-gain split.
- Hash-table layout, component-join order, lattice tie-breaking, numeric
  safeguards, and memory management are not specified. They are deterministic
  here and are reported where they affect observable counters.
- The paper does not define a cross-pattern flattening/deduplication rule. The
  repository adapter is deliberately separated from `Global-Pattern Rows`.
- Original query-generation seeds, source code, raw measurements, and the
  authors' SAPPER binary are unavailable. Exact 2012 chart values therefore
  cannot be claimed.
- The repeated driver deliberately starts a cold process per observation, so
  every algorithm rebuilds its own data index. The paper's 100-query runs shared
  a SAPPER index; the two response-time conventions are not numerically
  interchangeable.
- The default `DECQ_INTERMEDIATE_MATCH_LIMIT` is one million simultaneously
  live rows as an out-of-memory guard. The paper's plots can exceed that scale;
  any run marked `(reached)` is incomplete and the benchmark rejects it. Use a
  dedicated build with, for example,
  `-DSSM_GED_COMPILE_DEFINITIONS=DECQ_INTERMEDIATE_MATCH_LIMIT=5000000`, or `0`
  for an explicitly unlimited run after checking available memory.

The paper used GCC 4.3.2 with `-O3` on one 2.40 GHz Xeon machine with 4 GB RAM
and Debian 4.1.1-21. Its HPRD graph had 9,460 vertices, 37,081 edges, and 307
vertex labels. Synthetic defaults were 5,000 data vertices, average data degree
12, 100 labels, 40 query edges, average query degree 4, and `delta = 2`; each
query set contained 100 induced-subgraph queries with one to three noisy edges.
Current-repository measurements are suitable for a controlled relative
comparison, not a numerical recreation of those old plots. The default
`decq,treespan` comparison also differs from the paper's RO-ND/EO-ND/SAPPER
baselines.

For a fair report, keep compiler flags, hardware, timeout, output/intermediate
limits, dataset, query, and process policy identical. Include unique-result
agreement, `Global-Pattern Rows`, peak live rows, decomposition/local/global
phase times, and whether any configured limit was reached.
