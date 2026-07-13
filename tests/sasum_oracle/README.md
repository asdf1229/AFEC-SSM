# SASUM micro-oracle

These fixtures exercise the semantics shared by SASUM and SSM-GED:

- exact matching (`theta = 0`);
- one missing cycle edge;
- connectivity of the preserved query edges;
- injectivity and distinct automorphic mappings;
- extension of a pruned leaf base graph.
- a two-level K4 lattice whose second occurrence needs exactly two deletions.
- the paper's Figure 2 query shape, including its `14/8/8/4`
  query-subgraph/terminal/base/seed structural counts at `theta = 2`.

The over-threshold triangle/path rows also verify the common cycle-rank cap:
increasing `theta` beyond the maximum number of deletions compatible with a
connected spanning match must not discard the exact and lower-distance answers.

Run the oracle after building the algorithms:

```bash
python3 tests/sasum_oracle/run_oracle.py \
  build/ssm_ged_sasum \
  build/ssm_ged_treespan
```

The runner computes an independent brute-force answer for each fixture, checks
it against `expected.tsv`, and then compares every supplied executable with the
oracle count.
