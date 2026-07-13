#ifndef MATCHING_ALGORITHMS_DECQ_H_
#define MATCHING_ALGORITHMS_DECQ_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"

#include <cstdint>
#include <memory>

#ifndef DECQ_INTERMEDIATE_MATCH_LIMIT
#define DECQ_INTERMEDIATE_MATCH_LIMIT 1000000
#endif

#if DECQ_INTERMEDIATE_MATCH_LIMIT < 0
#error "DECQ_INTERMEDIATE_MATCH_LIMIT must be non-negative."
#endif

namespace ssm_ged {

// A clean-room implementation of the DecQ execution plan described in
// "Efficient Subgraph Similarity All-Matching" (DASFAA 2012).  The public
// statistics intentionally distinguish the paper's per-global-pattern rows
// from this repository's canonical, unique full mappings.
class DecQSolver {
public:
    struct DecQStats {
        long long total_time = 0;
        long long data_index_time = 0;
        long long decomposition_time = 0;
        long long local_pattern_time = 0;
        long long local_matching_time = 0;
        long long lattice_time = 0;
        long long minimal_merge_time = 0;
        long long edge_validation_time = 0;
        long long canonical_dedup_time = 0;

        size_t query_vertices = 0;
        size_t query_edges = 0;
        size_t indexed_labels = 0;
        size_t decomposition_nodes = 0;
        size_t decomposition_leaves = 0;
        size_t accepted_splits = 0;
        size_t split_candidates = 0;
        size_t local_patterns = 0;
        size_t local_edge_cuts_skipped = 0;
        size_t disconnected_local_patterns = 0;
        size_t local_components = 0;
        size_t exact_match_executions = 0;
        size_t local_search_partial_rows = 0;
        size_t local_match_rows = 0;
        size_t lattice_candidates = 0;
        size_t global_patterns = 0;
        size_t global_edge_cuts_skipped = 0;
        size_t minimal_patterns = 0;
        size_t nonminimal_patterns = 0;
        size_t minimal_patterns_processed = 0;
        size_t nonminimal_patterns_processed = 0;
        size_t decomposition_cache_hits = 0;
        size_t decomposition_cache_misses = 0;
        size_t hash_joins = 0;
        size_t hash_join_input_rows = 0;
        size_t hash_join_output_rows = 0;
        size_t edge_validations = 0;
        size_t materialized_match_rows = 0;
        size_t peak_live_match_rows = 0;
        size_t global_pattern_rows = 0;
        size_t noncanonical_pattern_rows = 0;
        size_t duplicate_canonical_mappings = 0;
        size_t invalid_global_rows = 0;
        size_t result_count = 0;

        double minimum_edge_selectivity = 0.0;
        double average_edge_selectivity = 0.0;
        double maximum_edge_selectivity = 0.0;
        double unsplit_estimated_cost = 0.0;
        double decomposition_estimated_cost = 0.0;

        uint64_t result_digest_sum = 0;
        uint64_t result_digest_xor = 0;
        bool output_limit_reached = false;
        bool intermediate_limit_reached = false;
    } stats;

    DecQSolver();
    ~DecQSolver();

    bool init(const Graph *query_graph, const Graph *data_graph, ui threshold);
    void match(MatchResults &results);
    void printStats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    DecQSolver(const DecQSolver &);
    DecQSolver &operator=(const DecQSolver &);
};

void Approximate_DecQ(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold);

} // namespace ssm_ged

#endif // MATCHING_ALGORITHMS_DECQ_H_
