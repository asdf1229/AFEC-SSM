#ifndef MATCHING_ALGORITHMS_SASUM_H_
#define MATCHING_ALGORITHMS_SASUM_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"

#include <memory>

namespace ssm {

class SASUMSolver {
public:
    struct SASUMStats {
        long long total_time = 0;
        long long data_index_time = 0;
        long long lattice_time = 0;
        long long base_generation_time = 0;
        long long seed_selection_time = 0;
        long long exact_matching_time = 0;
        long long derivation_time = 0;
        size_t query_subgraph_count = 0;
        size_t terminal_graph_count = 0;
        size_t base_graph_count = 0;
        size_t seed_graph_count = 0;
        size_t exact_matching_executions = 0;
        size_t fallback_exact_executions = 0;
        size_t exact_embeddings = 0;
        size_t seed_match_rows = 0;
        size_t query_pattern_match_rows = 0;
        size_t duplicate_outputs = 0;
        size_t result_count = 0;
        size_t peak_live_match_rows = 0;
        bool output_limit_reached = false;
        bool intermediate_limit_reached = false;
    } stats;

    SASUMSolver();
    ~SASUMSolver();

    bool init(const Graph *query_graph, const Graph *data_graph, ui threshold);
    void match(MatchResults &results);
    void printStats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    SASUMSolver(const SASUMSolver &);
    SASUMSolver &operator=(const SASUMSolver &);
};

void Approximate_SASUM(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold);

} // namespace ssm

#endif // MATCHING_ALGORITHMS_SASUM_H_
