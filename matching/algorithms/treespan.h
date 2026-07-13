#ifndef MATCHING_ALGORITHMS_TREESPAN_H_
#define MATCHING_ALGORITHMS_TREESPAN_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"

#include <memory>

class TreeSpanSolver {
public:
    struct TreeStats {
        long long total_time = 0;
        long long init_time = 0;
        long long filter_time = 0;
        long long index_time = 0;
        long long search_time = 0;
        long long enum_time = 0;
        long long verify_time = 0;
        long long recursion_calls = 0;
        long long candidate_range_hits = 0;
        long long candidate_range_misses = 0;
        long long candidate_edge_check_calls = 0;
        long long replacement_calls = 0;
        long long sequences_count = 0;
        long long enum_call_count = 0;
        long long duplicate_results = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
    } stats;

    TreeSpanSolver();
    ~TreeSpanSolver();

    bool init(const Graph *query_graph, const Graph *data_graph, ui threshold);
    void match(ssm_ged::MatchResults &results);
    void printStats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    TreeSpanSolver(const TreeSpanSolver &);
    TreeSpanSolver &operator=(const TreeSpanSolver &);
};

void Approximate_TreeSpan(const Graph *query_graph, const Graph *data_graph,
    ssm_ged::MatchResults &results, ui threshold);

#endif // MATCHING_ALGORITHMS_TREESPAN_H_
