#ifndef MATCHING_ALGORITHMS_S3AND_H_
#define MATCHING_ALGORITHMS_S3AND_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"

#include <memory>

class S3ANDSolver {
public:
    struct S3ANDStats {
        long long total_time = 0;
        long long init_time = 0;
        long long filter_time = 0;
        long long index_time = 0;
        long long dfs_time = 0;
        unsigned long long entry_node_visits = 0;
        unsigned long long entry_node_prunes = 0;
        unsigned long long vertex_visits = 0;
        unsigned long long candidate_checks = 0;
        unsigned long long recursion_calls = 0;
        unsigned long long prune_calls = 0;
        unsigned long long final_checks = 0;
        ui filter_candidate_count = 0;
        ui index_entry_count = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
    } stats;

    S3ANDSolver();
    ~S3ANDSolver();

    bool init(const Graph *query_graph, const Graph *data_graph, ui threshold);
    void match(ssm_ged::MatchResults &results);
    void printStats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    S3ANDSolver(const S3ANDSolver &);
    S3ANDSolver &operator=(const S3ANDSolver &);
};

void Approximate_S3AND(const Graph *query_graph, const Graph *data_graph,
    ssm_ged::MatchResults &results, ui threshold);

#endif // MATCHING_ALGORITHMS_S3AND_H_
