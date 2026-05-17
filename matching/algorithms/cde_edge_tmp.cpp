#include "matching/run_matching.h"
#include "matching/algorithms/cde_edge_tmp.h"

namespace ssm_ged {

    namespace {

        void run_cde_edge_tmp(const Graph *query_graph, const Graph *data_graph,
            MatchResults &results, ui threshold)
        {
            Approximate_CDE_EdgeIE(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            "cde_edge_tmp",
            "CDE-Edge-TMP",
            &run_cde_edge_tmp
        };
        return definition;
    }

} // namespace ssm_ged
