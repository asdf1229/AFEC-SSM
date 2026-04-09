#include "matching/run_matching.h"
#include "matching/algorithms/cde_match.h"

namespace ssm_ged {

    namespace {

        void run_cde_match(const Graph *query_graph, const Graph *data_graph, MatchResults &results, ui threshold)
        {
            Approximate_Matching(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            "cde_match",
            "CDE-Match",
            &run_cde_match
        };
        return definition;
    }

} // namespace ssm_ged
