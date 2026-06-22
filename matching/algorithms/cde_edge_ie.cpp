#include "matching/run_matching.h"
#include "matching/algorithms/cde_edge_ie.h"

#ifndef CDE_EDGE_IE_ALGORITHM_KEY
#define CDE_EDGE_IE_ALGORITHM_KEY "cde_edge_ie"
#endif

#ifndef CDE_EDGE_IE_DISPLAY_NAME
#define CDE_EDGE_IE_DISPLAY_NAME "CDE-Edge-IE"
#endif

namespace ssm_ged {

    namespace {

        void run_cde_edge_ie(const Graph *query_graph, const Graph *data_graph, MatchResults &results, ui threshold)
        {
            Approximate_CDE_EdgeIE(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            CDE_EDGE_IE_ALGORITHM_KEY,
            CDE_EDGE_IE_DISPLAY_NAME,
            &run_cde_edge_ie
        };
        return definition;
    }

} // namespace ssm_ged
