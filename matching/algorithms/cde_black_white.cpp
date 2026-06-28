#include "matching/run_matching.h"
#include "matching/algorithms/cde_black_white.h"

#ifndef CDE_BLACK_WHITE_ALGORITHM_KEY
#define CDE_BLACK_WHITE_ALGORITHM_KEY "cde_black_white"
#endif

#ifndef CDE_BLACK_WHITE_DISPLAY_NAME
#define CDE_BLACK_WHITE_DISPLAY_NAME "CDE-Black-White"
#endif

namespace ssm_ged {

    namespace {

        void run_cde_black_white(const Graph *query_graph, const Graph *data_graph, MatchResults &results, ui threshold)
        {
            Approximate_CDE_BlackWhite(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            CDE_BLACK_WHITE_ALGORITHM_KEY,
            CDE_BLACK_WHITE_DISPLAY_NAME,
            &run_cde_black_white
        };
        return definition;
    }

} // namespace ssm_ged
