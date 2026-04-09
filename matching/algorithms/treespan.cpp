#include "matching/run_matching.h"
#include "matching/algorithms/treespan.h"

namespace ssm_ged {

    namespace {

        void run_treespan(const Graph *query_graph, const Graph *data_graph, MatchResults &results, ui threshold)
        {
            Approximate_TreeSpan(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            "treespan",
            "TreeSpan",
            &run_treespan
        };
        return definition;
    }

} // namespace ssm_ged
