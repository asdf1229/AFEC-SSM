#include "matching/run_matching.h"
#include "matching/algorithms/s3and.h"

namespace ssm_ged {

    namespace {

        void run_s3and(const Graph *query_graph, const Graph *data_graph,
            MatchResults &results, ui threshold)
        {
            Approximate_S3AND(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            "s3and",
            "S3AND",
            &run_s3and
        };
        return definition;
    }

} // namespace ssm_ged
