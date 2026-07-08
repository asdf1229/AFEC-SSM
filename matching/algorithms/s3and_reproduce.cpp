#include "matching/run_matching.h"
#include "matching/algorithms/s3and_reproduce.h"

namespace ssm_ged {

    namespace {

        void run_s3and_reproduce(const Graph *query_graph, const Graph *data_graph,
            MatchResults &results, ui threshold)
        {
            Approximate_S3AND_Reproduce(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            "s3and_reproduce",
            "S3AND-Reproduce",
            &run_s3and_reproduce
        };
        return definition;
    }

} // namespace ssm_ged
