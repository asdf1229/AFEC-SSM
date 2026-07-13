#include "matching/algorithms/sasum.h"

namespace ssm_ged {

namespace {

void run_sasum(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold)
{
    Approximate_SASUM(query_graph, data_graph, results, threshold);
}

} // namespace

const AlgorithmDefinition &create_algorithm_definition()
{
    static const AlgorithmDefinition definition = {
        "sasum",
        "SASUM",
        &run_sasum
    };
    return definition;
}

} // namespace ssm_ged
