#ifndef SSM_GED_MATCHING_RUN_MATCHING_H_
#define SSM_GED_MATCHING_RUN_MATCHING_H_

#include "graph/graph.h"

namespace ssm_ged {

    using MatchResults = std::vector<std::vector<std::pair<ui, ui> > >;
    using AlgorithmEntry = void (*)(const Graph *, const Graph *, MatchResults &, ui);

    struct AlgorithmDefinition {
        std::string key;
        std::string display_name;
        AlgorithmEntry entry;
    };

    const AlgorithmDefinition &create_algorithm_definition();
    int run_algorithm_main(int argc, char *argv[], const AlgorithmDefinition &algorithm);

} // namespace ssm_ged

#endif
