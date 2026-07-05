#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#include "matching/algorithms/cde_black_white/context.h"

inline void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph, std::vector<std::vector<std::pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    cde_black_white::Workspace workspace;
    if (workspace.init(query_graph, data_graph, threshold)) {
        workspace.match(M_ANS);
    }

    workspace.stats.total_time = t_total.elapsed();
    workspace.printStats();
    ssm_ged::set_reported_result_count(workspace.stats.result_count);
}

#endif
