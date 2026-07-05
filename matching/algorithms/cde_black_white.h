#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#include "matching/algorithms/cde_black_white/context.h"

inline void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph, std::vector<std::vector<std::pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    cde_black_white::MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
