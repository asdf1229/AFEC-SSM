#include "matching/algorithms/afec/afec.h"

namespace afec {

#if !AFEC_ENABLE_ANCHOR_FRONTIER

void MatchingSolver::buildNoAFOrder(ui root)
{
    no_af_order.clear();
    no_af_order.reserve(qn);

    vector<unsigned char> ordered(qn, 0);
    vector<ui> ordered_neighbor_count(qn, 0);

    no_af_order.push_back(root);
    ordered[root] = 1;
    for (ui nbr : q_neighbors[root]) {
        ordered_neighbor_count[nbr]++;
    }

    while (no_af_order.size() < qn) {
        ui best = qn;
        for (ui u = 0; u < qn; ++u) {
            if (ordered[u]) continue;
            if (best == qn ||
                candidates[u].size() < candidates[best].size() ||
                (candidates[u].size() == candidates[best].size() &&
                 ordered_neighbor_count[u] > ordered_neighbor_count[best]) ||
                (candidates[u].size() == candidates[best].size() &&
                 ordered_neighbor_count[u] == ordered_neighbor_count[best] &&
                 q_degree[u] > q_degree[best]) ||
                (candidates[u].size() == candidates[best].size() &&
                 ordered_neighbor_count[u] == ordered_neighbor_count[best] &&
                 q_degree[u] == q_degree[best] && u < best)) {
                best = u;
            }
        }

        assert(best != qn);
        no_af_order.push_back(best);
        ordered[best] = 1;
        for (ui nbr : q_neighbors[best]) {
            if (!ordered[nbr]) ordered_neighbor_count[nbr]++;
        }
    }
}

ui MatchingSolver::noAFMissingDelta(const NoAnchorFrontierState &state,
    ui u, ui v, ui limit)
{
    ui missing = 0;
    for (ui nbr : q_neighbors[u]) {
        int mapped_v = state.mapped_q[nbr];
        if (mapped_v == -1) continue;

        stats.graph_has_edge_checks++;
        if (!data_graph->hasEdge(v, (ui)mapped_v)) {
            ++missing;
            if (missing > limit) return missing;
        }
    }
    return missing;
}

bool MatchingSolver::noAFPreservedQueryGraphIsConnected(
    const NoAnchorFrontierState &state)
{
    if (qn <= 1) return true;

    no_af_visited.assign(qn, 0);
    no_af_queue.clear();
    no_af_queue.reserve(qn);
    no_af_queue.push_back(0);
    no_af_visited[0] = 1;

    for (size_t pos = 0; pos < no_af_queue.size(); ++pos) {
        ui u = no_af_queue[pos];
        assert(state.mapped_q[u] != -1);
        for (ui nbr : q_neighbors[u]) {
            if (no_af_visited[nbr]) continue;
            assert(state.mapped_q[nbr] != -1);
            stats.graph_has_edge_checks++;
            if (data_graph->hasEdge((ui)state.mapped_q[u],
                    (ui)state.mapped_q[nbr])) {
                no_af_visited[nbr] = 1;
                no_af_queue.push_back(nbr);
            }
        }
    }

    return no_af_queue.size() == qn;
}

void MatchingSolver::emitNoAFResult(const NoAnchorFrontierState &state)
{
    assert(state.part_M.size() == qn);
    stats.result_count++;
    ssm::report_result_progress(stats.result_count, stats.recursion_calls);
    if (outputLimitReached()) stats.output_limit_reached = true;
#ifndef NDEBUG
    results_ptr->push_back(state.part_M);
#endif
}

void MatchingSolver::searchNoAF(NoAnchorFrontierState &state, ui depth,
    ui cost)
{
    if (outputLimitReached()) return;
    stats.recursion_calls++;

    if (depth == qn) {
        if (noAFPreservedQueryGraphIsConnected(state)) {
            emitNoAFResult(state);
        }
        else {
            stats.prun_calls++;
        }
        return;
    }

    ui u = no_af_order[depth];
    bool branched = false;
    for (ui v : candidates[u]) {
        if (state.used_data_flag[v]) continue;

        ui remaining_budget = threshold - cost;
        ui delta = noAFMissingDelta(state, u, v, remaining_budget);
        if (delta > remaining_budget) {
            stats.prun_calls++;
            continue;
        }

        branched = true;
        state.mapped_q[u] = (int)v;
        state.used_data_flag[v] = 1;
        state.part_M.push_back({ u, v });

        searchNoAF(state, depth + 1, cost + delta);

        state.part_M.pop_back();
        state.used_data_flag[v] = 0;
        state.mapped_q[u] = -1;

        if (outputLimitReached()) break;
    }

    if (!branched) stats.prun_calls++;
}

#endif // !AFEC_ENABLE_ANCHOR_FRONTIER

} // namespace afec
