#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

bool MatchingSolver::buildTailWhite(const SearchState &state, ui cost, TailWhite &tail_white)
{
    ui white_u = tail_white.u;

    assert(cost <= threshold);
    assert(isWhite(state, white_u));
    assert(state.mapped_q[white_u] == -1);
#ifndef NDEBUG
    for (ui neighbor : q_neighbors[white_u]) assert(!isWhite(state, neighbor));
#endif

    ui remaining_budget = threshold - cost;
    tail_white.buckets.assign((size_t)remaining_budget + 1, vector<ui>());
    tail_white.feasible_count = 0;
    tail_white.min_delta = std::numeric_limits<ui>::max();

    WhiteCands white = state.white[white_u];
    assert(white.begin + white.count <= state.white_candidate_pool.size());
    assert(white.begin + white.count <= state.white_candidate_delta_pool.size());

    for (ui candidate_idx = 0; candidate_idx < white.count; ++candidate_idx) {
        size_t pool_idx = white.begin + candidate_idx;
        ui candidate = state.white_candidate_pool[pool_idx];
        if (isDataVertexUsed(state, candidate)) continue;
        ui delta = state.white_candidate_delta_pool[pool_idx];
        if (delta > remaining_budget) continue;

        tail_white.buckets[delta].push_back(candidate);
        tail_white.feasible_count++;

        if (delta < tail_white.min_delta) tail_white.min_delta = delta;
    }

    return tail_white.feasible_count > 0;
}

bool MatchingSolver::buildTailWhites(const SearchState &state, ui cost, vector<TailWhite> &tail_whites)
{
    tail_whites.clear();
    tail_whites.reserve(state.white_count);

    for (ui u = 0; u < qn; ++u) {
        if (!isWhite(state, u)) continue;

        TailWhite tail_white;
        tail_white.u = u;
        if (!buildTailWhite(state, cost, tail_white)) return false;
        tail_whites.push_back(std::move(tail_white));
    }

    assert(tail_whites.size() == (size_t)state.white_count);

    // 1. fewer feasible candidates first;
    // 2. larger minimum delta first;
    // 3. higher query degree first;
    // 4. smaller vertex id first.
    std::sort(tail_whites.begin(), tail_whites.end(),
        [this](const TailWhite &lhs, const TailWhite &rhs) {
            if (lhs.feasible_count != rhs.feasible_count) {
                return lhs.feasible_count < rhs.feasible_count;
            }
            if (lhs.min_delta != rhs.min_delta) {
                return lhs.min_delta > rhs.min_delta;
            }
            if (q_degree[lhs.u] != q_degree[rhs.u]) {
                return q_degree[lhs.u] > q_degree[rhs.u];
            }
            return lhs.u < rhs.u;
        });

    return true;
}

void MatchingSolver::enumTailWhites(SearchState &state, size_t pos, ui cost, vector<TailWhite> &tail_vertices)
{
    if (outputLimitReached()) return;
    stats.recursion_calls++;
    stats.terminal_tail_calls++;
    assert(cost <= threshold);

    if (pos == tail_vertices.size()) {
        emitResult(state);
        return;
    }

    TailWhite &tail_vertex = tail_vertices[pos];
    ui u = tail_vertex.u;
    assert(isWhite(state, u));
    assert(state.mapped_q[u] == -1);

    ui remaining_budget = threshold - cost;
    ui max_delta = std::min((ui)tail_vertex.buckets.size() - 1, remaining_budget);
    for (ui missing_delta = 0; missing_delta <= max_delta; ++missing_delta) {
        const vector<ui> &bucket = tail_vertex.buckets[missing_delta];
        for (ui v : bucket) {
            if (isDataVertexUsed(state, v)) continue;

            state.mapped_q[u] = (int)v;
            state.used_data_vertices.push_back(v);
            state.used_data_flag[v] = 1;
            state.part_M.push_back({ u, v });

            enumTailWhites(state, pos + 1, cost + missing_delta, tail_vertices);

            state.part_M.pop_back();
            state.used_data_flag[v] = 0;
            state.used_data_vertices.pop_back();
            state.mapped_q[u] = -1;

            if (outputLimitReached()) return;
        }
    }
}

} // namespace cde_black_white
