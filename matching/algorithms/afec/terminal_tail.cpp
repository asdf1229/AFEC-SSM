#include "matching/algorithms/afec/afec.h"

namespace afec {

#if AFEC_ENABLE_COST_SPLIT

void MatchingSolver::enumTailWhites(SearchState &state, size_t pos, ui cost,
    const vector<ui> &tail_vertices)
{
    if (outputLimitReached()) return;
    stats.recursion_calls++;
    stats.terminal_tail_calls++;
    assert(cost <= threshold);

    if (pos == tail_vertices.size()) {
        emitResult(state);
        return;
    }

    ui u = tail_vertices[pos];
    assert(isWhite(state, u));
    assert(state.mapped_q[u] == -1);

    const WhiteCands &white_bucket = state.white[u];
    assert(white_bucket.begin + white_bucket.count <= state.white_candidate_pool.size());
    assert(white_bucket.begin + white_bucket.count <= state.white_candidate_delta_pool.size());

    for (ui candidate_idx = 0; candidate_idx < white_bucket.count; ++candidate_idx) {
        size_t pool_idx = white_bucket.begin + candidate_idx;
        ui v = state.white_candidate_pool[pool_idx];
        if (isDataVertexUsed(state, v)) continue;

        ui candidate_delta = state.white_candidate_delta_pool[pool_idx];
        if (candidate_delta > threshold - cost) continue;
        assert(candidate_delta == 0);

        state.mapped_q[u] = (int)v;
        state.used_data_vertices.push_back(v);
        state.used_data_flag[v] = 1;
        state.part_M.push_back({ u, v });

        enumTailWhites(state, pos + 1, cost + candidate_delta, tail_vertices);

        state.part_M.pop_back();
        state.used_data_flag[v] = 0;
        state.used_data_vertices.pop_back();
        state.mapped_q[u] = -1;

        if (outputLimitReached()) return;
    }
}

#else

bool MatchingSolver::buildTerminalTail(const SearchState &state, ui cost,
    vector<TerminalTailVertex> &tail_vertices) const
{
    assert(cost <= threshold);

    tail_vertices.clear();
    tail_vertices.reserve(state.white_count);
    const ui remaining_budget = threshold - cost;

    for (ui u = 0; u < qn; ++u) {
        if (!isWhite(state, u)) continue;
#ifndef NDEBUG
        for (ui neighbor : q_neighbors[u]) assert(!isWhite(state, neighbor));
#endif

        TerminalTailVertex tail_vertex;
        tail_vertex.u = u;
        tail_vertex.candidates_by_delta.resize((size_t)remaining_budget + 1);

        const WhiteCands &white_bucket = state.white[u];
        assert(white_bucket.begin + white_bucket.count <= state.white_candidate_pool.size());
        assert(white_bucket.begin + white_bucket.count <= state.white_candidate_delta_pool.size());

        for (ui candidate_idx = 0; candidate_idx < white_bucket.count; ++candidate_idx) {
            size_t pool_idx = white_bucket.begin + candidate_idx;
            ui candidate = state.white_candidate_pool[pool_idx];
            if (isDataVertexUsed(state, candidate)) continue;

            ui candidate_delta = state.white_candidate_delta_pool[pool_idx];
            if (candidate_delta > remaining_budget) continue;

            tail_vertex.candidates_by_delta[candidate_delta].push_back(candidate);
            tail_vertex.feasible_count++;
            tail_vertex.min_delta = std::min(tail_vertex.min_delta, candidate_delta);
        }

        if (tail_vertex.feasible_count == 0) return false;
        tail_vertices.push_back(std::move(tail_vertex));
    }

    assert(tail_vertices.size() == (size_t)state.white_count);
    std::sort(tail_vertices.begin(), tail_vertices.end(),
        [this](const TerminalTailVertex &lhs, const TerminalTailVertex &rhs) {
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

void MatchingSolver::enumTailWhites(SearchState &state, size_t pos, ui cost,
    const vector<TerminalTailVertex> &tail_vertices)
{
    if (outputLimitReached()) return;
    stats.recursion_calls++;
    stats.terminal_tail_calls++;
    assert(cost <= threshold);

    if (pos == tail_vertices.size()) {
        emitResult(state);
        return;
    }

    const TerminalTailVertex &tail_vertex = tail_vertices[pos];
    ui u = tail_vertex.u;
    assert(isWhite(state, u));
    assert(state.mapped_q[u] == -1);

    ui remaining_budget = threshold - cost;
    ui max_delta = std::min(remaining_budget,
        (ui)tail_vertex.candidates_by_delta.size() - 1);
    for (ui candidate_delta = 0; candidate_delta <= max_delta; ++candidate_delta) {
        for (ui v : tail_vertex.candidates_by_delta[candidate_delta]) {
            if (isDataVertexUsed(state, v)) continue;

            state.mapped_q[u] = (int)v;
            state.used_data_vertices.push_back(v);
            state.used_data_flag[v] = 1;
            state.part_M.push_back({ u, v });

            enumTailWhites(state, pos + 1, cost + candidate_delta, tail_vertices);

            state.part_M.pop_back();
            state.used_data_flag[v] = 0;
            state.used_data_vertices.pop_back();
            state.mapped_q[u] = -1;

            if (outputLimitReached()) return;
        }
    }
}

#endif

} // namespace afec
