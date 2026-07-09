#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

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
#if CDE_BLACK_WHITE_ENABLE_SPLIT
        assert(candidate_delta == 0);
#endif

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

} // namespace cde_black_white
