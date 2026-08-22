#ifndef MATCHING_ALGORITHMS_AFEC_COST_SPLIT_H_
#define MATCHING_ALGORITHMS_AFEC_COST_SPLIT_H_

#include "matching/algorithms/afec/afec.h"

namespace afec {

#if AFEC_ENABLE_COST_SPLIT

inline bool MatchingSolver::chooseSplitWhite(const SearchState &state, ui cost,
    ui &split_u, bool &needs_split) const
{
    assert(cost <= threshold);

    needs_split = false;
    split_u = qn;

    const ui remaining_budget = threshold - cost;
    ui best_branch_count = std::numeric_limits<ui>::max();
    ui best_feasible_count = std::numeric_limits<ui>::max();
    ui best_min_delta = 0;

    for (ui u = 0; u < qn; ++u) {
        if (!isWhite(state, u)) continue;

        WhiteCands white_bucket = state.white[u];
        assert(white_bucket.begin + white_bucket.count <= state.white_candidate_pool.size());
        assert(white_bucket.begin + white_bucket.count <= state.white_candidate_delta_pool.size());

        vector<unsigned char> seen_delta((size_t)remaining_budget + 1, 0);
        ui feasible_count = 0;
        ui branch_count = 0;
        ui min_delta = std::numeric_limits<ui>::max();

        for (ui candidate_idx = 0; candidate_idx < white_bucket.count; ++candidate_idx) {
            size_t pool_idx = white_bucket.begin + candidate_idx;
            ui candidate = state.white_candidate_pool[pool_idx];
            if (isDataVertexUsed(state, candidate)) continue;

            ui delta = state.white_candidate_delta_pool[pool_idx];
            if (delta > remaining_budget) continue;

            feasible_count++;
            if (!seen_delta[delta]) {
                seen_delta[delta] = 1;
                branch_count++;
                if (delta < min_delta) min_delta = delta;
            }
        }

        if (feasible_count == 0) return false;

        // A singleton cost bucket has no deferred choice left, so materialize
        // it immediately as black instead of carrying it as a white vertex.
        bool white_needs_split = feasible_count == 1 ||
            branch_count > 1 || min_delta > 0;
        if (!white_needs_split) continue;

        bool better = !needs_split ||
            branch_count < best_branch_count ||
            (branch_count == best_branch_count && feasible_count < best_feasible_count) ||
            (branch_count == best_branch_count && feasible_count == best_feasible_count &&
                min_delta > best_min_delta) ||
            (branch_count == best_branch_count && feasible_count == best_feasible_count &&
                min_delta == best_min_delta && q_degree[u] > q_degree[split_u]) ||
            (branch_count == best_branch_count && feasible_count == best_feasible_count &&
                min_delta == best_min_delta && q_degree[u] == q_degree[split_u] &&
                u < split_u);

        if (!better) continue;

        needs_split = true;
        split_u = u;
        best_branch_count = branch_count;
        best_feasible_count = feasible_count;
        best_min_delta = min_delta;
    }

    return true;
}

inline void MatchingSolver::branchSplitWhite(SearchState &state, ui cost,
    ui white_u, const MatchingSolver::ContinueBranch &continue_branch)
{
    assert(isWhite(state, white_u));
    assert(cost <= threshold);

    stats.split_calls++;

    const ui remaining_budget = threshold - cost;
    vector<vector<ui>> buckets((size_t)remaining_budget + 1);

    WhiteCands white_bucket = state.white[white_u];
    assert(white_bucket.begin + white_bucket.count <= state.white_candidate_pool.size());
    assert(white_bucket.begin + white_bucket.count <= state.white_candidate_delta_pool.size());

    for (ui candidate_idx = 0; candidate_idx < white_bucket.count; ++candidate_idx) {
        size_t pool_idx = white_bucket.begin + candidate_idx;
        ui candidate = state.white_candidate_pool[pool_idx];
        if (isDataVertexUsed(state, candidate)) continue;

        ui delta = state.white_candidate_delta_pool[pool_idx];
        if (delta > remaining_budget) continue;
        buckets[delta].push_back(candidate);
    }

    for (ui delta = 0; delta <= remaining_budget; ++delta) {
        const vector<ui> &bucket = buckets[delta];
        if (bucket.empty()) continue;

        size_t undo_mark = mark();
        stats.split_branches++;

        if (bucket.size() == 1) {
            ui next_cost = cost;
            if (tryMapWhite(state, cost, white_u, bucket.front(), delta,
                    next_cost)) {
                assert(next_cost == cost + delta);
                branchSplitWhites(state, next_cost, continue_branch);
            }
        }
        else {
            vector<ui> zero_deltas(bucket.size(), 0);
            replaceBucket(state, white_u, bucket, zero_deltas);
            branchSplitWhites(state, cost + delta, continue_branch);
        }

        rollback(state, undo_mark);
        if (outputLimitReached()) return;
    }
}

inline void MatchingSolver::branchSplitWhites(SearchState &state, ui cost,
    const MatchingSolver::ContinueBranch &continue_branch)
{
    if (outputLimitReached()) return;
    if (cost > threshold) {
        stats.prun_calls++;
        return;
    }

    ui split_u = qn;
    bool needs_split = false;
    if (!chooseSplitWhite(state, cost, split_u, needs_split)) {
        stats.prun_calls++;
        return;
    }

    if (!needs_split) {
        continue_branch(state, cost);
        return;
    }

    branchSplitWhite(state, cost, split_u, continue_branch);
}

inline void MatchingSolver::continueAfterSplit(SearchState &state, ui cost)
{
    branchSplitWhites(state, cost,
        [this](SearchState &split_state, ui split_cost) {
            search(split_state, split_cost);
        });
}

#endif // AFEC_ENABLE_COST_SPLIT

} // namespace afec

#endif // MATCHING_ALGORITHMS_AFEC_COST_SPLIT_H_
