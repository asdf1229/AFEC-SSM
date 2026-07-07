#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

bool MatchingSolver::tryMapBlackWithDelta(SearchState &state, ui cost, ui u, ui v, ui delta, ui &next_cost)
{
    assert(u < qn && v < gn && candidates[u].contains(v));
    assert(state.color[u] == COLOR_UNSELECTED);
    assert(!isDataVertexUsed(state, v));

    next_cost = cost + delta;
    assert(next_cost <= threshold);

    setColor(state, u, COLOR_BLACK);
    setMap(state, u, (int)v);
    pushUsed(state, v);
    pushMatch(state, u, v);
    setSelectedCnt(state, state.selected_count + 1);

    for (ui neighbor : q_neighbors[u]) {
        if (!isBlack(state, neighbor)) continue;
        EdgeState state_uv = getEdge(state, u, neighbor);
        if (state_uv == EDGE_UNDECIDED) {
            ui mapped_neighbor = (ui)state.mapped_q[neighbor];
            bool adjacent = anchorAdjacent(neighbor, mapped_neighbor, u, v);
            if (!adjacent && isQueryBridgeEdge(u, neighbor)) return false;
            setEdge(state, u, neighbor, adjacent ? EDGE_PRESENT : EDGE_MISSING);
        }
    }

    for (ui neighbor : q_neighbors[u]) {
        if (!isWhite(state, neighbor)) continue;
        if (!refreshWhiteByBlack(state, neighbor, u, v, next_cost)) {
            return false;
        }
    }

    return true;
}

bool MatchingSolver::tryMapWhite(SearchState &state, ui cost, ui white_u, ui candidate, ui &next_cost)
{
    assert(isWhite(state, white_u));
    assert(candidate < gn);
    assert(!isDataVertexUsed(state, candidate));
    assert(candidates[white_u].contains(candidate));
    assert(cost <= threshold);

    setColor(state, white_u, COLOR_BLACK);
    assert(state.white_count > 0);
    setWhiteCnt(state, state.white_count - 1);
    setMap(state, white_u, (int)candidate);
    pushUsed(state, candidate);
    pushMatch(state, white_u, candidate);

    ui delta = 0;
    for (ui neighbor : q_neighbors[white_u]) {
        if (!isSelected(state, neighbor)) continue;
        if (isWhite(state, neighbor)) return false;

        ui mapped_neighbor = (ui)state.mapped_q[neighbor];
        bool adjacent = anchorAdjacent(neighbor, mapped_neighbor, white_u, candidate);
        EdgeState state_uv = getEdge(state, white_u, neighbor);
        if (state_uv == EDGE_PRESENT) {
            if (!adjacent) return false;
        }
        else if (state_uv == EDGE_MISSING) {
            if (adjacent) return false;
        }
        else { // EDGE_UNDECIDED
            if (!adjacent) {
                if (isQueryBridgeEdge(white_u, neighbor)) return false;
                delta++;
                if (cost + delta > threshold) return false;
            }
            setEdge(state, white_u, neighbor, adjacent ? EDGE_PRESENT : EDGE_MISSING);
        }
    }

    next_cost = cost + delta;
    return true;
}

void MatchingSolver::branchWhite(SearchState &state, ui cost, ui u)
{
    assert(u < qn);
    assert(state.color[u] == COLOR_UNSELECTED);
#ifndef NDEBUG
    for (ui neighbor : q_neighbors[u]) assert(!isWhite(state, neighbor));
#endif

    if (!buildWhiteCands(state, u, cost, nullptr)) return;
    if (candidate_result_buffer.empty()) return;

    size_t undo_mark = mark();
    setColor(state, u, COLOR_WHITE);
    replaceBucket(state, u, candidate_result_buffer);
    setSelectedCnt(state, state.selected_count + 1);
    setWhiteCnt(state, state.white_count + 1);
    search(state, cost);
    rollback(state, undo_mark);
    return;
}

void MatchingSolver::branchMatWhite(SearchState &state, ui cost,
    ui white_u, const MatchingSolver::ContinueBranch &continue_branch)
{
    assert(isWhite(state, white_u));

    WhiteCands white_bucket = state.white[white_u];
    assert(white_bucket.begin + white_bucket.count <= state.white_candidate_pool.size());
    for (ui candidate_idx = 0; candidate_idx < white_bucket.count; ++candidate_idx) {
        ui candidate = state.white_candidate_pool[white_bucket.begin + candidate_idx];
        if (isDataVertexUsed(state, candidate)) continue;

        size_t undo_mark = mark();
        ui next_cost = cost;
        if (!tryMapWhite(state, cost, white_u, candidate, next_cost)) {
            rollback(state, undo_mark);
            continue;
        }
        continue_branch(state, next_cost);
        rollback(state, undo_mark);
        if (outputLimitReached()) return;
    }
}

void MatchingSolver::branchMatWhites(SearchState &state, ui cost,
    const vector<ui> &white_vertices, size_t pos,
    const MatchingSolver::ContinueBranch &continue_branch)
{
    if (pos == white_vertices.size()) {
        continue_branch(state, cost);
        return;
    }

    ui white_u = white_vertices[pos];
    assert(isWhite(state, white_u));

    branchMatWhite(state, cost, white_u,
        [&](SearchState &next_state, ui next_cost) {
            branchMatWhites(next_state, next_cost, white_vertices, pos + 1, continue_branch);
        });
}

bool MatchingSolver::shouldExpandAsWhite(const SearchState &state, ui u,
    const vector<pair<ui, ui>> &anchor_candidates) const
{
    (void)state;
    (void)anchor_candidates;
#if CDE_BLACK_WHITE_STATIC_COLOR
    assert(u < static_color.size());
    assert(static_color[u] == COLOR_WHITE || static_color[u] == COLOR_BLACK);
    return static_color[u] == COLOR_WHITE;
#else
    // TODO
    return true;
#endif
}

void MatchingSolver::branchBlack(SearchState &state, ui cost, ui u,
    const vector<pair<ui, ui>> &anchor_candidates)
{
    auto map_candidate = [&](const pair<ui, ui> &candidate_delta) -> bool {
        size_t undo_mark = mark();

        ui next_cost = cost;
        if (!tryMapBlackWithDelta(state, cost, u, candidate_delta.first, candidate_delta.second, next_cost)) {
            rollback(state, undo_mark);
            return false;
        }

        search(state, next_cost);
        rollback(state, undo_mark);
        if (outputLimitReached()) return true;
        return false;
    };

    for (const auto &candidate : anchor_candidates) {
        if (map_candidate(candidate)) return;
    }
}

void MatchingSolver::branchBlackAnchor(SearchState &state, ui cost, ui u, ui anchor)
{
    assert(u < qn && anchor < qn);
    assert(isBlack(state, anchor));
    assert(state.color[u] == COLOR_UNSELECTED);
    assert(getEdge(state, u, anchor) == EDGE_PRESENT);

    // calculate candidates supported by the anchor
    ui mapped_anchor = (ui)state.mapped_q[anchor];
    const pair<size_t, ui> *range = findAdjRange(anchor, mapped_anchor, u);
    if (range == nullptr) return;
    vector<pair<ui, ui>> &anchor_candidates = branchCandsBuffer(state.selected_count);
    if (anchor_candidates.capacity() < range->second) anchor_candidates.reserve(range->second);
    for (const ui *it = rangeBegin(*range); it != rangeEnd(*range); ++it) {
        ui candidate = *it;
        if (isDataVertexUsed(state, candidate)) continue;
        ui delta = 0;
        if (!calcBlackDelta(state, u, candidate, cost, delta)) continue;
        anchor_candidates.push_back({ candidate, delta });
    }
    if (anchor_candidates.empty()) return;

    // Decide whether to delay u as white
    bool expand_as_white = shouldExpandAsWhite(state, u, anchor_candidates);

    if(expand_as_white) {
        vector<ui> &white_neighbors = whiteNbrsBuffer(state.selected_count);
        for (ui neighbor : q_neighbors[u]) {
            if (isWhite(state, neighbor)) white_neighbors.push_back(neighbor);
        }
        if (white_neighbors.empty()) branchWhite(state, cost, u);
        else {
            branchMatWhites(state, cost, white_neighbors, 0,
                [&](SearchState &materialized_state, ui materialized_cost) {
                    branchWhite(materialized_state, materialized_cost, u);
                });
        }
        return;
    }

    // Map u directly as black
    branchBlack(state, cost, u, anchor_candidates);
}

// branch all top_edges
void MatchingSolver::branchEdges(SearchState &state, ui cost, const vector<AnchorEdge> &top_edges, size_t edge_idx)
{
    if(outputLimitReached()) return;
    if(edge_idx >= top_edges.size()) return;
    assert(cost <= threshold);

    const AnchorEdge &edge = top_edges[edge_idx];
    ui u = edge.u;
    ui anchor = edge.anchor;
    assert(u < qn && anchor < qn);
    assert(state.color[u] == COLOR_UNSELECTED);
    assert(isSelected(state, anchor));
    assert(getEdge(state, u, anchor) == EDGE_UNDECIDED);

    // branch 1: include edge(u, anchor).
    size_t undo_mark = mark();
    setEdge(state, u, anchor, EDGE_PRESENT);
    if(isBlack(state, anchor)) {
        // black anchor
        branchBlackAnchor(state, cost, u, anchor);
    }
    else {
        // white anchor
        // matrialize all white backward neighbors of u, then branch u
        branchMatWhite(state, cost, anchor,
            [&](SearchState &materialized_state, ui materialized_cost) {
                branchBlackAnchor(materialized_state, materialized_cost, u, anchor);
            });
    }
    rollback(state, undo_mark);
    assert(state.color[u] == COLOR_UNSELECTED);
    assert(isSelected(state, anchor));
    assert(getEdge(state, u, anchor) == EDGE_UNDECIDED);

    if (outputLimitReached()) return;
    if (isQueryBridgeEdge(u, anchor)) return;

    // branch 2: exclude edge(u, anchor).
    if (cost + 1 > threshold) {
        stats.prun_calls++;
        return;
    }
    undo_mark = mark();
    setEdge(state, u, anchor, EDGE_MISSING);
    branchEdges(state, cost + 1, top_edges, edge_idx + 1);
    rollback(state, undo_mark);
}

void MatchingSolver::search(SearchState &state, ui cost)
{
    if (outputLimitReached()) return;

    if (cost > threshold) {
        stats.prun_calls++;
        return;
    }

    stats.recursion_calls++;

    if (state.selected_count == qn) {
        if (state.white_count == 0) {
            emitResult(state);
            return;
        }

#if CDE_BLACK_WHITE_USE_TERMINAL_TAIL
        vector<TailWhite> tail_vertices;
        if (!buildTailWhites(state, cost, tail_vertices)) {
            stats.prun_calls++;
            return;
        }
        enumTailWhites(state, 0, cost, tail_vertices);
#else
        vector<ui> &white_vertices = whiteNbrsBuffer(state.selected_count);
        for (ui u = 0; u < qn; ++u) {
            if (isWhite(state, u)) white_vertices.push_back(u);
        }
        assert(white_vertices.size() == (size_t)state.white_count);
        branchMatWhites(state, cost, white_vertices, 0,
            [&](SearchState &materialized_state, ui materialized_cost) {
                (void)materialized_cost;
                assert(materialized_state.white_count == 0);
                emitResult(materialized_state);
            });
#endif
        return;
    }

    vector<AnchorEdge> &top_edges = topEdgesBuffer(state.selected_count);
    if (!collectActiveEdges(state, threshold - cost + 1, top_edges)) {
        stats.prun_calls++;
        return;
    }

    branchEdges(state, cost, top_edges, 0);
}

} // namespace cde_black_white
