#include "matching/algorithms/afec/afec.h"

namespace afec {

void MatchingSolver::initState(SearchState &state) const
{
    state.mapped_q.assign(qn, -1);
    state.used_data_vertices.clear();
    state.used_data_vertices.reserve(qn);
    state.used_data_flag.assign(gn, 0);
    state.color.assign(qn, COLOR_UNSELECTED);
    state.edge_state.assign((size_t)qn * qn, EDGE_UNDECIDED);
    state.anchor_edges.clear();
    state.anchor_edge_pos.assign((size_t)qn * qn, -1);
    state.anchor_count.assign(qn, 0);
    state.white.clear();
    state.white.resize(qn);
    state.white_candidate_pool.clear();
    state.white_candidate_pool.reserve(stats.filter_candidate_count);
    state.white_candidate_delta_pool.clear();
    state.white_candidate_delta_pool.reserve(stats.filter_candidate_count);
    state.part_M.clear();
    state.part_M.reserve(qn);
    state.selected_count = 0;
    state.white_count = 0;
}

size_t MatchingSolver::edgeIdx(ui u, ui v) const
{
    return (size_t)u * qn + v;
}

EdgeState MatchingSolver::getEdge(const SearchState &state, ui u, ui v) const
{
    return state.edge_state[edgeIdx(u, v)];
}

void MatchingSolver::setEdgeRaw(SearchState &state, ui u, ui v, EdgeState edge_state_value) const
{
    state.edge_state[edgeIdx(u, v)] = edge_state_value;
}

void MatchingSolver::addAnchorEdgeRaw(SearchState &state, ui u, ui anchor) const
{
    assert(u < qn && anchor < qn);
    size_t idx = edgeIdx(u, anchor);
    if (state.anchor_edge_pos[idx] != -1) return;

    AnchorEdge edge;
    edge.u = u;
    edge.anchor = anchor;
    state.anchor_edge_pos[idx] = (int)state.anchor_edges.size();
    state.anchor_edges.push_back(edge);
    state.anchor_count[u]++;
}

void MatchingSolver::removeAnchorEdgeRaw(SearchState &state, ui u, ui anchor) const
{
    assert(u < qn && anchor < qn);
    size_t idx = edgeIdx(u, anchor);
    int pos = state.anchor_edge_pos[idx];
    if (pos == -1) return;

    size_t remove_pos = (size_t)pos;
    size_t last_pos = state.anchor_edges.size() - 1;
    if (remove_pos != last_pos) {
        AnchorEdge moved = state.anchor_edges[last_pos];
        state.anchor_edges[remove_pos] = moved;
        state.anchor_edge_pos[edgeIdx(moved.u, moved.anchor)] = pos;
    }
    state.anchor_edges.pop_back();
    state.anchor_edge_pos[idx] = -1;
    assert(state.anchor_count[u] > 0);
    state.anchor_count[u]--;
}

void MatchingSolver::refreshAnchorEdge(SearchState &state, ui u, ui v) const
{
    assert(u < qn && v < qn);

    removeAnchorEdgeRaw(state, u, v);
    removeAnchorEdgeRaw(state, v, u);

    if (getEdge(state, u, v) != EDGE_UNDECIDED) return;

    bool u_selected = isSelected(state, u);
    bool v_selected = isSelected(state, v);
    if (u_selected == v_selected) return;

    if (u_selected) addAnchorEdgeRaw(state, v, u);
    else addAnchorEdgeRaw(state, u, v);
}

bool MatchingSolver::isDataVertexUsed(const SearchState &state, ui v) const
{
    assert(v < gn);
    return state.used_data_flag[v] != 0;
}

bool MatchingSolver::isSelected(const SearchState &state, ui u) const
{
    assert(u < qn);
    return state.color[u] != COLOR_UNSELECTED;
}

bool MatchingSolver::isBlack(const SearchState &state, ui u) const
{
    assert(u < qn);
    return state.color[u] == COLOR_BLACK;
}

bool MatchingSolver::isWhite(const SearchState &state, ui u) const
{
    assert(u < qn);
    return state.color[u] == COLOR_WHITE;
}

size_t MatchingSolver::mark() const
{
    return undo_stack.size();
}

void MatchingSolver::rollback(SearchState &state, size_t mark)
{
    while (undo_stack.size() > mark) {
        UndoRecord undo = std::move(undo_stack.back());
        undo_stack.pop_back();

        switch (undo.kind) {
        case UNDO_MAPPED_Q:
            state.mapped_q[undo.u] = undo.old_mapped_q;
            break;
        case UNDO_COLOR:
            state.color[undo.u] = undo.old_color;
            for (ui nbr : q_neighbors[undo.u]) refreshAnchorEdge(state, undo.u, nbr);
            break;
        case UNDO_EDGE_STATE:
            setEdgeRaw(state, undo.u, undo.v, undo.old_edge_uv);
            setEdgeRaw(state, undo.v, undo.u, undo.old_edge_vu);
            refreshAnchorEdge(state, undo.u, undo.v);
            break;
        case UNDO_USED_DATA_SIZE:
            for (size_t i = undo.old_size; i < state.used_data_vertices.size(); ++i) {
                ui used_v = state.used_data_vertices[i];
                if (used_v < state.used_data_flag.size()) {
                    state.used_data_flag[used_v] = 0;
                }
            }
            state.used_data_vertices.resize(undo.old_size);
            break;
        case UNDO_PART_M_SIZE:
            state.part_M.resize(undo.old_size);
            break;
        case UNDO_SELECTED_COUNT:
            state.selected_count = undo.old_count;
            break;
        case UNDO_WHITE_COUNT:
            state.white_count = undo.old_count;
            break;
        case UNDO_WHITE_BUCKET:
            state.white_candidate_pool.resize(undo.old_size);
            state.white_candidate_delta_pool.resize(undo.old_size);
            state.white[undo.u] = undo.old_white;
            break;
        }
    }
}

void MatchingSolver::setMap(SearchState &state, ui u, int value)
{
    UndoRecord undo;
    undo.kind = UNDO_MAPPED_Q;
    undo.u = u;
    undo.old_mapped_q = state.mapped_q[u];
    undo_stack.push_back(std::move(undo));
    state.mapped_q[u] = value;
}

void MatchingSolver::setColor(SearchState &state, ui u, VertexColor value)
{
    UndoRecord undo;
    undo.kind = UNDO_COLOR;
    undo.u = u;
    undo.old_color = state.color[u];
    undo_stack.push_back(std::move(undo));
    state.color[u] = value;
    for (ui nbr : q_neighbors[u]) refreshAnchorEdge(state, u, nbr);
}

void MatchingSolver::setEdge(SearchState &state, ui u, ui v, EdgeState edge_state_value)
{
    UndoRecord undo;
    undo.kind = UNDO_EDGE_STATE;
    undo.u = u;
    undo.v = v;
    undo.old_edge_uv = getEdge(state, u, v);
    undo.old_edge_vu = getEdge(state, v, u);
    undo_stack.push_back(std::move(undo));
    setEdgeRaw(state, u, v, edge_state_value);
    setEdgeRaw(state, v, u, edge_state_value);
    refreshAnchorEdge(state, u, v);
}

void MatchingSolver::pushUsed(SearchState &state, ui v)
{
    UndoRecord undo;
    undo.kind = UNDO_USED_DATA_SIZE;
    undo.old_size = state.used_data_vertices.size();
    undo_stack.push_back(std::move(undo));
    state.used_data_vertices.push_back(v);
    state.used_data_flag[v] = 1;
}

void MatchingSolver::pushMatch(SearchState &state, ui u, ui v)
{
    UndoRecord undo;
    undo.kind = UNDO_PART_M_SIZE;
    undo.old_size = state.part_M.size();
    undo_stack.push_back(std::move(undo));
    state.part_M.push_back({ u, v });
}

void MatchingSolver::setSelectedCnt(SearchState &state, ui value)
{
    UndoRecord undo;
    undo.kind = UNDO_SELECTED_COUNT;
    undo.old_count = state.selected_count;
    undo_stack.push_back(std::move(undo));
    state.selected_count = value;
}

void MatchingSolver::setWhiteCnt(SearchState &state, ui value)
{
    UndoRecord undo;
    undo.kind = UNDO_WHITE_COUNT;
    undo.old_count = state.white_count;
    undo_stack.push_back(std::move(undo));
    state.white_count = value;
}

void MatchingSolver::replaceBucket(SearchState &state, ui u,
    const vector<ui> &candidates_to_store,
    const vector<ui> &candidate_deltas_to_store)
{
    assert(candidates_to_store.size() == candidate_deltas_to_store.size());
    assert(state.white_candidate_pool.size() ==
        state.white_candidate_delta_pool.size());

    UndoRecord undo;
    undo.kind = UNDO_WHITE_BUCKET;
    undo.u = u;
    undo.old_size = state.white_candidate_pool.size();
    undo.old_white = state.white[u];
    undo_stack.push_back(std::move(undo));

    WhiteCands bucket;
    bucket.begin = state.white_candidate_pool.size();
    bucket.count = (ui)candidates_to_store.size();
    bucket.feasible_count = bucket.count;
    state.white_candidate_pool.insert(state.white_candidate_pool.end(), candidates_to_store.begin(), candidates_to_store.end());
    state.white_candidate_delta_pool.insert(state.white_candidate_delta_pool.end(),
                                            candidate_deltas_to_store.begin(), candidate_deltas_to_store.end());
    state.white[u] = bucket;
}

} // namespace afec
