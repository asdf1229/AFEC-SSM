#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

ui MatchingSolver::chooseRoot()
{
    ui root = 0;
    for (ui u = 1; u < qn; ++u) {
        size_t cand_u = candidates[u].size();
        size_t cand_root = candidates[root].size();

        ui deg_u = q_degree[u];
        ui deg_root = q_degree[root];

        if (cand_u * deg_root < cand_root * deg_u) root = u;
    }
    return root;
}

void MatchingSolver::initColors()
{
    static_root = chooseRoot();
    static_color.assign(qn, COLOR_WHITE);
    static_color[static_root] = COLOR_BLACK;
}

ui MatchingSolver::chooseMatWhite(const SearchState &state) const
{
    // 选择候选数最少的 white 点作为优先具体化对象。
    ui chosen = qn;
    ui best_count = std::numeric_limits<ui>::max();
    for (ui u = 0; u < qn; ++u) {
        if (!isWhite(state, u)) {
            continue;
        }
        ui count = state.white[u].feasible_count;
        if (count < best_count || (count == best_count && u < chosen)) {
            chosen = u;
            best_count = count;
        }
    }
    return chosen;
}

#if CDE_BLACK_WHITE_FIXED_ORDER
struct MatchingSolver::FixedEdgePriorityEntry {
    ui u = 0;
    ui anchor = 0;
    unsigned long long pair_support = 0;
    ui u_candidate_count = 0;
    ui anchor_candidate_count = 0;
};

void MatchingSolver::initFixedEdgePriorities()
{
    // 根据静态候选边支持为查询边建立固定分支优先级。
    const ui invalid_priority = std::numeric_limits<ui>::max();
    static_edge_priority.assign(qn, vector<ui>(qn, invalid_priority));

    vector<FixedEdgePriorityEntry> entries;
    size_t directed_edge_count = 0;
    for (ui u = 0; u < qn; ++u) {
        directed_edge_count += q_neighbors[u].size();
    }
    entries.reserve(directed_edge_count);

    for (ui u = 0; u < qn; ++u) {
        for (ui anchor : q_neighbors[u]) {
            if (u >= anchor) {
                continue;
            }

            unsigned long long pair_support =
                static_edge_support[u][anchor];
            ui u_candidate_count = (ui)candidates[u].size();
            ui anchor_candidate_count = (ui)candidates[anchor].size();

            FixedEdgePriorityEntry forward;
            forward.u = u;
            forward.anchor = anchor;
            forward.pair_support = pair_support;
            forward.u_candidate_count = u_candidate_count;
            forward.anchor_candidate_count = anchor_candidate_count;
            entries.push_back(forward);

            FixedEdgePriorityEntry reverse;
            reverse.u = anchor;
            reverse.anchor = u;
            reverse.pair_support = pair_support;
            reverse.u_candidate_count = anchor_candidate_count;
            reverse.anchor_candidate_count = u_candidate_count;
            entries.push_back(reverse);
        }
    }

    std::sort(entries.begin(), entries.end(),
        [&](const FixedEdgePriorityEntry &lhs,
            const FixedEdgePriorityEntry &rhs) {
            __uint128_t lhs_scaled =
                (__uint128_t)lhs.pair_support *
                std::max((ui)1, rhs.anchor_candidate_count);
            __uint128_t rhs_scaled =
                (__uint128_t)rhs.pair_support *
                std::max((ui)1, lhs.anchor_candidate_count);
            if (lhs_scaled != rhs_scaled) {
                return lhs_scaled < rhs_scaled;
            }
            if (lhs.u_candidate_count != rhs.u_candidate_count) {
                return lhs.u_candidate_count < rhs.u_candidate_count;
            }
            if (lhs.pair_support != rhs.pair_support) {
                return lhs.pair_support < rhs.pair_support;
            }
            if (q_degree[lhs.u] != q_degree[rhs.u]) {
                return q_degree[lhs.u] > q_degree[rhs.u];
            }
            if (q_degree[lhs.anchor] != q_degree[rhs.anchor]) {
                return q_degree[lhs.anchor] > q_degree[rhs.anchor];
            }
            if (lhs.u != rhs.u) {
                return lhs.u < rhs.u;
            }
            return lhs.anchor < rhs.anchor;
        });

    assert(entries.size() <= (size_t)std::numeric_limits<ui>::max());
    for (size_t rank = 0; rank < entries.size(); ++rank) {
        const FixedEdgePriorityEntry &entry = entries[rank];
        static_edge_priority[entry.u][entry.anchor] = (ui)rank;
    }
}
#endif

void MatchingSolver::addFrontierEdgeRaw(SearchState &state, ui u,
    ui anchor) const
{
    if (u >= qn || anchor >= qn) {
        return;
    }
    size_t idx = edgeIdx(u, anchor);
    if (state.frontier_edge_pos[idx] != -1) {
        return;
    }

    ActiveEdge edge;
    edge.u = u;
    edge.anchor = anchor;
    state.frontier_edge_pos[idx] = (int)state.frontier_edges.size();
    state.frontier_edges.push_back(edge);
    state.frontier_count[u]++;
}

void MatchingSolver::removeFrontierEdgeRaw(SearchState &state, ui u,
    ui anchor) const
{
    if (u >= qn || anchor >= qn) {
        return;
    }
    size_t idx = edgeIdx(u, anchor);
    int pos = state.frontier_edge_pos[idx];
    if (pos == -1) {
        return;
    }

    size_t remove_pos = (size_t)pos;
    size_t last_pos = state.frontier_edges.size() - 1;
    if (remove_pos != last_pos) {
        ActiveEdge moved = state.frontier_edges[last_pos];
        state.frontier_edges[remove_pos] = moved;
        state.frontier_edge_pos[edgeIdx(moved.u, moved.anchor)] = pos;
    }
    state.frontier_edges.pop_back();
    state.frontier_edge_pos[idx] = -1;
    assert(state.frontier_count[u] > 0);
    state.frontier_count[u]--;
}

void MatchingSolver::refreshFrontierEdge(SearchState &state, ui u,
    ui v) const
{
    if (u >= qn || v >= qn) {
        return;
    }

    removeFrontierEdgeRaw(state, u, v);
    removeFrontierEdgeRaw(state, v, u);

    if (getEdge(state, u, v) != EDGE_UNDECIDED) {
        return;
    }

    bool u_selected = isSelected(state, u);
    bool v_selected = isSelected(state, v);
    if (u_selected == v_selected) {
        return;
    }

    if (u_selected) {
        addFrontierEdgeRaw(state, v, u);
    }
    else {
        addFrontierEdgeRaw(state, u, v);
    }
}

void MatchingSolver::refreshFrontierEdgesIncidentTo(SearchState &state,
    ui u) const
{
    if (u >= qn) {
        return;
    }
    for (ui neighbor : q_neighbors[u]) {
        refreshFrontierEdge(state, u, neighbor);
    }
}

double MatchingSolver::blackSupport(const SearchState &state, ui u,
    ui anchor) const
{
    // 估计 black anchor 对未选点 u 的剩余候选支持数量。
    if (!isBlack(state, anchor)) {
        return 0.0;
    }

    double support_count = 0.0;
    ui mapped_anchor = (ui)state.mapped_q[anchor];
#if CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SUPPORT
    const pair<size_t, ui> *range =
        findAdjRange(anchor, mapped_anchor, u);
    if (range == nullptr) {
        return 0.0;
    }

    for (const ui *it = rangeBegin(*range);
        it != rangeEnd(*range); ++it) {
        ui v = *it;
        if (isDataVertexUsed(state, v)) {
            continue;
        }
        support_count += 1.0;
    }
#else
    for (int candidate : candidates[u]) {
        ui v = (ui)candidate;
        if (isDataVertexUsed(state, v)) {
            continue;
        }
        if (data_graph->hasEdge(mapped_anchor, v)) {
            support_count += 1.0;
        }
    }
#endif
    return support_count;
}

double MatchingSolver::whiteSupport(const SearchState &state, ui anchor) const
{
    // 估计 white anchor 的分支支持，使用其候选桶可行数量。
    if (!isWhite(state, anchor)) {
        return 0.0;
    }
    return (double)std::max((ui)1, state.white[anchor].feasible_count);
}

bool MatchingSolver::betterEdge(const ActiveEdge &lhs,
    const ActiveEdge &rhs) const
{
    // 比较两条活跃边的分支优先级。
#if CDE_BLACK_WHITE_FIXED_ORDER
    ui lhs_priority = static_edge_priority[lhs.u][lhs.anchor];
    ui rhs_priority = static_edge_priority[rhs.u][rhs.anchor];
    if (lhs_priority != rhs_priority) {
        return lhs_priority < rhs_priority;
    }
    if (lhs.u != rhs.u) {
        return lhs.u < rhs.u;
    }
    return lhs.anchor < rhs.anchor;
#else
    double lhs_scaled = lhs.rank_support *
        (double)std::max((ui)1, rhs.live_anchor_count);
    double rhs_scaled = rhs.rank_support *
        (double)std::max((ui)1, lhs.live_anchor_count);
    double scale = std::max(1.0,
        std::max(std::fabs(lhs_scaled), std::fabs(rhs_scaled)));
    if (std::fabs(lhs_scaled - rhs_scaled) > 1e-12 * scale) {
        return lhs_scaled < rhs_scaled;
    }
    if (lhs.live_anchor_count != rhs.live_anchor_count) {
        return lhs.live_anchor_count > rhs.live_anchor_count;
    }
    if (lhs.query_degree != rhs.query_degree) {
        return lhs.query_degree > rhs.query_degree;
    }
    if (lhs.u != rhs.u) {
        return lhs.u < rhs.u;
    }
    return lhs.anchor < rhs.anchor;
#endif
}

void MatchingSolver::selectTopEdges(ui max_count,
    vector<ActiveEdge> &top_edges)
{
    // 从活跃边集合中按启发式选择前 max_count 条。
    size_t selected_limit = top_edges.size();
    if ((size_t)max_count < selected_limit) {
        selected_limit = (size_t)max_count;
    }

#if CDE_BLACK_WHITE_FIXED_ORDER
    auto better_edge = [&](const ActiveEdge &lhs,
        const ActiveEdge &rhs) {
        return betterEdge(lhs, rhs);
        };
    if (top_edges.size() > selected_limit) {
        partial_sort(top_edges.begin(), top_edges.begin() + selected_limit,
            top_edges.end(), better_edge);
    }
    else {
        sort(top_edges.begin(), top_edges.end(), better_edge);
    }
#elif CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
    const double gamma = (double)CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY_GAMMA;
    for (size_t selected_idx = 0; selected_idx < selected_limit; ++selected_idx) {
        size_t best_idx = selected_idx;
        for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
            if (betterEdge(top_edges[i], top_edges[best_idx])) {
                best_idx = i;
            }
        }

        if (best_idx != selected_idx) {
            std::swap(top_edges[selected_idx], top_edges[best_idx]);
        }

        const ActiveEdge &selected = top_edges[selected_idx];
        double candidate_count = (double)candidates[selected.u].size();
        if (candidate_count <= 0.0) {
            continue;
        }

        double factor = 1.0 - gamma * selected.rank_support / candidate_count;
        if (factor < 0.0) {
            factor = 0.0;
        }
        else if (factor > 1.0) {
            factor = 1.0;
        }

        if (factor == 1.0) {
            continue;
        }

        for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
            if (top_edges[i].u == selected.u) {
                top_edges[i].rank_support *= factor;
            }
        }
    }
#else
    auto better_edge = [&](const ActiveEdge &lhs,
        const ActiveEdge &rhs) {
        return betterEdge(lhs, rhs);
        };
    if (top_edges.size() > selected_limit) {
        partial_sort(top_edges.begin(), top_edges.begin() + selected_limit,
            top_edges.end(), better_edge);
    }
    else {
        sort(top_edges.begin(), top_edges.end(), better_edge);
    }
#endif
    top_edges.resize(selected_limit);
}

bool MatchingSolver::isActiveFrontierEdge(const SearchState &state,
    const ActiveEdge &edge) const
{
    ui u = edge.u;
    ui anchor = edge.anchor;
    return u < qn && anchor < qn &&
        state.color[u] == COLOR_UNSELECTED &&
        isSelected(state, anchor) &&
        getEdge(state, u, anchor) == EDGE_UNDECIDED &&
        state.frontier_count[u] > 0;
}

size_t MatchingSolver::labelFrontierComponents(const SearchState &state,
    vector<int> &component_id, vector<vector<ui>> &component_frontiers,
    vector<ui> &queue)
{
    // 从已有活跃 anchor 边的未选端出发，标记其所在未选连通分量。
    component_id.assign(qn, -1);
    for (vector<ui> &frontier : component_frontiers) {
        frontier.clear();
    }
    queue.clear();

    size_t component_count = 0;
    for (const ActiveEdge &edge : state.frontier_edges) {
        if (!isActiveFrontierEdge(state, edge)) {
            continue;
        }

        ui start = edge.u;
        if (component_id[start] != -1) {
            continue;
        }

        if (component_frontiers.size() == component_count) {
            component_frontiers.emplace_back();
        }
        vector<ui> &frontier = component_frontiers[component_count];

        int component = (int)component_count;
        component_id[start] = component;
        queue.clear();
        queue.push_back(start);

        for (size_t head = 0; head < queue.size(); ++head) {
            ui u = queue[head];
            if (state.frontier_count[u] > 0) {
                frontier.push_back(u);
            }

            for (ui neighbor : q_neighbors[u]) {
                if (state.color[neighbor] != COLOR_UNSELECTED ||
                    component_id[neighbor] != -1) {
                    continue;
                }
                component_id[neighbor] = component;
                queue.push_back(neighbor);
            }
        }

        component_count++;
    }

    return component_count;
}

void MatchingSolver::restrictTopEdgesToCoveredComponent(
    const SearchState &state, vector<ActiveEdge> &top_edges)
{
    if (top_edges.empty()) {
        return;
    }

    size_t component_count = labelFrontierComponents(state,
        component_id_buffer, component_frontiers_buffer,
        component_queue_buffer);
    if (component_count == 0) {
        return;
    }

    component_edge_counts_buffer.assign(component_count, 0);
    component_seen_counts_buffer.assign(component_count, 0);

    for (const ActiveEdge &edge : state.frontier_edges) {
        if (!isActiveFrontierEdge(state, edge)) {
            continue;
        }
        int component = component_id_buffer[edge.u];
        if (component < 0) {
            continue;
        }
        component_edge_counts_buffer[(size_t)component]++;
    }

    for (size_t edge_idx = 0; edge_idx < top_edges.size(); ++edge_idx) {
        const ActiveEdge &edge = top_edges[edge_idx];
        int component = component_id_buffer[edge.u];
        if (component < 0) {
            continue;
        }

        size_t id = (size_t)component;
        component_seen_counts_buffer[id]++;
        if (component_seen_counts_buffer[id] !=
            component_edge_counts_buffer[id]) {
            continue;
        }

        top_edges.resize(edge_idx + 1);
        return;
    }
}

bool MatchingSolver::collectActiveEdges(const SearchState &state,
    ui max_count, vector<ActiveEdge> &top_edges)
{
    // 收集当前状态下可分支的活跃边，并截取 top 边。
    top_edges.clear();
    if (max_count == 0) {
        return false;
    }

    for (const ActiveEdge &frontier_edge : state.frontier_edges) {
        ui u = frontier_edge.u;
        ui anchor = frontier_edge.anchor;
        if (!isActiveFrontierEdge(state, frontier_edge)) {
            continue;
        }

        ActiveEdge edge;
        edge.u = u;
        edge.anchor = anchor;
        edge.live_anchor_count = state.frontier_count[u];
        edge.query_degree = q_degree[u];
#if !CDE_BLACK_WHITE_FIXED_ORDER
        if (isBlack(state, anchor)) {
            edge.rank_support = blackSupport(state, u, anchor);
        }
        else {
            edge.rank_support = whiteSupport(state, anchor);
        }
#endif
        top_edges.push_back(edge);
    }

    if (top_edges.empty()) {
        return false;
    }
    selectTopEdges(max_count, top_edges);
    restrictTopEdgesToCoveredComponent(state, top_edges);
    return true;
}

} // namespace cde_black_white
