#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

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
    // Rank every directed query edge once from the filtered candidate graph.
    // Search states reuse this rank; only edge availability changes at runtime.
    const ui invalid_priority = std::numeric_limits<ui>::max();
    static_edge_support.assign(qn, vector<unsigned long long>(qn, 0));
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

            unsigned long long pair_support = 0;
            for (ui v : candidates[u]) {
                const pair<size_t, ui> *range = findAdjRange(u, v, anchor);
                if (range != nullptr) {
                    pair_support += range->second;
                }
            }
            static_edge_support[u][anchor] = pair_support;
            static_edge_support[anchor][u] = pair_support;

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

// black-anchor edge: support(u,a) = min(|C(u)|, deg_label(u)(v_a)).
double MatchingSolver::blackSupport(const SearchState &state, ui u, ui anchor) const
{
    assert(isBlack(state, anchor));

    ui mapped_anchor = (ui)state.mapped_q[anchor];
    const pair<size_t, ui> *range = findAdjRange(anchor, mapped_anchor, u);
    if (range == nullptr) return 0.0;

    return (double)range->second;
}

// white-anchor edge: support(u,a) = |C(a)|, i.e., the materialization branch count of the white anchor.
double MatchingSolver::whiteSupport(const SearchState &state, ui anchor) const
{
    assert(isWhite(state, anchor));
    return (double)state.white[anchor].feasible_count;
}

bool MatchingSolver::betterEdge(const AnchorEdge &lhs, const AnchorEdge &rhs) const
{
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
#else // CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
    double lhs_scaled = lhs.support * (double)std::max((ui)1, rhs.live_anchor_count);
    double rhs_scaled = rhs.support * (double)std::max((ui)1, lhs.live_anchor_count);
    double scale = std::max(1.0, std::max(std::fabs(lhs_scaled), std::fabs(rhs_scaled)));
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

void MatchingSolver::selectTopEdges(ui max_count, vector<AnchorEdge> &top_edges)
{
    size_t selected_limit = min((size_t)max_count, top_edges.size());

#if CDE_BLACK_WHITE_FIXED_ORDER
    auto better_edge = [&](const AnchorEdge &lhs,
        const AnchorEdge &rhs) {
        return betterEdge(lhs, rhs);
        };
    if (top_edges.size() > selected_limit) {
        partial_sort(top_edges.begin(), top_edges.begin() + selected_limit, top_edges.end(), better_edge);
    }
    else {
        sort(top_edges.begin(), top_edges.end(), better_edge);
    }
#else // CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
    const double gamma = (double)CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY_GAMMA;
    assert(gamma >= 0.0 && gamma <= 1.0);
    for (size_t selected_idx = 0; selected_idx < selected_limit; ++selected_idx) {
        size_t best_idx = selected_idx;
        for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
            if (betterEdge(top_edges[i], top_edges[best_idx])) best_idx = i;
        }
        if (best_idx != selected_idx) std::swap(top_edges[selected_idx], top_edges[best_idx]);

        const AnchorEdge &selected = top_edges[selected_idx];
        if (selected.anchor_color == COLOR_WHITE) continue;

        double candidate_count = (double)candidates[selected.u].size();
        assert(candidate_count > 0.0 && selected.support >= 0.0);
        assert(selected.support <= candidate_count);
        double factor = 1.0 - gamma * selected.support / candidate_count;
        assert(factor >= 0.0 && factor <= 1.0);

        for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
            if (top_edges[i].u == selected.u) {
                assert(top_edges[i].support >= 0.0);
                top_edges[i].support *= factor;
            }
        }
    }
#endif
    top_edges.resize(selected_limit);
}

bool MatchingSolver::isActiveAnchorEdge(const SearchState &state, const AnchorEdge &edge) const
{
    ui u = edge.u;
    ui anchor = edge.anchor;
    return u < qn && anchor < qn &&
        state.color[u] == COLOR_UNSELECTED &&
        isSelected(state, anchor) &&
        getEdge(state, u, anchor) == EDGE_UNDECIDED &&
        state.anchor_count[u] > 0;
}

size_t MatchingSolver::labelFrontierComponents(const SearchState &state, vector<int> &component_id, vector<ui> &queue)
{
    component_id.assign(qn, -1);
    queue.clear();

    size_t component_count = 0;
    for (const AnchorEdge &edge : state.anchor_edges) {
        assert(isActiveAnchorEdge(state, edge));
        ui start = edge.u;
        if (component_id[start] != -1) continue;

        int component = (int)component_count;
        component_id[start] = component;
        queue.clear();
        queue.push_back(start);

        for (size_t head = 0; head < queue.size(); ++head) {
            ui u = queue[head];
            for (ui neighbor : q_neighbors[u]) {
                if (state.color[neighbor] != COLOR_UNSELECTED) continue;
                if (component_id[neighbor] != -1) continue;
                component_id[neighbor] = component;
                queue.push_back(neighbor);
            }
        }

        component_count++;
    }

    return component_count;
}

void MatchingSolver::trimToCompleteComponent(const SearchState &state, vector<AnchorEdge> &top_edges)
{
    assert(!top_edges.empty());

    size_t component_count = labelFrontierComponents(state, component_id_buffer, component_queue_buffer);
    assert(component_count > 0);
    component_edge_counts_buffer.assign(component_count, 0);

    for (const AnchorEdge &edge : state.anchor_edges) {
        assert(isActiveAnchorEdge(state, edge));
        int component = component_id_buffer[edge.u];
        assert(component >= 0);
        component_edge_counts_buffer[(size_t)component]++;
    }

    for (size_t edge_idx = 0; edge_idx < top_edges.size(); ++edge_idx) {
        const AnchorEdge &edge = top_edges[edge_idx];
        int component = component_id_buffer[edge.u];
        assert(component >= 0);

        size_t id = (size_t)component;
        component_edge_counts_buffer[id]--;
        if (component_edge_counts_buffer[id] != 0) continue;

        top_edges.resize(edge_idx + 1);
        return;
    }
}

bool MatchingSolver::collectActiveEdges(const SearchState &state, ui max_count, vector<AnchorEdge> &top_edges)
{
    top_edges.clear();
    assert(max_count > 0);

    for (const AnchorEdge &anchor_edge : state.anchor_edges) {
        assert(isActiveAnchorEdge(state, anchor_edge));

        ui u = anchor_edge.u;
        ui anchor = anchor_edge.anchor;

        AnchorEdge edge;
        edge.u = u;
        edge.anchor = anchor;
        edge.anchor_color = state.color[anchor];
        edge.live_anchor_count = state.anchor_count[u];
        edge.query_degree = q_degree[u];
#if !CDE_BLACK_WHITE_FIXED_ORDER
        if (edge.anchor_color == COLOR_BLACK) {
            edge.support = blackSupport(state, u, anchor);
        }
        else {
            assert(edge.anchor_color == COLOR_WHITE);
            edge.support = whiteSupport(state, anchor);
        }
#endif
        top_edges.push_back(edge);
    }

    if (top_edges.empty()) return false;

    selectTopEdges(max_count, top_edges);
    trimToCompleteComponent(state, top_edges);
    return true;
}

} // namespace cde_black_white
