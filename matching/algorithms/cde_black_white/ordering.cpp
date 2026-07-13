#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

#if !CDE_BLACK_WHITE_FIXED_ORDER && !CDE_BLACK_WHITE_RANDOM_ORDER
void MatchingSolver::initStaticEdgeSupports()
{
    // pairSupport(a,u) is the number of filtered candidate-graph edges
    // between C(a) and C(u). The data and query graphs are undirected, so one
    // traversal per undirected query edge supplies both directed lookups.
    static_candidate_count.assign(qn, 0);
    for (ui u = 0; u < qn; ++u) {
        static_candidate_count[u] = (ui)candidates[u].size();
    }
    static_edge_support.assign(qn, vector<unsigned long long>(qn, 0));

    for (ui anchor = 0; anchor < qn; ++anchor) {
        for (ui u : q_neighbors[anchor]) {
            if (anchor >= u) continue;

            unsigned long long pair_support = 0;
            for (ui anchor_candidate : candidates[anchor]) {
                const pair<size_t, ui> *range =
                    findAdjRange(anchor, anchor_candidate, u);
                if (range != nullptr) pair_support += range->second;
            }
            static_edge_support[anchor][u] = pair_support;
            static_edge_support[u][anchor] = pair_support;
        }
    }
}
#endif

#if CDE_BLACK_WHITE_FIXED_ORDER
struct MatchingSolver::FixedEdgePriorityEntry {
    ui u = 0;
    ui anchor = 0;
    unsigned long long numerator = 0;
    unsigned long long denominator = 1;
    ui u_candidate_count = 0;
    ui edge_label_frequency = 0;
};

void MatchingSolver::initFixedEdgePriorities()
{
    // Match TreeSpan's directed edge priority:
    // |C(u)| * phi(u, anchor) / phi(u), where u is the new vertex.
    const ui invalid_priority = std::numeric_limits<ui>::max();
    static_edge_priority.assign(qn, vector<ui>(qn, invalid_priority));

    vector<ui> data_label_frequency(label_count, 0);
    vector<vector<ui>> data_edge_label_frequency(
        label_count, vector<ui>(label_count, 0));

    for (ui v = 0; v < gn; ++v) {
        LabelID label = data_graph->getVertexLabel(v);
        if (label >= 0 && (ui)label < label_count) {
            data_label_frequency[(ui)label]++;
        }
    }

    for (ui v = 0; v < gn; ++v) {
        LabelID v_label = data_graph->getVertexLabel(v);
        if (v_label < 0 || (ui)v_label >= label_count) continue;

        ui degree = 0;
        const ui *neighbors = data_graph->getVertexNeighbors(v, degree);
        for (ui i = 0; i < degree; ++i) {
            ui neighbor = neighbors[i];
            if (v >= neighbor) continue;

            LabelID neighbor_label = data_graph->getVertexLabel(neighbor);
            if (neighbor_label < 0 || (ui)neighbor_label >= label_count) continue;

            data_edge_label_frequency[(ui)v_label][(ui)neighbor_label]++;
            data_edge_label_frequency[(ui)neighbor_label][(ui)v_label]++;
        }
    }

    vector<FixedEdgePriorityEntry> entries;
    size_t directed_edge_count = 0;
    for (ui u = 0; u < qn; ++u) {
        directed_edge_count += q_neighbors[u].size();
    }
    entries.reserve(directed_edge_count);

    for (ui u = 0; u < qn; ++u) {
        for (ui anchor : q_neighbors[u]) {
            LabelID u_label = query_graph->getVertexLabel(u);
            LabelID anchor_label = query_graph->getVertexLabel(anchor);

            ui vertex_frequency = 1;
            ui edge_frequency = 0;
            if (u_label >= 0 && (ui)u_label < label_count) {
                vertex_frequency = std::max(
                    (ui)1, data_label_frequency[(ui)u_label]);
                if (anchor_label >= 0 && (ui)anchor_label < label_count) {
                    edge_frequency = data_edge_label_frequency
                        [(ui)u_label][(ui)anchor_label];
                }
            }

            FixedEdgePriorityEntry entry;
            entry.u = u;
            entry.anchor = anchor;
            entry.u_candidate_count = (ui)candidates[u].size();
            entry.edge_label_frequency = edge_frequency;
            entry.numerator =
                (unsigned long long)entry.u_candidate_count * edge_frequency;
            entry.denominator = vertex_frequency;
            entries.push_back(entry);
        }
    }

    std::sort(entries.begin(), entries.end(),
        [&](const FixedEdgePriorityEntry &lhs,
            const FixedEdgePriorityEntry &rhs) {
            __uint128_t lhs_scaled = (__uint128_t)lhs.numerator * rhs.denominator;
            __uint128_t rhs_scaled = (__uint128_t)rhs.numerator * lhs.denominator;
            if (lhs_scaled != rhs_scaled) {
                return lhs_scaled < rhs_scaled;
            }
            if (lhs.u_candidate_count != rhs.u_candidate_count) {
                return lhs.u_candidate_count < rhs.u_candidate_count;
            }
            if (lhs.edge_label_frequency != rhs.edge_label_frequency) {
                return lhs.edge_label_frequency < rhs.edge_label_frequency;
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

// black-anchor edge: support(u,a) = |C(u) intersect N_G(M(a))|.
double MatchingSolver::blackSupport(const SearchState &state, ui u, ui anchor) const
{
    assert(isBlack(state, anchor));

    ui mapped_anchor = (ui)state.mapped_q[anchor];
    const pair<size_t, ui> *range = findAdjRange(anchor, mapped_anchor, u);
    if (range == nullptr) return 0.0;

    return (double)range->second;
}

// white-anchor edge: support(u,a) = |B(a)| * pairSupport(a,u) / |C(a)|.
double MatchingSolver::whiteSupport(const SearchState &state, ui u, ui anchor) const
{
    assert(isWhite(state, anchor));
    assert(u < qn && anchor < qn);

    double anchor_candidate_count = (double)static_candidate_count[anchor];
    assert(anchor_candidate_count > 0.0);
    return (double)state.white[anchor].feasible_count *
        (double)static_edge_support[anchor][u] / anchor_candidate_count;
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
#else
    // Compare support / live-anchor-count without performing division.
    double lhs_scaled = lhs.support * (double)std::max((ui)1, rhs.live_anchor_count);
    double rhs_scaled = rhs.support * (double)std::max((ui)1, lhs.live_anchor_count);
    if (lhs_scaled != rhs_scaled) {
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

#if CDE_BLACK_WHITE_RANDOM_ORDER
    // Shuffle once per search state, then retain the random Top-(b+1) prefix.
    std::shuffle(top_edges.begin(), top_edges.end(), random_order_rng);
#else
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
#if !CDE_BLACK_WHITE_FIXED_ORDER && !CDE_BLACK_WHITE_RANDOM_ORDER
        if (edge.anchor_color == COLOR_BLACK) {
            edge.support = blackSupport(state, u, anchor);
        }
        else {
            assert(edge.anchor_color == COLOR_WHITE);
            edge.support = whiteSupport(state, u, anchor);
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
