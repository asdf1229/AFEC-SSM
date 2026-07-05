#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_ORDERING_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_ORDERING_H_

#ifndef CDE_BLACK_WHITE_INSIDE_MATCHING_SOLVER
#error "This internal header must be included from cde_black_white/context.h inside MatchingSolver."
#endif

#if CDE_BLACK_WHITE_FIXED_ORDER
    struct FixedEdgePriorityEntry {
        ui u = 0;
        ui anchor = 0;
        unsigned long long pair_support = 0;
        ui u_candidate_count = 0;
        ui anchor_candidate_count = 0;
    };

    void initFixedEdgePriorities()
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

    ui chooseRoot()
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

    void initColors()
    {
        static_root = chooseRoot();
        static_color.assign(qn, COLOR_WHITE);
        static_color[static_root] = COLOR_BLACK;
    }

    double blackSupport(const SearchState &state, ui u,
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

    double whiteSupport(const SearchState &state, ui anchor) const
    {
        // 估计 white anchor 的分支支持，使用其候选桶可行数量。
        if (!isWhite(state, anchor)) {
            return 0.0;
        }
        return (double)std::max((ui)1, state.white[anchor].feasible_count);
    }

    bool betterEdge(const ActiveEdge &lhs,
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

    void selectTopEdges(ui max_count,
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

    bool collectActiveEdges(const SearchState &state,
        ui max_count, vector<ActiveEdge> &top_edges)
    {
        // 收集当前状态下可分支的活跃边，并截取 top 边。
        top_edges.clear();
        if (max_count == 0) {
            return false;
        }

        for (ui u = 0; u < qn; ++u) {
            if (state.color[u] != COLOR_UNSELECTED) {
                continue;
            }

            ui live_anchor_count = 0;
            for (ui anchor : q_neighbors[u]) {
                if (isSelected(state, anchor) &&
                    getEdge(state, u, anchor) == EDGE_UNDECIDED) {
                    live_anchor_count++;
                }
            }
            if (live_anchor_count == 0) {
                continue;
            }

            for (ui anchor : q_neighbors[u]) {
                if (!isSelected(state, anchor) ||
                    getEdge(state, u, anchor) != EDGE_UNDECIDED) {
                    continue;
                }

                ActiveEdge edge;
                edge.u = u;
                edge.anchor = anchor;
                edge.live_anchor_count = live_anchor_count;
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
        }

        if (top_edges.empty()) {
            return false;
        }
        selectTopEdges(max_count, top_edges);
        return true;
    }

    ui chooseMatWhite(const SearchState &state) const
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


#endif
