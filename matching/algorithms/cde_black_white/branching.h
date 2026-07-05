#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_BRANCHING_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_BRANCHING_H_

#ifndef CDE_BLACK_WHITE_INSIDE_WORKSPACE
#error "This internal header must be included from cde_black_white/context.h inside Workspace."
#endif

    bool tryBindBlack(SearchState &state, ui cost, ui u, ui v,
        ui &next_cost)
    {
        // 尝试将未选查询点 u 绑定为 black 映射到 v。
        if (u >= qn || v >= gn || !candidates[u].contains(v)) {
            return false;
        }
        if (state.color[u] != COLOR_UNSELECTED || isDataVertexUsed(state, v)) {
            return false;
        }

        ui delta = 0;
        if (!calcBlackDelta(state, u, v, cost, delta)) {
            return false;
        }
        next_cost = cost + delta;
        if (next_cost > threshold) {
            return false;
        }

        setColor(state, u, COLOR_BLACK);
        setMap(state, u, (int)v);
        pushUsed(state, v);
        pushMatch(state, u, v);
        setSelectedCnt(state, state.selected_count + 1);

        for (ui neighbor : q_neighbors[u]) {
            if (!isBlack(state, neighbor)) {
                continue;
            }
            EdgeState state_uv = getEdge(state, u, neighbor);
            if (state_uv != EDGE_UNDECIDED) {
                continue;
            }
            bool adjacent = anchorAdjacent(
                neighbor, (ui)state.mapped_q[neighbor], u, v);
            setEdge(state, u, neighbor,
                adjacent ? EDGE_PRESENT : EDGE_MISSING);
        }

        for (ui neighbor : q_neighbors[u]) {
            if (!isWhite(state, neighbor)) {
                continue;
            }
            if (!refreshWhiteByBlack(state, neighbor, u, v, next_cost)) {
                return false;
            }
        }
        return true;
    }

    bool tryMaterializeWhite(SearchState &state, ui cost, ui white_u,
        ui candidate, ui bucket_delta, ui &next_cost)
    {
        // 尝试把已选 white 点具体化为 black 映射。
        if (!isWhite(state, white_u) || candidate >= gn ||
            isDataVertexUsed(state, candidate) ||
            !candidates[white_u].contains(candidate)) {
            return false;
        }

        next_cost = cost + bucket_delta;
        if (next_cost > threshold) {
            return false;
        }

        setColor(state, white_u, COLOR_BLACK);
        assert(state.white_count > 0);
        setWhiteCnt(state, state.white_count - 1);
        setMap(state, white_u, (int)candidate);
        pushUsed(state, candidate);
        pushMatch(state, white_u, candidate);

        for (ui neighbor : q_neighbors[white_u]) {
            if (!isSelected(state, neighbor)) {
                continue;
            }
            if (isWhite(state, neighbor)) {
                return false;
            }

            ui mapped_neighbor = (ui)state.mapped_q[neighbor];
            bool adjacent = anchorAdjacent(
                neighbor, mapped_neighbor, white_u, candidate);
            EdgeState state_uv = getEdge(state, white_u, neighbor);
            if (state_uv == EDGE_PRESENT) {
                if (!adjacent) return false;
            }
            else if (state_uv == EDGE_MISSING) {
                if (adjacent) return false;
            }
            else {
                setEdge(state, white_u, neighbor,
                    adjacent ? EDGE_PRESENT : EDGE_MISSING);
            }
        }
        return true;
    }

    bool branchWhite(SearchState &state, ui cost, ui u)
    {
        // 分支：将未选查询点 u 设为 white，并递归继续搜索。
        if (state.color[u] != COLOR_UNSELECTED) {
            return false;
        }
        if (u >= static_color.size() || static_color[u] != COLOR_WHITE) {
            return false;
        }

        for (ui neighbor : q_neighbors[u]) {
            if (isWhite(state, neighbor)) {
                return false;
            }
        }

        if (!initWhiteCands(state, u, cost)) {
            return false;
        }
        if (candidate_result_buffer.empty()) {
            return false;
        }

        size_t undo_mark = mark();
        setColor(state, u, COLOR_WHITE);
        replaceBucket(state, u, candidate_result_buffer);
        setSelectedCnt(state, state.selected_count + 1);
        setWhiteCnt(state, state.white_count + 1);
        search(state, cost);
        rollback(state, undo_mark);
        return true;
    }

    bool branchBlack(SearchState &state, ui cost, ui u,
        ui required_anchor = std::numeric_limits<ui>::max())
    {
        // 分支：枚举未选查询点 u 的 black 映射候选。
        bool emitted_branch = false;
        auto try_candidate = [&](ui candidate) -> bool {
            size_t undo_mark = mark();
            ui next_cost = cost;
            if (!tryBindBlack(state, cost, u, candidate, next_cost)) {
                rollback(state, undo_mark);
                return false;
            }
            emitted_branch = true;
            search(state, next_cost);
            rollback(state, undo_mark);
            if (outputLimitReached()) {
                return true;
            }
            return false;
        };

#if CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_ANCHOR_BRANCH
        if (required_anchor < qn && isBlack(state, required_anchor) &&
            getEdge(state, u, required_anchor) == EDGE_PRESENT) {
            ui mapped_anchor = (ui)state.mapped_q[required_anchor];
            const pair<size_t, ui> *range = findAdjRange(
                required_anchor, mapped_anchor, u);
            if (range == nullptr) {
                return false;
            }

            for (const ui *it = rangeBegin(*range);
                it != rangeEnd(*range); ++it) {
                ui candidate = *it;
                if (try_candidate(candidate)) {
                    return true;
                }
            }
            return emitted_branch;
        }
#endif

        for (int candidate : candidates[u]) {
            if (try_candidate((ui)candidate)) {
                return true;
            }
        }
        return emitted_branch;
    }

    template <typename Continue>
    bool branchMatWhite(SearchState &state, ui cost,
        ui white_u, Continue continue_branch)
    {
        // 分支：枚举一个 white 点的具体映射，再交给后续回调继续。
        if (!isWhite(state, white_u)) {
            return false;
        }

        bool emitted_branch = false;
        WhiteCandidateBuckets white_bucket = state.white[white_u];
        assert(white_bucket.begin + white_bucket.count <=
            state.white_candidate_pool.size());
        for (ui candidate_idx = 0; candidate_idx < white_bucket.count;
            ++candidate_idx) {
            ui candidate =
                state.white_candidate_pool[white_bucket.begin + candidate_idx];
            if (isDataVertexUsed(state, candidate)) {
                continue;
            }

            ui delta = 0;
            if (!calcBlackDelta(state, white_u, candidate, cost, delta)) {
                continue;
            }

            size_t undo_mark = mark();
            ui next_cost = cost;
            if (!tryMaterializeWhite(state, cost, white_u,
                candidate, delta, next_cost)) {
                rollback(state, undo_mark);
                continue;
            }
            if (continue_branch(state, next_cost)) {
                emitted_branch = true;
            }
            rollback(state, undo_mark);
            if (outputLimitReached()) {
                return true;
            }
        }
        return emitted_branch;
    }

    template <typename Continue>
    bool branchMatWhites(SearchState &state, ui cost,
        const vector<ui> &white_vertices, size_t pos, Continue continue_branch)
    {
        // 分支：按顺序枚举一组 white 点的具体映射。
        if (pos == white_vertices.size()) {
            return continue_branch(state, cost);
        }

        ui white_u = white_vertices[pos];
        if (!isWhite(state, white_u)) {
            return branchMatWhites(state, cost, white_vertices,
                pos + 1, continue_branch);
        }

        return branchMatWhite(state, cost, white_u,
            [&](SearchState &next_state, ui next_cost) -> bool {
                return branchMatWhites(next_state, next_cost,
                    white_vertices, pos + 1, continue_branch);
            });
    }

    bool branchBlackAnchor(SearchState &state, ui cost, ui u,
        ui anchor)
    {
        // 在 anchor 已是 black 且边存在时，决定 u 走 white 或 black 分支。
        if (!isBlack(state, anchor) || state.color[u] != COLOR_UNSELECTED ||
            getEdge(state, u, anchor) != EDGE_PRESENT) {
            return false;
        }

        vector<ui> &white_neighbors =
            whiteNbrsBuffer(state.selected_count);
        collectWhiteNbrs(state, u, white_neighbors);
        bool emitted_white_branch = false;
        if (white_neighbors.empty()) {
            emitted_white_branch = branchWhite(state, cost, u);
        }
        else {
            emitted_white_branch = branchMatWhites(state, cost,
                white_neighbors, 0,
                [&](SearchState &materialized_state, ui materialized_cost) -> bool {
                    return branchWhite(materialized_state,
                        materialized_cost, u);
                });
        }

        if (emitted_white_branch) {
            return true;
        }
        return branchBlack(state, cost, u, anchor);
    }

    bool branchPresentEdge(SearchState &state, ui cost,
        const ActiveEdge &edge)
    {
        // 存在边分支：先标记活跃边存在，再扩展相关顶点。
        ui u = edge.u;
        ui anchor = edge.anchor;
        if (state.color[u] != COLOR_UNSELECTED ||
            !isSelected(state, anchor) ||
            getEdge(state, u, anchor) != EDGE_UNDECIDED) {
            return false;
        }

        size_t undo_mark = mark();
        setEdge(state, u, anchor, EDGE_PRESENT);

        bool emitted_branch = false;
        if (isBlack(state, anchor)) {
            emitted_branch = branchBlackAnchor(state, cost, u, anchor);
        }
        else {
            emitted_branch = branchMatWhite(state, cost, anchor,
            [&](SearchState &materialized_state, ui materialized_cost) -> bool {
                return branchBlackAnchor(materialized_state,
                    materialized_cost, u, anchor);
            });
        }
        rollback(state, undo_mark);
        return emitted_branch;
    }

    void branchEdges(SearchState &state, ui cost,
        const vector<ActiveEdge> &top_edges, size_t edge_idx)
    {
        // 对 top_edges 依次执行存在边/缺失边分支。
        if (outputLimitReached() || cost > threshold ||
            edge_idx >= top_edges.size()) {
            return;
        }

        const ActiveEdge &edge = top_edges[edge_idx];
        branchPresentEdge(state, cost, edge);
        if (outputLimitReached()) {
            return;
        }

        if (cost + 1 > threshold) {
            stats.prun_calls++;
            return;
        }

        size_t undo_mark = mark();
        if (state.color[edge.u] == COLOR_UNSELECTED &&
            isSelected(state, edge.anchor) &&
            getEdge(state, edge.u, edge.anchor) == EDGE_UNDECIDED) {
            setEdge(state, edge.u, edge.anchor, EDGE_MISSING);
            branchEdges(state, cost + 1, top_edges, edge_idx + 1);
        }
        rollback(state, undo_mark);
    }

    void search(SearchState &state, ui cost)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }
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

            vector<TerminalTailVertex> tail_vertices;
            if (!buildTailWhites(state, cost,
                tail_vertices)) {
                stats.prun_calls++;
                return;
            }
            enumTailWhites(state, 0, cost, tail_vertices);
            return;
        }

        vector<ActiveEdge> &top_edges =
            topEdgesBuffer(state.selected_count);
        ui max_branch_edges = threshold - cost + 1;
        if (!collectActiveEdges(state, max_branch_edges, top_edges)) {
            stats.prun_calls++;
            return;
        }

        branchEdges(state, cost, top_edges, 0);
    }

#endif
