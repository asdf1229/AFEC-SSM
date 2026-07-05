#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_STATE_UNDO_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_STATE_UNDO_H_

#ifndef CDE_BLACK_WHITE_INSIDE_MATCHING_SOLVER
#error "This internal header must be included from cde_black_white/context.h inside MatchingSolver."
#endif

    void initState(SearchState &state) const
    {
        state.mapped_q.assign(qn, -1);
        state.used_data_vertices.clear();
        state.used_data_vertices.reserve(qn);
        state.used_data_flag.assign(gn, 0);
        state.color.assign(qn, COLOR_UNSELECTED);
        state.edge_state.assign((size_t)qn * qn, EDGE_UNDECIDED);
        state.white.clear();
        state.white.resize(qn);
        state.white_candidate_pool.clear();
        state.white_candidate_pool.reserve(stats.filter_candidate_count);
        state.part_M.clear();
        state.part_M.reserve(qn);
        state.selected_count = 0;
        state.white_count = 0;
    }

    size_t edgeIdx(ui u, ui v) const
    {
        return (size_t)u * qn + v;
    }

    EdgeState getEdge(const SearchState &state, ui u, ui v) const
    {
        return state.edge_state[edgeIdx(u, v)];
    }

    void setEdgeRaw(SearchState &state, ui u, ui v,
        EdgeState edge_state_value) const
    {
        // 不记录 undo，直接写入有向查询边 (u, v) 的状态。
        state.edge_state[edgeIdx(u, v)] = edge_state_value;
    }

    vector<ActiveEdge> &topEdgesBuffer(ui depth)
    {
        if (top_edges_buffer_by_depth.size() <= depth) {
            top_edges_buffer_by_depth.resize((size_t)depth + 1);
        }
        vector<ActiveEdge> &buffer = top_edges_buffer_by_depth[depth];
        buffer.clear();
        if (buffer.capacity() < qn) {
            buffer.reserve(qn);
        }
        return buffer;
    }

    vector<ui> &whiteNbrsBuffer(ui depth)
    {
        // 获取指定搜索深度复用的 white 邻居缓冲区。
        if (white_neighbors_buffer_by_depth.size() <= depth) {
            white_neighbors_buffer_by_depth.resize((size_t)depth + 1);
        }
        vector<ui> &buffer = white_neighbors_buffer_by_depth[depth];
        buffer.clear();
        if (buffer.capacity() < qn) {
            buffer.reserve(qn);
        }
        return buffer;
    }

    bool tryBindRoot(SearchState &state, ui root, ui v) const
    {
        // 尝试把根查询点绑定到数据点 v，作为搜索初始 black 映射。
        if (root >= qn || v >= gn || !candidates[root].contains(v)) {
            return false;
        }
        state.color[root] = COLOR_BLACK;
        state.mapped_q[root] = (int)v;
        state.used_data_vertices.push_back(v);
        state.used_data_flag[v] = 1;
        state.part_M.push_back({ root, v });
        state.selected_count = 1;
        return true;
    }

    bool isDataVertexUsed(const SearchState &state, ui v) const
    {
        return v < state.used_data_flag.size() && state.used_data_flag[v] != 0;
    }

    bool isSelected(const SearchState &state, ui u) const
    {
        return u < state.color.size() && state.color[u] != COLOR_UNSELECTED;
    }

    bool isBlack(const SearchState &state, ui u) const
    {
        return u < state.color.size() && state.color[u] == COLOR_BLACK;
    }

    bool isWhite(const SearchState &state, ui u) const
    {
        return u < state.color.size() && state.color[u] == COLOR_WHITE;
    }

    size_t mark() const
    {
        // 返回当前 undo 栈大小，作为后续回滚标记。
        return undo_stack.size();
    }

    void rollback(SearchState &state, size_t mark)
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
                break;
            case UNDO_EDGE_STATE:
                setEdgeRaw(state, undo.u, undo.v, undo.old_edge_uv);
                setEdgeRaw(state, undo.v, undo.u, undo.old_edge_vu);
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
                state.white[undo.u] = undo.old_white;
                break;
            }
        }
    }

    void setMap(SearchState &state, ui u, int value)
    {
        // 设置查询点 u 的映射值，并记录 undo。
        UndoRecord undo;
        undo.kind = UNDO_MAPPED_Q;
        undo.u = u;
        undo.old_mapped_q = state.mapped_q[u];
        undo_stack.push_back(std::move(undo));
        state.mapped_q[u] = value;
    }

    void setColor(SearchState &state, ui u, VertexColor value)
    {
        // 设置查询点 u 的颜色，并记录 undo。
        UndoRecord undo;
        undo.kind = UNDO_COLOR;
        undo.u = u;
        undo.old_color = state.color[u];
        undo_stack.push_back(std::move(undo));
        state.color[u] = value;
    }

    void setEdge(SearchState &state, ui u, ui v,
        EdgeState edge_state_value)
    {
        // 对称设置查询边 (u, v) 的状态，并记录 undo。
        UndoRecord undo;
        undo.kind = UNDO_EDGE_STATE;
        undo.u = u;
        undo.v = v;
        undo.old_edge_uv = getEdge(state, u, v);
        undo.old_edge_vu = getEdge(state, v, u);
        undo_stack.push_back(std::move(undo));
        setEdgeRaw(state, u, v, edge_state_value);
        setEdgeRaw(state, v, u, edge_state_value);
    }

    void pushUsed(SearchState &state, ui v)
    {
        // 将数据点 v 标记为已使用，并记录 used 列表大小以便回滚。
        UndoRecord undo;
        undo.kind = UNDO_USED_DATA_SIZE;
        undo.old_size = state.used_data_vertices.size();
        undo_stack.push_back(std::move(undo));
        state.used_data_vertices.push_back(v);
        state.used_data_flag[v] = 1;
    }

    void pushMatch(SearchState &state, ui u, ui v)
    {
        // 将匹配对 (u, v) 加入部分匹配，并记录大小以便回滚。
        UndoRecord undo;
        undo.kind = UNDO_PART_M_SIZE;
        undo.old_size = state.part_M.size();
        undo_stack.push_back(std::move(undo));
        state.part_M.push_back({ u, v });
    }

    void setSelectedCnt(SearchState &state, ui value)
    {
        // 更新已选查询点数量，并记录 undo。
        UndoRecord undo;
        undo.kind = UNDO_SELECTED_COUNT;
        undo.old_count = state.selected_count;
        undo_stack.push_back(std::move(undo));
        state.selected_count = value;
    }

    void setWhiteCnt(SearchState &state, ui value)
    {
        // 更新 white 查询点数量，并记录 undo。
        UndoRecord undo;
        undo.kind = UNDO_WHITE_COUNT;
        undo.old_count = state.white_count;
        undo_stack.push_back(std::move(undo));
        state.white_count = value;
    }

    void replaceBucket(SearchState &state, ui u,
        const vector<ui> &candidates_to_store)
    {
        // 用新的候选列表替换 u 的 white bucket，并记录旧 bucket。
        UndoRecord undo;
        undo.kind = UNDO_WHITE_BUCKET;
        undo.u = u;
        undo.old_size = state.white_candidate_pool.size();
        undo.old_white = state.white[u];
        undo_stack.push_back(std::move(undo));

        WhiteCandidateBuckets bucket;
        bucket.begin = state.white_candidate_pool.size();
        assert(candidates_to_store.size() <=
            (size_t)std::numeric_limits<ui>::max());
        bucket.count = (ui)candidates_to_store.size();
        bucket.feasible_count = bucket.count;
        state.white_candidate_pool.insert(state.white_candidate_pool.end(),
            candidates_to_store.begin(), candidates_to_store.end());
        state.white[u] = bucket;
    }


#endif
