#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_TERMINAL_TAIL_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_TERMINAL_TAIL_H_

#ifndef CDE_BLACK_WHITE_INSIDE_MATCHING_SOLVER
#error "This internal header must be included from cde_black_white/context.h inside MatchingSolver."
#endif

    bool buildTailBuckets(
        const SearchState &state, ui white_u, ui cost,
        vector<vector<ui>> &buckets, ui &feasible_count, ui &min_delta)
    {
        // 为终端阶段的一个 white 点按缺边增量构建候选分桶。
        assert(cost <= threshold);
        if (!isWhite(state, white_u) || state.mapped_q[white_u] != -1) {
            return false;
        }

        for (ui neighbor : q_neighbors[white_u]) {
            if (isWhite(state, neighbor)) {
                return false;
            }
        }

        ui remaining_budget = threshold - cost;
        buckets.assign((size_t)remaining_budget + 1, vector<ui>());
        feasible_count = 0;
        min_delta = std::numeric_limits<ui>::max();

        WhiteCandidateBuckets white = state.white[white_u];
        assert(white.begin + white.count <= state.white_candidate_pool.size());
        for (ui candidate_idx = 0; candidate_idx < white.count;
            ++candidate_idx) {
            ui candidate =
                state.white_candidate_pool[white.begin + candidate_idx];
            if (isDataVertexUsed(state, candidate)) {
                continue;
            }

            ui delta = 0;
            if (!calcBlackDelta(state, white_u, candidate, cost, delta)) {
                continue;
            }
            if (delta > remaining_budget) {
                continue;
            }

            buckets[delta].push_back(candidate);
            feasible_count++;
            if (delta < min_delta) {
                min_delta = delta;
            }
        }

        return feasible_count > 0;
    }

    bool buildTailWhites(const SearchState &state,
        ui cost, vector<TerminalTailVertex> &tail_vertices)
    {
        // 为所有剩余 white 点构建终端 tail 枚举结构。
        tail_vertices.clear();
        tail_vertices.reserve(state.white_count);

        for (ui u = 0; u < qn; ++u) {
            if (!isWhite(state, u)) {
                continue;
            }

            TerminalTailVertex tail_vertex;
            tail_vertex.u = u;
            if (!buildTailBuckets(state, u, cost,
                tail_vertex.buckets, tail_vertex.feasible_count,
                tail_vertex.min_delta)) {
                return false;
            }
            tail_vertices.push_back(std::move(tail_vertex));
        }

        if (tail_vertices.size() != (size_t)state.white_count) {
            return false;
        }

        std::sort(tail_vertices.begin(), tail_vertices.end(),
            [this](const TerminalTailVertex &lhs, const TerminalTailVertex &rhs) {
                if (lhs.feasible_count != rhs.feasible_count) {
                    return lhs.feasible_count < rhs.feasible_count;
                }
                if (lhs.min_delta != rhs.min_delta) {
                    return lhs.min_delta < rhs.min_delta;
                }
                if (q_degree[lhs.u] != q_degree[rhs.u]) {
                    return q_degree[lhs.u] > q_degree[rhs.u];
                }
                return lhs.u < rhs.u;
            });
        return true;
    }

    void enumTailWhites(SearchState &state,
        size_t pos, ui cost, vector<TerminalTailVertex> &tail_vertices)
    {
        // 递归枚举终端阶段所有 white 点的具体映射。
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }
        stats.recursion_calls++;
        stats.terminal_tail_calls++;
        assert(cost <= threshold);

        if (pos == tail_vertices.size()) {
            emitResult(state);
            return;
        }

        TerminalTailVertex &tail_vertex = tail_vertices[pos];
        ui u = tail_vertex.u;
        assert(isWhite(state, u));
        assert(state.mapped_q[u] == -1);

        ui remaining_budget = threshold - cost;
        ui max_delta = std::min((ui)tail_vertex.buckets.size() - 1,
            remaining_budget);
        for (ui missing_delta = 0; missing_delta <= max_delta; ++missing_delta) {
            const vector<ui> &bucket = tail_vertex.buckets[missing_delta];
            for (ui v : bucket) {
                if (isDataVertexUsed(state, v)) {
                    continue;
                }

                state.mapped_q[u] = (int)v;
                state.used_data_vertices.push_back(v);
                state.used_data_flag[v] = 1;
                state.part_M.push_back({ u, v });

                enumTailWhites(state, pos + 1,
                    cost + missing_delta, tail_vertices);

                state.part_M.pop_back();
                state.used_data_flag[v] = 0;
                state.used_data_vertices.pop_back();
                state.mapped_q[u] = -1;

                if (outputLimitReached()) {
                    return;
                }
            }
        }
    }

    void emitResult(const SearchState &state)
    {
        // 输出一个完整匹配，并更新结果计数和输出上限状态。
        assert(state.part_M.size() == qn);
        stats.result_count++;
        noteOutputLimitIfReached();
#ifndef NDEBUG
        results_ptr->push_back(state.part_M);
#endif
    }


#endif
