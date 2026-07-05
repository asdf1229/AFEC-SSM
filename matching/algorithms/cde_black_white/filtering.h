#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_FILTERING_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_FILTERING_H_

#ifndef CDE_BLACK_WHITE_INSIDE_WORKSPACE
#error "This internal header must be included from cde_black_white/context.h inside Workspace."
#endif

    // ========================================================================
    // Filtering
    // ========================================================================
    struct CandidateFilter {
    private:
        Workspace &solver;

        struct LabelFreq {
            LabelID label = 0;
            ui count = 0;
        };

        struct QueryLabelReq {
            LabelID label = 0;
            ui bridge = 0;
            ui normal = 0;
        };

        struct BridgeArc {
            ui from = 0;
            ui to = 0;
        };

        struct BridgeNbr {
            ui to = 0;
            ui support_arc = 0;
        };

        // static graph/profile cache
        vector<LabelID> query_label;
        vector<LabelID> data_label;
        vector<ui> data_degree;

        // NLF filtering
        vector<vector<LabelFreq>> data_label_freqs;
        vector<vector<QueryLabelReq>> query_label_reqs;
        vector<vector<ui>> data_by_label;
        vector<ui> query_bridge_degree;

        // Bridge filtering
        vector<BridgeArc> bridge_arcs;
        vector<vector<BridgeNbr>> bridge_nbrs;
        vector<ui> bridge_support;
        queue<pair<ui, ui>> removed;

        // spoke filtering
        vector<vector<ui>>  spoke_adj;
        vector<int>         match_right;
        vector<ui>          seen_right;
        ui                  seen_token = 1;
        vector<char>        left_is_bridge;
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
        queue<ui>           pending_spokes;
        vector<char>        queued_spoke;
#endif

    public:
        explicit CandidateFilter(Workspace &solver)
            : solver(solver),
            spoke_adj(solver.qn, vector<ui>()),
            match_right(solver.max_g_deg, -1),
            seen_right(solver.max_g_deg, 0),
            left_is_bridge(solver.qn, 0)
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
            , queued_spoke(solver.qn, 0)
#endif
        {
            // 绑定外层求解器，并初始化过滤阶段需要的临时匹配结构。
        }

        bool run()
        {
            // 依次执行桥边索引、NLF、桥边闭包和可选 spoke 过滤。
            if (!timed(&TimeStats::filter_bridge_time, [&] {
                buildBridgeIndex();
                return true;
            })) return false;

            if (!timed(&TimeStats::filter_nlf_time, [&] {
                return filterByNLF();
            })) return false;

#if CDE_BLACK_WHITE_ENABLE_BRIDGE_FILTERING
            if (!timed(&TimeStats::filter_bridge_time, [&] {
                return initBridgeSupport() && propFilter();
            })) return false;
#endif

#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
            if (!timed(&TimeStats::filter_spoke_time, [&] {
                pushAllSpokes();
                return propFilter();
            })) return false;
#endif

            updateCandidateCount();
            return true;
        }

    private:
        template <typename Fn>
        bool timed(long long TimeStats::*field, Fn &&fn)
        {
            // 执行一个过滤步骤，并在调试构建中累计该步骤耗时。
#ifndef NDEBUG
            Timer t;
#endif
            bool ok = fn();
#ifndef NDEBUG
            solver.stats.*field += t.elapsed();
#endif
            return ok;
        }

        void updateCandidateCount()
        {
            // 统计过滤后所有查询点的候选总数。
            solver.stats.filter_candidate_count = 0;
            for (ui u = 0; u < solver.qn; ++u) {
                solver.stats.filter_candidate_count += (ui)solver.candidates[u].size();
            }
        }

        ui addBridgeArc(ui from, ui to)
        {
            // 添加一条有向桥边弧，返回其编号。
            bridge_arcs.push_back({ from, to });
            return (ui)bridge_arcs.size() - 1;
        }

        void markBridgeNeighbor(ui from, ui to)
        {
            // 在查询邻接缓存中标记 from 到 to 这条边为桥边。
            bool marked = false;
            const vector<ui> &neighbors = solver.q_neighbors[from];
            for (ui i = 0; i < solver.q_degree[from]; ++i) {
                if (neighbors[i] == to) {
                    solver.q_neighbor_is_bridge[from][i] = 1;
                    marked = true;
                    break;
                }
            }
            assert(marked);
            (void)marked;
        }

        void addBridge(ui a, ui b)
        {
            // 记录一条无向桥边，并建立两个方向的支持依赖。
            ui ab = addBridgeArc(a, b);
            ui ba = addBridgeArc(b, a);
            bridge_nbrs[a].push_back({ b, ba });
            bridge_nbrs[b].push_back({ a, ab });
            markBridgeNeighbor(a, b);
            markBridgeNeighbor(b, a);
        }

        void tarjan(ui u, ui parent, vector<int> &dfn, vector<int> &low, int &time)
        {
            // 使用 Tarjan DFS 发现查询图中的桥边。
            dfn[u] = low[u] = ++time;
            for (ui v : solver.q_neighbors[u]) {
                if (dfn[v] == 0) {
                    tarjan(v, u, dfn, low, time);
                    low[u] = std::min(low[u], low[v]);
                    if (low[v] > dfn[u]) {
                        addBridge(u, v);
                    }
                }
                else if (v != parent) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
        }

        void buildBridgeIndex()
        {
            // 构建桥边邻接索引，并初始化每个查询邻接位置的桥边标记。
            bridge_arcs.clear();
            bridge_nbrs.assign(solver.qn, vector<BridgeNbr>());
            solver.q_neighbor_is_bridge.assign(solver.qn, vector<char>());

            for (ui u = 0; u < solver.qn; ++u) {
                solver.q_neighbor_is_bridge[u].assign(solver.q_degree[u], 0);
            }

#if CDE_BLACK_WHITE_ENABLE_BRIDGE_FILTERING
            vector<int> dfn(solver.qn, 0);
            vector<int> low(solver.qn, 0);
            int tim = 0;

            for (ui u = 0; u < solver.qn; ++u) {
                if (dfn[u] == 0) {
                    tarjan(u, solver.qn, dfn, low, tim);
                }
            }
#endif
        }

        void buildLabelCache()
        {
            // 缓存查询点标签、数据点标签和数据点度数。
            query_label.assign(solver.qn, 0);
            for (ui u = 0; u < solver.qn; ++u) {
                query_label[u] = solver.query_graph->getVertexLabel(u);
            }

            data_label.assign(solver.gn, 0);
            data_degree.assign(solver.gn, 0);
            for (ui v = 0; v < solver.gn; ++v) {
                data_label[v] = solver.data_graph->getVertexLabel(v);
                data_degree[v] = solver.data_graph->getVertexDegree(v);
            }
        }

        void buildQReqs()
        {
            // 按标签汇总每个查询点的桥边和非桥边邻居需求。
            query_label_reqs.assign(solver.qn, vector<QueryLabelReq>());
            query_bridge_degree.assign(solver.qn, 0);
            vector<ui> bridge_counts(solver.label_count, 0);
            vector<ui> non_bridge_counts(solver.label_count, 0);
            vector<ui> touched_labels;

            for (ui u = 0; u < solver.qn; ++u) {
                query_label_reqs[u].reserve(std::min(solver.q_degree[u], solver.label_count));
                touched_labels.clear();

                for (ui i = 0; i < solver.q_degree[u]; ++i) {
                    ui u1 = solver.q_neighbors[u][i];
                    LabelID label = query_label[u1];
                    assert(label >= 0 && (ui)label < solver.label_count);
                    ui label_idx = (ui)label;

                    if (bridge_counts[label_idx] == 0 && non_bridge_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    if (solver.q_neighbor_is_bridge[u][i]) {
                        query_bridge_degree[u]++;
                        bridge_counts[label_idx]++;
                    }
                    else {
                        non_bridge_counts[label_idx]++;
                    }
                }

                sort(touched_labels.begin(), touched_labels.end());
                for (ui label : touched_labels) {
                    ui bridge_count = bridge_counts[label];
                    ui non_bridge_count = non_bridge_counts[label];
                    query_label_reqs[u].push_back({
                        (LabelID)label, bridge_count, non_bridge_count
                    });
                    bridge_counts[label] = 0;
                    non_bridge_counts[label] = 0;
                }
            }
        }

        void buildGFreqs()
        {
            // 按标签汇总每个数据点的邻居频次，并建立数据点标签分桶。
            data_label_freqs.assign(solver.gn, vector<LabelFreq>());
            data_by_label.assign(solver.label_count, vector<ui>());
            vector<ui> label_counts(solver.label_count, 0);
            vector<ui> touched_labels;

            for (ui v = 0; v < solver.gn; ++v) {
                ui deg = 0;
                const ui *neighbors = solver.data_graph->getVertexNeighbors(v, deg);
                data_label_freqs[v].reserve(std::min(deg, solver.label_count));
                touched_labels.clear();

                LabelID vertex_label = data_label[v];
                assert(vertex_label >= 0 && (ui)vertex_label < solver.label_count);
                if (vertex_label >= 0 && (ui)vertex_label < solver.label_count) {
                    data_by_label[(ui)vertex_label].push_back(v);
                }

                for (ui i = 0; i < deg; ++i) {
                    LabelID label = data_label[neighbors[i]];
                    assert(label >= 0 && (ui)label < solver.label_count);
                    ui label_idx = (ui)label;
                    if (label_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    label_counts[label_idx]++;
                }

                sort(touched_labels.begin(), touched_labels.end());
                for (ui label : touched_labels) {
                    data_label_freqs[v].push_back({ (LabelID)label, label_counts[label] });
                    label_counts[label] = 0;
                }
            }
        }

        ui nlfDiff(ui u, ui v) const
        {
            // 计算查询点 u 映射到数据点 v 时的邻域标签频次缺口。
            ui diff = 0;
            const vector<QueryLabelReq> &query_reqs = query_label_reqs[u];
            const vector<LabelFreq> &data_freqs = data_label_freqs[v];
            size_t data_idx = 0;

            for (const auto &req : query_reqs) {
                while (data_idx < data_freqs.size() && data_freqs[data_idx].label < req.label) {
                    data_idx++;
                }

                ui data_count = 0;
                if (data_idx < data_freqs.size() && data_freqs[data_idx].label == req.label) {
                    data_count = data_freqs[data_idx].count;
                }

                ui bridge_need = req.bridge;
                if (bridge_need > data_count) {
                    return solver.threshold + 1;
                }

                ui non_bridge_need = req.normal;
                ui non_bridge_available = data_count - bridge_need;
                if (non_bridge_need > non_bridge_available) {
                    diff += (non_bridge_need - non_bridge_available);
                }
                if (diff > solver.threshold) {
                    return diff;
                }
            }
            return diff;
        }

        bool filterByNLF()
        {
            // 使用标签和邻域标签频次过滤初始候选集。
            buildLabelCache();
            buildQReqs();
            buildGFreqs();

            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = query_label[u];
                assert(lu >= 0 && (ui)lu < data_by_label.size());
                if (lu < 0 || (ui)lu >= data_by_label.size()) return false;

                for (ui v : data_by_label[(ui)lu]) {
                    if (query_bridge_degree[u] > data_degree[v]) continue;
                    if (solver.q_degree[u] > data_degree[v] + solver.threshold) continue;
                    if (nlfDiff(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        ui &support(ui arc_id, ui v)
        {
            // 返回桥边弧 arc_id 在数据点 v 上的支持计数引用。
            return bridge_support[(size_t)arc_id * solver.gn + v];
        }

        bool pruneCandidate(ui u, ui v)
        {
            // 删除候选 (u, v)，并把相关传播任务加入队列。
            if (!solver.candidates[u].contains(v)) {
                return true;
            }

            solver.candidates[u].remove(v);
            if (solver.candidates[u].empty()) {
                return false;
            }
            removed.push({ u, v });
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
            for (ui nbr_u : solver.q_neighbors[u]) {
                pushSpoke(nbr_u);
            }
#endif
            return true;
        }

        bool initBridgeSupport()
        {
            // 初始化桥边候选支持，并删除零支持候选。
            bridge_support.assign((size_t)bridge_arcs.size() * solver.gn, 0);
            vector<pair<ui, ui>> zero_support_candidates;

            for (ui arc_id = 0; arc_id < (ui)bridge_arcs.size(); ++arc_id) {
                const BridgeArc &arc = bridge_arcs[arc_id];
                for (ui v : solver.candidates[arc.from]) {
                    ui deg = 0;
                    const ui *nbrs = solver.data_graph->getVertexNeighbors(v, deg);
                    ui support = 0;
                    for (ui i = 0; i < deg; ++i) {
                        if (solver.candidates[arc.to].contains(nbrs[i])) {
                            support++;
                        }
                    }
                    this->support(arc_id, v) = support;
                    if (support == 0) {
                        zero_support_candidates.push_back({ arc.from, v });
                    }
                }
            }

            for (const auto &candidate : zero_support_candidates) {
                if (!pruneCandidate(candidate.first, candidate.second)) {
                    return false;
                }
            }
            return true;
        }

        bool propBridge()
        {
            // 沿桥边支持关系传播候选删除，直到删除队列清空。
            while (!removed.empty()) {
                ui removed_u = removed.front().first;
                ui removed_v = removed.front().second;
                removed.pop();

                ui deg = 0;
                const ui *nbrs = solver.data_graph->getVertexNeighbors(removed_v, deg);

                for (const BridgeNbr &bridge_nbr : bridge_nbrs[removed_u]) {
                    ui affected_u = bridge_nbr.to;
                    ui arc_id = bridge_nbr.support_arc;
                    for (ui i = 0; i < deg; ++i) {
                        ui v = nbrs[i];
                        if (!solver.candidates[affected_u].contains(v)) {
                            continue;
                        }
                        ui &support_count = support(arc_id, v);
                        if (support_count == 0) {
                            continue;
                        }
                        support_count--;
                        if (support_count == 0 &&
                            !pruneCandidate(affected_u, v)) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        bool propFilter()
        {
            // 交替执行桥边传播和 spoke 传播，直到达到过滤闭包。
            while (true) {
                if (!propBridge()) {
                    return false;
                }
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
                if (pending_spokes.empty()) {
                    break;
                }

                ui u = pending_spokes.front();
                pending_spokes.pop();
                queued_spoke[u] = 0;
                if (!filterSpoke(u)) {
                    return false;
                }
#else
                break;
#endif
            }
            return true;
        }

        bool augmentSpoke(ui left_idx)
        {
            // 在 spoke 二分图中为指定左点寻找增广路径。
            for (ui right_idx : spoke_adj[left_idx]) {
                if (seen_right[right_idx] == seen_token) continue;
                seen_right[right_idx] = seen_token;
                if (match_right[right_idx] < 0 ||
                    augmentSpoke((ui)match_right[right_idx])) {
                    match_right[right_idx] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        bool tryAugmentSpoke(ui left_idx)
        {
            // 刷新访问 token 后尝试为一个 spoke 左点增广。
            seen_token++;
            if (seen_token == 0) {
                std::fill(seen_right.begin(), seen_right.end(), 0);
                seen_token = 1;
            }
            return augmentSpoke(left_idx);
        }

        void buildSpokeAdj(ui u, ui v, ui &deg_u, ui &deg_v)
        {
            // 为候选 (u, v) 建立查询邻居到数据邻居的 spoke 二分图。
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            deg_u = solver.q_degree[u];
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            for (ui i = 0; i < deg_u; ++i) {
                spoke_adj[i].clear();
                ui u1 = u_neighbors[i];
                left_is_bridge[i] = solver.q_neighbor_is_bridge[u][i];
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        spoke_adj[i].push_back(j);
                    }
                }
            }
        }

        bool checkSpoke(ui u, ui v, ui budget)
        {
            // 检查候选 (u, v) 的 spoke 匹配是否满足缺边预算。
            ui deg_u = 0;
            ui deg_v = 0;
            buildSpokeAdj(u, v, deg_u, deg_v);

            std::fill(match_right.begin(), match_right.begin() + deg_v, -1);

            ui optional_count = 0;
            for (ui i = 0; i < deg_u; ++i) {
                if (left_is_bridge[i]) {
                    if (!tryAugmentSpoke(i)) {
                        return false;
                    }
                }
                else {
                    optional_count++;
                }
            }

            if (optional_count <= budget) {
                return true;
            }

            ui required = optional_count - budget;
            ui matched = 0;
            ui processed = 0;

            for (ui i = 0; i < deg_u; ++i) {
                if (left_is_bridge[i]) {
                    continue;
                }

                processed++;
                if (tryAugmentSpoke(i)) {
                    matched++;
                    if (matched >= required) {
                        return true;
                    }
                }

                ui remaining = optional_count - processed;
                if (matched + remaining < required) {
                    return false;
                }
            }

            return matched >= required;
        }

        ui edgeBudget(ui u) const
        {
            // 返回查询点 u 的关联边最多可缺失数量。
            if (solver.q_degree[u] == 0) return 0;
            return std::min(solver.threshold, solver.q_degree[u] - 1);
        }

#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
        void pushSpoke(ui u)
        {
            // 将查询点加入 spoke 待处理队列，避免重复入队。
            if (queued_spoke[u]) {
                return;
            }
            pending_spokes.push(u);
            queued_spoke[u] = 1;
        }

        void pushAllSpokes()
        {
            // 将所有查询点加入 spoke 队列，作为初始传播入口。
            for (ui u = 0; u < solver.qn; ++u) {
                pushSpoke(u);
            }
        }

        bool filterSpoke(ui u)
        {
            // 删除查询点 u 下所有不满足 spoke 约束的候选。
            ui budget = edgeBudget(u);

            vector<ui> to_remove;
            for (ui v : solver.candidates[u]) {
                if (!checkSpoke(u, v, budget)) {
                    to_remove.push_back(v);
                }
            }
            if (to_remove.empty()) return true;

            for (ui v : to_remove) {
                if (!pruneCandidate(u, v)) {
                    return false;
                }
            }
            return true;
        }
#endif
    };

    bool runCandidateFiltering()
    {
        return CandidateFilter(*this).run();
    }


#endif
