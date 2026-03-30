#include "utility/utility.h"
// #include "graph/graph.h"

using namespace std;

const int INT_MAX = 0x7fffffff;

// ============================================================================
// MatchingSolver Implementation (CDE + Dynamic Kernel-and-Shell Strategy)
// ============================================================================

class MatchingSolver {
public:
    MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

    bool init(const Graph *q, const Graph *g, ui match_threshold)
    {
        Timer t_init;
        t_init.restart();

        query_graph = q;
        data_graph = g;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();

        resetState();

        // 初始化每个查询顶点的未映射邻居数 (初始即为其度数)
        for (ui i = 0; i < qn; ++i) {
            unmapped_degrees[i] = query_graph->getVertexDegree(i);
        }

        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

        Timer t_filter;
        bool res = calVerticesFilter();
        stats.filter_time = t_filter.elapsed();
        if (!res) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        Timer t_order;
        generateMatchingOrder();
        stats.order_time = t_order.elapsed();

        stats.init_time = t_init.elapsed();
        return true;
    }

    void match(vector<vector<pair<ui, ui>>> &results)
    {
        Timer t_search;
        t_search.restart();

        results_ptr = &results;
        results_ptr->clear();

        ui v = _order[0];
        ui deg; const ui *nbrs = query_graph->getVertexNeighbors(v, deg);

        for (ui v0 : candidates[v]) {
            mapped_q[v] = v0;
            mapped_g[v0] = v;
            in_Mq[v] = 1;
            matched_count = 1;
            part_M.push_back({ v, v0 }); // v -> v0

            updateFrontierStatus(v);

            for (ui i = 0; i < deg; ++i) {
                ui nbr = nbrs[i];
                unmapped_degrees[nbr]--;
                if (!is_excluded[nbr][v]) {
                    anchor_count[nbr]++;
                }
                updateFrontierStatus(nbr);
            }

            dfs(0, (int)v);

            for (ui i = 0; i < deg; ++i) {
                ui nbr = nbrs[i];
                unmapped_degrees[nbr]++;
                if (!is_excluded[nbr][v]) {
                    anchor_count[nbr]--;
                }
                updateFrontierStatus(nbr);
            }

            mapped_q[v] = -1;
            mapped_g[v0] = -1;
            in_Mq[v] = 0;
            matched_count = 0;
            part_M.pop_back();

            updateFrontierStatus(v);
        }

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct TimeStats {
        long long total_time = 0;
        // init breakdown
        long long init_time = 0;
        long long filter_time = 0;
        long long order_time = 0;
        // search breakdown
        long long dfs_time = 0;
        long long branch_time = 0;   // candidate enumeration & matching in dfs
        long long lb_time = 0;       // computeLowerBound
        long long frontier_time = 0; // building U_frontier
        long long shell_time = 0;    // NEW: shell batch processing time
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        printf("\n--- CDE-Match with Dynamic KSS Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("  - Order Time:      %.4lf ms (%.2f%% of Init)\n", stats.order_time / 1000.0, pct(stats.order_time, stats.init_time));
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("  - Shell Time:      %.4lf ms (%.2f%% of Search)\n", stats.shell_time / 1000.0, pct(stats.shell_time, stats.dfs_time));
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
        printf("Results Found:       %zu\n", results_ptr ? results_ptr->size() : 0);
        printf("-----------------------------------------------------------\n");
    }

private:
    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;
    ui threshold;
    ui qn, gn;

    vector<vector<ui>> candidates;
    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;

    // --- CDE & Dynamic KSS Core States ---
    vector<int> mapped_q;            // f: query vertex → data vertex (-1 if unmapped)
    vector<int> mapped_g;            // reverse: data vertex → query vertex (-1 if unused)
    vector<char> in_Mq;              // 当前已匹配集合 M_q
    vector<vector<int>> is_excluded; // X: 排除边集合 (symmetric, is_excluded[u][v] > 0 means (u,v) ∈ X)
    vector<char> in_P;               // P: 延期顶点集合
    ui matched_count;
    vector<pair<ui, ui>> part_M;

    // 动态维护每个顶点有多少个邻居还尚未被映射
    vector<ui> unmapped_degrees;

    vector<ui> _order;               // 固定搜索顺序 π
    vector<ui> order_rank;           // 用于还原原始匹配顺序的排名

    // --- Incremental U_frontier States ---
    vector<ui> anchor_count;         // 记录每个顶点当前有多少个未排斥(un-excluded)的已匹配邻居
    vector<int> frontier_pos;        // 记录顶点在 active_frontier 中的位置，-1表示不在其中
    vector<ui> active_frontier;      // 动态维护的 U_frontier 集合


    // O(1) 增量更新顶点 u 的 Frontier 状态
    inline void updateFrontierStatus(ui u)
    {
        bool should_be = (!in_Mq[u] && !in_P[u] && unmapped_degrees[u] > 0 && anchor_count[u] > 0);
        bool is_in = (frontier_pos[u] != -1);

        if (should_be && !is_in) {
            frontier_pos[u] = active_frontier.size();
            active_frontier.push_back(u);
        }
        else if (!should_be && is_in) {
            ui idx = frontier_pos[u];
            ui last_u = active_frontier.back();
            active_frontier[idx] = last_u;
            frontier_pos[last_u] = idx;
            active_frontier.pop_back();
            frontier_pos[u] = -1;
        }
    }

    void resetState()
    {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        in_Mq.assign(qn, 0);
        is_excluded.assign(qn, vector<int>(qn, 0));
        in_P.assign(qn, 0);
        matched_count = 0;
        part_M.clear();
        part_M.reserve(qn);
        candidates.clear();
        candidates.resize(qn);
        unmapped_degrees.assign(qn, 0);
        _order.clear();

        anchor_count.assign(qn, 0);
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
        order_rank.assign(qn, 0);

        stats = TimeStats();
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees)
    {
        ui n = g->getVerticesCount();
        ui num_labels = g->getLabelsCount();
        counts.assign(n, vector<ui>(num_labels, 0));
        degrees.assign(n, 0);
        for (ui u = 0; u < n; ++u) {
            ui deg; const ui *neighbors = g->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                counts[u][g->getVertexLabel(neighbors[i])]++;
                degrees[u]++;
            }
        }
    }

    ui computeNLF(ui u, ui v)
    {
        ui diff = 0;
        size_t sz = Lq_counts[u].size();
        for (size_t i = 0; i < sz; ++i) {
            if (Lq_counts[u][i] > Lg_counts[v][i]) {
                diff += (Lq_counts[u][i] - Lg_counts[v][i]);
            }
        }
        return diff;
    }

    bool calVerticesFilter()
    {
        for (ui u = 0; u < qn; ++u) {
            LabelID label_u = query_graph->getVertexLabel(u);
            for (ui v = 0; v < gn; ++v) {
                if (label_u != data_graph->getVertexLabel(v)) continue;
                if (Lq_degrees[u] > Lg_degrees[v] + threshold) continue;
                if (computeNLF(u, v) <= threshold) {
                    candidates[u].push_back(v);
                }
            }
            if (candidates[u].empty()) return false;
            sort(candidates[u].begin(), candidates[u].end());
        }

#ifndef NDEBUG
        printf("candidates nums:\n");
        for (ui u = 0; u < qn; ++u) {
            printf("u = %u: %zu candidates\n", u, candidates[u].size());
        }
#endif

        return true;
    }

    void generateMatchingOrder()
    {
        _order.clear();
        vector<bool> visited(qn, false);

        ui root = 0;
        size_t min_cand = candidates[0].size();
        for (ui u = 1; u < qn; ++u) {
            if (candidates[u].size() < min_cand) {
                min_cand = candidates[u].size();
                root = u;
            }
        }

        _order.push_back(root);
        visited[root] = true;

        for (ui step = 1; step < qn; ++step) {
            bool found = false;
            ui best_u = 0;
            ui best_anchor_sz = qn;
            ui best_cand_sz = qn;

            for (ui u = 0; u < qn; ++u) {
                if (visited[u]) continue;

                // anchor size
                ui current_anchor_sz = 0;
                ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
                for (ui i = 0; i < deg; ++i) if (visited[nbrs[i]]) current_anchor_sz++;

                if (current_anchor_sz == 0) continue;

                ui current_cand_sz = candidates[u].size();

                bool is_better = false;

                if (found == false) {
                    found = true;
                    is_better = true;
                }
                else {
                    if (current_cand_sz < best_cand_sz) {
                        is_better = true;
                    }
                    else if (current_cand_sz == best_cand_sz) {
                        if (current_anchor_sz > best_anchor_sz) { // modified
                            is_better = true;
                        }
                        else if (current_anchor_sz == best_anchor_sz) {
                            if (u < best_u) is_better = true;
                        }
                    }
                }

                if (is_better) {
                    best_u = u;
                    best_anchor_sz = current_anchor_sz;
                    best_cand_sz = current_cand_sz;
                }
            }

            assert(found);

            visited[best_u] = true;
            _order.push_back(best_u);
        }

        // 缓存排名以便于后续 O(F log F) 的 Frontier 快速定序
        order_rank.assign(qn, 0);
        for (ui i = 0; i < qn; ++i) {
            order_rank[_order[i]] = i;
        }

#ifndef NDEBUG
        printf("matching order: ");
        for (auto u : _order) {
            printf("%u ", u);
        }
        printf("\n");
#endif
    }

    // TODO
    ui computeLowerBound(ui current_miss)
    {
        return current_miss;
    }

    // 在查询图中，仅使用 (E_q \ X \ E_cut(u)) 的边，
    // 检查 u 是否还能到达 M_q 中任意点。
    bool isConnectionMandatory(ui u)
    {
        vector<char> visited(qn, 0);
        queue<ui> bfs_q;
        visited[u] = 1;
        bfs_q.push(u);

        while (!bfs_q.empty()) {
            ui cur = bfs_q.front();
            bfs_q.pop();

            ui deg; const ui *nbrs = query_graph->getVertexNeighbors(cur, deg);
            for (ui i = 0; i < deg; ++i) {
                ui nxt = nbrs[i];

                if (visited[nxt]) continue;

                // 去掉 X 中的边
                if (is_excluded[cur][nxt]) continue;

                // 去掉 E_cut(u): u 与当前 M_q 之间的边
                if ((cur == u && in_Mq[nxt]) || (nxt == u && in_Mq[cur])) continue;

                visited[nxt] = 1;

                if (in_Mq[nxt]) return false;

                bfs_q.push(nxt);
            }
        }

        return true;
    }

    // =====================================================
    // DFS 搜索过程
    // 包含 CDE 原有逻辑及 Dynamic KSS (Kernel-and-Shell) 机制
    // 对应伪代码: Procedure DFS(M_part, cost, X, P, u_new)
    //
    // cost:  当前部分映射中已确认缺失的查询边数量
    // u_new: 最近一次新加入映射的查询顶点 (-1 表示 undefined)
    // =====================================================
    void dfs(ui cost, int u_new)
    {
        assert(matched_count > 0);
        assert(matched_count <= qn);
        assert(cost <= threshold);

        Timer t_lb;
        ui lb = computeLowerBound(cost);
        stats.lb_time += t_lb.elapsed();

        if (lb > threshold) {
            stats.prun_calls++;
            return;
        }

        stats.recursion_calls++;

        if (matched_count == qn) {
            results_ptr->push_back(part_M);
            return;
        }

        vector<ui> reEnabledP;
        if (u_new >= 0) {
            ui deg; const ui *nbrs = query_graph->getVertexNeighbors((ui)u_new, deg);
            for (ui i = 0; i < deg; ++i) {
                ui w = nbrs[i];
                if (in_P[w]) {
                    reEnabledP.push_back(w);
                    in_P[w] = 0;
                    updateFrontierStatus(w);
                }
            }
        }

        // ====================================================================
        // NEW: 动态 Kernel & Shell 判定
        // 当所有未分配顶点的所有邻居都已被访问（即未映射度数均为 0）时，
        // 说明它们互不相连构成绝对的独立集 (Shell)。此时拦截下来批量枚举。
        // ====================================================================
        bool all_shells = true;
        vector<ui> local_shells;
        for (ui u : _order) {
            if (!in_Mq[u]) {
                if (unmapped_degrees[u] == 0) {
                    local_shells.push_back(u);
                }
                else {
                    all_shells = false;
                    break;
                }
            }
        }

        if (all_shells && !local_shells.empty()) {
            Timer t_shell;
            processShell(0, cost, local_shells);
            stats.shell_time += t_shell.elapsed();

            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }

        Timer t_frontier;
        vector<ui> U_frontier = active_frontier;
        // TODO
        if (U_frontier.size() > 1) {
            // 通过极小规模排序来保证原本的遍历顺序
            sort(U_frontier.begin(), U_frontier.end(), [&](ui a, ui b) {
                return order_rank[a] < order_rank[b];
                });
        }
        stats.frontier_time += t_frontier.elapsed();

        if (U_frontier.empty()) {
            stats.prun_calls++;
            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }

        Timer t_branch;
        long long child_dfs_time = 0;
        vector<ui> local_P;
        vector<pair<ui, ui>> local_X;
        ui current_cost = cost;

        for (ui u : U_frontier) {
            ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            vector<ui> U_anchor;
            for (ui i = 0; i < deg; ++i) if (in_Mq[nbrs[i]]) U_anchor.push_back(nbrs[i]);

            // --- matching u (Kernel 扩展) ---
            for (ui v : candidates[u]) {
                if (mapped_g[v] != -1) continue;

                bool conflict = false;
                bool connects = false;
                ui delta = 0;

                for (ui ua : U_anchor) {
                    bool has_edge = data_graph->hasEdge(v, mapped_q[ua]);
                    if (is_excluded[u][ua]) {
                        if (has_edge) { conflict = true; break; }
                    }
                    else {
                        if (has_edge) connects = true;
                        else delta++;
                    }
                }

                if (conflict || !connects || current_cost + delta > threshold) continue;

                mapped_q[u] = v;
                mapped_g[v] = u;
                in_Mq[u] = 1;
                matched_count++;
                part_M.push_back({ u, v });

                updateFrontierStatus(u);

                for (ui i = 0; i < deg; ++i) {
                    ui nbr = nbrs[i];
                    unmapped_degrees[nbr]--;
                    if (!is_excluded[nbr][u]) {
                        anchor_count[nbr]++;
                    }
                    updateFrontierStatus(nbr);
                }

                Timer t_child;
                dfs(current_cost + delta, (int)u);
                child_dfs_time += t_child.elapsed();


                for (ui i = 0; i < deg; ++i) {
                    ui nbr = nbrs[i];
                    unmapped_degrees[nbr]++;
                    if (!is_excluded[nbr][u]) {
                        anchor_count[nbr]--;
                    }
                    updateFrontierStatus(nbr);
                }

                part_M.pop_back();
                matched_count--;
                in_Mq[u] = 0;
                mapped_g[v] = -1;
                mapped_q[u] = -1;

                updateFrontierStatus(u);
            }

            if (isConnectionMandatory(u)) break;

            // --- delay u ---
            ui delta = 0;
            for (ui ua : U_anchor) if (!is_excluded[u][ua]) if (++current_cost + delta > threshold) break;

            in_P[u] = 1;
            local_P.push_back(u);
            updateFrontierStatus(u);

            for (ui ua : U_anchor) {
                if (!is_excluded[u][ua]) {
                    is_excluded[u][ua] = 1;
                    is_excluded[ua][u] = 1;
                    local_X.push_back({ u, ua });

                    anchor_count[u]--;
                }
            }
        }

        stats.branch_time += t_branch.elapsed() - child_dfs_time;

        // backtracking
        for (ui u : local_P) {
            in_P[u] = 0;
            updateFrontierStatus(u);
        }
        for (auto &e : local_X) {
            is_excluded[e.first][e.second] = 0;
            is_excluded[e.second][e.first] = 0;

            if (in_Mq[e.second]) anchor_count[e.first]++;
            if (in_Mq[e.first]) anchor_count[e.second]++;
            updateFrontierStatus(e.first);
            updateFrontierStatus(e.second);
        }

        for (ui w : reEnabledP) {
            in_P[w] = 1;
            updateFrontierStatus(w);
        }
    }

    // 独立集 (Shell) 批量枚举处理
    // 注意：因为 Shell 之间无边，所以在此函数中不需要再维护 unmapped_degrees
    void processShell(size_t idx, ui cost, const vector<ui> &shells)
    {
        stats.recursion_calls++;

        if (idx == shells.size()) {
            results_ptr->push_back(part_M);
            return;
        }

        ui u = shells[idx];
        ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);

        vector<ui> U_anchor;
        for (ui i = 0; i < deg; ++i) {
            // Shell 此时的定义保证了它的所有邻居必定均在 M_q 中
            U_anchor.push_back(nbrs[i]);
        }

        for (ui v : candidates[u]) {
            if (mapped_g[v] != -1) continue;

            bool conflict = false;
            bool connects = false;
            ui delta = 0;

            for (ui ua : U_anchor) {
                bool has_edge = data_graph->hasEdge(v, mapped_q[ua]);
                // Shell 节点之前可能由于 CDE 被置入排斥集合 X
                if (is_excluded[u][ua]) {
                    if (has_edge) { conflict = true; break; }
                }
                else {
                    if (has_edge) connects = true;
                    else delta++;
                }
            }

            if (conflict || !connects || cost + delta > threshold) continue;

            mapped_q[u] = v;
            mapped_g[v] = u;
            in_Mq[u] = 1;
            matched_count++;
            part_M.push_back({ u, v });

            processShell(idx + 1, cost + delta, shells);

            part_M.pop_back();
            matched_count--;
            in_Mq[u] = 0;
            mapped_g[v] = -1;
            mapped_q[u] = -1;
        }
    }
};


// ============================================================================
// TreeSpanSolver Implementation (EnumerateOnDemand Strategy)
// ============================================================================

class TreeSpanSolver {
public:
    TreeSpanSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

    bool init(const Graph *q, const Graph *d, ui match_threshold)
    {
        Timer t_init; // Start Init Timer
        t_init.restart();

        // Reset Stats
        stats = TreeStats();

        query_graph = q;
        data_graph = d;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();

        // 1. 初始化全局映射和候选集
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);

        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

        candidates.assign(qn, vector<ui>());
        if (!calVerticesFilter()) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        // 2. 预处理所有查询边（用于生成树时的边选择）
        all_q_edges = getAllQueryEdges();

        stats.init_time = t_init.elapsed();
        return true;
    }

    void match(vector<vector<pair<ui, ui>>> &results)
    {
        Timer t_search; // Start Search Timer
        t_search.restart();

        results_ptr = &results;
        results_ptr->clear();

        // 3. 生成初始 QISequence (基于 Prim 算法)
        // 选择候选集最小的点作为根
        ui root = selectRoot();
        QISequence initial_seq = generateInitialSequence(root);

        // 4. 开始 SimSearch (Algorithm 2)
        // h=1 因为 seq[0] 是根节点，不需要匹配边，只需要匹配点
        // 我们在 h=0 处做特殊处理，循环根节点的候选集，然后从 h=1 开始递归

        ui start_u = initial_seq.S[0];

        for (ui v : candidates[start_u]) {
            mapped_q[start_u] = v;
            mapped_g[v] = start_u;

            // 参数: depth h, current sequence, missing_edges gamma
            SimSearch(1, initial_seq, 0);

            mapped_q[start_u] = -1;
            mapped_g[v] = -1;
        }

        stats.search_time = t_search.elapsed();
    }

    // --- Added Stats Structure ---
    struct TreeStats {
        long long total_time = 0;
        long long init_time = 0;
        long long search_time = 0;      // Total time in match() excluding init
        long long verify_time = 0;      // Time spent checking edge existence/gamma
        long long reorder_time = 0;     // Time spent finding replacements and re-running Prim
        long long recursion_calls = 0;
        long long reorder_calls = 0;    // How many times the tree structure was changed
    } stats;

    void printStats() const
    {
        printf("\n--- TreeSpan Matching Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, (double)stats.init_time / stats.total_time * 100);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, (double)stats.search_time / stats.total_time * 100);

        // Percentages below are relative to Search Time
        printf("- Verify Time:       %.4lf ms (%.2f%% of Search)\n", stats.verify_time / 1000.0, (stats.search_time > 0 ? (double)stats.verify_time / stats.search_time * 100 : 0));
        printf("- Reorder Time:      %.4lf ms (%.2f%% of Search)\n", stats.reorder_time / 1000.0, (stats.search_time > 0 ? (double)stats.reorder_time / stats.search_time * 100 : 0));

        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Reorder Ops:         %lld\n", stats.reorder_calls);
        printf("Results Found:       %zu\n", results_ptr ? results_ptr->size() : 0);
        printf("---------------------------------------\n");
    }

private:
    // --- Data Structures ---

    struct QEdge {
        ui u, v;
        bool operator==(const QEdge &other) const
        {
            return (u == other.u && v == other.v) || (u == other.v && v == other.u);
        }
        // 为了放入 set 或比较，标准化为 u < v
        QEdge canonical() const
        {
            return (u < v) ? *this : QEdge{ v, u };
        }
        bool operator<(const QEdge &other) const
        {
            if (u != other.u) return u < other.u;
            return v < other.v;
        }
    };

    // QISequence: 对应论文中的 "seq"
    // S: 顶点访问顺序
    // sEdge: S[i] 对应的生成树边 (连接 S[i] 和 S[0...i-1] 中的某点)
    // bEdges: S[i] 的反向边集合
    // R: Edge Exclusion Set (排斥边集)
    struct QISequence {
        vector<ui> S;             // Vertices in order
        vector<QEdge> sEdge;      // Spanning edges. sEdge[i] is edge for S[i]
        vector<vector<QEdge>> bEdges; // Backward edges for S[i]
        set<QEdge> R;             // Exclusion set
    };

    const Graph *query_graph;
    const Graph *data_graph;
    ui threshold;
    ui qn, gn;
    vector<vector<pair<ui, ui>>> *results_ptr;

    vector<vector<ui>> candidates;
    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;
    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<QEdge> all_q_edges;

    // --- Core Algorithm: SimSearch (Algorithm 2) ---

    /**
     * @brief 递归搜索 + 按需生成生成树
     *
     * @param h 当前匹配深度 (对应 seq.S[h])
     * @param seq 当前使用的 QISequence
     * @param gamma 当前缺失边数量
     */
    void SimSearch(ui h, QISequence seq, ui gamma)
    {
        stats.recursion_calls++;

        // Base Case: 成功匹配所有点
        if (h == qn) {
            vector<pair<ui, ui>> res;
            res.reserve(qn);
            for (ui u = 0; u < qn; ++u) res.push_back({ u, (ui)mapped_q[u] });
            results_ptr->push_back(res);
            return;
        }

        // --- Phase 1: Go-Down (纵向扩展) ---
        // 尝试为 seq.S[h] 寻找匹配点
        ui u_curr = seq.S[h];
        QEdge tree_edge = seq.sEdge[h];

        // 确定父节点 (树边是 (parent, u_curr))
        ui u_parent = (tree_edge.u == u_curr) ? tree_edge.v : tree_edge.u;
        ui v_parent = mapped_q[u_parent];

        // 获取父节点在数据图中的邻居作为候选
        ui deg_g; const ui *neighbors_g = data_graph->getVertexNeighbors(v_parent, deg_g);

        // 优化: 先将 u_curr 的候选集放入 bloom filter 或 hash set 加速查找? 
        // 这里直接使用 binary_search (前提: candidates 已排序)

        for (ui k = 0; k < deg_g; ++k) {
            ui v_curr = neighbors_g[k];

            if (mapped_g[v_curr] != -1) continue; // 必须未被匹配
            if (!binary_search(candidates[u_curr].begin(), candidates[u_curr].end(), v_curr)) continue; // 必须是候选点

            Timer t_check; // Measure verification overhead
            t_check.restart();

            // 验证并计算新的 gamma
            ui new_gamma = gamma;
            bool possible = true;

            // 检查反向边 (Backward Edges)
            for (const auto &bedge : seq.bEdges[h]) {
                ui u_target = (bedge.u == u_curr) ? bedge.v : bedge.u;
                ui v_target = mapped_q[u_target]; // u_target 在 S[0...h-1] 中，必定已匹配

                bool edge_exists = data_graph->hasEdge(v_curr, v_target);

                // R check: 如果边在排除集 R 中，且在数据图中存在，则非法
                // (根据论文定义，R 中的边不应被诱导出来，用于避免重复)
                bool in_R = seq.R.count(bedge.canonical());

                if (in_R && edge_exists) {
                    possible = false;
                    break;
                }

                if (!edge_exists) {
                    if (!in_R) {
                        // 正常缺失，计入 gamma
                        new_gamma++;
                    }
                    else {
                        // 在 R 中且不存在：这是符合预期的，因为我们排除了这条边。
                        // 这里不需要增加 gamma，因为树的生成逻辑已经考虑了这种情况?
                        // 论文 Line 4: if alpha + gamma <= theta. alpha 是 missing back edges.
                        // R 中的边是被“切断”的边，理论上不算做“匹配缺失”，而是“结构变更”。
                        // 但是 Algorithm 2 Line 4 明确计算 alpha。
                        // 如果我们在 R 中，我们显式地不希望这条边匹配。
                        // 如果它确实没匹配，这不算做“missing error”。
                    }
                }
            }
            stats.verify_time += t_check.elapsed();

            if (possible && new_gamma <= threshold) {
                mapped_q[u_curr] = v_curr;
                mapped_g[v_curr] = u_curr;

                SimSearch(h + 1, seq, new_gamma);

                mapped_q[u_curr] = -1;
                mapped_g[v_curr] = -1;
            }
        }

        // --- Phase 2: Alternating-Reordering (横向重组) ---
        // 对应 Algorithm 2 Lines 10-16
        // 尝试替换当前的树边 seq.sEdge[h]，生成新的生成树并继续搜索
        // 条件: 只有当我们还有“额度”去切断一条边时 (gamma < threshold)
        // 注意: 这里的逻辑是将原本的树边放入 R (视为缺失)，然后寻找替代路径

        if (gamma < threshold) {
            Timer t_reorder; // Measure Reordering Overhead
            t_reorder.restart();
            // 尝试替换 sEdge[h]
            QEdge current_tree_edge = seq.sEdge[h];
            QEdge new_edge;

            // 查找是否有一条非树边可以替代 current_tree_edge
            // 并且不在当前的 R 中
            if (findReplacement(seq, h, new_edge)) {
                // 构造新的 seq'
                // 1. 复制当前前缀 S[0...h-1] (因为 h 之前的树结构不变)
                // 2. 将当前树边加入 R
                // 3. 使用 new_edge 替换，并重新运行 Prim 算法生成后缀 S[h...n-1]

                QISequence next_seq;
                next_seq.R = seq.R;
                next_seq.R.insert(current_tree_edge.canonical());

                // 执行 reOrdering 生成新序列
                // 注意：传入 gamma + 1，因为我们刚刚人为切断了一条边 (current_tree_edge)
                // 实际上，Algorithm 2 Line 16 递归调用的是 SimSearch(h, seq', gamma + 1)
                // 这意味着我们在当前层 h，换了一棵树，重新尝试匹配 S'[h]

                if (reorderSequence(seq, h, new_edge, next_seq)) {
                    stats.reorder_time += t_reorder.elapsed(); // Stop timer before recursion

                    SimSearch(h, next_seq, gamma + 1);

                    // Resume timer isn't strictly necessary as we are back in this scope,
                    // but usually we count logic time, not recursion wait time for this block.
                }
                else {
                    stats.reorder_time += t_reorder.elapsed();
                }
            }
            else {
                stats.reorder_time += t_reorder.elapsed();
            }
        }
    }

    // --- Helper Functions for Tree Generation ---

    // 寻找替代边：连接 S[0...h-1] 集合 和 (V - S[0...h-1]) 集合的边
    // 且该边不能是原本的树边，也不能在 R 中
    bool findReplacement(const QISequence &seq, ui h, QEdge &out_edge)
    {
        // 构建已访问集合 (Prefix)
        vector<bool> visited(qn, false);
        for (ui i = 0; i < h; ++i) visited[seq.S[i]] = true;

        // 当前被移除的树边
        QEdge removed = seq.sEdge[h].canonical();

        // 寻找最优替代边 (遵循 Prim 序: 最小权/ID)
        bool found = false;
        QEdge best_e;

        for (const auto &e : all_q_edges) {
            bool u_vis = visited[e.u];
            bool v_vis = visited[e.v];

            // 必须是一个在 visited 中，一个不在 (Cut property)
            if (u_vis != v_vis) {
                QEdge cand = e.canonical();
                // 不能是刚被移除的边
                if (cand == removed) continue;
                // 不能在 R 中
                if (seq.R.count(cand)) continue;

                // 简单的最小字典序策略作为权重
                if (!found) {
                    best_e = cand;
                    found = true;
                }
                else {
                    if (cand < best_e) best_e = cand;
                }
            }
        }

        if (found) out_edge = best_e;
        return found;
    }

    // 基于前缀和新边，重新生成后续序列 (reOrdering)
    bool reorderSequence(const QISequence &old_seq, ui h, QEdge new_edge, QISequence &new_seq)
    {
        stats.reorder_calls++; // Count operations
        // 1. 复制前缀
        new_seq.S.resize(h);
        new_seq.sEdge.resize(h);
        new_seq.bEdges.resize(h); // bEdges 需要重新计算吗？前缀内部的 bEdges 不变

        vector<bool> visited(qn, false);
        for (ui i = 0; i < h; ++i) {
            new_seq.S[i] = old_seq.S[i];
            new_seq.sEdge[i] = old_seq.sEdge[i];
            new_seq.bEdges[i] = old_seq.bEdges[i];
            visited[old_seq.S[i]] = true;
        }

        // 2. 添加第 h 个节点 (由 new_edge 引入)
        ui u_next = (visited[new_edge.u]) ? new_edge.v : new_edge.u;
        if (visited[u_next]) return false; // Should not happen if findReplacement logic is correct

        new_seq.S.push_back(u_next);
        new_seq.sEdge.push_back(new_edge);
        visited[u_next] = true;

        // 计算 h 的 bEdges
        vector<QEdge> current_bEdges;
        ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u_next, deg);
        for (ui k = 0; k < deg; ++k) {
            ui neighbor = nbrs[k];
            if (visited[neighbor]) {
                QEdge e = { u_next, neighbor };
                if (!(e.canonical() == new_edge.canonical())) {
                    current_bEdges.push_back(e);
                }
            }
        }
        new_seq.bEdges.push_back(current_bEdges);

        // 3. 运行 Prim 算法生成剩余部分 (h+1 ... qn-1)
        // 允许的边集合: all_edges - R
        // 注意：已被加入 S 的点都在 visited 中

        for (ui step = h + 1; step < qn; ++step) {
            QEdge best_prim_edge;
            bool found = false;

            // 扫描所有连接 visited 和 unvisited 的边
            // 同样使用简单的遍历 (对于小图足够，大图可用优先队列优化)
            for (const auto &e : all_q_edges) {
                if (new_seq.R.count(e.canonical())) continue; // Skip R

                bool u_vis = visited[e.u];
                bool v_vis = visited[e.v];

                if (u_vis != v_vis) {
                    QEdge cand = e.canonical();
                    if (!found || cand < best_prim_edge) {
                        best_prim_edge = cand;
                        found = true;
                    }
                }
            }

            if (!found) return false; // 图不连通 (由于 R 的切割)

            ui v_new = visited[best_prim_edge.u] ? best_prim_edge.v : best_prim_edge.u;
            new_seq.S.push_back(v_new);
            new_seq.sEdge.push_back(best_prim_edge);
            visited[v_new] = true;

            // 计算 bEdges
            vector<QEdge> b_edges;
            const ui *v_nbrs = query_graph->getVertexNeighbors(v_new, deg);
            for (ui k = 0; k < deg; ++k) {
                ui neighbor = v_nbrs[k];
                if (visited[neighbor]) {
                    QEdge e = { v_new, neighbor };
                    if (!(e.canonical() == best_prim_edge.canonical())) {
                        b_edges.push_back(e);
                    }
                }
            }
            new_seq.bEdges.push_back(b_edges);
        }

        return true;
    }

    // 生成初始序列 (Prim starting from root)
    QISequence generateInitialSequence(ui root)
    {
        QISequence seq;
        seq.S.push_back(root);
        seq.sEdge.push_back({ (ui)-1, (ui)-1 }); // Root has no parent edge
        seq.bEdges.push_back({});

        // 调用 reorderSequence 的一部分逻辑来填充 (实际上就是 Prim)
        // 造一个假的 "0" 长度前缀，利用 reorderSequence 填充整个列表
        // 但这里手写 Prim 更清晰

        vector<bool> visited(qn, false);
        visited[root] = true;

        for (ui i = 1; i < qn; ++i) {
            QEdge best_edge;
            bool found = false;

            for (const auto &e : all_q_edges) {
                bool u_vis = visited[e.u];
                bool v_vis = visited[e.v];
                if (u_vis != v_vis) {
                    if (!found || e.canonical() < best_edge) {
                        best_edge = e.canonical();
                        found = true;
                    }
                }
            }

            if (found) {
                ui next_u = visited[best_edge.u] ? best_edge.v : best_edge.u;
                seq.S.push_back(next_u);
                seq.sEdge.push_back(best_edge);
                visited[next_u] = true;

                vector<QEdge> bes;
                ui deg; const ui *nbrs = query_graph->getVertexNeighbors(next_u, deg);
                for (ui k = 0; k < deg; ++k) {
                    ui nbr = nbrs[k];
                    if (visited[nbr]) {
                        QEdge be = { next_u, nbr };
                        if (!(be.canonical() == best_edge)) bes.push_back(be);
                    }
                }
                seq.bEdges.push_back(bes);
            }
        }
        return seq;
    }

    // --- Utilities ---

    ui selectRoot()
    {
        ui best_u = 0;
        size_t min_cand = candidates[0].size();
        for (ui u = 1; u < qn; ++u) {
            if (candidates[u].size() < min_cand) {
                min_cand = candidates[u].size();
                best_u = u;
            }
        }
        return best_u;
    }

    vector<QEdge> getAllQueryEdges() const
    {
        vector<QEdge> edges;
        edges.reserve(query_graph->getEdgesCount());
        for (ui u = 0; u < qn; ++u) {
            ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (u < v) edges.push_back({ u, v });
            }
        }
        // 排序确保确定性
        sort(edges.begin(), edges.end());
        return edges;
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees)
    {
        ui n = g->getVerticesCount();
        ui num_labels = g->getLabelsCount();
        counts.assign(n, vector<ui>(num_labels, 0));
        degrees.assign(n, 0);
        for (ui u = 0; u < n; ++u) {
            ui deg; const ui *neighbors = g->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                counts[u][g->getVertexLabel(neighbors[i])]++;
                degrees[u]++;
            }
        }
    }

    ui computeDelta(ui u, ui v)
    {
        ui diff = 0;
        size_t sz = Lq_counts[u].size();
        for (size_t i = 0; i < sz; ++i) {
            if (Lq_counts[u][i] > Lg_counts[v][i]) diff += (Lq_counts[u][i] - Lg_counts[v][i]);
        }
        return diff;
    }

    bool calVerticesFilter()
    {
        for (ui u = 0; u < qn; ++u) {
            LabelID label_u = query_graph->getVertexLabel(u);
            for (ui v = 0; v < gn; ++v) {
                if (label_u != data_graph->getVertexLabel(v)) continue;
                if (computeDelta(u, v) <= threshold) {
                    candidates[u].push_back(v);
                }
            }
            if (candidates[u].empty()) return false;
            sort(candidates[u].begin(), candidates[u].end());
        }
        return true;
    }
};

// ============================================================
// Top-level function: Approximate_Matching
// ============================================================
void Approximate_Matching(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
}

void Approximate_TreeSpan(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    TreeSpanSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
}