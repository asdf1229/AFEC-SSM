#ifndef MATCHING_ALGORITHMS_CDE_MATCH_H_
#define MATCHING_ALGORITHMS_CDE_MATCH_H_

#include "graph/graph.h"
#include "utility/utility.h"
#include "utility/mybitset.h"

using namespace std;

// ============================================================================
// MatchingSolver Implementation
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

        // 获取最大数据图度数，用于预分配内存
        max_g_deg = 0;
        for (ui i = 0; i < gn; ++i) {
            ui d; data_graph->getVertexNeighbors(i, d);
            if (d > max_g_deg) max_g_deg = d;
        }

        resetState();

        // 预处理 q_matrix 和 q_neighbors
        q_matrix.assign(qn, vector<bool>(qn, false));
        for (ui u = 0; u < qn; ++u) {
            ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                q_matrix[u][nbrs[i]] = true;
                q_neighbors[u].push_back(nbrs[i]);
            }
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

        // MyBitset 支持基于 Iterator 的 range-based for 循环
        for (ui v0 : candidates[v]) {
            mapped_q[v] = v0;
            mapped_g[v0] = v;
            in_Mq[v] = 1;
            matched_count = 1;
            part_M.push_back({ v, v0 }); // v -> v0

            updateFrontierStatus(v);

            // 用邻接矩阵 q_matrix 替代 getVertexNeighbors 遍历
            for (ui i = 0; i < qn; ++i) {
                if (q_matrix[v][i]) {
                    if (!is_excluded[i][v]) {
                        anchor_count[i]++;
                    }
                    updateFrontierStatus(i);
                }
            }

            dfs(0, (int)v);

            for (ui i = 0; i < qn; ++i) {
                if (q_matrix[v][i]) {
                    if (!is_excluded[i][v]) {
                        anchor_count[i]--;
                    }
                    updateFrontierStatus(i);
                }
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
        long long shell_time = 0;    // shell batch processing time
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        printf("\n--- CDE-Match Time Analysis ---\n");
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
    ui max_g_deg;

    vector<vector<bool>> q_matrix;
    vector<vector<ui>> q_neighbors;

    // 修改处 1：将 candidates 声明为 vector<MyBitset>
    vector<MyBitset> candidates;
    // 全局 cand 排除集，使用 MyBitset 保证 O(1) 插入和查询
    vector<MyBitset> x_cand;

    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;

    // --- 过滤用辅助存储 ---
    vector<vector<ui>> spoke_lb;
    vector<vector<ui>> adjL;       // 用于二分图匹配
    vector<int> match_arr;
    vector<bool> vis_arr;
    MyBitset used_1hop;            // 用于1Hop搜索

    // --- CDE & Dynamic KSS Core States ---
    vector<int> mapped_q;            // f: query vertex → data vertex (-1 if unmapped)
    vector<int> mapped_g;            // reverse: data vertex → query vertex (-1 if unused)
    vector<char> in_Mq;              // 当前已匹配集合 M_q
    vector<vector<int>> is_excluded; // X: 排除边集合 (symmetric, is_excluded[u][v] > 0 means (u,v) ∈ X)
    vector<char> in_P;               // P: 延期顶点集合
    ui matched_count;
    vector<pair<ui, ui>> part_M;

    vector<ui> _order;               // 固定搜索顺序 π
    vector<ui> order_rank;           // 用于还原原始匹配顺序的排名

    // --- Incremental U_frontier States ---
    vector<ui> anchor_count;         // 记录每个顶点当前有多少个未排斥(un-excluded)的已匹配邻居
    vector<int> frontier_pos;        // 记录顶点在 active_frontier 中的位置，-1表示不在其中
    vector<ui> active_frontier;      // 动态维护的 U_frontier 集合

    // O(1) 增量更新顶点 u 的 Frontier 状态
    inline void updateFrontierStatus(ui u)
    {
        bool should_be = (!in_Mq[u] && !in_P[u] && anchor_count[u] > 0);
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

        // 修改处 2：初始化 candidates，为其分配 gn（数据图顶点数）范围的 MyBitset
        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        x_cand.clear();
        x_cand.assign(qn, MyBitset(gn));

        q_neighbors.assign(qn, vector<ui>());
        spoke_lb.assign(qn, vector<ui>(gn, 0));
        adjL.assign(qn, vector<ui>());
        match_arr.assign(max_g_deg, -1);
        vis_arr.assign(max_g_deg, false);
        used_1hop = MyBitset(gn);

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

    // ========================================================================
    // 过滤模块重构：三段式过滤
    // ========================================================================
    bool calVerticesFilter()
    {
        // Step 0: coarse filter
        if (!filterNLF()) return false;

        // Step 1: LB_spoke
        if (!filterSpoke()) return false;

        // Step 2: LB_1hop
        if (!filterOneHop()) return false;

#ifndef NDEBUG
        // statistics
        std::vector<ui> missing_edges_dist(threshold + 1, 0);
        std::vector<std::vector<ui>> vertex_missing_edges_dist(qn, std::vector<ui>(threshold + 1, 0));
        ui total_candidates_count = 0;

        for (ui u = 0; u < qn; ++u) {
            total_candidates_count += candidates[u].size();

            for (ui v : candidates[u]) {
                ui min_missing_edges = computeLBSpoke(u, v);

                if (min_missing_edges <= threshold) {
                    missing_edges_dist[min_missing_edges]++;
                    vertex_missing_edges_dist[u][min_missing_edges]++;
                }
            }
        }

        printf("\n================ Candidate Missing Edges Statistics ================\n");
        printf("Total valid candidates across all query vertices: %u\n", total_candidates_count);
        for (ui i = 0; i <= threshold; ++i) {
            double percent = (total_candidates_count == 0) ? 0.0 : (double)missing_edges_dist[i] / total_candidates_count * 100.0;
            printf("Missing edges = %u: %6u candidates (%6.2f %%)\n", i, missing_edges_dist[i], percent);
        }
        printf("====================================================================\n\n");

        printf("candidates nums and missing edges distribution:\n");
        for (ui u = 0; u < qn; ++u) {
            ui cand_size = candidates[u].size();
            ui deg_u = q_neighbors[u].size();
            ui max_miss_allowed = (deg_u > 0) ? std::min(threshold, deg_u - 1) : threshold;

            printf("u = %u: %6d candidates", u, cand_size);

            if (cand_size > 0) {
                printf("  [ ");
                for (ui i = 0; i <= max_miss_allowed; ++i) {
                    double percent = (double)vertex_missing_edges_dist[u][i] / cand_size * 100.0;
                    printf("%u-miss: %5.1f%%", i, percent);
                    if (i < max_miss_allowed) printf(" | ");
                }
                printf(" ]");
            }
            printf("\n");
        }
#endif

        return true;
    }

    // --- NLF filter ---
    bool filterNLF()
    {
        for (ui u = 0; u < qn; ++u) {
            LabelID lu = query_graph->getVertexLabel(u);
            for (ui v = 0; v < gn; ++v) {
                if (lu != data_graph->getVertexLabel(v)) continue;
                if (Lq_degrees[u] > Lg_degrees[v] + threshold) continue;
                if (computeNLF(u, v) > threshold) continue;
                candidates[u].insert(v);
            }
            if (candidates[u].empty()) return false;
        }
        return true;
    }

    // --- Spoke filter ---
    bool dfs_match(ui u, const vector<vector<ui>> &adj, vector<int> &match, vector<bool> &vis)
    {
        for (ui right_node : adj[u]) {
            if (vis[right_node]) continue;
            vis[right_node] = true;
            if (match[right_node] < 0 || dfs_match((ui)match[right_node], adj, match, vis)) {
                match[right_node] = (int)u;
                return true;
            }
        }
        return false;
    }

    ui hungarian(const vector<vector<ui>> &adj, ui left_size, ui right_size)
    {
        ui mu = 0;
        std::fill(match_arr.begin(), match_arr.begin() + right_size, -1);
        for (ui i = 0; i < left_size; ++i) {
            std::fill(vis_arr.begin(), vis_arr.begin() + right_size, false);
            if (dfs_match(i, adj, match_arr, vis_arr)) mu++;
        }
        return mu;
    }

    ui computeLBSpoke(ui u, ui v)
    {
        const vector<ui> &S = q_neighbors[u];
        assert(!S.empty());
        ui degv; const ui *nbrs = data_graph->getVertexNeighbors(v, degv);
        for (ui i = 0; i < S.size(); ++i) {
            adjL[i].clear();
            ui u1 = S[i];
            for (ui j = 0; j < degv; ++j) {
                ui v1 = nbrs[j];
                if (candidates[u1].contains(v1)) adjL[i].push_back(j);
            }
        }
        ui mu = hungarian(adjL, (ui)S.size(), degv);
        return (ui)S.size() - mu;
    }

    bool filterSpoke()
    {
        bool changed = true;
        while (changed) {
            changed = false;
            for (ui u = 0; u < qn; ++u) {
                vector<ui> to_remove;
                for (ui v : candidates[u]) {
                    ui lb = computeLBSpoke(u, v);
                    spoke_lb[u][v] = lb;
                    if (lb > threshold) to_remove.push_back(v);
                }
                if (!to_remove.empty()) {
                    for (ui v : to_remove) candidates[u].remove(v);
                    changed = true;
                }
                if (candidates[u].empty()) return false;
            }
        }
        return true;
    }

    // --- One-hop filter ---
    ui countInnerEdgesAmongNeighbors(ui u)
    {
        ui count = 0;
        const vector<ui> &S = q_neighbors[u];
        for (ui i = 0; i < S.size(); ++i) {
            for (ui j = i + 1; j < S.size(); ++j) {
                if (q_matrix[S[i]][S[j]]) count++;
            }
        }
        return count;
    }

    vector<ui> orderDfsOneHop(const vector<ui> &S, const vector<vector<ui>> &S_cand)
    {
        vector<ui> ord(S.size());
        iota(ord.begin(), ord.end(), 0);

        sort(ord.begin(), ord.end(), [&](ui a, ui b) {
            // inside 域小的优先
            if (S_cand[a].size() != S_cand[b].size()) {
                return S_cand[a].size() < S_cand[b].size();
            }
            // 按邻居子图度数高的优先
            ui deg_a = 0, deg_b = 0;
            for (ui z : S) {
                if (q_matrix[S[a]][z]) deg_a++;
                if (q_matrix[S[b]][z]) deg_b++;
            }
            return deg_a > deg_b;
            });
        return ord;
    }

    ui residualSpokeLB(ui pos, const vector<ui> &ord, const vector<vector<ui>> &S_cand)
    {
        // 简单安全的剩余代价下界：如果后续点所有可能的 S_cand 元素都被 used，
        // 或者 S_cand 本来就为空，则必定走 OUT 分支，必定产生缺失。
        ui rem = 0;
        for (ui k = pos; k < ord.size(); ++k) {
            ui i = ord[k];
            bool has_free = false;
            for (ui y : S_cand[i]) {
                if (!used_1hop.contains(y)) {
                    has_free = true;
                    break;
                }
            }
            if (!has_free) rem++;
        }
        return rem;
    }

    bool dfsOneHop(ui pos, const vector<ui> &S, const vector<vector<ui>> &S_cand, vector<int> &part_S, const vector<ui> &ord, ui current_cost)
    {

        if (current_cost > threshold) return false;
        if (pos == ord.size()) return true;

        ui rem_lb = residualSpokeLB(pos, ord, S_cand);
        if (current_cost + rem_lb > threshold) return false;

        ui i = ord[pos];
        ui x = S[i];

        // branch 1：x -> OUT
        part_S[i] = -1;
        if (dfsOneHop(pos + 1, S, S_cand, part_S, ord, current_cost + 1)) {
            part_S[i] = -2;
            return true;
        }

        // branch 2：x -> inside y
        for (ui y : S_cand[i]) {
            if (used_1hop.contains(y)) continue;

            ui delta_inner = 0;

            for (ui j = 0; j < S.size(); ++j) {
                if (part_S[j] < 0) continue;
                ui z = S[j];
                if (!q_matrix[x][z]) continue;

                ui yz = (ui)part_S[j];
                if (!data_graph->hasEdge(y, yz)) {
                    delta_inner++;
                }
            }

            if (current_cost + delta_inner > threshold) continue;

            used_1hop.insert(y);
            part_S[i] = (int)y;

            if (dfsOneHop(pos + 1, S, S_cand, part_S, ord, current_cost + delta_inner)) {
                used_1hop.remove(y);
                part_S[i] = -2;
                return true;
            }

            used_1hop.remove(y);
            part_S[i] = -2;
        }

        part_S[i] = -2;
        return false;
    }

    bool checkLBOneHop(ui u, ui v)
    {
        const vector<ui> &S = q_neighbors[u];

        vector<vector<ui>> S_cand(S.size());
        ui degv; const ui *nbrs = data_graph->getVertexNeighbors(v, degv);

        for (ui i = 0; i < S.size(); ++i) {
            ui u1 = S[i];
            for (ui j = 0; j < degv; ++j) {
                ui y = nbrs[j];
                if (candidates[u1].contains(y)) S_cand[i].push_back(y);
            }
        }

        vector<ui> ord = orderDfsOneHop(S, S_cand);
        vector<int> part_S(S.size(), -2); // -2=未处理, -1=OUT, >=0 表示数据点

        return dfsOneHop(0, S, S_cand, part_S, ord, 0);
    }

    bool filterOneHop()
    {
        vector<pair<ui, ui>> del;

        for (ui u = 0; u < qn; ++u) {
            ui du = q_neighbors[u].size();
            if (du <= 1) continue;
            if (countInnerEdgesAmongNeighbors(u) == 0) continue;

            for (ui v : candidates[u]) {
                ui lb = spoke_lb[u][v];

                // TODO
                // if (lb != threshold) continue;

                if (!checkLBOneHop(u, v)) {
                    del.push_back({ u, v });
                }
            }
        }

        for (auto p : del) {
            candidates[p.first].remove(p.second);
        }

        for (ui u = 0; u < qn; ++u) {
            if (candidates[u].empty()) return false;
        }
        return true;
    }
    // ========================================================================

    void generateMatchingOrder()
    {
        _order.clear();
        vector<bool> visited(qn, false);

        ui root = 0;
        int min_cand = candidates[0].size();
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
                for (ui i = 0; i < qn; ++i) {
                    if (q_matrix[u][i] && visited[i]) current_anchor_sz++;
                }

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
                        if (current_anchor_sz > best_anchor_sz) {
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

    vector<ui> getBestComponentFrontier()
    {
        // 1. 遍历一遍，判断当前查询图被 Mq 分割为了几个不相连的部分
        vector<bool> visited_unmapped(qn, false);
        vector<vector<ui>> components;

        for (ui i = 0; i < qn; ++i) {
            if (!in_Mq[i] && !visited_unmapped[i]) {
                vector<ui> comp;
                queue<ui> q;
                q.push(i);
                visited_unmapped[i] = true;

                while (!q.empty()) {
                    ui curr = q.front();
                    q.pop();
                    comp.push_back(curr);

                    for (ui nbr = 0; nbr < qn; ++nbr) {
                        if (q_matrix[curr][nbr]) {
                            // 仅考虑未映射点之间的连通性
                            if (!in_Mq[nbr] && !visited_unmapped[nbr]) {
                                visited_unmapped[nbr] = true;
                                q.push(nbr);
                            }
                        }
                    }
                }
                components.push_back(comp);
            }
        }

        // 2. 选择其中（点数 * 每个点的待选集）最少的部分，以及对应部分的 U_frontier
        vector<ui> best_U_frontier;

        // 使用 double 替代 unsigned long long 来存储代价，彻底避免乘法溢出风险
        // 注意：需要包含头文件 <limits>
        double min_cost = std::numeric_limits<double>::max();

        for (const auto &comp : components) {
            vector<ui> comp_frontier;
            double sum_cand = 1.0;

            for (ui v : comp) {
                sum_cand *= candidates[v].size();
                // 检查 v 是否属于当前的 active_frontier (O(1)判断)
                if (frontier_pos[v] != -1) {
                    comp_frontier.push_back(v);
                }
            }

            // 仅对具有合法 Frontier 的连通分量进行代价评估
            assert(!comp_frontier.empty());
            if (!comp_frontier.empty()) {
                if (sum_cand < min_cost) {
                    min_cost = sum_cand;
                    best_U_frontier = comp_frontier;
                }
            }
        }

        // 3. 返回找到的最优 U_frontier
        return best_U_frontier;
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
            for (ui w = 0; w < qn; ++w) {
                if (q_matrix[u_new][w] && in_P[w]) {
                    reEnabledP.push_back(w);
                    in_P[w] = 0;
                    updateFrontierStatus(w);
                }
            }
        }

        Timer t_frontier;

        // 提前拦截：如果当前没有任何 Frontier 点，则直接剪枝回溯
        if (active_frontier.empty()) {
            stats.prun_calls++;
            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }

        vector<ui> U_frontier = getBestComponentFrontier();
        // TODO dynamic component-based frontier selection
        // TODO choose best u
        sort(U_frontier.begin(), U_frontier.end(), [&](ui a, ui b) {
            return order_rank[a] < order_rank[b];
            });
        stats.frontier_time += t_frontier.elapsed();

        Timer t_branch;
        long long child_dfs_time = 0;
        vector<ui> local_P;
        vector<pair<ui, ui>> local_X;

        // 用于在当前 DFS 状态结束时回溯 x_cand 排除集
        vector<pair<ui, ui>> local_x_cand;

        ui current_cost = cost;

        for (ui u : U_frontier) {
            vector<ui> U_anchor;

            for (ui i = 0; i < qn; ++i) {
                if (q_matrix[u][i] && in_Mq[i]) {
                    U_anchor.push_back(i);
                }
            }

            bool threshold_exceeded = false;

            // --- matching u 按 anchor 进行细分 ---
            for (ui ua : U_anchor) {
                if (is_excluded[u][ua]) continue;

                // 第一阶段：精准收集所有与该 anchor 连通的候选点
                vector<ui> anchor_v_list;
                ui deg; const ui *nbrs = data_graph->getVertexNeighbors(mapped_q[ua], deg);
                for (ui j = 0; j < deg; ++j) {
                    ui v = nbrs[j];
                    // O(1) 利用 Bitset 检查 v 是否为 u 的候选集
                    if (candidates[u].contains(v)) {
                        anchor_v_list.push_back(v);
                    }
                }

                // 第二阶段：第一分支 - 从这些候选点中尝试扩展
                for (ui v : anchor_v_list) {
                    // 如果已经被映射过，或者已经进入了全局/分支排除集，则直接跳过
                    if (mapped_g[v] != -1) continue;
                    if (x_cand[u].contains(v)) continue;

                    bool conflict = false;
                    ui delta = 0;

                    for (ui other_ua : U_anchor) {
                        if (other_ua == ua) continue;
                        if (is_excluded[u][other_ua]) continue;

                        bool has_edge = data_graph->hasEdge(v, mapped_q[other_ua]);

                        if (!has_edge) delta++;
                    }

                    if (conflict || current_cost + delta > threshold) continue;

                    mapped_q[u] = v;
                    mapped_g[v] = u;
                    in_Mq[u] = 1;
                    matched_count++;
                    part_M.push_back({ u, v });

                    updateFrontierStatus(u);

                    for (ui i = 0; i < qn; ++i) {
                        if (q_matrix[u][i]) {
                            if (!is_excluded[i][u]) {
                                anchor_count[i]++;
                            }
                            updateFrontierStatus(i);
                        }
                    }

                    Timer t_child;
                    dfs(current_cost + delta, (int)u);
                    child_dfs_time += t_child.elapsed();

                    for (ui i = 0; i < qn; ++i) {
                        if (q_matrix[u][i]) {
                            if (!is_excluded[i][u]) {
                                anchor_count[i]--;
                            }
                            updateFrontierStatus(i);
                        }
                    }

                    part_M.pop_back();
                    matched_count--;
                    in_Mq[u] = 0;
                    mapped_g[v] = -1;
                    mapped_q[u] = -1;

                    updateFrontierStatus(u);
                }

                // 第三阶段：第二分支 - 排斥该边 (u, ua)
                current_cost++;
                is_excluded[u][ua] = 1;
                is_excluded[ua][u] = 1;
                anchor_count[u]--;
                local_X.push_back({ u, ua });

                if (current_cost > threshold) {
                    threshold_exceeded = true;
                    break;
                }

                // 一旦该边被加入排斥集，意味着任何与 mapped_q[ua] 连通的候选点 v 都不可能作为 u 的映射。
                // 我们将上面完整收集的 anchor_v_list 全部加入 x_cand 排斥集。
                for (ui v : anchor_v_list) {
                    if (!x_cand[u].contains(v)) {
                        x_cand[u].insert(v);
                        local_x_cand.push_back({ u, v });
                    }
                }
            }

            // 如果还没能把 u 塞进 P 集合就已经超过阈值，整条线死路一条，直接 break 并回溯。
            if (threshold_exceeded) {
                break;
            }

            // 当所有 anchor 都断开（且容错尚未爆表），证明该 u 与当前任意锚点都没连通，延期 u
            in_P[u] = 1;
            local_P.push_back(u);
            updateFrontierStatus(u);
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

        // 精准清空当前 DFS 分支中加入的排斥集候选点（基于局部回溯记录）
        for (auto &p : local_x_cand) {
            x_cand[p.first].remove(p.second);
        }

        for (ui w : reEnabledP) {
            in_P[w] = 1;
            updateFrontierStatus(w);
        }
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

#endif
