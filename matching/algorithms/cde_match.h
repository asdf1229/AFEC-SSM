#ifndef MATCHING_ALGORITHMS_CDE_MATCH_H_
#define MATCHING_ALGORITHMS_CDE_MATCH_H_

#include "graph/graph.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
#include <limits>

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
        label_count = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());

        if (qn == 0 || gn == 0) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        max_g_deg = data_graph->getMaxDegree();

        resetState();

        q_matrix.assign(qn, vector<char>(qn, 0));
        q_neighbors.assign(qn, vector<ui>());

        for (ui u = 0; u < qn; ++u) {
            ui deg = 0;
            const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) {
                ui nbr = nbrs[i];
                q_matrix[u][nbr] = 1;
                q_neighbors[u].push_back(nbr);
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

        stats.init_time = t_init.elapsed();
        return true;
    }

    void match(vector<vector<pair<ui, ui>>> &results)
    {
        Timer t_search;
        t_search.restart();

        results_ptr = &results;
        results_ptr->clear();

        ui root = getInitialRoot();

        for (ui v0 : candidates[root]) {
            mapped_q[root] = (int)v0;
            mapped_g[v0] = (int)root;
            in_Mq[root] = 1;
            matched_count = 1;
            part_M.push_back({ root, v0 });

            onVertexMatchStateChanged(root, true);

            dfs(0, (int)root);

            onVertexMatchStateChanged(root, false);

            mapped_q[root] = -1;
            mapped_g[v0] = -1;
            in_Mq[root] = 0;
            matched_count = 0;
            part_M.pop_back();
        }

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct TimeStats {
        long long total_time = 0;
        // init breakdown
        long long init_time = 0;
        long long filter_time = 0;
        // search breakdown
        long long dfs_time = 0;
        long long branch_time = 0;   // candidate enumeration & matching in dfs
        long long lb_time = 0;       // computeLowerBound
        long long frontier_time = 0; // building U_frontier
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
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
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
    ui label_count;
    ui max_g_deg;

    vector<vector<char>> q_matrix;
    vector<vector<ui>> q_neighbors;

    vector<MyBitset> candidates;
    vector<MyBitset> x_cand;

    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;

    // --- filter ---
    vector<vector<ui>> spoke_lb;
    vector<vector<ui>> adjL;       // 用于二分图匹配
    vector<int> match_arr;
    vector<bool> vis_arr;
    MyBitset used_1hop;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<char> in_Mq;
    vector<vector<char>> is_excluded;
    vector<char> in_P;
    ui matched_count;
    vector<pair<ui, ui>> part_M;

    vector<ui> anchor_count;
    vector<int> frontier_pos;
    vector<ui> active_frontier;
    vector<int> lb_match_right;
    vector<ui> lb_seen_right;
    ui lb_seen_token;
    vector<ui> lb_data_frontier_mark;
    ui lb_data_frontier_token;

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

    void onVertexMatchStateChanged(ui u, bool matched)
    {
        ui label_u = (ui)query_graph->getVertexLabel(u);
        if (matched) {
            updateFrontierStatus(u);
        }

        for (ui nbr : q_neighbors[u]) {
            if (!is_excluded[nbr][u]) {
                if (matched) {
                    anchor_count[nbr]++;
                }
                else {
                    anchor_count[nbr]--;
                }
            }
            updateFrontierStatus(nbr);
        }

        if (!matched) {
            updateFrontierStatus(u);
        }
    }

    void resetState()
    {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        in_Mq.assign(qn, 0);
        is_excluded.assign(qn, vector<char>(qn, 0));
        in_P.assign(qn, 0);
        matched_count = 0;
        part_M.clear();
        part_M.reserve(qn);

        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        x_cand.clear();
        x_cand.assign(qn, MyBitset(gn));

        q_matrix.clear();
        q_neighbors.clear();
        spoke_lb.assign(qn, vector<ui>(gn, 0));
        adjL.assign(qn, vector<ui>());
        match_arr.assign(max_g_deg, -1);
        vis_arr.assign(max_g_deg, false);
        used_1hop = MyBitset(gn);

        anchor_count.assign(qn, 0);
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
        lb_match_right.assign(gn, -1);
        lb_seen_right.assign(gn, 0);
        lb_seen_token = 1;
        lb_data_frontier_mark.assign(gn, 0);
        lb_data_frontier_token = 1;

        stats = TimeStats();
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees)
    {
        ui n = g->getVerticesCount();
        counts.assign(n, vector<ui>(label_count, 0));
        degrees.assign(n, 0);

        for (ui u = 0; u < n; ++u) {
            ui deg = 0;
            const ui *neighbors = g->getVertexNeighbors(u, deg);
            degrees[u] = deg;
            for (ui i = 0; i < deg; ++i) {
                LabelID label = g->getVertexLabel(neighbors[i]);
                if (label >= 0 && (ui)label < label_count) {
                    counts[u][(ui)label]++;
                }
            }
        }
    }

    // ========================================================================
    // Filtering
    // ========================================================================
    bool calVerticesFilter()
    {
        // Step 0: LB_NLF
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
        ui degv; const ui *nbrs = data_graph->getVertexNeighbors(v, degv);
        for (ui i = 0; i < S.size(); ++i) {
            adjL[i].clear(); // adjL: bipartite graph
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
        queue<ui> q;
        vector<char> in_q(qn, 1);

        for (ui u = 0; u < qn; ++u) q.push(u);

        while (!q.empty()) {
            ui u = q.front(); q.pop();
            in_q[u] = 0;

            vector<ui> to_remove;
            for (ui v : candidates[u]) {
                ui lb = computeLBSpoke(u, v);
                spoke_lb[u][v] = lb;
                if (lb > threshold) to_remove.push_back(v);
            }

            if (to_remove.empty()) continue;

            for (ui v : to_remove) candidates[u].remove(v);
            if (candidates[u].empty()) return false;

            for (ui v : q_neighbors[u]) {
                if (!in_q[v]) {
                    q.push(v);
                    in_q[v] = 1;
                }
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
            if (S_cand[a].size() != S_cand[b].size()) {
                return S_cand[a].size() < S_cand[b].size();
            }
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
        vector<int> part_S(S.size(), -2); // -2=unprocessed, -1=OUT, (>=0)->data vertex

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

    ui getInitialRoot()
    {
        auto betterRoot = [&](ui a, ui b) -> bool {
            ui cand_a = candidates[a].size();
            ui cand_b = candidates[b].size();
            if (cand_a != cand_b) {
                return cand_a < cand_b;
            }

            ui deg_a = q_neighbors[a].size();
            ui deg_b = q_neighbors[b].size();
            if (deg_a != deg_b) {
                return deg_a > deg_b;
            }

            ui inner_a = countInnerEdgesAmongNeighbors(a);
            ui inner_b = countInnerEdgesAmongNeighbors(b);
            if (inner_a != inner_b) {
                return inner_a > inner_b;
            }

            return a < b;
            };

        ui root = 0;
        for (ui u = 1; u < qn; ++u) {
            if (betterRoot(u, root)) {
                root = u;
            }
        }
        return root;
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

                    for (ui nbr : q_neighbors[curr]) {
                        if (!in_Mq[nbr] && !visited_unmapped[nbr]) {
                            visited_unmapped[nbr] = true;
                            q.push(nbr);
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

    // ========================================================================
    // Lower Bound based Pruning
    // ========================================================================
    ui computeCrossDelta(ui v, const vector<ui> &anchors) const
    {
        ui delta = 0;
        for (ui a : anchors) {
            assert(mapped_q[a] >= 0);
            if (!data_graph->hasEdge(v, (ui)mapped_q[a])) {
                delta++;
            }
        }
        return delta;
    }

    bool findBudgetFeasibleAugment(ui left_idx, const vector<vector<ui>> &adj)
    {
        for (ui v : adj[left_idx]) {
            if (lb_seen_right[v] == lb_seen_token) {
                continue;
            }
            lb_seen_right[v] = lb_seen_token;

            if (lb_match_right[v] < 0 || findBudgetFeasibleAugment((ui)lb_match_right[v], adj)) {
                lb_match_right[v] = (int)left_idx;
                return true;
            }
        }
        return false;
    }

    ui computeBudgetLayerDeficit(const vector<ui> &U_frontier,
        const vector<vector<ui>> &anchors_by_idx,
        ui remaining_budget)
    {
        ui frontier_size = (ui)U_frontier.size();
        if (frontier_size <= 1 || remaining_budget == 0) {
            return 0;
        }

        vector<vector<vector<ui>>> exact_adj(frontier_size, vector<vector<ui>>(remaining_budget));
        vector<ui> right_vertices;

        if (++lb_data_frontier_token == 0) {
            std::fill(lb_data_frontier_mark.begin(), lb_data_frontier_mark.end(), 0);
            lb_data_frontier_token = 1;
        }

        for (ui i = 0; i < frontier_size; ++i) {
            ui u = U_frontier[i];
            const vector<ui> &anchors = anchors_by_idx[i];

            for (ui v : candidates[u]) {
                if (mapped_g[v] != -1 || x_cand[u].contains(v)) {
                    continue;
                }

                ui alpha = computeCrossDelta(v, anchors);
                if (alpha >= remaining_budget) {
                    continue;
                }

                exact_adj[i][alpha].push_back(v);
                if (lb_data_frontier_mark[v] != lb_data_frontier_token) {
                    lb_data_frontier_mark[v] = lb_data_frontier_token;
                    right_vertices.push_back(v);
                }
            }
        }

        vector<vector<ui>> cumulative_adj(frontier_size);
        ui extra_lb = 0;

        for (ui k = 0; k < remaining_budget; ++k) {
            for (ui i = 0; i < frontier_size; ++i) {
                const vector<ui> &delta = exact_adj[i][k];
                cumulative_adj[i].insert(cumulative_adj[i].end(), delta.begin(), delta.end());
            }

            for (ui v : right_vertices) {
                lb_match_right[v] = -1;
            }

            ui mu = 0;
            for (ui i = 0; i < frontier_size; ++i) {
                if (++lb_seen_token == 0) {
                    std::fill(lb_seen_right.begin(), lb_seen_right.end(), 0);
                    lb_seen_token = 1;
                }

                if (findBudgetFeasibleAugment(i, cumulative_adj)) {
                    mu++;
                }
            }

            ui deficit = frontier_size - mu;
            extra_lb += deficit;
            if (extra_lb > remaining_budget) {
                return extra_lb;
            }
        }

        return extra_lb;
    }

    bool shouldPruneByLowerBounds(ui current_miss, const vector<ui> &U_frontier)
    {
        ui remaining_budget = threshold - current_miss;
        if (U_frontier.empty()) {
            return false;
        }

        ui frontier_size = (ui)U_frontier.size();
        vector<vector<ui>> anchors_by_idx(frontier_size);
        ui sum_lb = current_miss;

        for (ui i = 0; i < frontier_size; ++i) {
            ui u = U_frontier[i];
            vector<ui> &anchors = anchors_by_idx[i];
            anchors.clear();
            for (ui nbr : q_neighbors[u]) {
                if (in_Mq[nbr] && !is_excluded[u][nbr]) {
                    anchors.push_back(nbr);
                }
            }

            ui best_boundary = threshold + 1;

            for (ui v : candidates[u]) {
                if (mapped_g[v] != -1 || x_cand[u].contains(v)) {
                    continue;
                }

                ui alpha = computeCrossDelta(v, anchors);
                best_boundary = std::min(best_boundary, alpha);
            }

            if (best_boundary == threshold + 1) {
                return true;
            }

            // Search-time lower bounds must be residual to the current partial match.
            // The static spoke_lb(u, v) is computed before DFS and still accounts for
            // incident query edges that may already have been paid via exclusions in
            // current_miss. Using it directly here can double count those edges and
            // prune valid solutions, so the cheap bound only uses the live anchor cut.
            sum_lb += best_boundary;
            if (sum_lb > threshold) {
                return true;
            }
        }

        ui competition_lb = computeBudgetLayerDeficit(U_frontier, anchors_by_idx, remaining_budget);
        return current_miss + competition_lb > threshold;
    }
    // ========================================================================

    // =====================================================
    // Procedure DFS(M_part, cost, X, P, u_new)
    //
    // cost:  current cost of partial match M_part
    // u_new: the most recently mapped query vertex (-1 means undefined)
    // =====================================================
    void dfs(ui cost, int u_new)
    {
        assert(matched_count > 0);
        assert(matched_count <= qn);
        assert(cost <= threshold);

        if (matched_count == qn) {
            stats.recursion_calls++;
            results_ptr->push_back(part_M);
            return;
        }

        stats.recursion_calls++;

        vector<ui> reEnabledP;
        if (u_new >= 0) {
            for (ui w : q_neighbors[(ui)u_new]) {
                if (in_P[w]) {
                    reEnabledP.push_back(w);
                    in_P[w] = 0;
                    updateFrontierStatus(w);
                }
            }
        }

        Timer t_frontier;

        if (active_frontier.empty()) {
            // stats.prun_calls++;
            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }

        vector<ui> U_frontier = getBestComponentFrontier();
        stats.frontier_time += t_frontier.elapsed();

#ifdef LOWER_BOUND
        Timer t_lb;
        if (shouldPruneByLowerBounds(cost, U_frontier)) {
            stats.lb_time += t_lb.elapsed();
            stats.prun_calls++;

            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }
        stats.lb_time += t_lb.elapsed();
#endif

        Timer t_branch;
        long long child_dfs_time = 0;
        vector<ui> local_P;
        vector<pair<ui, ui>> local_X;
        vector<pair<ui, ui>> local_x_cand;

        ui current_cost = cost;

        for (ui u : U_frontier) {
            vector<ui> U_anchor;
            U_anchor.clear();
            for (ui nbr : q_neighbors[u]) {
                if (in_Mq[nbr]) {
                    U_anchor.push_back(nbr);
                }
            }

            bool threshold_exceeded = false;

            // --- matching u 按 anchor 进行细分 ---
            for (ui ua : U_anchor) {
                if (is_excluded[u][ua]) continue;

                // branch 1: matching u by (u, ua) edge
                vector<ui> anchor_v_list;
                ui deg; const ui *nbrs = data_graph->getVertexNeighbors(mapped_q[ua], deg);
                for (ui j = 0; j < deg; ++j) {
                    ui v = nbrs[j];
                    if (candidates[u].contains(v)) anchor_v_list.push_back(v);
                }

                for (ui v : anchor_v_list) {
                    if (mapped_g[v] != -1) continue;
                    if (x_cand[u].contains(v)) continue;

                    ui delta = 0;

                    for (ui other_ua : U_anchor) {
                        if (other_ua == ua) continue;
                        if (is_excluded[u][other_ua]) continue;

                        bool has_edge = data_graph->hasEdge(v, mapped_q[other_ua]);

                        if (!has_edge) delta++;
                    }

                    if (current_cost + delta > threshold) continue;

                    mapped_q[u] = (int)v;
                    mapped_g[v] = (int)u;
                    in_Mq[u] = 1;
                    matched_count++;
                    part_M.push_back({ u, v });

                    onVertexMatchStateChanged(u, true);

                    Timer t_child;
                    dfs(current_cost + delta, (int)u);
                    child_dfs_time += t_child.elapsed();

                    onVertexMatchStateChanged(u, false);

                    part_M.pop_back();
                    matched_count--;
                    in_Mq[u] = 0;
                    mapped_g[v] = -1;
                    mapped_q[u] = -1;
                }

                // branch 2 : excluding (u, ua)
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
    // ========================================================================
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
