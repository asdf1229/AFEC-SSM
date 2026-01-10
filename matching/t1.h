#include "utility/utility.h"
// #include "graph/graph.h"

using namespace std;

class MatchingSolver{
public:
    MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}
    bool init(const Graph *q, const Graph *g, ui match_threshold)
    {
        // 1. Initialization
        query_graph = q;
        data_graph = g;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();
        resetState();

        // 2. Initial Setup
        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

        bool res = calVerticesFilter();
        if(res == false) {
            printf("!!! No candidates for some query vertex. No possible mapping. !!!\n");
            assert(false);
            return false;
        }
        return true;
    }
    void match(vector<vector<pair<ui, ui>>> &results) 
    {
        results_ptr = &results;
        vector<pair<ui, ui>> initial_Mcand; // Empty for root

        // 3. Start Recursion
        Timer t_dfs;
        dfs(0, (ui)-1, (ui)-1, initial_Mcand);
        stats.dfs_time += t_dfs.elapsed();
    }

    // --- Stats ---
    struct TimeStats {
        long long total_time = 0;
        long long dfs_time = 0;         // DFS 耗时
        long long get_cand_time = 0;    // ExpandableMappings 耗时
        long long lb_time = 0;          // calLowerBound 耗时
        long long branch_time = 0;      // 分支时间
        long long recursion_calls = 0;  // 递归调用次数计数
        long long prun_calls = 0;       // 剪枝次数

        // calLowerBound 内部细分时间
        long long lb_identify_sets_time = 0; // 1. Identify Nq, Rq, Ng, Rg
        long long lb_remain_time = 0; // 2. Calculate LB_remain
        long long lb_front_time = 0; // 3. Calculate LB_front
        long long lb_mwpm_time = 0; // MWPM 求解时间
    }stats;

    void printStats() const {
        printf("\n--- Matching Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("DFS Time:            %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, (double)stats.dfs_time/stats.total_time*100);
        printf("- Mcand Time:        %.4lf ms (%.2f%% of DFS)\n", stats.get_cand_time / 1000.0, (double)stats.get_cand_time/stats.dfs_time*100);
        printf("- Branch Time:       %.4lf ms (%.2f%% of DFS)\n", stats.branch_time / 1000.0, (double)stats.branch_time/stats.dfs_time*100);
        printf("- LowerBound Time:   %.4lf ms (%.2f%% of DFS)\n", stats.lb_time / 1000.0, (double)stats.lb_time/stats.dfs_time*100);
        printf("  - Identify Sets:   %.4lf ms (%.2f%% of LB)\n", stats.lb_identify_sets_time / 1000.0, (double)stats.lb_identify_sets_time/stats.lb_time*100);
        printf("  - LB Remain:       %.4lf ms (%.2f%% of LB)\n", stats.lb_remain_time / 1000.0, (double)stats.lb_remain_time/stats.lb_time*100);
        printf("  - LB Front (MWPM): %.4lf ms (%.2f%% of LB)\n", stats.lb_front_time / 1000.0, (double)stats.lb_front_time/stats.lb_time*100);
        printf("    - MWPM:          %.4lf ms (%.2f%% of LB Front)\n", stats.lb_mwpm_time / 1000.0, (double)stats.lb_mwpm_time/stats.lb_front_time*100);
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
        printf("Results Found:       %zu\n", results_ptr ? results_ptr->size() : 0);
        printf("---------------------------------\n");
    }
private:
    // --- Input Data ---
    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui> > > *results_ptr; // Results
    ui threshold;
    ui qn, gn;

    vector<vector<ui>> candidates; // Pointer to external candidates

    // --- Global State (Modified during DFS) ---
    vector<int> mapped_q; // Query vertex -> Data vertex
    vector<int> mapped_g; // Data vertex -> Query vertex
    vector<pair<ui, ui>> part_M; // Current partial mapping
    unordered_set<pair<ui, ui>, PairHash> X; // Exclusion set
    

    // --- Dynamic Lookahead / Filtering Stats ---
    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;

    // --- OPTIMIZATION: Incremental Cache State ---
    vector<int> q_last_update; // Timestamp when Lq_counts[u] changed
    vector<int> g_last_update; // Timestamp when Lg_counts[v] changed
    int global_time;           // Increments on every update
    vector<vector<ui>> delta_cache;        // Cache [u][v] -> cost
    vector<vector<int>> delta_cache_time;  // Cache [u][v] -> valid timestamp

    // --- Lower Bound Optimization & Memoization ---
    vector<ui> min_costs;
    vector<ui> best_cands;
    vector<int> g_mark;
    int mark_token;

    // --- Buffers for LB ---

    // Members for Identify Sets
    vector<int> lb_in_Nq, lb_in_Ng;
    vector<vector<ui>> lb_adj_Mcand;
    int lb_token;
    
    // Members for MWPM (KM Algorithm) to avoid malloc in loop
    vector<int> km_lx, km_ly, km_mx, km_my;
    vector<int> km_slack, km_slackmy;
    vector<int> km_prev, km_queue;
    vector<char> km_visX, km_visY;
    
    // Members for LB Front
    vector<int> lb_v_to_col;

    // State Helpers
    void resetState() {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        part_M.clear();
        part_M.reserve(qn);
        X.clear();
        if (results_ptr != nullptr) results_ptr->clear();

        // Reset Candidates
        candidates.clear();
        candidates.resize(qn);

        min_costs.assign(qn, (ui)-1);
        best_cands.assign(qn, (ui)-1);
        g_mark.assign(gn, 0);
        mark_token = 0;

        // Reset LB buffers
        lb_in_Nq.assign(qn, 0);
        lb_in_Ng.assign(gn, 0);
        lb_adj_Mcand.assign(qn, vector<ui>());
        lb_token = 0;
        lb_v_to_col.assign(gn, -1);

        // Reset Cache Structures
        q_last_update.assign(qn, 0);
        g_last_update.assign(gn, 0);
        global_time = 0;
        
        // Initialize cache matrix if needed (preserves memory across calls if object reused, but resets content logic)
        if (delta_cache.size() != qn) {
            delta_cache.assign(qn, vector<ui>(gn, 0));
            delta_cache_time.assign(qn, vector<int>(gn, 0));
        } else {
            // Fast clear by just resetting timestamps to 0
            for(auto &row : delta_cache_time) fill(row.begin(), row.end(), 0);
        }

        // Reset Stats
        stats = TimeStats();
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees)
    {
        ui n = g->getVerticesCount();
        ui num_labels = g->getLabelsCount();
        counts.assign(n, vector<ui>(num_labels, 0));
        degrees.assign(n, 0);

        for (ui u = 0; u < n; ++u) {
            ui deg; const ui* neighbors = g->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = neighbors[i];
                LabelID l_v = g->getVertexLabel(v);
                counts[u][l_v]++;
                degrees[u]++;
            }
        }
    }

    bool calVerticesFilter()
    {
        for (ui u = 0; u < qn; ++u) {
            LabelID label_u = query_graph->getVertexLabel(u);
            for (ui v = 0; v < gn; ++v) {
                if (label_u != data_graph->getVertexLabel(v)) continue;
                ui delta = computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
                assert(delta_cache.size() == qn && delta_cache[u].size() == gn && global_time == 0);
                delta_cache[u][v] = delta;
                delta_cache_time[u][v] = global_time;
                if (delta <= threshold) candidates[u].push_back(v);
            }
            if (candidates[u].empty())  return false;
        }
        return true;
    }

    void updateNeighborCounts(const Graph *g, ui vertex, const vector<int> &mapped_status, 
                              vector<vector<ui>> &counts, vector<ui> &degrees, bool is_add, vector<int> &timestamps)
    {
        LabelID label = g->getVertexLabel(vertex);
        ui deg; const ui* neighbors = g->getVertexNeighbors(vertex, deg);
        for (ui i = 0; i < deg; ++i) {
            ui neighbor = neighbors[i];
            if (mapped_status[neighbor] == -1) {
                if (is_add) {
                    counts[neighbor][label]++;
                    degrees[neighbor]++;
                }
                else {
                    counts[neighbor][label]--;
                    degrees[neighbor]--;
                }
                // [MOD] Mark this neighbor as modified at current global time
                timestamps[neighbor] = ++global_time;
            }
        }
    }

    // Missing Edges & Delta

    ui computeDeltaInnerDoubledVector(const vector<ui> &Lq, const vector<ui> &Lg)
    {
        ui diff = 0;
        for (size_t i = 0; i < Lq.size(); ++i) if (Lq[i] > Lg[i]) diff += (Lq[i] - Lg[i]);
        return diff;
    }

    // [MOD]: Memoized Delta Calculation
    inline ui getDelta(ui u, ui v) {
        // If the cache is newer than the last modification of either u or v, use it
        if (delta_cache_time[u][v] >= q_last_update[u] && 
            delta_cache_time[u][v] >= g_last_update[v]) {
            return delta_cache[u][v];
        }

        // Otherwise compute, store, and update timestamp
        ui val = computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
        delta_cache[u][v] = val;
        // The cache is valid for the current global_time
        delta_cache_time[u][v] = global_time; 
        return val;
    }
    
    ui calIncrementalMissing(ui new_u, ui new_v)
    {
        ui delta_missing = 0;
        for (const auto &p : part_M) {
            ui u_prev = p.first;
            ui v_prev = p.second;
            if (query_graph->hasEdge(new_u, u_prev) && !data_graph->hasEdge(new_v, v_prev)) delta_missing++;
        }
        return delta_missing;
    }

    // Candidate Generation

    void getIncrementalExpandableMappings(ui new_u, ui new_v, const vector<pair<ui, ui>> &fa_Mcand, vector<pair<ui,ui> > &Mcand)
    {
        Mcand.clear();

        if (part_M.empty()) {
            assert(new_u == (ui)-1 && new_v == (ui)-1);
            // select 0 as starting vertex
            ui u = 0;
            assert(mapped_q[u] == -1);
            for(ui v : candidates[u]) {
                assert(mapped_g[v] == -1 && !X.count({u, v}));
                Mcand.emplace_back(u, v);
            }
            return;
        }

        Mcand.reserve(fa_Mcand.size() + 64);

        // 1. Inherit from parent
        for (const auto &p : fa_Mcand) {
            if (p.first == new_u || p.second == new_v) continue;
            if (X.count(p)) continue;
            Mcand.push_back(p);
        }

        ui q_count; const ui* q_neighbors = query_graph->getVertexNeighbors(new_u, q_count);
        ui g_count; const ui* g_neighbors = data_graph->getVertexNeighbors(new_v, g_count);
        size_t inherit_end = Mcand.size();
        static vector<bool> local_g_mark;
        if (local_g_mark.size() < gn) local_g_mark.resize(gn, false);

        // 2. Generate new candidates based on neighbors
        for(ui i = 0; i < g_count; ++i) local_g_mark[g_neighbors[i]] = true;
        for(ui i = 0; i < q_count; ++i) {
            ui u2 = q_neighbors[i];
            if(mapped_q[u2] != -1) continue;

            for(ui v2 : candidates[u2]) {
                if(mapped_g[v2] != -1) continue;
                if(local_g_mark[v2]) {
                    if(X.count({u2, v2})) continue;
                    Mcand.emplace_back(u2, v2);
                }
            }
        }
        for(ui i = 0; i < g_count; ++i) local_g_mark[g_neighbors[i]] = false;

        // 3. Deduplicate
        if (Mcand.size() > inherit_end) {
            sort(Mcand.begin() + inherit_end, Mcand.end());
            if (inherit_end > 0) inplace_merge(Mcand.begin(), Mcand.begin() + inherit_end, Mcand.end());
            Mcand.erase(unique(Mcand.begin(), Mcand.end()), Mcand.end());
        }
    }

    /**
     * @brief 求解最小权重完美匹配 (Minimum Weight Perfect Matching)
     * 逻辑参考自 Application.txt 中的 Hungarian 实现。
     * 
     * @param cost_matrix 代价矩阵，行表示查询图节点，列表示数据图节点
     * @return 最小权重之和
     */
    ui solveMWPM(const vector<vector<ui>>& cost_matrix, ui limit) {
        if (cost_matrix.empty()) return 0;
        int n = (int)cost_matrix.size();
        int m = (int)cost_matrix[0].size();
        const int INF_INT = 1e9;

        // Resize buffers if necessary
        if (km_lx.size() < (size_t)n) { km_lx.resize(n); km_mx.resize(n); km_prev.resize(n); km_queue.resize(n); km_visX.resize(n); }
        if (km_ly.size() < (size_t)m) { km_ly.resize(m); km_my.resize(m); km_slack.resize(m); km_slackmy.resize(m); km_visY.resize(m); }

        // Reset buffers
        fill(km_lx.begin(), km_lx.begin() + n, 0);
        fill(km_ly.begin(), km_ly.begin() + m, 0);
        fill(km_mx.begin(), km_mx.begin() + n, -1);
        fill(km_my.begin(), km_my.begin() + m, -1);

        // 1. Initialization & Row Reduction
        long long lb_row = 0;
        for (int i = 0; i < n; i++) {
            int min_val = INF_INT;
            for (int j = 0; j < m; j++) {
                if ((int)cost_matrix[i][j] < min_val) min_val = (int)cost_matrix[i][j];
            }
            km_lx[i] = min_val;
            if (min_val == INF_INT) return limit + 1;
            lb_row += min_val;
        }
        if (lb_row > limit) return limit + 1;

        // Column Reduction check omitted for brevity/speed trade-off, or can be added back.
        
        // 2. Greedy Initialization
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (km_my[j] == -1 && (int)cost_matrix[i][j] == km_lx[i]) {
                    km_mx[i] = j;
                    km_my[j] = i;
                    break;
                }
            }
        }

        // 3. Augmentation
        for (int u = 0; u < n; u++) {
            if (km_mx[u] != -1) continue;

            fill(km_visX.begin(), km_visX.begin() + n, 0);
            fill(km_visY.begin(), km_visY.begin() + m, 0);
            
            for (int j = 0; j < m; j++) {
                km_slack[j] = (int)cost_matrix[u][j] - km_lx[u] - km_ly[j];
                km_slackmy[j] = u;
            }

            km_visX[u] = 1;
            int q_size = 0;
            km_queue[q_size++] = u;
            int target_y = -1, last_x = -1;

            while (true) {
                int q_ptr = 0;
                while (q_ptr < q_size) {
                    int v = km_queue[q_ptr++];
                    for (int j = 0; j < m; j++) {
                        if (!km_visY[j] && (int)cost_matrix[v][j] == km_lx[v] + km_ly[j]) {
                            if (km_my[j] == -1) { last_x = v; target_y = j; goto found_path; }
                            km_visY[j] = 1;
                            int next_x = km_my[j];
                            km_visX[next_x] = 1;
                            km_prev[next_x] = v;
                            km_queue[q_size++] = next_x;
                            for (int k = 0; k < m; k++) {
                                int s = (int)cost_matrix[next_x][k] - km_lx[next_x] - km_ly[k];
                                if (!km_visY[k] && s < km_slack[k]) {
                                    km_slack[k] = s; km_slackmy[k] = next_x;
                                }
                            }
                        }
                    }
                }

                int delta = INF_INT;
                for (int j = 0; j < m; j++) if (!km_visY[j] && km_slack[j] < delta) delta = km_slack[j];
                if (delta == INF_INT) return limit + 1;

                for (int i = 0; i < n; i++) if (km_visX[i]) km_lx[i] += delta;
                for (int j = 0; j < m; j++) {
                    if (km_visY[j]) km_ly[j] -= delta;
                    else km_slack[j] -= delta;
                }

                q_size = 0;
                for (int j = 0; j < m; j++) {
                    if (!km_visY[j] && km_slack[j] == 0) {
                        if (km_my[j] == -1) { last_x = km_slackmy[j]; target_y = j; goto found_path; }
                        km_visY[j] = 1;
                        if (!km_visX[km_my[j]]) {
                            int next_x = km_my[j];
                            km_visX[next_x] = 1;
                            km_prev[next_x] = km_slackmy[j];
                            km_queue[q_size++] = next_x;
                            for (int k = 0; k < m; k++) {
                                int s = (int)cost_matrix[next_x][k] - km_lx[next_x] - km_ly[k];
                                if (!km_visY[k] && s < km_slack[k]) { km_slack[k] = s; km_slackmy[k] = next_x; }
                            }
                        }
                    }
                }
            }

            found_path:
            while (true) {
                int ty = km_mx[last_x];
                km_mx[last_x] = target_y;
                km_my[target_y] = last_x;
                if (last_x == u) break;
                target_y = ty;
                last_x = km_prev[last_x];
            }
        }

        unsigned long long total_cost = 0;
        for (int i = 0; i < n; i++) if (km_mx[i] != -1) total_cost += cost_matrix[i][km_mx[i]];
        if (total_cost > limit) return limit + 1;
        return total_cost;
    }

    // Lower Bound
    ui calLowerBound(const vector<pair<ui, ui>> &Mcand, ui last_u, ui last_v, 
                     vector<tuple<ui, ui, ui>> &restore_log)
    {
        Timer t_part;
        t_part.restart();

        // 1. Identify Nq, Ng, Rq
        lb_token++;
        if (lb_token == 0) {
            fill(lb_in_Nq.begin(), lb_in_Nq.end(), 0);
            fill(lb_in_Ng.begin(), lb_in_Ng.end(), 0);
            lb_token = 1;
        }

        vector<ui> Nq, Ng, Rq; 
        // Clearing adj_Mcand
        for(ui i = 0; i < qn; ++i) lb_adj_Mcand[i].clear();

        for (const auto& p : Mcand) {
            ui u = p.first; ui v = p.second;
            lb_adj_Mcand[u].push_back(v);
            if (lb_in_Nq[u] != lb_token) { lb_in_Nq[u] = lb_token; Nq.push_back(u); }
            if (lb_in_Ng[v] != lb_token) { lb_in_Ng[v] = lb_token; Ng.push_back(v); }
        }
        for (ui u = 0; u < qn; ++u) if (mapped_q[u] == -1 && lb_in_Nq[u] != lb_token) Rq.push_back(u);
        stats.lb_identify_sets_time += t_part.elapsed();

        // 2. Calculate LB_remain (Greedy relaxation for Rq) [MOD](doubled)
        t_part.restart();
        ui LB_remain_doubled = 0;
        vector<ui> dirty_nodes;

        if(last_u == (ui)-1) {
            dirty_nodes = Rq;
        } else {
            mark_token++;
            ui g_deg; const ui* g_neighbors = data_graph->getVertexNeighbors(last_v, g_deg);
            for(ui i = 0; i < g_deg; i++) g_mark[g_neighbors[i]] = mark_token;
            g_mark[last_v] = mark_token;

            for (ui u : Rq) {
                bool is_dirty = false;
                if (query_graph->hasEdge(u, last_u)) is_dirty = true;
                else {
                    ui old_best = best_cands[u];
                    if (old_best != (ui)-1 && g_mark[old_best] == mark_token) is_dirty = true;
                }
                if (is_dirty) dirty_nodes.push_back(u);
            }
        }

        for (ui u : dirty_nodes) {
            restore_log.emplace_back(u, min_costs[u], best_cands[u]);

            ui min_u_cost = (ui)-1;
            ui best_v = (ui)-1;
            ui deg_u = Lq_degrees[u];

            for (ui v : candidates[u]) {
                if (mapped_g[v] != -1) continue;
                ui deg_v = Lg_degrees[v];
                if (min_u_cost != (ui)-1) if (deg_u > deg_v && (deg_u - deg_v) >= min_u_cost) continue;

                // ui cost = computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
                // [MOD] USE OPTIMIZED DELTA
                ui cost = getDelta(u, v);

                if (cost < min_u_cost) {
                    min_u_cost = cost;
                    best_v = v;
                    if (min_u_cost == 0) break;
                }
            }
            min_costs[u] = min_u_cost;
            best_cands[u] = best_v;
        }
        
        for (ui u : Rq) {
            if (min_costs[u] == (ui)-1) return threshold + 1;
            LB_remain_doubled += min_costs[u];
        }
        stats.lb_remain_time += t_part.elapsed();

        if (LB_remain_doubled > 2 * threshold) return threshold + 1;

        // 3. Calculate LB_front (MWPM)
        Timer t_mwpm;
        t_part.restart();
        ui LB_front_doubled = 0;
        
        vector<pair<LabelID, ui>> Nq_labeled;
        Nq_labeled.reserve(Nq.size());
        for (ui u : Nq) Nq_labeled.emplace_back(query_graph->getVertexLabel(u), u);
        sort(Nq_labeled.begin(), Nq_labeled.end());

        ui remaining_budget_doubled = 2 * threshold - LB_remain_doubled;
        vector<ui> cols_to_reset;

        size_t ptr = 0;
        while(ptr < Nq_labeled.size()) {
            LabelID label = Nq_labeled[ptr].first;
            size_t qtr = ptr;
            while(qtr < Nq_labeled.size() && Nq_labeled[qtr].first == label) qtr++;

            size_t u_count = qtr - ptr;
            vector<ui> v_list;
            v_list.reserve(Ng.size());
            cols_to_reset.clear();

            for (ui v : Ng) {
                if (data_graph->getVertexLabel(v) == label) {
                    lb_v_to_col[v] = v_list.size();
                    v_list.push_back(v);
                    cols_to_reset.push_back(v);
                }
            }

            vector<ui> w_ext_doubled(u_count);
            for (size_t i = 0; i < u_count; ++i) {
                ui u = Nq_labeled[ptr + i].second;
                ui cross_fixed = 0;
                ui deg; const ui* neighbors = query_graph->getVertexNeighbors(u, deg);
                for(ui k = 0; k < deg; ++k) if(mapped_q[neighbors[k]] != -1) cross_fixed++;

                ui min_inner_doubled = (ui)-1; 
                bool found_rg = false;
                ui deg_u = Lq_degrees[u];

                for (ui v : candidates[u]) {
                    if (mapped_g[v] != -1) continue;
                    if (lb_in_Ng[v] == lb_token) continue;
                    
                    ui deg_v = Lg_degrees[v];
                    if (found_rg && deg_u > deg_v && (deg_u - deg_v) >= min_inner_doubled) continue;

                    // ui d = computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
                    // [MOD] USE OPTIMIZED DELTA
                    ui d = getDelta(u, v);

                    if (!found_rg || d < min_inner_doubled) { min_inner_doubled = d; found_rg = true; }
                    if (min_inner_doubled == 0) break;
                }
                if(!found_rg) min_inner_doubled = 0;
                w_ext_doubled[i] = (2 * cross_fixed) + min_inner_doubled;
            }

            size_t row_n = u_count;
            size_t col_n = v_list.size() + row_n;
            const ui INF_COST = 1e9;
            vector<vector<ui>> cost_matrix(row_n, vector<ui>(col_n, INF_COST));

            for (size_t i = 0; i < row_n; ++i) {
                ui u = Nq_labeled[ptr + i].second;
                // Virtual
                for (size_t k = 0; k < row_n; ++k) cost_matrix[i][v_list.size() + k] = w_ext_doubled[i];
                // Real
                for (ui v : lb_adj_Mcand[u]) {
                    int col_idx = lb_v_to_col[v];
                    if (col_idx != -1) {
                        ui d_cross = 0;
                        ui u_deg; const ui* u_neighbors = query_graph->getVertexNeighbors(u, u_deg);
                        for (ui k = 0; k < u_deg; ++k) {
                            ui un = u_neighbors[k];
                            if (mapped_q[un] != -1 && !data_graph->hasEdge(v, (ui)mapped_q[un])) d_cross++;
                        }

                        // cost_matrix[i][col_idx] = (2 * d_cross) + computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
                        // [MOD] USE OPTIMIZED DELTA
                        cost_matrix[i][col_idx] = (2 * d_cross) + getDelta(u, v);
                    }
                }
            }

            t_mwpm.restart();
            ui local_limit = remaining_budget_doubled - LB_front_doubled;
            ui local_cost = solveMWPM(cost_matrix, local_limit);
            stats.lb_mwpm_time += t_mwpm.elapsed();

            if (local_cost > local_limit) {
                LB_front_doubled = remaining_budget_doubled + 1;
                for(ui v : cols_to_reset) lb_v_to_col[v] = -1;
                break;
            }
            LB_front_doubled += local_cost;
            for(ui v : cols_to_reset) lb_v_to_col[v] = -1;
            ptr = qtr;
        }
        stats.lb_front_time += t_part.elapsed();

        return (LB_front_doubled + LB_remain_doubled + 1) / 2;
    }

    void dfs(ui missing, ui last_u, ui last_v, const vector<pair<ui, ui>> &fa_Mcand) 
    {
        stats.recursion_calls++;

        if (part_M.size() == qn) {
            results_ptr->push_back(part_M);
            return;
        }

        Timer t_cand;
        vector<pair<ui,ui> > Mcand;
        getIncrementalExpandableMappings(last_u, last_v, fa_Mcand, Mcand);
        stats.get_cand_time += t_cand.elapsed();
        if(Mcand.empty()) return;

        Timer t_lb;
        ui lb = 0;
        vector<tuple<ui, ui, ui>> restore_log;

#ifdef LOWER_BOUND
        lb = calLowerBound(Mcand, last_u, last_v, restore_log);
#endif
        stats.lb_time += t_lb.elapsed();

        Timer t_branch;
        if(missing + lb <= threshold) {
            vector<pair<ui, ui>> local_X_additions;
            local_X_additions.reserve(Mcand.size());
            for (auto p : Mcand) {
                t_branch.restart();
                ui u = p.first;
                ui v = p.second;
                ui delta = calIncrementalMissing(u, v);

                if(missing + delta <= threshold) {
                    // [MOD] Update Counts and Timestamps
                    updateNeighborCounts(query_graph, u, mapped_q, Lq_counts, Lq_degrees, false, q_last_update);
                    updateNeighborCounts(data_graph, v, mapped_g, Lg_counts, Lg_degrees, false, g_last_update);
                    mapped_q[u] = v;
                    mapped_g[v] = u;
                    part_M.emplace_back(u, v);
                    stats.branch_time += t_branch.elapsed();
                    
                    dfs(missing + delta, u, v, Mcand);

                    t_branch.restart();
                    mapped_q[u] = -1;
                    mapped_g[v] = -1;
                    part_M.pop_back();
                    // [MOD]
                    updateNeighborCounts(query_graph, u, mapped_q, Lq_counts, Lq_degrees, true, q_last_update);
                    updateNeighborCounts(data_graph, v, mapped_g, Lg_counts, Lg_degrees, true, g_last_update);
                }
                X.insert({u, v});
                local_X_additions.push_back({u, v});
                stats.branch_time += t_branch.elapsed();
            }
            t_branch.restart();
            for (const auto& p : local_X_additions) X.erase(p);
            stats.branch_time += t_branch.elapsed();
        } else {
            stats.prun_calls++;
        }

        // Restore LB state
        t_branch.restart();
        for (const auto &rec : restore_log) {
            min_costs[get<0>(rec)] = get<1>(rec);
            best_cands[get<0>(rec)] = get<2>(rec);
        }
        stats.branch_time += t_branch.elapsed();
    }
};

// ============================================================================
// TreeSpanSolver Implementation (EnumerateOnDemand Strategy)
// ============================================================================

class TreeSpanSolver {
public:
    TreeSpanSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

    bool init(const Graph *q, const Graph *d, ui match_threshold) {
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

    void match(vector<vector<pair<ui, ui>>> &results) {
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

    void printStats() const {
        printf("\n--- TreeSpan Matching Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, (double)stats.init_time/stats.total_time*100);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, (double)stats.search_time/stats.total_time*100);
        
        // Percentages below are relative to Search Time
        printf("- Verify Time:       %.4lf ms (%.2f%% of Search)\n", stats.verify_time / 1000.0, (stats.search_time > 0 ? (double)stats.verify_time/stats.search_time*100 : 0));
        printf("- Reorder Time:      %.4lf ms (%.2f%% of Search)\n", stats.reorder_time / 1000.0, (stats.search_time > 0 ? (double)stats.reorder_time/stats.search_time*100 : 0));
        
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Reorder Ops:         %lld\n", stats.reorder_calls);
        printf("Results Found:       %zu\n", results_ptr ? results_ptr->size() : 0);
        printf("---------------------------------------\n");
    }

private:
    // --- Data Structures ---
    
    struct QEdge {
        ui u, v;
        bool operator==(const QEdge& other) const {
            return (u == other.u && v == other.v) || (u == other.v && v == other.u);
        }
        // 为了放入 set 或比较，标准化为 u < v
        QEdge canonical() const {
            return (u < v) ? *this : QEdge{v, u};
        }
        bool operator<(const QEdge& other) const {
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
    void SimSearch(ui h, QISequence seq, ui gamma) {
        // Base Case: 成功匹配所有点
        if (h == qn) {
            vector<pair<ui, ui>> res;
            res.reserve(qn);
            for(ui u=0; u<qn; ++u) res.push_back({u, (ui)mapped_q[u]});
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
        ui deg_g; const ui* neighbors_g = data_graph->getVertexNeighbors(v_parent, deg_g);
        
        // 优化: 先将 u_curr 的候选集放入 bloom filter 或 hash set 加速查找? 
        // 这里直接使用 binary_search (前提: candidates 已排序)
        
        for (ui k = 0; k < deg_g; ++k) {
            ui v_curr = neighbors_g[k];
            
            if (mapped_g[v_curr] != -1) continue; // 必须未被匹配
            if (!binary_search(candidates[u_curr].begin(), candidates[u_curr].end(), v_curr)) continue; // 必须是候选点

            // 验证并计算新的 gamma
            ui new_gamma = gamma;
            bool possible = true;

            // 检查反向边 (Backward Edges)
            for (const auto& bedge : seq.bEdges[h]) {
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
                    } else {
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
                    SimSearch(h, next_seq, gamma + 1);
                }
            }
        }
    }

    // --- Helper Functions for Tree Generation ---

    // 寻找替代边：连接 S[0...h-1] 集合 和 (V - S[0...h-1]) 集合的边
    // 且该边不能是原本的树边，也不能在 R 中
    bool findReplacement(const QISequence& seq, ui h, QEdge& out_edge) {
        // 构建已访问集合 (Prefix)
        vector<bool> visited(qn, false);
        for(ui i=0; i<h; ++i) visited[seq.S[i]] = true;

        // 当前被移除的树边
        QEdge removed = seq.sEdge[h].canonical();

        // 寻找最优替代边 (遵循 Prim 序: 最小权/ID)
        bool found = false;
        QEdge best_e;

        for (const auto& e : all_q_edges) {
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
                } else {
                    if (cand < best_e) best_e = cand;
                }
            }
        }

        if (found) out_edge = best_e;
        return found;
    }

    // 基于前缀和新边，重新生成后续序列 (reOrdering)
    bool reorderSequence(const QISequence& old_seq, ui h, QEdge new_edge, QISequence& new_seq) {
        // 1. 复制前缀
        new_seq.S.resize(h);
        new_seq.sEdge.resize(h);
        new_seq.bEdges.resize(h); // bEdges 需要重新计算吗？前缀内部的 bEdges 不变
        
        vector<bool> visited(qn, false);
        for(ui i=0; i<h; ++i) {
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
        ui deg; const ui* nbrs = query_graph->getVertexNeighbors(u_next, deg);
        for(ui k=0; k<deg; ++k) {
            ui neighbor = nbrs[k];
            if (visited[neighbor]) {
                QEdge e = {u_next, neighbor};
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
            for (const auto& e : all_q_edges) {
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
            const ui* v_nbrs = query_graph->getVertexNeighbors(v_new, deg);
            for(ui k=0; k<deg; ++k) {
                ui neighbor = v_nbrs[k];
                if (visited[neighbor]) {
                    QEdge e = {v_new, neighbor};
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
    QISequence generateInitialSequence(ui root) {
        QISequence seq;
        seq.S.push_back(root);
        seq.sEdge.push_back({(ui)-1, (ui)-1}); // Root has no parent edge
        seq.bEdges.push_back({});
        
        // 调用 reorderSequence 的一部分逻辑来填充 (实际上就是 Prim)
        // 造一个假的 "0" 长度前缀，利用 reorderSequence 填充整个列表
        // 但这里手写 Prim 更清晰
        
        vector<bool> visited(qn, false);
        visited[root] = true;

        for (ui i = 1; i < qn; ++i) {
            QEdge best_edge;
            bool found = false;

            for (const auto& e : all_q_edges) {
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
                ui deg; const ui* nbrs = query_graph->getVertexNeighbors(next_u, deg);
                for(ui k=0; k<deg; ++k) {
                    ui nbr = nbrs[k];
                    if (visited[nbr]) {
                        QEdge be = {next_u, nbr};
                        if (!(be.canonical() == best_edge)) bes.push_back(be);
                    }
                }
                seq.bEdges.push_back(bes);
            }
        }
        return seq;
    }

    // --- Utilities ---

    ui selectRoot() {
        ui best_u = 0;
        size_t min_cand = candidates[0].size();
        for(ui u=1; u<qn; ++u) {
            if (candidates[u].size() < min_cand) {
                min_cand = candidates[u].size();
                best_u = u;
            }
        }
        return best_u;
    }

    vector<QEdge> getAllQueryEdges() const {
        vector<QEdge> edges;
        edges.reserve(query_graph->getEdgesCount());
        for (ui u = 0; u < qn; ++u) {
            ui deg; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (u < v) edges.push_back({u, v});
            }
        }
        // 排序确保确定性
        sort(edges.begin(), edges.end());
        return edges;
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees) {
        ui n = g->getVerticesCount();
        ui num_labels = g->getLabelsCount();
        counts.assign(n, vector<ui>(num_labels, 0));
        degrees.assign(n, 0);
        for (ui u = 0; u < n; ++u) {
            ui deg; const ui* neighbors = g->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                counts[u][g->getVertexLabel(neighbors[i])]++;
                degrees[u]++;
            }
        }
    }

    ui computeDelta(ui u, ui v) {
        ui diff = 0;
        size_t sz = Lq_counts[u].size();
        for (size_t i = 0; i < sz; ++i) {
            if (Lq_counts[u][i] > Lg_counts[v][i]) diff += (Lq_counts[u][i] - Lg_counts[v][i]);
        }
        return diff;
    }

    bool calVerticesFilter() {
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
    
    long long total_time = t_total.elapsed();
    printf("\n--- Total Time ---\n");
    printf("Total Time:          %.4lf ms\n", total_time / 1000.0);
    // solver.printStats();
}