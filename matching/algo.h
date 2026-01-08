#include "utility/utility.h"
// #include "graph/graph.h"

using namespace std;

class MatchingSolver{
public:
    MatchingSolver() : query_graph(nullptr), data_graph(nullptr), candidates_ptr(nullptr), results_ptr(nullptr) {}
    void match(const Graph *q, const Graph *d, const vector<vector<ui>> &cands,
               vector<vector<pair<ui, ui>>> &results, ui match_threshold) 
    {
        // 1. Initialization
        query_graph = q;
        data_graph = d;
        candidates_ptr = &cands;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();
        results_ptr = &results;
        resetState();

        Timer t_total;
        t_total.restart();

        // 2. Initial Setup
        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

        vector<pair<ui, ui>> initial_Mcand; // Empty for root

        // 3. Start Recursion
        Timer t_dfs;
        dfs(0, (ui)-1, (ui)-1, initial_Mcand);
        stats.dfs_time += t_dfs.elapsed();

        stats.total_time = t_total.elapsed();

        printStats();
    }
private:
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
    };

    // --- Input Data ---
    const Graph *query_graph;
    const Graph *data_graph;
    const vector<vector<ui>> *candidates_ptr; // Pointer to external candidates
    vector<vector<pair<ui, ui>>> *results_ptr; // Results
    ui threshold;
    ui qn, gn;

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

    // --- Stats ---
    TimeStats stats;

    // State Helpers
    void resetState() {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        part_M.clear();
        part_M.reserve(qn);
        X.clear();
        if (results_ptr != nullptr) {
            results_ptr->clear();
        }

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
        size_t sz = Lq.size();
        for (size_t i = 0; i < sz; ++i) if (Lq[i] > Lg[i]) diff += (Lq[i] - Lg[i]);
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
        const auto& candidates = *candidates_ptr;

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
        const auto& candidates = *candidates_ptr;

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
};

class TreeSpanSolver{
    
};

// ============================================================
// Top-level function: Approximate_Matching
// ============================================================
void Approximate_Matching(const Graph *query_graph, const Graph *data_graph,
                          vector<vector<ui> > &candidates, 
                          vector<vector<pair<ui, ui> > > &M_ANS,
                          ui threshold)
{
    MatchingSolver matcher;
    matcher.match(query_graph, data_graph, candidates, M_ANS, threshold);
}