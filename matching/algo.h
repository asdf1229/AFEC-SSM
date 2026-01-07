#include "utility/utility.h"
// #include "graph/graph.h"

using namespace std;

struct TimeStats {
    long long total_time = 0;
    long long get_cand_time = 0;   // ExpandableMappings 耗时
    long long lb_time = 0;         // calLowerBound 耗时
    long long missing_time = 0;    // calPartialMissing 耗时
    long long recursion_calls = 0; // 递归调用次数计数
    long long prun_calls = 0;      // 剪枝次数

    // calLowerBound 内部细分时间
    long long lb_identify_sets_time = 0; // 1. Identify Nq, Rq, Ng, Rg
    long long lb_remain_time = 0; // 2. Calculate LB_remain
    long long lb_front_time = 0; // 3. Calculate LB_front
};

LabelID label2int(const string str, map<string, LabelID> &M)
{
	if(M.find(str) == M.end()) M[str] = M.size();
	return M[str];
}

void load_graph(const string &input_graph, Graph *graph,
                map<string, LabelID> &vM, map<string, LabelID> &eM)
{
    ifstream fin(input_graph);
    string line;

    while (getline(fin, line)) if (!line.empty() && line[0] == 't') break;

    if(!fin) {
        printf("!!! Cannot open graph file %s !!!\n", input_graph.c_str());
        assert(false);
        exit(1);
    }

    istringstream head(line);
    char tchar; 
    string sharp, id;
    if (!(head >> tchar >> sharp >> id) || tchar != 't') {
        fprintf(stderr, "!!! Invalid graph header line: %s !!!\n", line.c_str());
        assert(false);
        exit(1);
    }

    vector<pair<ui, LabelID> > vertices;
    vector<pair<pair<ui, ui>, LabelID> > undirected_edges;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        char type = line[0];
        if (type == 't') break;

        istringstream iss(line);

        if(type == 'v') {
            char c;
            ui vid;
            string vlab;
            if (!(iss >> c >> vid >> vlab)) {
                fprintf(stderr, "!!! Invalid vertex line: %s !!!\n", line.c_str());
                assert(false);
                exit(1);
            }
            vertices.emplace_back(vid, label2int(vlab, vM));
        }
        else if (type == 'e') {
            char c;
            ui u, v;
            string elab;

            if (!(iss >> c >> u >> v)) {
                fprintf(stderr, "!!! Invalid edge line: %s !!!\n", line.c_str());
                assert(false);
                exit(1);
            }

            if (!(iss >> elab)) elab = "__NO_EDGE_LABEL__";

            if (u == v) continue;

            LabelID L = label2int(elab, eM);
            undirected_edges.emplace_back(make_pair(min(u, v), max(u, v)), L);
        }
        else {
            fprintf(stderr, "!!! Unknown line type: %s !!!\n", line.c_str());
            assert(false);
            exit(1);
        }
    }

    sort(vertices.begin(), vertices.end());
    sort(undirected_edges.begin(), undirected_edges.end());
    vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());
    undirected_edges.erase(unique(undirected_edges.begin(), undirected_edges.end()), undirected_edges.end());

    vector<pair<pair<ui, ui>, LabelID> > edges;
    edges.reserve(undirected_edges.size() * 2);

    for(auto &e: undirected_edges) {
        edges.emplace_back(make_pair(e.first.first, e.first.second), e.second);
        edges.emplace_back(make_pair(e.first.second, e.first.first), e.second);
    }

    sort(edges.begin(), edges.end());

    graph->build_graph(id, vertices, edges);
}

bool calVerticesFilter(const Graph *query_graph, const Graph *data_graph, vector<vector<ui> > &candidates)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();
    candidates.clear();
    candidates.resize(qn);

    for (ui i = 0; i < qn; ++i) {
        LabelID label_i = query_graph->getVertexLabel(i);
        for (ui j = 0; j < gn; ++j) {
            LabelID label_j = data_graph->getVertexLabel(j);
            if (label_i == label_j) candidates[i].push_back(j);
        }
        if (candidates[i].empty())  return false;
    }
    return true;
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
                          vector<vector<ui>> &counts, vector<ui> &degrees, bool is_add)
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
        }
    }
}

ui calPartialMissing(const Graph *query_graph, const Graph *data_graph,
                     const vector<pair<ui, ui>> &part_M)
{
    ui missing = 0;
    ui mn = part_M.size();
    for (ui i = 0; i < mn; ++i) {
        ui u1 = part_M[i].first;
        ui v1 = part_M[i].second;
        for (ui j = i + 1; j < mn; ++j) {
            ui u2 = part_M[j].first;
            ui v2 = part_M[j].second;
            if (query_graph->hasEdge(u1, u2)) {
                if (!data_graph->hasEdge(v1, v2))
                    missing++;
            }
        }
    }
    return missing;
}

ui calIncrementalMissing(const Graph *query_graph, const Graph *data_graph,
                         const vector<pair<ui, ui>> &part_M, ui new_u, ui new_v)
{
    ui delta_missing = 0;
    for (const auto &p : part_M) {
        ui u_prev = p.first;
        ui v_prev = p.second;
        
        if (query_graph->hasEdge(new_u, u_prev) && !data_graph->hasEdge(new_v, v_prev)) delta_missing++;
    }
    return delta_missing;
}

void getIncrementalExpandableMappings(const Graph *query_graph, const Graph *data_graph,
                                      const vector<int> &mapped_q, const vector<int> &mapped_g,
                                      const vector<vector<ui> > &candidates,
                                      const unordered_set<pair<ui, ui>, PairHash> &X, 
                                      const vector<pair<ui, ui>> &fa_Mcand, 
                                      const vector<pair<ui, ui>> &part_M,
                                      ui new_u, ui new_v,                         
                                      vector<pair<ui,ui> > &Mcand)
{
    Mcand.clear();

    if (part_M.empty()) {
        assert(new_u == (ui)-1 && new_v == (ui)-1);
        // select 0 as starting vertex
        ui u = 0;
        assert(mapped_q[u] == -1);
        for(ui v : candidates[u]) {
            assert(mapped_g[v] == -1);
            assert(!X.count({u, v}));
            Mcand.emplace_back(u, v);
        }
        return;
    }

    Mcand.reserve(fa_Mcand.size() + 64);

    for (const auto &p : fa_Mcand) {
        if (p.first == new_u || p.second == new_v) continue;
        if (X.count(p)) continue; 
        Mcand.push_back(p);
    }

    ui q_count;
    const ui* q_neighbors = query_graph->getVertexNeighbors(new_u, q_count);
    ui g_count;
    const ui* g_neighbors = data_graph->getVertexNeighbors(new_v, g_count);

    // 记录继承部分的结束位置，用于后续归并
    size_t inherit_end = Mcand.size();
    ui gn = data_graph->getVerticesCount();
    static vector<bool> g_mark;
    if (g_mark.size() < gn) g_mark.resize(gn, false);


    for(ui i = 0; i < g_count; ++i) g_mark[g_neighbors[i]] = true;

    for(ui i = 0; i < q_count; ++i) {
        ui u2 = q_neighbors[i];
        if(mapped_q[u2] != -1) continue;

        for(ui v2 : candidates[u2]) {
            if(mapped_g[v2] != -1) continue;
            if(g_mark[v2]) {
                if(X.count({u2, v2})) continue;
                Mcand.emplace_back(u2, v2);
            }
        }


        for(ui j = 0; j < g_count; ++j) {
            ui v2 = g_neighbors[j];
            if(mapped_g[v2] != -1) continue; 

            bool is_cand = false;
            for(ui c : candidates[u2]) if(c == v2) { is_cand = true; break; }
            if(!is_cand) continue;

            if(X.count({u2, v2})) continue;

            Mcand.emplace_back(u2, v2);
        }
    }

    for(ui i = 0; i < g_count; ++i) g_mark[g_neighbors[i]] = false;

    // 3. 优化去重
    // Mcand 前半部分 (0 ~ inherit_end-1) 已经有序
    // Mcand 后半部分 (inherit_end ~ end) 是新加入的，局部可能无序
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
ui solveMWPM(const vector<vector<ui>>& cost_matrix, ui limit) { //[MOD] Add limit parameter
    if (cost_matrix.empty()) return 0;

    int n = (int)cost_matrix.size();    // 左侧节点数 (Nq)
    int m = (int)cost_matrix[0].size(); // 右侧节点数 (Ng + Virtual)
    const int INF_INT = 1e9;

    // lx, ly 为可行顶标；mx, my 为匹配状态
    vector<int> lx(n, 0), ly(m, 0);
    vector<int> mx(n, -1), my(m, -1);
    vector<int> slack(m);       // KM 优化：降低修改顶标的时间复杂度
    vector<int> slackmy(m);    // 记录 slack 来源的左侧节点
    vector<int> prev(n);       // 记录交错树中的前驱
    vector<int> queue(n);      // BFS 队列
    vector<char> visX(n), visY(m);

    // 1. 初始化顶标 (Initialization) & 行规约检查 (Row Reduction)
    // [MOD] Strategy D & A: 计算行最小值之和
    long long lb_row = 0; 
    for (int i = 0; i < n; i++) {
        int min_val = INF_INT;
        for (int j = 0; j < m; j++) {
            if ((int)cost_matrix[i][j] < min_val) min_val = (int)cost_matrix[i][j];
        }
        lx[i] = min_val; // 左侧顶标初始化为行最小值
        if (min_val == INF_INT) return limit + 1; // 某行全为 INF，无解
        lb_row += min_val; 
    }

    // [MOD] 如果仅行最小值就已经超过限制，直接返回
    if (lb_row > limit) return limit + 1;

    // [MOD] Strategy A: Column Reduction (列规约) 下界检查
    // 计算 reduce 后的矩阵 C'[i][j] = C[i][j] - lx[i] 的列最小值
    // 真正的下界 >= lb_row + sum(前 n 小的列最小值)
    {
        vector<int> col_mins(m, INF_INT);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                int reduced_val = (int)cost_matrix[i][j] - lx[i];
                if (reduced_val < col_mins[j]) col_mins[j] = reduced_val;
            }
        }
        
        // 我们需要找到最小的 n 个列规约值。
        // 由于 n 通常较小，且 m 也较小，直接排序或 nth_element 均可。
        // 这里使用 sort，因为 m 通常很小 (degree size)
        sort(col_mins.begin(), col_mins.end());
        
        long long lb_col = 0;
        for (int k = 0; k < n; ++k) {
            if (col_mins[k] == INF_INT) return limit + 1;
            lb_col += col_mins[k];
        }

        if (lb_row + lb_col > limit) return limit + 1;
    }

    // 2. 贪心匹配初始化 (Greedy Match Initialization)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (my[j] == -1 && (int)cost_matrix[i][j] == lx[i]) {
                mx[i] = j;
                my[j] = i;
                break;
            }
        }
    }

    // 3. 增广 (Augmentation)
    for (int u = 0; u < n; u++) {
        if (mx[u] != -1) continue; // 已匹配则跳过

        fill(visX.begin(), visX.end(), 0);
        fill(visY.begin(), visY.end(), 0);
        
        // 初始化 slack 数组
        for (int j = 0; j < m; j++) {
            slack[j] = (int)cost_matrix[u][j] - lx[u] - ly[j];
            slackmy[j] = u;
        }

        visX[u] = 1;
        int q_size = 0;
        queue[q_size++] = u;

        int target_y = -1;
        int last_x = -1;

        while (true) {
            // 尝试在当前相等子图中寻找增广路
            int q_ptr = 0;
            while (q_ptr < q_size) {
                int v = queue[q_ptr++];
                for (int j = 0; j < m; j++) {
                    if (!visY[j] && (int)cost_matrix[v][j] == lx[v] + ly[j]) {
                        if (my[j] == -1) { // 找到未匹配的右侧节点
                            last_x = v; target_y = j; goto found_path;
                        }
                        visY[j] = 1;
                        int next_x = my[j];
                        visX[next_x] = 1;
                        prev[next_x] = v;
                        queue[q_size++] = next_x;
                        // 更新新加入 X 集合的节点带来的 slack
                        for (int k = 0; k < m; k++) {
                            int s = (int)cost_matrix[next_x][k] - lx[next_x] - ly[k];
                            if (!visY[k] && s < slack[k]) {
                                slack[k] = s; slackmy[k] = next_x;
                            }
                        }
                    }
                }
            }

            // 找不到增广路，修改顶标以扩大相等子图
            int delta = INF_INT;
            for (int j = 0; j < m; j++) if (!visY[j] && slack[j] < delta) delta = slack[j];
            
            // [MOD] 如果无法继续增广，或者 delta 还是无穷大，说明无解
            if (delta == INF_INT) return limit + 1;
            
            // [MOD] Strategy D (Lazy Check):
            // 注意：严格来说 sum(lx)+sum(ly) 是单调递减直到最优解的。
            // 但如果当前的 delta 过大，虽然不一定意味着立刻超标（因为可能有负调整），
            // 在整数权重的 KM 中，如果 delta 极大通常意味着不可行。
            // 更安全的检测是每轮增广结束后检测总 cost。
            
            for (int i = 0; i < n; i++) if (visX[i]) lx[i] += delta;
            for (int j = 0; j < m; j++) {
                if (visY[j]) ly[j] -= delta;
                else slack[j] -= delta;
            }

            // 顶标修改后，检查是否有新的边进入相等子图
            q_size = 0;
            for (int j = 0; j < m; j++) {
                if (!visY[j] && slack[j] == 0) {
                    if (my[j] == -1) {
                        last_x = slackmy[j]; target_y = j; goto found_path;
                    }
                    visY[j] = 1;
                    if (!visX[my[j]]) {
                        int next_x = my[j];
                        visX[next_x] = 1;
                        prev[next_x] = slackmy[j];
                        queue[q_size++] = next_x;
                        for (int k = 0; k < m; k++) {
                            int s = (int)cost_matrix[next_x][k] - lx[next_x] - ly[k];
                            if (!visY[k] && s < slack[k]) {
                                slack[k] = s; slackmy[k] = next_x;
                            }
                        }
                    }
                }
            }
        }

        found_path:
        // 沿增广路回溯更新匹配
        while (true) {
            int ty = mx[last_x];
            mx[last_x] = target_y;
            my[target_y] = last_x;
            if (last_x == u) break;
            target_y = ty;
            last_x = prev[last_x];
        }
    }

    // 统计总代价
    unsigned long long total_cost = 0;
    for (int i = 0; i < n; i++) {
        if (mx[i] != -1) total_cost += cost_matrix[i][mx[i]];
    }
    
    // [MOD] Strategy D: 最后再次检查总代价
    if (total_cost > limit) return limit + 1;

    return total_cost;
}

ui computeDeltaInnerDoubledVector(const vector<ui> &Lq, const vector<ui> &Lg)
{
    ui diff = 0;
    size_t sz = Lq.size();
    for (size_t i = 0; i < sz; ++i) {
        if (Lq[i] > Lg[i]) diff += (Lq[i] - Lg[i]);
    }
    return diff;
}

void getUnmappedNeighborLabelsVector(const Graph *g, ui u, const vector<int> &mapped, vector<ui> &label_counts) {
    ui count;
    const ui* neighbors = g->getVertexNeighbors(u, count);
    for (ui i = 0; i < count; ++i) {
        ui neighbor = neighbors[i];
        if (mapped[neighbor] == -1) {
            label_counts[g->getVertexLabel(neighbor)]++;
        }
    }
}

ui calLowerBound(const Graph *query_graph, const Graph *data_graph,
                 const vector<int> &mapped_q, const vector<int> &mapped_g,
                 const vector<vector<ui> > &candidates, const vector<pair<ui, ui>> &part_M, 
                 vector<pair<ui,ui> > &Mcand, TimeStats &stats, ui threshold, 
                 const vector<vector<ui>> &Lq_counts, const vector<vector<ui>> &Lg_counts,
                 const vector<ui> &Lq_degrees, const vector<ui> &Lg_degrees,
                 vector<ui> &min_costs, vector<ui> &best_cands, 
                 vector<tuple<ui, ui, ui>> &restore_log, // log: <u, old_cost, old_cand>
                 ui last_u, ui last_v, vector<int> &g_mark, int &mark_token)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();

    Timer t_part;
    t_part.restart();

    // 1. Identify Nq, Ng, Rq, Rg
    static vector<int> in_Nq;
    static vector<int> in_Ng;
    static int lb_token = 0;
    
    static vector<ui> Nq;
    static vector<ui> Ng;
    static vector<ui> Rq;
    // vector<ui> Rg;
    static vector<vector<ui>> adj_Mcand; // Adjacency list for bipartite matching

    if (in_Nq.size() < qn) {
        in_Nq.resize(qn, 0);
        adj_Mcand.resize(qn);
    }
    if (in_Ng.size() < gn) {
        in_Ng.resize(gn, 0);
    }

    lb_token++;
    if (lb_token == 0) { 
        fill(in_Nq.begin(), in_Nq.end(), 0);
        fill(in_Ng.begin(), in_Ng.end(), 0);
        lb_token = 1;
    }

    Nq.clear();
    Ng.clear();
    Rq.clear();
    for(ui i = 0; i < qn; ++i) adj_Mcand[i].clear();

    for (const auto& p : Mcand) {
        ui u = p.first; ui v = p.second;
        adj_Mcand[u].push_back(v);
        if (in_Nq[u] != lb_token) { 
            in_Nq[u] = lb_token; 
            Nq.push_back(u); 
        }
        if (in_Ng[v] != lb_token) { 
            in_Ng[v] = lb_token; 
            Ng.push_back(v); 
        }
    }
    for (ui u = 0; u < qn; ++u) if (mapped_q[u] == -1 && in_Nq[u] != lb_token) Rq.push_back(u);
    // for (ui v = 0; v < gn; ++v) if (mapped_g[v] == -1 && !in_Ng[v]) Rg.push_back(v);
    stats.lb_identify_sets_time += t_part.elapsed();

    // 2. Calculate LB_remain (Greedy relaxation for Rq) [MOD](doubled)
    t_part.restart();
    ui LB_remain_doubled = 0;
    vector<ui> dirty_nodes;
    if(last_u == (ui)-1) {
        dirty_nodes = Rq;
    }
    else {
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

    // update dirty nodes
    for (ui u : dirty_nodes) {
        restore_log.emplace_back(u, min_costs[u], best_cands[u]);

        ui min_u_cost = (ui)-1;
        ui best_v = (ui)-1;
        ui deg_u = Lq_degrees[u];

        for (ui v : candidates[u]) {
            if (mapped_g[v] != -1) continue;
            ui deg_v = Lg_degrees[v];
            if (min_u_cost != (ui)-1) if (deg_u > deg_v && (deg_u - deg_v) >= min_u_cost) continue;
            ui cost = computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
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

    // early exit
    if (LB_remain_doubled > 2 * threshold) return threshold + 1;

    // 3. Calculate LB_front (MWPM for Nq)
    t_part.restart();
    ui LB_front_doubled = 0;
    
    vector<pair<LabelID, ui>> Nq_labeled;
    Nq_labeled.reserve(Nq.size());
    for (ui u : Nq) Nq_labeled.emplace_back(query_graph->getVertexLabel(u), u);
    sort(Nq_labeled.begin(), Nq_labeled.end());
    
    // Optimization: Avoid O(|V_G|) initialization for v_to_col
    static vector<int> v_to_col;
    if (v_to_col.size() < gn) v_to_col.resize(gn, -1);
    // Use a list to track modified indices for fast reset
    vector<ui> cols_to_reset; 

    assert(LB_remain_doubled <= 2 * threshold);
    ui remaining_budget_doubled = 2 * threshold - LB_remain_doubled;

    size_t ptr = 0;
    while(ptr < Nq_labeled.size()) {
        LabelID label = Nq_labeled[ptr].first;
        size_t qtr = ptr;
        while(qtr < Nq_labeled.size() && Nq_labeled[qtr].first == label) qtr++;

        size_t u_count = qtr - ptr;
        
        // Build v_list from Ng (subset of Mcand targets)
        vector<ui> v_list;
        v_list.reserve(Ng.size()); // Heuristic reserve
        cols_to_reset.clear();

        // Note: iterating Ng is O(|Ng|), much smaller than O(|V_G|)
        for (ui v : Ng) {
            if (data_graph->getVertexLabel(v) == label) {
                v_to_col[v] = v_list.size();
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

            // Optimization 2: Implicit Rg. 
            // Instead of iterating Rg (large), iterate candidates[u] (small)
            // and check conditions: unmapped AND not in Ng.
            const auto& cands = candidates[u];
            for (ui v : cands) {
                // Filter: must be unmapped
                if (mapped_g[v] != -1) continue;
                // Filter: must NOT be in Ng (this makes it effectively in Rg)
                if (in_Ng[v] == lb_token) continue;
                
                // Now v is effectively in Rg
                ui deg_v = Lg_degrees[v];
                if (found_rg && deg_u > deg_v && (deg_u - deg_v) >= min_inner_doubled) continue;
                ui d = computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
                if (!found_rg || d < min_inner_doubled) { min_inner_doubled = d; found_rg = true; }
                if (min_inner_doubled == 0) break; 
            }

            if(!found_rg) min_inner_doubled = 0; 
            w_ext_doubled[i] = (2 * cross_fixed) + min_inner_doubled;
        }

        size_t row_n = u_count;
        size_t col_n = v_list.size() + row_n;
        const ui INF_COST = 1e9;
        // cost_matrix size is small (based on Mcand size), re-allocation is acceptable or can be optimized further
        vector<vector<ui>> cost_matrix(row_n, vector<ui>(col_n, INF_COST));

        for (size_t i = 0; i < row_n; ++i) {
            ui u = Nq_labeled[ptr + i].second;
            // Virtual
            for (size_t k = 0; k < row_n; ++k) cost_matrix[i][v_list.size() + k] = w_ext_doubled[i];
            // Real
            for (ui v : adj_Mcand[u]) {
                // v is in Ng, check if it has a column index
                int col_idx = v_to_col[v];
                if (col_idx != -1) {
                    ui d_cross = 0;
                    ui u_deg; const ui* u_neighbors = query_graph->getVertexNeighbors(u, u_deg);
                    for (ui k = 0; k < u_deg; ++k) {
                        ui un = u_neighbors[k];
                        if (mapped_q[un] != -1 && !data_graph->hasEdge(v, (ui)mapped_q[un])) d_cross++;
                    }
                    cost_matrix[i][col_idx] = (2 * d_cross) + computeDeltaInnerDoubledVector(Lq_counts[u], Lg_counts[v]);
                }
            }
        }

        assert(LB_front_doubled <= remaining_budget_doubled);

        ui local_limit = remaining_budget_doubled - LB_front_doubled;
        ui local_cost = solveMWPM(cost_matrix, local_limit);
        
        if (local_cost > local_limit) {
            LB_front_doubled = remaining_budget_doubled + 1; // 标记剪枝
            for(ui v : cols_to_reset) v_to_col[v] = -1;
            break;
        }
        
        LB_front_doubled += local_cost;
        
        // Fast reset for v_to_col
        for(ui v : cols_to_reset) v_to_col[v] = -1;
        
        ptr = qtr;
    }
    stats.lb_front_time += t_part.elapsed();

    return (LB_front_doubled + LB_remain_doubled + 1) / 2;
}

void DFS_Approximate(const Graph *query_graph, const Graph *data_graph,
                     const vector<vector<ui>> &candidates, unordered_set<pair<ui, ui>, PairHash> &X,
                     vector<int> &mapped_q, vector<int> &mapped_g,
                     vector<pair<ui, ui>> &part_M, vector<vector<pair<ui,ui>>> &M_ANS,
                     ui threshold, ui missing, TimeStats &stats, 
                     vector<vector<ui>> &Lq_counts, vector<vector<ui>> &Lg_counts,
                     vector<ui> &Lq_degrees, vector<ui> &Lg_degrees,
                     ui last_u, ui last_v, const vector<pair<ui, ui>> &fa_Mcand, 
                     vector<ui> &min_costs, vector<ui> &best_cands, 
                     vector<int> &g_mark, int &mark_token)
{
#ifndef NDEBUG
    for(ui i = 0; i < part_M.size(); ++i) {
        ui u = part_M[i].first;
        ui v = part_M[i].second;
        assert(u < query_graph->getVerticesCount());
        assert(v < data_graph->getVerticesCount());
        assert(mapped_q[u] == static_cast<int>(v));
        assert(mapped_g[v] == static_cast<int>(u));
    }
    assert(missing == calPartialMissing(query_graph, data_graph, part_M));
#endif
    stats.recursion_calls++;

    if (part_M.size() == query_graph->getVerticesCount()) {
        assert(missing <= threshold);
        M_ANS.push_back(part_M);
        return;
    }

    Timer t_cand;
    vector<pair<ui,ui> > Mcand;
    getIncrementalExpandableMappings(query_graph, data_graph, mapped_q, mapped_g, 
                                     candidates, X, fa_Mcand, part_M, 
                                     last_u, last_v, Mcand);
    stats.get_cand_time += t_cand.elapsed();
    if(Mcand.empty()) return;

    Timer t_lb;
    ui lb = 0;

    vector<tuple<ui, ui, ui>> restore_log;

#ifdef LOWER_BOUND
    lb = calLowerBound(query_graph, data_graph, mapped_q, mapped_g, candidates, part_M, Mcand, stats, threshold,
                       Lq_counts, Lg_counts, Lq_degrees, Lg_degrees,
                       min_costs, best_cands, restore_log, last_u, last_v, g_mark, mark_token);
#endif
    stats.lb_time += t_lb.elapsed();

    if(missing + lb <= threshold) {
        vector<pair<ui, ui>> local_X_additions;
        local_X_additions.reserve(Mcand.size());

        for (auto p : Mcand) {
            ui u = p.first;
            ui v = p.second;

            Timer t_miss;
            ui delta = calIncrementalMissing(query_graph, data_graph, part_M, u, v);
            stats.missing_time += t_miss.elapsed();

            if(missing + delta <= threshold) {
                assert(X.count({u, v}) == 0);
                assert(mapped_q[u] == -1);
                assert(mapped_g[v] == -1);

                updateNeighborCounts(query_graph, u, mapped_q, Lq_counts, Lq_degrees, false);
                updateNeighborCounts(data_graph, v, mapped_g, Lg_counts, Lg_degrees, false);
                mapped_q[u] = v;
                mapped_g[v] = u;
                part_M.emplace_back(u, v);

                DFS_Approximate(query_graph, data_graph, candidates, X, mapped_q, mapped_g, 
                                part_M, M_ANS, threshold, missing + delta, stats, 
                                Lq_counts, Lg_counts, Lq_degrees, Lg_degrees, u, v, Mcand,
                                min_costs, best_cands, g_mark, mark_token);

                mapped_q[u] = -1;
                mapped_g[v] = -1;
                part_M.pop_back();
                updateNeighborCounts(query_graph, u, mapped_q, Lq_counts, Lq_degrees, true);
                updateNeighborCounts(data_graph, v, mapped_g, Lg_counts, Lg_degrees, true);
            }
            X.insert({u, v});
            local_X_additions.push_back({u, v});
        }
        for (const auto& p : local_X_additions) X.erase(p);
    } 
    else {
        stats.prun_calls++;
    }

    for (const auto &rec : restore_log) {
        min_costs[get<0>(rec)] = get<1>(rec);
        best_cands[get<0>(rec)] = get<2>(rec);
    }
}

// ============================================================
// Top-level function: Approximate_Matching_v2
// ============================================================
void Approximate_Matching_v2(const Graph *query_graph, const Graph *data_graph,
                             vector<vector<ui> > &candidates, 
                             vector<vector<pair<ui, ui> > > &M_ANS,
                             ui threshold)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();

    M_ANS.clear();
    assert(qn && gn);

    vector<int> mapped_q(qn, -1);
    vector<int> mapped_g(gn, -1);
    unordered_set<pair<ui, ui>, PairHash> X;
    vector<pair<ui, ui>> part_M; // (Mq, Mg)

    TimeStats stats;
    
    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;
    initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
    initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

    vector<ui> min_costs(qn, (ui)-1);
    vector<ui> best_cands(qn, (ui)-1);
    vector<int> g_mark(gn, 0); // 用于 O(1) 检查数据节点是否受影响
    int mark_token = 0;

    vector<pair<ui, ui>> initial_Mcand;

    Timer t_total;
    DFS_Approximate(query_graph, data_graph, candidates, X, mapped_q, mapped_g,
                    part_M, M_ANS, threshold, 0, stats, 
                    Lq_counts, Lg_counts, Lq_degrees, Lg_degrees, (ui)-1, (ui)-1, initial_Mcand,
                    min_costs, best_cands, g_mark, mark_token);
    stats.total_time = t_total.elapsed();

    printf("\n--- Matching v2 Time Analysis ---\n");
    printf("Total Time:          %lld us\n", stats.total_time);
    printf("Mcand Generation:    %lld us (%.2f%%)\n", stats.get_cand_time, (double)stats.get_cand_time/stats.total_time*100);
    printf("LowerBound Total:    %lld us (%.2f%%)\n", stats.lb_time, (double)stats.lb_time/stats.total_time*100);
    printf("  - Identify Sets:   %lld us (%.2f%% of LB)\n", stats.lb_identify_sets_time, (double)stats.lb_identify_sets_time/stats.lb_time*100);
    printf("  - LB Remain:       %lld us (%.2f%% of LB)\n", stats.lb_remain_time, (double)stats.lb_remain_time/stats.lb_time*100);
    printf("  - LB Front (MWPM): %lld us (%.2f%% of LB)\n", stats.lb_front_time, (double)stats.lb_front_time/stats.lb_time*100);
    printf("Missing Edge Calc:   %lld us (%.2f%%)\n", stats.missing_time, (double)stats.missing_time/stats.total_time*100);
    printf("Recursion Calls:     %lld\n", stats.recursion_calls);
    printf("Pruning Calls:       %lld\n", stats.prun_calls);
    printf("Results Found:       %zu\n", M_ANS.size());
    printf("---------------------------------\n");
}