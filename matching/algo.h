#include "utility/utility.h"
// #include "graph/graph.h"

using namespace std;

struct TimeStats {
    long long total_time = 0;
    long long get_cand_time = 0;   // getExpandableMappings 耗时
    long long lb_time = 0;         // calLowerBound 耗时
    long long missing_time = 0;    // calPartialMissing 耗时
    long long recursion_calls = 0; // 递归调用次数计数
    long long prun_calls = 0;

    // calLowerBound 内部细分时间
    long long lb_identify_sets_time = 0; // 1. Identify Nq, Rq, Ng, Rg
    long long lb_precompute_delta_time = 0; // 3. Precompute delta_inner
    long long lb_remain_time = 0; // 4. Calculate LB_remain
    long long lb_front_time = 0; // 5. Calculate LB_front
    long long lb_inner_time = 0;
};

LabelID label2int(const std::string str, std::map<std::string, LabelID> &M) {
	if(M.find(str) == M.end()) M[str] = M.size();
	return M[str];
}

void load_graph(const std::string &input_graph, Graph *graph,
                std::map<std::string, LabelID> &vM, std::map<std::string, LabelID> &eM) {

    std::ifstream fin(input_graph);
    std::string line;

    while (std::getline(fin, line)) if (!line.empty() && line[0] == 't') break;

    if(!fin) {
        printf("!!! Cannot open graph file %s !!!\n", input_graph.c_str());
        assert(false);
        return;
    }

    std::istringstream head(line);
    char tchar; 
    std::string sharp, id;
    head >> tchar >> sharp >> id;
    assert(tchar == 't');

    std::vector<std::pair<ui, LabelID> > vertices;
    std::vector<std::pair<std::pair<ui, ui>, LabelID> > undirected_edges;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        char type = line[0];
        if (type == 't') break;

        std::istringstream iss(line);

        if(type == 'v') {
            char c;
            ui vid;
            std::string vlab;
            iss >> c >> vid >> vlab;
            vertices.emplace_back(vid, label2int(vlab, vM));
        }
        else if (type == 'e') {
            char c;
            ui u, v;
            std::string elab;
            iss >> c >> u >> v >> elab;
            if(u == v) continue;
            LabelID L = label2int(elab, eM);
            undirected_edges.emplace_back(std::make_pair(std::min(u, v), std::max(u, v)), L);
        }
    }

    std::sort(vertices.begin(), vertices.end());
    std::sort(undirected_edges.begin(), undirected_edges.end());

    vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());
    undirected_edges.erase(unique(undirected_edges.begin(), undirected_edges.end()), undirected_edges.end());

    std::vector<std::pair<std::pair<ui, ui>, LabelID> > edges;

    for(auto &e: undirected_edges) {
        edges.emplace_back(std::make_pair(e.first.first, e.first.second), e.second);
        edges.emplace_back(std::make_pair(e.first.second, e.first.first), e.second);
    }

    std::sort(edges.begin(), edges.end());

    graph->build_graph(id, vertices, edges);

// #ifndef NDEBUG
//     graph.print_graph();
// #endif

	return;
}

bool calVerticesFilter(const Graph *query_graph, const Graph *data_graph,
                       vector<vector<ui> > &candidates) {
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();

    candidates.clear();
    candidates.resize(qn);

    for (ui i = 0; i < qn; ++i) {
        int label_i = query_graph->getVertexLabel(i);

        for (ui j = 0; j < gn; ++j) {
            int label_j = data_graph->getVertexLabel(j);
            if (label_i == label_j) candidates[i].push_back(j);
        }

        if (candidates[i].empty())  return false;
    }

    return true;
}

void getExpandableMappings(const Graph *query_graph, const Graph *data_graph,
                     const vector<int> &mapped_q, const vector<int> &mapped_g,
                     const vector<vector<ui> > &candidates, const vector<pair<ui, ui>> &part_M, 
                     const unordered_set<pair<ui, ui>, PairHash> X, vector<pair<ui,ui> > &Mcand)
{
    Mcand.clear();

    // part_M = {}
    if(part_M.empty()) {
        // select 0 as starting vertex
        ui u = 0;
        assert(mapped_q[u] == -1);
        for(ui v : candidates[u]) {
            assert(mapped_g[v] == -1);
            if(X.count({u, v})) continue;
            Mcand.emplace_back(u, v);
        }
        return;
    }

    unordered_set<uint64_t> seen;
    for(auto p : part_M) {
        ui u = p.first;
        ui v = p.second;

        ui q_count;
        const ui* q_neighbors = query_graph->getVertexNeighbors(u, q_count);
        for(ui i = 0; i < q_count; ++i) {
            ui u2 = q_neighbors[i];
            if(mapped_q[u2] != -1) continue;

            ui g_count;
            const ui* g_neighbors = data_graph->getVertexNeighbors(v, g_count);
            for(ui j = 0; j < g_count; ++j) {
                ui v2 = g_neighbors[j];

                if(mapped_g[v2] != -1) continue;

                if(find(candidates[u2].begin(), candidates[u2].end(), v2) == candidates[u2].end())
                    continue;
                
                if(X.count({u2, v2})) continue;

                uint64_t key = (uint64_t(u2) << 32) | v2;
                if(seen.insert(key).second) {
                    Mcand.emplace_back(u2, v2);
                }
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

/**
 * @brief 求解最小权重完美匹配 (Minimum Weight Perfect Matching)
 * 逻辑参考自 Application.txt 中的 Hungarian 实现。
 * 
 * @param cost_matrix 代价矩阵，行表示查询图节点，列表示数据图节点
 * @return 最小权重之和
 */
ui solveMWPM(const vector<vector<ui>>& cost_matrix) {
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

    // 1. 初始化 (Initialization - Greedy Match)
    for (int i = 0; i < n; i++) {
        int min_val = INF_INT;
        for (int j = 0; j < m; j++) {
            if ((int)cost_matrix[i][j] < min_val) min_val = (int)cost_matrix[i][j];
        }
        lx[i] = min_val; // 左侧顶标初始化为行最小值
        for (int j = 0; j < m; j++) {
            if (my[j] == -1 && (int)cost_matrix[i][j] == lx[i]) {
                mx[i] = j;
                my[j] = i;
                break;
            }
        }
    }

    // 2. 增广 (Augmentation)
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
    ui total_cost = 0;
    for (int i = 0; i < n; i++) {
        if (mx[i] != -1) total_cost += cost_matrix[i][mx[i]];
    }
    return total_cost;
}

//  

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

ui computeDeltaInnerDoubledVector(const vector<ui> &Lq, const vector<ui> &Lg) {
    ui diff = 0;
    size_t sz = Lq.size();
    for (size_t i = 0; i < sz; ++i) {
        if (Lq[i] > Lg[i]) {
            diff += (Lq[i] - Lg[i]);
        }
    }
    return diff;
}

ui calLowerBound(const Graph *query_graph, const Graph *data_graph,
                 const vector<int> &mapped_q, const vector<int> &mapped_g,
                 const vector<vector<ui> > &candidates, const vector<pair<ui, ui>> &part_M, 
                 vector<pair<ui,ui> > &Mcand, TimeStats &stats, ui threshold)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();
    ui ln = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());

    Timer t_part;
    t_part.restart();

    // 1. Identify Nq, Ng, Rq, Rg
    // build adj_Mcand
    vector<vector<ui>> adj_Mcand(qn);
    vector<uint8_t> in_Nq(qn, 0);
    vector<uint8_t> in_Ng(gn, 0);
    vector<ui> Nq, Ng, Rq, Rg;
    Nq.reserve(Mcand.size());
    Ng.reserve(Mcand.size());
    Rq.reserve(qn);
    Rg.reserve(gn);

    for (const auto& p : Mcand) {
        ui u = p.first;
        ui v = p.second;
        adj_Mcand[u].push_back(v);

        if (!in_Nq[u]) {
            in_Nq[u] = true;
            Nq.push_back(u);
        }
        if (!in_Ng[v]) {
            in_Ng[v] = true;
            Ng.push_back(v);
        }
    }

    for (ui u = 0; u < qn; ++u) if (mapped_q[u] == -1 && !in_Nq[u]) Rq.push_back(u);
    for (ui v = 0; v < gn; ++v) if (mapped_g[v] == -1 && !in_Ng[v]) Rg.push_back(v);
    stats.lb_identify_sets_time += t_part.elapsed();

    // 2. Precompute delta_inner components for all unmapped
    t_part.restart();
    vector<vector<ui>> Lq_sets(qn, vector<ui>(ln, 0));
    vector<vector<ui>> Lg_sets(gn, vector<ui>(ln, 0));

    for(ui u : Nq) getUnmappedNeighborLabelsVector(query_graph, u, mapped_q, Lq_sets[u]);
    for(ui u : Rq) getUnmappedNeighborLabelsVector(query_graph, u, mapped_q, Lq_sets[u]);
    for(ui v : Ng) getUnmappedNeighborLabelsVector(data_graph, v, mapped_g, Lg_sets[v]);
    for(ui v : Rg) getUnmappedNeighborLabelsVector(data_graph, v, mapped_g, Lg_sets[v]);
    stats.lb_precompute_delta_time += t_part.elapsed();

    // 3. Calculate LB_remain (Greedy relaxation for Rq) [MOD](double)
    // ui LB_remain = 0;
    t_part.restart();
    ui LB_remain_doubled = 0;
    for (ui u : Rq) {
        ui min_u_cost = (ui)-1;
        // LabelID lu = query_graph->getVertexLabel(u);

        for (ui v : candidates[u]) if(mapped_g[v] == -1) {
            min_u_cost = min(min_u_cost, computeDeltaInnerDoubledVector(Lq_sets[u], Lg_sets[v]));
        }
        if (min_u_cost != (ui)-1) LB_remain_doubled += min_u_cost;
        // else {
        //     LB_remain_doubled += threshold * 2;
        // }
    }
    stats.lb_remain_time += t_part.elapsed();

    // 4. Calculate LB_front (MWPM for Nq)
    t_part.restart();
    ui LB_front_doubled = 0;
    vector<pair<LabelID, ui>> Nq_labeled;
    for (ui u : Nq) Nq_labeled.emplace_back(query_graph->getVertexLabel(u), u);
    sort(Nq_labeled.begin(), Nq_labeled.end());
    
    size_t ptr = 0;
    while(ptr < Nq_labeled.size()) {
        LabelID label = Nq_labeled[ptr].first;
        size_t qtr = ptr;
        while(qtr < Nq_labeled.size() && Nq_labeled[qtr].first == Nq_labeled[ptr].first) qtr++;

        vector<ui> u_list;
        for(size_t i = ptr; i < qtr; i++) u_list.emplace_back(Nq_labeled[i].second);
        vector<ui> v_list;
        for (ui v : Ng) if (data_graph->getVertexLabel(v) == label) v_list.push_back(v);

        // calculate w_ext(u) (doubled): 2 * |Nq(u) \cap Mq| + min delta_inner_doubled
        vector<ui> w_ext_doubled(u_list.size());
        for (size_t i = 0; i < u_list.size(); ++i) {
            ui u = u_list[i];
            ui cross_fixed = 0;
            ui deg;
            const ui* neighbors = query_graph->getVertexNeighbors(u, deg);
            for (ui j = 0; j < deg; ++j) if (mapped_q[neighbors[j]] != -1) cross_fixed++;

            ui min_inner_doubled = 0; 
            bool found_rg = false;
            for (ui v : Rg) {
                if (data_graph->getVertexLabel(v) == label) {
                    ui d = computeDeltaInnerDoubledVector(Lq_sets[u], Lg_sets[v]);
                    if (!found_rg || d < min_inner_doubled) { min_inner_doubled = d; found_rg = true; }
                }
            }
            w_ext_doubled[i] = (2 * cross_fixed) + min_inner_doubled;
        }

        // 构建代价矩阵: 行 = u_list, 列 = v_list (Ng) + 虚拟节点 (代表 Rg)
        size_t row_n = u_list.size();
        size_t col_n = v_list.size() + row_n;
        vector<vector<ui>> cost_matrix(row_n, vector<ui>(col_n));

        for (size_t i = 0; i < row_n; ++i) {
            ui u = u_list[i];
            // 真实节点 Ng 的代价: 2 * delta_cross + delta_inner_doubled
            for (size_t j = 0; j < v_list.size(); ++j) {
                ui v = v_list[j];
                ui d_cross = 0;
                ui u_deg;
                const ui* u_neighbors = query_graph->getVertexNeighbors(u, u_deg);
                for (ui k = 0; k < u_deg; ++k) {
                    ui un = u_neighbors[k];
                    if (mapped_q[un] != -1) {
                        if (!data_graph->hasEdge(v, (ui)mapped_q[un])) d_cross++;
                    }
                }
                cost_matrix[i][j] = (2 * d_cross) + computeDeltaInnerDoubledVector(Lq_sets[u], Lg_sets[v]);
            }
            // 虚拟节点的代价 (映射到 Rg)
            for (size_t j = 0; j < row_n; ++j) {
                cost_matrix[i][v_list.size() + j] = w_ext_doubled[i];
            }
        }
        LB_front_doubled += solveMWPM(cost_matrix);
        ptr = qtr;
    }
    stats.lb_front_time += t_part.elapsed();

    return (LB_front_doubled + LB_remain_doubled + 1) / 2;
}

void DFS_Approximate(const Graph *query_graph, const Graph *data_graph,
                     const vector<vector<ui>> &candidates, unordered_set<pair<ui, ui>, PairHash> X,
                     vector<int> &mapped_q, vector<int> &mapped_g,
                     vector<pair<ui, ui>> &part_M, vector<vector<pair<ui,ui>>> &M_ANS,
                     ui threshold, ui missing, TimeStats &stats)
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
#endif
    stats.recursion_calls++;

    ui qn = query_graph->getVerticesCount();
    assert(missing == calPartialMissing(query_graph, data_graph, part_M));

    if (part_M.size() == qn) {
        assert(missing <= threshold);
        M_ANS.push_back(part_M);
        return;
    }

    Timer t_cand;
    vector<pair<ui,ui> > Mcand;
    getExpandableMappings(query_graph, data_graph, mapped_q, mapped_g, candidates, part_M, X, Mcand);
    stats.get_cand_time += t_cand.elapsed();

    if(Mcand.empty()) return;

    Timer t_lb;
    ui lb = 0;
#ifdef LOWER_BOUND
    lb = calLowerBound(query_graph, data_graph, mapped_q, mapped_g, candidates, part_M, Mcand, stats, threshold);
#endif
    stats.lb_time += t_lb.elapsed();
    if(missing + lb > threshold) {
        // pruning
        stats.prun_calls++;
        return;
    }

    for (auto p : Mcand) {
        ui u = p.first;
        ui v = p.second;

        assert(X.count({u, v}) == 0);
        assert(mapped_q[u] == -1);
        assert(mapped_g[v] == -1);
        mapped_q[u] = v;
        mapped_g[v] = u;
        part_M.emplace_back(u, v);

        Timer t_miss;
        ui new_missing = calPartialMissing(query_graph, data_graph, part_M);
        stats.missing_time += t_miss.elapsed();
        if (new_missing <= threshold) {
            DFS_Approximate(query_graph, data_graph, candidates, X,
                            mapped_q, mapped_g, part_M, M_ANS, threshold, new_missing, stats);
        }

        mapped_q[u] = -1;
        mapped_g[v] = -1;
        part_M.pop_back();

        X.insert({u, v});
    }
}

// ============================================================
// Top-level function: Approximate_Matching_v2
// ============================================================
void Approximate_Matching_v2(const Graph *query_graph, const Graph *data_graph,
                             vector<vector<ui> > &candidates, 
                             vector<ui> &order, // unused in v2
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
    
    Timer t_total;
    DFS_Approximate(query_graph, data_graph, candidates, X,
                    mapped_q, mapped_g, part_M, M_ANS, threshold, 0, stats);
    stats.total_time = t_total.elapsed();

    printf("\n--- Matching v2 Time Analysis ---\n");
    printf("Total Time:          %lld us\n", stats.total_time);
    printf("Mcand Generation:    %lld us (%.2f%%)\n", stats.get_cand_time, (double)stats.get_cand_time/stats.total_time*100);
    printf("LowerBound Pruning:  %lld us (%.2f%%)\n", stats.lb_time, (double)stats.lb_time/stats.total_time*100);
    printf("  - Identify Sets:   %lld us (%.2f%% of LB)\n", stats.lb_identify_sets_time, (double)stats.lb_identify_sets_time/stats.lb_time*100);
    printf("  - Precompute Delta:%lld us (%.2f%% of LB)\n", stats.lb_precompute_delta_time, (double)stats.lb_precompute_delta_time/stats.lb_time*100);
    printf("  - LB Remain:       %lld us (%.2f%% of LB)\n", stats.lb_remain_time, (double)stats.lb_remain_time/stats.lb_time*100);
    printf("  - LB Front (MWPM): %lld us (%.2f%% of LB)\n", stats.lb_front_time, (double)stats.lb_front_time/stats.lb_time*100);
    printf("  - LB Inner:        %lld us (%.2f%% of LB)\n", stats.lb_inner_time, (double)stats.lb_inner_time/stats.lb_time*100);
    printf("Missing Edge Calc:   %lld us (%.2f%%)\n", stats.missing_time, (double)stats.missing_time/stats.total_time*100);
    printf("Recursion Calls:     %lld\n", stats.recursion_calls);
    printf("Results Found:       %zu\n", M_ANS.size());
    printf("---------------------------------\n");
}