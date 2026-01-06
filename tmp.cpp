    // // 1. Identify Nq, Rq (Query side)
    // vector<ui> Nq, Rq;
    // for (ui u = 0; u < qn; ++u) {
    //     if (mapped_q[u] != -1) continue;
        
    //     bool is_frontier = false;
    //     ui deg;
    //     const ui* neighbors = query_graph->getVertexNeighbors(u, deg);
    //     for (ui i = 0; i < deg; ++i) {
    //         if (mapped_q[neighbors[i]] != -1) {
    //             is_frontier = true;
    //             break;
    //         }
    //     }
    //     if (is_frontier) Nq.push_back(u);
    //     else Rq.push_back(u);
    // }

    // // 2. Identify Ng, Rg (Data side)
    // vector<ui> Ng, Rg;
    // for (ui v = 0; v < gn; ++v) {
    //     if (mapped_g[v] != -1) continue;

    //     bool is_frontier = false;
    //     ui deg;
    //     const ui* neighbors = data_graph->getVertexNeighbors(v, deg);
    //     for (ui i = 0; i < deg; ++i) {
    //         if (mapped_g[neighbors[i]] != -1) {
    //             is_frontier = true;
    //             break;
    //         }
    //     }
    //     if (is_frontier) Ng.push_back(v);
    //     else Rg.push_back(v);
    // }
    // 4. Calculate LB_front (MWPM for Nq)
    // t_part.restart();
    // ui LB_front_doubled = 0;
    // map<LabelID, vector<ui>> Nq_by_label;
    // for (ui u : Nq) Nq_by_label[query_graph->getVertexLabel(u)].push_back(u);
    
    //     for (auto const& [label, u_list] : Nq_by_label) {
    //     vector<ui> v_list;
    //     for (ui v : Ng) if (data_graph->getVertexLabel(v) == label) v_list.push_back(v);

    //     // 计算 w_ext(u) 的翻倍值: 2 * |Nq(u) \cap Mq| + min delta_inner_doubled
    //     vector<ui> w_ext_doubled(u_list.size());
    //     for (size_t i = 0; i < u_list.size(); ++i) {
    //         ui u = u_list[i];
    //         ui cross_fixed = 0;
    //         ui deg;
    //         const ui* neighbors = query_graph->getVertexNeighbors(u, deg);
    //         for (ui j = 0; j < deg; ++j) if (mapped_q[neighbors[j]] != -1) cross_fixed++;

    //         ui min_inner_doubled = 0; 
    //         bool found_rg = false;
    //         for (ui v : Rg) {
    //             if (data_graph->getVertexLabel(v) == label) {
    //                 ui d = computeDeltaInnerDoubled(Lq_sets[u], Lg_sets[v]);
    //                 if (!found_rg || d < min_inner_doubled) { min_inner_doubled = d; found_rg = true; }
    //             }
    //         }
    //         w_ext_doubled[i] = (2 * cross_fixed) + min_inner_doubled;
    //     }

    //     // 构建代价矩阵: 行 = u_list, 列 = v_list (Ng) + 虚拟节点 (代表 Rg)
    //     size_t row_n = u_list.size();
    //     size_t col_n = v_list.size() + row_n;
    //     vector<vector<ui>> cost_matrix(row_n, vector<ui>(col_n));

    //     for (size_t i = 0; i < row_n; ++i) {
    //         ui u = u_list[i];
    //         // 真实节点 Ng 的代价: 2 * delta_cross + delta_inner_doubled
    //         for (size_t j = 0; j < v_list.size(); ++j) {
    //             ui v = v_list[j];
    //             ui d_cross = 0;
    //             ui u_deg;
    //             const ui* u_neighbors = query_graph->getVertexNeighbors(u, u_deg);
    //             for (ui k = 0; k < u_deg; ++k) {
    //                 ui un = u_neighbors[k];
    //                 if (mapped_q[un] != -1) {
    //                     if (!data_graph->hasEdge(v, (ui)mapped_q[un])) d_cross++;
    //                 }
    //             }
    //             cost_matrix[i][j] = (2 * d_cross) + computeDeltaInnerDoubled(Lq_sets[u], Lg_sets[v]);
    //         }
    //         // 虚拟节点的代价 (映射到 Rg)
    //         for (size_t j = 0; j < row_n; ++j) {
    //             cost_matrix[i][v_list.size() + j] = w_ext_doubled[i];
    //         }
    //     }
    //     LB_front_doubled += solveMWPM(cost_matrix);
    // }

    // // delta_inner: 1/2 * |Lq \ Lg|
// double computeDeltaInner(const map<LabelID, ui> &Lq, const map<LabelID, ui> &Lg) {
//     ui diff = 0;
//     for (auto const& [label, count] : Lq) {
//         ui g_count = 0;
//         if (Lg.find(label) != Lg.end()) g_count = Lg.at(label);
//         if (count > g_count) diff += (count - g_count);
//     }
//     return diff / 2.0; // Factor 1/2
// }

// ui computeDeltaInnerDoubled(const map<LabelID, ui> &Lq, const map<LabelID, ui> &Lg) {
//     ui diff = 0;
//     for (auto const& [label, count] : Lq) {
//         ui g_count = 0;
//         auto it = Lg.find(label);
//         if (it != Lg.end()) g_count = it->second;
//         if (count > g_count) diff += (count - g_count);
//     }
//     return diff; 
// }