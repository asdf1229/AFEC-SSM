#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_CANDIDATE_INDEX_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_CANDIDATE_INDEX_H_

#ifndef CDE_BLACK_WHITE_INSIDE_WORKSPACE
#error "This internal header must be included from cde_black_white/context.h inside Workspace."
#endif

    unsigned long long adjKey(ui data_vertex, ui query_neighbor) const
    {
        return ((unsigned long long)data_vertex << 32) | (unsigned long long)query_neighbor;
    }

    void buildAdjIndex()
    {
        candidate_adj_index.clear();
        candidate_adj_index.resize(qn);
        candidate_adj_pool.clear();

        size_t total_candidate_count = 0;
        for (ui u = 0; u < qn; ++u) {
            total_candidate_count += (size_t)candidates[u].size();
            candidate_adj_index[u].reserve((size_t)candidates[u].size() * q_neighbors[u].size());
        }
        candidate_source_buffer.reserve(total_candidate_count);
        candidate_result_buffer.reserve(total_candidate_count);
        candidate_intersection_buffer.reserve(total_candidate_count);
        candidate_range_buffer.reserve(qn);

        for (ui u = 0; u < qn; ++u) {
            for (ui v : candidates[u]) {
                ui degree = 0; const ui *neighbors = data_graph->getVertexNeighbors(v, degree);
                for (ui query_neighbor : q_neighbors[u]) {
                    size_t begin = candidate_adj_pool.size();
                    for (ui i = 0; i < degree; ++i) {
                        ui data_neighbor = neighbors[i];
                        if (candidates[query_neighbor].contains(data_neighbor)) {
                            candidate_adj_pool.push_back(data_neighbor);
                        }
                    }

                    ui len = (ui)(candidate_adj_pool.size() - begin);
                    if (len == 0) continue;

                    pair<size_t, ui> range(begin, len);
                    candidate_adj_index[u].emplace(adjKey(v, query_neighbor), range);
                }
            }
        }
    }

    // return to_data set
    const pair<size_t, ui> *findAdjRange(ui from_query, ui from_data, ui to_query) const
    {
        if (from_query >= candidate_adj_index.size()) return nullptr;
        const auto &index = candidate_adj_index[from_query];
        auto it = index.find(adjKey(from_data, to_query));
        if (it == index.end()) return nullptr;
        return &it->second;
    }

    const ui *rangeBegin(const pair<size_t, ui> &range) const
    {
        return candidate_adj_pool.data() + range.first;
    }

    const ui *rangeEnd(const pair<size_t, ui> &range) const
    {
        return rangeBegin(range) + range.second;
    }

    bool rangeHas(const pair<size_t, ui> &range, ui value) const
    {
        return std::binary_search(rangeBegin(range), rangeEnd(range), value);
    }

    bool candAdjacent(ui from_query, ui from_data, ui to_query, ui to_data)
    {
        stats.candidate_edge_check_calls++;
        const pair<size_t, ui> *range = findAdjRange(from_query, from_data, to_query);
        if (range == nullptr) {
            stats.candidate_range_misses++;
            return false;
        }
        stats.candidate_range_hits++;
        return rangeHas(*range, to_data);
    }

    bool hasDataEdge(ui u, ui v)
    {
        stats.graph_has_edge_checks++;
        return data_graph->hasEdge(u, v);
    }

    bool anchorAdjacent(ui anchor_query, ui anchor_data, ui target_query, ui target_data)
    {
        stats.candidate_edge_check_calls++;
        const pair<size_t, ui> *range = findAdjRange(anchor_query, anchor_data, target_query);
        if (range == nullptr) {
            stats.candidate_range_misses++;
            return false;
        }
        stats.candidate_range_hits++;
        return rangeHas(*range, target_data);
    }


#endif
