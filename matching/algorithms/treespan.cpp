#include "matching/algorithms/treespan.h"

#include "utility/utility.h"
#include "utility/mybitset.h"

using namespace std;

struct TreeSpanSolver::Impl {
    TreeStats &stats;

    explicit Impl(TreeStats &stats_ref)
        : stats(stats_ref), query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr)
    {}

    bool init(const Graph *q, const Graph *d, ui match_threshold)
    {
        Timer t_init;
        t_init.restart();

        resetAll();

        query_graph = q;
        data_graph = d;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();
        label_count = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());

        if (qn == 0 || gn == 0) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        mapped_q.assign(qn, -1);
        used_data_vertices.clear();
        used_data_vertices.reserve(qn);
        used_data_flag.assign(gn, 0);
        undo_stack.clear();
        undo_stack.reserve((size_t)qn * 2);
        candidates.assign(qn, MyBitset(gn));

        buildQueryMatrix();
        Timer t_filter;
        t_filter.restart();
        initGlobalLabelCounts(query_graph, Lq_counts);
        initGlobalLabelCounts(data_graph, Lg_counts);
        initDataLabelBuckets();

        bool filter_ok = calVerticesFilter();
        stats.filter_time = t_filter.elapsed();
        if (!filter_ok) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        cacheCandidateCounts();

        Timer t_index;
        t_index.restart();
        buildAdjIndex();
        stats.index_time = t_index.elapsed();

        initOrderingStatistics();
        cacheDirectedEdgeScores();

        root = selectRoot();

        if (!buildInitialTree(initial_tree)) {
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
        result_keys.clear();
        std::fill(mapped_q.begin(), mapped_q.end(), -1);
        std::fill(used_data_flag.begin(), used_data_flag.end(), 0);
        used_data_vertices.clear();
        undo_stack.clear();

        const QISequence *seq = getQISequence(initial_tree);
        if (seq == nullptr) {
            stats.search_time = t_search.elapsed();
            return;
        }

        ui start_u = seq->S[0];
        size_t root_mark = mark();
        for (ui v : candidates[start_u]) {
            if (outputLimitReached()) break;

            setMap(start_u, (int)v);
            pushUsed(v);

            SimSearchOnDemand(1, initial_tree, 0);

            rollback(root_mark);
            assert(mark() == root_mark);

            if (outputLimitReached()) break;
        }

        stats.search_time = t_search.elapsed();
        stats.sequences_count = (long long)sequence_cache.size();
    }

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        long long search_accounted_time = stats.enum_time + stats.verify_time;
        long long search_other_time = stats.search_time > search_accounted_time
            ? stats.search_time - search_accounted_time : 0;

        printf("\n--- TreeSpan Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("  - Adj Index Time:  %.4lf ms (%.2f%% of Init)\n", stats.index_time / 1000.0, pct(stats.index_time, stats.init_time));
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
        printf("  - Q_T OD Time:     %.4lf ms (%.2f%% of Search)\n", stats.enum_time / 1000.0, pct(stats.enum_time, stats.search_time));
        printf("  - Verify Time:     %.4lf ms (%.2f%% of Search)\n", stats.verify_time / 1000.0, pct(stats.verify_time, stats.search_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.search_time));
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Replacement Calls:   %lld\n", stats.replacement_calls);
        printf("QISequences Cached:  %lld\n", stats.sequences_count);
        printf("OD Calls:            %lld\n", stats.enum_call_count);
        printf("Range Hits/Misses:   %lld / %lld\n",
            stats.candidate_range_hits, stats.candidate_range_misses);
        printf("Range Edge Checks:   %lld\n", stats.candidate_edge_check_calls);
        printf("Duplicate Results:   %lld\n", stats.duplicate_results);
        printf("Results Found:       %zu\n", stats.result_count);
#if MATCH_OUTPUT_LIMIT > 0
        printf("Output Limit:        %zu%s\n",
            (size_t)MATCH_OUTPUT_LIMIT,
            stats.output_limit_reached ? " (reached)" : "");
#endif
        printf("-----------------------------------------------------------\n");
    }

    static const ui INVALID_EDGE = (ui)-1;

    struct DirectedScore {
        unsigned long long numerator = 0;
        unsigned long long denominator = 1;
    };

    struct QEdge {
        ui u = 0;
        ui v = 0;
        ui id = 0;
        DirectedScore u_to_v; // u is the new vertex and v is the anchor.
        DirectedScore v_to_u; // v is the new vertex and u is the anchor.
    };

    struct LabelCount {
        LabelID label = 0;
        ui count = 0;
    };

    struct TreeState {
        vector<ui> edges; // Ordered tree edges, each edge adding one new vertex.
        set<ui> R;       // Edges excluded by previous replacements.
    };

    struct QISequence {
        vector<ui> S;                  // Vertex visit order.
        vector<ui> sEdgeIds;           // sEdgeIds[i] connects S[i] to the prefix.
        vector<vector<ui>> bEdgeIds;   // Non-tree edges from S[i] to the prefix.
        set<ui> R;
    };

    enum UndoKind : unsigned char {
        UNDO_MAPPED_Q = 0,
        UNDO_USED_DATA_SIZE = 1
    };

    struct UndoRecord {
        UndoKind kind = UNDO_MAPPED_Q;
        ui u = 0;
        int old_mapped_q = -1;
        size_t old_size = 0;
    };

    struct DisjointSet {
        vector<int> parent;
        vector<int> rank;

        explicit DisjointSet(ui n) : parent(n), rank(n, 0)
        {
            for (ui i = 0; i < n; ++i) parent[i] = (int)i;
        }

        int find(int x)
        {
            if (parent[x] != x) parent[x] = find(parent[x]);
            return parent[x];
        }

        bool unite(ui a, ui b)
        {
            int ra = find((int)a);
            int rb = find((int)b);
            if (ra == rb) return false;
            if (rank[ra] < rank[rb]) swap(ra, rb);
            parent[rb] = ra;
            if (rank[ra] == rank[rb]) rank[ra]++;
            return true;
        }
    };

    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;

    ui threshold = 0;
    ui qn = 0;
    ui gn = 0;
    ui label_count = 0;
    ui root = 0;

    vector<MyBitset> candidates;
    vector<ui> candidate_counts;
    vector<vector<LabelCount>> Lq_counts;
    vector<vector<LabelCount>> Lg_counts;
    vector<vector<ui>> data_vertices_by_label;
    vector<ui> data_label_frequency;
    vector<vector<ui>> data_edge_label_frequency;

    vector<vector<char>> q_matrix;
    vector<vector<int>> q_edge_id;
    vector<QEdge> all_q_edges;

    vector<int> mapped_q;
    vector<ui> used_data_vertices;
    vector<unsigned char> used_data_flag;
    vector<UndoRecord> undo_stack;

    vector<absl::flat_hash_map<unsigned long long, pair<size_t, ui>>> candidate_adj_index;
    vector<ui> candidate_adj_pool;

    TreeState initial_tree;
    map<string, QISequence> sequence_cache;
    unordered_set<string> result_keys;

    bool outputLimitReached() const
    {
        return (size_t)MATCH_OUTPUT_LIMIT > 0 &&
            stats.result_count >= (size_t)MATCH_OUTPUT_LIMIT;
    }

    void noteOutputLimitIfReached()
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
        }
    }

    void resetAll()
    {
        stats = TreeStats();
        query_graph = nullptr;
        data_graph = nullptr;
        results_ptr = nullptr;
        threshold = 0;
        qn = 0;
        gn = 0;
        label_count = 0;
        root = 0;

        candidates.clear();
        candidate_counts.clear();
        Lq_counts.clear();
        Lg_counts.clear();
        data_vertices_by_label.clear();
        data_label_frequency.clear();
        data_edge_label_frequency.clear();
        q_matrix.clear();
        q_edge_id.clear();
        all_q_edges.clear();
        mapped_q.clear();
        used_data_vertices.clear();
        used_data_flag.clear();
        undo_stack.clear();
        candidate_adj_index.clear();
        candidate_adj_pool.clear();
        initial_tree = TreeState();
        sequence_cache.clear();
        result_keys.clear();
    }

    unsigned long long adjKey(ui data_vertex, ui query_neighbor) const
    {
        return ((unsigned long long)data_vertex << 32) |
            (unsigned long long)query_neighbor;
    }

    void buildAdjIndex()
    {
        candidate_adj_index.clear();
        candidate_adj_index.resize(qn);
        candidate_adj_pool.clear();

        for (ui u = 0; u < qn; ++u) {
            ui query_degree = 0;
            query_graph->getVertexNeighbors(u, query_degree);
            candidate_adj_index[u].reserve(
                (size_t)candidate_counts[u] * query_degree);
        }

        for (ui u = 0; u < qn; ++u) {
            ui query_degree = 0;
            const ui *query_neighbors =
                query_graph->getVertexNeighbors(u, query_degree);

            for (ui v : candidates[u]) {
                ui data_degree = 0;
                const ui *data_neighbors =
                    data_graph->getVertexNeighbors(v, data_degree);

                for (ui i = 0; i < query_degree; ++i) {
                    ui query_neighbor = query_neighbors[i];
                    size_t begin = candidate_adj_pool.size();

                    for (ui j = 0; j < data_degree; ++j) {
                        ui data_neighbor = data_neighbors[j];
                        if (candidates[query_neighbor].contains(data_neighbor)) {
                            candidate_adj_pool.push_back(data_neighbor);
                        }
                    }

                    ui len = (ui)(candidate_adj_pool.size() - begin);
                    if (len == 0) continue;

                    candidate_adj_index[u].emplace(
                        adjKey(v, query_neighbor), pair<size_t, ui>(begin, len));
                }
            }
        }
    }

    const pair<size_t, ui> *findAdjRange(ui from_query, ui from_data,
        ui to_query) const
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

    bool anchorAdjacent(ui anchor_query, ui anchor_data, ui target_query,
        ui target_data)
    {
        stats.candidate_edge_check_calls++;
        const pair<size_t, ui> *range =
            findAdjRange(anchor_query, anchor_data, target_query);
        if (range == nullptr) {
            stats.candidate_range_misses++;
            return false;
        }
        stats.candidate_range_hits++;
        return rangeHas(*range, target_data);
    }

    bool isDataVertexUsed(ui v) const
    {
        assert(v < used_data_flag.size());
        return used_data_flag[v] != 0;
    }

    size_t mark() const
    {
        return undo_stack.size();
    }

    void rollback(size_t undo_mark)
    {
        assert(undo_mark <= undo_stack.size());
        while (undo_stack.size() > undo_mark) {
            UndoRecord undo = std::move(undo_stack.back());
            undo_stack.pop_back();

            switch (undo.kind) {
            case UNDO_MAPPED_Q:
                mapped_q[undo.u] = undo.old_mapped_q;
                break;
            case UNDO_USED_DATA_SIZE:
                for (size_t i = undo.old_size;
                    i < used_data_vertices.size(); ++i) {
                    ui used_v = used_data_vertices[i];
                    assert(used_v < used_data_flag.size());
                    used_data_flag[used_v] = 0;
                }
                used_data_vertices.resize(undo.old_size);
                break;
            }
        }
    }

    void setMap(ui u, int value)
    {
        assert(u < mapped_q.size());
        assert(value >= 0 && (ui)value < gn);
        assert(mapped_q[u] == -1);
        UndoRecord undo;
        undo.kind = UNDO_MAPPED_Q;
        undo.u = u;
        undo.old_mapped_q = mapped_q[u];
        undo_stack.push_back(std::move(undo));
        mapped_q[u] = value;
    }

    void pushUsed(ui v)
    {
        assert(v < used_data_flag.size());
        assert(!isDataVertexUsed(v));

        UndoRecord undo;
        undo.kind = UNDO_USED_DATA_SIZE;
        undo.old_size = used_data_vertices.size();
        undo_stack.push_back(std::move(undo));

        used_data_vertices.push_back(v);
        used_data_flag[v] = 1;
    }

    void buildQueryMatrix()
    {
        q_matrix.assign(qn, vector<char>(qn, 0));
        q_edge_id.assign(qn, vector<int>(qn, -1));
        all_q_edges.clear();

        for (ui u = 0; u < qn; ++u) {
            ui deg;
            const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                q_matrix[u][v] = 1;
                if (u < v) {
                    QEdge e;
                    e.u = u;
                    e.v = v;
                    e.id = (ui)all_q_edges.size();
                    all_q_edges.push_back(e);
                    q_edge_id[u][v] = (int)e.id;
                    q_edge_id[v][u] = (int)e.id;
                }
            }
        }
    }

    void initOrderingStatistics()
    {
        data_label_frequency.assign(label_count, 0);
        data_edge_label_frequency.assign(label_count, vector<ui>(label_count, 0));

        for (ui u = 0; u < gn; ++u) {
            LabelID label = data_graph->getVertexLabel(u);
            if (label >= 0 && (ui)label < label_count) {
                data_label_frequency[(ui)label]++;
            }
        }

        for (ui u = 0; u < gn; ++u) {
            LabelID lu = data_graph->getVertexLabel(u);
            if (lu < 0 || (ui)lu >= label_count) continue;

            ui deg;
            const ui *nbrs = data_graph->getVertexNeighbors(u, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (u >= v) continue;

                LabelID lv = data_graph->getVertexLabel(v);
                if (lv < 0 || (ui)lv >= label_count) continue;

                data_edge_label_frequency[(ui)lu][(ui)lv]++;
                data_edge_label_frequency[(ui)lv][(ui)lu]++;
            }
        }
    }

    DirectedScore makeDirectedScore(ui new_u, ui anchor) const
    {
        DirectedScore score;
        score.numerator =
            (unsigned long long)candidate_counts[new_u] *
            (unsigned long long)dataEdgeLabelFrequency(new_u, anchor);
        score.denominator = (unsigned long long)dataLabelFrequency(new_u);
        return score;
    }

    void cacheDirectedEdgeScores()
    {
        for (QEdge &edge : all_q_edges) {
            edge.u_to_v = makeDirectedScore(edge.u, edge.v);
            edge.v_to_u = makeDirectedScore(edge.v, edge.u);
        }
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<LabelCount>> &counts)
    {
        ui n = g->getVerticesCount();
        counts.assign(n, vector<LabelCount>());
        vector<ui> label_counts(label_count, 0);
        vector<ui> touched_labels;

        for (ui u = 0; u < n; ++u) {
            ui deg;
            const ui *neighbors = g->getVertexNeighbors(u, deg);
            counts[u].reserve(std::min(deg, label_count));
            touched_labels.clear();

            for (ui i = 0; i < deg; ++i) {
                LabelID label = g->getVertexLabel(neighbors[i]);
                if (label >= 0 && (ui)label < label_count) {
                    ui label_idx = (ui)label;
                    if (label_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    label_counts[label_idx]++;
                }
            }

            sort(touched_labels.begin(), touched_labels.end());
            for (ui label : touched_labels) {
                counts[u].push_back({ (LabelID)label, label_counts[label] });
                label_counts[label] = 0;
            }
        }
    }

    void initDataLabelBuckets()
    {
        data_vertices_by_label.assign(label_count, vector<ui>());
        for (ui v = 0; v < gn; ++v) {
            LabelID label = data_graph->getVertexLabel(v);
            if (label >= 0 && (ui)label < label_count) {
                data_vertices_by_label[(ui)label].push_back(v);
            }
        }
    }

    ui computeDelta(ui u, ui v) const
    {
        ui diff = 0;
        const vector<LabelCount> &query_counts = Lq_counts[u];
        const vector<LabelCount> &data_counts = Lg_counts[v];
        size_t data_idx = 0;

        for (const auto &query_count : query_counts) {
            while (data_idx < data_counts.size() && data_counts[data_idx].label < query_count.label) {
                data_idx++;
            }

            ui data_count = 0;
            if (data_idx < data_counts.size() && data_counts[data_idx].label == query_count.label) {
                data_count = data_counts[data_idx].count;
            }

            if (query_count.count > data_count) {
                diff += query_count.count - data_count;
                if (diff > threshold) {
                    return diff;
                }
            }
        }
        return diff;
    }

    bool calVerticesFilter()
    {
        for (ui u = 0; u < qn; ++u) {
            LabelID label_u = query_graph->getVertexLabel(u);
            if (label_u < 0 || (ui)label_u >= data_vertices_by_label.size()) {
                return false;
            }
            for (ui v : data_vertices_by_label[(ui)label_u]) {
                if (query_graph->getVertexDegree(u) > data_graph->getVertexDegree(v) + threshold) continue;
                if (computeDelta(u, v) <= threshold) candidates[u].insert(v);
            }
            if (candidates[u].empty()) return false;
        }
        return true;
    }

    void cacheCandidateCounts()
    {
        candidate_counts.assign(qn, 0);
        for (ui u = 0; u < qn; ++u) {
            candidate_counts[u] = (ui)candidates[u].size();
        }
    }

    ui selectRoot()
    {
        ui best_u = 0;
        ui best_cands = candidate_counts[0];
        ui best_degree = query_graph->getVertexDegree(0);

        for (ui u = 1; u < qn; ++u) {
            ui cand_size = candidate_counts[u];
            ui degree = query_graph->getVertexDegree(u);
            if (cand_size < best_cands || (cand_size == best_cands && degree > best_degree)) {
                best_u = u;
                best_cands = cand_size;
                best_degree = degree;
            }
        }
        return best_u;
    }

    ui dataLabelFrequency(ui u) const
    {
        LabelID label = query_graph->getVertexLabel(u);
        if (label < 0 || (ui)label >= data_label_frequency.size()) {
            return 1;
        }
        return max((ui)1, data_label_frequency[(ui)label]);
    }

    ui dataEdgeLabelFrequency(ui u, ui anchor) const
    {
        LabelID lu = query_graph->getVertexLabel(u);
        LabelID la = query_graph->getVertexLabel(anchor);
        if (lu < 0 || la < 0 ||
            (ui)lu >= data_edge_label_frequency.size() ||
            (ui)la >= data_edge_label_frequency[(ui)lu].size()) {
            return 0;
        }
        return data_edge_label_frequency[(ui)lu][(ui)la];
    }

    bool edgeCrossesVisited(ui edge_id, const vector<char> &visited,
        ui &new_u, ui &anchor) const
    {
        const QEdge &edge = all_q_edges[edge_id];
        bool u_vis = visited[edge.u] != 0;
        bool v_vis = visited[edge.v] != 0;
        if (u_vis == v_vis) return false;

        if (u_vis) {
            anchor = edge.u;
            new_u = edge.v;
        }
        else {
            anchor = edge.v;
            new_u = edge.u;
        }
        return true;
    }

    bool isBetterDirectedEdge(ui edge_id, ui new_u, ui anchor,
        ui best_edge_id, ui best_new_u, ui best_anchor) const
    {
        if (best_edge_id == INVALID_EDGE) return true;

        if (isDirectedScoreLess(new_u, anchor, best_new_u, best_anchor)) {
            return true;
        }
        if (isDirectedScoreLess(best_new_u, best_anchor, new_u, anchor)) {
            return false;
        }

        if (new_u != best_new_u) return new_u < best_new_u;
        if (anchor != best_anchor) return anchor < best_anchor;
        return edge_id < best_edge_id;
    }

    bool isDirectedScoreLess(ui lhs_new_u, ui lhs_anchor,
        ui rhs_new_u, ui rhs_anchor) const
    {
        int lhs_edge_id = q_edge_id[lhs_new_u][lhs_anchor];
        int rhs_edge_id = q_edge_id[rhs_new_u][rhs_anchor];
        assert(lhs_edge_id >= 0 && rhs_edge_id >= 0);

        const QEdge &lhs_edge = all_q_edges[(ui)lhs_edge_id];
        const QEdge &rhs_edge = all_q_edges[(ui)rhs_edge_id];
        const DirectedScore &lhs_score = lhs_new_u == lhs_edge.u
            ? lhs_edge.u_to_v : lhs_edge.v_to_u;
        const DirectedScore &rhs_score = rhs_new_u == rhs_edge.u
            ? rhs_edge.u_to_v : rhs_edge.v_to_u;

        __uint128_t lhs =
            (__uint128_t)lhs_score.numerator * rhs_score.denominator;
        __uint128_t rhs =
            (__uint128_t)rhs_score.numerator * lhs_score.denominator;
        if (lhs != rhs) return lhs < rhs;

        if (candidate_counts[lhs_new_u] != candidate_counts[rhs_new_u]) {
            return candidate_counts[lhs_new_u] < candidate_counts[rhs_new_u];
        }

        ui lhs_edge_freq = dataEdgeLabelFrequency(lhs_new_u, lhs_anchor);
        ui rhs_edge_freq = dataEdgeLabelFrequency(rhs_new_u, rhs_anchor);
        if (lhs_edge_freq != rhs_edge_freq) return lhs_edge_freq < rhs_edge_freq;
        return false;
    }

    ui chooseBestCrossingEdge(const vector<char> &visited,
        const vector<ui> &edge_ids, ui &best_new_u) const
    {
        ui best_edge_id = INVALID_EDGE;
        ui best_anchor = 0;
        best_new_u = 0;

        for (ui edge_id : edge_ids) {
            ui new_u = 0;
            ui anchor = 0;
            if (!edgeCrossesVisited(edge_id, visited, new_u, anchor)) continue;
            if (isBetterDirectedEdge(edge_id, new_u, anchor,
                best_edge_id, best_new_u, best_anchor)) {
                best_edge_id = edge_id;
                best_new_u = new_u;
                best_anchor = anchor;
            }
        }

        return best_edge_id;
    }

    ui chooseBestCrossingEdge(const vector<char> &visited, ui &best_new_u) const
    {
        ui best_edge_id = INVALID_EDGE;
        ui best_anchor = 0;
        best_new_u = 0;

        for (const QEdge &edge : all_q_edges) {
            ui new_u = 0;
            ui anchor = 0;
            if (!edgeCrossesVisited(edge.id, visited, new_u, anchor)) continue;
            if (isBetterDirectedEdge(edge.id, new_u, anchor,
                best_edge_id, best_new_u, best_anchor)) {
                best_edge_id = edge.id;
                best_new_u = new_u;
                best_anchor = anchor;
            }
        }

        return best_edge_id;
    }

    bool markPrefixVertices(const TreeState &state, ui edge_count,
        vector<char> &visited) const
    {
        visited.assign(qn, 0);
        visited[root] = 1;

        for (ui i = 0; i < edge_count; ++i) {
            const QEdge &edge = all_q_edges[state.edges[i]];
            bool u_vis = visited[edge.u] != 0;
            bool v_vis = visited[edge.v] != 0;
            if (u_vis == v_vis) return false;
            visited[u_vis ? edge.v : edge.u] = 1;
        }

        return true;
    }

    bool buildInitialTree(TreeState &state) const
    {
        state.edges.clear();
        state.R.clear();

        if (qn == 1) return true;
        assert(!all_q_edges.empty());

        vector<char> visited(qn, 0);
        visited[root] = 1;

        for (ui step = 0; step + 1 < qn; ++step) {
            ui next_vertex = 0;
            ui best = chooseBestCrossingEdge(visited, next_vertex);

            if (best == INVALID_EDGE) return false;

            visited[next_vertex] = 1;
            state.edges.push_back(best);
        }

        return state.edges.size() + 1 == qn;
    }

    bool replaceAndReorder(const TreeState &state, ui h, TreeState &next_state)
    {
        if (h >= state.edges.size()) return false;

        stats.replacement_calls++;

        ui removed_edge_id = state.edges[h];
        vector<char> in_tree(all_q_edges.size(), 0);
        for (ui edge_id : state.edges) in_tree[edge_id] = 1;

        DisjointSet dsu(qn);
        for (ui i = 0; i < state.edges.size(); ++i) {
            if (i == h) continue;
            const QEdge &edge = all_q_edges[state.edges[i]];
            dsu.unite(edge.u, edge.v);
        }

        ui replacement = INVALID_EDGE;
        ui replacement_new_u = 0;
        ui replacement_anchor = 0;
        int root_component = dsu.find((int)root);
        for (const QEdge &edge : all_q_edges) {
            if (in_tree[edge.id]) continue;
            if (state.R.find(edge.id) != state.R.end()) continue;

            int u_component = dsu.find((int)edge.u);
            int v_component = dsu.find((int)edge.v);
            if (u_component == v_component) continue;

            ui new_u = 0;
            ui anchor = 0;
            // The endpoint in the root component is reached first in the
            // reordered rooted tree, so it is the anchor of this direction.
            if (u_component == root_component) {
                anchor = edge.u;
                new_u = edge.v;
            }
            else {
                assert(v_component == root_component);
                anchor = edge.v;
                new_u = edge.u;
            }

            if (isBetterDirectedEdge(edge.id, new_u, anchor,
                replacement, replacement_new_u, replacement_anchor)) {
                replacement = edge.id;
                replacement_new_u = new_u;
                replacement_anchor = anchor;
            }
        }

        if (replacement == INVALID_EDGE) return false;

        next_state.R = state.R;
        next_state.R.insert(removed_edge_id);

        return reorderAfterReplacement(state, h, replacement, next_state);
    }

    bool reorderAfterReplacement(const TreeState &state, ui h, ui replacement, TreeState &next_state) const
    {
        vector<char> visited;
        next_state.edges.clear();

        if (!markPrefixVertices(state, h, visited)) return false;
        for (ui i = 0; i < h; ++i) next_state.edges.push_back(state.edges[i]);

        vector<ui> remaining;
        remaining.reserve(state.edges.size() - h);
        remaining.push_back(replacement);
        for (ui i = h + 1; i < state.edges.size(); ++i) {
            if (next_state.R.find(state.edges[i]) == next_state.R.end()) {
                remaining.push_back(state.edges[i]);
            }
        }

        while (!remaining.empty()) {
            ui best_pos = INVALID_EDGE;
            ui best_new_u = 0;
            ui best_edge_id = chooseBestCrossingEdge(visited, remaining, best_new_u);

            for (ui i = 0; i < remaining.size(); ++i) {
                if (remaining[i] == best_edge_id) {
                    best_pos = i;
                    break;
                }
            }

            if (best_edge_id == INVALID_EDGE) return false;

            visited[best_new_u] = 1;
            next_state.edges.push_back(best_edge_id);
            remaining.erase(remaining.begin() + best_pos);
        }

        return next_state.edges.size() + 1 == qn;
    }

    bool buildQISequence(const TreeState &state, QISequence &seq) const
    {
        seq.S.clear();
        seq.sEdgeIds.clear();
        seq.bEdgeIds.clear();
        seq.R = state.R;

        vector<char> visited(qn, 0);
        visited[root] = 1;

        seq.S.push_back(root);
        seq.sEdgeIds.push_back(INVALID_EDGE);
        seq.bEdgeIds.push_back(vector<ui>());

        for (ui edge_id : state.edges) {
            const QEdge &edge = all_q_edges[edge_id];
            bool u_vis = visited[edge.u] != 0;
            bool v_vis = visited[edge.v] != 0;
            if (u_vis == v_vis) return false;

            ui next_u = u_vis ? edge.v : edge.u;

            vector<ui> backward_edges;
            ui deg;
            const ui *nbrs = query_graph->getVertexNeighbors(next_u, deg);
            for (ui i = 0; i < deg; ++i) {
                ui nbr = nbrs[i];
                if (!visited[nbr]) continue;

                int back_id = q_edge_id[next_u][nbr];
                if (back_id >= 0 && (ui)back_id != edge_id) {
                    backward_edges.push_back((ui)back_id);
                }
            }

            visited[next_u] = 1;
            seq.S.push_back(next_u);
            seq.sEdgeIds.push_back(edge_id);
            seq.bEdgeIds.push_back(backward_edges);
        }

        return seq.S.size() == qn;
    }

    const QISequence *getQISequence(const TreeState &state)
    {
        Timer t_enum;
        t_enum.restart();

        string key = stateKey(state);
        auto it = sequence_cache.find(key);
        if (it != sequence_cache.end()) {
            stats.enum_time += t_enum.elapsed();
            return &it->second;
        }

        QISequence seq;
        if (!buildQISequence(state, seq)) {
            stats.enum_time += t_enum.elapsed();
            return nullptr;
        }

        auto inserted = sequence_cache.emplace(key, seq);
        stats.sequences_count = (long long)sequence_cache.size();
        stats.enum_time += t_enum.elapsed();
        return &inserted.first->second;
    }

    void SimSearchOnDemand(ui h, const TreeState &state, ui gamma)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }

        stats.recursion_calls++;
        stats.enum_call_count++;

        const QISequence *seq_ptr = getQISequence(state);
        if (seq_ptr == nullptr) return;
        const QISequence &seq = *seq_ptr;

        if (h == qn) {
            emitIfValid(seq);
            return;
        }

        ui u_curr = seq.S[h];
        ui tree_edge_id = seq.sEdgeIds[h];
        const QEdge &tree_edge = all_q_edges[tree_edge_id];
        ui u_parent = (tree_edge.u == u_curr) ? tree_edge.v : tree_edge.u;
        int v_parent_int = mapped_q[u_parent];
        if (v_parent_int < 0) return;

        ui v_parent = (ui)v_parent_int;
        const pair<size_t, ui> *tree_range =
            findAdjRange(u_parent, v_parent, u_curr);

        const ui *tree_begin = tree_range == nullptr
            ? nullptr : rangeBegin(*tree_range);
        const ui *tree_end = tree_range == nullptr
            ? nullptr : rangeEnd(*tree_range);

        for (const ui *it = tree_begin; it != tree_end; ++it) {
            ui v_curr = *it;
            if (isDataVertexUsed(v_curr)) continue;

            ui new_gamma = gamma;
            bool possible = true;

            for (ui back_edge_id : seq.bEdgeIds[h]) {
                const QEdge &back_edge = all_q_edges[back_edge_id];
                ui u_target = (back_edge.u == u_curr) ? back_edge.v : back_edge.u;
                int v_target_int = mapped_q[u_target];
                if (v_target_int < 0) {
                    possible = false;
                    break;
                }

                bool edge_exists = anchorAdjacent(
                    u_target, (ui)v_target_int, u_curr, v_curr);

                bool in_R = seq.R.find(back_edge_id) != seq.R.end();
                if (in_R) {
                    if (edge_exists) {
                        possible = false;
                        break;
                    }
                }
                else if (!edge_exists) {
                    new_gamma++;
                    if (new_gamma > threshold) {
                        possible = false;
                        break;
                    }
                }
            }

            if (!possible || new_gamma > threshold) continue;

            size_t branch_mark = mark();
            setMap(u_curr, (int)v_curr);
            pushUsed(v_curr);

            SimSearchOnDemand(h + 1, state, new_gamma);

            rollback(branch_mark);
            assert(mark() == branch_mark);

            if (outputLimitReached()) break;
        }

        if (!outputLimitReached() && gamma < threshold && h > 0) {
            TreeState next_state;
            ui current_tree_edge_pos = h - 1;
            if (replaceAndReorder(state, current_tree_edge_pos, next_state)) {
                SimSearchOnDemand(h, next_state, gamma + 1);
            }
        }
    }

    void emitIfValid(const QISequence &seq)
    {
        Timer t_verify;
        t_verify.restart();

        ui missing_edges = 0;

        for (const QEdge &edge : all_q_edges) {
            int vu_int = mapped_q[edge.u];
            int vv_int = mapped_q[edge.v];
            if (vu_int < 0 || vv_int < 0) {
                stats.verify_time += t_verify.elapsed();
                return;
            }

            bool exists = anchorAdjacent(
                edge.u, (ui)vu_int, edge.v, (ui)vv_int);

            if (seq.R.find(edge.id) != seq.R.end() && exists) {
                stats.verify_time += t_verify.elapsed();
                return;
            }
            if (!exists) missing_edges++;
            if (missing_edges > threshold) {
                stats.verify_time += t_verify.elapsed();
                return;
            }
        }

        string key = mappingKey();
        if (!result_keys.insert(key).second) {
            stats.duplicate_results++;
            stats.verify_time += t_verify.elapsed();
            return;
        }

#ifndef NDEBUG
        vector<pair<ui, ui>> res;
        res.reserve(qn);
        for (ui u = 0; u < qn; ++u) {
            res.push_back({ u, (ui)mapped_q[u] });
        }
        results_ptr->push_back(res);
#endif
        stats.result_count++;
        noteOutputLimitIfReached();

        stats.verify_time += t_verify.elapsed();
    }

    string stateKey(const TreeState &state) const
    {
        string key;
        key.reserve(state.edges.size() * 6 + state.R.size() * 6 + 4);
        key += "T:";
        for (ui edge_id : state.edges) {
            key += to_string(edge_id);
            key += ',';
        }
        key += "|R:";
        for (ui edge_id : state.R) {
            key += to_string(edge_id);
            key += ',';
        }
        return key;
    }

    string mappingKey() const
    {
        string key;
        key.reserve(qn * 8);
        for (ui u = 0; u < qn; ++u) {
            key += to_string(mapped_q[u]);
            key += ',';
        }
        return key;
    }
};

const ui TreeSpanSolver::Impl::INVALID_EDGE;

TreeSpanSolver::TreeSpanSolver() : impl_(new Impl(stats)) {}

TreeSpanSolver::~TreeSpanSolver() = default;

bool TreeSpanSolver::init(const Graph *query_graph, const Graph *data_graph,
    ui threshold)
{
    return impl_->init(query_graph, data_graph, threshold);
}

void TreeSpanSolver::match(ssm_ged::MatchResults &results)
{
    impl_->match(results);
}

void TreeSpanSolver::printStats() const
{
    impl_->printStats();
}

void Approximate_TreeSpan(const Graph *query_graph, const Graph *data_graph,
    ssm_ged::MatchResults &results, ui threshold)
{
    Timer t_total;
    t_total.restart();

    TreeSpanSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(results);
    }

    solver.stats.total_time = t_total.elapsed();
    ssm_ged::set_reported_phase_times(solver.stats.init_time,
        solver.stats.search_time);
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

namespace ssm_ged {

namespace {

void run_treespan(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold)
{
    ::Approximate_TreeSpan(query_graph, data_graph, results, threshold);
}

} // namespace

const AlgorithmDefinition &create_algorithm_definition()
{
    static const AlgorithmDefinition definition = {
        "treespan",
        "TreeSpan",
        &run_treespan
    };
    return definition;
}

} // namespace ssm_ged
