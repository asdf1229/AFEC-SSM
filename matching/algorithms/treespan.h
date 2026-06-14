#ifndef MATCHING_ALGORITHMS_TREESPAN_H_
#define MATCHING_ALGORITHMS_TREESPAN_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"

using namespace std;

class TreeSpanSolver {
public:
    TreeSpanSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

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
        mapped_g.assign(gn, -1);
        candidates.assign(qn, MyBitset(gn));

        buildQueryMatrix();
        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

        if (!calVerticesFilter()) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        cacheCandidateCounts();
        initOrderingStatistics();
        orderQueryEdges();

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

        const QISequence *seq = getQISequence(initial_tree);
        if (seq == nullptr) {
            stats.search_time = t_search.elapsed();
            return;
        }

        ui start_u = seq->S[0];
        for (ui v : candidates[start_u]) {
            if (outputLimitReached()) break;

            mapped_q[start_u] = (int)v;
            mapped_g[v] = (int)start_u;

            SimSearchOnDemand(1, initial_tree, 0);

            mapped_q[start_u] = -1;
            mapped_g[v] = -1;

            if (outputLimitReached()) break;
        }

        stats.search_time = t_search.elapsed();
        stats.sequences_count = (long long)sequence_cache.size();
    }

    struct TreeStats {
        long long total_time = 0;
        long long init_time = 0;
        long long search_time = 0;
        long long enum_time = 0;
        long long verify_time = 0;
        long long recursion_calls = 0;
        long long has_edge_calls = 0;
        long long replacement_calls = 0;
        long long sequences_count = 0;
        long long enum_call_count = 0;
        long long duplicate_results = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
    } stats;

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
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
        printf("  - Q_T OD Time:     %.4lf ms (%.2f%% of Search)\n", stats.enum_time / 1000.0, pct(stats.enum_time, stats.search_time));
        printf("  - Verify Time:     %.4lf ms (%.2f%% of Search)\n", stats.verify_time / 1000.0, pct(stats.verify_time, stats.search_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.search_time));
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Replacement Calls:   %lld\n", stats.replacement_calls);
        printf("QISequences Cached:  %lld\n", stats.sequences_count);
        printf("OD Calls:            %lld\n", stats.enum_call_count);
        printf("[hasEdge] Calls:     %lld\n", stats.has_edge_calls);
        printf("Duplicate Results:   %lld\n", stats.duplicate_results);
        printf("Results Found:       %zu\n", stats.result_count);
#if MATCH_OUTPUT_LIMIT > 0
        printf("Output Limit:        %zu%s\n",
            (size_t)MATCH_OUTPUT_LIMIT,
            stats.output_limit_reached ? " (reached)" : "");
#endif
        printf("-----------------------------------------------------------\n");
    }

private:
    static const ui INVALID_EDGE = (ui)-1;

    struct QEdge {
        ui u = 0;
        ui v = 0;
        ui id = 0;
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
    vector<vector<ui>> Lq_counts;
    vector<vector<ui>> Lg_counts;
    vector<ui> Lq_degrees;
    vector<ui> Lg_degrees;
    vector<ui> data_label_frequency;
    vector<vector<ui>> data_edge_label_frequency;

    vector<vector<char>> q_matrix;
    vector<vector<int>> q_edge_id;
    vector<QEdge> all_q_edges;

    vector<int> mapped_q;
    vector<int> mapped_g;

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
        Lq_degrees.clear();
        Lg_degrees.clear();
        data_label_frequency.clear();
        data_edge_label_frequency.clear();
        q_matrix.clear();
        q_edge_id.clear();
        all_q_edges.clear();
        mapped_q.clear();
        mapped_g.clear();
        initial_tree = TreeState();
        sequence_cache.clear();
        result_keys.clear();
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
                    all_q_edges.push_back(e);
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

    void orderQueryEdges()
    {
        sort(all_q_edges.begin(), all_q_edges.end(),
            [this](const QEdge &lhs, const QEdge &rhs) {
                ui lhs_new = lhs.u;
                ui lhs_anchor = lhs.v;
                if (isDirectedScoreLess(lhs.v, lhs.u, lhs.u, lhs.v)) {
                    lhs_new = lhs.v;
                    lhs_anchor = lhs.u;
                }

                ui rhs_new = rhs.u;
                ui rhs_anchor = rhs.v;
                if (isDirectedScoreLess(rhs.v, rhs.u, rhs.u, rhs.v)) {
                    rhs_new = rhs.v;
                    rhs_anchor = rhs.u;
                }

                if (isDirectedScoreLess(lhs_new, lhs_anchor, rhs_new, rhs_anchor)) {
                    return true;
                }
                if (isDirectedScoreLess(rhs_new, rhs_anchor, lhs_new, lhs_anchor)) {
                    return false;
                }
                if (lhs.u != rhs.u) return lhs.u < rhs.u;
                return lhs.v < rhs.v;
            });

        q_edge_id.assign(qn, vector<int>(qn, -1));
        for (ui i = 0; i < all_q_edges.size(); ++i) {
            all_q_edges[i].id = i;
            ui u = all_q_edges[i].u;
            ui v = all_q_edges[i].v;
            q_edge_id[u][v] = (int)i;
            q_edge_id[v][u] = (int)i;
        }
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees)
    {
        ui n = g->getVerticesCount();
        counts.assign(n, vector<ui>(label_count, 0));
        degrees.assign(n, 0);

        for (ui u = 0; u < n; ++u) {
            ui deg;
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

    ui computeDelta(ui u, ui v) const
    {
        ui diff = 0;
        for (ui i = 0; i < label_count; ++i) {
            if (Lq_counts[u][i] > Lg_counts[v][i]) {
                diff += Lq_counts[u][i] - Lg_counts[v][i];
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
        unsigned long long lhs_num =
            (unsigned long long)candidate_counts[lhs_new_u] *
            (unsigned long long)dataEdgeLabelFrequency(lhs_new_u, lhs_anchor);
        unsigned long long rhs_num =
            (unsigned long long)candidate_counts[rhs_new_u] *
            (unsigned long long)dataEdgeLabelFrequency(rhs_new_u, rhs_anchor);
        unsigned long long lhs_den = (unsigned long long)dataLabelFrequency(lhs_new_u);
        unsigned long long rhs_den = (unsigned long long)dataLabelFrequency(rhs_new_u);

        __uint128_t lhs = (__uint128_t)lhs_num * rhs_den;
        __uint128_t rhs = (__uint128_t)rhs_num * lhs_den;
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
        for (const QEdge &edge : all_q_edges) {
            if (in_tree[edge.id]) continue;
            if (state.R.find(edge.id) != state.R.end()) continue;
            if (dsu.find((int)edge.u) != dsu.find((int)edge.v)) {
                replacement = edge.id;
                break;
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
        ui deg_g;
        const ui *neighbors_g = data_graph->getVertexNeighbors(v_parent, deg_g);

        for (ui i = 0; i < deg_g; ++i) {
            ui v_curr = neighbors_g[i];
            if (mapped_g[v_curr] != -1) continue;
            if (!candidates[u_curr].contains(v_curr)) continue;

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

                bool edge_exists = data_graph->hasEdge(v_curr, (ui)v_target_int);
                stats.has_edge_calls++;

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

            mapped_q[u_curr] = (int)v_curr;
            mapped_g[v_curr] = (int)u_curr;

            SimSearchOnDemand(h + 1, state, new_gamma);

            mapped_q[u_curr] = -1;
            mapped_g[v_curr] = -1;

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

            bool exists = data_graph->hasEdge((ui)vu_int, (ui)vv_int);
            stats.has_edge_calls++;

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

void Approximate_TreeSpan(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui>>> &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    TreeSpanSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
