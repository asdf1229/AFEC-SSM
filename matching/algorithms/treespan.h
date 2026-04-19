#ifndef MATCHING_ALGORITHMS_TREESPAN_H_
#define MATCHING_ALGORITHMS_TREESPAN_H_

#include "graph/graph.h"
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
            mapped_q[start_u] = (int)v;
            mapped_g[v] = (int)start_u;

            SimSearchOnDemand(1, initial_tree, 0);

            mapped_q[start_u] = -1;
            mapped_g[v] = -1;
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
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        printf("\n--- TreeSpan Matching Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Q_T OD Time:     %.4lf ms\n", stats.enum_time / 1000.0);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
        printf("  - Verify Time:     %.4lf ms (%.2f%% of Search)\n", stats.verify_time / 1000.0, pct(stats.verify_time, stats.search_time));
        printf("QISequences Cached:  %lld\n", stats.sequences_count);
        printf("OD Calls:            %lld\n", stats.enum_call_count);
        printf("Replacement Calls:   %lld\n", stats.replacement_calls);
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("[hasEdge] Calls:     %lld\n", stats.has_edge_calls);
        printf("Duplicate Results:   %lld\n", stats.duplicate_results);
        printf("Results Found:       %zu\n", results_ptr ? results_ptr->size() : 0);
        printf("---------------------------------------\n");
    }

private:
    static const ui INVALID_EDGE = (ui)-1;

    struct QEdge {
        ui u = 0;
        ui v = 0;
        ui weight = 0;
        ui id = 0;

        bool operator<(const QEdge &other) const
        {
            if (weight != other.weight) return weight < other.weight;
            if (u != other.u) return u < other.u;
            return v < other.v;
        }
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
    vector<vector<ui>> Lq_counts;
    vector<vector<ui>> Lg_counts;
    vector<ui> Lq_degrees;
    vector<ui> Lg_degrees;

    vector<vector<char>> q_matrix;
    vector<vector<int>> q_edge_id;
    vector<QEdge> all_q_edges;

    vector<int> mapped_q;
    vector<int> mapped_g;

    TreeState initial_tree;
    map<string, QISequence> sequence_cache;
    unordered_set<string> result_keys;

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
        Lq_counts.clear();
        Lg_counts.clear();
        Lq_degrees.clear();
        Lg_degrees.clear();
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

    void orderQueryEdges()
    {
        for (QEdge &e : all_q_edges) {
            ui cu = (ui)candidates[e.u].size();
            ui cv = (ui)candidates[e.v].size();
            ui du = query_graph->getVertexDegree(e.u);
            ui dv = query_graph->getVertexDegree(e.v);
            e.weight = cu + cv + du + dv;
        }

        sort(all_q_edges.begin(), all_q_edges.end());

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

    ui selectRoot()
    {
        ui best_u = 0;
        int best_cands = candidates[0].size();
        ui best_degree = query_graph->getVertexDegree(0);

        for (ui u = 1; u < qn; ++u) {
            int cand_size = candidates[u].size();
            ui degree = query_graph->getVertexDegree(u);
            if (cand_size < best_cands || (cand_size == best_cands && degree > best_degree)) {
                best_u = u;
                best_cands = cand_size;
                best_degree = degree;
            }
        }
        return best_u;
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
            ui best = INVALID_EDGE;

            for (const QEdge &e : all_q_edges) {
                bool u_vis = visited[e.u] != 0;
                bool v_vis = visited[e.v] != 0;
                if (u_vis != v_vis) {
                    best = e.id;
                    break;
                }
            }

            if (best == INVALID_EDGE) return false;

            const QEdge &edge = all_q_edges[best];
            ui next_vertex = visited[edge.u] ? edge.v : edge.u;
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
        vector<char> visited(qn, 0);
        visited[root] = 1;
        next_state.edges.clear();

        for (ui i = 0; i < h; ++i) {
            ui edge_id = state.edges[i];
            const QEdge &edge = all_q_edges[edge_id];
            bool u_vis = visited[edge.u] != 0;
            bool v_vis = visited[edge.v] != 0;
            if (u_vis == v_vis) return false;

            visited[u_vis ? edge.v : edge.u] = 1;
            next_state.edges.push_back(edge_id);
        }

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
            ui best_edge_id = INVALID_EDGE;

            for (ui i = 0; i < remaining.size(); ++i) {
                const QEdge &edge = all_q_edges[remaining[i]];
                bool u_vis = visited[edge.u] != 0;
                bool v_vis = visited[edge.v] != 0;
                if (u_vis == v_vis) continue;
                if (best_edge_id == INVALID_EDGE || remaining[i] < best_edge_id) {
                    best_edge_id = remaining[i];
                    best_pos = i;
                }
            }

            if (best_edge_id == INVALID_EDGE) return false;

            const QEdge &edge = all_q_edges[best_edge_id];
            visited[visited[edge.u] ? edge.v : edge.u] = 1;
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
        }

        if (gamma < threshold && h > 0) {
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

        vector<pair<ui, ui>> res;
        res.reserve(qn);
        for (ui u = 0; u < qn; ++u) {
            res.push_back({ u, (ui)mapped_q[u] });
        }
        results_ptr->push_back(res);

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
}

#endif
