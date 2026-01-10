Here is the updated code for the `TreeSpanSolver` class and the `Approximate_TreeSpan` function. I have added a `TreeStats` structure, a `printStats` function, and inserted `Timer` calls to measure initialization, search, verification, and reordering times.

```cpp
class TreeSpanSolver {
public:
    TreeSpanSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

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
        ui root = selectRoot();
        QISequence initial_seq = generateInitialSequence(root);

        // 4. 开始 SimSearch
        ui start_u = initial_seq.S[0];
        
        for (ui v : candidates[start_u]) {
            mapped_q[start_u] = v;
            mapped_g[v] = start_u;
            
            SimSearch(1, initial_seq, 0);

            mapped_q[start_u] = -1;
            mapped_g[v] = -1;
        }

        stats.search_time = t_search.elapsed();
    }

private:
    // ... (Data Structures struct QEdge, QISequence remain unchanged) ...
    struct QEdge {
        ui u, v;
        bool operator==(const QEdge& other) const {
            return (u == other.u && v == other.v) || (u == other.v && v == other.u);
        }
        QEdge canonical() const {
            return (u < v) ? *this : QEdge{v, u};
        }
        bool operator<(const QEdge& other) const {
            if (u != other.u) return u < other.u;
            return v < other.v;
        }
    };

    struct QISequence {
        vector<ui> S;             
        vector<QEdge> sEdge;      
        vector<vector<QEdge>> bEdges; 
        set<QEdge> R;             
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

    // --- Core Algorithm: SimSearch ---

    void SimSearch(ui h, QISequence seq, ui gamma) {
        stats.recursion_calls++;

        if (h == qn) {
            vector<pair<ui, ui>> res;
            res.reserve(qn);
            for(ui u=0; u<qn; ++u) res.push_back({u, (ui)mapped_q[u]});
            results_ptr->push_back(res);
            return;
        }

        // --- Phase 1: Go-Down (纵向扩展) ---
        ui u_curr = seq.S[h];
        QEdge tree_edge = seq.sEdge[h];
        ui u_parent = (tree_edge.u == u_curr) ? tree_edge.v : tree_edge.u;
        ui v_parent = mapped_q[u_parent];

        ui deg_g; const ui* neighbors_g = data_graph->getVertexNeighbors(v_parent, deg_g);
        
        for (ui k = 0; k < deg_g; ++k) {
            ui v_curr = neighbors_g[k];
            
            if (mapped_g[v_curr] != -1) continue; 
            if (!binary_search(candidates[u_curr].begin(), candidates[u_curr].end(), v_curr)) continue;

            Timer t_check; // Measure verification overhead
            t_check.restart();

            ui new_gamma = gamma;
            bool possible = true;

            for (const auto& bedge : seq.bEdges[h]) {
                ui u_target = (bedge.u == u_curr) ? bedge.v : bedge.u;
                ui v_target = mapped_q[u_target]; 
                
                bool edge_exists = data_graph->hasEdge(v_curr, v_target);
                bool in_R = seq.R.count(bedge.canonical());
                
                if (in_R && edge_exists) {
                    possible = false; 
                    break;
                }
                
                if (!edge_exists) {
                    if (!in_R) new_gamma++;
                }
            }
            stats.verify_time += t_check.elapsed();

            if (possible && new_gamma <= threshold) {
                mapped_q[u_curr] = v_curr;
                mapped_g[v_curr] = u_curr;
                
                SimSearch(h + 1, seq, new_gamma);
                
                mapped_q[u_curr] = -1;
                mapped_g[v_curr] = -1;
            }
        }

        // --- Phase 2: Alternating-Reordering (横向重组) ---
        if (gamma < threshold) {
            Timer t_reorder; // Measure Reordering Overhead
            t_reorder.restart();

            QEdge current_tree_edge = seq.sEdge[h];
            QEdge new_edge;
            
            if (findReplacement(seq, h, new_edge)) {
                QISequence next_seq;
                next_seq.R = seq.R;
                next_seq.R.insert(current_tree_edge.canonical());
                
                if (reorderSequence(seq, h, new_edge, next_seq)) {
                    stats.reorder_time += t_reorder.elapsed(); // Stop timer before recursion
                    
                    SimSearch(h, next_seq, gamma + 1);
                    
                    // Resume timer isn't strictly necessary as we are back in this scope,
                    // but usually we count logic time, not recursion wait time for this block.
                } else {
                    stats.reorder_time += t_reorder.elapsed();
                }
            } else {
                stats.reorder_time += t_reorder.elapsed();
            }
        }
    }

    // --- Helper Functions ---

    bool findReplacement(const QISequence& seq, ui h, QEdge& out_edge) {
        vector<bool> visited(qn, false);
        for(ui i=0; i<h; ++i) visited[seq.S[i]] = true;

        QEdge removed = seq.sEdge[h].canonical();
        bool found = false;
        QEdge best_e;

        for (const auto& e : all_q_edges) {
            bool u_vis = visited[e.u];
            bool v_vis = visited[e.v];
            
            if (u_vis != v_vis) {
                QEdge cand = e.canonical();
                if (cand == removed) continue;
                if (seq.R.count(cand)) continue;

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

    bool reorderSequence(const QISequence& old_seq, ui h, QEdge new_edge, QISequence& new_seq) {
        stats.reorder_calls++; // Count operations

        new_seq.S.resize(h);
        new_seq.sEdge.resize(h);
        new_seq.bEdges.resize(h);
        
        vector<bool> visited(qn, false);
        for(ui i=0; i<h; ++i) {
            new_seq.S[i] = old_seq.S[i];
            new_seq.sEdge[i] = old_seq.sEdge[i];
            new_seq.bEdges[i] = old_seq.bEdges[i];
            visited[old_seq.S[i]] = true;
        }

        ui u_next = (visited[new_edge.u]) ? new_edge.v : new_edge.u;
        if (visited[u_next]) return false;

        new_seq.S.push_back(u_next);
        new_seq.sEdge.push_back(new_edge);
        visited[u_next] = true;
        
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

        for (ui step = h + 1; step < qn; ++step) {
            QEdge best_prim_edge;
            bool found = false;
            
            for (const auto& e : all_q_edges) {
                if (new_seq.R.count(e.canonical())) continue;

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

            if (!found) return false;

            ui v_new = visited[best_prim_edge.u] ? best_prim_edge.v : best_prim_edge.u;
            new_seq.S.push_back(v_new);
            new_seq.sEdge.push_back(best_prim_edge);
            visited[v_new] = true;

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

    QISequence generateInitialSequence(ui root) {
        QISequence seq;
        seq.S.push_back(root);
        seq.sEdge.push_back({(ui)-1, (ui)-1});
        seq.bEdges.push_back({});
        
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

void Approximate_TreeSpan(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
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
```