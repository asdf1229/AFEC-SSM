#ifndef MATCHING_ALGORITHMS_CDE_EDGE_IE_H_
#define MATCHING_ALGORITHMS_CDE_EDGE_IE_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
using namespace std;

// ============================================================================
// MatchingSolver Implementation
// ============================================================================
class MatchingSolver {
public:
    MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

    bool init(const Graph *q, const Graph *g, ui match_threshold)
    {
        Timer t_init;
        t_init.restart();

        query_graph = q;
        data_graph = g;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();
        if(qn == 0 || gn == 0) return false;
        max_g_deg = data_graph->getMaxDegree();
        label_count = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());

        resetState();

        // init q_matrix and q_neighbors
        q_matrix.assign(qn, vector<char>(qn, 0));
        q_neighbors.assign(qn, vector<ui>());
        q_degree.assign(qn, 0);
        for (ui u = 0; u < qn; ++u) {
            ui deg = 0; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) {
                ui u1 = nbrs[i];
                q_matrix[u][u1] = 1;
                q_neighbors[u].push_back(u1);
            }
            q_degree[u] = (ui)q_neighbors[u].size();
        }

        Timer t_filter;
        bool res = runCandidateFiltering();
        stats.filter_time = t_filter.elapsed();
        if (!res) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        Timer t_label_degree;
        initDataLabelDegreeIndex();
        stats.init_label_degree_time = t_label_degree.elapsed();
#if CDE_EDGE_IE_FIXED_ORDER
        initFixedEdgePriorities();
#endif
        stats.init_time = t_init.elapsed();
        return true;
    }

    void match(vector<vector<pair<ui, ui>>> &results)
    {
        Timer t_search;
        t_search.restart();

        results_ptr = &results;
        results_ptr->clear();
        initDfsBuffer();

        BranchSelector branch_selector(*this);
        ui root = branch_selector.selectInitialRoot();

        for (ui v0 : candidates[root]) {
            mapped_q[root] = (int)v0;
            mapped_g[v0] = (int)root;
            part_M.push_back({ root, v0 });

            updateFrontier(root, true);

            dfs(0);

            part_M.pop_back();

            mapped_g[v0] = -1;
            mapped_q[root] = -1;
            updateFrontier(root, false);

            if (outputLimitReached()) break;
        }

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct TimeStats {
        long long total_time = 0;
        // init breakdown
        long long init_time = 0;
        long long filter_time = 0;
        long long filter_nlf_time = 0;
        long long filter_bridge_time = 0;
        long long filter_spoke_time = 0;
        long long init_label_degree_time = 0;
        ui filter_candidate_count = 0;
        // search breakdown
        long long dfs_time = 0;
        long long dfs_buffer_clear_time = 0;
        long long terminal_time = 0;
        long long terminal_scan_time = 0;
        long long terminal_tail_build_time = 0;
        long long terminal_bucket_time = 0;
        long long terminal_candidate_collect_time = 0;
        long long terminal_min_delta_time = 0;
        long long terminal_tail_enum_time = 0;
        long long frontier_time = 0; // building and ordering U_frontier
        long long frontier_select_time = 0;
        long long frontier_component_time = 0;
        long long frontier_sort_time = 0;
        long long frontier_score_time = 0;
        long long frontier_score_live_anchor_time = 0;
        long long frontier_score_live_candidate_time = 0;
        long long frontier_score_query_degree_time = 0;
        long long frontier_score_cheap_support_time = 0;
        long long branch_time = 0;   // candidate enumeration & matching in dfs
        long long branch_cal_edge_support_time = 0;
        long long branch_count_anchors_time = 0;
        long long branch_delta_bucket_time = 0;
        long long support_update_time = 0;
        long long candidate_loop_time = 0;
        long long exclude_update_time = 0;
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
        unsigned long long terminal_scan_calls = 0;
        unsigned long long terminal_tail_build_calls = 0;
        unsigned long long terminal_tail_calls = 0;
        unsigned long long terminal_prune_calls = 0;
        unsigned long long terminal_delayed_vertices = 0;
        unsigned long long terminal_bucket_candidate_checks = 0;
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        long long search_accounted_time = stats.dfs_buffer_clear_time +
            stats.terminal_time + stats.frontier_time + stats.branch_time;
        long long search_other_time = stats.dfs_time > search_accounted_time
            ? stats.dfs_time - search_accounted_time : 0;

        printf("\n--- CDE-Edge-IE Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("    - NLF:           %.4lf ms (%.2f%% of Filter)\n", stats.filter_nlf_time / 1000.0, pct(stats.filter_nlf_time, stats.filter_time));
        printf("    - Bridge/Index:  %.4lf ms (%.2f%% of Filter)\n", stats.filter_bridge_time / 1000.0, pct(stats.filter_bridge_time, stats.filter_time));
#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
        printf("    - Spoke:         %.4lf ms (%.2f%% of Filter)\n", stats.filter_spoke_time / 1000.0, pct(stats.filter_spoke_time, stats.filter_time));
#endif
        printf("  - Label Deg Index: %.4lf ms (%.2f%% of Init)\n", stats.init_label_degree_time / 1000.0, pct(stats.init_label_degree_time, stats.init_time));
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - DFS Buffer Clear:%.4lf ms (%.2f%% of Search)\n", stats.dfs_buffer_clear_time / 1000.0, pct(stats.dfs_buffer_clear_time, stats.dfs_time));
        printf("  - Terminal Time:   %.4lf ms (%.2f%% of Search)\n", stats.terminal_time / 1000.0, pct(stats.terminal_time, stats.dfs_time));
        printf("    - Scan:          %.4lf ms (%.2f%% of Terminal)\n", stats.terminal_scan_time / 1000.0, pct(stats.terminal_scan_time, stats.terminal_time));
        printf("    - Tail Build:    %.4lf ms (%.2f%% of Terminal)\n", stats.terminal_tail_build_time / 1000.0, pct(stats.terminal_tail_build_time, stats.terminal_time));
        printf("    - Cand Buckets:  %.4lf ms (%.2f%% of Terminal)\n", stats.terminal_bucket_time / 1000.0, pct(stats.terminal_bucket_time, stats.terminal_time));
        printf("      - Collect Supp:%.4lf ms (%.2f%% of Cand Buckets)\n", stats.terminal_candidate_collect_time / 1000.0, pct(stats.terminal_candidate_collect_time, stats.terminal_bucket_time));
        printf("    - Min Delta:     %.4lf ms (%.2f%% of Terminal)\n", stats.terminal_min_delta_time / 1000.0, pct(stats.terminal_min_delta_time, stats.terminal_time));
        printf("    - Tail Enum:     %.4lf ms (%.2f%% of Terminal)\n", stats.terminal_tail_enum_time / 1000.0, pct(stats.terminal_tail_enum_time, stats.terminal_time));
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("    - Select Best:   %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_select_time / 1000.0, pct(stats.frontier_select_time, stats.frontier_time));
        printf("    - Component:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_component_time / 1000.0, pct(stats.frontier_component_time, stats.frontier_time));
#ifdef ENABLE_FRONTIER_ORDERING
        printf("    - Sort Hook:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_sort_time / 1000.0, pct(stats.frontier_sort_time, stats.frontier_time));
#endif
        printf("    - Score Total:   %.4lf ms (%.2f%% of Select+Sort)\n", stats.frontier_score_time / 1000.0, pct(stats.frontier_score_time, stats.frontier_select_time + stats.frontier_sort_time));
        printf("      - Live Anchors:   %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_live_anchor_time / 1000.0, pct(stats.frontier_score_live_anchor_time, stats.frontier_score_time));
        printf("      - Live Candidates:%.4lf ms (%.2f%% of Score)\n", stats.frontier_score_live_candidate_time / 1000.0, pct(stats.frontier_score_live_candidate_time, stats.frontier_score_time));
        printf("      - Query Degree:   %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_query_degree_time / 1000.0, pct(stats.frontier_score_query_degree_time, stats.frontier_score_time));
        printf("      - Cheap Support:  %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_cheap_support_time / 1000.0, pct(stats.frontier_score_cheap_support_time, stats.frontier_score_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("    - branch_cal_edge_support_ms: %.4lf (%.2f%% of Branch)\n", stats.branch_cal_edge_support_time / 1000.0, pct(stats.branch_cal_edge_support_time, stats.branch_time));
        printf("    - branch_count_anchors_ms:    %.4lf (%.2f%% of Branch)\n", stats.branch_count_anchors_time / 1000.0, pct(stats.branch_count_anchors_time, stats.branch_time));
        printf("    - branch_delta_bucket_ms:     %.4lf (%.2f%% of Branch)\n", stats.branch_delta_bucket_time / 1000.0, pct(stats.branch_delta_bucket_time, stats.branch_time));
        printf("    - support_update_ms:          %.4lf (%.2f%% of Branch)\n", stats.support_update_time / 1000.0, pct(stats.support_update_time, stats.branch_time));
        printf("    - candidate_loop_ms:          %.4lf (%.2f%% of Branch)\n", stats.candidate_loop_time / 1000.0, pct(stats.candidate_loop_time, stats.branch_time));
        printf("    - exclude_update_ms:          %.4lf (%.2f%% of Branch)\n", stats.exclude_update_time / 1000.0, pct(stats.exclude_update_time, stats.branch_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.dfs_time));
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
        printf("Terminal Scan Calls: %llu\n", stats.terminal_scan_calls);
        printf("Terminal Build Calls:%llu\n", stats.terminal_tail_build_calls);
        printf("Terminal Tail Calls: %llu\n", stats.terminal_tail_calls);
        printf("Terminal Prunes:     %llu\n", stats.terminal_prune_calls);
        printf("Terminal Delayed:    %llu\n", stats.terminal_delayed_vertices);
        printf("Terminal Cand Checks:%llu\n", stats.terminal_bucket_candidate_checks);
        printf("Results Found:       %zu\n", stats.result_count);
#if MATCH_OUTPUT_LIMIT > 0
        printf("Output Limit:        %zu%s\n",
            (size_t)MATCH_OUTPUT_LIMIT,
            stats.output_limit_reached ? " (reached)" : "");
#endif
        printf("-----------------------------------------------------------\n");
        fflush(stdout);
    }

private:
    struct DataLabelDegreeCount {
        LabelID label = 0;
        ui count = 0;
    };

    struct ActiveEdge {
        ui u = 0;       // unmatched endpoint
        ui anchor = 0;  // matched endpoint
        ui anchor_support = std::numeric_limits<ui>::max();
        ui live_anchor_count = 0;
        ui query_degree = 0;
    };

    struct EdgeScoreCache {
        vector<vector<ActiveEdge>> active_edges_by_vertex;
        vector<char> active_edges_cached;
    };

    struct SupportSnapshot {
        ui u = 0;
        ui anchor = 0;
        ui value = 0;
        char dirty = 0;
    };

    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;
    ui threshold;
    ui qn, gn;
    ui label_count;
    ui max_g_deg;
    vector<vector<char>> q_matrix;
    vector<vector<ui>> q_neighbors;
    vector<ui> q_degree;
    vector<vector<char>> q_neighbor_is_bridge;
#if CDE_EDGE_IE_FIXED_ORDER
    vector<vector<ui>> fixed_edge_priority;
#endif
    vector<vector<DataLabelDegreeCount>> data_label_degrees;

    vector<MyBitset> candidates;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<vector<char>> excluded_edges;
    vector<MyBitset> excluded_cands;
    vector<pair<ui, ui>> part_M;

    vector<ui> anchor_count;
    vector<vector<ui>> support;
    vector<vector<char>> support_dirty;
    vector<SupportSnapshot> *active_support_snapshots = nullptr;
    vector<int> frontier_pos;
    vector<ui> active_frontier;
    vector<ui> data_vertex_mark;
    vector<ui> data_vertex_mark_pos;
    ui data_vertex_mark_token = 0;

    struct TerminalScan {
        ui unmatched_count = 0;
        ui terminal_count = 0;
        ui terminal_frontier_count = 0;
        ui nonterminal_frontier_count = 0;

        bool allRemainingTerminal() const
        {
            return unmatched_count > 0 && unmatched_count == terminal_count;
        }

        bool hasNonterminalFrontier() const
        {
            return nonterminal_frontier_count > 0;
        }
    };

    struct TerminalTailVertex {
        ui u = 0;
        ui feasible_count = 0;
        ui min_delta = std::numeric_limits<ui>::max();
        vector<vector<ui>> buckets;
    };

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

    void resetState()
    {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        excluded_edges.assign(qn, vector<char>(qn, 0));
        part_M.clear();
        part_M.reserve(qn);

        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        excluded_cands.clear();
        excluded_cands.assign(qn, MyBitset(gn));

        q_matrix.clear();
        q_neighbors.clear();
        q_degree.clear();
        q_neighbor_is_bridge.clear();
#if CDE_EDGE_IE_FIXED_ORDER
        fixed_edge_priority.clear();
#endif
        data_label_degrees.clear();

        anchor_count.assign(qn, 0);
        support.assign(qn, vector<ui>(qn, 0));
        support_dirty.assign(qn, vector<char>(qn, 1));
        active_support_snapshots = nullptr;
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
        data_vertex_mark.assign(gn, 0);
        data_vertex_mark_pos.assign(gn, 0);
        data_vertex_mark_token = 0;
        stats = TimeStats();
    }

    // ========================================================================
    // Filtering
    // ========================================================================
    struct CandidateFilter {
    private:
        MatchingSolver &solver;

        struct LabelFreq {
            LabelID label = 0;
            ui count = 0;
        };

        struct QueryLabelReq {
            LabelID label = 0;
            ui bridge = 0;
            ui normal = 0;
        };

        struct BridgeArc {
            ui from = 0;
            ui to = 0;
        };

        struct BridgeNbr {
            ui to = 0;
            ui support_arc = 0;
        };

        // static graph/profile cache
        vector<LabelID> query_label;
        vector<LabelID> data_label;
        vector<ui> data_degree;

        // NLF filtering
        vector<vector<LabelFreq>> data_label_freqs;
        vector<vector<QueryLabelReq>> query_label_reqs;
        vector<vector<ui>> data_by_label;
        vector<ui> query_bridge_degree;

        // Bridge filtering
        vector<BridgeArc> bridge_arcs;
        vector<vector<BridgeNbr>> bridge_nbrs;
        vector<ui> bridge_support;
        queue<pair<ui, ui>> removed;

        // spoke filtering
        vector<vector<ui>>  spoke_adj;
        vector<int>         match_right;
        vector<ui>          seen_right;
        ui                  seen_token = 1;
        vector<char>        left_is_bridge;
#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
        queue<ui>           pending_spokes;
        vector<char>        queued_spoke;
#endif

    public:
        explicit CandidateFilter(MatchingSolver &solver)
            : solver(solver),
            spoke_adj(solver.qn, vector<ui>()),
            match_right(solver.max_g_deg, -1),
            seen_right(solver.max_g_deg, 0),
            left_is_bridge(solver.qn, 0)
#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
            , queued_spoke(solver.qn, 0)
#endif
        {}

        bool run()
        {
            if (!timed(&MatchingSolver::TimeStats::filter_bridge_time, [&] {
                buildBridgeIndex();
                return true;
            })) return false;

            if (!timed(&MatchingSolver::TimeStats::filter_nlf_time, [&] {
                return filterNLF();
            })) return false;

#if CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING
            if (!timed(&MatchingSolver::TimeStats::filter_bridge_time, [&] {
                return initBridgeSupports() && propagateFilterClosure();
            })) return false;
#endif

#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
            if (!timed(&MatchingSolver::TimeStats::filter_spoke_time, [&] {
                enqueueAllSpokeVertices();
                return propagateFilterClosure();
            })) return false;
#endif

            updateCandidateCount();
            return true;
        }

    private:
        template <typename Fn>
        bool timed(long long MatchingSolver::TimeStats::*field, Fn &&fn)
        {
            Timer t;
            bool ok = fn();
            solver.stats.*field += t.elapsed();
            return ok;
        }

        void updateCandidateCount()
        {
            solver.stats.filter_candidate_count = 0;
            for (ui u = 0; u < solver.qn; ++u) {
                solver.stats.filter_candidate_count += (ui)solver.candidates[u].size();
            }
        }

        ui addBridgeArc(ui from, ui to)
        {
            bridge_arcs.push_back({ from, to });
            return (ui)bridge_arcs.size() - 1;
        }

        void markBridgeNeighbor(ui from, ui to)
        {
            bool marked = false;
            const vector<ui> &neighbors = solver.q_neighbors[from];
            for (ui i = 0; i < solver.q_degree[from]; ++i) {
                if (neighbors[i] == to) {
                    solver.q_neighbor_is_bridge[from][i] = 1;
                    marked = true;
                    break;
                }
            }
            assert(marked);
            (void)marked;
        }

        void addBridge(ui a, ui b)
        {
            ui ab = addBridgeArc(a, b);
            ui ba = addBridgeArc(b, a);
            bridge_nbrs[a].push_back({ b, ba });
            bridge_nbrs[b].push_back({ a, ab });
            markBridgeNeighbor(a, b);
            markBridgeNeighbor(b, a);
        }

        void tarjan(ui u, ui parent, vector<int> &dfn, vector<int> &low, int &time)
        {
            dfn[u] = low[u] = ++time;
            for (ui v : solver.q_neighbors[u]) {
                if (dfn[v] == 0) {
                    tarjan(v, u, dfn, low, time);
                    low[u] = std::min(low[u], low[v]);
                    if (low[v] > dfn[u]) {
                        addBridge(u, v);
                    }
                }
                else if (v != parent) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
        }

        void buildBridgeIndex()
        {
            bridge_arcs.clear();
            bridge_nbrs.assign(solver.qn, vector<BridgeNbr>());
            solver.q_neighbor_is_bridge.assign(solver.qn, vector<char>());

            for (ui u = 0; u < solver.qn; ++u) {
                solver.q_neighbor_is_bridge[u].assign(solver.q_degree[u], 0);
            }

#if CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING
            vector<int> dfn(solver.qn, 0);
            vector<int> low(solver.qn, 0);
            int tim = 0;

            for (ui u = 0; u < solver.qn; ++u) {
                if (dfn[u] == 0) {
                    tarjan(u, solver.qn, dfn, low, tim);
                }
            }
#endif
        }

        void buildVertexCache()
        {
            query_label.assign(solver.qn, 0);
            for (ui u = 0; u < solver.qn; ++u) {
                query_label[u] = solver.query_graph->getVertexLabel(u);
            }

            data_label.assign(solver.gn, 0);
            data_degree.assign(solver.gn, 0);
            for (ui v = 0; v < solver.gn; ++v) {
                data_label[v] = solver.data_graph->getVertexLabel(v);
                data_degree[v] = solver.data_graph->getVertexDegree(v);
            }
        }

        void buildQueryLabelReqs()
        {
            query_label_reqs.assign(solver.qn, vector<QueryLabelReq>());
            query_bridge_degree.assign(solver.qn, 0);
            vector<ui> bridge_counts(solver.label_count, 0);
            vector<ui> non_bridge_counts(solver.label_count, 0);
            vector<ui> touched_labels;

            for (ui u = 0; u < solver.qn; ++u) {
                query_label_reqs[u].reserve(std::min(solver.q_degree[u], solver.label_count));
                touched_labels.clear();

                for (ui i = 0; i < solver.q_degree[u]; ++i) {
                    ui u1 = solver.q_neighbors[u][i];
                    LabelID label = query_label[u1];
                    assert(label >= 0 && (ui)label < solver.label_count);
                    ui label_idx = (ui)label;

                    if (bridge_counts[label_idx] == 0 && non_bridge_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    if (solver.q_neighbor_is_bridge[u][i]) {
                        query_bridge_degree[u]++;
                        bridge_counts[label_idx]++;
                    }
                    else {
                        non_bridge_counts[label_idx]++;
                    }
                }

                sort(touched_labels.begin(), touched_labels.end());
                for (ui label : touched_labels) {
                    ui bridge_count = bridge_counts[label];
                    ui non_bridge_count = non_bridge_counts[label];
                    query_label_reqs[u].push_back({
                        (LabelID)label, bridge_count, non_bridge_count
                    });
                    bridge_counts[label] = 0;
                    non_bridge_counts[label] = 0;
                }
            }
        }

        void buildDataLabelFreqs()
        {
            data_label_freqs.assign(solver.gn, vector<LabelFreq>());
            data_by_label.assign(solver.label_count, vector<ui>());
            vector<ui> label_counts(solver.label_count, 0);
            vector<ui> touched_labels;

            for (ui v = 0; v < solver.gn; ++v) {
                ui deg = 0;
                const ui *neighbors = solver.data_graph->getVertexNeighbors(v, deg);
                data_label_freqs[v].reserve(std::min(deg, solver.label_count));
                touched_labels.clear();

                LabelID vertex_label = data_label[v];
                assert(vertex_label >= 0 && (ui)vertex_label < solver.label_count);
                if (vertex_label >= 0 && (ui)vertex_label < solver.label_count) {
                    data_by_label[(ui)vertex_label].push_back(v);
                }

                for (ui i = 0; i < deg; ++i) {
                    LabelID label = data_label[neighbors[i]];
                    assert(label >= 0 && (ui)label < solver.label_count);
                    ui label_idx = (ui)label;
                    if (label_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    label_counts[label_idx]++;
                }

                sort(touched_labels.begin(), touched_labels.end());
                for (ui label : touched_labels) {
                    data_label_freqs[v].push_back({ (LabelID)label, label_counts[label] });
                    label_counts[label] = 0;
                }
            }
        }

        ui nlfDiff(ui u, ui v) const
        {
            ui diff = 0;
            const vector<QueryLabelReq> &query_reqs = query_label_reqs[u];
            const vector<LabelFreq> &data_freqs = data_label_freqs[v];
            size_t data_idx = 0;

            for (const auto &req : query_reqs) {
                while (data_idx < data_freqs.size() && data_freqs[data_idx].label < req.label) {
                    data_idx++;
                }

                ui data_count = 0;
                if (data_idx < data_freqs.size() && data_freqs[data_idx].label == req.label) {
                    data_count = data_freqs[data_idx].count;
                }

                ui bridge_need = req.bridge;
                if (bridge_need > data_count) {
                    return solver.threshold + 1;
                }

                ui non_bridge_need = req.normal;
                ui non_bridge_available = data_count - bridge_need;
                if (non_bridge_need > non_bridge_available) {
                    diff += (non_bridge_need - non_bridge_available);
                }
                if (diff > solver.threshold) {
                    return diff;
                }
            }
            return diff;
        }

        bool filterNLF()
        {
            buildVertexCache();
            buildQueryLabelReqs();
            buildDataLabelFreqs();

            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = query_label[u];
                assert(lu >= 0 && (ui)lu < data_by_label.size());
                if (lu < 0 || (ui)lu >= data_by_label.size()) return false;

                for (ui v : data_by_label[(ui)lu]) {
                    if (query_bridge_degree[u] > data_degree[v]) continue;
                    if (solver.q_degree[u] > data_degree[v] + solver.threshold) continue;
                    if (nlfDiff(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        ui &support(ui arc_id, ui v)
        {
            return bridge_support[(size_t)arc_id * solver.gn + v];
        }

        bool pruneCandidate(ui u, ui v)
        {
            if (!solver.candidates[u].contains(v)) {
                return true;
            }

            solver.candidates[u].remove(v);
            if (solver.candidates[u].empty()) {
                return false;
            }
            removed.push({ u, v });
#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
            for (ui nbr_u : solver.q_neighbors[u]) {
                enqueueSpokeVertex(nbr_u);
            }
#endif
            return true;
        }

        bool initBridgeSupports()
        {
            bridge_support.assign((size_t)bridge_arcs.size() * solver.gn, 0);
            vector<pair<ui, ui>> zero_support_candidates;

            for (ui arc_id = 0; arc_id < (ui)bridge_arcs.size(); ++arc_id) {
                const BridgeArc &arc = bridge_arcs[arc_id];
                for (ui v : solver.candidates[arc.from]) {
                    ui deg = 0;
                    const ui *nbrs = solver.data_graph->getVertexNeighbors(v, deg);
                    ui support = 0;
                    for (ui i = 0; i < deg; ++i) {
                        if (solver.candidates[arc.to].contains(nbrs[i])) {
                            support++;
                        }
                    }
                    this->support(arc_id, v) = support;
                    if (support == 0) {
                        zero_support_candidates.push_back({ arc.from, v });
                    }
                }
            }

            for (const auto &candidate : zero_support_candidates) {
                if (!pruneCandidate(candidate.first, candidate.second)) {
                    return false;
                }
            }
            return true;
        }

        bool propagateBridgeRemovals()
        {
            while (!removed.empty()) {
                ui removed_u = removed.front().first;
                ui removed_v = removed.front().second;
                removed.pop();

                ui deg = 0;
                const ui *nbrs = solver.data_graph->getVertexNeighbors(removed_v, deg);

                for (const BridgeNbr &bridge_nbr : bridge_nbrs[removed_u]) {
                    ui affected_u = bridge_nbr.to;
                    ui arc_id = bridge_nbr.support_arc;
                    for (ui i = 0; i < deg; ++i) {
                        ui v = nbrs[i];
                        if (!solver.candidates[affected_u].contains(v)) {
                            continue;
                        }
                        ui &support_count = support(arc_id, v);
                        if (support_count == 0) {
                            continue;
                        }
                        support_count--;
                        if (support_count == 0 &&
                            !pruneCandidate(affected_u, v)) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        bool propagateFilterClosure()
        {
            while (true) {
                if (!propagateBridgeRemovals()) {
                    return false;
                }
#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
                if (pending_spokes.empty()) {
                    break;
                }

                ui u = pending_spokes.front();
                pending_spokes.pop();
                queued_spoke[u] = 0;
                if (!processSpokeVertex(u)) {
                    return false;
                }
#else
                break;
#endif
            }
            return true;
        }

        bool augmentSpoke(ui left_idx)
        {
            for (ui right_idx : spoke_adj[left_idx]) {
                if (seen_right[right_idx] == seen_token) continue;
                seen_right[right_idx] = seen_token;
                if (match_right[right_idx] < 0 ||
                    augmentSpoke((ui)match_right[right_idx])) {
                    match_right[right_idx] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        bool tryAugmentSpoke(ui left_idx)
        {
            seen_token++;
            if (seen_token == 0) {
                std::fill(seen_right.begin(), seen_right.end(), 0);
                seen_token = 1;
            }
            return augmentSpoke(left_idx);
        }

        void buildSpokeAdj(ui u, ui v, ui &deg_u, ui &deg_v)
        {
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            deg_u = solver.q_degree[u];
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            for (ui i = 0; i < deg_u; ++i) {
                spoke_adj[i].clear();
                ui u1 = u_neighbors[i];
                left_is_bridge[i] = solver.q_neighbor_is_bridge[u][i];
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        spoke_adj[i].push_back(j);
                    }
                }
            }
        }

        bool spokeFeasible(ui u, ui v, ui budget)
        {
            ui deg_u = 0;
            ui deg_v = 0;
            buildSpokeAdj(u, v, deg_u, deg_v);

            std::fill(match_right.begin(), match_right.begin() + deg_v, -1);

            ui optional_count = 0;
            for (ui i = 0; i < deg_u; ++i) {
                if (left_is_bridge[i]) {
                    if (!tryAugmentSpoke(i)) {
                        return false;
                    }
                }
                else {
                    optional_count++;
                }
            }

            if (optional_count <= budget) {
                return true;
            }

            ui required = optional_count - budget;
            ui matched = 0;
            ui processed = 0;

            for (ui i = 0; i < deg_u; ++i) {
                if (left_is_bridge[i]) {
                    continue;
                }

                processed++;
                if (tryAugmentSpoke(i)) {
                    matched++;
                    if (matched >= required) {
                        return true;
                    }
                }

                ui remaining = optional_count - processed;
                if (matched + remaining < required) {
                    return false;
                }
            }

            return matched >= required;
        }

        ui maxMissingIncidentEdges(ui u) const
        {
            if (solver.q_degree[u] == 0) return 0;
            return std::min(solver.threshold, solver.q_degree[u] - 1);
        }

#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
        void enqueueSpokeVertex(ui u)
        {
            if (queued_spoke[u]) {
                return;
            }
            pending_spokes.push(u);
            queued_spoke[u] = 1;
        }

        void enqueueAllSpokeVertices()
        {
            for (ui u = 0; u < solver.qn; ++u) {
                enqueueSpokeVertex(u);
            }
        }

        bool processSpokeVertex(ui u)
        {
            ui budget = maxMissingIncidentEdges(u);

            vector<ui> to_remove;
            for (ui v : solver.candidates[u]) {
                if (!spokeFeasible(u, v, budget)) {
                    to_remove.push_back(v);
                }
            }
            if (to_remove.empty()) return true;

            for (ui v : to_remove) {
                if (!pruneCandidate(u, v)) {
                    return false;
                }
            }
            return true;
        }
#endif
    };

    bool runCandidateFiltering()
    {
        return CandidateFilter(*this).run();
    }

    void initDataLabelDegreeIndex()
    {
        data_label_degrees.assign(gn, vector<DataLabelDegreeCount>());
        vector<ui> label_counts(label_count, 0);

        for (ui v = 0; v < gn; ++v) {
            ui deg = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(v, deg);
            data_label_degrees[v].reserve(std::min(deg, label_count));

            for (ui i = 0; i < deg; ++i) {
                LabelID label = data_graph->getVertexLabel(neighbors[i]);
                assert(label >= 0 && (ui)label < label_count);
                label_counts[(ui)label]++;
            }

            for (ui label = 0; label < label_count; ++label) {
                ui count = label_counts[label];
                if (count != 0) {
                    data_label_degrees[v].push_back({ (LabelID)label, count });
                    label_counts[label] = 0;
                }
            }
        }
    }

    ui dataLabelDegree(ui v, LabelID label) const
    {
        if (v >= gn || label < 0 || (ui)label >= label_count) {
            return 0;
        }

        const vector<DataLabelDegreeCount> &counts = data_label_degrees[v];
        auto it = std::lower_bound(counts.begin(), counts.end(), label,
            [](const DataLabelDegreeCount &item, LabelID target_label) {
                return item.label < target_label;
            });
        if (it == counts.end() || it->label != label) {
            return 0;
        }
        return it->count;
    }

    ui cheapEdgeSupport(ui u, ui anchor)
    {
        if (u >= qn || anchor >= qn || mapped_q[anchor] == -1 ||
            excluded_edges[u][anchor]) {
            return 0;
        }

        LabelID label = query_graph->getVertexLabel(u);
        ui label_degree = dataLabelDegree((ui)mapped_q[anchor], label);
        ui candidate_count = (ui)candidates[u].size();
        return std::min(candidate_count, label_degree);
    }

#if CDE_EDGE_IE_FIXED_ORDER
    struct FixedEdgePriorityEntry {
        ui u = 0;
        ui anchor = 0;
        unsigned long long pair_support = 0;
        ui u_candidate_count = 0;
        ui anchor_candidate_count = 0;
    };

    unsigned long long countStaticCandidateEdgePairs(ui u, ui u1) const
    {
        unsigned long long count = 0;
        for (ui v : candidates[u]) {
            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(v, degree);
            for (ui i = 0; i < degree; ++i) {
                if (candidates[u1].contains(neighbors[i])) {
                    count++;
                }
            }
        }
        return count;
    }

    void initFixedEdgePriorities()
    {
        const ui invalid_priority = std::numeric_limits<ui>::max();
        fixed_edge_priority.assign(qn, vector<ui>(qn, invalid_priority));

        vector<FixedEdgePriorityEntry> entries;
        size_t directed_edge_count = 0;
        for (ui u = 0; u < qn; ++u) {
            directed_edge_count += q_neighbors[u].size();
        }
        entries.reserve(directed_edge_count);

        for (ui u = 0; u < qn; ++u) {
            for (ui u1 : q_neighbors[u]) {
                if (u >= u1) {
                    continue;
                }

                ui u_candidate_count = (ui)candidates[u].size();
                ui u1_candidate_count = (ui)candidates[u1].size();
                ui scan_u = u_candidate_count <= u1_candidate_count ? u : u1;
                ui target_u = scan_u == u ? u1 : u;
                unsigned long long pair_support =
                    countStaticCandidateEdgePairs(scan_u, target_u);

                FixedEdgePriorityEntry forward;
                forward.u = u;
                forward.anchor = u1;
                forward.pair_support = pair_support;
                forward.u_candidate_count = u_candidate_count;
                forward.anchor_candidate_count = u1_candidate_count;
                entries.push_back(forward);

                FixedEdgePriorityEntry reverse;
                reverse.u = u1;
                reverse.anchor = u;
                reverse.pair_support = pair_support;
                reverse.u_candidate_count = u1_candidate_count;
                reverse.anchor_candidate_count = u_candidate_count;
                entries.push_back(reverse);
            }
        }

        std::sort(entries.begin(), entries.end(),
            [&](const FixedEdgePriorityEntry &lhs,
                const FixedEdgePriorityEntry &rhs) {
                // pair_support / |Cand(anchor)| estimates the support seen after
                // fixing one candidate for the anchor. Compare the ratios
                // exactly so the priority is deterministic.
                __uint128_t lhs_scaled =
                    (__uint128_t)lhs.pair_support *
                    std::max((ui)1, rhs.anchor_candidate_count);
                __uint128_t rhs_scaled =
                    (__uint128_t)rhs.pair_support *
                    std::max((ui)1, lhs.anchor_candidate_count);
                if (lhs_scaled != rhs_scaled) {
                    return lhs_scaled < rhs_scaled;
                }
                if (lhs.u_candidate_count != rhs.u_candidate_count) {
                    return lhs.u_candidate_count < rhs.u_candidate_count;
                }
                if (lhs.pair_support != rhs.pair_support) {
                    return lhs.pair_support < rhs.pair_support;
                }
                if (q_degree[lhs.u] != q_degree[rhs.u]) {
                    return q_degree[lhs.u] > q_degree[rhs.u];
                }
                if (q_degree[lhs.anchor] != q_degree[rhs.anchor]) {
                    return q_degree[lhs.anchor] > q_degree[rhs.anchor];
                }
                if (lhs.u != rhs.u) {
                    return lhs.u < rhs.u;
                }
                return lhs.anchor < rhs.anchor;
            });

        assert(entries.size() <= (size_t)std::numeric_limits<ui>::max());
        for (size_t rank = 0; rank < entries.size(); ++rank) {
            const FixedEdgePriorityEntry &entry = entries[rank];
            fixed_edge_priority[entry.u][entry.anchor] = (ui)rank;
        }
    }
#endif
    // ========================================================================

    // ========================================================================
    // Branching
    // ========================================================================
    struct BranchSelector {
    private:
        MatchingSolver &solver;

    public:
        explicit BranchSelector(MatchingSolver &solver) : solver(solver) {}

        ui selectInitialRoot() const
        {
            ui root = 0;
            for (ui u = 1; u < solver.qn; ++u) {
                if (solver.candidates[u].size() < solver.candidates[root].size() ||
                    (solver.candidates[u].size() == solver.candidates[root].size() &&
                        solver.q_degree[u] > solver.q_degree[root])) {
                    root = u;
                }
            }
            return root;
        }

        bool collectTopActiveEdges(const vector<ui> &component_frontier, ui max_count,
            vector<ActiveEdge> &top_edges, EdgeScoreCache *edge_score_cache = nullptr,
            const vector<char> *skip_query_vertices = nullptr) const
        {
            top_edges.clear();
            if (max_count == 0) {
                return false;
            }

            vector<ActiveEdge> uncached_edges;
            for (ui u : component_frontier) {
                if (u >= solver.qn) {
                    continue;
                }
                if (shouldSkipQueryVertex(u, skip_query_vertices)) {
                    continue;
                }
                if (edge_score_cache != nullptr) {
                    const vector<ActiveEdge> &cached_edges = cachedActiveEdgesForVertex(u, *edge_score_cache);
                    top_edges.insert(top_edges.end(), cached_edges.begin(), cached_edges.end());
                }
                else {
                    collectActiveEdgesForVertex(u, uncached_edges);
                    top_edges.insert(top_edges.end(), uncached_edges.begin(), uncached_edges.end());
                }
            }

            if (top_edges.empty()) {
                return false;
            }

            auto better_edge = [&](const ActiveEdge &lhs, const ActiveEdge &rhs) {
                return isBetterActiveEdge(lhs, rhs);
                };
            if (top_edges.size() > max_count) {
                partial_sort(top_edges.begin(), top_edges.begin() + max_count,
                    top_edges.end(), better_edge);
                top_edges.resize(max_count);
            }
            else {
                sort(top_edges.begin(), top_edges.end(), better_edge);
            }
            return true;
        }

        bool restrictTopEdgesToCoveredComponent(vector<ActiveEdge> &top_edges,
            EdgeScoreCache &edge_score_cache, vector<int> &component_id,
            vector<vector<ui>> &component_frontiers,
            vector<ui> &component_edge_counts, vector<ui> &component_seen_counts,
            vector<double> &component_support_sums,
            double &covered_component_support_sum, bool &has_zero_support_component,
            const vector<char> *skip_query_vertices = nullptr) const
        {
            covered_component_support_sum = std::numeric_limits<double>::max();
            has_zero_support_component = false;
            if (top_edges.empty()) {
                return false;
            }

            size_t component_count = labelUnmatchedComponents(component_id,
                component_frontiers, skip_query_vertices);
            component_edge_counts.assign(component_count, 0);
            component_seen_counts.assign(component_count, 0);
            component_support_sums.assign(component_count, 0.0);

            for (size_t id = 0; id < component_count; ++id) {
                for (ui component_u : component_frontiers[id]) {
                    const vector<ActiveEdge> &edges =
                        cachedActiveEdgesForVertex(component_u, edge_score_cache);
                    component_edge_counts[id] += (ui)edges.size();
#if !CDE_EDGE_IE_FIXED_ORDER
                    for (const ActiveEdge &edge : edges) {
                        component_support_sums[id] += edge.anchor_support;
                    }
#endif
                }
            }

            for (size_t edge_idx = 0; edge_idx < top_edges.size(); ++edge_idx) {
                const ActiveEdge &edge = top_edges[edge_idx];
                if (edge.u >= component_id.size() || component_id[edge.u] < 0) {
                    continue;
                }

                size_t id = (size_t)component_id[edge.u];
                if (id >= component_edge_counts.size() ||
                    component_edge_counts[id] == 0) {
                    continue;
                }

                component_seen_counts[id]++;
                if (component_seen_counts[id] != component_edge_counts[id]) {
                    continue;
                }

#if !CDE_EDGE_IE_FIXED_ORDER
                if (component_support_sums[id] == 0.0) {
                    has_zero_support_component = true;
                    return false;
                }
#endif

                covered_component_support_sum = component_support_sums[id];
                top_edges.resize(edge_idx + 1);
                return true;
            }

            return false;
        }

    private:
        bool shouldSkipQueryVertex(ui u, const vector<char> *skip_query_vertices) const
        {
            return skip_query_vertices != nullptr &&
                u < skip_query_vertices->size() && (*skip_query_vertices)[u];
        }

        ui countLiveAnchors(ui u) const
        {
            return u < solver.qn ? solver.anchor_count[u] : 0;
        }

        ui estimateEdgeSupport(ui u, ui anchor) const
        {
            Timer t_score;
            ui estimate = solver.cheapEdgeSupport(u, anchor);
            long long elapsed = t_score.elapsed();
            solver.stats.frontier_score_cheap_support_time += elapsed;
            solver.stats.frontier_score_time += elapsed;
            return estimate;
        }

        bool isBetterActiveEdge(const ActiveEdge &lhs, const ActiveEdge &rhs) const
        {
#if CDE_EDGE_IE_FIXED_ORDER
            ui lhs_priority = solver.fixed_edge_priority[lhs.u][lhs.anchor];
            ui rhs_priority = solver.fixed_edge_priority[rhs.u][rhs.anchor];
            if (lhs_priority != rhs_priority) {
                return lhs_priority < rhs_priority;
            }
            if (lhs.u != rhs.u) {
                return lhs.u < rhs.u;
            }
            return lhs.anchor < rhs.anchor;
#else
            if (lhs.anchor_support != rhs.anchor_support) {
                return lhs.anchor_support < rhs.anchor_support;
            }
            if (lhs.live_anchor_count != rhs.live_anchor_count) {
                return lhs.live_anchor_count > rhs.live_anchor_count;
            }
            if (lhs.query_degree != rhs.query_degree) {
                return lhs.query_degree > rhs.query_degree;
            }
            if (lhs.u != rhs.u) {
                return lhs.u < rhs.u;
            }
            return lhs.anchor < rhs.anchor;
#endif
        }

        void collectActiveEdgesForVertex(ui u, vector<ActiveEdge> &edges) const
        {
            edges.clear();
            if (u >= solver.qn || solver.mapped_q[u] != -1 || solver.frontier_pos[u] == -1) {
                return;
            }

            ui live_anchor_count = countLiveAnchors(u);
            if (live_anchor_count == 0) {
                return;
            }

            for (ui anchor : solver.q_neighbors[u]) {
                if (solver.mapped_q[anchor] == -1 || solver.excluded_edges[u][anchor]) {
                    continue;
                }

                ActiveEdge edge;
                edge.u = u;
                edge.anchor = anchor;
#if !CDE_EDGE_IE_FIXED_ORDER
                edge.anchor_support = estimateEdgeSupport(u, anchor);
#endif
                edge.live_anchor_count = live_anchor_count;
                edge.query_degree = solver.q_degree[u];
                edges.push_back(edge);
            }
        }

        const vector<ActiveEdge> &cachedActiveEdgesForVertex(ui u, EdgeScoreCache &cache) const
        {
            if (cache.active_edges_by_vertex.empty()) {
                cache.active_edges_by_vertex.resize(solver.qn);
                cache.active_edges_cached.assign(solver.qn, 0);
            }
            if (!cache.active_edges_cached[u]) {
                collectActiveEdgesForVertex(u, cache.active_edges_by_vertex[u]);
                cache.active_edges_cached[u] = 1;
            }
            return cache.active_edges_by_vertex[u];
        }

        size_t labelUnmatchedComponents(vector<int> &component_id,
            vector<vector<ui>> &component_frontiers,
            const vector<char> *skip_query_vertices) const
        {
            // Label every unmatched component in one graph sweep, then reuse the
            // labels for all top edges in this DFS state.
            component_id.assign(solver.qn, -1);
            for (vector<ui> &frontier : component_frontiers) {
                frontier.clear();
            }

            size_t component_count = 0;
            queue<ui> q;
            for (ui start = 0; start < solver.qn; ++start) {
                if (solver.mapped_q[start] != -1 ||
                    shouldSkipQueryVertex(start, skip_query_vertices) ||
                    component_id[start] != -1) {
                    continue;
                }

                int id = (int)component_count++;
                if ((size_t)id == component_frontiers.size()) {
                    component_frontiers.emplace_back();
                }
                component_id[start] = id;
                q.push(start);

                while (!q.empty()) {
                    ui curr = q.front();
                    q.pop();

                    if (solver.frontier_pos[curr] != -1) {
                        component_frontiers[(size_t)id].push_back(curr);
                    }

                    for (ui nbr : solver.q_neighbors[curr]) {
                        if (solver.mapped_q[nbr] == -1 &&
                            !shouldSkipQueryVertex(nbr, skip_query_vertices) &&
                            component_id[nbr] == -1) {
                            component_id[nbr] = id;
                            q.push(nbr);
                        }
                    }
                }
            }
            return component_count;
        }

    };

    ui nextDataVertexMarkToken()
    {
        if (++data_vertex_mark_token == 0) {
            std::fill(data_vertex_mark.begin(), data_vertex_mark.end(), 0);
            data_vertex_mark_token = 1;
        }
        return data_vertex_mark_token;
    }

    ui countAnchorsByMark(ui u, ui selected_anchor, const vector<ui> &candidate_vertices, vector<ui> &anchor_counts)
    {
        anchor_counts.assign(candidate_vertices.size(), 0);
        ui token = nextDataVertexMarkToken();
        for (ui i = 0; i < candidate_vertices.size(); ++i) {
            ui v = candidate_vertices[i];
            assert(v < gn);
            data_vertex_mark[v] = token;
            data_vertex_mark_pos[v] = i;
        }
        ui anchor_num = 0;

        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;

            anchor_num++;
            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (data_vertex_mark[v] == token) {
                    anchor_counts[data_vertex_mark_pos[v]]++;
                }
            }
        }

        return anchor_num;
    }

    ui countAnchorsByHasEdge(ui u, ui selected_anchor, const vector<ui> &candidate_vertices, vector<ui> &anchor_counts) const
    {
        anchor_counts.assign(candidate_vertices.size(), 0);
        ui anchor_num = 0;
        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
            anchor_num++;
        }

        for (ui i = 0; i < candidate_vertices.size(); ++i) {
            ui v = candidate_vertices[i];
            for (ui anchor : q_neighbors[u]) {
                if (anchor == selected_anchor) continue;
                if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
                if (data_graph->hasEdge(v, (ui)mapped_q[anchor])) {
                    anchor_counts[i]++;
                }
            }
        }

        return anchor_num;
    }

    // heuristic to decide whether to use mark-based counting or has-edge checking.
    bool useMarkForAnchorCount(ui u, ui selected_anchor, size_t vertex_count) const
    {
        auto binary_search_cost = [](ui degree) -> size_t {
            size_t cost = 1;
            for (; degree > 1; degree >>= 1) ++cost;
            return cost;
        };

        size_t mark_cost = vertex_count * 2;
        size_t has_edge_cost = 0;

        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;

            ui mapped_anchor = (ui)mapped_q[anchor];
            ui degree = data_graph->getVertexDegree(mapped_anchor);

            mark_cost += degree;
            has_edge_cost += vertex_count * binary_search_cost(degree);
        }

        return mark_cost <= has_edge_cost;
    }

    // count u's candidate vertices that are adjacent to each anchor's mapped data vertex
    // and store the counts in anchor_counts
    // return the number of anchors in query graph excluding selected_anchor and excluded ones.
    ui countAnchors(ui u, ui ua, const vector<ui> &candidate_vertices, vector<ui> &anchor_counts)
    {
        if (candidate_vertices.empty()) return 0;
        if (useMarkForAnchorCount(u, ua, candidate_vertices.size())) {
            return countAnchorsByMark(u, ua, candidate_vertices, anchor_counts);
        }
        return countAnchorsByHasEdge(u, ua, candidate_vertices, anchor_counts);
    }

    // ========================================================================
    // Maintain Frontier and Anchor Support
    // ========================================================================
    // calculate support[u][anchor] or calculate all candidates of u by (u, anchor)
    template <typename Visitor>
    ui calEdgeSupport(ui u, ui anchor, Visitor visit) const
    {
        if (u >= qn || anchor >= qn || mapped_q[anchor] == -1 || excluded_edges[u][anchor]) {
            return 0;
        }

        ui count = 0;
        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
        for (ui i = 0; i < deg; ++i) {
            ui v = nbrs[i];
            if (!candidates[u].contains(v)) continue;
            if (mapped_g[v] != -1) continue;
            if (excluded_cands[u].contains(v)) continue;
            count++;
            visit(v);
        }
        return count;
    }

    inline void recordSupportSnapshot(ui u, ui anchor)
    {
        if (active_support_snapshots == nullptr) {
            return;
        }
        SupportSnapshot snapshot;
        snapshot.u = u;
        snapshot.anchor = anchor;
        snapshot.value = support[u][anchor];
        snapshot.dirty = support_dirty[u][anchor];
        active_support_snapshots->push_back(snapshot);
    }

    inline void markSupportDirty(ui u, ui anchor)
    {
        if (u >= qn || anchor >= qn || support_dirty[u][anchor]) {
            return;
        }
        recordSupportSnapshot(u, anchor);
        support_dirty[u][anchor] = 1;
    }

    inline void markLiveAnchorSupportDirty(ui u)
    {
        if (u >= qn || mapped_q[u] != -1 || frontier_pos[u] == -1) {
            return;
        }
        for (ui anchor : q_neighbors[u]) {
            if (mapped_q[anchor] != -1 && !excluded_edges[u][anchor]) {
                markSupportDirty(u, anchor);
            }
        }
    }

    // update active_frontier, frontier_pos, anchor_support
    inline void updateFrontierStatus(ui u)
    {
        bool should_be = (mapped_q[u] == -1 && anchor_count[u] > 0);
        bool is_in = (frontier_pos[u] != -1);

        if (should_be && !is_in) {
            // add u to frontier
            frontier_pos[u] = active_frontier.size();
            active_frontier.push_back(u);
        }
        else if (!should_be && is_in) {
            // remove u from frontier
            ui idx = frontier_pos[u];
            ui last_u = active_frontier.back();
            active_frontier[idx] = last_u;
            frontier_pos[last_u] = idx;
            active_frontier.pop_back();
            frontier_pos[u] = -1;
        }
    }

    // when u becomes matched/unmatched, update anchor_count, active_frontier and anchor_support
    void updateFrontier(ui u, bool matched)
    {
        Timer t_support_update;
        if (matched) updateFrontierStatus(u);

        // update u's neighbors in frontier
        for (ui u1 : q_neighbors[u]) {
            if(excluded_edges[u][u1]) continue;

            if (matched) {
                anchor_count[u1]++;
                updateFrontierStatus(u1);

                if (mapped_q[u1] == -1 && frontier_pos[u1] != -1) {
                    markSupportDirty(u1, u);
                }
            }
            else {
                anchor_count[u1]--;
                updateFrontierStatus(u1);
            }
        }

        if (!matched) updateFrontierStatus(u);
        stats.support_update_time += t_support_update.elapsed();
    }

    // mark (u, anchor) as excluded, update anchor_count[u], active_frontier and anchor_support
    void excludeFrontierEdge(ui u, ui anchor)
    {
        excluded_edges[u][anchor] = 1;
        excluded_edges[anchor][u] = 1;
        assert(anchor_count[u] > 0);
        anchor_count[u]--;
        if(anchor_count[u] == 0) updateFrontierStatus(u);
        markLiveAnchorSupportDirty(u);
    }

    // restore (u, anchor), update anchor_count[u], active_frontier and anchor_support
    void restoreFrontierEdge(ui u, ui anchor)
    {
        excluded_edges[u][anchor] = 0;
        excluded_edges[anchor][u] = 0;
        anchor_count[u]++;
        if(anchor_count[u] == 1) updateFrontierStatus(u);
    }

    struct SupportUndoScope {
        MatchingSolver &solver;
        vector<SupportSnapshot> *previous_snapshots;
        vector<SupportSnapshot> &snapshots;

        SupportUndoScope(MatchingSolver &solver, vector<SupportSnapshot> &snapshots)
            : solver(solver),
            previous_snapshots(solver.active_support_snapshots),
            snapshots(snapshots)
        {
            solver.active_support_snapshots = &snapshots;
        }

        ~SupportUndoScope()
        {
            for (auto it = snapshots.rbegin(); it != snapshots.rend(); ++it) {
                solver.support[it->u][it->anchor] = it->value;
                solver.support_dirty[it->u][it->anchor] = it->dirty;
            }
            snapshots.clear();
            solver.active_support_snapshots = previous_snapshots;
        }
    };
    // ========================================================================

    // ========================================================================
    // DFS Buffer
    // ========================================================================
    struct DfsBuffer {
        vector<ActiveEdge> top_edges;
        vector<ui> candidate_vertices;
        vector<ui> candidate_anchor_counts;
        vector<int> component_id;
        vector<vector<ui>> component_frontiers;
        vector<ui> component_edge_counts;
        vector<ui> component_seen_counts;
        vector<double> component_support_sums;
        vector<char> terminal_skip;
        vector<ui> terminal_vertices;
        vector<ui> active_terminal_vertices;
        vector<TerminalTailVertex> terminal_tail_vertices;
        vector<ui> terminal_candidate_vertices;
        vector<ui> terminal_candidate_support_counts;
        vector<vector<ui>> candidate_delta_buckets;
        vector<ui> candidate_delta_touched;
        EdgeScoreCache edge_score_cache;
        BranchSelector branch_selector;
        vector<SupportSnapshot> local_support_snapshots;
        explicit DfsBuffer(MatchingSolver &solver)
            : branch_selector(solver)
        {}

        // reserve space for dfs buffers
        void reserve(ui threshold, ui max_g_deg, ui qn, ui gn)
        {
            top_edges.reserve((size_t)threshold + 1);
            candidate_vertices.reserve(max_g_deg);
            candidate_anchor_counts.reserve(max_g_deg);
            component_id.reserve(qn);
            component_frontiers.reserve(qn);
            component_edge_counts.reserve(qn);
            component_seen_counts.reserve(qn);
            component_support_sums.reserve(qn);
            terminal_skip.reserve(qn);
            terminal_vertices.reserve(qn);
            active_terminal_vertices.reserve(qn);
            terminal_tail_vertices.reserve(qn);
            terminal_candidate_vertices.reserve(max_g_deg);
            terminal_candidate_support_counts.reserve(max_g_deg);
            candidate_delta_buckets.resize((size_t)threshold + 1);
            candidate_delta_touched.reserve((size_t)threshold + 1);
            local_excluded_edges.reserve((size_t)threshold + 1);
            size_t cand_size = std::min((size_t)gn, ((size_t)threshold + 1) * (size_t)max_g_deg);
            local_excluded_cands.reserve(cand_size);
        }

        void clearLocal()
        {
            clearCandidateDeltaBuckets();
            for (vector<ActiveEdge> &edges : edge_score_cache.active_edges_by_vertex) {
                edges.clear();
            }
            std::fill(edge_score_cache.active_edges_cached.begin(),
                edge_score_cache.active_edges_cached.end(), 0);
            component_id.clear();
            component_edge_counts.clear();
            component_seen_counts.clear();
            component_support_sums.clear();
            terminal_vertices.clear();
            active_terminal_vertices.clear();
            terminal_tail_vertices.clear();
            terminal_candidate_vertices.clear();
            terminal_candidate_support_counts.clear();
            local_excluded_edges.clear();
            local_excluded_cands.clear();
            local_support_snapshots.clear();
        }

        void clearCandidateDeltaBuckets()
        {
            for (ui delta : candidate_delta_touched) {
                if (delta < candidate_delta_buckets.size()) {
                    candidate_delta_buckets[delta].clear();
                }
            }
            candidate_delta_touched.clear();
        }

        void addCandidateDelta(ui delta, ui candidate)
        {
            assert(delta < candidate_delta_buckets.size());
            vector<ui> &bucket = candidate_delta_buckets[delta];
            if (bucket.empty()) {
                candidate_delta_touched.push_back(delta);
            }
            bucket.push_back(candidate);
        }

        void recordExcludedEdge(ui u, ui anchor)
        {
            local_excluded_edges.push_back({ u, anchor });
        }

        void recordExcludedCands(ui u, ui v)
        {
            local_excluded_cands.push_back({ u, v });
        }

        void restoreLocalChanges(MatchingSolver &solver)
        {
            for (auto &p : local_excluded_cands) {
                solver.excluded_cands[p.first].remove(p.second);
            }

            for (auto it = local_excluded_edges.rbegin(); it != local_excluded_edges.rend(); ++it) {
                const auto &e = *it;
                ui u = e.first;
                ui ua = e.second;
                assert(solver.mapped_q[u] == -1);
                assert(solver.mapped_q[ua] != -1);

                solver.restoreFrontierEdge(u, ua);
            }
        }

    private:
        vector<pair<ui, ui>> local_excluded_edges;
        vector<pair<ui, ui>> local_excluded_cands;
    };

    vector<DfsBuffer> dfs_buffers;

    void reserveDfsBuffer(DfsBuffer &buf) const
    {
        buf.reserve(threshold, max_g_deg, qn, gn);
    }

    void initDfsBuffer()
    {
        dfs_buffers.clear();
        dfs_buffers.reserve((size_t)qn + 1);
    }

    DfsBuffer &dfsBufferForDepth(size_t depth)
    {
        assert(depth <= qn);
        assert(dfs_buffers.capacity() >= (size_t)qn + 1);
        while (dfs_buffers.size() <= depth) {
            dfs_buffers.emplace_back(*this);
            reserveDfsBuffer(dfs_buffers.back());
        }
        return dfs_buffers[depth];
    }
    // ========================================================================

    // ========================================================================
    // Terminal-tail enumeration
    // ========================================================================
    bool isTerminalQueryVertex(ui u) const
    {
        if (u >= qn || mapped_q[u] != -1 || anchor_count[u] == 0) {
            return false;
        }

        for (ui nbr : q_neighbors[u]) {
            if (excluded_edges[u][nbr]) {
                continue;
            }
            if (mapped_q[nbr] == -1) {
                return false;
            }
        }
        return true;
    }

    TerminalScan markTerminalVertices(DfsBuffer &buf)
    {
        TerminalScan scan;
        if (buf.terminal_skip.size() != qn) {
            buf.terminal_skip.assign(qn, 0);
        }
        else {
            std::fill(buf.terminal_skip.begin(), buf.terminal_skip.end(), 0);
        }
        buf.terminal_vertices.clear();
        buf.active_terminal_vertices.clear();

        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] != -1) {
                continue;
            }

            scan.unmatched_count++;
            if (isTerminalQueryVertex(u)) {
                scan.terminal_count++;
                buf.terminal_vertices.push_back(u);
                if (frontier_pos[u] != -1) {
                    buf.terminal_skip[u] = 1;
                    buf.active_terminal_vertices.push_back(u);
                    scan.terminal_frontier_count++;
                }
            }
            else if (frontier_pos[u] != -1) {
                scan.nonterminal_frontier_count++;
            }
        }
        return scan;
    }

    ui collectTerminalCandidateSupports(ui u, vector<ui> &candidate_vertices,
        vector<ui> &candidate_support_counts)
    {
        candidate_vertices.clear();
        candidate_support_counts.clear();

        ui token = nextDataVertexMarkToken();
        ui live_anchor_count = 0;

        for (ui anchor : q_neighbors[u]) {
            if (excluded_edges[u][anchor] || mapped_q[anchor] == -1) {
                continue;
            }
            live_anchor_count++;

            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];

                if (!candidates[u].contains(v)) {
                    continue;
                }
                if (mapped_g[v] != -1) {
                    continue;
                }
                if (excluded_cands[u].contains(v)) {
                    continue;
                }

                if (data_vertex_mark[v] == token) {
                    candidate_support_counts[data_vertex_mark_pos[v]]++;
                    continue;
                }

                data_vertex_mark[v] = token;
                data_vertex_mark_pos[v] = (ui)candidate_vertices.size();
                candidate_vertices.push_back(v);
                candidate_support_counts.push_back(1);
                stats.terminal_bucket_candidate_checks++;
            }
        }
        return live_anchor_count;
    }

    bool buildTerminalCandidateBuckets(ui u, ui cost, vector<vector<ui>> &buckets,
        ui &feasible_count, ui &min_delta, vector<ui> &candidate_vertices,
        vector<ui> &candidate_support_counts)
    {
        Timer t_terminal_bucket;
        assert(cost <= threshold);
        ui remaining_budget = threshold - cost;
        buckets.assign((size_t)remaining_budget + 1, vector<ui>());
        feasible_count = 0;
        min_delta = std::numeric_limits<ui>::max();

        Timer t_collect_support;
        ui live_anchor_count = collectTerminalCandidateSupports(u,
            candidate_vertices, candidate_support_counts);
        stats.terminal_candidate_collect_time += t_collect_support.elapsed();
        if (live_anchor_count == 0) {
            stats.terminal_bucket_time += t_terminal_bucket.elapsed();
            return false;
        }

        for (size_t i = 0; i < candidate_vertices.size(); ++i) {
            ui support_count = candidate_support_counts[i];
            if (support_count == 0 || support_count > live_anchor_count) {
                continue;
            }
            ui missing_delta = live_anchor_count - support_count;
            if (missing_delta > remaining_budget) {
                continue;
            }

            buckets[missing_delta].push_back(candidate_vertices[i]);
            feasible_count++;
            if (missing_delta < min_delta) {
                min_delta = missing_delta;
            }
        }

        stats.terminal_bucket_time += t_terminal_bucket.elapsed();
        return feasible_count > 0;
    }

    bool buildTerminalTailVertices(const vector<ui> &terminal_vertices, ui cost,
        DfsBuffer &buf)
    {
        vector<TerminalTailVertex> &tail_vertices = buf.terminal_tail_vertices;
        tail_vertices.clear();
        tail_vertices.reserve(terminal_vertices.size());

        for (ui u : terminal_vertices) {
            TerminalTailVertex tail_vertex;
            tail_vertex.u = u;
            if (!buildTerminalCandidateBuckets(u, cost, tail_vertex.buckets,
                tail_vertex.feasible_count, tail_vertex.min_delta,
                buf.terminal_candidate_vertices,
                buf.terminal_candidate_support_counts)) {
                return false;
            }
            tail_vertices.push_back(std::move(tail_vertex));
        }

        std::sort(tail_vertices.begin(), tail_vertices.end(),
            [this](const TerminalTailVertex &lhs, const TerminalTailVertex &rhs) {
                if (lhs.feasible_count != rhs.feasible_count) {
                    return lhs.feasible_count < rhs.feasible_count;
                }
                if (lhs.min_delta != rhs.min_delta) {
                    return lhs.min_delta < rhs.min_delta;
                }
                if (q_degree[lhs.u] != q_degree[rhs.u]) {
                    return q_degree[lhs.u] > q_degree[rhs.u];
                }
                return lhs.u < rhs.u;
            });
        return true;
    }

    ui terminalMinDeltaSum(const vector<TerminalTailVertex> &tail_vertices,
        ui limit) const
    {
        ui sum = 0;
        for (const TerminalTailVertex &tail_vertex : tail_vertices) {
            if (sum > limit || tail_vertex.min_delta > limit - sum) {
                return limit + 1;
            }
            sum += tail_vertex.min_delta;
        }
        return sum;
    }

    void recordTerminalPrune()
    {
        stats.prun_calls++;
        stats.terminal_prune_calls++;
    }

    void enumerateTerminalTail(size_t pos, ui cost,
        vector<TerminalTailVertex> &tail_vertices)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }

        stats.terminal_tail_calls++;
        assert(cost <= threshold);

        if (pos == tail_vertices.size()) {
            assert(part_M.size() == qn);
            stats.result_count++;
            noteOutputLimitIfReached();
#ifndef NDEBUG
            results_ptr->push_back(part_M);
#endif
            return;
        }

        TerminalTailVertex &tail_vertex = tail_vertices[pos];
        ui u = tail_vertex.u;
        assert(mapped_q[u] == -1);

        ui remaining_budget = threshold - cost;
        ui max_delta = std::min((ui)tail_vertex.buckets.size() - 1, remaining_budget);
        for (ui missing_delta = 0; missing_delta <= max_delta; ++missing_delta) {
            const vector<ui> &bucket = tail_vertex.buckets[missing_delta];
            for (ui v : bucket) {
                if (mapped_g[v] != -1) {
                    continue;
                }

                mapped_q[u] = (int)v;
                mapped_g[v] = (int)u;
                part_M.push_back({ u, v });

                enumerateTerminalTail(pos + 1, cost + missing_delta, tail_vertices);

                part_M.pop_back();
                mapped_g[v] = -1;
                mapped_q[u] = -1;

                if (outputLimitReached()) {
                    return;
                }
            }
        }
    }
    // ========================================================================

    // =====================================================
    // Procedure DFS(M_part, cost, X)
    //
    // cost:  current cost of partial match M_part
    // X:     the set of excluded edges (u, ua)
    // =====================================================
    void dfs(ui cost)
    {
        // printf("part_M.size() = %zu, cost = %u\n", part_M.size(), cost);
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }

        assert(part_M.size() <= qn);
        assert(cost <= threshold);

        if (part_M.size() == qn) {
            stats.recursion_calls++;
            stats.result_count++;
            noteOutputLimitIfReached();
#ifndef NDEBUG
            results_ptr->push_back(part_M);
#endif
            return;
        }

        stats.recursion_calls++;

        DfsBuffer &buf = dfsBufferForDepth(part_M.size());
        Timer t_buffer_clear;
        buf.clearLocal();
        stats.dfs_buffer_clear_time += t_buffer_clear.elapsed();
        SupportUndoScope support_undo_scope(*this, buf.local_support_snapshots);
        vector<ui> &candidate_vertices = buf.candidate_vertices;
        vector<ui> &candidate_anchor_counts = buf.candidate_anchor_counts;
        vector<ActiveEdge> &top_edges = buf.top_edges;
        ui current_cost = cost;
        double selected_component_support_sum = std::numeric_limits<double>::max();
#if !CDE_EDGE_IE_FIXED_ORDER
        bool selected_covered_component = false;
#endif
        bool has_zero_support_component = false;
        const vector<char> *terminal_skip_vertices = nullptr;

        Timer t_terminal;
        Timer t_terminal_scan;
        TerminalScan terminal_scan = markTerminalVertices(buf);
        stats.terminal_scan_time += t_terminal_scan.elapsed();
        stats.terminal_scan_calls++;

        if (terminal_scan.allRemainingTerminal()) {
            Timer t_terminal_build;
            stats.terminal_tail_build_calls++;
            if (!buildTerminalTailVertices(buf.terminal_vertices,
                current_cost, buf)) {
                stats.terminal_tail_build_time += t_terminal_build.elapsed();
                stats.terminal_time += t_terminal.elapsed();
                recordTerminalPrune();
                return;
            }
            stats.terminal_tail_build_time += t_terminal_build.elapsed();

            Timer t_terminal_enum;
            enumerateTerminalTail(0, current_cost, buf.terminal_tail_vertices);
            stats.terminal_tail_enum_time += t_terminal_enum.elapsed();
            stats.terminal_time += t_terminal.elapsed();
            return;
        }

        if (terminal_scan.terminal_frontier_count > 0) {
            Timer t_terminal_build;
            stats.terminal_tail_build_calls++;
            if (!buildTerminalTailVertices(buf.active_terminal_vertices,
                current_cost, buf)) {
                stats.terminal_tail_build_time += t_terminal_build.elapsed();
                stats.terminal_time += t_terminal.elapsed();
                recordTerminalPrune();
                return;
            }
            stats.terminal_tail_build_time += t_terminal_build.elapsed();

            ui remaining_budget = threshold - current_cost;
            Timer t_min_delta;
            if (terminalMinDeltaSum(buf.terminal_tail_vertices,
                remaining_budget) > remaining_budget) {
                stats.terminal_min_delta_time += t_min_delta.elapsed();
                stats.terminal_time += t_terminal.elapsed();
                recordTerminalPrune();
                return;
            }
            stats.terminal_min_delta_time += t_min_delta.elapsed();

            if (terminal_scan.hasNonterminalFrontier()) {
                terminal_skip_vertices = &buf.terminal_skip;
                stats.terminal_delayed_vertices += terminal_scan.terminal_frontier_count;
            }
        }
        stats.terminal_time += t_terminal.elapsed();

        Timer t_frontier;
        ui max_branch_edges = threshold - current_cost + 1;
        {
            Timer t_select;
            if (!buf.branch_selector.collectTopActiveEdges(active_frontier,
                max_branch_edges, top_edges, &buf.edge_score_cache,
                terminal_skip_vertices)) {
                stats.frontier_time += t_frontier.elapsed();
                return;
            }
            stats.frontier_select_time += t_select.elapsed();
        }

#if !CDE_EDGE_IE_FIXED_ORDER
        bool all_selected_edges_have_zero_support = std::all_of(
            top_edges.begin(), top_edges.end(),
            [](const ActiveEdge &edge) {
                return edge.anchor_support == 0;
            });
        if (all_selected_edges_have_zero_support) {
            stats.frontier_time += t_frontier.elapsed();
            stats.prun_calls++;
            return;
        }
#endif

        {
            Timer t_component;
#if CDE_EDGE_IE_FIXED_ORDER
            (void)buf.branch_selector.restrictTopEdgesToCoveredComponent(top_edges,
                buf.edge_score_cache, buf.component_id, buf.component_frontiers,
                buf.component_edge_counts, buf.component_seen_counts,
                buf.component_support_sums, selected_component_support_sum,
                has_zero_support_component, terminal_skip_vertices);
#else
            selected_covered_component = buf.branch_selector.restrictTopEdgesToCoveredComponent(top_edges,
                buf.edge_score_cache, buf.component_id, buf.component_frontiers,
                buf.component_edge_counts, buf.component_seen_counts,
                buf.component_support_sums, selected_component_support_sum,
                has_zero_support_component, terminal_skip_vertices);
#endif
            stats.frontier_component_time += t_component.elapsed();
        }
        stats.frontier_time += t_frontier.elapsed();

#if !CDE_EDGE_IE_FIXED_ORDER
        if (has_zero_support_component ||
            (selected_covered_component && selected_component_support_sum == 0)) {
            stats.prun_calls++;
            return;
        }
#endif

        Timer t_branch;
        long long child_dfs_time = 0;

        ui first_branch_edge = 0;
        bool pruned_by_forced_zero = false;
#if !CDE_EDGE_IE_FIXED_ORDER
        for (; first_branch_edge < top_edges.size(); ++first_branch_edge) {
            const ActiveEdge &edge = top_edges[first_branch_edge];
            if (edge.anchor_support != 0) break;

            ui u = edge.u;
            ui ua = edge.anchor;
            assert(u < qn && ua < qn);
            assert(mapped_q[ua] != -1 && mapped_q[u] == -1);
            assert(!excluded_edges[u][ua] && frontier_pos[u] != -1);
            assert(q_matrix[u][ua]);

            Timer t_exclude_update;
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);
            stats.exclude_update_time += t_exclude_update.elapsed();

            if (current_cost > threshold) {
                stats.prun_calls++;
                pruned_by_forced_zero = true;
                break;
            }
        }
#endif

        for (ui edge_idx = first_branch_edge; !pruned_by_forced_zero && edge_idx < top_edges.size(); ++edge_idx) {
            const ActiveEdge &edge = top_edges[edge_idx];
            if (current_cost > threshold) break;

            ui u = edge.u;
            ui ua = edge.anchor;
            assert(u < qn && ua < qn);
            assert(mapped_q[ua] != -1 && mapped_q[u] == -1);
            assert(!excluded_edges[u][ua] && frontier_pos[u] != -1);
            assert(q_matrix[u][ua]);

            // calculate the candidate vertices of u supported by anchor ua
            candidate_vertices.clear();
            {
                Timer t_cal_edge_support;
                calEdgeSupport(u, ua, [&](ui v) {candidate_vertices.push_back(v);});
                stats.branch_cal_edge_support_time += t_cal_edge_support.elapsed();
            }
            // number of live anchors of u excluding anchor ua
            ui anchor_num = 0;
            {
                Timer t_count_anchors;
                anchor_num = countAnchors(u, ua, candidate_vertices, candidate_anchor_counts);
                stats.branch_count_anchors_time += t_count_anchors.elapsed();
            }

            ui remaining_budget = threshold - current_cost;
            buf.clearCandidateDeltaBuckets();
            {
                Timer t_delta_bucket;
                for (ui i = 0; i < candidate_vertices.size(); ++i) {
                    ui delta = anchor_num - candidate_anchor_counts[i];
                    if (delta <= remaining_budget) {
                        buf.addCandidateDelta(delta, candidate_vertices[i]);
                    }
                }
                stats.branch_delta_bucket_time += t_delta_bucket.elapsed();
            }

            // Include branch: use this frontier-anchor edge to add exactly one new query vertex.
            Timer t_candidate_loop;
            long long candidate_child_time_before = child_dfs_time;
            long long candidate_support_update_time = 0;
            for (ui delta = 0; delta <= remaining_budget; ++delta) {
                const vector<ui> &delta_bucket = buf.candidate_delta_buckets[delta];
                for (ui v : delta_bucket) {
                    ui next_cost = current_cost + delta;

                    mapped_q[u] = (int)v;
                    mapped_g[v] = (int)u;
                    part_M.push_back({ u, v });

                    long long support_update_before = stats.support_update_time;
                    updateFrontier(u, true);
                    candidate_support_update_time += stats.support_update_time - support_update_before;

                    Timer t_child;
                    dfs(next_cost);
                    child_dfs_time += t_child.elapsed();

                    part_M.pop_back();
                    mapped_g[v] = -1;
                    mapped_q[u] = -1;

                    support_update_before = stats.support_update_time;
                    updateFrontier(u, false);
                    candidate_support_update_time += stats.support_update_time - support_update_before;

                    if (outputLimitReached()) break;
                }
                if (outputLimitReached()) break;
            }
            {
                long long candidate_elapsed = t_candidate_loop.elapsed();
                long long candidate_excluded_time =
                    (child_dfs_time - candidate_child_time_before) +
                    candidate_support_update_time;
                stats.candidate_loop_time += candidate_elapsed > candidate_excluded_time
                    ? candidate_elapsed - candidate_excluded_time : 0;
            }

            if (outputLimitReached()) break;

            // Exclude branch: skip this edge, no recursion, just update cost and state.
            Timer t_exclude_update;
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);

            if (current_cost > threshold) {
                stats.exclude_update_time += t_exclude_update.elapsed();
                break;
            }

            for (ui v : candidate_vertices) {
                if (!excluded_cands[u].contains(v)) {
                    excluded_cands[u].insert(v);
                    buf.recordExcludedCands(u, v);
                }
            }
            markLiveAnchorSupportDirty(u);
            stats.exclude_update_time += t_exclude_update.elapsed();
        }

        {
            long long branch_elapsed = t_branch.elapsed();
            long long excluded_time = child_dfs_time;
            stats.branch_time += branch_elapsed > excluded_time
                ? branch_elapsed - excluded_time : 0;
        }

        buf.restoreLocalChanges(*this);
    }
    // ========================================================================
};

// ============================================================
// Top-level function: Approximate_CDE_EdgeIE
// ============================================================
void Approximate_CDE_EdgeIE(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
