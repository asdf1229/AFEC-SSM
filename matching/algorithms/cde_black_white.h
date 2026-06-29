#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#define CDE_BLACK_WHITE_TERMINAL_BUCKETS_DEFAULT 0
#define CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING 1
#define CDE_BLACK_WHITE_ENABLE_BRIDGE_FILTERING 1
#define CDE_BLACK_WHITE_FIXED_ORDER 0
#define CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY 1
#define CDE_BLACK_WHITE_USE_DYNAMIC_SEARCH 1
#define CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY_GAMMA 0.9

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
using namespace std;

// ============================================================================
// CDEBlackWhiteSolver Implementation
// ============================================================================
class CDEBlackWhiteSolver {
    enum VertexColor : unsigned char {
        COLOR_UNSELECTED = 0,
        COLOR_BLACK = 1,
        COLOR_WHITE = 2
    };

    enum EdgeState : unsigned char {
        EDGE_UNDECIDED = 0,
        EDGE_PRESENT = 1,
        EDGE_MISSING = 2
    };

    struct WhiteCandidateBuckets {
        vector<ui> candidates;
        ui feasible_count = 0;

        void clear()
        {
            candidates.clear();
            feasible_count = 0;
        }

        bool empty() const
        {
            return feasible_count == 0;
        }

        void reserve(size_t size)
        {
            candidates.reserve(size);
        }

        void addCandidate(ui candidate)
        {
            candidates.push_back(candidate);
            feasible_count++;
        }
    };

    struct BlackWhiteState {
        vector<int> mapped_q;
        vector<ui> used_data_vertices;
        vector<unsigned char> used_data_flag;
        vector<VertexColor> color;
        vector<EdgeState> edge_state;
        vector<WhiteCandidateBuckets> white;
        vector<pair<ui, ui>> part_M;
        ui selected_count = 0;
        ui white_count = 0;
    };

    struct BlackWhiteActiveEdge {
        ui u = 0;
        ui anchor = 0;
        double rank_support = std::numeric_limits<double>::max();
        ui live_anchor_count = 0;
        ui query_degree = 0;
    };

    enum BlackWhiteUndoKind : unsigned char {
        BW_UNDO_MAPPED_Q = 0,
        BW_UNDO_COLOR = 1,
        BW_UNDO_EDGE_STATE = 2,
        BW_UNDO_USED_DATA_SIZE = 3,
        BW_UNDO_PART_M_SIZE = 4,
        BW_UNDO_SELECTED_COUNT = 5,
        BW_UNDO_WHITE_COUNT = 6,
        BW_UNDO_WHITE_BUCKET = 7
    };

    struct BlackWhiteUndo {
        BlackWhiteUndoKind kind = BW_UNDO_MAPPED_Q;
        ui u = 0;
        ui v = 0;
        int old_mapped_q = -1;
        VertexColor old_color = COLOR_UNSELECTED;
        EdgeState old_edge_uv = EDGE_UNDECIDED;
        EdgeState old_edge_vu = EDGE_UNDECIDED;
        size_t old_size = 0;
        ui old_count = 0;
        WhiteCandidateBuckets old_white;
    };

public:
    CDEBlackWhiteSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

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

        initBlackWhiteStaticColors();

#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
        initDataLabelDegreeIndex();
#endif
#if CDE_BLACK_WHITE_FIXED_ORDER
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

#if CDE_BLACK_WHITE_USE_DYNAMIC_SEARCH
        ui root = bw_static_root < qn ? bw_static_root : selectInitialRootBlackWhite();
        for (ui v0 : candidates[root]) {
            BlackWhiteState state;
            initBlackWhiteState(state);
            if (!commitRootBlack(state, root, v0)) {
                continue;
            }
            bwSearch(state, 0);
            if (outputLimitReached()) break;
        }
#else
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
#endif

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
        ui filter_candidate_count = 0;
        // search breakdown
        long long dfs_time = 0;
        long long frontier_time = 0; // building and ordering U_frontier
        long long frontier_select_time = 0;
        long long frontier_component_time = 0;
        long long frontier_sort_time = 0;
        long long frontier_score_time = 0;
        long long frontier_score_live_anchor_time = 0;
        long long frontier_score_live_candidate_time = 0;
        long long frontier_score_query_degree_time = 0;
        long long frontier_score_anchor_support_time = 0;
        long long branch_time = 0;   // candidate enumeration & matching in dfs
        long long branch_cal_edge_support_time = 0;
        long long branch_count_anchors_time = 0;
        long long support_update_time = 0;
        long long candidate_loop_time = 0;
        long long exclude_update_time = 0;
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
#ifndef NDEBUG
        unsigned long long terminal_tail_calls = 0;
        unsigned long long terminal_prune_calls = 0;
        unsigned long long terminal_delayed_vertices = 0;
        unsigned long long terminal_bucket_candidate_checks = 0;
#endif
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        long long search_accounted_time = stats.frontier_time + stats.branch_time;
        long long search_other_time = stats.dfs_time > search_accounted_time
            ? stats.dfs_time - search_accounted_time : 0;

        printf("\n--- CDE-Black-White Time Analysis ---\n");
#ifdef NDEBUG
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms\n", stats.filter_time / 1000.0);
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.dfs_time));
#else
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("    - NLF:           %.4lf ms (%.2f%% of Filter)\n", stats.filter_nlf_time / 1000.0, pct(stats.filter_nlf_time, stats.filter_time));
        printf("    - Bridge:        %.4lf ms (%.2f%% of Filter)\n", stats.filter_bridge_time / 1000.0, pct(stats.filter_bridge_time, stats.filter_time));
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
        printf("    - Spoke:         %.4lf ms (%.2f%% of Filter)\n", stats.filter_spoke_time / 1000.0, pct(stats.filter_spoke_time, stats.filter_time));
#endif
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
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
        printf("      - Anchor Support: %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_anchor_support_time / 1000.0, pct(stats.frontier_score_anchor_support_time, stats.frontier_score_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("    - branch_cal_edge_support_ms: %.4lf (%.2f%% of Branch)\n", stats.branch_cal_edge_support_time / 1000.0, pct(stats.branch_cal_edge_support_time, stats.branch_time));
        printf("    - branch_count_anchors_ms:    %.4lf (%.2f%% of Branch)\n", stats.branch_count_anchors_time / 1000.0, pct(stats.branch_count_anchors_time, stats.branch_time));
        printf("    - support_update_ms:          %.4lf (%.2f%% of Branch)\n", stats.support_update_time / 1000.0, pct(stats.support_update_time, stats.branch_time));
        printf("    - candidate_loop_ms:          %.4lf (%.2f%% of Branch)\n", stats.candidate_loop_time / 1000.0, pct(stats.candidate_loop_time, stats.branch_time));
        printf("    - exclude_update_ms:          %.4lf (%.2f%% of Branch)\n", stats.exclude_update_time / 1000.0, pct(stats.exclude_update_time, stats.branch_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.dfs_time));
#endif
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
#ifndef NDEBUG
        printf("Terminal Tail Calls: %llu\n", stats.terminal_tail_calls);
        printf("Terminal Prunes:     %llu\n", stats.terminal_prune_calls);
        printf("Terminal Delayed:    %llu\n", stats.terminal_delayed_vertices);
        printf("Terminal Cand Checks:%llu\n", stats.terminal_bucket_candidate_checks);
#endif
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
#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
    struct DataLabelDegreeCount {
        LabelID label = 0;
        ui count = 0;
    };
#endif

    struct ActiveEdge {
        ui u = 0;       // unmatched endpoint
        ui anchor = 0;  // matched endpoint
        ui anchor_support = std::numeric_limits<ui>::max();
        double rank_support = std::numeric_limits<double>::max();
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
#if CDE_BLACK_WHITE_FIXED_ORDER
    vector<vector<ui>> fixed_edge_priority;
#endif
#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
    vector<vector<DataLabelDegreeCount>> data_label_degrees;
#endif

    vector<MyBitset> candidates;
    ui bw_static_root = 0;
    vector<VertexColor> bw_static_color;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<vector<char>> excluded_edges;
    vector<MyBitset> excluded_cands;
    vector<pair<ui, ui>> part_M;
    vector<BlackWhiteUndo> black_white_undo;

    vector<ui> anchor_count;
    vector<vector<ui>> support;
    vector<vector<char>> support_dirty;
    vector<SupportSnapshot> *active_support_snapshots = nullptr;
    vector<vector<BlackWhiteActiveEdge>> bw_top_edges_buffer_by_depth;
    vector<vector<ui>> bw_white_neighbors_buffer_by_depth;
    vector<int> frontier_pos;
    vector<ui> active_frontier;
    vector<ui> data_vertex_mark;
    vector<ui> data_vertex_mark_pos;
    bool terminal_buckets_enabled = true;
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
#if CDE_BLACK_WHITE_FIXED_ORDER
        fixed_edge_priority.clear();
#endif
#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
        data_label_degrees.clear();
#endif

        anchor_count.assign(qn, 0);
        support.assign(qn, vector<ui>(qn, 0));
        support_dirty.assign(qn, vector<char>(qn, 1));
        active_support_snapshots = nullptr;
        bw_top_edges_buffer_by_depth.clear();
        bw_top_edges_buffer_by_depth.resize((size_t)qn + 1);
        bw_white_neighbors_buffer_by_depth.clear();
        bw_white_neighbors_buffer_by_depth.resize((size_t)qn + 1);
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
        data_vertex_mark.assign(gn, 0);
        data_vertex_mark_pos.assign(gn, 0);
        data_vertex_mark_token = 0;
        black_white_undo.clear();
        stats = TimeStats();
        terminal_buckets_enabled = CDE_BLACK_WHITE_TERMINAL_BUCKETS_DEFAULT != 0;
        bw_static_root = 0;
        bw_static_color.clear();
    }

    // ========================================================================
    // Filtering
    // ========================================================================
    struct CandidateFilter {
    private:
        CDEBlackWhiteSolver &solver;

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
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
        queue<ui>           pending_spokes;
        vector<char>        queued_spoke;
#endif

    public:
        explicit CandidateFilter(CDEBlackWhiteSolver &solver)
            : solver(solver),
            spoke_adj(solver.qn, vector<ui>()),
            match_right(solver.max_g_deg, -1),
            seen_right(solver.max_g_deg, 0),
            left_is_bridge(solver.qn, 0)
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
            , queued_spoke(solver.qn, 0)
#endif
        {}

        bool run()
        {
            if (!timed(&CDEBlackWhiteSolver::TimeStats::filter_bridge_time, [&] {
                buildBridgeIndex();
                return true;
            })) return false;

            if (!timed(&CDEBlackWhiteSolver::TimeStats::filter_nlf_time, [&] {
                return filterNLF();
            })) return false;

#if CDE_BLACK_WHITE_ENABLE_BRIDGE_FILTERING
            if (!timed(&CDEBlackWhiteSolver::TimeStats::filter_bridge_time, [&] {
                return initBridgeSupports() && propagateFilterClosure();
            })) return false;
#endif

#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
            if (!timed(&CDEBlackWhiteSolver::TimeStats::filter_spoke_time, [&] {
                enqueueAllSpokeVertices();
                return propagateFilterClosure();
            })) return false;
#endif

            updateCandidateCount();
            return true;
        }

    private:
        template <typename Fn>
        bool timed(long long CDEBlackWhiteSolver::TimeStats::*field, Fn &&fn)
        {
#ifndef NDEBUG
            Timer t;
#endif
            bool ok = fn();
#ifndef NDEBUG
            solver.stats.*field += t.elapsed();
#endif
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

#if CDE_BLACK_WHITE_ENABLE_BRIDGE_FILTERING
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
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
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
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
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

#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
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

#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
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

    double initialTopkDecayRankSupport(ui u, ui anchor)
    {
        if (u >= qn || anchor >= qn || mapped_q[anchor] == -1) {
            return 0.0;
        }

        LabelID label = query_graph->getVertexLabel(u);
        ui label_degree = dataLabelDegree((ui)mapped_q[anchor], label);
        ui candidate_count = (ui)candidates[u].size();
        return (double)std::min(candidate_count, label_degree);
    }
#endif

#if CDE_BLACK_WHITE_FIXED_ORDER
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
        CDEBlackWhiteSolver &solver;

    public:
        explicit BranchSelector(CDEBlackWhiteSolver &solver) : solver(solver) {}

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

#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
            selectTopActiveEdgesWithDecay(max_count, top_edges);
#else
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
#endif
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
#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
                    for (const ActiveEdge &edge : edges) {
                        component_support_sums[id] += edge.rank_support;
                    }
#elif !CDE_BLACK_WHITE_FIXED_ORDER
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

#if !CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY && !CDE_BLACK_WHITE_FIXED_ORDER
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
            ui count = 0;
            for (ui anchor : solver.q_neighbors[u]) {
                if (solver.mapped_q[anchor] != -1 && !solver.excluded_edges[u][anchor]) {
                    count++;
                }
            }
            return count;
        }

        ui countEdgeSupport(ui u, ui anchor) const
        {
            if (u >= solver.qn || anchor >= solver.qn ||
                solver.mapped_q[anchor] == -1 || solver.excluded_edges[u][anchor]) {
                return 0;
            }

            if (solver.support_dirty[u][anchor]) {
#ifndef NDEBUG
                Timer t_score;
#endif
                solver.recordSupportSnapshot(u, anchor);
                solver.support[u][anchor] = solver.calEdgeSupport(u, anchor, [](ui) {});
                solver.support_dirty[u][anchor] = 0;
#ifndef NDEBUG
                long long elapsed = t_score.elapsed();
                solver.stats.frontier_score_anchor_support_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
#endif
            }
#ifndef NDEBUG
            else {
                Timer t_score;
                ui exact_count = solver.calEdgeSupport(u, anchor, [](ui) {});
                long long elapsed = t_score.elapsed();
                solver.stats.frontier_score_anchor_support_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
                assert(solver.support[u][anchor] >= exact_count);
            }
#endif
            return solver.support[u][anchor];
        }

        bool isBetterActiveEdge(const ActiveEdge &lhs, const ActiveEdge &rhs) const
        {
#if CDE_BLACK_WHITE_FIXED_ORDER
            ui lhs_priority = solver.fixed_edge_priority[lhs.u][lhs.anchor];
            ui rhs_priority = solver.fixed_edge_priority[rhs.u][rhs.anchor];
            if (lhs_priority != rhs_priority) {
                return lhs_priority < rhs_priority;
            }
            if (lhs.u != rhs.u) {
                return lhs.u < rhs.u;
            }
            return lhs.anchor < rhs.anchor;
#elif CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
            double lhs_scaled =
                lhs.rank_support * (double)std::max((ui)1, rhs.live_anchor_count);
            double rhs_scaled =
                rhs.rank_support * (double)std::max((ui)1, lhs.live_anchor_count);
            double scale = std::max(1.0,
                std::max(std::fabs(lhs_scaled), std::fabs(rhs_scaled)));
            if (std::fabs(lhs_scaled - rhs_scaled) > 1e-12 * scale) {
                return lhs_scaled < rhs_scaled;
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

#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
        void selectTopActiveEdgesWithDecay(ui max_count, vector<ActiveEdge> &top_edges) const
        {
            size_t selected_limit = top_edges.size();
            if ((size_t)max_count < selected_limit) {
                selected_limit = (size_t)max_count;
            }

            const double gamma = (double)CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY_GAMMA;
            for (size_t selected_idx = 0; selected_idx < selected_limit; ++selected_idx) {
                size_t best_idx = selected_idx;
                for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
                    if (isBetterActiveEdge(top_edges[i], top_edges[best_idx])) {
                        best_idx = i;
                    }
                }

                if (best_idx != selected_idx) {
                    std::swap(top_edges[selected_idx], top_edges[best_idx]);
                }

                const ActiveEdge &selected = top_edges[selected_idx];
                double candidate_count = (double)solver.candidates[selected.u].size();
                if (candidate_count <= 0.0) {
                    continue;
                }

                double factor = 1.0 - gamma * selected.rank_support / candidate_count;
                if (factor < 0.0) {
                    factor = 0.0;
                }
                else if (factor > 1.0) {
                    factor = 1.0;
                }

                if (factor == 1.0) {
                    continue;
                }

                for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
                    if (top_edges[i].u == selected.u) {
                        top_edges[i].rank_support *= factor;
                    }
                }
            }

            top_edges.resize(selected_limit);
        }
#endif

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
#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
                edge.rank_support = solver.initialTopkDecayRankSupport(u, anchor);
#elif !CDE_BLACK_WHITE_FIXED_ORDER
                edge.anchor_support = countEdgeSupport(u, anchor);
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

    ui markDataVertexNeighbors(ui v)
    {
        ui token = nextDataVertexMarkToken();
        if (v >= gn) {
            return token;
        }

        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors(v, deg);
        for (ui i = 0; i < deg; ++i) {
            data_vertex_mark[nbrs[i]] = token;
        }
        return token;
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
        CDEBlackWhiteSolver &solver;
        vector<SupportSnapshot> *previous_snapshots;
        vector<SupportSnapshot> &snapshots;

        SupportUndoScope(CDEBlackWhiteSolver &solver, vector<SupportSnapshot> &snapshots)
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
        EdgeScoreCache edge_score_cache;
        BranchSelector branch_selector;
        vector<SupportSnapshot> local_support_snapshots;
        explicit DfsBuffer(CDEBlackWhiteSolver &solver)
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
            local_excluded_edges.reserve((size_t)threshold + 1);
            size_t cand_size = std::min((size_t)gn, ((size_t)threshold + 1) * (size_t)max_g_deg);
            local_excluded_cands.reserve(cand_size);
        }

        void clearLocal()
        {
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
            local_excluded_edges.clear();
            local_excluded_cands.clear();
            local_support_snapshots.clear();
        }

        void recordExcludedEdge(ui u, ui anchor)
        {
            local_excluded_edges.push_back({ u, anchor });
        }

        void recordExcludedCands(ui u, ui v)
        {
            local_excluded_cands.push_back({ u, v });
        }

        void restoreLocalChanges(CDEBlackWhiteSolver &solver)
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
    // Dynamic black/white search
    // ========================================================================
    void initBlackWhiteState(BlackWhiteState &state) const
    {
        state.mapped_q.assign(qn, -1);
        state.used_data_vertices.clear();
        state.used_data_vertices.reserve(qn);
        state.used_data_flag.assign(gn, 0);
        state.color.assign(qn, COLOR_UNSELECTED);
        state.edge_state.assign((size_t)qn * qn, EDGE_UNDECIDED);
        state.white.clear();
        state.white.resize(qn);
        state.part_M.clear();
        state.part_M.reserve(qn);
        state.selected_count = 0;
        state.white_count = 0;
    }

    size_t blackWhiteEdgeIndex(ui u, ui v) const
    {
        return (size_t)u * qn + v;
    }

    EdgeState getBlackWhiteEdgeState(const BlackWhiteState &state, ui u, ui v) const
    {
        return state.edge_state[blackWhiteEdgeIndex(u, v)];
    }

    void setBlackWhiteEdgeStateRaw(BlackWhiteState &state, ui u, ui v,
        EdgeState edge_state_value) const
    {
        state.edge_state[blackWhiteEdgeIndex(u, v)] = edge_state_value;
    }

    vector<BlackWhiteActiveEdge> &blackWhiteTopEdgesBuffer(ui depth)
    {
        if (bw_top_edges_buffer_by_depth.size() <= depth) {
            bw_top_edges_buffer_by_depth.resize((size_t)depth + 1);
        }
        vector<BlackWhiteActiveEdge> &buffer = bw_top_edges_buffer_by_depth[depth];
        buffer.clear();
        if (buffer.capacity() < qn) {
            buffer.reserve(qn);
        }
        return buffer;
    }

    vector<ui> &blackWhiteWhiteNeighborsBuffer(ui depth)
    {
        if (bw_white_neighbors_buffer_by_depth.size() <= depth) {
            bw_white_neighbors_buffer_by_depth.resize((size_t)depth + 1);
        }
        vector<ui> &buffer = bw_white_neighbors_buffer_by_depth[depth];
        buffer.clear();
        if (buffer.capacity() < qn) {
            buffer.reserve(qn);
        }
        return buffer;
    }

    ui selectInitialRootBlackWhite()
    {
        ui root = 0;
        for (ui u = 1; u < qn; ++u) {
            if (candidates[u].size() < candidates[root].size() ||
                (candidates[u].size() == candidates[root].size() &&
                    q_degree[u] > q_degree[root])) {
                root = u;
            }
        }
        return root;
    }

    void initBlackWhiteStaticColors()
    {
        bw_static_root = selectInitialRootBlackWhite();
        bw_static_color.assign(qn, COLOR_WHITE);
        if (bw_static_root < qn) {
            bw_static_color[bw_static_root] = COLOR_BLACK;
        }
    }

    bool shouldPreferStaticWhite(ui u) const
    {
        return u < bw_static_color.size() && bw_static_color[u] == COLOR_WHITE;
    }

    bool commitRootBlack(BlackWhiteState &state, ui root, ui v) const
    {
        if (root >= qn || v >= gn || !candidates[root].contains(v)) {
            return false;
        }
        state.color[root] = COLOR_BLACK;
        state.mapped_q[root] = (int)v;
        state.used_data_vertices.push_back(v);
        state.used_data_flag[v] = 1;
        state.part_M.push_back({ root, v });
        state.selected_count = 1;
        return true;
    }

    bool isDataVertexUsed(const BlackWhiteState &state, ui v) const
    {
        return v < state.used_data_flag.size() && state.used_data_flag[v] != 0;
    }

    bool isSelected(const BlackWhiteState &state, ui u) const
    {
        return u < state.color.size() && state.color[u] != COLOR_UNSELECTED;
    }

    bool isBlack(const BlackWhiteState &state, ui u) const
    {
        return u < state.color.size() && state.color[u] == COLOR_BLACK;
    }

    bool isWhite(const BlackWhiteState &state, ui u) const
    {
        return u < state.color.size() && state.color[u] == COLOR_WHITE;
    }

    size_t markBlackWhiteState() const
    {
        return black_white_undo.size();
    }

    void rollbackBlackWhiteState(BlackWhiteState &state, size_t mark)
    {
        while (black_white_undo.size() > mark) {
            BlackWhiteUndo undo = std::move(black_white_undo.back());
            black_white_undo.pop_back();

            switch (undo.kind) {
            case BW_UNDO_MAPPED_Q:
                state.mapped_q[undo.u] = undo.old_mapped_q;
                break;
            case BW_UNDO_COLOR:
                state.color[undo.u] = undo.old_color;
                break;
            case BW_UNDO_EDGE_STATE:
                setBlackWhiteEdgeStateRaw(state, undo.u, undo.v, undo.old_edge_uv);
                setBlackWhiteEdgeStateRaw(state, undo.v, undo.u, undo.old_edge_vu);
                break;
            case BW_UNDO_USED_DATA_SIZE:
                for (size_t i = undo.old_size; i < state.used_data_vertices.size(); ++i) {
                    ui used_v = state.used_data_vertices[i];
                    if (used_v < state.used_data_flag.size()) {
                        state.used_data_flag[used_v] = 0;
                    }
                }
                state.used_data_vertices.resize(undo.old_size);
                break;
            case BW_UNDO_PART_M_SIZE:
                state.part_M.resize(undo.old_size);
                break;
            case BW_UNDO_SELECTED_COUNT:
                state.selected_count = undo.old_count;
                break;
            case BW_UNDO_WHITE_COUNT:
                state.white_count = undo.old_count;
                break;
            case BW_UNDO_WHITE_BUCKET:
                state.white[undo.u] = std::move(undo.old_white);
                break;
            }
        }
    }

    void setBlackWhiteMappedQ(BlackWhiteState &state, ui u, int value)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_MAPPED_Q;
        undo.u = u;
        undo.old_mapped_q = state.mapped_q[u];
        black_white_undo.push_back(std::move(undo));
        state.mapped_q[u] = value;
    }

    void setBlackWhiteColor(BlackWhiteState &state, ui u, VertexColor value)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_COLOR;
        undo.u = u;
        undo.old_color = state.color[u];
        black_white_undo.push_back(std::move(undo));
        state.color[u] = value;
    }

    void setBlackWhiteEdgeState(BlackWhiteState &state, ui u, ui v,
        EdgeState edge_state_value)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_EDGE_STATE;
        undo.u = u;
        undo.v = v;
        undo.old_edge_uv = getBlackWhiteEdgeState(state, u, v);
        undo.old_edge_vu = getBlackWhiteEdgeState(state, v, u);
        black_white_undo.push_back(std::move(undo));
        setBlackWhiteEdgeStateRaw(state, u, v, edge_state_value);
        setBlackWhiteEdgeStateRaw(state, v, u, edge_state_value);
    }

    void pushBlackWhiteUsedDataVertex(BlackWhiteState &state, ui v)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_USED_DATA_SIZE;
        undo.old_size = state.used_data_vertices.size();
        black_white_undo.push_back(std::move(undo));
        state.used_data_vertices.push_back(v);
        state.used_data_flag[v] = 1;
    }

    void pushBlackWhitePartM(BlackWhiteState &state, ui u, ui v)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_PART_M_SIZE;
        undo.old_size = state.part_M.size();
        black_white_undo.push_back(std::move(undo));
        state.part_M.push_back({ u, v });
    }

    void setBlackWhiteSelectedCount(BlackWhiteState &state, ui value)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_SELECTED_COUNT;
        undo.old_count = state.selected_count;
        black_white_undo.push_back(std::move(undo));
        state.selected_count = value;
    }

    void setBlackWhiteWhiteCount(BlackWhiteState &state, ui value)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_WHITE_COUNT;
        undo.old_count = state.white_count;
        black_white_undo.push_back(std::move(undo));
        state.white_count = value;
    }

    void replaceBlackWhiteBucket(BlackWhiteState &state, ui u,
        WhiteCandidateBuckets &&bucket)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_WHITE_BUCKET;
        undo.u = u;
        undo.old_white = std::move(state.white[u]);
        black_white_undo.push_back(std::move(undo));
        state.white[u] = std::move(bucket);
    }

    bool computeBlackNeighborDelta(const BlackWhiteState &state, ui u, ui v,
        ui cost, ui &delta) const
    {
        delta = 0;
        for (ui neighbor : q_neighbors[u]) {
            if (!isBlack(state, neighbor)) {
                continue;
            }

            ui mapped_neighbor = (ui)state.mapped_q[neighbor];
            bool adjacent = data_graph->hasEdge(v, mapped_neighbor);
            EdgeState state_uv = getBlackWhiteEdgeState(state, u, neighbor);
            if (state_uv == EDGE_PRESENT) {
                if (!adjacent) return false;
            }
            else if (state_uv == EDGE_MISSING) {
                if (adjacent) return false;
            }
            else if (!adjacent) {
                delta++;
                if (cost + delta > threshold) {
                    return false;
                }
            }
        }
        return true;
    }

    bool rebuildWhiteCandidatesForCurrentState(BlackWhiteState &state,
        ui white_u, ui cost)
    {
        const WhiteCandidateBuckets &white = state.white[white_u];
        WhiteCandidateBuckets next_white;
        next_white.reserve(white.candidates.size());

        for (ui candidate : white.candidates) {
            if (isDataVertexUsed(state, candidate)) {
                continue;
            }

            ui delta = 0;
            if (computeBlackNeighborDelta(state, white_u, candidate, cost, delta)) {
                next_white.addCandidate(candidate);
            }
        }

        bool has_candidate = !next_white.empty();
        replaceBlackWhiteBucket(state, white_u, std::move(next_white));
        return has_candidate;
    }

    bool buildWhiteCandidateBuckets(const BlackWhiteState &state, ui u, ui cost,
        WhiteCandidateBuckets &white_candidate_buckets) const
    {
        white_candidate_buckets.clear();
        if (cost > threshold || !isSelectedByBlackNeighbor(state, u)) {
            return false;
        }

        for (int candidate : candidates[u]) {
            ui v = (ui)candidate;
            if (isDataVertexUsed(state, v)) {
                continue;
            }

            ui delta = 0;
            if (!computeBlackNeighborDelta(state, u, v, cost, delta)) {
                continue;
            }

            white_candidate_buckets.addCandidate(v);
        }
        return !white_candidate_buckets.empty();
    }

    bool isSelectedByBlackNeighbor(const BlackWhiteState &state, ui u) const
    {
        for (ui neighbor : q_neighbors[u]) {
            if (isBlack(state, neighbor) &&
                getBlackWhiteEdgeState(state, u, neighbor) != EDGE_MISSING) {
                return true;
            }
        }
        return false;
    }

    void collectSelectedWhiteNeighbors(const BlackWhiteState &state, ui u,
        vector<ui> &white_neighbors) const
    {
        white_neighbors.clear();
        for (ui neighbor : q_neighbors[u]) {
            if (isWhite(state, neighbor)) {
                white_neighbors.push_back(neighbor);
            }
        }
    }

    bool chooseBlackWhite(const BlackWhiteState &state, ui u,
        const WhiteCandidateBuckets &white_candidate_buckets) const
    {
        if (!shouldPreferStaticWhite(u) || white_candidate_buckets.empty()) {
            return false;
        }

        for (ui neighbor : q_neighbors[u]) {
            if (isWhite(state, neighbor)) {
                return false;
            }
        }
        return true;
    }

    bool updateWhiteBucketForBlackNeighbor(BlackWhiteState &state, ui white_u,
        ui black_u, ui black_v, ui cost)
    {
        assert(isWhite(state, white_u));
        assert(isBlack(state, black_u));
        (void)black_u;
        (void)black_v;
        if (cost > threshold) {
            return false;
        }

        return rebuildWhiteCandidatesForCurrentState(state, white_u, cost);
    }

    bool commitNewBlackVertex(BlackWhiteState &state, ui cost, ui u, ui v,
        ui &next_cost)
    {
        if (u >= qn || v >= gn || !candidates[u].contains(v)) {
            return false;
        }
        if (state.color[u] != COLOR_UNSELECTED || isDataVertexUsed(state, v)) {
            return false;
        }

        ui delta = 0;
        if (!computeBlackNeighborDelta(state, u, v, cost, delta)) {
            return false;
        }
        next_cost = cost + delta;
        if (next_cost > threshold) {
            return false;
        }

        setBlackWhiteColor(state, u, COLOR_BLACK);
        setBlackWhiteMappedQ(state, u, (int)v);
        pushBlackWhiteUsedDataVertex(state, v);
        pushBlackWhitePartM(state, u, v);
        setBlackWhiteSelectedCount(state, state.selected_count + 1);

        for (ui neighbor : q_neighbors[u]) {
            if (!isBlack(state, neighbor)) {
                continue;
            }
            EdgeState state_uv = getBlackWhiteEdgeState(state, u, neighbor);
            if (state_uv != EDGE_UNDECIDED) {
                continue;
            }
            bool adjacent = data_graph->hasEdge(v, (ui)state.mapped_q[neighbor]);
            setBlackWhiteEdgeState(state, u, neighbor,
                adjacent ? EDGE_PRESENT : EDGE_MISSING);
        }

        for (ui neighbor : q_neighbors[u]) {
            if (!isWhite(state, neighbor)) {
                continue;
            }
            if (!updateWhiteBucketForBlackNeighbor(state, neighbor, u, v, next_cost)) {
                return false;
            }
        }
        return true;
    }

    bool commitMaterializedWhiteVertex(BlackWhiteState &state, ui cost, ui white_u,
        ui candidate, ui bucket_delta, ui &next_cost)
    {
        if (!isWhite(state, white_u) || candidate >= gn ||
            isDataVertexUsed(state, candidate) ||
            !candidates[white_u].contains(candidate)) {
            return false;
        }

        next_cost = cost + bucket_delta;
        if (next_cost > threshold) {
            return false;
        }

        setBlackWhiteColor(state, white_u, COLOR_BLACK);
        assert(state.white_count > 0);
        setBlackWhiteWhiteCount(state, state.white_count - 1);
        setBlackWhiteMappedQ(state, white_u, (int)candidate);
        pushBlackWhiteUsedDataVertex(state, candidate);
        pushBlackWhitePartM(state, white_u, candidate);

        for (ui neighbor : q_neighbors[white_u]) {
            if (!isSelected(state, neighbor)) {
                continue;
            }
            if (isWhite(state, neighbor)) {
                return false;
            }

            ui mapped_neighbor = (ui)state.mapped_q[neighbor];
            bool adjacent = data_graph->hasEdge(candidate, mapped_neighbor);
            EdgeState state_uv = getBlackWhiteEdgeState(state, white_u, neighbor);
            if (state_uv == EDGE_PRESENT) {
                if (!adjacent) return false;
            }
            else if (state_uv == EDGE_MISSING) {
                if (adjacent) return false;
            }
            else {
                setBlackWhiteEdgeState(state, white_u, neighbor,
                    adjacent ? EDGE_PRESENT : EDGE_MISSING);
            }
        }
        return true;
    }

    bool addWhiteVertexBranch(BlackWhiteState &state, ui cost, ui u)
    {
        if (state.color[u] != COLOR_UNSELECTED) {
            return false;
        }

        for (ui neighbor : q_neighbors[u]) {
            if (isWhite(state, neighbor)) {
                return false;
            }
        }

        WhiteCandidateBuckets white_candidate_buckets;
        if (!buildWhiteCandidateBuckets(state, u, cost, white_candidate_buckets)) {
            return false;
        }
        if (!chooseBlackWhite(state, u, white_candidate_buckets)) {
            return false;
        }

        size_t mark = markBlackWhiteState();
        setBlackWhiteColor(state, u, COLOR_WHITE);
        replaceBlackWhiteBucket(state, u, std::move(white_candidate_buckets));
        setBlackWhiteSelectedCount(state, state.selected_count + 1);
        setBlackWhiteWhiteCount(state, state.white_count + 1);
        bwSearch(state, cost);
        rollbackBlackWhiteState(state, mark);
        return true;
    }

    bool addBlackVertexBranches(BlackWhiteState &state, ui cost, ui u,
        ui required_anchor = std::numeric_limits<ui>::max())
    {
        bool emitted_branch = false;
        auto try_candidate = [&](ui candidate) -> bool {
            size_t mark = markBlackWhiteState();
            ui next_cost = cost;
            if (!commitNewBlackVertex(state, cost, u, candidate, next_cost)) {
                rollbackBlackWhiteState(state, mark);
                return false;
            }
            emitted_branch = true;
            bwSearch(state, next_cost);
            rollbackBlackWhiteState(state, mark);
            if (outputLimitReached()) {
                return true;
            }
            return false;
        };

        if (required_anchor < qn && isBlack(state, required_anchor) &&
            getBlackWhiteEdgeState(state, u, required_anchor) == EDGE_PRESENT) {
            ui mapped_anchor = (ui)state.mapped_q[required_anchor];
            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(mapped_anchor, degree);
            for (ui i = 0; i < degree; ++i) {
                ui candidate = neighbors[i];
                if (!candidates[u].contains(candidate)) {
                    continue;
                }
                if (try_candidate(candidate)) {
                    return true;
                }
            }
            return emitted_branch;
        }

        for (int candidate : candidates[u]) {
            if (try_candidate((ui)candidate)) {
                return true;
            }
        }
        return emitted_branch;
    }

    template <typename Continue>
    bool materializeWhiteVertexBranches(BlackWhiteState &state, ui cost,
        ui white_u, Continue continue_branch)
    {
        if (!isWhite(state, white_u)) {
            return false;
        }

        bool emitted_branch = false;
        const vector<ui> &white_candidates = state.white[white_u].candidates;
        for (ui candidate : white_candidates) {
            if (isDataVertexUsed(state, candidate)) {
                continue;
            }

            ui delta = 0;
            if (!computeBlackNeighborDelta(state, white_u, candidate, cost, delta)) {
                continue;
            }

            size_t mark = markBlackWhiteState();
            ui next_cost = cost;
            if (!commitMaterializedWhiteVertex(state, cost, white_u,
                candidate, delta, next_cost)) {
                rollbackBlackWhiteState(state, mark);
                continue;
            }
            if (continue_branch(state, next_cost)) {
                emitted_branch = true;
            }
            rollbackBlackWhiteState(state, mark);
            if (outputLimitReached()) {
                return true;
            }
        }
        return emitted_branch;
    }

    template <typename Continue>
    bool materializeWhiteSetBranches(BlackWhiteState &state, ui cost,
        const vector<ui> &white_vertices, size_t pos, Continue continue_branch)
    {
        if (pos == white_vertices.size()) {
            return continue_branch(state, cost);
        }

        ui white_u = white_vertices[pos];
        if (!isWhite(state, white_u)) {
            return materializeWhiteSetBranches(state, cost, white_vertices,
                pos + 1, continue_branch);
        }

        return materializeWhiteVertexBranches(state, cost, white_u,
            [&](BlackWhiteState &next_state, ui next_cost) -> bool {
                return materializeWhiteSetBranches(next_state, next_cost,
                    white_vertices, pos + 1, continue_branch);
            });
    }

    bool forcedIncludeBlackAnchor(BlackWhiteState &state, ui cost, ui u,
        ui anchor)
    {
        if (!isBlack(state, anchor) || state.color[u] != COLOR_UNSELECTED ||
            getBlackWhiteEdgeState(state, u, anchor) != EDGE_PRESENT) {
            return false;
        }

        vector<ui> &white_neighbors =
            blackWhiteWhiteNeighborsBuffer(state.selected_count);
        collectSelectedWhiteNeighbors(state, u, white_neighbors);
        bool emitted_white_branch = false;
        if (white_neighbors.empty()) {
            emitted_white_branch = addWhiteVertexBranch(state, cost, u);
        }
        else {
            emitted_white_branch = materializeWhiteSetBranches(state, cost,
                white_neighbors, 0,
                [&](BlackWhiteState &materialized_state, ui materialized_cost) -> bool {
                    return addWhiteVertexBranch(materialized_state,
                        materialized_cost, u);
                });
        }

        if (emitted_white_branch) {
            return true;
        }
        return addBlackVertexBranches(state, cost, u, anchor);
    }

    bool includeActiveEdgeBranch(BlackWhiteState &state, ui cost,
        const BlackWhiteActiveEdge &edge)
    {
        ui u = edge.u;
        ui anchor = edge.anchor;
        if (state.color[u] != COLOR_UNSELECTED ||
            !isSelected(state, anchor) ||
            getBlackWhiteEdgeState(state, u, anchor) != EDGE_UNDECIDED) {
            return false;
        }

        size_t mark = markBlackWhiteState();
        setBlackWhiteEdgeState(state, u, anchor, EDGE_PRESENT);

        bool emitted_branch = false;
        if (isBlack(state, anchor)) {
            emitted_branch = forcedIncludeBlackAnchor(state, cost, u, anchor);
        }
        else {
            emitted_branch = materializeWhiteVertexBranches(state, cost, anchor,
            [&](BlackWhiteState &materialized_state, ui materialized_cost) -> bool {
                return forcedIncludeBlackAnchor(materialized_state,
                    materialized_cost, u, anchor);
            });
        }
        rollbackBlackWhiteState(state, mark);
        return emitted_branch;
    }

    double estimateBlackAnchorSupport(const BlackWhiteState &state, ui u,
        ui anchor) const
    {
        if (!isBlack(state, anchor)) {
            return 0.0;
        }

        double support_count = 0.0;
        ui mapped_anchor = (ui)state.mapped_q[anchor];
        ui degree = 0;
        const ui *neighbors = data_graph->getVertexNeighbors(mapped_anchor, degree);
        for (ui i = 0; i < degree; ++i) {
            ui v = neighbors[i];
            if (!candidates[u].contains(v)) {
                continue;
            }
            if (isDataVertexUsed(state, v)) {
                continue;
            }
            support_count += 1.0;
        }
        return support_count;
    }

    double estimateWhiteAnchorSupport(const BlackWhiteState &state, ui anchor) const
    {
        if (!isWhite(state, anchor)) {
            return 0.0;
        }
        return (double)std::max((ui)1, state.white[anchor].feasible_count);
    }

    bool isBetterBlackWhiteActiveEdge(const BlackWhiteActiveEdge &lhs,
        const BlackWhiteActiveEdge &rhs) const
    {
        double lhs_scaled = lhs.rank_support *
            (double)std::max((ui)1, rhs.live_anchor_count);
        double rhs_scaled = rhs.rank_support *
            (double)std::max((ui)1, lhs.live_anchor_count);
        double scale = std::max(1.0,
            std::max(std::fabs(lhs_scaled), std::fabs(rhs_scaled)));
        if (std::fabs(lhs_scaled - rhs_scaled) > 1e-12 * scale) {
            return lhs_scaled < rhs_scaled;
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
    }

    void selectTopBlackWhiteActiveEdges(ui max_count,
        vector<BlackWhiteActiveEdge> &top_edges)
    {
        size_t selected_limit = top_edges.size();
        if ((size_t)max_count < selected_limit) {
            selected_limit = (size_t)max_count;
        }

#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
        const double gamma = (double)CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY_GAMMA;
        for (size_t selected_idx = 0; selected_idx < selected_limit; ++selected_idx) {
            size_t best_idx = selected_idx;
            for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
                if (isBetterBlackWhiteActiveEdge(top_edges[i], top_edges[best_idx])) {
                    best_idx = i;
                }
            }

            if (best_idx != selected_idx) {
                std::swap(top_edges[selected_idx], top_edges[best_idx]);
            }

            const BlackWhiteActiveEdge &selected = top_edges[selected_idx];
            double candidate_count = (double)candidates[selected.u].size();
            if (candidate_count <= 0.0) {
                continue;
            }

            double factor = 1.0 - gamma * selected.rank_support / candidate_count;
            if (factor < 0.0) {
                factor = 0.0;
            }
            else if (factor > 1.0) {
                factor = 1.0;
            }

            if (factor == 1.0) {
                continue;
            }

            for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
                if (top_edges[i].u == selected.u) {
                    top_edges[i].rank_support *= factor;
                }
            }
        }
#else
        auto better_edge = [&](const BlackWhiteActiveEdge &lhs,
            const BlackWhiteActiveEdge &rhs) {
            return isBetterBlackWhiteActiveEdge(lhs, rhs);
            };
        if (top_edges.size() > selected_limit) {
            partial_sort(top_edges.begin(), top_edges.begin() + selected_limit,
                top_edges.end(), better_edge);
        }
        else {
            sort(top_edges.begin(), top_edges.end(), better_edge);
        }
#endif
        top_edges.resize(selected_limit);
    }

    bool collectTopBlackWhiteActiveEdges(const BlackWhiteState &state,
        ui max_count, vector<BlackWhiteActiveEdge> &top_edges)
    {
        top_edges.clear();
        if (max_count == 0) {
            return false;
        }

        for (ui u = 0; u < qn; ++u) {
            if (state.color[u] != COLOR_UNSELECTED) {
                continue;
            }

            ui live_anchor_count = 0;
            for (ui anchor : q_neighbors[u]) {
                if (isSelected(state, anchor) &&
                    getBlackWhiteEdgeState(state, u, anchor) == EDGE_UNDECIDED) {
                    live_anchor_count++;
                }
            }
            if (live_anchor_count == 0) {
                continue;
            }

            for (ui anchor : q_neighbors[u]) {
                if (!isSelected(state, anchor) ||
                    getBlackWhiteEdgeState(state, u, anchor) != EDGE_UNDECIDED) {
                    continue;
                }

                BlackWhiteActiveEdge edge;
                edge.u = u;
                edge.anchor = anchor;
                edge.live_anchor_count = live_anchor_count;
                edge.query_degree = q_degree[u];
                if (isBlack(state, anchor)) {
                    edge.rank_support = estimateBlackAnchorSupport(state, u, anchor);
                }
                else {
                    edge.rank_support = estimateWhiteAnchorSupport(state, anchor);
                }
                top_edges.push_back(edge);
            }
        }

        if (top_edges.empty()) {
            return false;
        }
        selectTopBlackWhiteActiveEdges(max_count, top_edges);
        return true;
    }

    ui chooseWhiteVertexToMaterialize(const BlackWhiteState &state) const
    {
        ui chosen = qn;
        ui best_count = std::numeric_limits<ui>::max();
        for (ui u = 0; u < qn; ++u) {
            if (!isWhite(state, u)) {
                continue;
            }
            ui count = state.white[u].feasible_count;
            if (count < best_count || (count == best_count && u < chosen)) {
                chosen = u;
                best_count = count;
            }
        }
        return chosen;
    }

    bool buildBlackWhiteTerminalWhiteCandidateBuckets(
        const BlackWhiteState &state, ui white_u, ui cost,
        vector<vector<ui>> &buckets, ui &feasible_count, ui &min_delta) const
    {
        assert(cost <= threshold);
        if (!isWhite(state, white_u) || state.mapped_q[white_u] != -1) {
            return false;
        }

        for (ui neighbor : q_neighbors[white_u]) {
            if (isWhite(state, neighbor)) {
                return false;
            }
        }

        ui remaining_budget = threshold - cost;
        buckets.assign((size_t)remaining_budget + 1, vector<ui>());
        feasible_count = 0;
        min_delta = std::numeric_limits<ui>::max();

        const WhiteCandidateBuckets &white = state.white[white_u];
        for (ui candidate : white.candidates) {
            if (isDataVertexUsed(state, candidate)) {
                continue;
            }

            ui delta = 0;
            if (!computeBlackNeighborDelta(state, white_u, candidate, cost, delta)) {
                continue;
            }
            if (delta > remaining_budget) {
                continue;
            }

            buckets[delta].push_back(candidate);
            feasible_count++;
            if (delta < min_delta) {
                min_delta = delta;
            }
        }

        return feasible_count > 0;
    }

    bool buildBlackWhiteTerminalWhiteTailVertices(const BlackWhiteState &state,
        ui cost, vector<TerminalTailVertex> &tail_vertices) const
    {
        tail_vertices.clear();
        tail_vertices.reserve(state.white_count);

        for (ui u = 0; u < qn; ++u) {
            if (!isWhite(state, u)) {
                continue;
            }

            TerminalTailVertex tail_vertex;
            tail_vertex.u = u;
            if (!buildBlackWhiteTerminalWhiteCandidateBuckets(state, u, cost,
                tail_vertex.buckets, tail_vertex.feasible_count,
                tail_vertex.min_delta)) {
                return false;
            }
            tail_vertices.push_back(std::move(tail_vertex));
        }

        if (tail_vertices.size() != (size_t)state.white_count) {
            return false;
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

    void enumerateBlackWhiteTerminalWhiteTail(BlackWhiteState &state,
        size_t pos, ui cost, vector<TerminalTailVertex> &tail_vertices)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }
        assert(cost <= threshold);

        if (pos == tail_vertices.size()) {
            emitBlackWhiteResult(state);
            return;
        }

        TerminalTailVertex &tail_vertex = tail_vertices[pos];
        ui u = tail_vertex.u;
        assert(isWhite(state, u));
        assert(state.mapped_q[u] == -1);

        ui remaining_budget = threshold - cost;
        ui max_delta = std::min((ui)tail_vertex.buckets.size() - 1,
            remaining_budget);
        for (ui missing_delta = 0; missing_delta <= max_delta; ++missing_delta) {
            const vector<ui> &bucket = tail_vertex.buckets[missing_delta];
            for (ui v : bucket) {
                if (isDataVertexUsed(state, v)) {
                    continue;
                }

                state.mapped_q[u] = (int)v;
                state.used_data_vertices.push_back(v);
                state.used_data_flag[v] = 1;
                state.part_M.push_back({ u, v });

                enumerateBlackWhiteTerminalWhiteTail(state, pos + 1,
                    cost + missing_delta, tail_vertices);

                state.part_M.pop_back();
                state.used_data_flag[v] = 0;
                state.used_data_vertices.pop_back();
                state.mapped_q[u] = -1;

                if (outputLimitReached()) {
                    return;
                }
            }
        }
    }

    void emitBlackWhiteResult(const BlackWhiteState &state)
    {
        assert(state.part_M.size() == qn);
        stats.result_count++;
        noteOutputLimitIfReached();
#ifndef NDEBUG
        results_ptr->push_back(state.part_M);
#endif
    }

    void executeBlackWhiteTopEdges(BlackWhiteState &state, ui cost,
        const vector<BlackWhiteActiveEdge> &top_edges, size_t edge_idx)
    {
        if (outputLimitReached() || cost > threshold ||
            edge_idx >= top_edges.size()) {
            return;
        }

        const BlackWhiteActiveEdge &edge = top_edges[edge_idx];
        includeActiveEdgeBranch(state, cost, edge);
        if (outputLimitReached()) {
            return;
        }

        if (cost + 1 > threshold) {
            stats.prun_calls++;
            return;
        }

        size_t mark = markBlackWhiteState();
        if (state.color[edge.u] == COLOR_UNSELECTED &&
            isSelected(state, edge.anchor) &&
            getBlackWhiteEdgeState(state, edge.u, edge.anchor) == EDGE_UNDECIDED) {
            setBlackWhiteEdgeState(state, edge.u, edge.anchor, EDGE_MISSING);
            executeBlackWhiteTopEdges(state, cost + 1, top_edges, edge_idx + 1);
        }
        rollbackBlackWhiteState(state, mark);
    }

    void bwSearch(BlackWhiteState &state, ui cost)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }
        if (cost > threshold) {
            stats.prun_calls++;
            return;
        }

        stats.recursion_calls++;

        if (state.selected_count == qn) {
            if (state.white_count == 0) {
                emitBlackWhiteResult(state);
                return;
            }

            vector<TerminalTailVertex> tail_vertices;
            if (!buildBlackWhiteTerminalWhiteTailVertices(state, cost,
                tail_vertices)) {
                stats.prun_calls++;
                return;
            }
            enumerateBlackWhiteTerminalWhiteTail(state, 0, cost, tail_vertices);
            return;
        }

        vector<BlackWhiteActiveEdge> &top_edges =
            blackWhiteTopEdgesBuffer(state.selected_count);
        ui max_branch_edges = threshold - cost + 1;
        if (!collectTopBlackWhiteActiveEdges(state, max_branch_edges, top_edges)) {
            stats.prun_calls++;
            return;
        }

        executeBlackWhiteTopEdges(state, cost, top_edges, 0);
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

    ui terminalLiveAnchorCount(ui u) const
    {
        ui count = 0;
        for (ui anchor : q_neighbors[u]) {
            if (!excluded_edges[u][anchor] && mapped_q[anchor] != -1) {
                count++;
            }
        }
        return count;
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

    ui terminalMissingDelta(ui u, ui v, ui limit = std::numeric_limits<ui>::max()) const
    {
        ui delta = 0;
        for (ui anchor : q_neighbors[u]) {
            if (excluded_edges[u][anchor] || mapped_q[anchor] == -1) {
                continue;
            }
            if (!data_graph->hasEdge(v, (ui)mapped_q[anchor])) {
                delta++;
                if (delta > limit) {
                    return delta;
                }
            }
        }
        return delta;
    }

    template <typename Visitor>
    bool visitTerminalSupportedCandidates(ui u, Visitor visit)
    {
        if (++data_vertex_mark_token == 0) {
            std::fill(data_vertex_mark.begin(), data_vertex_mark.end(), 0);
            data_vertex_mark_token = 1;
        }
        ui token = data_vertex_mark_token;

        for (ui anchor : q_neighbors[u]) {
            if (excluded_edges[u][anchor] || mapped_q[anchor] == -1) {
                continue;
            }

            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (data_vertex_mark[v] == token) {
                    continue;
                }
                data_vertex_mark[v] = token;
#ifndef NDEBUG
                stats.terminal_bucket_candidate_checks++;
#endif

                if (!candidates[u].contains(v)) {
                    continue;
                }
                if (mapped_g[v] != -1) {
                    continue;
                }
                if (excluded_cands[u].contains(v)) {
                    continue;
                }
                if (!visit(v)) {
                    return false;
                }
            }
        }
        return true;
    }

    bool buildTerminalCandidateBuckets(ui u, ui cost, vector<vector<ui>> &buckets,
        ui &feasible_count, ui &min_delta)
    {
        assert(cost <= threshold);
        ui remaining_budget = threshold - cost;
        buckets.assign((size_t)remaining_budget + 1, vector<ui>());
        feasible_count = 0;
        min_delta = std::numeric_limits<ui>::max();
        ui live_anchor_count = terminalLiveAnchorCount(u);
        if (live_anchor_count == 0) {
            return false;
        }

        visitTerminalSupportedCandidates(u, [&](ui v) -> bool {
            ui missing_delta = terminalMissingDelta(u, v, remaining_budget);
            if (missing_delta > remaining_budget) {
                return true;
            }
            if (missing_delta >= live_anchor_count) {
                return true;
            }

            buckets[missing_delta].push_back(v);
            feasible_count++;
            if (missing_delta < min_delta) {
                min_delta = missing_delta;
            }
            return true;
        });

        return feasible_count > 0;
    }

    bool terminalVertexHasFeasibleCandidate(ui u, ui cost)
    {
        assert(cost <= threshold);
        ui remaining_budget = threshold - cost;
        ui live_anchor_count = terminalLiveAnchorCount(u);
        if (live_anchor_count == 0) {
            return false;
        }

        bool found = false;
        visitTerminalSupportedCandidates(u, [&](ui v) -> bool {
            ui missing_delta = terminalMissingDelta(u, v, remaining_budget);
            if (missing_delta <= remaining_budget &&
                missing_delta < live_anchor_count) {
                found = true;
                return false;
            }
            return true;
        });
        return found;
    }

    bool activeTerminalVerticesHaveCandidate(const vector<ui> &terminal_vertices,
        ui cost)
    {
        for (ui u : terminal_vertices) {
            if (!terminalVertexHasFeasibleCandidate(u, cost)) {
                return false;
            }
        }
        return true;
    }

    bool buildTerminalTailVertices(const vector<ui> &terminal_vertices, ui cost,
        vector<TerminalTailVertex> &tail_vertices)
    {
        tail_vertices.clear();
        tail_vertices.reserve(terminal_vertices.size());

        for (ui u : terminal_vertices) {
            TerminalTailVertex tail_vertex;
            tail_vertex.u = u;
            if (!buildTerminalCandidateBuckets(u, cost, tail_vertex.buckets,
                tail_vertex.feasible_count, tail_vertex.min_delta)) {
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

    void recordTerminalPrune()
    {
        stats.prun_calls++;
#ifndef NDEBUG
        stats.terminal_prune_calls++;
#endif
    }

    void enumerateTerminalTail(size_t pos, ui cost,
        vector<TerminalTailVertex> &tail_vertices)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }

#ifndef NDEBUG
        stats.terminal_tail_calls++;
#endif
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
        buf.clearLocal();
        SupportUndoScope support_undo_scope(*this, buf.local_support_snapshots);
        vector<ui> &candidate_vertices = buf.candidate_vertices;
        vector<ui> &candidate_anchor_counts = buf.candidate_anchor_counts;
        vector<ActiveEdge> &top_edges = buf.top_edges;
        ui current_cost = cost;
        double selected_component_support_sum = std::numeric_limits<double>::max();
#if !CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY && !CDE_BLACK_WHITE_FIXED_ORDER
        bool selected_covered_component = false;
#endif
        bool has_zero_support_component = false;
        const vector<char> *terminal_skip_vertices = nullptr;

        if (terminal_buckets_enabled) {
            TerminalScan terminal_scan = markTerminalVertices(buf);

            if (terminal_scan.allRemainingTerminal()) {
                if (!buildTerminalTailVertices(buf.terminal_vertices, current_cost,
                    buf.terminal_tail_vertices)) {
                    recordTerminalPrune();
                    return;
                }

                enumerateTerminalTail(0, current_cost, buf.terminal_tail_vertices);
                return;
            }

            if (terminal_scan.terminal_frontier_count > 0) {
                if (!activeTerminalVerticesHaveCandidate(buf.active_terminal_vertices,
                    current_cost)) {
                    recordTerminalPrune();
                    return;
                }

                if (terminal_scan.hasNonterminalFrontier()) {
                    terminal_skip_vertices = &buf.terminal_skip;
#ifndef NDEBUG
                    stats.terminal_delayed_vertices += terminal_scan.terminal_frontier_count;
#endif
                }
            }
        }

        Timer t_frontier;
        ui max_branch_edges = threshold - current_cost + 1;
        {
#ifndef NDEBUG
            Timer t_select;
#endif
            if (!buf.branch_selector.collectTopActiveEdges(active_frontier,
                max_branch_edges, top_edges, &buf.edge_score_cache,
                terminal_skip_vertices)) {
                stats.frontier_time += t_frontier.elapsed();
                return;
            }
#ifndef NDEBUG
            stats.frontier_select_time += t_select.elapsed();
#endif
        }

#if !CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY && !CDE_BLACK_WHITE_FIXED_ORDER
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
#ifndef NDEBUG
            Timer t_component;
#endif
#if CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY || CDE_BLACK_WHITE_FIXED_ORDER
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
#ifndef NDEBUG
            stats.frontier_component_time += t_component.elapsed();
#endif
        }
        stats.frontier_time += t_frontier.elapsed();

#if !CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY && !CDE_BLACK_WHITE_FIXED_ORDER
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
#if !CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY && !CDE_BLACK_WHITE_FIXED_ORDER
        for (; first_branch_edge < top_edges.size(); ++first_branch_edge) {
            const ActiveEdge &edge = top_edges[first_branch_edge];
            if (edge.anchor_support != 0) break;

            ui u = edge.u;
            ui ua = edge.anchor;
            assert(u < qn && ua < qn);
            assert(mapped_q[ua] != -1 && mapped_q[u] == -1);
            assert(!excluded_edges[u][ua] && frontier_pos[u] != -1);
            assert(q_matrix[u][ua]);

#ifndef NDEBUG
            Timer t_exclude_update;
#endif
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);
#ifndef NDEBUG
            stats.exclude_update_time += t_exclude_update.elapsed();
#endif

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
#ifndef NDEBUG
            {
                Timer t_cal_edge_support;
                calEdgeSupport(u, ua, [&](ui v) {candidate_vertices.push_back(v);});
                stats.branch_cal_edge_support_time += t_cal_edge_support.elapsed();
            }
#else
            calEdgeSupport(u, ua, [&](ui v) {candidate_vertices.push_back(v);});
#endif
            // number of live anchors of u excluding anchor ua
            ui anchor_num = 0;
#ifndef NDEBUG
            {
                Timer t_count_anchors;
                anchor_num = countAnchors(u, ua, candidate_vertices, candidate_anchor_counts);
                stats.branch_count_anchors_time += t_count_anchors.elapsed();
            }
#else
            anchor_num = countAnchors(u, ua, candidate_vertices, candidate_anchor_counts);
#endif

            // Include branch: use this frontier-anchor edge to add exactly one new query vertex.
#ifndef NDEBUG
            Timer t_candidate_loop;
            long long candidate_child_time_before = child_dfs_time;
            long long candidate_support_update_time = 0;
#endif
            for (ui i = 0; i < candidate_vertices.size(); ++i) {
                ui v = candidate_vertices[i];
                ui delta = anchor_num - candidate_anchor_counts[i];
                ui next_cost = current_cost + delta;
                if (next_cost > threshold) {
                    continue;
                }

                mapped_q[u] = (int)v;
                mapped_g[v] = (int)u;
                part_M.push_back({ u, v });

#ifndef NDEBUG
                long long support_update_before = stats.support_update_time;
#endif
                updateFrontier(u, true);
#ifndef NDEBUG
                candidate_support_update_time += stats.support_update_time - support_update_before;
#endif

                Timer t_child;
                dfs(next_cost);
                child_dfs_time += t_child.elapsed();

                part_M.pop_back();
                mapped_g[v] = -1;
                mapped_q[u] = -1;

#ifndef NDEBUG
                support_update_before = stats.support_update_time;
#endif
                updateFrontier(u, false);
#ifndef NDEBUG
                candidate_support_update_time += stats.support_update_time - support_update_before;
#endif

                if (outputLimitReached()) break;
            }
#ifndef NDEBUG
            {
                long long candidate_elapsed = t_candidate_loop.elapsed();
                long long candidate_excluded_time =
                    (child_dfs_time - candidate_child_time_before) +
                    candidate_support_update_time;
                stats.candidate_loop_time += candidate_elapsed > candidate_excluded_time
                    ? candidate_elapsed - candidate_excluded_time : 0;
            }
#endif

            if (outputLimitReached()) break;

            // Exclude branch: skip this edge, no recursion, just update cost and state.
#ifndef NDEBUG
            Timer t_exclude_update;
#endif
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);

            if (current_cost > threshold) {
#ifndef NDEBUG
                stats.exclude_update_time += t_exclude_update.elapsed();
#endif
                break;
            }

            for (ui v : candidate_vertices) {
                if (!excluded_cands[u].contains(v)) {
                    excluded_cands[u].insert(v);
                    buf.recordExcludedCands(u, v);
                }
            }
#ifndef NDEBUG
            stats.exclude_update_time += t_exclude_update.elapsed();
#endif
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
// Top-level function: Approximate_CDE_BlackWhite
// ============================================================
void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    CDEBlackWhiteSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
