#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#define CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING 1
#define CDE_BLACK_WHITE_ENABLE_BRIDGE_FILTERING 1
#define CDE_BLACK_WHITE_FIXED_ORDER 0
#define CDE_BLACK_WHITE_STATIC_COLOR 1
#define CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY 1
#define CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY_GAMMA 0.9
#ifndef CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SOURCE
#define CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SOURCE 1
#endif
#ifndef CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_ANCHOR_BRANCH
#define CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_ANCHOR_BRANCH 1
#endif
#ifndef CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SUPPORT
#define CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SUPPORT 1
#endif

#if CDE_BLACK_WHITE_FIXED_ORDER && CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
#error "CDE_BLACK_WHITE_FIXED_ORDER and CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY are mutually exclusive."
#endif

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
#include <iterator>
#include <unordered_map>
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
        size_t begin = 0;
        ui count = 0;
        ui feasible_count = 0;

        void clear()
        {
            begin = 0;
            count = 0;
            feasible_count = 0;
        }

        bool empty() const
        {
            return feasible_count == 0;
        }
    };

    struct BlackWhiteState {
        vector<int> mapped_q;
        vector<ui> used_data_vertices;
        vector<unsigned char> used_data_flag;
        vector<VertexColor> color;
        vector<EdgeState> edge_state;
        vector<WhiteCandidateBuckets> white;
        vector<ui> white_candidate_pool;
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

    struct CandidateAdjRange {
        size_t begin = 0;
        ui len = 0;
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

        // init q_neighbors
        q_neighbors.assign(qn, vector<ui>());
        q_degree.assign(qn, 0);
        for (ui u = 0; u < qn; ++u) {
            ui deg = 0; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) {
                ui u1 = nbrs[i];
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

        initBlackWhiteCandidateAdjIndex();

#if CDE_BLACK_WHITE_STATIC_COLOR
        initBWColor();
#endif

#if CDE_BLACK_WHITE_FIXED_ORDER
        initBlackWhiteFixedEdgePriorities();
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

        ui root = bw_static_root < qn ? bw_static_root : selectBWRoot();
        for (ui v0 : candidates[root]) {
            BlackWhiteState state;
            initBlackWhiteState(state);
            if (!commitRootBlack(state, root, v0)) {
                continue;
            }
            bwSearch(state, 0);
            if (outputLimitReached()) break;
        }

        stats.search_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.search_time;
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
        long long search_time = 0;
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
        long long white_bucket_rebuilds = 0;
        long long graph_has_edge_checks = 0;
        long long candidate_range_hits = 0;
        long long candidate_range_misses = 0;
        long long candidate_intersection_calls = 0;
        long long candidate_edge_check_calls = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        printf("\n--- CDE-Black-White Time Analysis ---\n");
#ifdef NDEBUG
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms\n", stats.filter_time / 1000.0);
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
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
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
#endif
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
        printf("White Rebuilds:      %lld\n", stats.white_bucket_rebuilds);
        printf("Graph hasEdge Calls: %lld\n", stats.graph_has_edge_checks);
        printf("Range Hits/Misses:   %lld / %lld\n",
            stats.candidate_range_hits, stats.candidate_range_misses);
        printf("Range Intersections: %lld\n", stats.candidate_intersection_calls);
        printf("Range Edge Checks:   %lld\n", stats.candidate_edge_check_calls);
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
    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;
    ui threshold;
    ui qn, gn;
    ui label_count;
    ui max_g_deg;
    vector<vector<ui>> q_neighbors;
    vector<ui> q_degree;
    vector<vector<char>> q_neighbor_is_bridge;

    vector<MyBitset> candidates;
    ui bw_static_root = 0;
    vector<VertexColor> bw_static_color;
    vector<vector<unsigned long long>> bw_static_edge_support;
#if CDE_BLACK_WHITE_FIXED_ORDER
    vector<vector<ui>> bw_static_edge_priority;
#endif

    vector<BlackWhiteUndo> black_white_undo;

    vector<unordered_map<unsigned long long, CandidateAdjRange>> bw_candidate_adj_index;
    vector<ui> bw_candidate_adj_pool;
    vector<CandidateAdjRange> bw_candidate_range_buffer;
    vector<ui> bw_candidate_source_buffer;
    vector<ui> bw_candidate_result_buffer;
    vector<ui> bw_candidate_intersection_buffer;
    vector<ui> bw_candidate_batch_mark;
    vector<ui> bw_candidate_batch_pos;
    vector<ui> bw_candidate_batch_present_hits;
    vector<ui> bw_candidate_batch_undecided_hits;
    vector<unsigned char> bw_candidate_batch_valid;
    ui bw_candidate_batch_token = 0;

    vector<vector<BlackWhiteActiveEdge>> bw_top_edges_buffer_by_depth;
    vector<vector<ui>> bw_white_neighbors_buffer_by_depth;

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
        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        q_neighbors.clear();
        q_degree.clear();
        q_neighbor_is_bridge.clear();

        bw_top_edges_buffer_by_depth.clear();
        bw_top_edges_buffer_by_depth.resize((size_t)qn + 1);
        bw_white_neighbors_buffer_by_depth.clear();
        bw_white_neighbors_buffer_by_depth.resize((size_t)qn + 1);
        black_white_undo.clear();
        bw_candidate_adj_index.clear();
        bw_candidate_adj_pool.clear();
        bw_candidate_range_buffer.clear();
        bw_candidate_source_buffer.clear();
        bw_candidate_result_buffer.clear();
        bw_candidate_intersection_buffer.clear();
        bw_candidate_batch_mark.assign(gn, 0);
        bw_candidate_batch_pos.assign(gn, 0);
        bw_candidate_batch_present_hits.clear();
        bw_candidate_batch_undecided_hits.clear();
        bw_candidate_batch_valid.clear();
        bw_candidate_batch_token = 0;
        stats = TimeStats();
        bw_static_root = 0;
        bw_static_color.clear();
        bw_static_edge_support.clear();
#if CDE_BLACK_WHITE_FIXED_ORDER
        bw_static_edge_priority.clear();
#endif
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

    unsigned long long blackWhiteCandidateAdjKey(ui data_vertex, ui query_neighbor) const
    {
        return ((unsigned long long)data_vertex << 32) |
            (unsigned long long)query_neighbor;
    }

    void initBlackWhiteCandidateAdjIndex()
    {
        bw_candidate_adj_index.clear();
        bw_candidate_adj_index.resize(qn);
        bw_candidate_adj_pool.clear();

        size_t total_candidate_count = 0;
        for (ui u = 0; u < qn; ++u) {
            total_candidate_count += (size_t)candidates[u].size();
            bw_candidate_adj_index[u].reserve(
                (size_t)candidates[u].size() * q_neighbors[u].size());
        }
        bw_candidate_source_buffer.reserve(total_candidate_count);
        bw_candidate_result_buffer.reserve(total_candidate_count);
        bw_candidate_intersection_buffer.reserve(total_candidate_count);
        bw_candidate_range_buffer.reserve(qn);

        for (ui u = 0; u < qn; ++u) {
            for (ui v : candidates[u]) {
                ui degree = 0;
                const ui *neighbors = data_graph->getVertexNeighbors(v, degree);

                for (ui query_neighbor : q_neighbors[u]) {
                    size_t begin = bw_candidate_adj_pool.size();
                    for (ui i = 0; i < degree; ++i) {
                        ui data_neighbor = neighbors[i];
                        if (candidates[query_neighbor].contains(data_neighbor)) {
                            bw_candidate_adj_pool.push_back(data_neighbor);
                        }
                    }

                    ui len = (ui)(bw_candidate_adj_pool.size() - begin);
                    if (len == 0) {
                        continue;
                    }

                    CandidateAdjRange range;
                    range.begin = begin;
                    range.len = len;
                    bw_candidate_adj_index[u].emplace(
                        blackWhiteCandidateAdjKey(v, query_neighbor), range);
                }
            }
        }
    }

    const CandidateAdjRange *findBlackWhiteCandidateAdjRange(ui from_query,
        ui from_data, ui to_query) const
    {
        if (from_query >= bw_candidate_adj_index.size()) {
            return nullptr;
        }

        const auto &index = bw_candidate_adj_index[from_query];
        auto it = index.find(blackWhiteCandidateAdjKey(from_data, to_query));
        if (it == index.end()) {
            return nullptr;
        }
        return &it->second;
    }

    const ui *blackWhiteRangeBegin(const CandidateAdjRange &range) const
    {
        return bw_candidate_adj_pool.data() + range.begin;
    }

    const ui *blackWhiteRangeEnd(const CandidateAdjRange &range) const
    {
        return blackWhiteRangeBegin(range) + range.len;
    }

    bool blackWhiteRangeContains(const CandidateAdjRange &range, ui value) const
    {
        return std::binary_search(blackWhiteRangeBegin(range),
            blackWhiteRangeEnd(range), value);
    }

    bool blackWhiteCandidateAdjacent(ui from_query, ui from_data,
        ui to_query, ui to_data)
    {
        stats.candidate_edge_check_calls++;
        const CandidateAdjRange *range =
            findBlackWhiteCandidateAdjRange(from_query, from_data, to_query);
        if (range == nullptr) {
            stats.candidate_range_misses++;
            return false;
        }

        stats.candidate_range_hits++;
        return blackWhiteRangeContains(*range, to_data);
    }

    bool blackWhiteGraphHasEdge(ui u, ui v)
    {
        stats.graph_has_edge_checks++;
        return data_graph->hasEdge(u, v);
    }

    bool blackWhiteCandidateAdjacentFromAnchor(ui anchor_query, ui anchor_data,
        ui target_query, ui target_data)
    {
        stats.candidate_edge_check_calls++;
        const CandidateAdjRange *range =
            findBlackWhiteCandidateAdjRange(anchor_query, anchor_data, target_query);
        if (range == nullptr) {
            stats.candidate_range_misses++;
            return false;
        }

        stats.candidate_range_hits++;
        return blackWhiteRangeContains(*range, target_data);
    }

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
        state.white_candidate_pool.clear();
        state.white_candidate_pool.reserve(stats.filter_candidate_count);
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

#if CDE_BLACK_WHITE_FIXED_ORDER
    struct BlackWhiteFixedEdgePriorityEntry {
        ui u = 0;
        ui anchor = 0;
        unsigned long long pair_support = 0;
        ui u_candidate_count = 0;
        ui anchor_candidate_count = 0;
    };

    void initBlackWhiteFixedEdgePriorities()
    {
        const ui invalid_priority = std::numeric_limits<ui>::max();
        bw_static_edge_priority.assign(qn, vector<ui>(qn, invalid_priority));

        vector<BlackWhiteFixedEdgePriorityEntry> entries;
        size_t directed_edge_count = 0;
        for (ui u = 0; u < qn; ++u) {
            directed_edge_count += q_neighbors[u].size();
        }
        entries.reserve(directed_edge_count);

        for (ui u = 0; u < qn; ++u) {
            for (ui anchor : q_neighbors[u]) {
                if (u >= anchor) {
                    continue;
                }

                unsigned long long pair_support =
                    bw_static_edge_support[u][anchor];
                ui u_candidate_count = (ui)candidates[u].size();
                ui anchor_candidate_count = (ui)candidates[anchor].size();

                BlackWhiteFixedEdgePriorityEntry forward;
                forward.u = u;
                forward.anchor = anchor;
                forward.pair_support = pair_support;
                forward.u_candidate_count = u_candidate_count;
                forward.anchor_candidate_count = anchor_candidate_count;
                entries.push_back(forward);

                BlackWhiteFixedEdgePriorityEntry reverse;
                reverse.u = anchor;
                reverse.anchor = u;
                reverse.pair_support = pair_support;
                reverse.u_candidate_count = anchor_candidate_count;
                reverse.anchor_candidate_count = u_candidate_count;
                entries.push_back(reverse);
            }
        }

        std::sort(entries.begin(), entries.end(),
            [&](const BlackWhiteFixedEdgePriorityEntry &lhs,
                const BlackWhiteFixedEdgePriorityEntry &rhs) {
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
            const BlackWhiteFixedEdgePriorityEntry &entry = entries[rank];
            bw_static_edge_priority[entry.u][entry.anchor] = (ui)rank;
        }
    }
#endif

    ui selectBWRoot()
    {
        ui root = 0;
        for (ui u = 1; u < qn; ++u) {
            size_t cand_u = candidates[u].size();
            size_t cand_root = candidates[root].size();

            ui deg_u = q_degree[u];
            ui deg_root = q_degree[root];

            if (cand_u * deg_root > cand_root * deg_u) root = u;
        }
        return root;
    }

    void initBWColor()
    {
        bw_static_root = selectBWRoot();
        bw_static_color.assign(qn, COLOR_WHITE);
        bw_static_color[bw_static_root] = COLOR_BLACK;
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
                state.white_candidate_pool.resize(undo.old_size);
                state.white[undo.u] = undo.old_white;
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
        const vector<ui> &candidates_to_store)
    {
        BlackWhiteUndo undo;
        undo.kind = BW_UNDO_WHITE_BUCKET;
        undo.u = u;
        undo.old_size = state.white_candidate_pool.size();
        undo.old_white = state.white[u];
        black_white_undo.push_back(std::move(undo));

        WhiteCandidateBuckets bucket;
        bucket.begin = state.white_candidate_pool.size();
        assert(candidates_to_store.size() <=
            (size_t)std::numeric_limits<ui>::max());
        bucket.count = (ui)candidates_to_store.size();
        bucket.feasible_count = bucket.count;
        state.white_candidate_pool.insert(state.white_candidate_pool.end(),
            candidates_to_store.begin(), candidates_to_store.end());
        state.white[u] = bucket;
    }

    bool computeBlackNeighborDelta(const BlackWhiteState &state, ui u, ui v,
        ui cost, ui &delta)
    {
        delta = 0;
        for (ui neighbor : q_neighbors[u]) {
            if (!isBlack(state, neighbor)) {
                continue;
            }

            ui mapped_neighbor = (ui)state.mapped_q[neighbor];
            bool adjacent = blackWhiteCandidateAdjacentFromAnchor(
                neighbor, mapped_neighbor, u, v);
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

    bool collectBlackWhitePositiveRanges(const BlackWhiteState &state, ui u,
        vector<CandidateAdjRange> &ranges)
    {
        ranges.clear();
        for (ui neighbor : q_neighbors[u]) {
            if (!isBlack(state, neighbor) ||
                getBlackWhiteEdgeState(state, u, neighbor) != EDGE_PRESENT) {
                continue;
            }

            const CandidateAdjRange *range = findBlackWhiteCandidateAdjRange(
                neighbor, (ui)state.mapped_q[neighbor], u);
            if (range == nullptr) {
                stats.candidate_range_misses++;
                return false;
            }
            stats.candidate_range_hits++;
            ranges.push_back(*range);
        }
        return true;
    }

    void buildBlackWhiteCandidateSourceFromRanges(
        vector<CandidateAdjRange> &ranges, vector<ui> &source)
    {
        source.clear();
        if (ranges.empty()) {
            return;
        }

        std::sort(ranges.begin(), ranges.end(),
            [](const CandidateAdjRange &lhs, const CandidateAdjRange &rhs) {
                return lhs.len < rhs.len;
            });

        if (ranges.size() == 1) {
            const CandidateAdjRange &range = ranges.front();
            source.assign(blackWhiteRangeBegin(range), blackWhiteRangeEnd(range));
            return;
        }

        size_t total_len = 0;
        for (const CandidateAdjRange &range : ranges) {
            total_len += range.len;
        }

        size_t edge_check_threshold = (size_t)ranges.front().len *
            ranges.size() * 8;
        if (total_len <= edge_check_threshold) {
            stats.candidate_intersection_calls++;
            const CandidateAdjRange &first = ranges.front();
            source.assign(blackWhiteRangeBegin(first), blackWhiteRangeEnd(first));

            for (size_t i = 1; i < ranges.size() && !source.empty(); ++i) {
                bw_candidate_intersection_buffer.clear();
                const CandidateAdjRange &range = ranges[i];
                std::set_intersection(source.begin(), source.end(),
                    blackWhiteRangeBegin(range), blackWhiteRangeEnd(range),
                    std::back_inserter(bw_candidate_intersection_buffer));
                source.swap(bw_candidate_intersection_buffer);
            }
            return;
        }

        const CandidateAdjRange &shortest = ranges.front();
        for (const ui *it = blackWhiteRangeBegin(shortest);
            it != blackWhiteRangeEnd(shortest); ++it) {
            ui candidate = *it;
            bool supported = true;
            for (size_t i = 1; i < ranges.size(); ++i) {
                stats.candidate_edge_check_calls++;
                if (!blackWhiteRangeContains(ranges[i], candidate)) {
                    supported = false;
                    break;
                }
            }
            if (supported) {
                source.push_back(candidate);
            }
        }
    }

    void appendFeasibleWhiteCandidate(const BlackWhiteState &state, ui u,
        ui candidate, ui cost, vector<ui> &result)
    {
        if (isDataVertexUsed(state, candidate)) {
            return;
        }

        ui delta = 0;
        if (computeBlackNeighborDelta(state, u, candidate, cost, delta)) {
            result.push_back(candidate);
        }
    }

    bool whiteBucketContains(const BlackWhiteState &state,
        const WhiteCandidateBuckets &bucket, ui candidate) const
    {
        assert(bucket.begin + bucket.count <= state.white_candidate_pool.size());
        const ui *begin = state.white_candidate_pool.data() + bucket.begin;
        const ui *end = begin + bucket.count;
        return std::binary_search(begin, end, candidate);
    }

    ui nextBlackWhiteCandidateBatchToken()
    {
        if (bw_candidate_batch_mark.size() < gn) {
            bw_candidate_batch_mark.assign(gn, 0);
            bw_candidate_batch_pos.assign(gn, 0);
            bw_candidate_batch_token = 0;
        }

        bw_candidate_batch_token++;
        if (bw_candidate_batch_token == 0) {
            std::fill(bw_candidate_batch_mark.begin(),
                bw_candidate_batch_mark.end(), 0);
            bw_candidate_batch_token = 1;
        }
        return bw_candidate_batch_token;
    }

    void addBlackWhiteBatchRangeHits(const CandidateAdjRange *range, ui token,
        vector<ui> &hits)
    {
        if (range == nullptr) {
            return;
        }

        for (const ui *it = blackWhiteRangeBegin(*range);
            it != blackWhiteRangeEnd(*range); ++it) {
            ui candidate = *it;
            if (candidate >= bw_candidate_batch_mark.size() ||
                bw_candidate_batch_mark[candidate] != token) {
                continue;
            }

            ui pos = bw_candidate_batch_pos[candidate];
            if (pos < hits.size()) {
                hits[pos]++;
            }
        }
    }

    void invalidateBlackWhiteBatchRange(const CandidateAdjRange *range, ui token)
    {
        if (range == nullptr) {
            return;
        }

        for (const ui *it = blackWhiteRangeBegin(*range);
            it != blackWhiteRangeEnd(*range); ++it) {
            ui candidate = *it;
            if (candidate >= bw_candidate_batch_mark.size() ||
                bw_candidate_batch_mark[candidate] != token) {
                continue;
            }

            ui pos = bw_candidate_batch_pos[candidate];
            if (pos < bw_candidate_batch_valid.size()) {
                bw_candidate_batch_valid[pos] = 0;
            }
        }
    }

    void appendFeasibleCandidatesBatch(const BlackWhiteState &state, ui u,
        ui cost, const vector<ui> &source, vector<ui> &result)
    {
        if (source.empty()) {
            return;
        }

        ui token = nextBlackWhiteCandidateBatchToken();
        size_t source_count = source.size();
        bw_candidate_batch_present_hits.assign(source_count, 0);
        bw_candidate_batch_undecided_hits.assign(source_count, 0);
        bw_candidate_batch_valid.assign(source_count, 0);

        for (size_t i = 0; i < source_count; ++i) {
            ui candidate = source[i];
            if (candidate >= gn) {
                continue;
            }
            bw_candidate_batch_mark[candidate] = token;
            bw_candidate_batch_pos[candidate] = (ui)i;
            bw_candidate_batch_valid[i] =
                isDataVertexUsed(state, candidate) ? 0 : 1;
        }

        ui present_count = 0;
        ui undecided_count = 0;
        for (ui neighbor : q_neighbors[u]) {
            if (!isBlack(state, neighbor)) {
                continue;
            }

            ui mapped_neighbor = (ui)state.mapped_q[neighbor];
            EdgeState state_uv = getBlackWhiteEdgeState(state, u, neighbor);
            const CandidateAdjRange *range =
                findBlackWhiteCandidateAdjRange(neighbor, mapped_neighbor, u);

            stats.candidate_edge_check_calls += (long long)source_count;
            if (range == nullptr) {
                stats.candidate_range_misses++;
            }
            else {
                stats.candidate_range_hits++;
            }

            if (state_uv == EDGE_PRESENT) {
                present_count++;
                addBlackWhiteBatchRangeHits(range, token,
                    bw_candidate_batch_present_hits);
            }
            else if (state_uv == EDGE_MISSING) {
                invalidateBlackWhiteBatchRange(range, token);
            }
            else {
                undecided_count++;
                addBlackWhiteBatchRangeHits(range, token,
                    bw_candidate_batch_undecided_hits);
            }
        }

        for (size_t i = 0; i < source_count; ++i) {
            if (!bw_candidate_batch_valid[i]) {
                continue;
            }
            if (bw_candidate_batch_present_hits[i] != present_count) {
                continue;
            }

            ui adjacent_undecided = bw_candidate_batch_undecided_hits[i];
            ui delta = undecided_count > adjacent_undecided
                ? undecided_count - adjacent_undecided : 0;
            if (cost + delta <= threshold) {
                result.push_back(source[i]);
            }
        }
    }

    void copyWhiteBucketCandidates(const BlackWhiteState &state,
        const WhiteCandidateBuckets &bucket, vector<ui> &target) const
    {
        assert(bucket.begin + bucket.count <= state.white_candidate_pool.size());
        target.assign(state.white_candidate_pool.begin() + bucket.begin,
            state.white_candidate_pool.begin() + bucket.begin + bucket.count);
    }

    void filterSourceByWhiteBucket(const BlackWhiteState &state,
        const WhiteCandidateBuckets &bucket, const vector<ui> &source,
        vector<ui> &target) const
    {
        target.clear();
        for (ui candidate : source) {
            if (whiteBucketContains(state, bucket, candidate)) {
                target.push_back(candidate);
            }
        }
    }

    void collectAllBlackWhiteCandidates(ui u, vector<ui> &target)
    {
        target.clear();
        for (int candidate : candidates[u]) {
            target.push_back((ui)candidate);
        }
    }

    void appendFeasibleCandidatesFromWhiteBucket(const BlackWhiteState &state,
        ui u, ui cost, const WhiteCandidateBuckets &bucket, vector<ui> &result)
    {
        assert(bucket.begin + bucket.count <= state.white_candidate_pool.size());
        for (ui i = 0; i < bucket.count; ++i) {
            ui candidate = state.white_candidate_pool[bucket.begin + i];
            appendFeasibleWhiteCandidate(state, u, candidate, cost, result);
        }
    }

    bool buildWhiteCandidateBuffer(BlackWhiteState &state, ui u, ui cost,
        const WhiteCandidateBuckets *existing_bucket)
    {
        bw_candidate_result_buffer.clear();
        if (cost > threshold) {
            return false;
        }
        stats.white_bucket_rebuilds++;

#if CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SOURCE
        if (!collectBlackWhitePositiveRanges(state, u,
            bw_candidate_range_buffer)) {
            return false;
        }

        if (!bw_candidate_range_buffer.empty()) {
            buildBlackWhiteCandidateSourceFromRanges(bw_candidate_range_buffer,
                bw_candidate_source_buffer);

            if (existing_bucket != nullptr &&
                existing_bucket->count <= bw_candidate_source_buffer.size()) {
                copyWhiteBucketCandidates(state, *existing_bucket,
                    bw_candidate_intersection_buffer);
                appendFeasibleCandidatesBatch(state, u, cost,
                    bw_candidate_intersection_buffer, bw_candidate_result_buffer);
            }
            else {
                const vector<ui> *source = &bw_candidate_source_buffer;
                if (existing_bucket != nullptr) {
                    filterSourceByWhiteBucket(state, *existing_bucket,
                        bw_candidate_source_buffer,
                        bw_candidate_intersection_buffer);
                    source = &bw_candidate_intersection_buffer;
                }
                appendFeasibleCandidatesBatch(state, u, cost, *source,
                    bw_candidate_result_buffer);
            }
        }
        else if (existing_bucket != nullptr) {
            copyWhiteBucketCandidates(state, *existing_bucket,
                bw_candidate_intersection_buffer);
            appendFeasibleCandidatesBatch(state, u, cost,
                bw_candidate_intersection_buffer, bw_candidate_result_buffer);
        }
        else {
#endif
        if (existing_bucket != nullptr) {
            copyWhiteBucketCandidates(state, *existing_bucket,
                bw_candidate_intersection_buffer);
            appendFeasibleCandidatesBatch(state, u, cost,
                bw_candidate_intersection_buffer, bw_candidate_result_buffer);
        }
        else {
            collectAllBlackWhiteCandidates(u, bw_candidate_source_buffer);
            appendFeasibleCandidatesBatch(state, u, cost,
                bw_candidate_source_buffer, bw_candidate_result_buffer);
        }
#if CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SOURCE
        }
#endif

        return !bw_candidate_result_buffer.empty();
    }

    bool rebuildWhiteCandidatesForCurrentState(BlackWhiteState &state,
        ui white_u, ui cost)
    {
        WhiteCandidateBuckets old_bucket = state.white[white_u];
        if (!buildWhiteCandidateBuffer(state, white_u, cost, &old_bucket)) {
            return false;
        }

        replaceBlackWhiteBucket(state, white_u, bw_candidate_result_buffer);
        return true;
    }

    bool buildWhiteCandidateBuckets(BlackWhiteState &state, ui u, ui cost)
    {
        if (cost > threshold || !isSelectedByBlackNeighbor(state, u)) {
            return false;
        }
        return buildWhiteCandidateBuffer(state, u, cost, nullptr);
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
        ui white_candidate_count) const
    {
        if (!shouldPreferStaticWhite(u) || white_candidate_count == 0) {
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
            bool adjacent = blackWhiteCandidateAdjacentFromAnchor(
                neighbor, (ui)state.mapped_q[neighbor], u, v);
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
            bool adjacent = blackWhiteCandidateAdjacentFromAnchor(
                neighbor, mapped_neighbor, white_u, candidate);
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

        if (!buildWhiteCandidateBuckets(state, u, cost)) {
            return false;
        }
        if (!chooseBlackWhite(state, u,
            (ui)bw_candidate_result_buffer.size())) {
            return false;
        }

        size_t mark = markBlackWhiteState();
        setBlackWhiteColor(state, u, COLOR_WHITE);
        replaceBlackWhiteBucket(state, u, bw_candidate_result_buffer);
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

#if CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_ANCHOR_BRANCH
        if (required_anchor < qn && isBlack(state, required_anchor) &&
            getBlackWhiteEdgeState(state, u, required_anchor) == EDGE_PRESENT) {
            ui mapped_anchor = (ui)state.mapped_q[required_anchor];
            const CandidateAdjRange *range = findBlackWhiteCandidateAdjRange(
                required_anchor, mapped_anchor, u);
            if (range == nullptr) {
                return false;
            }

            for (const ui *it = blackWhiteRangeBegin(*range);
                it != blackWhiteRangeEnd(*range); ++it) {
                ui candidate = *it;
                if (try_candidate(candidate)) {
                    return true;
                }
            }
            return emitted_branch;
        }
#endif

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
        WhiteCandidateBuckets white_bucket = state.white[white_u];
        assert(white_bucket.begin + white_bucket.count <=
            state.white_candidate_pool.size());
        for (ui candidate_idx = 0; candidate_idx < white_bucket.count;
            ++candidate_idx) {
            ui candidate =
                state.white_candidate_pool[white_bucket.begin + candidate_idx];
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
#if CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SUPPORT
        const CandidateAdjRange *range =
            findBlackWhiteCandidateAdjRange(anchor, mapped_anchor, u);
        if (range == nullptr) {
            return 0.0;
        }

        for (const ui *it = blackWhiteRangeBegin(*range);
            it != blackWhiteRangeEnd(*range); ++it) {
            ui v = *it;
            if (isDataVertexUsed(state, v)) {
                continue;
            }
            support_count += 1.0;
        }
#else
        for (int candidate : candidates[u]) {
            ui v = (ui)candidate;
            if (isDataVertexUsed(state, v)) {
                continue;
            }
            if (data_graph->hasEdge(mapped_anchor, v)) {
                support_count += 1.0;
            }
        }
#endif
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
#if CDE_BLACK_WHITE_FIXED_ORDER
        ui lhs_priority = bw_static_edge_priority[lhs.u][lhs.anchor];
        ui rhs_priority = bw_static_edge_priority[rhs.u][rhs.anchor];
        if (lhs_priority != rhs_priority) {
            return lhs_priority < rhs_priority;
        }
        if (lhs.u != rhs.u) {
            return lhs.u < rhs.u;
        }
        return lhs.anchor < rhs.anchor;
#else
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
#endif
    }

    void selectTopBlackWhiteActiveEdges(ui max_count,
        vector<BlackWhiteActiveEdge> &top_edges)
    {
        size_t selected_limit = top_edges.size();
        if ((size_t)max_count < selected_limit) {
            selected_limit = (size_t)max_count;
        }

#if CDE_BLACK_WHITE_FIXED_ORDER
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
#elif CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
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
#if !CDE_BLACK_WHITE_FIXED_ORDER
                if (isBlack(state, anchor)) {
                    edge.rank_support = estimateBlackAnchorSupport(state, u, anchor);
                }
                else {
                    edge.rank_support = estimateWhiteAnchorSupport(state, anchor);
                }
#endif
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
        vector<vector<ui>> &buckets, ui &feasible_count, ui &min_delta)
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

        WhiteCandidateBuckets white = state.white[white_u];
        assert(white.begin + white.count <= state.white_candidate_pool.size());
        for (ui candidate_idx = 0; candidate_idx < white.count;
            ++candidate_idx) {
            ui candidate =
                state.white_candidate_pool[white.begin + candidate_idx];
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
        ui cost, vector<TerminalTailVertex> &tail_vertices)
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
