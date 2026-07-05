#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#include "matching/algorithms/cde_black_white/config.h"
#include "matching/algorithms/cde_black_white/types.h"
#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/mybitset.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <queue>
#include <unordered_map>

namespace cde_black_white {

class MatchingSolver {
public:
    MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

    bool init(const Graph *q, const Graph *g, ui match_threshold)
    {
        // init info, filter candidates, build candidate adjacency index, and initialize static colors/edge priorities
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

        // init q_matrix
        q_neighbors.assign(qn, vector<ui>());
        q_degree.assign(qn, 0);
        for (ui u = 0; u < qn; ++u) {
            ui deg = 0; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) q_neighbors[u].push_back(nbrs[i]);
            q_degree[u] = (ui)q_neighbors[u].size();
        }

        Timer t_filter;
        bool res = runCandidateFiltering();
        stats.filter_time = t_filter.elapsed();
        if (!res) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        buildAdjIndex();

#if CDE_BLACK_WHITE_STATIC_COLOR
        initColors();
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

        ui root = static_root < qn ? static_root : chooseRoot();
        for (ui v0 : candidates[root]) {
            SearchState state;
            initState(state);
            if (!tryBindRoot(state, root, v0)) continue;
            search(state, 0);
            if (outputLimitReached()) break;
        }

        stats.search_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.search_time;
    }


    TimeStats stats;

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
        printf("Terminal Tail Calls: %lld\n", stats.terminal_tail_calls);
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
    struct CandidateFilter;
#if CDE_BLACK_WHITE_FIXED_ORDER
    struct FixedEdgePriorityEntry;
#endif

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
    ui static_root = 0;
    vector<VertexColor> static_color;
    vector<vector<unsigned long long>> static_edge_support;
#if CDE_BLACK_WHITE_FIXED_ORDER
    vector<vector<ui>> static_edge_priority;
#endif

    vector<UndoRecord> undo_stack;

    vector<unordered_map<unsigned long long, pair<size_t, ui>>> candidate_adj_index;
    vector<ui> candidate_adj_pool;
    vector<pair<size_t, ui>> candidate_range_buffer;
    vector<ui> candidate_source_buffer;
    vector<ui> candidate_result_buffer;
    vector<ui> candidate_intersection_buffer;
    vector<ui> candidate_batch_mark;
    vector<ui> candidate_batch_pos;
    vector<ui> candidate_batch_present_hits;
    vector<ui> candidate_batch_undecided_hits;
    vector<unsigned char> candidate_batch_valid;
    ui candidate_batch_token = 0;

    vector<vector<ActiveEdge>> top_edges_buffer_by_depth;
    vector<vector<ui>> white_neighbors_buffer_by_depth;


    bool outputLimitReached() const
    {
        return (size_t)MATCH_OUTPUT_LIMIT > 0 && stats.result_count >= (size_t)MATCH_OUTPUT_LIMIT;
    }

    void noteOutputLimitIfReached()
    {
        if (outputLimitReached()) stats.output_limit_reached = true;
    }

    void resetState()
    {
        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        q_neighbors.clear();
        q_degree.clear();
        q_neighbor_is_bridge.clear();

        top_edges_buffer_by_depth.clear();
        top_edges_buffer_by_depth.resize((size_t)qn + 1);
        white_neighbors_buffer_by_depth.clear();
        white_neighbors_buffer_by_depth.resize((size_t)qn + 1);
        undo_stack.clear();
        candidate_adj_index.clear();
        candidate_adj_pool.clear();
        candidate_range_buffer.clear();
        candidate_source_buffer.clear();
        candidate_result_buffer.clear();
        candidate_intersection_buffer.clear();
        candidate_batch_mark.assign(gn, 0);
        candidate_batch_pos.assign(gn, 0);
        candidate_batch_present_hits.clear();
        candidate_batch_undecided_hits.clear();
        candidate_batch_valid.clear();
        candidate_batch_token = 0;
        stats = TimeStats();
        static_root = 0;
        static_color.clear();
        static_edge_support.clear();
#if CDE_BLACK_WHITE_FIXED_ORDER
        static_edge_priority.clear();
#endif
    }

    bool runCandidateFiltering();

    unsigned long long adjKey(ui data_vertex, ui query_neighbor) const;
    void buildAdjIndex();
    const pair<size_t, ui> *findAdjRange(ui from_query, ui from_data,
        ui to_query) const;
    const ui *rangeBegin(const pair<size_t, ui> &range) const;
    const ui *rangeEnd(const pair<size_t, ui> &range) const;
    bool rangeHas(const pair<size_t, ui> &range, ui value) const;
    bool candAdjacent(ui from_query, ui from_data, ui to_query, ui to_data);
    bool hasDataEdge(ui u, ui v);
    bool anchorAdjacent(ui anchor_query, ui anchor_data, ui target_query,
        ui target_data);

    void initState(SearchState &state) const;
    size_t edgeIdx(ui u, ui v) const;
    EdgeState getEdge(const SearchState &state, ui u, ui v) const;
    void setEdgeRaw(SearchState &state, ui u, ui v,
        EdgeState edge_state_value) const;
    vector<ActiveEdge> &topEdgesBuffer(ui depth);
    vector<ui> &whiteNbrsBuffer(ui depth);
    bool tryBindRoot(SearchState &state, ui root, ui v) const;
    bool isDataVertexUsed(const SearchState &state, ui v) const;
    bool isSelected(const SearchState &state, ui u) const;
    bool isBlack(const SearchState &state, ui u) const;
    bool isWhite(const SearchState &state, ui u) const;
    size_t mark() const;
    void rollback(SearchState &state, size_t mark);
    void setMap(SearchState &state, ui u, int value);
    void setColor(SearchState &state, ui u, VertexColor value);
    void setEdge(SearchState &state, ui u, ui v, EdgeState edge_state_value);
    void pushUsed(SearchState &state, ui v);
    void pushMatch(SearchState &state, ui u, ui v);
    void setSelectedCnt(SearchState &state, ui value);
    void setWhiteCnt(SearchState &state, ui value);
    void replaceBucket(SearchState &state, ui u,
        const vector<ui> &candidates_to_store);

    bool calcBlackDelta(const SearchState &state, ui u, ui v, ui cost,
        ui &delta);
    bool collectPosRanges(const SearchState &state, ui u,
        vector<pair<size_t, ui>> &ranges);
    void buildRangeSource(vector<pair<size_t, ui>> &ranges,
        vector<ui> &source);
    void addFeasibleCand(const SearchState &state, ui u, ui candidate,
        ui cost, vector<ui> &result);
    bool bucketHas(const SearchState &state,
        const WhiteCandidateBuckets &bucket, ui candidate) const;
    ui nextBatchToken();
    void addRangeHits(const pair<size_t, ui> *range, ui token,
        vector<ui> &hits);
    void invalidateRange(const pair<size_t, ui> *range, ui token);
    void addFeasibleBatch(const SearchState &state, ui u, ui cost,
        const vector<ui> &source, vector<ui> &result);
    void copyBucketCands(const SearchState &state,
        const WhiteCandidateBuckets &bucket, vector<ui> &target) const;
    void filterByBucket(const SearchState &state,
        const WhiteCandidateBuckets &bucket, const vector<ui> &source,
        vector<ui> &target) const;
    void collectAllCands(ui u, vector<ui> &target);
    void addBucketCands(const SearchState &state, ui u, ui cost,
        const WhiteCandidateBuckets &bucket, vector<ui> &result);
    bool buildWhiteCands(SearchState &state, ui u, ui cost,
        const WhiteCandidateBuckets *existing_bucket);
    bool refreshWhiteCands(SearchState &state, ui white_u, ui cost);
    bool initWhiteCands(SearchState &state, ui u, ui cost);
    bool isSelectedByBlackNeighbor(const SearchState &state, ui u) const;
    void collectWhiteNbrs(const SearchState &state, ui u,
        vector<ui> &white_neighbors) const;
    bool refreshWhiteByBlack(SearchState &state, ui white_u, ui black_u,
        ui black_v, ui cost);

#if CDE_BLACK_WHITE_FIXED_ORDER
    void initFixedEdgePriorities();
#endif
    ui chooseRoot();
    void initColors();
    double blackSupport(const SearchState &state, ui u, ui anchor) const;
    double whiteSupport(const SearchState &state, ui anchor) const;
    bool betterEdge(const ActiveEdge &lhs, const ActiveEdge &rhs) const;
    void selectTopEdges(ui max_count, vector<ActiveEdge> &top_edges);
    bool collectActiveEdges(const SearchState &state, ui max_count,
        vector<ActiveEdge> &top_edges);
    ui chooseMatWhite(const SearchState &state) const;

    bool buildTailBuckets(const SearchState &state, ui white_u, ui cost,
        vector<vector<ui>> &buckets, ui &feasible_count, ui &min_delta);
    bool buildTailWhites(const SearchState &state, ui cost,
        vector<TerminalTailVertex> &tail_vertices);
    void enumTailWhites(SearchState &state, size_t pos, ui cost,
        vector<TerminalTailVertex> &tail_vertices);
    void emitResult(const SearchState &state);

    bool tryBindBlack(SearchState &state, ui cost, ui u, ui v,
        ui &next_cost);
    bool tryMaterializeWhite(SearchState &state, ui cost, ui white_u,
        ui candidate, ui bucket_delta, ui &next_cost);
    bool branchWhite(SearchState &state, ui cost, ui u);
    bool branchBlack(SearchState &state, ui cost, ui u,
        ui required_anchor = std::numeric_limits<ui>::max());
    template <typename Continue>
    bool branchMatWhite(SearchState &state, ui cost, ui white_u,
        Continue continue_branch);
    template <typename Continue>
    bool branchMatWhites(SearchState &state, ui cost,
        const vector<ui> &white_vertices, size_t pos,
        Continue continue_branch);
    bool branchBlackAnchor(SearchState &state, ui cost, ui u, ui anchor);
    bool branchPresentEdge(SearchState &state, ui cost,
        const ActiveEdge &edge);
    void branchEdges(SearchState &state, ui cost,
        const vector<ActiveEdge> &top_edges, size_t edge_idx);
    void search(SearchState &state, ui cost);
};

} // namespace cde_black_white

#include "matching/algorithms/cde_black_white/filtering.h"
#include "matching/algorithms/cde_black_white/candidate_index.h"
#include "matching/algorithms/cde_black_white/state_undo.h"
#include "matching/algorithms/cde_black_white/white_candidates.h"
#include "matching/algorithms/cde_black_white/ordering.h"
#include "matching/algorithms/cde_black_white/terminal_tail.h"
#include "matching/algorithms/cde_black_white/branching.h"

inline void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph, std::vector<std::vector<std::pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    cde_black_white::MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
