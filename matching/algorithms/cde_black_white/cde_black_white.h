#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#include "configuration/types.h"
#include "graph/graph.h"
#include "matching/algorithms/cde_black_white/config.h"
#include "matching/algorithms/cde_black_white/types.h"
#include "matching/run_matching.h"
#include "utility/mybitset.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <queue>
#include <unordered_map>

namespace cde_black_white {

class MatchingSolver {
public:
    MatchingSolver();

    bool init(const Graph *q, const Graph *g, ui match_threshold);
    void match(vector<vector<pair<ui, ui>>> &results);
    void printStats() const;

    TimeStats stats;

private:
    struct CandidateFilter;
#if CDE_BLACK_WHITE_FIXED_ORDER
    struct FixedEdgePriorityEntry;
#endif
    using ContinueBranch = std::function<bool(SearchState &, ui)>;

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
    vector<int> component_id_buffer;
    vector<vector<ui>> component_frontiers_buffer;
    vector<ui> component_queue_buffer;
    vector<ui> component_edge_counts_buffer;
    vector<ui> component_seen_counts_buffer;

    bool outputLimitReached() const;
    void noteOutputLimitIfReached();
    void resetState();
    void resetBuffers();

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
    void addFrontierEdgeRaw(SearchState &state, ui u, ui anchor) const;
    void removeFrontierEdgeRaw(SearchState &state, ui u, ui anchor) const;
    void refreshFrontierEdge(SearchState &state, ui u, ui v) const;
    void refreshFrontierEdgesIncidentTo(SearchState &state, ui u) const;
    double blackSupport(const SearchState &state, ui u, ui anchor) const;
    double whiteSupport(const SearchState &state, ui anchor) const;
    bool betterEdge(const ActiveEdge &lhs, const ActiveEdge &rhs) const;
    void selectTopEdges(ui max_count, vector<ActiveEdge> &top_edges);
    bool isActiveFrontierEdge(const SearchState &state,
        const ActiveEdge &edge) const;
    size_t labelFrontierComponents(const SearchState &state,
        vector<int> &component_id, vector<vector<ui>> &component_frontiers,
        vector<ui> &queue);
    void restrictTopEdgesToCoveredComponent(const SearchState &state,
        vector<ActiveEdge> &top_edges);
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
    bool branchMatWhite(SearchState &state, ui cost, ui white_u,
        const ContinueBranch &continue_branch);
    bool branchMatWhites(SearchState &state, ui cost,
        const vector<ui> &white_vertices, size_t pos,
        const ContinueBranch &continue_branch);
    bool branchBlackAnchor(SearchState &state, ui cost, ui u, ui anchor);
    bool branchPresentEdge(SearchState &state, ui cost,
        const ActiveEdge &edge);
    void branchEdges(SearchState &state, ui cost,
        const vector<ActiveEdge> &top_edges, size_t edge_idx);
    void search(SearchState &state, ui cost);
};

} // namespace cde_black_white

void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph,
    std::vector<std::vector<std::pair<ui, ui>>> &M_ANS, ui threshold);

#endif
