#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_H_

#include "configuration/config.h"
#include "configuration/types.h"
#include "graph/graph.h"
#include "matching/algorithms/cde_black_white/types.h"
#include "matching/run_matching.h"
#include "utility/mybitset.h"

#ifndef CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SOURCE
#define CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_SOURCE 1
#endif
#ifndef CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_ANCHOR_BRANCH
#define CDE_BLACK_WHITE_USE_CANDIDATE_RANGE_ANCHOR_BRANCH 1
#endif

#define CDE_BLACK_WHITE_USE_FLAT_HASH_MAP 1

#if CDE_BLACK_WHITE_FIXED_ORDER && CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY
#error "CDE_BLACK_WHITE_FIXED_ORDER and CDE_BLACK_WHITE_TOPK_SUPPORT_DECAY are mutually exclusive."
#endif

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
    using ContinueBranch = std::function<void(SearchState &, ui)>;

    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;
    ui threshold;
    ui qn, gn;
    ui label_count;
    ui max_g_deg;
    vector<vector<ui>> q_neighbors;
    vector<ui> q_degree;
    vector<vector<char>> q_bridge_matrix;

    vector<MyBitset> candidates;
    ui static_root = 0;
    vector<VertexColor> static_color;
    vector<vector<unsigned long long>> static_edge_support;
#if CDE_BLACK_WHITE_FIXED_ORDER
    vector<vector<ui>> static_edge_priority;
#endif

    vector<UndoRecord> undo_stack;

#if CDE_BLACK_WHITE_USE_FLAT_HASH_MAP
    vector<absl::flat_hash_map<unsigned long long, pair<size_t, ui>>> candidate_adj_index;
#else
    vector<unordered_map<unsigned long long, pair<size_t, ui>>> candidate_adj_index;
#endif
    vector<ui> candidate_adj_pool;
    vector<pair<size_t, ui>> candidate_range_buffer;
    vector<ui> candidate_source_buffer;
    vector<ui> candidate_result_buffer;
    vector<ui> candidate_result_delta_buffer;
    vector<ui> candidate_intersection_buffer;

    vector<ui> candidate_batch_mark;
    vector<ui> candidate_batch_pos;
    vector<ui> candidate_batch_present_hits;
    vector<ui> candidate_batch_undecided_hits;
    vector<unsigned char> candidate_batch_valid;
    ui candidate_batch_token = 0;

    vector<vector<AnchorEdge>> top_edges_buffer;
    vector<vector<pair<ui, ui>>> branch_cands_buffer;
    vector<vector<ui>> white_neighbors_buffer;
    vector<int> component_id_buffer;
    vector<ui> component_queue_buffer;
    vector<ui> component_edge_counts_buffer;

    bool outputLimitReached() const;
    void resetState();
    void resetBuffers();
    ui chooseRoot();
    void initColors();
    void initQueryBridge();
    void tarjan(ui u, ui parent, vector<int> &dfn, vector<int> &low, int &tim);
    bool isQueryBridgeEdge(ui u, ui v) const;

    bool runCandidateFiltering();

    unsigned long long adjKey(ui data_vertex, ui query_neighbor) const;
    void buildAdjIndex();
    const pair<size_t, ui> *findAdjRange(ui from_query, ui from_data,
        ui to_query) const;
    const ui *rangeBegin(const pair<size_t, ui> &range) const;
    const ui *rangeEnd(const pair<size_t, ui> &range) const;
    bool rangeHas(const pair<size_t, ui> &range, ui value) const;
    bool anchorAdjacent(ui anchor_query, ui anchor_data, ui target_query,
        ui target_data);

    void initState(SearchState &state) const;
    size_t edgeIdx(ui u, ui v) const;
    EdgeState getEdge(const SearchState &state, ui u, ui v) const;
    void setEdgeRaw(SearchState &state, ui u, ui v,
        EdgeState edge_state_value) const;
    void addAnchorEdgeRaw(SearchState &state, ui u, ui anchor) const;
    void removeAnchorEdgeRaw(SearchState &state, ui u, ui anchor) const;
    void refreshAnchorEdge(SearchState &state, ui u, ui v) const;
    vector<AnchorEdge> &topEdgesBuffer(ui depth);
    vector<pair<ui, ui>> &branchCandsBuffer(ui depth);
    vector<ui> &whiteNbrsBuffer(ui depth);
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
        const vector<ui> &candidates_to_store,
        const vector<ui> &candidate_deltas_to_store);

    bool calcBlackDelta(const SearchState &state, ui u, ui v, ui cost,
        ui &delta);
    bool collectRequiredRanges(const SearchState &state, ui u,
        vector<pair<size_t, ui>> &ranges);
    void intersectRequiredRanges(vector<pair<size_t, ui>> &ranges,
        vector<ui> &source);
    bool bucketHas(const SearchState &state,
        const WhiteCands &bucket, ui candidate) const;
    ui nextBatchToken();
    void addRangeHits(const pair<size_t, ui> *range, ui token,
        vector<ui> &hits);
    void invalidateRange(const pair<size_t, ui> *range, ui token);
    void addFeasibleBatch(const SearchState &state, ui u, ui cost,
        const vector<ui> &source, vector<ui> &result,
        vector<ui> &result_deltas);
    void copyBucketCands(const SearchState &state,
        const WhiteCands &bucket, vector<ui> &target) const;
    void filterByBucket(const SearchState &state,
        const WhiteCands &bucket, const vector<ui> &source,
        vector<ui> &target) const;
    void collectAllCands(ui u, vector<ui> &target);
    bool buildWhiteCands(SearchState &state, ui u, ui cost);
    bool refreshWhiteByBlack(SearchState &state, ui white_u, ui black_u,
        ui black_v, ui cost);

#if CDE_BLACK_WHITE_FIXED_ORDER
    void initFixedEdgePriorities();
#endif
    double blackSupport(const SearchState &state, ui u, ui anchor) const;
    double whiteSupport(const SearchState &state, ui anchor) const;
    bool betterEdge(const AnchorEdge &lhs, const AnchorEdge &rhs) const;
    void selectTopEdges(ui max_count, vector<AnchorEdge> &top_edges);
    bool isActiveAnchorEdge(const SearchState &state,
        const AnchorEdge &edge) const;
    size_t labelFrontierComponents(const SearchState &state, vector<int> &component_id, vector<ui> &queue);
    void trimToCompleteComponent(const SearchState &state,
        vector<AnchorEdge> &top_edges);
    bool collectActiveEdges(const SearchState &state, ui max_count,
        vector<AnchorEdge> &top_edges);

    bool buildTailWhite(const SearchState &state, ui cost, TailWhite &tail_white);
    bool buildTailWhites(const SearchState &state, ui cost,
        vector<TailWhite> &tail_vertices);
    void enumTailWhites(SearchState &state, size_t pos, ui cost,
        vector<TailWhite> &tail_vertices);
    void emitResult(const SearchState &state);

    bool tryMapBlackWithDelta(SearchState &state, ui cost, ui u, ui v,
        ui delta, ui &next_cost);
    bool tryMapWhite(SearchState &state, ui cost, ui white_u,
        ui candidate, ui &next_cost);
    bool shouldExpandAsWhite(const SearchState &state, ui u,
        const vector<pair<ui, ui>> &anchor_candidates) const;
    void branchWhite(SearchState &state, ui cost, ui u);
    void branchBlack(SearchState &state, ui cost, ui u,
        const vector<pair<ui, ui>> &anchor_candidates);
    void branchMatWhite(SearchState &state, ui cost, ui white_u,
        const ContinueBranch &continue_branch);
    void branchMatWhites(SearchState &state, ui cost,
        const vector<ui> &white_vertices, size_t pos,
        const ContinueBranch &continue_branch);
    void branchBlackAnchor(SearchState &state, ui cost, ui u, ui anchor);
    void branchEdges(SearchState &state, ui cost,
        const vector<AnchorEdge> &top_edges, size_t edge_idx);
    void search(SearchState &state, ui cost);
};

} // namespace cde_black_white

void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph,
    std::vector<std::vector<std::pair<ui, ui>>> &M_ANS, ui threshold);

#endif
