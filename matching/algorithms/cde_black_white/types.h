#ifndef MATCHING_ALGORITHMS_CDE_BLACK_WHITE_TYPES_H_
#define MATCHING_ALGORITHMS_CDE_BLACK_WHITE_TYPES_H_

#include "utility/utility.h"
#include <limits>
#include <utility>
#include <vector>

namespace cde_black_white {

using namespace std;

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

struct WhiteCands {
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

struct AnchorEdge {
    ui u = 0;
    ui anchor = 0;
    VertexColor anchor_color = COLOR_UNSELECTED;
    double support = std::numeric_limits<double>::max();
    ui live_anchor_count = 0;
    ui query_degree = 0;
};

struct SearchState {
    vector<int> mapped_q;
    vector<ui> used_data_vertices;
    vector<unsigned char> used_data_flag;
    vector<VertexColor> color;
    vector<EdgeState> edge_state;
    vector<AnchorEdge> anchor_edges;
    vector<int> anchor_edge_pos;
    vector<ui> anchor_count;
    vector<WhiteCands> white;
    vector<ui> white_candidate_pool;
    vector<pair<ui, ui>> part_M;
    ui selected_count = 0;
    ui white_count = 0;
};

enum UndoKind : unsigned char {
    UNDO_MAPPED_Q = 0,
    UNDO_COLOR = 1,
    UNDO_EDGE_STATE = 2,
    UNDO_USED_DATA_SIZE = 3,
    UNDO_PART_M_SIZE = 4,
    UNDO_SELECTED_COUNT = 5,
    UNDO_WHITE_COUNT = 6,
    UNDO_WHITE_BUCKET = 7
};

struct UndoRecord {
    UndoKind kind = UNDO_MAPPED_Q;
    ui u = 0;
    ui v = 0;
    int old_mapped_q = -1;
    VertexColor old_color = COLOR_UNSELECTED;
    EdgeState old_edge_uv = EDGE_UNDECIDED;
    EdgeState old_edge_vu = EDGE_UNDECIDED;
    size_t old_size = 0;
    ui old_count = 0;
    WhiteCands old_white;
};

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
    long long terminal_tail_calls = 0;
    long long white_bucket_rebuilds = 0;
    long long graph_has_edge_checks = 0;
    long long candidate_range_hits = 0;
    long long candidate_range_misses = 0;
    long long candidate_intersection_calls = 0;
    long long candidate_edge_check_calls = 0;
    size_t result_count = 0;
    bool output_limit_reached = false;
};

struct TailWhite {
    ui u = 0;
    ui feasible_count = 0;
    ui min_delta = std::numeric_limits<ui>::max();
    vector<vector<ui>> buckets;
};


} // namespace cde_black_white

#endif
