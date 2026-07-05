#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

void MatchingSolver::resetBuffers()
{
    top_edges_buffer_by_depth.clear();
    top_edges_buffer_by_depth.resize((size_t)qn + 1);
    white_neighbors_buffer_by_depth.clear();
    white_neighbors_buffer_by_depth.resize((size_t)qn + 1);

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
}

vector<ActiveEdge> &MatchingSolver::topEdgesBuffer(ui depth)
{
    if (top_edges_buffer_by_depth.size() <= depth) {
        top_edges_buffer_by_depth.resize((size_t)depth + 1);
    }
    vector<ActiveEdge> &buffer = top_edges_buffer_by_depth[depth];
    buffer.clear();
    if (buffer.capacity() < qn) {
        buffer.reserve(qn);
    }
    return buffer;
}

vector<ui> &MatchingSolver::whiteNbrsBuffer(ui depth)
{
    if (white_neighbors_buffer_by_depth.size() <= depth) {
        white_neighbors_buffer_by_depth.resize((size_t)depth + 1);
    }
    vector<ui> &buffer = white_neighbors_buffer_by_depth[depth];
    buffer.clear();
    if (buffer.capacity() < qn) {
        buffer.reserve(qn);
    }
    return buffer;
}

} // namespace cde_black_white
