#include "matching/algorithms/afee/afee.h"

namespace afee {

void MatchingSolver::resetBuffers()
{
    top_edges_buffer.clear();
    top_edges_buffer.resize((size_t)qn + 1);
    branch_cands_buffer.clear();
    branch_cands_buffer.resize((size_t)qn + 1);
    white_neighbors_buffer.clear();
    white_neighbors_buffer.resize((size_t)qn + 1);
    component_id_buffer.clear();
    component_id_buffer.reserve(qn);
    component_queue_buffer.clear();
    component_queue_buffer.reserve(qn);
    component_edge_counts_buffer.clear();
    component_edge_counts_buffer.reserve(qn);

    candidate_range_buffer.clear();
    candidate_source_buffer.clear();
    candidate_result_buffer.clear();
    candidate_result_delta_buffer.clear();
    candidate_intersection_buffer.clear();
    candidate_batch_mark.assign(gn, 0);
    candidate_batch_pos.assign(gn, 0);
    candidate_batch_present_hits.clear();
    candidate_batch_undecided_hits.clear();
    candidate_batch_valid.clear();
    candidate_batch_token = 0;
}

vector<AnchorEdge> &MatchingSolver::topEdgesBuffer(ui depth)
{
    if (top_edges_buffer.size() <= depth) top_edges_buffer.resize((size_t)depth + 1);
    vector<AnchorEdge> &buffer = top_edges_buffer[depth];
    buffer.clear();
    if (buffer.capacity() < qn) buffer.reserve(qn);
    return buffer;
}

vector<ui> &MatchingSolver::whiteNbrsBuffer(ui depth)
{
    if (white_neighbors_buffer.size() <= depth) white_neighbors_buffer.resize((size_t)depth + 1);
    vector<ui> &buffer = white_neighbors_buffer[depth];
    buffer.clear();
    if (buffer.capacity() < qn) buffer.reserve(qn);
    return buffer;
}

vector<pair<ui, ui>> &MatchingSolver::branchCandsBuffer(ui depth)
{
    if (branch_cands_buffer.size() <= depth) branch_cands_buffer.resize((size_t)depth + 1);
    vector<pair<ui, ui>> &buffer = branch_cands_buffer[depth];
    buffer.clear();
    return buffer;
}

} // namespace afee
