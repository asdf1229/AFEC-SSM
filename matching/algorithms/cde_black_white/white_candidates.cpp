#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

bool MatchingSolver::calcBlackDelta(const SearchState &state, ui u, ui v, ui cost, ui &delta)
{
    delta = 0;
    for (ui neighbor : q_neighbors[u]) {
        if (!isBlack(state, neighbor)) continue;

        ui mapped_neighbor = (ui)state.mapped_q[neighbor];
        bool adjacent = anchorAdjacent(neighbor, mapped_neighbor, u, v);
        EdgeState state_uv = getEdge(state, u, neighbor);
        if (state_uv == EDGE_PRESENT) {
            if (!adjacent) return false;
        }
        else if (state_uv == EDGE_MISSING) {
            if (adjacent) return false;
        }
        else if (!adjacent) { // EDGE_UNDECIDED
            if (isQueryBridgeEdge(u, neighbor)) return false;
            delta++;
            if (cost + delta > threshold) {
                return false;
            }
        }
    }
    return true;
}

// Collect candidate ranges forced by black neighbors whose query edges must be present.
bool MatchingSolver::collectRequiredRanges(const SearchState &state, ui u, vector<pair<size_t, ui>> &ranges)
{
    ranges.clear();
    for (ui neighbor : q_neighbors[u]) {
        if (!isBlack(state, neighbor)) continue;

        EdgeState state_uv = getEdge(state, u, neighbor);
        if (state_uv == EDGE_MISSING) continue;
        if (state_uv == EDGE_UNDECIDED && !isQueryBridgeEdge(u, neighbor)) continue;
        // EDGE_PRESENT or EDGE_UNDECIDED with query bridge edge

        const pair<size_t, ui> *range = findAdjRange(neighbor, (ui)state.mapped_q[neighbor], u);
        if (range == nullptr) {
            stats.candidate_range_misses++;
            return false;
        }
        stats.candidate_range_hits++;
        ranges.push_back(*range);
    }
    return true;
}

// intersection of multiple candidate ranges, result stored in source
void MatchingSolver::intersectRequiredRanges(vector<pair<size_t, ui>> &ranges, vector<ui> &source)
{
    source.clear();
    assert(!ranges.empty());

    if (ranges.size() == 1) {
        const pair<size_t, ui> &range = ranges.front();
        source.assign(rangeBegin(range), rangeEnd(range));
        return;
    }

    std::sort(ranges.begin(), ranges.end(),
        [](const pair<size_t, ui> &lhs, const pair<size_t, ui> &rhs) {
            return lhs.second < rhs.second;
        });

    size_t total_len = 0;
    for (const pair<size_t, ui> &range : ranges) total_len += range.second;
    size_t edge_check_threshold = (size_t)ranges.front().second * ranges.size() * 8;

    // hybrid intersection strategy
    if (total_len <= edge_check_threshold) {
        stats.candidate_intersection_calls++;
        const pair<size_t, ui> &first = ranges.front();
        source.assign(rangeBegin(first), rangeEnd(first));

        for (size_t i = 1; i < ranges.size() && !source.empty(); ++i) {
            candidate_intersection_buffer.clear();
            const pair<size_t, ui> &range = ranges[i];
            std::set_intersection(source.begin(), source.end(),
                                rangeBegin(range), rangeEnd(range),
                            std::back_inserter(candidate_intersection_buffer));
            source.swap(candidate_intersection_buffer);
        }
        return;
    }

    // total_len > edge_check_threshold
    const pair<size_t, ui> &shortest = ranges.front();
    for (const ui *it = rangeBegin(shortest); it != rangeEnd(shortest); ++it) {
        ui candidate = *it;
        bool supported = true;
        for (size_t i = 1; i < ranges.size(); ++i) {
            stats.candidate_edge_check_calls++;
            if (!rangeHas(ranges[i], candidate)) {
                supported = false;
                break;
            }
        }
        if (supported) source.push_back(candidate);
    }
}

ui MatchingSolver::nextBatchToken()
{
    if (candidate_batch_mark.size() < gn) {
        candidate_batch_mark.assign(gn, 0);
        candidate_batch_pos.assign(gn, 0);
        candidate_batch_token = 0;
    }

    candidate_batch_token++;
    if (candidate_batch_token == 0) {
        std::fill(candidate_batch_mark.begin(), candidate_batch_mark.end(), 0);
        candidate_batch_token = 1;
    }
    return candidate_batch_token;
}

void MatchingSolver::addRangeHits(const pair<size_t, ui> *range, ui token, vector<ui> &hits)
{
    if (range == nullptr) {
        return;
    }

    for (const ui *it = rangeBegin(*range); it != rangeEnd(*range); ++it) {
        ui candidate = *it;
        if (candidate >= candidate_batch_mark.size() || candidate_batch_mark[candidate] != token) {
            continue;
        }

        ui pos = candidate_batch_pos[candidate];
        if (pos < hits.size()) {
            hits[pos]++;
        }
    }
}

void MatchingSolver::invalidateRange(const pair<size_t, ui> *range, ui token)
{
    if (range == nullptr) {
        return;
    }

    for (const ui *it = rangeBegin(*range); it != rangeEnd(*range); ++it) {
        ui candidate = *it;
        if (candidate >= candidate_batch_mark.size() || candidate_batch_mark[candidate] != token) {
            continue;
        }

        ui pos = candidate_batch_pos[candidate];
        if (pos < candidate_batch_valid.size()) {
            candidate_batch_valid[pos] = 0;
        }
    }
}

void MatchingSolver::computeWhiteCandCosts(const SearchState &state, ui u, ui cost,
    const vector<ui> &source, vector<ui> &result, vector<ui> &result_deltas)
{
    if (source.empty()) return;

    ui token = nextBatchToken();
    size_t source_count = source.size();
    candidate_batch_present_hits.assign(source_count, 0);
    candidate_batch_undecided_hits.assign(source_count, 0);
    candidate_batch_valid.assign(source_count, 0);

    for (size_t i = 0; i < source_count; ++i) {
        ui candidate = source[i];
        assert(candidate < gn);
        candidate_batch_mark[candidate] = token;
        candidate_batch_pos[candidate] = (ui)i;
        candidate_batch_valid[i] = isDataVertexUsed(state, candidate) ? 0 : 1;
    }

    ui present_count = 0;
    ui undecided_count = 0;
    for (ui neighbor : q_neighbors[u]) {
        if (!isBlack(state, neighbor)) continue;

        ui mapped_neighbor = (ui)state.mapped_q[neighbor];
        EdgeState state_uv = getEdge(state, u, neighbor);
        const pair<size_t, ui> *range = findAdjRange(neighbor, mapped_neighbor, u);

        stats.candidate_edge_check_calls += (long long)source_count;
        if (range == nullptr) stats.candidate_range_misses++;
        else stats.candidate_range_hits++;

        if (state_uv == EDGE_PRESENT) {
            present_count++;
            addRangeHits(range, token, candidate_batch_present_hits);
        }
        else if (state_uv == EDGE_MISSING) {
            invalidateRange(range, token);
        }
        else if (isQueryBridgeEdge(u, neighbor)) {
            present_count++;
            addRangeHits(range, token, candidate_batch_present_hits);
        }
        else {
            undecided_count++;
            addRangeHits(range, token, candidate_batch_undecided_hits);
        }
    }

    for (size_t i = 0; i < source_count; ++i) {
        if (!candidate_batch_valid[i]) continue;
        if (candidate_batch_present_hits[i] != present_count) continue;

        ui adjacent_undecided = candidate_batch_undecided_hits[i];
        assert(undecided_count >= adjacent_undecided);
        ui delta = undecided_count - adjacent_undecided;
        if (cost + delta <= threshold) {
            result.push_back(source[i]);
            result_deltas.push_back(delta);
        }
    }
}

bool MatchingSolver::buildWhiteCands(SearchState &state, ui u, ui cost)
{
    assert(cost <= threshold);
    candidate_result_buffer.clear();
    candidate_result_delta_buffer.clear();
    stats.white_bucket_rebuilds++;

    if (!collectRequiredRanges(state, u, candidate_range_buffer)) return false;
    assert(!candidate_range_buffer.empty());

    intersectRequiredRanges(candidate_range_buffer, candidate_source_buffer);
    if(candidate_source_buffer.empty()) return false;

    computeWhiteCandCosts(state, u, cost, candidate_source_buffer, candidate_result_buffer, candidate_result_delta_buffer);
    assert(candidate_result_buffer.size() == candidate_result_delta_buffer.size());
    return !candidate_result_buffer.empty();
}

bool MatchingSolver::refreshWhiteByBlack(SearchState &state, ui white_u, ui black_u, ui black_v, ui cost)
{
    assert(isWhite(state, white_u));
    assert(isBlack(state, black_u));
    assert(state.mapped_q[black_u] == (int)black_v);
    assert(cost <= threshold);

    WhiteCands old_bucket = state.white[white_u];
    assert(old_bucket.begin + old_bucket.count <= state.white_candidate_pool.size());
    assert(old_bucket.begin + old_bucket.count <= state.white_candidate_delta_pool.size());

    EdgeState state_uv = getEdge(state, white_u, black_u);
    assert(state_uv == EDGE_MISSING || state_uv == EDGE_UNDECIDED);
    const pair<size_t, ui> *range = findAdjRange(black_u, black_v, white_u);

    stats.candidate_edge_check_calls += (long long)old_bucket.count;
    if (range == nullptr) stats.candidate_range_misses++;
    else stats.candidate_range_hits++;

    ui remaining_budget = threshold - cost;
    candidate_result_buffer.clear();
    candidate_result_delta_buffer.clear();

    for (ui candidate_idx = 0; candidate_idx < old_bucket.count; ++candidate_idx) {
        size_t pool_idx = old_bucket.begin + candidate_idx;
        ui candidate = state.white_candidate_pool[pool_idx];
        if (isDataVertexUsed(state, candidate)) continue;

        bool adjacent = (range != nullptr) && rangeHas(*range, candidate);
        ui delta = state.white_candidate_delta_pool[pool_idx];
        bool keep = true;

        if(state_uv == EDGE_MISSING) {
            keep = !adjacent;
        }
        else if (!adjacent) { // EDGE_UNDECIDED
            if (isQueryBridgeEdge(white_u, black_u)) keep = false;
            else delta++;
        }

        if (!keep || delta > remaining_budget) continue;
        candidate_result_buffer.push_back(candidate);
        candidate_result_delta_buffer.push_back(delta);
    }

    if (candidate_result_buffer.empty()) return false;
    replaceBucket(state, white_u, candidate_result_buffer, candidate_result_delta_buffer);
    return true;
}

} // namespace cde_black_white
