#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

bool MatchingSolver::calcBlackDelta(const SearchState &state, ui u, ui v, ui cost, ui &delta)
{
    delta = 0;
    for (ui neighbor_idx = 0; neighbor_idx < q_degree[u]; ++neighbor_idx) {
        ui neighbor = q_neighbors[u][neighbor_idx];
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
            if (q_neighbor_is_bridge[u][neighbor_idx]) return false;
            delta++;
            if (cost + delta > threshold) {
                return false;
            }
        }
    }
    return true;
}

bool MatchingSolver::collectPosRanges(const SearchState &state, ui u,
    vector<pair<size_t, ui>> &ranges)
{
    // 收集 u 相对所有已确定存在 black 边的正向候选范围。
    ranges.clear();
    for (ui neighbor : q_neighbors[u]) {
        if (!isBlack(state, neighbor) ||
            getEdge(state, u, neighbor) != EDGE_PRESENT) {
            continue;
        }

        const pair<size_t, ui> *range = findAdjRange(
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

void MatchingSolver::buildRangeSource(
    vector<pair<size_t, ui>> &ranges, vector<ui> &source)
{
    // 将多个正向候选范围合并成后续过滤使用的候选源。
    source.clear();
    if (ranges.empty()) {
        return;
    }

    std::sort(ranges.begin(), ranges.end(),
        [](const pair<size_t, ui> &lhs, const pair<size_t, ui> &rhs) {
            return lhs.second < rhs.second;
        });

    if (ranges.size() == 1) {
        const pair<size_t, ui> &range = ranges.front();
        source.assign(rangeBegin(range), rangeEnd(range));
        return;
    }

    size_t total_len = 0;
    for (const pair<size_t, ui> &range : ranges) {
        total_len += range.second;
    }

    size_t edge_check_threshold = (size_t)ranges.front().second *
        ranges.size() * 8;
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

    const pair<size_t, ui> &shortest = ranges.front();
    for (const ui *it = rangeBegin(shortest);
        it != rangeEnd(shortest); ++it) {
        ui candidate = *it;
        bool supported = true;
        for (size_t i = 1; i < ranges.size(); ++i) {
            stats.candidate_edge_check_calls++;
            if (!rangeHas(ranges[i], candidate)) {
                supported = false;
                break;
            }
        }
        if (supported) {
            source.push_back(candidate);
        }
    }
}

bool MatchingSolver::bucketHas(const SearchState &state,
    const WhiteCands &bucket, ui candidate) const
{
    // 判断 white bucket 中是否包含指定候选。
    assert(bucket.begin + bucket.count <= state.white_candidate_pool.size());
    const ui *begin = state.white_candidate_pool.data() + bucket.begin;
    const ui *end = begin + bucket.count;
    return std::binary_search(begin, end, candidate);
}

ui MatchingSolver::nextBatchToken()
{
    // 生成候选批处理标记 token，溢出时重置标记数组。
    if (candidate_batch_mark.size() < gn) {
        candidate_batch_mark.assign(gn, 0);
        candidate_batch_pos.assign(gn, 0);
        candidate_batch_token = 0;
    }

    candidate_batch_token++;
    if (candidate_batch_token == 0) {
        std::fill(candidate_batch_mark.begin(),
            candidate_batch_mark.end(), 0);
        candidate_batch_token = 1;
    }
    return candidate_batch_token;
}

void MatchingSolver::addRangeHits(const pair<size_t, ui> *range, ui token,
    vector<ui> &hits)
{
    // 对批处理候选统计一个候选范围内的命中次数。
    if (range == nullptr) {
        return;
    }

    for (const ui *it = rangeBegin(*range);
        it != rangeEnd(*range); ++it) {
        ui candidate = *it;
        if (candidate >= candidate_batch_mark.size() ||
            candidate_batch_mark[candidate] != token) {
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
    // 将批处理候选中落入缺失边范围的候选标记为无效。
    if (range == nullptr) {
        return;
    }

    for (const ui *it = rangeBegin(*range);
        it != rangeEnd(*range); ++it) {
        ui candidate = *it;
        if (candidate >= candidate_batch_mark.size() ||
            candidate_batch_mark[candidate] != token) {
            continue;
        }

        ui pos = candidate_batch_pos[candidate];
        if (pos < candidate_batch_valid.size()) {
            candidate_batch_valid[pos] = 0;
        }
    }
}

void MatchingSolver::addFeasibleBatch(const SearchState &state, ui u,
    ui cost, const vector<ui> &source, vector<ui> &result)
{
    // 批量检查候选源中哪些候选满足当前 black 邻居约束。
    if (source.empty()) {
        return;
    }

    ui token = nextBatchToken();
    size_t source_count = source.size();
    candidate_batch_present_hits.assign(source_count, 0);
    candidate_batch_undecided_hits.assign(source_count, 0);
    candidate_batch_valid.assign(source_count, 0);

    for (size_t i = 0; i < source_count; ++i) {
        ui candidate = source[i];
        if (candidate >= gn) {
            continue;
        }
        candidate_batch_mark[candidate] = token;
        candidate_batch_pos[candidate] = (ui)i;
        candidate_batch_valid[i] =
            isDataVertexUsed(state, candidate) ? 0 : 1;
    }

    ui present_count = 0;
    ui undecided_count = 0;
    for (ui neighbor_idx = 0; neighbor_idx < q_degree[u]; ++neighbor_idx) {
        ui neighbor = q_neighbors[u][neighbor_idx];
        if (!isBlack(state, neighbor)) {
            continue;
        }

        ui mapped_neighbor = (ui)state.mapped_q[neighbor];
        EdgeState state_uv = getEdge(state, u, neighbor);
        const pair<size_t, ui> *range =
            findAdjRange(neighbor, mapped_neighbor, u);

        stats.candidate_edge_check_calls += (long long)source_count;
        if (range == nullptr) {
            stats.candidate_range_misses++;
        }
        else {
            stats.candidate_range_hits++;
        }

        if (state_uv == EDGE_PRESENT) {
            present_count++;
            addRangeHits(range, token,
                candidate_batch_present_hits);
        }
        else if (state_uv == EDGE_MISSING) {
            invalidateRange(range, token);
        }
        else if (q_neighbor_is_bridge[u][neighbor_idx]) {
            present_count++;
            addRangeHits(range, token,
                candidate_batch_present_hits);
        }
        else {
            undecided_count++;
            addRangeHits(range, token,
                candidate_batch_undecided_hits);
        }
    }

    for (size_t i = 0; i < source_count; ++i) {
        if (!candidate_batch_valid[i]) {
            continue;
        }
        if (candidate_batch_present_hits[i] != present_count) {
            continue;
        }

        ui adjacent_undecided = candidate_batch_undecided_hits[i];
        ui delta = undecided_count > adjacent_undecided
            ? undecided_count - adjacent_undecided : 0;
        if (cost + delta <= threshold) {
            result.push_back(source[i]);
        }
    }
}

void MatchingSolver::copyBucketCands(const SearchState &state,
    const WhiteCands &bucket, vector<ui> &target) const
{
    // 将 white bucket 中保存的候选复制到目标缓冲区。
    assert(bucket.begin + bucket.count <= state.white_candidate_pool.size());
    target.assign(state.white_candidate_pool.begin() + bucket.begin,
        state.white_candidate_pool.begin() + bucket.begin + bucket.count);
}

void MatchingSolver::filterByBucket(const SearchState &state,
    const WhiteCands &bucket, const vector<ui> &source,
    vector<ui> &target) const
{
    // 用已有 white bucket 过滤候选源，保留仍在 bucket 中的候选。
    target.clear();
    for (ui candidate : source) {
        if (bucketHas(state, bucket, candidate)) {
            target.push_back(candidate);
        }
    }
}

void MatchingSolver::collectAllCands(ui u, vector<ui> &target)
{
    // 收集查询点 u 的全部静态候选。
    target.clear();
    for (int candidate : candidates[u]) {
        target.push_back((ui)candidate);
    }
}

bool MatchingSolver::buildWhiteCands(SearchState &state, ui u, ui cost,
    const WhiteCands *existing_bucket)
{
    assert(cost <= threshold);
    candidate_result_buffer.clear();
    stats.white_bucket_rebuilds++;

    if (!collectPosRanges(state, u, candidate_range_buffer)) return false;

    if (!candidate_range_buffer.empty()) {
        buildRangeSource(candidate_range_buffer, candidate_source_buffer);

        if (existing_bucket != nullptr &&
            existing_bucket->count <= candidate_source_buffer.size()) {
            copyBucketCands(state, *existing_bucket,
                candidate_intersection_buffer);
            addFeasibleBatch(state, u, cost,
                candidate_intersection_buffer, candidate_result_buffer);
        }
        else {
            const vector<ui> *source = &candidate_source_buffer;
            if (existing_bucket != nullptr) {
                filterByBucket(state, *existing_bucket,
                    candidate_source_buffer,
                    candidate_intersection_buffer);
                source = &candidate_intersection_buffer;
            }
            addFeasibleBatch(state, u, cost, *source,
                candidate_result_buffer);
        }
    }
    else if (existing_bucket != nullptr) {
        copyBucketCands(state, *existing_bucket,
            candidate_intersection_buffer);
        addFeasibleBatch(state, u, cost,
            candidate_intersection_buffer, candidate_result_buffer);
    }
    else {
        if (existing_bucket != nullptr) {
            copyBucketCands(state, *existing_bucket,
                candidate_intersection_buffer);
            addFeasibleBatch(state, u, cost,
                candidate_intersection_buffer, candidate_result_buffer);
        }
        else {
            collectAllCands(u, candidate_source_buffer);
            addFeasibleBatch(state, u, cost,
                candidate_source_buffer, candidate_result_buffer);
        }
    }

    return !candidate_result_buffer.empty();
}

bool MatchingSolver::refreshWhiteCands(SearchState &state,
    ui white_u, ui cost)
{
    // 基于当前状态重建已有 white 点的候选桶。
    WhiteCands old_bucket = state.white[white_u];
    if (!buildWhiteCands(state, white_u, cost, &old_bucket)) {
        return false;
    }

    replaceBucket(state, white_u, candidate_result_buffer);
    return true;
}

bool MatchingSolver::refreshWhiteByBlack(SearchState &state, ui white_u,
    ui black_u, ui black_v, ui cost)
{
    // 新增 black 邻居后刷新 white_u 的候选桶。
    assert(isWhite(state, white_u));
    assert(isBlack(state, black_u));
    (void)black_u;
    (void)black_v;
    if (cost > threshold) {
        return false;
    }

    return refreshWhiteCands(state, white_u, cost);
}

} // namespace cde_black_white
