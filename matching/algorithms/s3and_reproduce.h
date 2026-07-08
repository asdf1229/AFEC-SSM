#ifndef MATCHING_ALGORITHMS_S3AND_REPRODUCE_H_
#define MATCHING_ALGORITHMS_S3AND_REPRODUCE_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"

#include <random>

using namespace std;

class S3ANDReproduceSolver {
public:
    S3ANDReproduceSolver()
        : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr),
          threshold(0), qn(0), gn(0), rng(42)
    {}

    bool init(const Graph *q, const Graph *g, ui match_threshold)
    {
        Timer t_init;
        query_graph = q;
        data_graph = g;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();

        resetState();

        if (qn == 0 || gn == 0 || qn > gn) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        buildGraphCaches();

        Timer t_index;
        constructIndexBalanced();
        stats.index_time = t_index.elapsed();

        Timer t_filter;
        bool ok = buildCandidates();
        stats.filter_time = t_filter.elapsed();
        stats.init_time = t_init.elapsed();
        return ok;
    }

    void match(vector<vector<pair<ui, ui> > > &results)
    {
        Timer t_search;
        results_ptr = &results;
        results_ptr->clear();

        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        CCSandR(0, 0);

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct S3ANDReproduceStats {
        long long total_time = 0;
        long long init_time = 0;
        long long filter_time = 0;
        long long index_time = 0;
        long long dfs_time = 0;
        unsigned long long entry_node_visits = 0;
        unsigned long long entry_node_prunes = 0;
        unsigned long long vertex_visits = 0;
        unsigned long long candidate_checks = 0;
        unsigned long long recursion_calls = 0;
        unsigned long long prune_calls = 0;
        unsigned long long final_checks = 0;
        ui filter_candidate_count = 0;
        ui index_entry_count = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
    } stats;

    void printStats() const
    {
        printf("\n--- S3AND Reproduce Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms\n", stats.init_time / 1000.0);
        printf("Filter Time:         %.4lf ms\n", stats.filter_time / 1000.0);
        printf("Index Time:          %.4lf ms\n", stats.index_time / 1000.0);
        printf("CCSandR Time:        %.4lf ms\n", stats.dfs_time / 1000.0);
        printf("Index Entries:       %u\n", stats.index_entry_count);
        printf("Entry Visits:        %llu\n", stats.entry_node_visits);
        printf("Entry Prunes:        %llu\n", stats.entry_node_prunes);
        printf("Leaf Vertex Visits:  %llu\n", stats.vertex_visits);
        printf("Candidate Checks:    %llu\n", stats.candidate_checks);
        printf("Candidates:          %u\n", stats.filter_candidate_count);
        printf("Recursion Calls:     %llu\n", stats.recursion_calls);
        printf("Prune Calls:         %llu\n", stats.prune_calls);
        printf("Final Checks:        %llu\n", stats.final_checks);
        printf("Results Found:       %zu\n", stats.result_count);
        if (stats.output_limit_reached) {
            printf("Output Limit:        reached (%d)\n", MATCH_OUTPUT_LIMIT);
        }
    }

private:
    struct AuxSynopsis {
        vector<LabelID> labels;          // S3AND BV under exact one-label semantics.
        vector<LabelID> neighbor_labels; // S3AND NBV under exact one-label semantics.
        ui nk = 0;                       // Maximum incident-edge capacity in this entry.
    };

    struct IndexEntry {
        bool leaf = false;
        ui vertex = 0;
        ui level = 0;
        vector<ui> children;
        vector<ui> query_vertices;       // The source-code Q list carried by an entry.
        AuxSynopsis aux;
    };

    struct HeapEntry {
        ui entry = 0;
        ui key = 0;
        ui level = 0;

        bool operator<(const HeapEntry &rhs) const
        {
            if (level != rhs.level) return level > rhs.level;
            if (key != rhs.key) return key < rhs.key;
            return entry > rhs.entry;
        }
    };

    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui> > > *results_ptr;
    ui threshold;
    ui qn;
    ui gn;

    vector<vector<ui> > q_neighbors;
    vector<ui> q_degree;
    vector<ui> g_degree;
    vector<vector<ui> > candidates;
    vector<IndexEntry> index_entries;
    vector<ui> root_entries;
    vector<ui> order;
    vector<ui> pivot_by_order;
    vector<vector<ui> > backward_neighbors;
    vector<int> mapped_q;
    vector<int> mapped_g;
    mt19937 rng;

    enum {
        kPartitionNumber = 16,
        kCostModelIterations = 5
    };

    void resetState()
    {
        stats = S3ANDReproduceStats();
        q_neighbors.clear();
        q_degree.clear();
        g_degree.clear();
        candidates.clear();
        index_entries.clear();
        root_entries.clear();
        order.clear();
        pivot_by_order.clear();
        backward_neighbors.clear();
        mapped_q.clear();
        mapped_g.clear();
        rng.seed(42);
    }

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

    void buildGraphCaches()
    {
        q_neighbors.assign(qn, vector<ui>());
        q_degree.assign(qn, 0);
        for (ui u = 0; u < qn; ++u) {
            ui deg = 0;
            const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_degree[u] = deg;
            q_neighbors[u].assign(nbrs, nbrs + deg);
        }

        g_degree.assign(gn, 0);
        for (ui v = 0; v < gn; ++v) {
            g_degree[v] = data_graph->getVertexDegree(v);
        }
    }

    static void appendAndUnique(vector<LabelID> &dst,
        const vector<LabelID> &src)
    {
        dst.insert(dst.end(), src.begin(), src.end());
        sort(dst.begin(), dst.end());
        dst.erase(unique(dst.begin(), dst.end()), dst.end());
    }

    static bool hasLabel(const vector<LabelID> &labels, LabelID label)
    {
        return binary_search(labels.begin(), labels.end(), label);
    }

    AuxSynopsis makeVertexAux(ui v) const
    {
        AuxSynopsis aux;
        aux.labels.push_back(data_graph->getVertexLabel(v));

        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors(v, deg);
        aux.neighbor_labels.reserve(deg);
        for (ui i = 0; i < deg; ++i) {
            aux.neighbor_labels.push_back(data_graph->getVertexLabel(nbrs[i]));
        }
        sort(aux.neighbor_labels.begin(), aux.neighbor_labels.end());
        aux.neighbor_labels.erase(unique(aux.neighbor_labels.begin(),
            aux.neighbor_labels.end()), aux.neighbor_labels.end());

        aux.nk = deg;
        return aux;
    }

    AuxSynopsis aggregateAux(const vector<ui> &children) const
    {
        AuxSynopsis aux;
        for (ui child_id : children) {
            const AuxSynopsis &child = index_entries[child_id].aux;
            appendAndUnique(aux.labels, child.labels);
            appendAndUnique(aux.neighbor_labels, child.neighbor_labels);
            aux.nk = max(aux.nk, child.nk);
        }
        return aux;
    }

    bool nodeFeatureLess(ui lhs, ui rhs) const
    {
        LabelID ll = data_graph->getVertexLabel(lhs);
        LabelID lr = data_graph->getVertexLabel(rhs);
        if (ll != lr) return ll < lr;

        ui dl = g_degree[lhs];
        ui dr = g_degree[rhs];
        if (dl != dr) return dl > dr;
        return lhs < rhs;
    }

    ui labelDistance(LabelID lhs, LabelID rhs) const
    {
        return lhs == rhs ? 0 : 2;
    }

    vector<vector<ui> > initializePartition(const vector<ui> &nodes,
        ui partition_count)
    {
        vector<ui> shuffled = nodes;
        shuffle(shuffled.begin(), shuffled.end(), rng);

        vector<ui> centers(shuffled.begin(), shuffled.begin() + partition_count);
        vector<vector<ui> > partition(partition_count);
        for (ui i = 0; i < partition_count; ++i) {
            partition[i].push_back(centers[i]);
        }

        ui soft_cap = (ui)((nodes.size() + partition_count - 1) / partition_count);
        for (ui node : shuffled) {
            if (find(centers.begin(), centers.end(), node) != centers.end()) {
                continue;
            }

            ui best = partition_count;
            ui best_distance = numeric_limits<ui>::max();
            for (ui i = 0; i < partition_count; ++i) {
                if (partition[i].size() >= soft_cap) continue;
                ui dist = labelDistance(data_graph->getVertexLabel(node),
                    data_graph->getVertexLabel(centers[i]));
                dist = dist * 1000000u +
                    (g_degree[node] > g_degree[centers[i]]
                        ? g_degree[node] - g_degree[centers[i]]
                        : g_degree[centers[i]] - g_degree[node]);
                if (dist < best_distance ||
                    (dist == best_distance && partition[i].size() < partition[best].size())) {
                    best = i;
                    best_distance = dist;
                }
            }

            if (best == partition_count) {
                best = 0;
                for (ui i = 1; i < partition_count; ++i) {
                    if (partition[i].size() < partition[best].size()) best = i;
                }
            }
            partition[best].push_back(node);
        }

        return partition;
    }

    LabelID partitionCenterLabel(const vector<ui> &part) const
    {
        vector<LabelID> labels;
        labels.reserve(part.size());
        for (ui v : part) labels.push_back(data_graph->getVertexLabel(v));
        sort(labels.begin(), labels.end());

        LabelID best_label = labels[0];
        ui best_count = 0;
        for (size_t i = 0; i < labels.size();) {
            size_t j = i + 1;
            while (j < labels.size() && labels[j] == labels[i]) ++j;
            ui count = (ui)(j - i);
            if (count > best_count) {
                best_label = labels[i];
                best_count = count;
            }
            i = j;
        }
        return best_label;
    }

    unsigned long long partitionCost(const vector<vector<ui> > &partition,
        const vector<LabelID> &centers) const
    {
        unsigned long long cost = 0;
        for (ui i = 0; i < partition.size(); ++i) {
            for (ui v : partition[i]) {
                cost += labelDistance(data_graph->getVertexLabel(v), centers[i]);
            }
        }
        return cost;
    }

    vector<vector<ui> > costModel(vector<vector<ui> > partition,
        const vector<ui> &nodes, ui partition_count)
    {
        ui soft_cap = (ui)((nodes.size() + partition_count - 1) / partition_count);
        vector<LabelID> centers(partition_count, 0);
        for (ui i = 0; i < partition_count; ++i) {
            centers[i] = partitionCenterLabel(partition[i]);
        }
        unsigned long long best_cost = partitionCost(partition, centers);

        vector<ui> ordered = nodes;
        sort(ordered.begin(), ordered.end(),
            [this](ui lhs, ui rhs) { return nodeFeatureLess(lhs, rhs); });

        for (ui iter = 0; iter < kCostModelIterations; ++iter) {
            vector<vector<ui> > next(partition_count);
            for (ui v : ordered) {
                ui best = partition_count;
                ui best_distance = numeric_limits<ui>::max();
                for (ui i = 0; i < partition_count; ++i) {
                    if (next[i].size() >= soft_cap) continue;
                    ui dist = labelDistance(data_graph->getVertexLabel(v), centers[i]);
                    if (dist < best_distance ||
                        (dist == best_distance &&
                         (best == partition_count ||
                          next[i].size() < next[best].size()))) {
                        best = i;
                        best_distance = dist;
                    }
                }

                if (best == partition_count) {
                    best = 0;
                    for (ui i = 1; i < partition_count; ++i) {
                        if (next[i].size() < next[best].size()) best = i;
                    }
                }
                next[best].push_back(v);
            }

            bool has_empty = false;
            for (const auto &part : next) {
                if (part.empty()) {
                    has_empty = true;
                    break;
                }
            }
            if (has_empty) break;

            vector<LabelID> next_centers(partition_count, 0);
            for (ui i = 0; i < partition_count; ++i) {
                next_centers[i] = partitionCenterLabel(next[i]);
            }

            unsigned long long next_cost = partitionCost(next, next_centers);
            if (next_cost > best_cost && iter > 0) break;
            if (next == partition) break;

            partition.swap(next);
            centers.swap(next_centers);
            best_cost = next_cost;
        }

        return partition;
    }

    vector<ui> rootTree(const vector<ui> &nodes, ui level)
    {
        if (nodes.size() <= (size_t)kPartitionNumber) {
            vector<ui> leaves;
            leaves.reserve(nodes.size());
            for (ui v : nodes) {
                ui entry_id = (ui)index_entries.size();
                index_entries.push_back(IndexEntry());
                IndexEntry &entry = index_entries.back();
                entry.leaf = true;
                entry.vertex = v;
                entry.level = level;
                entry.aux = makeVertexAux(v);
                leaves.push_back(entry_id);
            }
            return leaves;
        }

        ui partition_count = min((ui)kPartitionNumber, (ui)nodes.size());
        vector<vector<ui> > partition = initializePartition(nodes, partition_count);
        partition = costModel(partition, nodes, partition_count);

        vector<ui> entries;
        entries.reserve(partition.size());
        for (const auto &part : partition) {
            if (part.empty()) continue;

            vector<ui> child_entries = rootTree(part, level + 1);
            ui entry_id = (ui)index_entries.size();
            index_entries.push_back(IndexEntry());
            IndexEntry &entry = index_entries.back();
            entry.leaf = false;
            entry.level = level;
            entry.children.swap(child_entries);
            entry.aux = aggregateAux(entry.children);
            entries.push_back(entry_id);
        }

        return entries;
    }

    void constructIndexBalanced()
    {
        vector<ui> vertices;
        vertices.reserve(gn);
        for (ui v = 0; v < gn; ++v) vertices.push_back(v);
        root_entries = rootTree(vertices, 0);
        stats.index_entry_count = (ui)index_entries.size();
    }

    bool pruningKeywordSet(const IndexEntry &entry, ui q) const
    {
        return hasLabel(entry.aux.labels, query_graph->getVertexLabel(q));
    }

    bool pruningLBNDPhase1(const IndexEntry &entry, ui q) const
    {
        ui lb_nd = q_degree[q] > entry.aux.nk ? q_degree[q] - entry.aux.nk : 0;
        return lb_nd <= threshold;
    }

    bool pruningLBNDPhase2(const IndexEntry &entry, ui q) const
    {
        ui lb_nd = q_degree[q];
        for (ui q_nbr : q_neighbors[q]) {
            if (hasLabel(entry.aux.neighbor_labels,
                    query_graph->getVertexLabel(q_nbr))) {
                --lb_nd;
            }
        }
        return lb_nd <= threshold;
    }

    bool entryMayMatchQueryVertex(const IndexEntry &entry, ui q) const
    {
        return pruningKeywordSet(entry, q) &&
            pruningLBNDPhase1(entry, q) &&
            pruningLBNDPhase2(entry, q);
    }

    void retrieveCandidatesFromIndex()
    {
        candidates.assign(qn, vector<ui>());
        priority_queue<HeapEntry> heap;

        for (ui entry_id : root_entries) {
            IndexEntry &entry = index_entries[entry_id];
            entry.query_vertices.clear();
            for (ui q = 0; q < qn; ++q) {
                if (entryMayMatchQueryVertex(entry, q)) {
                    entry.query_vertices.push_back(q);
                }
            }
            if (entry.query_vertices.empty()) {
                stats.entry_node_prunes++;
            }
            else {
                heap.push({ entry_id, entry.aux.nk, entry.level });
            }
        }

        while (!heap.empty()) {
            HeapEntry heap_entry = heap.top();
            heap.pop();

            IndexEntry &entry = index_entries[heap_entry.entry];
            stats.entry_node_visits++;

            if (entry.leaf) {
                stats.vertex_visits++;
                for (ui q : entry.query_vertices) {
                    stats.candidate_checks++;
                    if (query_graph->getVertexLabel(q) !=
                        data_graph->getVertexLabel(entry.vertex)) {
                        continue;
                    }
                    if (!entryMayMatchQueryVertex(entry, q)) continue;
                    candidates[q].push_back(entry.vertex);
                }
                continue;
            }

            for (ui child_id : entry.children) {
                IndexEntry &child = index_entries[child_id];
                child.query_vertices.clear();
                for (ui q : entry.query_vertices) {
                    if (entryMayMatchQueryVertex(child, q)) {
                        child.query_vertices.push_back(q);
                    }
                }
                if (child.query_vertices.empty()) {
                    stats.entry_node_prunes++;
                }
                else {
                    heap.push({ child_id, child.aux.nk, child.level });
                }
            }
        }
    }

    bool buildCandidates()
    {
        retrieveCandidatesFromIndex();
        stats.filter_candidate_count = 0;

        for (ui q = 0; q < qn; ++q) {
            vector<ui> &cand = candidates[q];
            sort(cand.begin(), cand.end());
            cand.erase(unique(cand.begin(), cand.end()), cand.end());
            if (cand.empty()) return false;
            stats.filter_candidate_count += (ui)cand.size();
        }

        selectQueryOrderWithPivot();
        generateBackwardNeighbors();
        return true;
    }

    ui betterOrderVertex(ui lhs, ui rhs) const
    {
        if (rhs == qn) return lhs;
        if (candidates[lhs].size() != candidates[rhs].size()) {
            return candidates[lhs].size() < candidates[rhs].size() ? lhs : rhs;
        }
        if (q_degree[lhs] != q_degree[rhs]) {
            return q_degree[lhs] > q_degree[rhs] ? lhs : rhs;
        }
        return lhs < rhs ? lhs : rhs;
    }

    void selectQueryOrderWithPivot()
    {
        order.clear();
        pivot_by_order.clear();
        vector<char> processed(qn, 0);

        ui start = 0;
        for (ui q = 1; q < qn; ++q) {
            start = betterOrderVertex(q, start);
        }
        order.push_back(start);
        pivot_by_order.push_back(qn);
        processed[start] = 1;

        while (order.size() < qn) {
            ui next = qn;
            ui pivot = qn;

            for (ui prev : order) {
                for (ui nbr : q_neighbors[prev]) {
                    if (processed[nbr]) continue;
                    ui old_next = next;
                    ui chosen = betterOrderVertex(nbr, next);
                    if (chosen == nbr && old_next != nbr) {
                        next = nbr;
                        pivot = prev;
                    }
                }
            }

            if (next == qn) {
                for (ui q = 0; q < qn; ++q) {
                    if (processed[q]) continue;
                    ui old_next = next;
                    next = betterOrderVertex(q, next);
                    if (next == q && old_next != q) pivot = qn;
                }
            }

            assert(next != qn);
            order.push_back(next);
            pivot_by_order.push_back(pivot);
            processed[next] = 1;
        }
    }

    void generateBackwardNeighbors()
    {
        backward_neighbors.assign(qn, vector<ui>());
        vector<char> visited(qn, 0);
        visited[order[0]] = 1;

        for (ui i = 1; i < order.size(); ++i) {
            ui u = order[i];
            ui pivot = pivot_by_order[i];
            for (ui nbr : q_neighbors[u]) {
                if (visited[nbr] && nbr != pivot) {
                    backward_neighbors[u].push_back(nbr);
                }
            }
            visited[u] = 1;
        }
    }

    bool candidateSatisfiesBackwardNeighbors(ui u, ui v) const
    {
        for (ui u_nbr : backward_neighbors[u]) {
            int mapped_v = mapped_q[u_nbr];
            assert(mapped_v != -1);
            if (!data_graph->hasEdge(v, (ui)mapped_v)) return false;
        }
        return true;
    }

    ui missingEdgesToMappedNeighbors(ui u, ui v, ui limit) const
    {
        ui missing = 0;
        for (ui u_nbr : q_neighbors[u]) {
            int mapped_v = mapped_q[u_nbr];
            if (mapped_v == -1) continue;

            if (!data_graph->hasEdge(v, (ui)mapped_v)) {
                ++missing;
                if (missing > limit) return missing;
            }
        }
        return missing;
    }

    ui countMappedQueryNeighbors(ui u) const
    {
        ui count = 0;
        for (ui u_nbr : q_neighbors[u]) {
            if (mapped_q[u_nbr] != -1) ++count;
        }
        return count;
    }

    ui minMissingToMappedNeighbors(ui u, ui limit) const
    {
        ui best = limit + 1;
        for (ui v : candidates[u]) {
            if (mapped_g[v] != -1) continue;
            ui missing = missingEdgesToMappedNeighbors(u, v, best == 0 ? 0 : best - 1);
            if (missing < best) {
                best = missing;
                if (best == 0) return 0;
            }
        }
        return best;
    }

    ui boundaryMissingLowerBound(ui cost, ui limit) const
    {
        ui lower_bound = 0;
        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] != -1) continue;
            if (countMappedQueryNeighbors(u) == 0) continue;

            ui remaining = limit >= lower_bound ? limit - lower_bound : 0;
            ui min_missing = minMissingToMappedNeighbors(u, remaining);
            if (min_missing > remaining) return limit + 1;
            lower_bound += min_missing;
            if (lower_bound > limit || cost + lower_bound > threshold) {
                return lower_bound;
            }
        }
        return lower_bound;
    }

    bool partialFeasibilityCheck(ui cost) const
    {
        if (cost > threshold) return false;
        ui remaining_budget = threshold - cost;
        ui lower_bound = boundaryMissingLowerBound(cost, remaining_budget);
        return lower_bound <= remaining_budget && cost + lower_bound <= threshold;
    }

    ui missingEdgesInMappedQuery(ui limit) const
    {
        ui missing = 0;
        for (ui u = 0; u < qn; ++u) {
            for (ui nbr : q_neighbors[u]) {
                if (u >= nbr) continue;
                assert(mapped_q[u] != -1 && mapped_q[nbr] != -1);
                if (!data_graph->hasEdge((ui)mapped_q[u], (ui)mapped_q[nbr])) {
                    ++missing;
                    if (missing > limit) return missing;
                }
            }
        }
        return missing;
    }

    bool mappedImageIsConnected() const
    {
        if (qn <= 1) return true;

        vector<char> visited(qn, 0);
        queue<ui> bfs;
        bfs.push(0);
        ui visited_count = 0;

        while (!bfs.empty()) {
            ui u = bfs.front();
            bfs.pop();
            if (visited[u]) continue;

            visited[u] = 1;
            ++visited_count;

            for (ui nbr : q_neighbors[u]) {
                if (visited[nbr]) continue;
                assert(mapped_q[u] != -1 && mapped_q[nbr] != -1);
                if (data_graph->hasEdge((ui)mapped_q[u], (ui)mapped_q[nbr])) {
                    bfs.push(nbr);
                }
            }
        }

        return visited_count == qn;
    }

    bool ANDConstraintSUM()
    {
        stats.final_checks++;
        if (missingEdgesInMappedQuery(threshold) > threshold) return false;
        return mappedImageIsConnected();
    }

    void emitResult()
    {
        stats.result_count++;
#ifndef NDEBUG
        vector<pair<ui, ui> > result;
        result.reserve(qn);
        for (ui q = 0; q < qn; ++q) {
            result.push_back({ q, (ui)mapped_q[q] });
        }
        results_ptr->push_back(result);
#endif
        noteOutputLimitIfReached();
    }

    void CCSandR(ui depth, ui cost)
    {
        if (outputLimitReached()) return;
        stats.recursion_calls++;

        if (cost > threshold) {
            stats.prune_calls++;
            return;
        }

        if (depth == qn) {
            if (ANDConstraintSUM()) {
                emitResult();
            }
            else {
                stats.prune_calls++;
            }
            return;
        }

        ui u = order[depth];
        bool branched = false;
        for (ui v : candidates[u]) {
            if (mapped_g[v] != -1) continue;
            if (depth > 0 && !candidateSatisfiesBackwardNeighbors(u, v)) continue;

            ui remaining_budget = threshold - cost;
            ui delta = missingEdgesToMappedNeighbors(u, v, remaining_budget);
            if (delta > remaining_budget) continue;

            branched = true;
            mapped_q[u] = (int)v;
            mapped_g[v] = (int)u;

            if (partialFeasibilityCheck(cost + delta)) {
                CCSandR(depth + 1, cost + delta);
            }
            else {
                stats.prune_calls++;
            }

            mapped_g[v] = -1;
            mapped_q[u] = -1;

            if (outputLimitReached()) break;
        }

        if (!branched) {
            stats.prune_calls++;
        }
    }
};

void Approximate_S3AND_Reproduce(const Graph *query_graph, const Graph *data_graph,
    vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    S3ANDReproduceSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
