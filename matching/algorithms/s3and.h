#ifndef MATCHING_ALGORITHMS_S3AND_H_
#define MATCHING_ALGORITHMS_S3AND_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"

using namespace std;

class S3ANDSolver {
public:
    S3ANDSolver()
        : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr),
          threshold(0), qn(0), gn(0), branch_stamp(1), index_root(0)
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
        buildIndex();
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

        ui root = query_plan.empty() ? selectInitialRoot() : query_plan[0];
        for (ui v : candidates[root]) {
            if (outputLimitReached()) break;
            if (!candidateCanStillConnect(root, v)) {
                stats.connectivity_prunes++;
                continue;
            }

            mapped_q[root] = (int)v;
            mapped_g[v] = (int)root;
            if (partialFeasibilityCheck(0)) {
                dfs(1, 0);
            }
            else {
                stats.prune_calls++;
            }
            mapped_g[v] = -1;
            mapped_q[root] = -1;
        }

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct S3ANDStats {
        long long total_time = 0;
        long long init_time = 0;
        long long filter_time = 0;
        long long index_time = 0;
        long long dfs_time = 0;
        unsigned long long recursion_calls = 0;
        unsigned long long prune_calls = 0;
        unsigned long long index_node_visits = 0;
        unsigned long long index_node_prunes = 0;
        unsigned long long candidate_checks = 0;
        unsigned long long forced_lower_bound_prunes = 0;
        unsigned long long connectivity_prunes = 0;
        ui filter_candidate_count = 0;
        ui index_node_count = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
    } stats;

    void printStats() const
    {
        printf("\n--- S3AND Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms\n", stats.init_time / 1000.0);
        printf("Filter Time:         %.4lf ms\n", stats.filter_time / 1000.0);
        printf("Index Time:          %.4lf ms\n", stats.index_time / 1000.0);
        printf("DFS Time:            %.4lf ms\n", stats.dfs_time / 1000.0);
        printf("Index Nodes:         %u\n", stats.index_node_count);
        printf("Index Node Visits:   %llu\n", stats.index_node_visits);
        printf("Index Node Prunes:   %llu\n", stats.index_node_prunes);
        printf("Candidate Checks:    %llu\n", stats.candidate_checks);
        printf("Candidates:          %u\n", stats.filter_candidate_count);
        printf("Recursion Calls:     %llu\n", stats.recursion_calls);
        printf("Prune Calls:         %llu\n", stats.prune_calls);
        printf("LB Prunes:           %llu\n", stats.forced_lower_bound_prunes);
        printf("Connectivity Prunes: %llu\n", stats.connectivity_prunes);
        printf("Results Found:       %zu\n", stats.result_count);
        if (stats.output_limit_reached) {
            printf("Output Limit:        reached (%d)\n", MATCH_OUTPUT_LIMIT);
        }
    }

private:
    struct LabelCount {
        LabelID label;
        ui count;
    };

    struct IndexNode {
        bool leaf = false;
        ui max_degree = 0;
        vector<ui> vertices;
        vector<ui> children;
        vector<LabelID> labels;
        vector<LabelCount> max_neighbor_label_counts;
    };

    struct HeapEntry {
        ui node;
        ui key;

        bool operator<(const HeapEntry &rhs) const
        {
            if (key != rhs.key) return key < rhs.key;
            return node > rhs.node;
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
    vector<vector<LabelCount> > Lq_counts;
    vector<vector<LabelCount> > Lg_counts;
    vector<vector<ui> > candidates;
    vector<vector<char> > candidate_marks;
    vector<vector<char> > candidate_neighbor_marks;
    vector<IndexNode> index_nodes;
    vector<ui> query_plan;
    vector<int> branch_seen_stamp;
    vector<unsigned short> branch_hit_count;
    vector<ui> branch_touched;
    int branch_stamp;
    vector<int> mapped_q;
    vector<int> mapped_g;
    ui index_root;

    enum {
        kIndexFanout = 16,
        kIndexLeafSize = 512
    };

    void resetState()
    {
        stats = S3ANDStats();
        q_neighbors.clear();
        q_degree.clear();
        g_degree.clear();
        Lq_counts.clear();
        Lg_counts.clear();
        candidates.clear();
        candidate_marks.clear();
        candidate_neighbor_marks.clear();
        index_nodes.clear();
        query_plan.clear();
        branch_seen_stamp.clear();
        branch_hit_count.clear();
        branch_touched.clear();
        branch_stamp = 1;
        mapped_q.clear();
        mapped_g.clear();
        index_root = 0;
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

        buildLabelCounts(query_graph, qn, Lq_counts);
        buildLabelCounts(data_graph, gn, Lg_counts);
    }

    void buildLabelCounts(const Graph *graph, ui n,
        vector<vector<LabelCount> > &counts) const
    {
        counts.assign(n, vector<LabelCount>());
        vector<LabelID> labels;

        for (ui u = 0; u < n; ++u) {
            ui deg = 0;
            const ui *nbrs = graph->getVertexNeighbors(u, deg);
            labels.clear();
            labels.reserve(deg);

            for (ui i = 0; i < deg; ++i) {
                labels.push_back(graph->getVertexLabel(nbrs[i]));
            }
            sort(labels.begin(), labels.end());

            vector<LabelCount> &out = counts[u];
            out.reserve(labels.size());
            for (size_t i = 0; i < labels.size();) {
                size_t j = i + 1;
                while (j < labels.size() && labels[j] == labels[i]) ++j;
                out.push_back({ labels[i], (ui)(j - i) });
                i = j;
            }
        }
    }

    ui computeNeighborLabelDeficit(ui u, ui v) const
    {
        const vector<LabelCount> &qc = Lq_counts[u];
        const vector<LabelCount> &gc = Lg_counts[v];
        size_t gi = 0;
        ui deficit = 0;

        for (const auto &q_count : qc) {
            while (gi < gc.size() && gc[gi].label < q_count.label) {
                ++gi;
            }

            ui data_count = 0;
            if (gi < gc.size() && gc[gi].label == q_count.label) {
                data_count = gc[gi].count;
            }

            if (q_count.count > data_count) {
                deficit += q_count.count - data_count;
                if (deficit > threshold) return deficit;
            }
        }

        return deficit;
    }

    bool hasLabel(const vector<LabelID> &labels, LabelID label) const
    {
        return binary_search(labels.begin(), labels.end(), label);
    }

    bool nodeFeatureLess(ui lhs, ui rhs) const
    {
        LabelID ll = data_graph->getVertexLabel(lhs);
        LabelID lr = data_graph->getVertexLabel(rhs);
        if (ll != lr) return ll < lr;

        ui dl = g_degree[lhs];
        ui dr = g_degree[rhs];
        if (dl != dr) return dl > dr;

        const vector<LabelCount> &lc = Lg_counts[lhs];
        const vector<LabelCount> &rc = Lg_counts[rhs];
        size_t n = min(lc.size(), rc.size());
        for (size_t i = 0; i < n; ++i) {
            if (lc[i].label != rc[i].label) return lc[i].label < rc[i].label;
            if (lc[i].count != rc[i].count) return lc[i].count > rc[i].count;
        }
        if (lc.size() != rc.size()) return lc.size() > rc.size();
        return lhs < rhs;
    }

    void computeNodeAggregates(IndexNode &node, const vector<ui> &vertices)
    {
        node.max_degree = 0;
        node.labels.clear();
        node.max_neighbor_label_counts.clear();

        map<LabelID, ui> max_counts;
        node.labels.reserve(vertices.size());

        for (ui v : vertices) {
            node.max_degree = max(node.max_degree, g_degree[v]);
            node.labels.push_back(data_graph->getVertexLabel(v));

            const vector<LabelCount> &counts = Lg_counts[v];
            for (const auto &item : counts) {
                ui &current = max_counts[item.label];
                if (item.count > current) current = item.count;
            }
        }

        sort(node.labels.begin(), node.labels.end());
        node.labels.erase(unique(node.labels.begin(), node.labels.end()), node.labels.end());

        node.max_neighbor_label_counts.reserve(max_counts.size());
        for (const auto &item : max_counts) {
            node.max_neighbor_label_counts.push_back({ item.first, item.second });
        }
    }

    ui buildIndexNode(const vector<ui> &vertices)
    {
        ui node_id = (ui)index_nodes.size();
        index_nodes.push_back(IndexNode());
        computeNodeAggregates(index_nodes[node_id], vertices);

        if (vertices.size() <= (size_t)kIndexLeafSize) {
            index_nodes[node_id].leaf = true;
            index_nodes[node_id].vertices = vertices;
            return node_id;
        }

        vector<ui> ordered = vertices;
        sort(ordered.begin(), ordered.end(),
            [this](ui lhs, ui rhs) { return nodeFeatureLess(lhs, rhs); });

        size_t part_count = min((size_t)kIndexFanout,
            (ordered.size() + (size_t)kIndexLeafSize - 1) / (size_t)kIndexLeafSize);
        if (part_count < 2) part_count = 2;

        index_nodes[node_id].leaf = false;
        index_nodes[node_id].children.reserve(part_count);
        for (size_t p = 0; p < part_count; ++p) {
            size_t begin = (ordered.size() * p) / part_count;
            size_t end = (ordered.size() * (p + 1)) / part_count;
            if (begin >= end) continue;

            vector<ui> child_vertices(ordered.begin() + begin, ordered.begin() + end);
            ui child_id = buildIndexNode(child_vertices);
            index_nodes[node_id].children.push_back(child_id);
        }

        return node_id;
    }

    void buildIndex()
    {
        index_nodes.clear();
        vector<ui> vertices;
        vertices.reserve(gn);
        for (ui v = 0; v < gn; ++v) {
            vertices.push_back(v);
        }

        if (!vertices.empty()) {
            index_root = buildIndexNode(vertices);
        }
        stats.index_node_count = (ui)index_nodes.size();
    }

    ui computeNodeNeighborLabelDeficit(ui u, const IndexNode &node) const
    {
        const vector<LabelCount> &qc = Lq_counts[u];
        const vector<LabelCount> &nc = node.max_neighbor_label_counts;
        size_t ni = 0;
        ui deficit = 0;

        for (const auto &q_count : qc) {
            while (ni < nc.size() && nc[ni].label < q_count.label) {
                ++ni;
            }

            ui node_count = 0;
            if (ni < nc.size() && nc[ni].label == q_count.label) {
                node_count = nc[ni].count;
            }

            if (q_count.count > node_count) {
                deficit += q_count.count - node_count;
                if (deficit > threshold) return deficit;
            }
        }

        return deficit;
    }

    bool nodeMayContainCandidate(ui u, const IndexNode &node) const
    {
        if (!hasLabel(node.labels, query_graph->getVertexLabel(u))) return false;
        if (q_degree[u] > node.max_degree + threshold) return false;
        if (computeNodeNeighborLabelDeficit(u, node) > threshold) return false;
        return true;
    }

    ui nodeTraversalKey(ui u, const IndexNode &node) const
    {
        ui deficit = computeNodeNeighborLabelDeficit(u, node);
        ui coverage = q_degree[u] > deficit ? q_degree[u] - deficit : 0;
        return coverage * 1000000u + node.max_degree;
    }

    void retrieveCandidatesFromIndex(ui u, vector<ui> &cand)
    {
        cand.clear();
        if (index_nodes.empty()) return;

        priority_queue<HeapEntry> heap;
        const IndexNode &root = index_nodes[index_root];
        if (nodeMayContainCandidate(u, root)) {
            heap.push({ index_root, nodeTraversalKey(u, root) });
        }
        else {
            stats.index_node_prunes++;
        }

        while (!heap.empty()) {
            HeapEntry entry = heap.top();
            heap.pop();
            stats.index_node_visits++;

            const IndexNode &node = index_nodes[entry.node];
            if (!nodeMayContainCandidate(u, node)) {
                stats.index_node_prunes++;
                continue;
            }

            if (node.leaf) {
                for (ui v : node.vertices) {
                    stats.candidate_checks++;
                    if (query_graph->getVertexLabel(u) != data_graph->getVertexLabel(v)) continue;
                    if (q_degree[u] > g_degree[v] + threshold) continue;
                    if (computeNeighborLabelDeficit(u, v) > threshold) continue;
                    cand.push_back(v);
                }
            }
            else {
                for (ui child_id : node.children) {
                    const IndexNode &child = index_nodes[child_id];
                    if (nodeMayContainCandidate(u, child)) {
                        heap.push({ child_id, nodeTraversalKey(u, child) });
                    }
                    else {
                        stats.index_node_prunes++;
                    }
                }
            }
        }
    }

    bool buildCandidates()
    {
        candidates.assign(qn, vector<ui>());
        candidate_marks.assign(qn, vector<char>(gn, 0));
        stats.filter_candidate_count = 0;

        for (ui u = 0; u < qn; ++u) {
            vector<ui> &cand = candidates[u];
            retrieveCandidatesFromIndex(u, cand);
            sort(cand.begin(), cand.end());
            cand.erase(unique(cand.begin(), cand.end()), cand.end());

            if (cand.empty()) {
                return false;
            }
            for (ui v : cand) {
                candidate_marks[u][v] = 1;
            }
            stats.filter_candidate_count += (ui)cand.size();
        }

        buildQueryPlan();
        buildCandidateNeighborMarks();
        branch_seen_stamp.assign(gn, 0);
        branch_hit_count.assign(gn, 0);
        branch_touched.clear();
        branch_touched.reserve(data_graph->getMaxDegree());
        branch_stamp = 1;
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        return true;
    }

    void buildCandidateNeighborMarks()
    {
        candidate_neighbor_marks.assign(qn, vector<char>(gn, 0));
        for (ui u = 0; u < qn; ++u) {
            for (ui v : candidates[u]) {
                ui deg = 0;
                const ui *nbrs = data_graph->getVertexNeighbors(v, deg);
                for (ui i = 0; i < deg; ++i) {
                    candidate_neighbor_marks[u][nbrs[i]] = 1;
                }
            }
        }
    }

    ui selectInitialRoot() const
    {
        ui best_u = 0;
        for (ui u = 1; u < qn; ++u) {
            if (candidates[u].size() < candidates[best_u].size() ||
                (candidates[u].size() == candidates[best_u].size() &&
                 q_degree[u] > q_degree[best_u])) {
                best_u = u;
            }
        }
        return best_u;
    }

    void buildQueryPlan()
    {
        query_plan.clear();
        if (qn == 0) return;

        vector<char> used(qn, 0);
        ui root = selectInitialRoot();
        query_plan.push_back(root);
        used[root] = 1;

        while (query_plan.size() < qn) {
            ui best_u = qn;
            for (ui planned_u : query_plan) {
                for (ui nbr : q_neighbors[planned_u]) {
                    if (used[nbr]) continue;
                    if (best_u == qn ||
                        candidates[nbr].size() < candidates[best_u].size() ||
                        (candidates[nbr].size() == candidates[best_u].size() &&
                         (q_degree[nbr] > q_degree[best_u] ||
                          (q_degree[nbr] == q_degree[best_u] && nbr < best_u)))) {
                        best_u = nbr;
                    }
                }
            }

            if (best_u == qn) {
                for (ui u = 0; u < qn; ++u) {
                    if (used[u]) continue;
                    if (best_u == qn ||
                        candidates[u].size() < candidates[best_u].size() ||
                        (candidates[u].size() == candidates[best_u].size() &&
                         (q_degree[u] > q_degree[best_u] ||
                          (q_degree[u] == q_degree[best_u] && u < best_u)))) {
                        best_u = u;
                    }
                }
            }

            assert(best_u != qn);
            used[best_u] = 1;
            query_plan.push_back(best_u);
        }
    }

    ui queryPlanRank(ui u) const
    {
        for (ui i = 0; i < query_plan.size(); ++i) {
            if (query_plan[i] == u) return i;
        }
        return qn + u;
    }

    ui countMappedQueryNeighbors(ui u) const
    {
        ui count = 0;
        for (ui nbr : q_neighbors[u]) {
            if (mapped_q[nbr] != -1) ++count;
        }
        return count;
    }

    ui selectNextQueryVertex(ui)
    {
        ui best_u = qn;
        ui best_live = 0;

        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] != -1) continue;

            ui live = countMappedQueryNeighbors(u);
            if (best_u == qn) {
                best_u = u;
                best_live = live;
                continue;
            }

            bool best_connected = best_live > 0;
            bool curr_connected = live > 0;
            if (curr_connected != best_connected) {
                if (curr_connected) {
                    best_u = u;
                    best_live = live;
                }
                continue;
            }

            if (curr_connected && live != best_live) {
                if (live > best_live) {
                    best_u = u;
                    best_live = live;
                }
                continue;
            }

            if (candidates[u].size() < candidates[best_u].size() ||
                (candidates[u].size() == candidates[best_u].size() &&
                 (q_degree[u] > q_degree[best_u] ||
                  (q_degree[u] == q_degree[best_u] &&
                   queryPlanRank(u) < queryPlanRank(best_u))))) {
                best_u = u;
                best_live = live;
            }
        }

        return best_u;
    }

    ui missingEdgesToMappedNeighbors(ui u, ui v, ui limit) const
    {
        ui missing = 0;
        for (ui nbr : q_neighbors[u]) {
            int mapped_v = mapped_q[nbr];
            if (mapped_v == -1) continue;

            if (!data_graph->hasEdge(v, (ui)mapped_v)) {
                ++missing;
                if (missing > limit) return missing;
            }
        }
        return missing;
    }

    bool existsUnusedCandidateAdjacent(ui u, ui data_v) const
    {
        if (u >= qn) return false;
        if (!candidate_neighbor_marks.empty() &&
            !candidate_neighbor_marks[u][data_v]) {
            return false;
        }

        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors(data_v, deg);
        if (deg <= candidates[u].size()) {
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (mapped_g[v] == -1 && candidate_marks[u][v]) return true;
            }
            return false;
        }

        for (ui v : candidates[u]) {
            if (mapped_g[v] != -1) continue;
            if (data_graph->hasEdge(data_v, v)) return true;
        }
        return false;
    }

    ui minMissingToMappedNeighbors(ui u)
    {
        ui live = 0;
        for (ui nbr : q_neighbors[u]) {
            if (mapped_q[nbr] != -1) {
                ++live;
            }
        }
        if (live == 0) return 0;

        advanceBranchStamp();
        ui best_hits = 0;

        for (ui nbr : q_neighbors[u]) {
            if (mapped_q[nbr] == -1) continue;

            ui anchor_v = (ui)mapped_q[nbr];
            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors(anchor_v, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (mapped_g[v] != -1 || !candidate_marks[u][v]) continue;

                if (branch_seen_stamp[v] != branch_stamp) {
                    branch_seen_stamp[v] = branch_stamp;
                    branch_hit_count[v] = 0;
                }
                ++branch_hit_count[v];
                if ((ui)branch_hit_count[v] > best_hits) {
                    best_hits = (ui)branch_hit_count[v];
                    if (best_hits == live) return 0;
                }
            }
        }

        return live - best_hits;
    }

    ui boundaryMissingLowerBound(ui cost, ui limit)
    {
        ui lower_bound = 0;
        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] != -1) continue;
            if (countMappedQueryNeighbors(u) == 0) continue;

            ui remaining = limit >= lower_bound ? limit - lower_bound : 0;
            ui min_missing = minMissingToMappedNeighbors(u);
            if (min_missing > remaining) return limit + 1;
            lower_bound += min_missing;
            if (lower_bound > limit || cost + lower_bound > threshold) {
                return lower_bound;
            }
        }
        return lower_bound;
    }

    bool candidateCanStillConnect(ui u, ui v) const
    {
        if (qn <= 1) return true;

        for (ui nbr : q_neighbors[u]) {
            int mapped_v = mapped_q[nbr];
            if (mapped_v != -1) {
                if (data_graph->hasEdge(v, (ui)mapped_v)) return true;
            }
            else if (existsUnusedCandidateAdjacent(nbr, v)) {
                return true;
            }
        }

        return false;
    }

    bool retainedComponentsCanStillConnect() const
    {
        if (qn <= 1) return true;

        ui mapped_count = 0;
        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] != -1) ++mapped_count;
        }
        if (mapped_count == 0) return true;

        vector<int> component(qn, -1);
        ui component_count = 0;
        for (ui start = 0; start < qn; ++start) {
            if (mapped_q[start] == -1 || component[start] != -1) continue;

            queue<ui> q;
            q.push(start);
            component[start] = (int)component_count;

            while (!q.empty()) {
                ui u = q.front();
                q.pop();
                for (ui nbr : q_neighbors[u]) {
                    if (mapped_q[nbr] == -1 || component[nbr] != -1) continue;
                    if (!data_graph->hasEdge((ui)mapped_q[u], (ui)mapped_q[nbr])) continue;
                    component[nbr] = (int)component_count;
                    q.push(nbr);
                }
            }

            ++component_count;
        }

        if (mapped_count == qn) return component_count == 1;

        vector<char> has_future_bridge(component_count, 0);
        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] == -1) continue;
            int comp = component[u];
            if (comp < 0) continue;

            for (ui nbr : q_neighbors[u]) {
                if (mapped_q[nbr] != -1) continue;
                if (existsUnusedCandidateAdjacent(nbr, (ui)mapped_q[u])) {
                    has_future_bridge[(ui)comp] = 1;
                    break;
                }
            }
        }

        for (ui c = 0; c < component_count; ++c) {
            if (!has_future_bridge[c]) return false;
        }
        return true;
    }

    bool partialFeasibilityCheck(ui cost)
    {
        if (cost > threshold) return false;

        ui remaining_budget = threshold - cost;
        ui boundary_lb = boundaryMissingLowerBound(cost, remaining_budget);
        if (boundary_lb > remaining_budget || cost + boundary_lb > threshold) {
            stats.forced_lower_bound_prunes++;
            return false;
        }

        if (!retainedComponentsCanStillConnect()) {
            stats.connectivity_prunes++;
            return false;
        }

        return true;
    }

    void advanceBranchStamp()
    {
        ++branch_stamp;
        if (branch_stamp == numeric_limits<int>::max()) {
            fill(branch_seen_stamp.begin(), branch_seen_stamp.end(), 0);
            branch_stamp = 1;
        }
    }

    const vector<ui> *collectBranchCandidates(ui u, ui remaining_budget, vector<ui> &out)
    {
        out.clear();

        ui live = 0;
        for (ui nbr : q_neighbors[u]) {
            if (mapped_q[nbr] != -1) {
                ++live;
            }
        }

        if (live == 0 || live <= remaining_budget) {
            return &candidates[u];
        }

        ui required_hits = live - remaining_budget;
        branch_touched.clear();
        advanceBranchStamp();

        for (ui nbr : q_neighbors[u]) {
            if (mapped_q[nbr] == -1) continue;

            ui anchor_v = (ui)mapped_q[nbr];
            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors(anchor_v, deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (mapped_g[v] != -1 || !candidate_marks[u][v]) continue;

                if (branch_seen_stamp[v] != branch_stamp) {
                    branch_seen_stamp[v] = branch_stamp;
                    branch_hit_count[v] = 0;
                    branch_touched.push_back(v);
                }
                ++branch_hit_count[v];
            }
        }

        out.reserve(branch_touched.size());
        for (ui v : branch_touched) {
            if ((ui)branch_hit_count[v] >= required_hits) {
                out.push_back(v);
            }
        }
        return &out;
    }

    bool mappedImageIsConnected() const
    {
        vector<char> visited(qn, 0);
        queue<ui> q;
        q.push(0);
        ui visited_count = 0;

        while (!q.empty()) {
            ui u = q.front();
            q.pop();
            if (visited[u]) continue;

            visited[u] = 1;
            ++visited_count;

            for (ui nbr : q_neighbors[u]) {
                if (visited[nbr]) continue;
                assert(mapped_q[u] != -1 && mapped_q[nbr] != -1);
                if (data_graph->hasEdge((ui)mapped_q[u], (ui)mapped_q[nbr])) {
                    q.push(nbr);
                }
            }
        }

        return visited_count == qn;
    }

    void emitResult()
    {
        stats.result_count++;
#ifndef NDEBUG
        vector<pair<ui, ui> > result;
        result.reserve(qn);
        for (ui u = 0; u < qn; ++u) {
            result.push_back({ u, (ui)mapped_q[u] });
        }
        results_ptr->push_back(result);
#endif
        noteOutputLimitIfReached();
    }

    void dfs(ui depth, ui cost)
    {
        if (outputLimitReached()) return;
        if (cost > threshold) {
            stats.prune_calls++;
            return;
        }

        stats.recursion_calls++;

        if (depth == qn) {
            if (mappedImageIsConnected()) {
                emitResult();
            }
            else {
                stats.prune_calls++;
            }
            return;
        }

        ui remaining_budget = threshold - cost;
        ui u = selectNextQueryVertex(remaining_budget);
        if (u == qn) {
            stats.prune_calls++;
            return;
        }

        bool took_branch = false;
        vector<ui> branch_candidates;
        const vector<ui> *branch_source =
            collectBranchCandidates(u, remaining_budget, branch_candidates);

        for (ui v : *branch_source) {
            if (mapped_g[v] != -1) continue;
            if (!candidateCanStillConnect(u, v)) {
                stats.connectivity_prunes++;
                continue;
            }

            ui delta = missingEdgesToMappedNeighbors(u, v, remaining_budget);
            if (delta > remaining_budget) continue;

            took_branch = true;
            mapped_q[u] = (int)v;
            mapped_g[v] = (int)u;

            if (partialFeasibilityCheck(cost + delta)) {
                dfs(depth + 1, cost + delta);
            }
            else {
                stats.prune_calls++;
            }

            mapped_g[v] = -1;
            mapped_q[u] = -1;

            if (outputLimitReached()) break;
        }

        if (!took_branch) {
            stats.prune_calls++;
        }
    }
};

void Approximate_S3AND(const Graph *query_graph, const Graph *data_graph,
    vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    S3ANDSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
