#include "matching/algorithms/decq.h"

#include <deque>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>

namespace ssm {

struct DecQSolver::Impl {
    typedef std::vector<int> Mapping;
    typedef std::vector<Mapping> MatchSet;

    struct QueryEdge {
        ui u = 0;
        ui v = 0;
        long double selectivity = 0.0L;
    };

    struct Relation {
        std::vector<unsigned char> vertices;
        std::shared_ptr<MatchSet> rows;
    };

    struct Bundle {
        std::vector<Relation> components;
    };

    struct DecompositionNode {
        std::vector<size_t> edge_ids;
        int left = -1;
        int right = -1;
        std::unordered_map<std::string, std::shared_ptr<Bundle> > table;

        bool isLeaf() const { return left < 0 && right < 0; }
    };

    struct GlobalPattern {
        std::vector<unsigned char> missing;
        std::shared_ptr<MatchSet> matches;
        bool owns_live_rows = false;
        bool minimal = false;
    };

    DecQStats &stats;
    const Graph *query_graph = nullptr;
    const Graph *data_graph = nullptr;
    MatchResults *results = nullptr;
    ui threshold = 0;
    ui qn = 0;
    ui gn = 0;

    std::vector<QueryEdge> query_edges;
    std::vector<std::vector<size_t> > query_incident_edges;
    std::vector<std::vector<ui> > data_vertices_by_label;
    std::vector<std::unordered_map<LabelID, ui> > data_nlf;
    std::vector<DecompositionNode> decomposition;
    std::vector<std::vector<GlobalPattern> > lattice;
    std::vector<std::unordered_map<std::string, size_t> > lattice_indices;
    std::unordered_set<std::string> emitted_mappings;
    size_t live_match_rows = 0;

    explicit Impl(DecQStats &stats_ref) : stats(stats_ref) {}

    void reset()
    {
        stats = DecQStats();
        query_graph = nullptr;
        data_graph = nullptr;
        results = nullptr;
        threshold = 0;
        qn = 0;
        gn = 0;
        query_edges.clear();
        query_incident_edges.clear();
        data_vertices_by_label.clear();
        data_nlf.clear();
        decomposition.clear();
        lattice.clear();
        lattice_indices.clear();
        emitted_mappings.clear();
        live_match_rows = 0;
    }

    static bool sameUndirectedEdgeEndpoint(const QueryEdge &a, const QueryEdge &b)
    {
        return a.u == b.u || a.u == b.v || a.v == b.u || a.v == b.v;
    }

    static std::string binaryKey(const std::vector<unsigned char> &bits)
    {
        if (bits.empty()) return std::string();
        return std::string(reinterpret_cast<const char *>(bits.data()), bits.size());
    }

    std::string nodeKey(size_t node_id,
        const std::vector<unsigned char> &global_missing) const
    {
        const std::vector<size_t> &edge_ids = decomposition[node_id].edge_ids;
        std::string key;
        key.reserve(edge_ids.size());
        for (size_t edge_id : edge_ids) {
            key.push_back(static_cast<char>(global_missing[edge_id]));
        }
        return key;
    }

    void buildQueryEdges()
    {
        query_edges.clear();
        query_incident_edges.assign(qn, std::vector<size_t>());
        for (ui u = 0; u < qn; ++u) {
            ui degree = 0;
            const ui *neighbors = query_graph->getVertexNeighbors(u, degree);
            for (ui i = 0; i < degree; ++i) {
                ui v = neighbors[i];
                if (u >= v) continue;
                QueryEdge edge;
                edge.u = u;
                edge.v = v;
                size_t edge_id = query_edges.size();
                query_edges.push_back(edge);
                query_incident_edges[u].push_back(edge_id);
                query_incident_edges[v].push_back(edge_id);
            }
        }
        stats.query_vertices = qn;
        stats.query_edges = query_edges.size();
    }

    void buildDataIndex()
    {
        LabelID maximum_label = 0;
        for (ui u = 0; u < qn; ++u) {
            maximum_label = std::max(maximum_label, query_graph->getVertexLabel(u));
        }
        for (ui v = 0; v < gn; ++v) {
            maximum_label = std::max(maximum_label, data_graph->getVertexLabel(v));
        }

        data_vertices_by_label.assign(static_cast<size_t>(maximum_label) + 1,
            std::vector<ui>());
        for (ui v = 0; v < gn; ++v) {
            data_vertices_by_label[data_graph->getVertexLabel(v)].push_back(v);
        }
        for (const std::vector<ui> &vertices : data_vertices_by_label) {
            if (!vertices.empty()) stats.indexed_labels++;
        }

        data_nlf.assign(gn, std::unordered_map<LabelID, ui>());
        for (ui v = 0; v < gn; ++v) {
            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(v, degree);
            std::unordered_map<LabelID, ui> &counts = data_nlf[v];
            for (ui i = 0; i < degree; ++i) {
                counts[data_graph->getVertexLabel(neighbors[i])]++;
            }
        }

        long double selectivity_sum = 0.0L;
        long double selectivity_min = std::numeric_limits<long double>::max();
        long double selectivity_max = 0.0L;
        std::map<std::pair<LabelID, LabelID>, long double> selectivity_cache;
        for (QueryEdge &edge : query_edges) {
            LabelID left_label = query_graph->getVertexLabel(edge.u);
            LabelID right_label = query_graph->getVertexLabel(edge.v);
            std::pair<LabelID, LabelID> label_pair = std::minmax(left_label,
                right_label);
            std::map<std::pair<LabelID, LabelID>, long double>::const_iterator cached =
                selectivity_cache.find(label_pair);
            if (cached != selectivity_cache.end()) {
                edge.selectivity = cached->second;
                selectivity_sum += edge.selectivity;
                selectivity_min = std::min(selectivity_min, edge.selectivity);
                selectivity_max = std::max(selectivity_max, edge.selectivity);
                continue;
            }

            left_label = label_pair.first;
            right_label = label_pair.second;
            const std::vector<ui> &left_vertices = data_vertices_by_label[left_label];
            const std::vector<ui> &right_vertices = data_vertices_by_label[right_label];

            size_t edge_mapping_count = 0;
            for (ui data_u : left_vertices) {
                ui degree = 0;
                const ui *neighbors = data_graph->getVertexNeighbors(data_u, degree);
                for (ui i = 0; i < degree; ++i) {
                    if (data_graph->getVertexLabel(neighbors[i]) == right_label) {
                        edge_mapping_count++;
                    }
                }
            }

            long double denominator = 0.0L;
            if (left_label == right_label) {
                size_t count = left_vertices.size();
                if (count > 1) {
                    denominator = static_cast<long double>(count) * (count - 1);
                }
            }
            else {
                denominator = static_cast<long double>(left_vertices.size()) *
                    static_cast<long double>(right_vertices.size());
            }
            edge.selectivity = denominator > 0.0L
                ? static_cast<long double>(edge_mapping_count) / denominator
                : 0.0L;
            selectivity_cache.insert(std::make_pair(label_pair, edge.selectivity));
            selectivity_sum += edge.selectivity;
            selectivity_min = std::min(selectivity_min, edge.selectivity);
            selectivity_max = std::max(selectivity_max, edge.selectivity);
        }

        if (!query_edges.empty()) {
            stats.minimum_edge_selectivity = static_cast<double>(selectivity_min);
            stats.average_edge_selectivity = static_cast<double>(
                selectivity_sum / query_edges.size());
            stats.maximum_edge_selectivity = static_cast<double>(selectivity_max);
        }
    }

    bool queryConnectedAfterMissing(
        const std::vector<unsigned char> &missing) const
    {
        if (qn == 0) return false;
        if (qn == 1) return true;
        std::vector<unsigned char> visited(qn, 0);
        std::queue<ui> pending;
        visited[0] = 1;
        pending.push(0);
        size_t reached = 0;
        while (!pending.empty()) {
            ui u = pending.front();
            pending.pop();
            reached++;
            for (size_t edge_id : query_incident_edges[u]) {
                if (missing[edge_id]) continue;
                const QueryEdge &edge = query_edges[edge_id];
                ui v = edge.u == u ? edge.v : edge.u;
                if (!visited[v]) {
                    visited[v] = 1;
                    pending.push(v);
                }
            }
        }
        return reached == qn;
    }

    std::vector<unsigned char> verticesForEdges(
        const std::vector<size_t> &edge_ids) const
    {
        std::vector<unsigned char> vertices(qn, 0);
        for (size_t edge_id : edge_ids) {
            vertices[query_edges[edge_id].u] = 1;
            vertices[query_edges[edge_id].v] = 1;
        }
        if (qn == 1 && edge_ids.empty()) vertices[0] = 1;
        return vertices;
    }

    long double labelCandidateCount(ui query_vertex) const
    {
        LabelID label = query_graph->getVertexLabel(query_vertex);
        if (static_cast<size_t>(label) >= data_vertices_by_label.size()) return 0.0L;
        return static_cast<long double>(data_vertices_by_label[label].size());
    }

    long double estimateMatches(const std::vector<unsigned char> &vertices,
        const std::vector<unsigned char> &active_edges) const
    {
        long double estimate = 1.0L;
        bool has_vertex = false;
        for (ui u = 0; u < qn; ++u) {
            if (!vertices[u]) continue;
            has_vertex = true;
            estimate *= labelCandidateCount(u);
        }
        if (!has_vertex) return 1.0L;
        for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
            if (active_edges[edge_id] && vertices[query_edges[edge_id].u] &&
                vertices[query_edges[edge_id].v]) {
                estimate *= query_edges[edge_id].selectivity;
            }
        }
        return estimate;
    }

    std::vector<ui> genOrder(const std::vector<unsigned char> &vertices,
        const std::vector<unsigned char> &active_edges) const
    {
        std::vector<ui> order;
        std::vector<unsigned char> selected(qn, 0);
        size_t vertex_count = 0;
        for (unsigned char active : vertices) vertex_count += active != 0;
        order.reserve(vertex_count);

        long double prefix_estimate = 1.0L;
        while (order.size() < vertex_count) {
            ui best = qn;
            long double best_estimate = std::numeric_limits<long double>::max();
            for (ui u = 0; u < qn; ++u) {
                if (!vertices[u] || selected[u]) continue;
                long double estimate = prefix_estimate * labelCandidateCount(u);
                for (size_t edge_id : query_incident_edges[u]) {
                    if (!active_edges[edge_id]) continue;
                    const QueryEdge &edge = query_edges[edge_id];
                    ui other = edge.u == u ? edge.v : edge.u;
                    if (selected[other]) estimate *= edge.selectivity;
                }
                if (best == qn || estimate < best_estimate ||
                    (estimate == best_estimate && u < best)) {
                    best = u;
                    best_estimate = estimate;
                }
            }
            if (best == qn) break;
            selected[best] = 1;
            order.push_back(best);
            prefix_estimate = best_estimate;
        }
        return order;
    }

    void masksForEdgeSet(const std::vector<size_t> &edge_ids,
        std::vector<unsigned char> &vertices,
        std::vector<unsigned char> &active_edges) const
    {
        vertices = verticesForEdges(edge_ids);
        active_edges.assign(query_edges.size(), 0);
        for (size_t edge_id : edge_ids) active_edges[edge_id] = 1;
    }

    long double estimateFullMatches(const std::vector<size_t> &edge_ids) const
    {
        std::vector<unsigned char> vertices;
        std::vector<unsigned char> active_edges;
        masksForEdgeSet(edge_ids, vertices, active_edges);
        return estimateMatches(vertices, active_edges);
    }

    long double estimateIntermediateMatches(
        const std::vector<size_t> &edge_ids) const
    {
        std::vector<unsigned char> vertices;
        std::vector<unsigned char> active_edges;
        masksForEdgeSet(edge_ids, vertices, active_edges);
        std::vector<ui> order = genOrder(vertices, active_edges);
        std::vector<unsigned char> prefix(qn, 0);
        long double total = 0.0L;
        for (size_t i = 0; i + 1 < order.size(); ++i) {
            prefix[order[i]] = 1;
            total += estimateMatches(prefix, active_edges);
        }
        return total;
    }

    long double localPatternMultiplier(size_t edge_count) const
    {
        size_t maximum_deletions = std::min<size_t>(threshold, edge_count);
        long double total = 1.0L;
        long double combination = 1.0L;
        for (size_t i = 1; i <= maximum_deletions; ++i) {
            combination *= static_cast<long double>(edge_count - i + 1) /
                static_cast<long double>(i);
            total += combination;
        }
        return total;
    }

    long double splitGain(const std::vector<size_t> &parent,
        const std::vector<size_t> &left,
        const std::vector<size_t> &right) const
    {
        long double parent_search = localPatternMultiplier(parent.size()) *
            estimateIntermediateMatches(parent);
        long double left_search = localPatternMultiplier(left.size()) *
            (estimateFullMatches(left) + estimateIntermediateMatches(left));
        long double right_search = localPatternMultiplier(right.size()) *
            (estimateFullMatches(right) + estimateIntermediateMatches(right));
        return parent_search - left_search - right_search;
    }

    void collectTreeSubtree(size_t root,
        const std::vector<std::vector<size_t> > &children,
        std::vector<unsigned char> &in_subtree) const
    {
        in_subtree[root] = 1;
        for (size_t child : children[root]) {
            collectTreeSubtree(child, children, in_subtree);
        }
    }

    void splitNode(size_t node_id)
    {
        const std::vector<size_t> parent_edges = decomposition[node_id].edge_ids;
        const size_t minimum_fragment_size = static_cast<size_t>(threshold) + 1;
        if (parent_edges.size() < 2 * minimum_fragment_size) return;

        const size_t local_edge_count = parent_edges.size();
        std::vector<std::vector<size_t> > line_adjacency(local_edge_count);
        for (size_t i = 0; i < local_edge_count; ++i) {
            for (size_t j = i + 1; j < local_edge_count; ++j) {
                if (sameUndirectedEdgeEndpoint(query_edges[parent_edges[i]],
                    query_edges[parent_edges[j]])) {
                    line_adjacency[i].push_back(j);
                    line_adjacency[j].push_back(i);
                }
            }
        }

        struct CandidateSplit {
            std::vector<size_t> left;
            std::vector<size_t> right;
            size_t imbalance = 0;
            long double gain = 0.0L;
            std::string key;
        };
        std::vector<CandidateSplit> candidates;
        std::unordered_set<std::string> seen;

        for (size_t bfs_root = 0; bfs_root < local_edge_count; ++bfs_root) {
            std::vector<size_t> parent(local_edge_count, local_edge_count);
            std::vector<std::vector<size_t> > children(local_edge_count);
            std::queue<size_t> pending;
            parent[bfs_root] = bfs_root;
            pending.push(bfs_root);
            while (!pending.empty()) {
                size_t u = pending.front();
                pending.pop();
                for (size_t v : line_adjacency[u]) {
                    if (parent[v] != local_edge_count) continue;
                    parent[v] = u;
                    children[u].push_back(v);
                    pending.push(v);
                }
            }

            for (size_t cut_root = 0; cut_root < local_edge_count; ++cut_root) {
                if (cut_root == bfs_root || parent[cut_root] == local_edge_count) continue;
                std::vector<unsigned char> in_left(local_edge_count, 0);
                collectTreeSubtree(cut_root, children, in_left);
                size_t left_size = std::count(in_left.begin(), in_left.end(), 1);
                size_t right_size = local_edge_count - left_size;
                if (left_size < minimum_fragment_size ||
                    right_size < minimum_fragment_size) continue;

                std::string left_key;
                std::string right_key;
                left_key.reserve(local_edge_count);
                right_key.reserve(local_edge_count);
                for (unsigned char value : in_left) {
                    left_key.push_back(static_cast<char>(value));
                    right_key.push_back(static_cast<char>(!value));
                }
                bool swap_sides = right_key < left_key;
                const std::string canonical = swap_sides ? right_key : left_key;
                if (!seen.insert(canonical).second) continue;

                CandidateSplit candidate;
                candidate.key = canonical;
                for (size_t i = 0; i < local_edge_count; ++i) {
                    bool goes_left = static_cast<bool>(in_left[i]) != swap_sides;
                    (goes_left ? candidate.left : candidate.right).push_back(parent_edges[i]);
                }
                candidate.imbalance = candidate.left.size() > candidate.right.size()
                    ? candidate.left.size() - candidate.right.size()
                    : candidate.right.size() - candidate.left.size();
                candidates.push_back(candidate);
            }
        }

        stats.split_candidates += candidates.size();
        if (candidates.empty()) return;
        size_t best_balance = std::numeric_limits<size_t>::max();
        for (const CandidateSplit &candidate : candidates) {
            best_balance = std::min(best_balance, candidate.imbalance);
        }

        CandidateSplit *best = nullptr;
        for (CandidateSplit &candidate : candidates) {
            if (candidate.imbalance != best_balance) continue;
            candidate.gain = splitGain(parent_edges, candidate.left, candidate.right);
            if (candidate.gain <= 0.0L) continue;
            if (best == nullptr || candidate.gain > best->gain ||
                (candidate.gain == best->gain && candidate.key < best->key)) {
                best = &candidate;
            }
        }
        if (best == nullptr) return;

        size_t left_id = decomposition.size();
        decomposition.push_back(DecompositionNode());
        decomposition.back().edge_ids = best->left;
        size_t right_id = decomposition.size();
        decomposition.push_back(DecompositionNode());
        decomposition.back().edge_ids = best->right;
        decomposition[node_id].left = static_cast<int>(left_id);
        decomposition[node_id].right = static_cast<int>(right_id);
        stats.accepted_splits++;

        splitNode(left_id);
        splitNode(right_id);
    }

    long double estimatedPlanCost(size_t node_id) const
    {
        const DecompositionNode &node = decomposition[node_id];
        if (node.isLeaf()) {
            return localPatternMultiplier(node.edge_ids.size()) *
                (estimateFullMatches(node.edge_ids) +
                    estimateIntermediateMatches(node.edge_ids));
        }
        return estimatedPlanCost(static_cast<size_t>(node.left)) +
            estimatedPlanCost(static_cast<size_t>(node.right)) +
            localPatternMultiplier(node.edge_ids.size()) *
                estimateFullMatches(node.edge_ids);
    }

    void buildDecomposition()
    {
        decomposition.clear();
        decomposition.push_back(DecompositionNode());
        decomposition[0].edge_ids.resize(query_edges.size());
        std::iota(decomposition[0].edge_ids.begin(),
            decomposition[0].edge_ids.end(), static_cast<size_t>(0));
        stats.unsplit_estimated_cost = static_cast<double>(
            localPatternMultiplier(query_edges.size()) *
            (estimateFullMatches(decomposition[0].edge_ids) +
                estimateIntermediateMatches(decomposition[0].edge_ids)));
        splitNode(0);
        stats.decomposition_nodes = decomposition.size();
        for (const DecompositionNode &node : decomposition) {
            if (node.isLeaf()) stats.decomposition_leaves++;
        }
        stats.decomposition_estimated_cost = static_cast<double>(estimatedPlanCost(0));
    }

    void chooseRecursively(const std::vector<size_t> &items, size_t choose,
        size_t position, std::vector<size_t> &chosen,
        const std::function<void(const std::vector<size_t> &)> &callback) const
    {
        if (chosen.size() == choose) {
            callback(chosen);
            return;
        }
        size_t still_needed = choose - chosen.size();
        if (items.size() - position < still_needed) return;
        for (size_t i = position; i + still_needed <= items.size(); ++i) {
            chosen.push_back(items[i]);
            chooseRecursively(items, choose, i + 1, chosen, callback);
            chosen.pop_back();
        }
    }

    void forEachCombination(const std::vector<size_t> &items, size_t choose,
        const std::function<void(const std::vector<size_t> &)> &callback) const
    {
        std::vector<size_t> chosen;
        chosen.reserve(choose);
        chooseRecursively(items, choose, 0, chosen, callback);
    }

    void subsetLexicographicRec(const std::vector<size_t> &items,
        size_t maximum_size, size_t position, std::vector<size_t> &chosen,
        const std::function<void(const std::vector<size_t> &)> &callback) const
    {
        if (chosen.size() == maximum_size) return;
        for (size_t i = position; i < items.size(); ++i) {
            chosen.push_back(items[i]);
            callback(chosen);
            subsetLexicographicRec(items, maximum_size, i + 1, chosen, callback);
            chosen.pop_back();
        }
    }

    void forEachSubsetLexicographic(const std::vector<size_t> &items,
        size_t maximum_size,
        const std::function<void(const std::vector<size_t> &)> &callback) const
    {
        std::vector<size_t> chosen;
        callback(chosen);
        subsetLexicographicRec(items, maximum_size, 0, chosen, callback);
    }

    void noteAllocatedRows(size_t rows)
    {
        stats.materialized_match_rows += rows;
        live_match_rows += rows;
        stats.peak_live_match_rows = std::max(stats.peak_live_match_rows,
            live_match_rows);
#if DECQ_INTERMEDIATE_MATCH_LIMIT > 0
        if (live_match_rows > static_cast<size_t>(DECQ_INTERMEDIATE_MATCH_LIMIT)) {
            stats.intermediate_limit_reached = true;
        }
#endif
    }

    void releaseOwnedRows(GlobalPattern &pattern)
    {
        if (pattern.owns_live_rows && pattern.matches) {
            assert(live_match_rows >= pattern.matches->size());
            live_match_rows -= pattern.matches->size();
        }
        pattern.matches.reset();
        pattern.owns_live_rows = false;
    }

    void releaseTemporaryRelationRows(const Relation &relation)
    {
        if (relation.rows && relation.rows.use_count() == 1) {
            assert(live_match_rows >= relation.rows->size());
            live_match_rows -= relation.rows->size();
        }
    }

    bool stopped() const
    {
        return stats.intermediate_limit_reached || stats.output_limit_reached;
    }

    bool pendingRowsExceedLimit(size_t pending_rows)
    {
#if DECQ_INTERMEDIATE_MATCH_LIMIT > 0
        if (live_match_rows + pending_rows >
            static_cast<size_t>(DECQ_INTERMEDIATE_MATCH_LIMIT)) {
            stats.intermediate_limit_reached = true;
            return true;
        }
#else
        (void)pending_rows;
#endif
        return false;
    }

    Relation exactMatch(const std::vector<unsigned char> &vertices,
        const std::vector<size_t> &component_edges)
    {
        Timer timer;
        Relation relation;
        relation.vertices = vertices;
        relation.rows.reset(new MatchSet());
        stats.exact_match_executions++;

        std::vector<unsigned char> active_edges(query_edges.size(), 0);
        std::vector<std::vector<ui> > active_neighbors(qn);
        for (size_t edge_id : component_edges) {
            active_edges[edge_id] = 1;
            const QueryEdge &edge = query_edges[edge_id];
            active_neighbors[edge.u].push_back(edge.v);
            active_neighbors[edge.v].push_back(edge.u);
        }

        std::vector<std::vector<ui> > candidates(qn);
        for (ui u = 0; u < qn; ++u) {
            if (!vertices[u]) continue;
            LabelID label = query_graph->getVertexLabel(u);
            if (static_cast<size_t>(label) >= data_vertices_by_label.size()) continue;

            std::unordered_map<LabelID, ui> required_nlf;
            for (ui neighbor : active_neighbors[u]) {
                required_nlf[query_graph->getVertexLabel(neighbor)]++;
            }
            for (ui v : data_vertices_by_label[label]) {
                if (data_graph->getVertexDegree(v) < active_neighbors[u].size()) continue;
                bool nlf_ok = true;
                for (const std::pair<const LabelID, ui> &requirement : required_nlf) {
                    std::unordered_map<LabelID, ui>::const_iterator it =
                        data_nlf[v].find(requirement.first);
                    if (it == data_nlf[v].end() || it->second < requirement.second) {
                        nlf_ok = false;
                        break;
                    }
                }
                if (nlf_ok) candidates[u].push_back(v);
            }
        }

        std::vector<ui> order = genOrder(vertices, active_edges);
        Mapping mapping(qn, -1);
        std::vector<unsigned char> used_data(gn, 0);
        std::function<void(size_t)> search = [&](size_t depth) {
            if (stats.intermediate_limit_reached) return;
            if (depth == order.size()) {
                if (pendingRowsExceedLimit(relation.rows->size() + 1)) return;
                relation.rows->push_back(mapping);
                return;
            }

            ui u = order[depth];
            const auto explore_candidate = [&](ui v) {
                if (used_data[v]) return;
                bool compatible = true;
                for (ui query_neighbor : active_neighbors[u]) {
                    int mapped_neighbor = mapping[query_neighbor];
                    if (mapped_neighbor >= 0 &&
                        !data_graph->hasEdge(v, static_cast<ui>(mapped_neighbor))) {
                        compatible = false;
                        break;
                    }
                }
                if (!compatible) return;

                mapping[u] = static_cast<int>(v);
                used_data[v] = 1;
                if (depth + 1 < order.size()) stats.local_search_partial_rows++;
                search(depth + 1);
                used_data[v] = 0;
                mapping[u] = -1;
            };

            ui anchor_data = gn;
            ui anchor_degree = std::numeric_limits<ui>::max();
            for (ui query_neighbor : active_neighbors[u]) {
                if (mapping[query_neighbor] < 0) continue;
                ui mapped_neighbor = static_cast<ui>(mapping[query_neighbor]);
                ui degree = data_graph->getVertexDegree(mapped_neighbor);
                if (anchor_data == gn || degree < anchor_degree) {
                    anchor_data = mapped_neighbor;
                    anchor_degree = degree;
                }
            }
            if (anchor_data != gn) {
                ui degree = 0;
                const ui *neighbors = data_graph->getVertexNeighbors(anchor_data, degree);
                for (ui i = 0; i < degree; ++i) {
                    ui v = neighbors[i];
                    if (!std::binary_search(candidates[u].begin(), candidates[u].end(), v)) {
                        continue;
                    }
                    explore_candidate(v);
                    if (stats.intermediate_limit_reached) break;
                }
            }
            else {
                for (ui v : candidates[u]) {
                    explore_candidate(v);
                    if (stats.intermediate_limit_reached) break;
                }
            }
        };
        search(0);

        noteAllocatedRows(relation.rows->size());
        stats.local_matching_time += timer.elapsed();
        stats.local_match_rows += relation.rows->size();
        return relation;
    }

    std::shared_ptr<Bundle> buildLocalBundle(size_t node_id,
        const std::vector<unsigned char> &global_missing)
    {
        const DecompositionNode &node = decomposition[node_id];
        std::vector<size_t> active_edge_ids;
        for (size_t edge_id : node.edge_ids) {
            if (!global_missing[edge_id]) active_edge_ids.push_back(edge_id);
        }

        std::vector<std::vector<ui> > adjacency(qn);
        std::vector<unsigned char> active_vertices(qn, 0);
        for (size_t edge_id : active_edge_ids) {
            const QueryEdge &edge = query_edges[edge_id];
            active_vertices[edge.u] = 1;
            active_vertices[edge.v] = 1;
            adjacency[edge.u].push_back(edge.v);
            adjacency[edge.v].push_back(edge.u);
        }
        if (qn == 1 && node.edge_ids.empty()) active_vertices[0] = 1;
        std::vector<unsigned char> fragment_vertices = verticesForEdges(node.edge_ids);
        bool removed_isolated_vertex = false;
        for (ui u = 0; u < qn; ++u) {
            if (fragment_vertices[u] && !active_vertices[u]) {
                removed_isolated_vertex = true;
                break;
            }
        }

        std::shared_ptr<Bundle> bundle(new Bundle());
        std::vector<unsigned char> visited(qn, 0);
        for (ui start = 0; start < qn; ++start) {
            if (!active_vertices[start] || visited[start]) continue;
            std::vector<unsigned char> component_vertices(qn, 0);
            std::queue<ui> pending;
            pending.push(start);
            visited[start] = 1;
            while (!pending.empty()) {
                ui u = pending.front();
                pending.pop();
                component_vertices[u] = 1;
                for (ui v : adjacency[u]) {
                    if (!visited[v]) {
                        visited[v] = 1;
                        pending.push(v);
                    }
                }
            }

            std::vector<size_t> component_edges;
            for (size_t edge_id : active_edge_ids) {
                const QueryEdge &edge = query_edges[edge_id];
                if (component_vertices[edge.u] && component_vertices[edge.v]) {
                    component_edges.push_back(edge_id);
                }
            }
            bundle->components.push_back(exactMatch(component_vertices,
                component_edges));
            stats.local_components++;
            if (stats.intermediate_limit_reached) break;
        }
        if (removed_isolated_vertex || bundle->components.size() > 1) {
            stats.disconnected_local_patterns++;
        }
        return bundle;
    }

    void enumerateLocalPatterns()
    {
        Timer phase_timer;
        for (size_t node_id = 0; node_id < decomposition.size(); ++node_id) {
            DecompositionNode &node = decomposition[node_id];
            if (!node.isLeaf()) continue;
            size_t maximum_deletions = std::min<size_t>(threshold,
                node.edge_ids.size());
            forEachSubsetLexicographic(node.edge_ids, maximum_deletions,
                [&](const std::vector<size_t> &deleted_edges) {
                    if (stats.intermediate_limit_reached) return;
                    std::vector<unsigned char> missing(query_edges.size(), 0);
                    for (size_t edge_id : deleted_edges) missing[edge_id] = 1;
                    if (!queryConnectedAfterMissing(missing)) {
                        stats.local_edge_cuts_skipped++;
                        return;
                    }
                    std::string key = nodeKey(node_id, missing);
                    std::shared_ptr<Bundle> bundle = buildLocalBundle(node_id,
                        missing);
                    decomposition[node_id].table.insert(std::make_pair(key, bundle));
                    stats.local_patterns++;
                });
            if (stats.intermediate_limit_reached) break;
        }
        long long elapsed = phase_timer.elapsed();
        stats.local_pattern_time = elapsed > stats.local_matching_time
            ? elapsed - stats.local_matching_time : 0;
    }

    static bool relationsOverlap(const Relation &left, const Relation &right)
    {
        for (size_t u = 0; u < left.vertices.size(); ++u) {
            if (left.vertices[u] && right.vertices[u]) return true;
        }
        return false;
    }

    static std::string joinKey(const Mapping &mapping,
        const std::vector<ui> &join_vertices)
    {
        std::string key;
        key.reserve(join_vertices.size() * sizeof(int));
        for (ui u : join_vertices) {
            int value = mapping[u];
            key.append(reinterpret_cast<const char *>(&value), sizeof(value));
        }
        return key;
    }

    Relation hashJoin(const Relation &left, const Relation &right)
    {
        stats.hash_joins++;
        stats.hash_join_input_rows += left.rows->size() + right.rows->size();
        Relation output;
        output.vertices.resize(qn, 0);
        for (ui u = 0; u < qn; ++u) {
            output.vertices[u] = left.vertices[u] || right.vertices[u];
        }
        output.rows.reset(new MatchSet());

        std::vector<ui> join_vertices;
        for (ui u = 0; u < qn; ++u) {
            if (left.vertices[u] && right.vertices[u]) join_vertices.push_back(u);
        }
        assert(!join_vertices.empty());

        const Relation *build = &left;
        const Relation *probe = &right;
        if (build->rows->size() > probe->rows->size()) std::swap(build, probe);

        std::unordered_map<std::string, std::vector<size_t> > hash_table;
        hash_table.reserve(build->rows->size());
        for (size_t i = 0; i < build->rows->size(); ++i) {
            hash_table[joinKey((*build->rows)[i], join_vertices)].push_back(i);
        }

        for (const Mapping &probe_mapping : *probe->rows) {
            std::unordered_map<std::string, std::vector<size_t> >::const_iterator it =
                hash_table.find(joinKey(probe_mapping, join_vertices));
            if (it == hash_table.end()) continue;
            for (size_t build_row_id : it->second) {
                const Mapping &build_mapping = (*build->rows)[build_row_id];
                Mapping merged = build_mapping;
                std::vector<int> used_data_vertices;
                used_data_vertices.reserve(qn);
                bool compatible = true;
                for (ui u = 0; u < qn; ++u) {
                    if (build_mapping[u] >= 0) {
                        int data_vertex = build_mapping[u];
                        if (std::find(used_data_vertices.begin(),
                            used_data_vertices.end(), data_vertex) !=
                            used_data_vertices.end()) {
                            compatible = false;
                            break;
                        }
                        used_data_vertices.push_back(data_vertex);
                    }
                }
                if (!compatible) continue;
                for (ui u = 0; u < qn; ++u) {
                    if (probe_mapping[u] < 0) continue;
                    if (merged[u] >= 0) {
                        if (merged[u] != probe_mapping[u]) {
                            compatible = false;
                            break;
                        }
                    }
                    else {
                        int data_vertex = probe_mapping[u];
                        if (std::find(used_data_vertices.begin(),
                            used_data_vertices.end(), data_vertex) !=
                            used_data_vertices.end()) {
                            compatible = false;
                            break;
                        }
                        used_data_vertices.push_back(data_vertex);
                        merged[u] = probe_mapping[u];
                    }
                }
                if (!compatible) continue;
                if (pendingRowsExceedLimit(output.rows->size() + 1)) break;
                output.rows->push_back(merged);
            }
            if (stats.intermediate_limit_reached) break;
        }

        stats.hash_join_output_rows += output.rows->size();
        noteAllocatedRows(output.rows->size());
        return output;
    }

    std::shared_ptr<Bundle> mergeBundles(const Bundle &left,
        const Bundle &right)
    {
        std::shared_ptr<Bundle> merged(new Bundle());
        merged->components = left.components;
        merged->components.insert(merged->components.end(),
            right.components.begin(), right.components.end());

        while (!stats.intermediate_limit_reached) {
            size_t best_i = merged->components.size();
            size_t best_j = merged->components.size();
            long double best_product = std::numeric_limits<long double>::max();
            for (size_t i = 0; i < merged->components.size(); ++i) {
                for (size_t j = i + 1; j < merged->components.size(); ++j) {
                    if (!relationsOverlap(merged->components[i],
                        merged->components[j])) continue;
                    long double product = static_cast<long double>(
                        merged->components[i].rows->size()) *
                        static_cast<long double>(merged->components[j].rows->size());
                    if (best_i == merged->components.size() ||
                        product < best_product) {
                        best_i = i;
                        best_j = j;
                        best_product = product;
                    }
                }
            }
            if (best_i == merged->components.size()) break;
            Relation joined = hashJoin(merged->components[best_i],
                merged->components[best_j]);
            releaseTemporaryRelationRows(merged->components[best_i]);
            releaseTemporaryRelationRows(merged->components[best_j]);
            merged->components.erase(merged->components.begin() + best_j);
            merged->components.erase(merged->components.begin() + best_i);
            merged->components.push_back(joined);
        }
        return merged;
    }

    std::shared_ptr<Bundle> materializeNode(size_t node_id,
        const std::vector<unsigned char> &missing)
    {
        DecompositionNode &node = decomposition[node_id];
        std::string key = nodeKey(node_id, missing);
        std::unordered_map<std::string, std::shared_ptr<Bundle> >::iterator found =
            node.table.find(key);
        if (found != node.table.end()) {
            stats.decomposition_cache_hits++;
            return found->second;
        }
        stats.decomposition_cache_misses++;
        if (node.isLeaf()) {
            return std::shared_ptr<Bundle>();
        }

        std::shared_ptr<Bundle> left = materializeNode(
            static_cast<size_t>(node.left), missing);
        std::shared_ptr<Bundle> right = materializeNode(
            static_cast<size_t>(node.right), missing);
        if (!left || !right || stats.intermediate_limit_reached) {
            return std::shared_ptr<Bundle>();
        }
        std::shared_ptr<Bundle> result = mergeBundles(*left, *right);
        node.table.insert(std::make_pair(key, result));
        return result;
    }

    void buildGlobalLattice()
    {
        lattice.assign(static_cast<size_t>(threshold) + 1,
            std::vector<GlobalPattern>());
        lattice_indices.assign(static_cast<size_t>(threshold) + 1,
            std::unordered_map<std::string, size_t>());
        std::vector<size_t> all_edges(query_edges.size());
        std::iota(all_edges.begin(), all_edges.end(), static_cast<size_t>(0));

        for (size_t deleted_count = 0; deleted_count <= threshold;
            ++deleted_count) {
            forEachCombination(all_edges, deleted_count,
                [&](const std::vector<size_t> &deleted_edges) {
                    stats.lattice_candidates++;
                    std::vector<unsigned char> missing(query_edges.size(), 0);
                    for (size_t edge_id : deleted_edges) missing[edge_id] = 1;
                    if (!queryConnectedAfterMissing(missing)) {
                        stats.global_edge_cuts_skipped++;
                        return;
                    }
                    GlobalPattern pattern;
                    pattern.missing.swap(missing);
                    size_t pattern_id = lattice[deleted_count].size();
                    lattice[deleted_count].push_back(pattern);
                    lattice_indices[deleted_count].insert(std::make_pair(
                        binaryKey(lattice[deleted_count].back().missing), pattern_id));
                    stats.global_patterns++;
                });
        }

        for (size_t level = 0; level < lattice.size(); ++level) {
            for (GlobalPattern &pattern : lattice[level]) {
                bool has_connected_child = false;
                if (level + 1 < lattice.size()) {
                    for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
                        if (pattern.missing[edge_id]) continue;
                        std::vector<unsigned char> child_missing = pattern.missing;
                        child_missing[edge_id] = 1;
                        if (lattice_indices[level + 1].find(binaryKey(child_missing)) !=
                            lattice_indices[level + 1].end()) {
                            has_connected_child = true;
                            break;
                        }
                    }
                }
                pattern.minimal = !has_connected_child;
                if (pattern.minimal) stats.minimal_patterns++;
                else stats.nonminimal_patterns++;
            }
        }
    }

    std::shared_ptr<MatchSet> validateRestoredEdge(const MatchSet &child_matches,
        size_t restored_edge)
    {
        std::shared_ptr<MatchSet> output(new MatchSet());
        const QueryEdge &edge = query_edges[restored_edge];
        for (const Mapping &mapping : child_matches) {
            stats.edge_validations++;
            if (mapping[edge.u] >= 0 && mapping[edge.v] >= 0 &&
                data_graph->hasEdge(static_cast<ui>(mapping[edge.u]),
                    static_cast<ui>(mapping[edge.v]))) {
                if (pendingRowsExceedLimit(output->size() + 1)) break;
                output->push_back(mapping);
            }
        }
        noteAllocatedRows(output->size());
        return output;
    }

    static std::string mappingKey(const Mapping &mapping)
    {
        if (mapping.empty()) return std::string();
        return std::string(reinterpret_cast<const char *>(mapping.data()),
            mapping.size() * sizeof(int));
    }

    bool emitCanonicalMapping(const Mapping &mapping)
    {
        std::string key = mappingKey(mapping);
        if (!emitted_mappings.insert(key).second) {
            stats.duplicate_canonical_mappings++;
            return true;
        }
#if MATCH_OUTPUT_LIMIT > 0
        if (stats.result_count >= static_cast<size_t>(MATCH_OUTPUT_LIMIT)) {
            stats.output_limit_reached = true;
            return false;
        }
#endif
        stats.result_count++;
#ifndef NDEBUG
        std::vector<std::pair<ui, ui> > result;
        result.reserve(qn);
        for (ui u = 0; u < qn; ++u) {
            result.push_back(std::make_pair(u, static_cast<ui>(mapping[u])));
        }
        results->push_back(result);
#endif
        return true;
    }

    void canonicalizePattern(const GlobalPattern &pattern)
    {
        if (!pattern.matches) return;
        for (const Mapping &mapping : *pattern.matches) {
            bool valid = mapping.size() == qn;
            std::vector<unsigned char> used_data(gn, 0);
            std::vector<unsigned char> actual_missing(query_edges.size(), 0);
            if (valid) {
                for (ui u = 0; u < qn; ++u) {
                    if (mapping[u] < 0 || static_cast<ui>(mapping[u]) >= gn ||
                        used_data[static_cast<ui>(mapping[u])] ||
                        query_graph->getVertexLabel(u) != data_graph->getVertexLabel(
                            static_cast<ui>(mapping[u]))) {
                        valid = false;
                        break;
                    }
                    used_data[static_cast<ui>(mapping[u])] = 1;
                }
            }
            if (!valid) {
                stats.invalid_global_rows++;
                continue;
            }
            for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
                const QueryEdge &edge = query_edges[edge_id];
                if (!data_graph->hasEdge(static_cast<ui>(mapping[edge.u]),
                    static_cast<ui>(mapping[edge.v]))) {
                    actual_missing[edge_id] = 1;
                }
            }
            if (actual_missing != pattern.missing) {
                stats.noncanonical_pattern_rows++;
                continue;
            }
            if (!emitCanonicalMapping(mapping)) return;
        }
    }

    bool materializeMinimalPattern(GlobalPattern &pattern)
    {
        std::shared_ptr<Bundle> bundle = materializeNode(0, pattern.missing);
        if (!bundle || stats.intermediate_limit_reached) return false;
        if (bundle->components.empty()) {
            pattern.matches.reset(new MatchSet());
            return true;
        }
        if (bundle->components.size() != 1) {
            // A connected global pattern makes the overlap graph of all local
            // components connected. More than one component signals an
            // implementation or malformed-input error, not an empty result.
            return false;
        }
        pattern.matches = bundle->components[0].rows;
        return true;
    }

    bool materializeNonminimalPattern(size_t level, GlobalPattern &pattern)
    {
        size_t best_child = lattice[level + 1].size();
        size_t best_size = std::numeric_limits<size_t>::max();
        size_t restored_edge = query_edges.size();
        for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
            if (pattern.missing[edge_id]) continue;
            std::vector<unsigned char> child_missing = pattern.missing;
            child_missing[edge_id] = 1;
            std::unordered_map<std::string, size_t>::const_iterator found =
                lattice_indices[level + 1].find(binaryKey(child_missing));
            if (found == lattice_indices[level + 1].end()) continue;
            const GlobalPattern &child = lattice[level + 1][found->second];
            if (!child.matches) continue;
            if (best_child == lattice[level + 1].size() ||
                child.matches->size() < best_size) {
                best_child = found->second;
                best_size = child.matches->size();
                restored_edge = edge_id;
            }
        }
        if (best_child == lattice[level + 1].size()) return false;
        pattern.matches = validateRestoredEdge(
            *lattice[level + 1][best_child].matches, restored_edge);
        pattern.owns_live_rows = true;
        return true;
    }

    void runGlobalMatching()
    {
        for (int signed_level = static_cast<int>(threshold);
            signed_level >= 0 && !stopped(); --signed_level) {
            size_t level = static_cast<size_t>(signed_level);
            for (GlobalPattern &pattern : lattice[level]) {
                bool materialized = false;
                if (pattern.minimal) {
                    Timer timer;
                    materialized = materializeMinimalPattern(pattern);
                    stats.minimal_merge_time += timer.elapsed();
                    stats.minimal_patterns_processed++;
                }
                else {
                    Timer timer;
                    materialized = materializeNonminimalPattern(level, pattern);
                    stats.edge_validation_time += timer.elapsed();
                    stats.nonminimal_patterns_processed++;
                }
                if (!materialized || stats.intermediate_limit_reached) {
                    stats.intermediate_limit_reached = true;
                    break;
                }
                stats.global_pattern_rows += pattern.matches->size();
                Timer canonical_timer;
                canonicalizePattern(pattern);
                stats.canonical_dedup_time += canonical_timer.elapsed();
                if (stopped()) break;
            }

            if (!stopped() && level + 1 < lattice.size()) {
                for (GlobalPattern &child : lattice[level + 1]) {
                    releaseOwnedRows(child);
                }
            }
        }
        for (std::vector<GlobalPattern> &level : lattice) {
            for (GlobalPattern &pattern : level) releaseOwnedRows(pattern);
        }
    }

    bool initialize(const Graph *query, const Graph *data, ui match_threshold)
    {
        reset();
        Timer timer;
        query_graph = query;
        data_graph = data;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();
        if (qn == 0 || gn == 0 || qn > gn) {
            stats.data_index_time = timer.elapsed();
            stats.total_time = stats.data_index_time;
            return false;
        }

        buildQueryEdges();
        if (qn > 1 && query_edges.size() < static_cast<size_t>(qn - 1)) {
            stats.data_index_time = timer.elapsed();
            stats.total_time = stats.data_index_time;
            return false;
        }
        std::vector<unsigned char> no_missing(query_edges.size(), 0);
        if (!queryConnectedAfterMissing(no_missing)) {
            stats.data_index_time = timer.elapsed();
            stats.total_time = stats.data_index_time;
            return false;
        }
        size_t minimum_connected_edges = qn > 0 ? static_cast<size_t>(qn - 1) : 0;
        size_t cycle_rank = query_edges.size() - minimum_connected_edges;
        threshold = static_cast<ui>(std::min<size_t>(threshold, cycle_rank));
        buildDataIndex();
        stats.data_index_time = timer.elapsed();
        return true;
    }

    void execute(MatchResults &output)
    {
        results = &output;
        results->clear();
        Timer total_timer;

        Timer phase_timer;
        buildDecomposition();
        stats.decomposition_time = phase_timer.elapsed();

        enumerateLocalPatterns();
        if (!stats.intermediate_limit_reached) {
            phase_timer.restart();
            buildGlobalLattice();
            stats.lattice_time = phase_timer.elapsed();
            runGlobalMatching();
        }
        stats.total_time = stats.data_index_time + total_timer.elapsed();
    }
};

DecQSolver::DecQSolver() : impl_(new Impl(stats)) {}

DecQSolver::~DecQSolver() = default;

bool DecQSolver::init(const Graph *query_graph, const Graph *data_graph,
    ui match_threshold)
{
    return impl_->initialize(query_graph, data_graph, match_threshold);
}

void DecQSolver::match(MatchResults &results)
{
    impl_->execute(results);
}

void DecQSolver::printStats() const
{
    std::printf("\n--- DecQ Time Analysis ---\n");
    std::printf("Total Time:                 %.4lf ms\n", stats.total_time / 1000.0);
    std::printf("Data Label/NLF Index:       %.4lf ms\n", stats.data_index_time / 1000.0);
    std::printf("Query Decomposition:        %.4lf ms\n", stats.decomposition_time / 1000.0);
    std::printf("Local Pattern Generation:   %.4lf ms\n", stats.local_pattern_time / 1000.0);
    std::printf("Local Exact Matching:       %.4lf ms\n", stats.local_matching_time / 1000.0);
    std::printf("Global Lattice Generation:  %.4lf ms\n", stats.lattice_time / 1000.0);
    std::printf("Minimal Sharing Merge:      %.4lf ms\n", stats.minimal_merge_time / 1000.0);
    std::printf("Non-minimal Edge Validation: %.4lf ms\n", stats.edge_validation_time / 1000.0);
    std::printf("Canonical Output Adapter:   %.4lf ms\n", stats.canonical_dedup_time / 1000.0);
    std::printf("Query Vertices:             %zu\n", stats.query_vertices);
    std::printf("Query Edges:                %zu\n", stats.query_edges);
    std::printf("Indexed Labels:             %zu\n", stats.indexed_labels);
    std::printf("Decomposition Nodes:        %zu\n", stats.decomposition_nodes);
    std::printf("Fragments:                  %zu\n", stats.decomposition_leaves);
    std::printf("Accepted Splits:            %zu\n", stats.accepted_splits);
    std::printf("Split Candidates:           %zu\n", stats.split_candidates);
    std::printf("Local Patterns:             %zu\n", stats.local_patterns);
    std::printf("Local Edge Cuts Skipped:    %zu\n", stats.local_edge_cuts_skipped);
    std::printf("Disconnected Local Patterns: %zu\n", stats.disconnected_local_patterns);
    std::printf("Local Components:           %zu\n", stats.local_components);
    std::printf("Exact Match Executions:     %zu\n", stats.exact_match_executions);
    std::printf("Local Search Partial Rows:  %zu\n", stats.local_search_partial_rows);
    std::printf("Local Match Rows:           %zu\n", stats.local_match_rows);
    std::printf("Lattice Candidates:         %zu\n", stats.lattice_candidates);
    std::printf("Global Patterns:            %zu\n", stats.global_patterns);
    std::printf("Global Edge Cuts Skipped:   %zu\n", stats.global_edge_cuts_skipped);
    std::printf("Minimal Patterns:           %zu\n", stats.minimal_patterns);
    std::printf("Non-minimal Patterns:       %zu\n", stats.nonminimal_patterns);
    std::printf("Minimal Patterns Processed: %zu\n", stats.minimal_patterns_processed);
    std::printf("Non-minimal Patterns Processed: %zu\n", stats.nonminimal_patterns_processed);
    std::printf("Decomposition Cache Hits:   %zu\n", stats.decomposition_cache_hits);
    std::printf("Decomposition Cache Misses: %zu\n", stats.decomposition_cache_misses);
    std::printf("Hash Joins:                 %zu\n", stats.hash_joins);
    std::printf("Hash Join Input Rows:       %zu\n", stats.hash_join_input_rows);
    std::printf("Hash Join Output Rows:      %zu\n", stats.hash_join_output_rows);
    std::printf("Edge Validations:           %zu\n", stats.edge_validations);
    std::printf("Materialized Match Rows:    %zu\n", stats.materialized_match_rows);
    std::printf("Peak Live Match Rows:       %zu\n", stats.peak_live_match_rows);
    std::printf("Global-Pattern Rows:        %zu\n", stats.global_pattern_rows);
    std::printf("Non-canonical Pattern Rows: %zu\n", stats.noncanonical_pattern_rows);
    std::printf("Duplicate Canonical Mappings: %zu\n", stats.duplicate_canonical_mappings);
    std::printf("Invalid Global Rows:        %zu\n", stats.invalid_global_rows);
    std::printf("Unique Results:             %zu\n", stats.result_count);
    std::printf("Edge Selectivity (min/avg/max): %.6g / %.6g / %.6g\n",
        stats.minimum_edge_selectivity, stats.average_edge_selectivity,
        stats.maximum_edge_selectivity);
    std::printf("Estimated Plan Cost (unsplit/DecQ): %.6g / %.6g\n",
        stats.unsplit_estimated_cost, stats.decomposition_estimated_cost);
#if MATCH_OUTPUT_LIMIT > 0
    std::printf("Output Limit:               %zu%s\n",
        static_cast<size_t>(MATCH_OUTPUT_LIMIT),
        stats.output_limit_reached ? " (reached)" : "");
#endif
#if DECQ_INTERMEDIATE_MATCH_LIMIT > 0
    std::printf("Intermediate Match Limit:   %zu%s\n",
        static_cast<size_t>(DECQ_INTERMEDIATE_MATCH_LIMIT),
        stats.intermediate_limit_reached ? " (reached)" : "");
#endif
    std::printf("------------------------------------------\n");
}

void Approximate_DecQ(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold)
{
    DecQSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(results);
    }
    long long preprocessing_time = solver.stats.data_index_time +
        solver.stats.decomposition_time + solver.stats.local_pattern_time +
        solver.stats.lattice_time;
    long long search_time = solver.stats.total_time > preprocessing_time
        ? solver.stats.total_time - preprocessing_time : 0;
    set_reported_phase_times(preprocessing_time, search_time);
    set_reported_result_count(solver.stats.result_count);
    solver.printStats();
}

namespace {

void run_decq(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold)
{
    Approximate_DecQ(query_graph, data_graph, results, threshold);
}

} // namespace

const AlgorithmDefinition &create_algorithm_definition()
{
    static const AlgorithmDefinition definition = {
        "decq",
        "DecQ",
        &run_decq
    };
    return definition;
}

} // namespace ssm
