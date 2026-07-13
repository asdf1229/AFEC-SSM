#include "matching/algorithms/sasum.h"

#include "utility/utility.h"

#include <functional>
#include <unordered_map>

namespace ssm {

struct SASUMSolver::Impl {
    SASUMStats &stats;

    explicit Impl(SASUMStats &stats_ref)
        : stats(stats_ref), query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr),
          threshold(0), qn(0), gn(0)
    {}

    bool init(const Graph *query, const Graph *data, ui match_threshold)
    {
        reset();
        Timer index_timer;
        query_graph = query;
        data_graph = data;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();

        if (qn == 0 || gn == 0 || qn > gn) {
            stats.data_index_time = index_timer.elapsed();
            stats.total_time = stats.data_index_time;
            return false;
        }

        buildQueryEdges();
        Pattern full_query;
        full_query.vertices.assign(qn, 1);
        full_query.edges.assign(query_edges.size(), 1);
        if (!isConnected(full_query)) {
            stats.data_index_time = index_timer.elapsed();
            stats.total_time = stats.data_index_time;
            return false;
        }
        size_t minimum_connected_edges = qn > 0 ? static_cast<size_t>(qn - 1) : 0;
        if (query_edges.size() >= minimum_connected_edges) {
            size_t cycle_rank = query_edges.size() - minimum_connected_edges;
            if (static_cast<size_t>(threshold) > cycle_rank) {
                threshold = static_cast<ui>(cycle_rank);
            }
        }
        buildDataIndex();
        stats.data_index_time = index_timer.elapsed();
        return true;
    }

    void match(MatchResults &results)
    {
        Timer total_timer;
        results_ptr = &results;
        results_ptr->clear();

        Timer phase_timer;
        buildQueryLattice();
        stats.lattice_time = phase_timer.elapsed();

        if (lattice_levels.empty() || lattice_levels[0].empty()) {
            stats.total_time = stats.data_index_time + total_timer.elapsed();
            return;
        }

        size_t terminal_distance = lattice_levels.size() - 1;
        while (terminal_distance > 0 && lattice_levels[terminal_distance].empty()) {
            terminal_distance--;
        }

        if (terminal_distance == 0) {
            stats.terminal_graph_count = 1;
            phase_timer.restart();
            MatchSet query_matches;
            exactMatch(lattice_levels[0][0], query_matches);
            stats.exact_matching_executions = 1;
            stats.exact_matching_time = phase_timer.elapsed();
            stats.query_pattern_match_rows = query_matches.size();
            if (!intermediateLimitReached()) addOutputs(query_matches);
            stats.total_time = stats.data_index_time + total_timer.elapsed();
            return;
        }

        const std::vector<Pattern> &terminals = lattice_levels[terminal_distance];
        stats.terminal_graph_count = terminals.size();
        if (terminals.empty()) {
            stats.total_time = stats.data_index_time + total_timer.elapsed();
            return;
        }

        phase_timer.restart();
        generateBaseGraphs(terminals);
        stats.base_generation_time = phase_timer.elapsed();

        phase_timer.restart();
        selectSeeds(terminals.size());
        stats.seed_selection_time = phase_timer.elapsed();

        phase_timer.restart();
        std::vector<MatchSet> seed_matches(selected_base_ids.size());
        for (size_t i = 0; i < selected_base_ids.size(); ++i) {
            exactMatch(base_graphs[selected_base_ids[i]].pattern, seed_matches[i]);
            stats.exact_matching_executions++;
            stats.seed_match_rows += seed_matches[i].size();
            if (intermediateLimitReached()) {
                stats.exact_matching_time = phase_timer.elapsed();
                stats.total_time = stats.data_index_time + total_timer.elapsed();
                return;
            }
        }
        stats.exact_matching_time = phase_timer.elapsed();

        phase_timer.restart();
        long long fallback_exact_time = 0;
        const std::function<void()> finish_derivation = [&]() {
            long long gross_derivation_time = phase_timer.elapsed();
            stats.derivation_time = gross_derivation_time > fallback_exact_time
                ? gross_derivation_time - fallback_exact_time : 0;
            stats.exact_matching_time += fallback_exact_time;
            stats.total_time = stats.data_index_time + total_timer.elapsed();
        };
        std::vector<MatchSet> current_matches(terminals.size());
        for (size_t terminal_id = 0; terminal_id < terminals.size(); ++terminal_id) {
            size_t best_seed = selected_base_ids.size();
            size_t best_size = std::numeric_limits<size_t>::max();

            for (size_t seed_id = 0; seed_id < selected_base_ids.size(); ++seed_id) {
                const BaseGraph &base = base_graphs[selected_base_ids[seed_id]];
                if (!coversTerminal(base, terminal_id)) continue;
                if (seed_matches[seed_id].size() < best_size) {
                    best_size = seed_matches[seed_id].size();
                    best_seed = seed_id;
                }
            }

            if (best_seed == selected_base_ids.size()) {
                // The greedy cover should make this path unreachable.  Keeping a
                // correctness fallback makes malformed or degenerate inputs safe.
                Timer fallback_timer;
                exactMatch(terminals[terminal_id], current_matches[terminal_id]);
                fallback_exact_time += fallback_timer.elapsed();
                stats.exact_matching_executions++;
                stats.fallback_exact_executions++;
            }
            else {
                const Pattern &base = base_graphs[selected_base_ids[best_seed]].pattern;
                deriveCoveredMatches(base, terminals[terminal_id],
                    seed_matches[best_seed], current_matches[terminal_id]);
            }
            if (intermediateLimitReached()) {
                finish_derivation();
                return;
            }
            stats.query_pattern_match_rows += current_matches[terminal_id].size();
            addOutputs(current_matches[terminal_id]);
            if (outputLimitReached()) {
                finish_derivation();
                return;
            }
        }
        releaseMatchSets(seed_matches);

        for (size_t distance = terminal_distance; distance > 0; --distance) {
            const std::vector<Pattern> &parents = lattice_levels[distance - 1];
            const std::vector<Pattern> &children = lattice_levels[distance];
            const std::unordered_map<std::string, size_t> &child_index =
                lattice_indices[distance];
            std::vector<MatchSet> parent_matches(parents.size());

            for (size_t parent_id = 0; parent_id < parents.size(); ++parent_id) {
                const Pattern &parent = parents[parent_id];
                size_t best_child = children.size();
                size_t best_size = std::numeric_limits<size_t>::max();
                size_t added_edge = query_edges.size();

                for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
                    if (!parent.edges[edge_id]) continue;
                    Pattern child = parent;
                    child.edges[edge_id] = 0;
                    std::unordered_map<std::string, size_t>::const_iterator it =
                        child_index.find(patternKey(child));
                    if (it == child_index.end()) continue;

                    size_t child_id = it->second;
                    if (current_matches[child_id].size() < best_size) {
                        best_size = current_matches[child_id].size();
                        best_child = child_id;
                        added_edge = edge_id;
                    }
                }

                if (best_child == children.size()) {
                    Timer fallback_timer;
                    exactMatch(parent, parent_matches[parent_id]);
                    fallback_exact_time += fallback_timer.elapsed();
                    stats.exact_matching_executions++;
                    stats.fallback_exact_executions++;
                }
                else {
                    filterByEdge(current_matches[best_child], query_edges[added_edge],
                        parent_matches[parent_id]);
                }

                if (intermediateLimitReached()) {
                    finish_derivation();
                    return;
                }

                stats.query_pattern_match_rows += parent_matches[parent_id].size();
                addOutputs(parent_matches[parent_id]);
                if (outputLimitReached()) {
                    finish_derivation();
                    return;
                }
            }

            current_matches.swap(parent_matches);
            releaseMatchSets(parent_matches);
        }

        finish_derivation();
    }

    void printStats() const
    {
        std::printf("\n--- SASUM Time Analysis ---\n");
        std::printf("Total Time:             %.4lf ms\n", stats.total_time / 1000.0);
        std::printf("Data Label/NLF Index:   %.4lf ms\n", stats.data_index_time / 1000.0);
        std::printf("Lattice Generation:     %.4lf ms\n", stats.lattice_time / 1000.0);
        std::printf("Base Graph Generation:  %.4lf ms\n", stats.base_generation_time / 1000.0);
        std::printf("Seed Selection:         %.4lf ms\n", stats.seed_selection_time / 1000.0);
        std::printf("Exact Matching:         %.4lf ms\n", stats.exact_matching_time / 1000.0);
        std::printf("Derivation + Dedup:     %.4lf ms\n", stats.derivation_time / 1000.0);
        std::printf("Query Subgraphs:        %zu\n", stats.query_subgraph_count);
        std::printf("Terminal Graphs:        %zu\n", stats.terminal_graph_count);
        std::printf("Base Graphs:            %zu\n", stats.base_graph_count);
        std::printf("Seed Graphs:            %zu\n", stats.seed_graph_count);
        std::printf("Exact Match Executions: %zu\n", stats.exact_matching_executions);
        std::printf("Fallback Exact Runs:    %zu\n", stats.fallback_exact_executions);
        std::printf("Exact Embeddings:       %zu\n", stats.exact_embeddings);
        std::printf("Seed Match Rows:        %zu\n", stats.seed_match_rows);
        std::printf("Query-Pattern Rows:     %zu\n", stats.query_pattern_match_rows);
        std::printf("Peak Live Match Rows:   %zu\n", stats.peak_live_match_rows);
        std::printf("Cross-Pattern Dups:     %zu\n", stats.duplicate_outputs);
        std::printf("Results Found:          %zu\n", stats.result_count);
#if MATCH_OUTPUT_LIMIT > 0
        std::printf("Output Limit:           %zu%s\n",
            static_cast<size_t>(MATCH_OUTPUT_LIMIT),
            stats.output_limit_reached ? " (reached)" : "");
#endif
#if SASUM_INTERMEDIATE_MATCH_LIMIT > 0
        std::printf("Intermediate Limit:     %zu%s\n",
            static_cast<size_t>(SASUM_INTERMEDIATE_MATCH_LIMIT),
            stats.intermediate_limit_reached ? " (reached)" : "");
#endif
        std::printf("------------------------------------------\n");
    }

    struct QueryEdge {
        ui u = 0;
        ui v = 0;
    };

    struct Pattern {
        std::vector<unsigned char> vertices;
        std::vector<unsigned char> edges;
    };

    typedef std::vector<int> Mapping;
    typedef std::vector<Mapping> MatchSet;

    struct BaseGraph {
        Pattern pattern;
        std::vector<size_t> covered_terminals;
    };

    struct ExactSearchState {
        std::vector<std::vector<ui> > adjacency;
        std::vector<std::vector<ui> > candidates;
        std::vector<ui> order;
        Mapping mapping;
        std::vector<unsigned char> used_data;
        MatchSet *output = nullptr;
    };

    const Graph *query_graph;
    const Graph *data_graph;
    MatchResults *results_ptr;
    ui threshold;
    ui qn;
    ui gn;

    std::vector<QueryEdge> query_edges;
    std::vector<std::vector<ui> > data_vertices_by_label;
    std::vector<std::vector<std::pair<LabelID, ui> > > data_nlf;
    std::vector<std::vector<Pattern> > lattice_levels;
    std::vector<std::unordered_map<std::string, size_t> > lattice_indices;
    std::vector<BaseGraph> base_graphs;
    std::vector<size_t> selected_base_ids;
    std::unordered_set<std::string> result_keys;
    size_t live_intermediate_rows = 0;

    void reset()
    {
        stats = SASUMStats();
        query_graph = nullptr;
        data_graph = nullptr;
        results_ptr = nullptr;
        threshold = 0;
        qn = 0;
        gn = 0;
        query_edges.clear();
        data_vertices_by_label.clear();
        data_nlf.clear();
        lattice_levels.clear();
        lattice_indices.clear();
        base_graphs.clear();
        selected_base_ids.clear();
        result_keys.clear();
        live_intermediate_rows = 0;
    }

    void buildQueryEdges()
    {
        query_edges.clear();
        for (ui u = 0; u < qn; ++u) {
            ui degree = 0;
            const ui *neighbors = query_graph->getVertexNeighbors(u, degree);
            for (ui i = 0; i < degree; ++i) {
                ui v = neighbors[i];
                if (u < v) query_edges.push_back(QueryEdge{u, v});
            }
        }
    }

    void buildDataIndex()
    {
        ui label_count = std::max(query_graph->getLabelsCount(),
            data_graph->getLabelsCount());
        data_vertices_by_label.assign(label_count, std::vector<ui>());
        data_nlf.assign(gn, std::vector<std::pair<LabelID, ui> >());

        for (ui v = 0; v < gn; ++v) {
            LabelID label = data_graph->getVertexLabel(v);
            if (label >= 0 && static_cast<ui>(label) < label_count) {
                data_vertices_by_label[static_cast<ui>(label)].push_back(v);
            }

            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(v, degree);
            std::vector<LabelID> labels;
            labels.reserve(degree);
            for (ui i = 0; i < degree; ++i) {
                labels.push_back(data_graph->getVertexLabel(neighbors[i]));
            }
            std::sort(labels.begin(), labels.end());
            for (size_t i = 0; i < labels.size();) {
                size_t j = i + 1;
                while (j < labels.size() && labels[j] == labels[i]) ++j;
                data_nlf[v].push_back(std::make_pair(labels[i], static_cast<ui>(j - i)));
                i = j;
            }
        }
    }

    std::string patternKey(const Pattern &pattern) const
    {
        std::string key;
        key.reserve(pattern.vertices.size() + pattern.edges.size());
        for (size_t i = 0; i < pattern.vertices.size(); ++i) {
            key.push_back(static_cast<char>(pattern.vertices[i]));
        }
        for (size_t i = 0; i < pattern.edges.size(); ++i) {
            key.push_back(static_cast<char>(pattern.edges[i]));
        }
        return key;
    }

    size_t activeVertexCount(const Pattern &pattern) const
    {
        size_t count = 0;
        for (size_t i = 0; i < pattern.vertices.size(); ++i) {
            count += pattern.vertices[i] != 0;
        }
        return count;
    }

    bool isConnected(const Pattern &pattern) const
    {
        size_t active_count = activeVertexCount(pattern);
        if (active_count == 0) return false;

        ui start = qn;
        for (ui u = 0; u < qn; ++u) {
            if (pattern.vertices[u]) {
                start = u;
                break;
            }
        }
        if (active_count == 1) return start != qn;

        std::vector<unsigned char> visited(qn, 0);
        std::queue<ui> queue;
        queue.push(start);
        visited[start] = 1;
        size_t visited_count = 0;

        while (!queue.empty()) {
            ui u = queue.front();
            queue.pop();
            visited_count++;

            for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
                if (!pattern.edges[edge_id]) continue;
                const QueryEdge &edge = query_edges[edge_id];
                ui v = qn;
                if (edge.u == u) v = edge.v;
                else if (edge.v == u) v = edge.u;
                if (v == qn || !pattern.vertices[v] || visited[v]) continue;
                visited[v] = 1;
                queue.push(v);
            }
        }

        return visited_count == active_count;
    }

    void buildQueryLattice()
    {
        lattice_levels.assign(static_cast<size_t>(threshold) + 1,
            std::vector<Pattern>());
        lattice_indices.assign(static_cast<size_t>(threshold) + 1,
            std::unordered_map<std::string, size_t>());

        Pattern query;
        query.vertices.assign(qn, 1);
        query.edges.assign(query_edges.size(), 1);
        lattice_levels[0].push_back(query);
        lattice_indices[0][patternKey(query)] = 0;

        for (ui distance = 1; distance <= threshold; ++distance) {
            const std::vector<Pattern> &previous = lattice_levels[distance - 1];
            std::vector<Pattern> &level = lattice_levels[distance];
            std::unordered_map<std::string, size_t> &index = lattice_indices[distance];

            for (size_t pattern_id = 0; pattern_id < previous.size(); ++pattern_id) {
                const Pattern &parent = previous[pattern_id];
                for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
                    if (!parent.edges[edge_id]) continue;
                    Pattern child = parent;
                    child.edges[edge_id] = 0;
                    if (!isConnected(child)) continue;
                    std::string key = patternKey(child);
                    if (index.find(key) != index.end()) continue;
                    index[key] = level.size();
                    level.push_back(child);
                }
            }
        }

        stats.query_subgraph_count = 0;
        for (size_t i = 0; i < lattice_levels.size(); ++i) {
            stats.query_subgraph_count += lattice_levels[i].size();
        }
    }

    bool pruneEdge(const Pattern &terminal, size_t edge_id, Pattern &base) const
    {
        base = terminal;
        base.edges[edge_id] = 0;
        if (isConnected(base)) return true;

        std::vector<ui> degree(qn, 0);
        for (size_t i = 0; i < query_edges.size(); ++i) {
            if (!base.edges[i]) continue;
            degree[query_edges[i].u]++;
            degree[query_edges[i].v]++;
        }

        std::vector<ui> isolated;
        for (ui u = 0; u < qn; ++u) {
            if (base.vertices[u] && degree[u] == 0) isolated.push_back(u);
        }

        if (isolated.size() == 1) {
            base.vertices[isolated[0]] = 0;
            return isConnected(base);
        }

        if (activeVertexCount(terminal) == 2 && isolated.size() == 2) {
            // Definition 6 permits an arbitrary choice for a one-edge graph.
            base.vertices[isolated[1]] = 0;
            return isConnected(base);
        }

        return false;
    }

    void generateBaseGraphs(const std::vector<Pattern> &terminals)
    {
        base_graphs.clear();
        std::unordered_map<std::string, size_t> index;

        // All terminals are in one edge-count layer and every base has exactly
        // one fewer edge.  Therefore recording every pruning provenance is
        // equivalent to scanning all terminals for the paper's b subset-of t
        // cover relation.

        for (size_t terminal_id = 0; terminal_id < terminals.size(); ++terminal_id) {
            const Pattern &terminal = terminals[terminal_id];
            for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
                if (!terminal.edges[edge_id]) continue;
                Pattern base;
                if (!pruneEdge(terminal, edge_id, base)) continue;

                std::string key = patternKey(base);
                std::unordered_map<std::string, size_t>::iterator found = index.find(key);
                size_t base_id;
                if (found == index.end()) {
                    base_id = base_graphs.size();
                    index[key] = base_id;
                    BaseGraph item;
                    item.pattern = base;
                    base_graphs.push_back(item);
                }
                else {
                    base_id = found->second;
                }

                std::vector<size_t> &covered = base_graphs[base_id].covered_terminals;
                if (covered.empty() || covered.back() != terminal_id) {
                    covered.push_back(terminal_id);
                }
            }
        }

        stats.base_graph_count = base_graphs.size();
    }

    void selectSeeds(size_t terminal_count)
    {
        selected_base_ids.clear();
        std::vector<unsigned char> covered(terminal_count, 0);
        size_t covered_count = 0;

        while (covered_count < terminal_count) {
            size_t best_base = base_graphs.size();
            size_t best_gain = 0;

            for (size_t base_id = 0; base_id < base_graphs.size(); ++base_id) {
                size_t gain = 0;
                const std::vector<size_t> &cover = base_graphs[base_id].covered_terminals;
                for (size_t i = 0; i < cover.size(); ++i) {
                    gain += covered[cover[i]] == 0;
                }
                if (gain > best_gain) {
                    best_gain = gain;
                    best_base = base_id;
                }
            }

            if (best_base == base_graphs.size() || best_gain == 0) break;
            selected_base_ids.push_back(best_base);
            const std::vector<size_t> &cover = base_graphs[best_base].covered_terminals;
            for (size_t i = 0; i < cover.size(); ++i) {
                if (!covered[cover[i]]) {
                    covered[cover[i]] = 1;
                    covered_count++;
                }
            }
        }

        stats.seed_graph_count = selected_base_ids.size();
    }

    bool coversTerminal(const BaseGraph &base, size_t terminal_id) const
    {
        return std::binary_search(base.covered_terminals.begin(),
            base.covered_terminals.end(), terminal_id);
    }

    void buildPatternAdjacency(const Pattern &pattern,
        std::vector<std::vector<ui> > &adjacency) const
    {
        adjacency.assign(qn, std::vector<ui>());
        for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
            if (!pattern.edges[edge_id]) continue;
            const QueryEdge &edge = query_edges[edge_id];
            adjacency[edge.u].push_back(edge.v);
            adjacency[edge.v].push_back(edge.u);
        }
    }

    bool nlfFits(const std::vector<std::pair<LabelID, ui> > &required,
        const std::vector<std::pair<LabelID, ui> > &available) const
    {
        size_t available_id = 0;
        for (size_t i = 0; i < required.size(); ++i) {
            while (available_id < available.size() &&
                available[available_id].first < required[i].first) {
                available_id++;
            }
            if (available_id == available.size() ||
                available[available_id].first != required[i].first ||
                available[available_id].second < required[i].second) {
                return false;
            }
        }
        return true;
    }

    bool buildExactCandidates(const Pattern &pattern, ExactSearchState &state) const
    {
        buildPatternAdjacency(pattern, state.adjacency);
        state.candidates.assign(qn, std::vector<ui>());

        for (ui u = 0; u < qn; ++u) {
            if (!pattern.vertices[u]) continue;
            LabelID label = query_graph->getVertexLabel(u);
            if (label < 0 || static_cast<size_t>(label) >= data_vertices_by_label.size()) {
                return false;
            }

            std::vector<LabelID> neighbor_labels;
            neighbor_labels.reserve(state.adjacency[u].size());
            for (size_t i = 0; i < state.adjacency[u].size(); ++i) {
                neighbor_labels.push_back(
                    query_graph->getVertexLabel(state.adjacency[u][i]));
            }
            std::sort(neighbor_labels.begin(), neighbor_labels.end());
            std::vector<std::pair<LabelID, ui> > required;
            for (size_t i = 0; i < neighbor_labels.size();) {
                size_t j = i + 1;
                while (j < neighbor_labels.size() && neighbor_labels[j] == neighbor_labels[i]) ++j;
                required.push_back(std::make_pair(neighbor_labels[i], static_cast<ui>(j - i)));
                i = j;
            }

            const std::vector<ui> &bucket = data_vertices_by_label[static_cast<ui>(label)];
            for (size_t i = 0; i < bucket.size(); ++i) {
                ui v = bucket[i];
                if (data_graph->getVertexDegree(v) < state.adjacency[u].size()) continue;
                if (!nlfFits(required, data_nlf[v])) continue;
                state.candidates[u].push_back(v);
            }
            if (state.candidates[u].empty()) return false;
        }

        return true;
    }

    void buildExactOrder(const Pattern &pattern, ExactSearchState &state) const
    {
        state.order.clear();
        std::vector<unsigned char> selected(qn, 0);
        size_t active_count = activeVertexCount(pattern);

        ui root = qn;
        for (ui u = 0; u < qn; ++u) {
            if (!pattern.vertices[u]) continue;
            if (root == qn || state.candidates[u].size() < state.candidates[root].size() ||
                (state.candidates[u].size() == state.candidates[root].size() &&
                 state.adjacency[u].size() > state.adjacency[root].size())) {
                root = u;
            }
        }
        if (root == qn) return;
        state.order.push_back(root);
        selected[root] = 1;

        while (state.order.size() < active_count) {
            ui best = qn;
            size_t best_mapped_neighbors = 0;
            for (ui u = 0; u < qn; ++u) {
                if (!pattern.vertices[u] || selected[u]) continue;
                size_t mapped_neighbors = 0;
                for (size_t i = 0; i < state.adjacency[u].size(); ++i) {
                    mapped_neighbors += selected[state.adjacency[u][i]] != 0;
                }
                if (best == qn || mapped_neighbors > best_mapped_neighbors ||
                    (mapped_neighbors == best_mapped_neighbors &&
                     state.candidates[u].size() < state.candidates[best].size()) ||
                    (mapped_neighbors == best_mapped_neighbors &&
                     state.candidates[u].size() == state.candidates[best].size() &&
                     state.adjacency[u].size() > state.adjacency[best].size())) {
                    best = u;
                    best_mapped_neighbors = mapped_neighbors;
                }
            }
            if (best == qn) break;
            selected[best] = 1;
            state.order.push_back(best);
        }
    }

    bool exactCandidateValid(ui u, ui v, const ExactSearchState &state) const
    {
        if (state.used_data[v]) return false;
        if (!std::binary_search(state.candidates[u].begin(), state.candidates[u].end(), v)) {
            return false;
        }

        for (size_t i = 0; i < state.adjacency[u].size(); ++i) {
            ui query_neighbor = state.adjacency[u][i];
            int data_neighbor = state.mapping[query_neighbor];
            if (data_neighbor >= 0 && !data_graph->hasEdge(v, static_cast<ui>(data_neighbor))) {
                return false;
            }
        }
        return true;
    }

    bool intermediateLimitReached() const
    {
#if SASUM_INTERMEDIATE_MATCH_LIMIT > 0
        return stats.intermediate_limit_reached;
#else
        return false;
#endif
    }

    bool reserveIntermediateRow()
    {
#if SASUM_INTERMEDIATE_MATCH_LIMIT > 0
        if (live_intermediate_rows >=
            static_cast<size_t>(SASUM_INTERMEDIATE_MATCH_LIMIT)) {
            stats.intermediate_limit_reached = true;
            return false;
        }
#endif
        live_intermediate_rows++;
        stats.peak_live_match_rows = std::max(
            stats.peak_live_match_rows, live_intermediate_rows);
        return true;
    }

    void releaseMatchSet(MatchSet &matches)
    {
        assert(matches.size() <= live_intermediate_rows);
        live_intermediate_rows -= matches.size();
        MatchSet().swap(matches);
    }

    void releaseMatchSets(std::vector<MatchSet> &match_sets)
    {
        for (size_t i = 0; i < match_sets.size(); ++i) {
            releaseMatchSet(match_sets[i]);
        }
        std::vector<MatchSet>().swap(match_sets);
    }

    void exactDfs(size_t depth, ExactSearchState &state)
    {
        if (intermediateLimitReached()) return;
        if (depth == state.order.size()) {
            if (!reserveIntermediateRow()) return;
            state.output->push_back(state.mapping);
            stats.exact_embeddings++;
            return;
        }

        ui u = state.order[depth];
        ui anchor = qn;
        ui anchor_data_degree = std::numeric_limits<ui>::max();
        for (size_t i = 0; i < state.adjacency[u].size(); ++i) {
            ui neighbor = state.adjacency[u][i];
            if (state.mapping[neighbor] < 0) continue;
            ui data_degree = data_graph->getVertexDegree(
                static_cast<ui>(state.mapping[neighbor]));
            if (anchor == qn || data_degree < anchor_data_degree) {
                anchor = neighbor;
                anchor_data_degree = data_degree;
            }
        }

        if (anchor == qn) {
            const std::vector<ui> &candidates = state.candidates[u];
            for (size_t i = 0; i < candidates.size(); ++i) {
                ui v = candidates[i];
                if (!exactCandidateValid(u, v, state)) continue;
                state.mapping[u] = static_cast<int>(v);
                state.used_data[v] = 1;
                exactDfs(depth + 1, state);
                state.used_data[v] = 0;
                state.mapping[u] = -1;
                if (intermediateLimitReached()) return;
            }
        }
        else {
            ui mapped_anchor = static_cast<ui>(state.mapping[anchor]);
            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(mapped_anchor, degree);
            for (ui i = 0; i < degree; ++i) {
                ui v = neighbors[i];
                if (!exactCandidateValid(u, v, state)) continue;
                state.mapping[u] = static_cast<int>(v);
                state.used_data[v] = 1;
                exactDfs(depth + 1, state);
                state.used_data[v] = 0;
                state.mapping[u] = -1;
                if (intermediateLimitReached()) return;
            }
        }
    }

    void exactMatch(const Pattern &pattern, MatchSet &output)
    {
        releaseMatchSet(output);
        ExactSearchState state;
        state.output = &output;
        if (!buildExactCandidates(pattern, state)) return;
        buildExactOrder(pattern, state);
        if (state.order.size() != activeVertexCount(pattern)) return;
        state.mapping.assign(qn, -1);
        state.used_data.assign(gn, 0);
        exactDfs(0, state);
    }

    void filterByEdge(const MatchSet &input, const QueryEdge &edge, MatchSet &output)
    {
        releaseMatchSet(output);
        for (size_t i = 0; i < input.size(); ++i) {
            const Mapping &mapping = input[i];
            if (mapping[edge.u] < 0 || mapping[edge.v] < 0) continue;
            if (data_graph->hasEdge(static_cast<ui>(mapping[edge.u]),
                    static_cast<ui>(mapping[edge.v]))) {
                if (!reserveIntermediateRow()) return;
                output.push_back(mapping);
            }
        }
    }

    void deriveCoveredMatches(const Pattern &base, const Pattern &terminal,
        const MatchSet &input, MatchSet &output)
    {
        size_t added_edge_id = query_edges.size();
        for (size_t edge_id = 0; edge_id < query_edges.size(); ++edge_id) {
            if (terminal.edges[edge_id] && !base.edges[edge_id]) {
                added_edge_id = edge_id;
                break;
            }
        }
        if (added_edge_id == query_edges.size()) {
            releaseMatchSet(output);
            return;
        }

        ui added_vertex = qn;
        for (ui u = 0; u < qn; ++u) {
            if (terminal.vertices[u] && !base.vertices[u]) {
                added_vertex = u;
                break;
            }
        }

        const QueryEdge &added_edge = query_edges[added_edge_id];
        if (added_vertex == qn) {
            filterByEdge(input, added_edge, output);
            return;
        }

        ui anchor = added_edge.u == added_vertex ? added_edge.v : added_edge.u;
        releaseMatchSet(output);
        for (size_t match_id = 0; match_id < input.size(); ++match_id) {
            const Mapping &mapping = input[match_id];
            if (mapping[anchor] < 0) continue;
            ui mapped_anchor = static_cast<ui>(mapping[anchor]);
            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(mapped_anchor, degree);
            for (ui i = 0; i < degree; ++i) {
                ui candidate = neighbors[i];
                if (data_graph->getVertexLabel(candidate) !=
                    query_graph->getVertexLabel(added_vertex)) continue;

                bool used = false;
                for (ui u = 0; u < qn; ++u) {
                    if (mapping[u] == static_cast<int>(candidate)) {
                        used = true;
                        break;
                    }
                }
                if (used) continue;

                Mapping extended = mapping;
                extended[added_vertex] = static_cast<int>(candidate);
                if (!reserveIntermediateRow()) return;
                output.push_back(extended);
            }
        }
    }

    std::string mappingKey(const Mapping &mapping) const
    {
        return std::string(reinterpret_cast<const char *>(&mapping[0]),
            mapping.size() * sizeof(mapping[0]));
    }

    bool outputLimitReached() const
    {
#if MATCH_OUTPUT_LIMIT > 0
        return stats.output_limit_reached;
#else
        return false;
#endif
    }

    void addOutputs(const MatchSet &matches)
    {
        for (size_t i = 0; i < matches.size(); ++i) {
            if (outputLimitReached()) return;
            const Mapping &mapping = matches[i];
            bool complete = true;
            for (ui u = 0; u < qn; ++u) {
                if (mapping[u] < 0) {
                    complete = false;
                    break;
                }
            }
            if (!complete) continue;

            std::string key = mappingKey(mapping);
            if (!result_keys.insert(key).second) {
                stats.duplicate_outputs++;
                continue;
            }

#ifndef NDEBUG
            std::vector<std::pair<ui, ui> > result;
            result.reserve(qn);
            for (ui u = 0; u < qn; ++u) {
                result.push_back(std::make_pair(u, static_cast<ui>(mapping[u])));
            }
            results_ptr->push_back(result);
#endif

            stats.result_count++;

#if MATCH_OUTPUT_LIMIT > 0
            if (stats.result_count >= static_cast<size_t>(MATCH_OUTPUT_LIMIT)) {
                stats.output_limit_reached = true;
                return;
            }
#endif
        }
    }
};

SASUMSolver::SASUMSolver() : impl_(new Impl(stats)) {}

SASUMSolver::~SASUMSolver() = default;

bool SASUMSolver::init(const Graph *query_graph, const Graph *data_graph,
    ui threshold)
{
    return impl_->init(query_graph, data_graph, threshold);
}

void SASUMSolver::match(MatchResults &results)
{
    impl_->match(results);
}

void SASUMSolver::printStats() const
{
    impl_->printStats();
}

void Approximate_SASUM(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold)
{
    SASUMSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(results);
    }
    long long preprocessing_time = solver.stats.data_index_time +
        solver.stats.lattice_time + solver.stats.base_generation_time +
        solver.stats.seed_selection_time;
    long long search_time = solver.stats.total_time > preprocessing_time
        ? solver.stats.total_time - preprocessing_time : 0;
    set_reported_phase_times(preprocessing_time, search_time);
    solver.printStats();
    set_reported_result_count(solver.stats.result_count);
}

namespace {

void run_sasum(const Graph *query_graph, const Graph *data_graph,
    MatchResults &results, ui threshold)
{
    Approximate_SASUM(query_graph, data_graph, results, threshold);
}

} // namespace

const AlgorithmDefinition &create_algorithm_definition()
{
    static const AlgorithmDefinition definition = {
        "sasum",
        "SASUM",
        &run_sasum
    };
    return definition;
}

} // namespace ssm
