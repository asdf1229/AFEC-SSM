#ifndef MATCHING_ALGORITHMS_CDE_MATCH_H_
#define MATCHING_ALGORITHMS_CDE_MATCH_H_

#include "graph/graph.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
#include <limits>

using namespace std;

// ============================================================================
// MatchingSolver Implementation
// ============================================================================

class MatchingSolver {
public:
    MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

    bool init(const Graph *q, const Graph *g, ui match_threshold)
    {
        Timer t_init;
        t_init.restart();

        query_graph = q;
        data_graph = g;
        threshold = match_threshold;
        qn = query_graph->getVerticesCount();
        gn = data_graph->getVerticesCount();
        label_count = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());

        if (qn == 0 || gn == 0) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        max_g_deg = data_graph->getMaxDegree();

        resetState();

        q_matrix.assign(qn, vector<char>(qn, 0));
        q_neighbors.assign(qn, vector<ui>());

        for (ui u = 0; u < qn; ++u) {
            ui deg = 0;
            const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) {
                ui nbr = nbrs[i];
                q_matrix[u][nbr] = 1;
                q_neighbors[u].push_back(nbr);
            }
        }

        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);

        Timer t_filter;
        bool res = runCandidateFiltering();
        stats.filter_time = t_filter.elapsed();
        if (!res) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        stats.init_time = t_init.elapsed();
        return true;
    }

    void match(vector<vector<pair<ui, ui>>> &results)
    {
        Timer t_search;
        t_search.restart();

        results_ptr = &results;
        results_ptr->clear();

        BranchSelector branch_selector(*this);
        ui root = branch_selector.selectInitialRoot();

        for (ui v0 : candidates[root]) {
            mapped_q[root] = (int)v0;
            mapped_g[v0] = (int)root;
            in_Mq[root] = 1;
            matched_count = 1;
            part_M.push_back({ root, v0 });

            onVertexMatchStateChanged(root, true);

            dfs(0, (int)root, nullptr);

            onVertexMatchStateChanged(root, false);

            mapped_q[root] = -1;
            mapped_g[v0] = -1;
            in_Mq[root] = 0;
            matched_count = 0;
            part_M.pop_back();
        }

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct TimeStats {
        long long total_time = 0;
        // init breakdown
        long long init_time = 0;
        long long filter_time = 0;
        // search breakdown
        long long dfs_time = 0;
        long long branch_time = 0;   // candidate enumeration & matching in dfs
        long long lb_time = 0;       // computeLowerBound
        long long frontier_time = 0; // building U_frontier
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        printf("\n--- CDE-Match Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
        printf("Results Found:       %zu\n", results_ptr ? results_ptr->size() : 0);
        printf("-----------------------------------------------------------\n");
    }

private:
    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;
    ui threshold;
    ui qn, gn;
    ui label_count;
    ui max_g_deg;

    vector<vector<char>> q_matrix;
    vector<vector<ui>> q_neighbors;

    vector<MyBitset> candidates;
    vector<MyBitset> x_cand;

    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui> Lq_degrees, Lg_degrees;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<char> in_Mq;
    vector<vector<char>> is_excluded;
    vector<char> in_P;
    ui matched_count;
    vector<pair<ui, ui>> part_M;

    vector<ui> anchor_count;
    vector<int> frontier_pos;
    vector<ui> active_frontier;
    struct SearchFrontierState {
        vector<ui> component_vertices;
        vector<ui> component_frontier;
    };
    vector<ui> frontier_visit_mark;
    ui frontier_visit_token;
    vector<int> lb_match_right;
    vector<ui> lb_seen_right;
    ui lb_seen_token;
    vector<ui> lb_data_frontier_mark;
    ui lb_data_frontier_token;

    inline void updateFrontierStatus(ui u)
    {
        bool should_be = (!in_Mq[u] && !in_P[u] && anchor_count[u] > 0);
        bool is_in = (frontier_pos[u] != -1);

        if (should_be && !is_in) {
            frontier_pos[u] = active_frontier.size();
            active_frontier.push_back(u);
        }
        else if (!should_be && is_in) {
            ui idx = frontier_pos[u];
            ui last_u = active_frontier.back();
            active_frontier[idx] = last_u;
            frontier_pos[last_u] = idx;
            active_frontier.pop_back();
            frontier_pos[u] = -1;
        }
    }

    void onVertexMatchStateChanged(ui u, bool matched)
    {
        if (matched) {
            updateFrontierStatus(u);
        }

        for (ui nbr : q_neighbors[u]) {
            if (!is_excluded[nbr][u]) {
                if (matched) {
                    anchor_count[nbr]++;
                }
                else {
                    anchor_count[nbr]--;
                }
            }
            updateFrontierStatus(nbr);
        }

        if (!matched) {
            updateFrontierStatus(u);
        }
    }

    void resetState()
    {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        in_Mq.assign(qn, 0);
        is_excluded.assign(qn, vector<char>(qn, 0));
        in_P.assign(qn, 0);
        matched_count = 0;
        part_M.clear();
        part_M.reserve(qn);

        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        x_cand.clear();
        x_cand.assign(qn, MyBitset(gn));

        q_matrix.clear();
        q_neighbors.clear();

        anchor_count.assign(qn, 0);
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
        frontier_visit_mark.assign(qn, 0);
        frontier_visit_token = 1;
        lb_match_right.assign(gn, -1);
        lb_seen_right.assign(gn, 0);
        lb_seen_token = 1;
        lb_data_frontier_mark.assign(gn, 0);
        lb_data_frontier_token = 1;

        stats = TimeStats();
    }

    void initGlobalLabelCounts(const Graph *g, vector<vector<ui>> &counts, vector<ui> &degrees)
    {
        ui n = g->getVerticesCount();
        counts.assign(n, vector<ui>(label_count, 0));
        degrees.assign(n, 0);

        for (ui u = 0; u < n; ++u) {
            ui deg = 0;
            const ui *neighbors = g->getVertexNeighbors(u, deg);
            degrees[u] = deg;
            for (ui i = 0; i < deg; ++i) {
                LabelID label = g->getVertexLabel(neighbors[i]);
                if (label >= 0 && (ui)label < label_count) {
                    counts[u][(ui)label]++;
                }
            }
        }
    }

    // ========================================================================
    // Filtering
    // ========================================================================
    struct CandidateFilter {
        MatchingSolver &solver;
        vector<vector<ui>> spoke_adj_by_query_nbr;
        vector<int> spoke_match_right;
        vector<bool> spoke_seen_right;
        MyBitset used_one_hop_vertices;

        explicit CandidateFilter(MatchingSolver &solver)
            : solver(solver), used_one_hop_vertices((int)solver.gn)
        {
            spoke_adj_by_query_nbr.assign(solver.qn, vector<ui>());
            spoke_match_right.assign(solver.max_g_deg, -1);
            spoke_seen_right.assign(solver.max_g_deg, false);
        }

        bool run()
        {
            if (!applyNeighborhoodLabelFilter()) return false;
            if (!applySpokeFilter()) return false;
            if (!applyOneHopFilter()) return false;

#ifndef NDEBUG
            printCandidateStatistics();
#endif
            return true;
        }

    private:
        ui computeNeighborhoodLabelGap(ui u, ui v) const
        {
            ui diff = 0;
            size_t sz = solver.Lq_counts[u].size();
            for (size_t i = 0; i < sz; ++i) {
                if (solver.Lq_counts[u][i] > solver.Lg_counts[v][i]) {
                    diff += (solver.Lq_counts[u][i] - solver.Lg_counts[v][i]);
                }
            }
            return diff;
        }

        bool applyNeighborhoodLabelFilter()
        {
            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = solver.query_graph->getVertexLabel(u);
                for (ui v = 0; v < solver.gn; ++v) {
                    if (lu != solver.data_graph->getVertexLabel(v)) continue;
                    if (solver.Lq_degrees[u] > solver.Lg_degrees[v] + solver.threshold) continue;
                    if (computeNeighborhoodLabelGap(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        bool findSpokeAugmentPath(ui left_idx, const vector<vector<ui>> &adj)
        {
            for (ui right_idx : adj[left_idx]) {
                if (spoke_seen_right[right_idx]) continue;
                spoke_seen_right[right_idx] = true;
                if (spoke_match_right[right_idx] < 0 ||
                    findSpokeAugmentPath((ui)spoke_match_right[right_idx], adj)) {
                    spoke_match_right[right_idx] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        ui computeMaximumSpokeMatching(const vector<vector<ui>> &adj, ui left_size, ui right_size)
        {
            ui mu = 0;
            std::fill(spoke_match_right.begin(), spoke_match_right.begin() + right_size, -1);
            for (ui i = 0; i < left_size; ++i) {
                std::fill(spoke_seen_right.begin(), spoke_seen_right.begin() + right_size, false);
                if (findSpokeAugmentPath(i, adj)) mu++;
            }
            return mu;
        }

        ui computeSpokeLowerBound(ui u, ui v)
        {
            const vector<ui> &query_neighbors = solver.q_neighbors[u];
            ui data_deg = 0;
            const ui *data_neighbors = solver.data_graph->getVertexNeighbors(v, data_deg);

            for (ui i = 0; i < query_neighbors.size(); ++i) {
                spoke_adj_by_query_nbr[i].clear();
                ui query_nbr = query_neighbors[i];
                for (ui j = 0; j < data_deg; ++j) {
                    ui data_nbr = data_neighbors[j];
                    if (solver.candidates[query_nbr].contains(data_nbr)) {
                        spoke_adj_by_query_nbr[i].push_back(j);
                    }
                }
            }

            ui match_size = computeMaximumSpokeMatching(spoke_adj_by_query_nbr, (ui)query_neighbors.size(), data_deg);
            return (ui)query_neighbors.size() - match_size;
        }

        bool applySpokeFilter()
        {
            queue<ui> pending_vertices;
            vector<char> in_queue(solver.qn, 1);

            for (ui u = 0; u < solver.qn; ++u) pending_vertices.push(u);

            while (!pending_vertices.empty()) {
                ui u = pending_vertices.front();
                pending_vertices.pop();
                in_queue[u] = 0;

                vector<ui> to_remove;
                for (ui v : solver.candidates[u]) {
                    if (computeSpokeLowerBound(u, v) > solver.threshold) {
                        to_remove.push_back(v);
                    }
                }

                if (to_remove.empty()) continue;

                for (ui v : to_remove) solver.candidates[u].remove(v);
                if (solver.candidates[u].empty()) return false;

                for (ui nbr : solver.q_neighbors[u]) {
                    if (!in_queue[nbr]) {
                        pending_vertices.push(nbr);
                        in_queue[nbr] = 1;
                    }
                }
            }
            return true;
        }

        ui countNeighborInnerEdges(ui u) const
        {
            ui count = 0;
            const vector<ui> &query_neighbors = solver.q_neighbors[u];
            for (ui i = 0; i < query_neighbors.size(); ++i) {
                for (ui j = i + 1; j < query_neighbors.size(); ++j) {
                    if (solver.q_matrix[query_neighbors[i]][query_neighbors[j]]) count++;
                }
            }
            return count;
        }

        vector<ui> buildOneHopSearchOrder(const vector<ui> &query_neighbors,
            const vector<vector<ui>> &neighbor_candidates) const
        {
            vector<ui> ord(query_neighbors.size());
            iota(ord.begin(), ord.end(), 0);

            sort(ord.begin(), ord.end(), [&](ui a, ui b) {
                if (neighbor_candidates[a].size() != neighbor_candidates[b].size()) {
                    return neighbor_candidates[a].size() < neighbor_candidates[b].size();
                }
                ui deg_a = 0, deg_b = 0;
                for (ui z : query_neighbors) {
                    if (solver.q_matrix[query_neighbors[a]][z]) deg_a++;
                    if (solver.q_matrix[query_neighbors[b]][z]) deg_b++;
                }
                return deg_a > deg_b;
                });
            return ord;
        }

        ui computeResidualOneHopLowerBound(ui pos, const vector<ui> &ord,
            const vector<vector<ui>> &neighbor_candidates) const
        {
            ui rem = 0;
            for (ui k = pos; k < ord.size(); ++k) {
                ui i = ord[k];
                bool has_free = false;
                for (ui y : neighbor_candidates[i]) {
                    if (!used_one_hop_vertices.contains(y)) {
                        has_free = true;
                        break;
                    }
                }
                if (!has_free) rem++;
            }
            return rem;
        }

        bool searchOneHopAssignment(ui pos, const vector<ui> &query_neighbors,
            const vector<vector<ui>> &neighbor_candidates,
            vector<int> &partial_assignment,
            const vector<ui> &ord,
            ui current_cost)
        {
            if (current_cost > solver.threshold) return false;
            if (pos == ord.size()) return true;

            ui rem_lb = computeResidualOneHopLowerBound(pos, ord, neighbor_candidates);
            if (current_cost + rem_lb > solver.threshold) return false;

            ui i = ord[pos];
            ui x = query_neighbors[i];

            partial_assignment[i] = -1;
            if (searchOneHopAssignment(pos + 1, query_neighbors, neighbor_candidates,
                partial_assignment, ord, current_cost + 1)) {
                partial_assignment[i] = -2;
                return true;
            }

            for (ui y : neighbor_candidates[i]) {
                if (used_one_hop_vertices.contains(y)) continue;

                ui delta_inner = 0;
                for (ui j = 0; j < query_neighbors.size(); ++j) {
                    if (partial_assignment[j] < 0) continue;
                    ui z = query_neighbors[j];
                    if (!solver.q_matrix[x][z]) continue;

                    ui yz = (ui)partial_assignment[j];
                    if (!solver.data_graph->hasEdge(y, yz)) {
                        delta_inner++;
                    }
                }

                if (current_cost + delta_inner > solver.threshold) continue;

                used_one_hop_vertices.insert(y);
                partial_assignment[i] = (int)y;

                if (searchOneHopAssignment(pos + 1, query_neighbors, neighbor_candidates,
                    partial_assignment, ord, current_cost + delta_inner)) {
                    used_one_hop_vertices.remove(y);
                    partial_assignment[i] = -2;
                    return true;
                }

                used_one_hop_vertices.remove(y);
                partial_assignment[i] = -2;
            }

            partial_assignment[i] = -2;
            return false;
        }

        bool passesOneHopFilter(ui u, ui v)
        {
            const vector<ui> &query_neighbors = solver.q_neighbors[u];
            vector<vector<ui>> neighbor_candidates(query_neighbors.size());

            ui data_deg = 0;
            const ui *data_neighbors = solver.data_graph->getVertexNeighbors(v, data_deg);
            for (ui i = 0; i < query_neighbors.size(); ++i) {
                ui query_nbr = query_neighbors[i];
                for (ui j = 0; j < data_deg; ++j) {
                    ui y = data_neighbors[j];
                    if (solver.candidates[query_nbr].contains(y)) {
                        neighbor_candidates[i].push_back(y);
                    }
                }
            }

            vector<ui> ord = buildOneHopSearchOrder(query_neighbors, neighbor_candidates);
            vector<int> partial_assignment(query_neighbors.size(), -2);
            return searchOneHopAssignment(0, query_neighbors, neighbor_candidates, partial_assignment, ord, 0);
        }

        bool applyOneHopFilter()
        {
            vector<pair<ui, ui>> to_delete;

            for (ui u = 0; u < solver.qn; ++u) {
                ui deg_u = solver.q_neighbors[u].size();
                if (deg_u <= 1) continue;
                if (countNeighborInnerEdges(u) == 0) continue;

                for (ui v : solver.candidates[u]) {
                    if (!passesOneHopFilter(u, v)) {
                        to_delete.push_back({ u, v });
                    }
                }
            }

            for (auto p : to_delete) {
                solver.candidates[p.first].remove(p.second);
            }

            for (ui u = 0; u < solver.qn; ++u) {
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

#ifndef NDEBUG
        void printCandidateStatistics()
        {
            vector<ui> missing_edges_dist(solver.threshold + 1, 0);
            vector<vector<ui>> vertex_missing_edges_dist(solver.qn, vector<ui>(solver.threshold + 1, 0));
            ui total_candidates_count = 0;

            for (ui u = 0; u < solver.qn; ++u) {
                total_candidates_count += solver.candidates[u].size();

                for (ui v : solver.candidates[u]) {
                    ui min_missing_edges = computeSpokeLowerBound(u, v);
                    if (min_missing_edges <= solver.threshold) {
                        missing_edges_dist[min_missing_edges]++;
                        vertex_missing_edges_dist[u][min_missing_edges]++;
                    }
                }
            }

            printf("\n================ Candidate Missing Edges Statistics ================\n");
            printf("Total valid candidates across all query vertices: %u\n", total_candidates_count);
            for (ui i = 0; i <= solver.threshold; ++i) {
                double percent = (total_candidates_count == 0) ? 0.0 :
                    (double)missing_edges_dist[i] / total_candidates_count * 100.0;
                printf("Missing edges = %u: %6u candidates (%6.2f %%)\n", i, missing_edges_dist[i], percent);
            }
            printf("====================================================================\n\n");

            printf("candidates nums and missing edges distribution:\n");
            for (ui u = 0; u < solver.qn; ++u) {
                ui cand_size = solver.candidates[u].size();
                ui deg_u = solver.q_neighbors[u].size();
                ui max_miss_allowed = (deg_u > 0) ? std::min(solver.threshold, deg_u - 1) : solver.threshold;

                printf("u = %u: %6d candidates", u, cand_size);
                if (cand_size > 0) {
                    printf("  [ ");
                    for (ui i = 0; i <= max_miss_allowed; ++i) {
                        double percent = (double)vertex_missing_edges_dist[u][i] / cand_size * 100.0;
                        printf("%u-miss: %5.1f%%", i, percent);
                        if (i < max_miss_allowed) printf(" | ");
                    }
                    printf(" ]");
                }
                printf("\n");
            }
        }
#endif
    };

    bool runCandidateFiltering()
    {
        return CandidateFilter(*this).run();
    }
    // ========================================================================

    struct BranchSelector {
        MatchingSolver &solver;

        explicit BranchSelector(MatchingSolver &solver) : solver(solver) {}

        ui selectInitialRoot() const
        {
            ui root = 0;
            for (ui u = 1; u < solver.qn; ++u) {
                if (solver.candidates[u].size() < solver.candidates[root].size() ||
                    (solver.candidates[u].size() == solver.candidates[root].size() &&
                        solver.q_neighbors[u].size() > solver.q_neighbors[root].size())) {
                    root = u;
                }
            }
            return root;
        }

        void collectLiveAnchors(ui u, vector<ui> &anchors) const
        {
            anchors.clear();
            for (ui nbr : solver.q_neighbors[u]) {
                if (solver.in_Mq[nbr] && !solver.is_excluded[u][nbr]) {
                    anchors.push_back(nbr);
                }
            }
        }

        SearchFrontierState buildSearchFrontierState(const SearchFrontierState *parent_state)
        {
            SearchFrontierState state;
            vector<ui> preferred_frontier;
            collectPreferredDescendantFrontier(parent_state, preferred_frontier);

            const vector<ui> &selection_pool =
                preferred_frontier.empty() ? solver.active_frontier : preferred_frontier;
            if (selection_pool.empty()) {
                return state;
            }

            ui best_u = selectBestFrontierPoint(selection_pool);
            collectComponentStateFromSeed(best_u, state);
            return state;
        }

    private:
        struct FrontierPointScore {
            ui best_anchor_support = std::numeric_limits<ui>::max();
            ui live_candidate_count = std::numeric_limits<ui>::max();
            ui live_anchor_count = 0;
            ui query_degree = 0;
        };

        void advanceFrontierVisitToken()
        {
            if (++solver.frontier_visit_token == 0) {
                std::fill(solver.frontier_visit_mark.begin(), solver.frontier_visit_mark.end(), 0);
                solver.frontier_visit_token = 1;
            }
        }

        ui countLiveCandidates(ui u) const
        {
            ui cnt = 0;
            for (ui v : solver.candidates[u]) {
                if (solver.mapped_g[v] == -1 && !solver.x_cand[u].contains(v)) {
                    cnt++;
                }
            }
            return cnt;
        }

        ui countAnchorSupport(ui u, ui anchor) const
        {
            ui support = 0;
            ui deg = 0;
            const ui *nbrs = solver.data_graph->getVertexNeighbors((ui)solver.mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (solver.mapped_g[v] != -1 || solver.x_cand[u].contains(v)) {
                    continue;
                }
                if (solver.candidates[u].contains(v)) {
                    support++;
                }
            }
            return support;
        }

        FrontierPointScore evaluateFrontierPoint(ui u, vector<ui> &anchors_buf) const
        {
            FrontierPointScore score;
            collectLiveAnchors(u, anchors_buf);

            score.live_anchor_count = (ui)anchors_buf.size();
            score.live_candidate_count = countLiveCandidates(u);
            score.query_degree = (ui)solver.q_neighbors[u].size();

            for (ui anchor : anchors_buf) {
                score.best_anchor_support = std::min(score.best_anchor_support, countAnchorSupport(u, anchor));
            }

            if (score.best_anchor_support == std::numeric_limits<ui>::max()) {
                score.best_anchor_support = score.live_candidate_count;
            }

            return score;
        }

        bool isBetterFrontierPoint(ui lhs_u, const FrontierPointScore &lhs,
            ui rhs_u, const FrontierPointScore &rhs) const
        {
            if (lhs.best_anchor_support != rhs.best_anchor_support) {
                return lhs.best_anchor_support < rhs.best_anchor_support;
            }
            if (lhs.live_candidate_count != rhs.live_candidate_count) {
                return lhs.live_candidate_count < rhs.live_candidate_count;
            }
            if (lhs.live_anchor_count != rhs.live_anchor_count) {
                return lhs.live_anchor_count > rhs.live_anchor_count;
            }
            if (lhs.query_degree != rhs.query_degree) {
                return lhs.query_degree > rhs.query_degree;
            }
            return lhs_u < rhs_u;
        }

        ui selectBestFrontierPoint(const vector<ui> &frontier_points) const
        {
            assert(!frontier_points.empty());

            vector<ui> anchors_buf;
            ui best_u = frontier_points.front();
            FrontierPointScore best_score = evaluateFrontierPoint(best_u, anchors_buf);

            for (ui i = 1; i < frontier_points.size(); ++i) {
                ui u = frontier_points[i];
                FrontierPointScore score = evaluateFrontierPoint(u, anchors_buf);
                if (isBetterFrontierPoint(u, score, best_u, best_score)) {
                    best_u = u;
                    best_score = score;
                }
            }

            return best_u;
        }

        void collectComponentStateFromSeed(ui seed, SearchFrontierState &state)
        {
            state.component_vertices.clear();
            state.component_frontier.clear();

            if (solver.in_Mq[seed]) {
                return;
            }

            advanceFrontierVisitToken();

            queue<ui> q;
            solver.frontier_visit_mark[seed] = solver.frontier_visit_token;
            q.push(seed);

            while (!q.empty()) {
                ui curr = q.front();
                q.pop();

                state.component_vertices.push_back(curr);
                if (solver.frontier_pos[curr] != -1) {
                    state.component_frontier.push_back(curr);
                }

                for (ui nbr : solver.q_neighbors[curr]) {
                    if (!solver.in_Mq[nbr] &&
                        solver.frontier_visit_mark[nbr] != solver.frontier_visit_token) {
                        solver.frontier_visit_mark[nbr] = solver.frontier_visit_token;
                        q.push(nbr);
                    }
                }
            }
        }

        void collectPreferredDescendantFrontier(const SearchFrontierState *parent_state,
            vector<ui> &preferred_frontier)
        {
            preferred_frontier.clear();
            if (parent_state == nullptr || parent_state->component_vertices.empty()) {
                return;
            }

            advanceFrontierVisitToken();

            queue<ui> q;
            for (ui seed : parent_state->component_vertices) {
                if (solver.in_Mq[seed] ||
                    solver.frontier_visit_mark[seed] == solver.frontier_visit_token) {
                    continue;
                }

                solver.frontier_visit_mark[seed] = solver.frontier_visit_token;
                q.push(seed);

                while (!q.empty()) {
                    ui curr = q.front();
                    q.pop();

                    if (solver.frontier_pos[curr] != -1) {
                        preferred_frontier.push_back(curr);
                    }

                    for (ui nbr : solver.q_neighbors[curr]) {
                        if (!solver.in_Mq[nbr] &&
                            solver.frontier_visit_mark[nbr] != solver.frontier_visit_token) {
                            solver.frontier_visit_mark[nbr] = solver.frontier_visit_token;
                            q.push(nbr);
                        }
                    }
                }
            }
        }
    };

    // ========================================================================
    // Lower Bound based Pruning
    // ========================================================================
    struct LowerBoundPruner {
        MatchingSolver &solver;

        explicit LowerBoundPruner(MatchingSolver &solver) : solver(solver) {}

        bool shouldPrune(ui current_miss, const vector<ui> &frontier_vertices)
        {
            ui remaining_budget = solver.threshold - current_miss;
            if (frontier_vertices.empty()) {
                return false;
            }

            BranchSelector branch_selector(solver);
            ui frontier_size = (ui)frontier_vertices.size();
            vector<vector<ui>> anchors_by_frontier_idx(frontier_size);
            ui sum_lb = current_miss;

            for (ui i = 0; i < frontier_size; ++i) {
                ui u = frontier_vertices[i];
                vector<ui> &anchors = anchors_by_frontier_idx[i];
                branch_selector.collectLiveAnchors(u, anchors);

                ui best_anchor_cut = solver.threshold + 1;
                for (ui v : solver.candidates[u]) {
                    if (solver.mapped_g[v] != -1 || solver.x_cand[u].contains(v)) {
                        continue;
                    }

                    ui alpha = computeAnchorCutDelta(v, anchors);
                    best_anchor_cut = std::min(best_anchor_cut, alpha);
                }

                if (best_anchor_cut == solver.threshold + 1) {
                    return true;
                }

                // Search-time lower bounds must stay residual to the current
                // partial match, so the cheap bound only uses the live anchor cut.
                sum_lb += best_anchor_cut;
                if (sum_lb > solver.threshold) {
                    return true;
                }
            }

            ui competition_lb = computeBudgetCompetitionLowerBound(
                frontier_vertices, anchors_by_frontier_idx, remaining_budget);
            return current_miss + competition_lb > solver.threshold;
        }

    private:
        ui computeAnchorCutDelta(ui v, const vector<ui> &anchors) const
        {
            ui delta = 0;
            for (ui a : anchors) {
                assert(solver.mapped_q[a] >= 0);
                if (!solver.data_graph->hasEdge(v, (ui)solver.mapped_q[a])) {
                    delta++;
                }
            }
            return delta;
        }

        bool findBudgetFeasibleAugment(ui left_idx, const vector<vector<ui>> &adj)
        {
            for (ui v : adj[left_idx]) {
                if (solver.lb_seen_right[v] == solver.lb_seen_token) {
                    continue;
                }
                solver.lb_seen_right[v] = solver.lb_seen_token;

                if (solver.lb_match_right[v] < 0 ||
                    findBudgetFeasibleAugment((ui)solver.lb_match_right[v], adj)) {
                    solver.lb_match_right[v] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        ui computeBudgetCompetitionLowerBound(const vector<ui> &frontier_vertices,
            const vector<vector<ui>> &anchors_by_frontier_idx,
            ui remaining_budget)
        {
            ui frontier_size = (ui)frontier_vertices.size();
            if (frontier_size <= 1 || remaining_budget == 0) {
                return 0;
            }

            vector<vector<vector<ui>>> exact_adj(frontier_size, vector<vector<ui>>(remaining_budget));
            vector<ui> right_vertices;

            if (++solver.lb_data_frontier_token == 0) {
                std::fill(solver.lb_data_frontier_mark.begin(), solver.lb_data_frontier_mark.end(), 0);
                solver.lb_data_frontier_token = 1;
            }

            for (ui i = 0; i < frontier_size; ++i) {
                ui u = frontier_vertices[i];
                const vector<ui> &anchors = anchors_by_frontier_idx[i];

                for (ui v : solver.candidates[u]) {
                    if (solver.mapped_g[v] != -1 || solver.x_cand[u].contains(v)) {
                        continue;
                    }

                    ui alpha = computeAnchorCutDelta(v, anchors);
                    if (alpha >= remaining_budget) {
                        continue;
                    }

                    exact_adj[i][alpha].push_back(v);
                    if (solver.lb_data_frontier_mark[v] != solver.lb_data_frontier_token) {
                        solver.lb_data_frontier_mark[v] = solver.lb_data_frontier_token;
                        right_vertices.push_back(v);
                    }
                }
            }

            vector<vector<ui>> cumulative_adj(frontier_size);
            ui extra_lb = 0;

            for (ui k = 0; k < remaining_budget; ++k) {
                for (ui i = 0; i < frontier_size; ++i) {
                    const vector<ui> &delta = exact_adj[i][k];
                    cumulative_adj[i].insert(cumulative_adj[i].end(), delta.begin(), delta.end());
                }

                for (ui v : right_vertices) {
                    solver.lb_match_right[v] = -1;
                }

                ui mu = 0;
                for (ui i = 0; i < frontier_size; ++i) {
                    if (++solver.lb_seen_token == 0) {
                        std::fill(solver.lb_seen_right.begin(), solver.lb_seen_right.end(), 0);
                        solver.lb_seen_token = 1;
                    }

                    if (findBudgetFeasibleAugment(i, cumulative_adj)) {
                        mu++;
                    }
                }

                ui deficit = frontier_size - mu;
                extra_lb += deficit;
                if (extra_lb > remaining_budget) {
                    return extra_lb;
                }
            }

            return extra_lb;
        }
    };
    // ========================================================================

    // =====================================================
    // Procedure DFS(M_part, cost, X, P, u_new)
    //
    // cost:  current cost of partial match M_part
    // u_new: the most recently mapped query vertex (-1 means undefined)
    // =====================================================
    void dfs(ui cost, int u_new, const SearchFrontierState *parent_state)
    {
        assert(matched_count > 0);
        assert(matched_count <= qn);
        assert(cost <= threshold);

        if (matched_count == qn) {
            stats.recursion_calls++;
            results_ptr->push_back(part_M);
            return;
        }

        stats.recursion_calls++;

        vector<ui> reEnabledP;
        if (u_new >= 0) {
            for (ui w : q_neighbors[(ui)u_new]) {
                if (in_P[w]) {
                    reEnabledP.push_back(w);
                    in_P[w] = 0;
                    updateFrontierStatus(w);
                }
            }
        }

        Timer t_frontier;

        if (active_frontier.empty()) {
            // stats.prun_calls++;
            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }

        BranchSelector branch_selector(*this);
        SearchFrontierState current_state = branch_selector.buildSearchFrontierState(parent_state);
        const vector<ui> &U_frontier = current_state.component_frontier;
        stats.frontier_time += t_frontier.elapsed();

        if (U_frontier.empty()) {
            for (ui w : reEnabledP) {
                in_P[w] = 1;
                updateFrontierStatus(w);
            }
            return;
        }

#ifdef LOWER_BOUND
        //         Timer t_lb;
        //         if (LowerBoundPruner(*this).shouldPrune(cost, U_frontier)) {
        //             stats.lb_time += t_lb.elapsed();
        //             stats.prun_calls++;

        //             for (ui w : reEnabledP) {
        //                 in_P[w] = 1;
        //                 updateFrontierStatus(w);
        //             }
        //             return;
        //         }
        //         stats.lb_time += t_lb.elapsed();
#endif

        Timer t_branch;
        long long child_dfs_time = 0;
        vector<ui> local_P;
        vector<pair<ui, ui>> local_X;
        vector<pair<ui, ui>> local_x_cand;

        ui current_cost = cost;

        for (ui u : U_frontier) {
            vector<ui> U_anchor;
            branch_selector.collectLiveAnchors(u, U_anchor);

            bool threshold_exceeded = false;

            // --- matching u 按 anchor 进行细分 ---
            for (ui ua : U_anchor) {
                if (is_excluded[u][ua]) continue;

                // branch 1: matching u by (u, ua) edge
                vector<ui> anchor_v_list;
                ui deg; const ui *nbrs = data_graph->getVertexNeighbors(mapped_q[ua], deg);
                for (ui j = 0; j < deg; ++j) {
                    ui v = nbrs[j];
                    if (candidates[u].contains(v)) anchor_v_list.push_back(v);
                }

                for (ui v : anchor_v_list) {
                    if (mapped_g[v] != -1) continue;
                    if (x_cand[u].contains(v)) continue;

                    ui delta = 0;

                    for (ui other_ua : U_anchor) {
                        if (other_ua == ua) continue;
                        if (is_excluded[u][other_ua]) continue;

                        bool has_edge = data_graph->hasEdge(v, mapped_q[other_ua]);

                        if (!has_edge) delta++;
                    }

                    if (current_cost + delta > threshold) continue;

                    mapped_q[u] = (int)v;
                    mapped_g[v] = (int)u;
                    in_Mq[u] = 1;
                    matched_count++;
                    part_M.push_back({ u, v });

                    onVertexMatchStateChanged(u, true);

                    Timer t_child;
                    dfs(current_cost + delta, (int)u, &current_state);
                    child_dfs_time += t_child.elapsed();

                    onVertexMatchStateChanged(u, false);

                    part_M.pop_back();
                    matched_count--;
                    in_Mq[u] = 0;
                    mapped_g[v] = -1;
                    mapped_q[u] = -1;
                }

                // branch 2 : excluding (u, ua)
                current_cost++;
                is_excluded[u][ua] = 1;
                is_excluded[ua][u] = 1;
                anchor_count[u]--;
                local_X.push_back({ u, ua });

                if (current_cost > threshold) {
                    threshold_exceeded = true;
                    break;
                }

                for (ui v : anchor_v_list) {
                    if (!x_cand[u].contains(v)) {
                        x_cand[u].insert(v);
                        local_x_cand.push_back({ u, v });
                    }
                }
            }

            if (threshold_exceeded) {
                break;
            }

            in_P[u] = 1;
            local_P.push_back(u);
            updateFrontierStatus(u);
        }

        stats.branch_time += t_branch.elapsed() - child_dfs_time;

        // backtracking
        for (ui u : local_P) {
            in_P[u] = 0;
            updateFrontierStatus(u);
        }

        for (auto &e : local_X) {
            is_excluded[e.first][e.second] = 0;
            is_excluded[e.second][e.first] = 0;

            if (in_Mq[e.second]) anchor_count[e.first]++;
            if (in_Mq[e.first]) anchor_count[e.second]++;
            updateFrontierStatus(e.first);
            updateFrontierStatus(e.second);
        }

        for (auto &p : local_x_cand) {
            x_cand[p.first].remove(p.second);
        }

        for (ui w : reEnabledP) {
            in_P[w] = 1;
            updateFrontierStatus(w);
        }
    }
    // ========================================================================
};

// ============================================================
// Top-level function: Approximate_Matching
// ============================================================
void Approximate_Matching(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
}

#endif
