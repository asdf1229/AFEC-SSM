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
        max_g_deg = data_graph->getMaxDegree();
        label_count = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());

        if (qn == 0 || gn == 0) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        resetState();

        q_matrix.assign(qn, vector<char>(qn, 0));
        q_neighbors.assign(qn, vector<ui>());

        for (ui u = 0; u < qn; ++u) {
            ui deg = 0; const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) {
                ui u1 = nbrs[i];
                q_matrix[u][u1] = 1;
                q_neighbors[u].push_back(u1);
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
            part_M.push_back({ root, v0 });

            updateFrontier(root, true);

            dfs(0, (int)root, nullptr);

            updateFrontier(root, false);

            mapped_q[root] = -1;
            mapped_g[v0] = -1;
            in_Mq[root] = 0;
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
        long long lb_time = 0;       // computeLowerBound
        long long frontier_time = 0; // building U_frontier
        long long branch_time = 0;   // candidate enumeration & matching in dfs
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
    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui>  Lq_degrees, Lg_degrees;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<char> in_Mq;
    vector<vector<char>> is_excluded;
    vector<MyBitset> x_cand;
    vector<pair<ui, ui>> part_M;

    vector<ui> anchor_count;
    vector<int> frontier_pos;
    vector<ui> active_frontier;

    // Frontier state for the selected unmatched query component in the current DFS step.
    struct FrontierState {
        // All unmatched query vertices in the selected connected component
        vector<ui> component_vertices;
        // Active frontier vertices within the selected component
        vector<ui> component_frontier;
    };
    vector<ui> frontier_visit;
    ui frontier_token;

    vector<int> lb_match_right;
    vector<ui>  lb_seen_right;
    ui          lb_seen_token;
    vector<ui>  lb_data_frontier_mark;
    ui          lb_data_frontier_token;

    // Update active_frontier
    inline void updateFrontierStatus(ui u)
    {
        // Should u be in frontier?
        bool should_be = (!in_Mq[u] && anchor_count[u] > 0);
        // Is u already in frontier?
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

    // Update anchor counts and maintain active_frontier when vertex u becomes matched/unmatched
    // matched = true: add, false: remove
    void updateFrontier(ui u, bool matched)
    {
        if (matched) {
            updateFrontierStatus(u);
        }

        for (ui u1 : q_neighbors[u]) {
            if (!is_excluded[u1][u]) {
                if (matched) {
                    anchor_count[u1]++;
                }
                else {
                    anchor_count[u1]--;
                }
            }
            updateFrontierStatus(u1);
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
        frontier_visit.assign(qn, 0);
        frontier_token = 1;

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
        vector<vector<ui>>  spoke_matrix;
        vector<int>         spoke_match;
        vector<bool>        spoke_vis;
        MyBitset            onehop_vis;

        explicit CandidateFilter(MatchingSolver &solver)
            : solver(solver),
            spoke_matrix(solver.qn, vector<ui>()),
            spoke_match(solver.max_g_deg, -1),
            spoke_vis(solver.max_g_deg, false),
            onehop_vis((int)solver.gn)
        {
        }

        bool run()
        {
            if (!filterNLF()) return false;
            if (!filterSpoke()) return false;
            if (!filterOneHop()) return false;

#ifndef NDEBUG
            printCandStats();
#endif
            return true;
        }

    private:
        ui computeNLF(ui u, ui v) const
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

        bool filterNLF()
        {
            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = solver.query_graph->getVertexLabel(u);
                for (ui v = 0; v < solver.gn; ++v) {
                    if (lu != solver.data_graph->getVertexLabel(v)) continue;
                    if (solver.Lq_degrees[u] > solver.Lg_degrees[v] + solver.threshold) continue;
                    if (computeNLF(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        bool dfsMatchSpoke(ui left_idx, const vector<vector<ui>> &adj)
        {
            for (ui right_idx : adj[left_idx]) {
                if (spoke_vis[right_idx]) continue;
                spoke_vis[right_idx] = true;
                if (spoke_match[right_idx] < 0 ||
                    dfsMatchSpoke((ui)spoke_match[right_idx], adj)) {
                    spoke_match[right_idx] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        ui computeMaxMatchSpoke(const vector<vector<ui>> &adj, ui left_size, ui right_size)
        {
            ui mu = 0;
            std::fill(spoke_match.begin(), spoke_match.begin() + right_size, -1);
            for (ui i = 0; i < left_size; ++i) {
                std::fill(spoke_vis.begin(), spoke_vis.begin() + right_size, false);
                if (dfsMatchSpoke(i, adj)) mu++;
            }
            return mu;
        }

        // If u is matched to v, compute the minimum number of spoke edges
        // from u to its neighbors that cannot be supported by neighbors of v.
        ui computeLBSpoke(ui u, ui v)
        {
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            ui deg_u = (ui)u_neighbors.size();
            ui deg_v = 0;
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            for (ui i = 0; i < deg_u; ++i) {
                spoke_matrix[i].clear();
                ui u1 = u_neighbors[i];
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        spoke_matrix[i].push_back(j);
                    }
                }
            }

            ui match_size = computeMaxMatchSpoke(spoke_matrix, deg_u, deg_v);
            return deg_u - match_size;
        }

        bool filterSpoke()
        {
            queue<ui> q;
            vector<char> in_q(solver.qn, 1);
            for (ui u = 0; u < solver.qn; ++u) q.push(u);

            while (!q.empty()) {
                ui u = q.front(); q.pop();
                in_q[u] = 0;

                vector<ui> to_remove;
                for (ui v : solver.candidates[u]) {
                    if (computeLBSpoke(u, v) > solver.threshold) to_remove.push_back(v);
                }

                if (to_remove.empty()) continue;
                for (ui v : to_remove) solver.candidates[u].remove(v);
                if (solver.candidates[u].empty()) return false;

                for (ui nbr_u : solver.q_neighbors[u]) {
                    if (!in_q[nbr_u]) {
                        q.push(nbr_u);
                        in_q[nbr_u] = 1;
                    }
                }
            }
            return true;
        }

        ui countInnerEdges(ui u) const
        {
            ui count = 0;
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            for (ui i = 0; i < u_neighbors.size(); ++i) {
                for (ui j = i + 1; j < u_neighbors.size(); ++j) {
                    if (solver.q_matrix[u_neighbors[i]][u_neighbors[j]]) count++;
                }
            }
            return count;
        }

        // Build the DFS order for one-hop matching.
        // Fewer candidates first.
        // More inner edges first if tied.
        vector<ui> buildOneHopOrder(const vector<ui> &u_neighbors, const vector<vector<ui>> &cand) const
        {
            ui deg_u = (ui)u_neighbors.size();
            vector<ui> ord(deg_u);
            iota(ord.begin(), ord.end(), 0);

            sort(ord.begin(), ord.end(), [&](ui a, ui b) {
                if (cand[a].size() != cand[b].size()) {
                    return cand[a].size() < cand[b].size();
                }
                ui deg_a = 0, deg_b = 0;
                for (ui z : u_neighbors) {
                    if (solver.q_matrix[u_neighbors[a]][z]) deg_a++;
                    if (solver.q_matrix[u_neighbors[b]][z]) deg_b++;
                }
                return deg_a > deg_b;
                });
            return ord;
        }

        ui computeRemainLBOneHop(ui pos, const vector<ui> &ord,
            const vector<vector<ui>> &cand) const
        {
            ui rem = 0;
            for (ui p = pos; p < ord.size(); ++p) {
                ui u1 = ord[p];
                bool has_free = false;
                for (ui v1 : cand[u1]) {
                    if (!onehop_vis.contains(v1)) {
                        has_free = true;
                        break;
                    }
                }
                if (!has_free) rem++;
            }
            return rem;
        }

        // DFS for one-hop matching
        // Branch 1: skip u1, adding one missing spoke edge
        // Branch 2: match u1 to v1, adding newly missing inner edges
        bool dfsOneHop(ui pos, vector<int> &state, ui cost,
            const vector<ui> &ord, const vector<ui> &u_neighbors, const vector<vector<ui>> &cand)
        {
            if (cost > solver.threshold) return false;
            if (pos == ord.size()) return true;

            ui rem_lb = computeRemainLBOneHop(pos, ord, cand);
            if (cost + rem_lb > solver.threshold) return false;

            ui i = ord[pos];
            ui u1 = u_neighbors[i];

            // branch 1: skip u1
            state[i] = -1;
            if (dfsOneHop(pos + 1, state, cost + 1, ord, u_neighbors, cand)) {
                state[i] = -2;
                return true;
            }

            // branch 2: match u1 to v1
            for (ui v1 : cand[i]) {
                if (onehop_vis.contains(v1)) continue;

                ui delta_inner = 0;
                for (ui j = 0; j < u_neighbors.size(); ++j) {
                    if (state[j] < 0) continue;
                    ui u2 = u_neighbors[j];
                    if (!solver.q_matrix[u1][u2]) continue;

                    ui v2 = (ui)state[j];
                    if (!solver.data_graph->hasEdge(v1, v2)) delta_inner++;
                }

                if (cost + delta_inner > solver.threshold) continue;

                onehop_vis.insert(v1);
                state[i] = (int)v1;

                if (dfsOneHop(pos + 1, state, cost + delta_inner, ord, u_neighbors, cand)) {
                    onehop_vis.remove(v1);
                    state[i] = -2;
                    return true;
                }

                onehop_vis.remove(v1);
                state[i] = -2;
            }

            state[i] = -2;
            return false;
        }

        bool checkOneHop(ui u, ui v)
        {
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            ui deg_u = (ui)u_neighbors.size();
            ui deg_v = 0;
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);
            vector<vector<ui>> cand(deg_u);

            for (ui i = 0; i < deg_u; ++i) {
                ui u1 = u_neighbors[i];
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        cand[i].push_back(v1);
                    }
                }
            }

            // matching order
            vector<ui> ord = buildOneHopOrder(u_neighbors, cand);

            // the current DFS state of u_neighbors[i]:
            // -2: unprocessed
            // -1: skipped / unmatched
            // >=0: matched to the corresponding data vertex
            vector<int> state(deg_u, -2);

            return dfsOneHop(0, state, 0, ord, u_neighbors, cand);
        }

        bool filterOneHop()
        {
            vector<pair<ui, ui>> to_remove;

            for (ui u = 0; u < solver.qn; ++u) {
                if (solver.q_neighbors[u].size() <= 1) continue;
                if (countInnerEdges(u) == 0) continue;

                for (ui v : solver.candidates[u]) {
                    if (!checkOneHop(u, v)) {
                        to_remove.push_back({ u, v });
                    }
                }
            }

            for (auto p : to_remove) {
                solver.candidates[p.first].remove(p.second);
            }

            for (ui u = 0; u < solver.qn; ++u) {
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        void printCandStats()
        {
            vector<ui> missing_edges_dist(solver.threshold + 1, 0);
            vector<vector<ui>> vertex_missing_edges_dist(solver.qn, vector<ui>(solver.threshold + 1, 0));
            ui total_candidates_count = 0;

            for (ui u = 0; u < solver.qn; ++u) {
                total_candidates_count += solver.candidates[u].size();

                for (ui v : solver.candidates[u]) {
                    ui min_missing_edges = computeLBSpoke(u, v);
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
            for (ui u1 : solver.q_neighbors[u]) {
                if (solver.in_Mq[u1] && !solver.is_excluded[u][u1]) {
                    anchors.push_back(u1);
                }
            }
        }

        FrontierState buildFrontierState(const FrontierState *parent_state)
        {
            FrontierState state;
            vector<ui> preferred_frontier;
            collectPreferredFrontier(parent_state, preferred_frontier);

            const vector<ui> &frontier_candidates = preferred_frontier.empty() ? solver.active_frontier : preferred_frontier;
            assert(!frontier_candidates.empty());

            ui best_u = selectBestFrontierVertex(frontier_candidates);
            collectComponent(best_u, state);
            return state;
        }

    private:
        struct FrontierScore {
            ui best_anchor_support = std::numeric_limits<ui>::max();
            ui live_candidate_count = std::numeric_limits<ui>::max();
            ui live_anchor_count = 0;
            ui query_degree = 0;
        };

        // Move to the next BFS visit token
        void nextToken()
        {
            if (++solver.frontier_token == 0) {
                std::fill(solver.frontier_visit.begin(), solver.frontier_visit.end(), 0);
                solver.frontier_token = 1;
            }
        }

        // Count usable candidate data vertices for query vertex u
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

        // Count usable candidates for u among the data neighbors of a matched anchor
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

        FrontierScore scoreFrontier(ui u, vector<ui> &anchors_buf) const
        {
            FrontierScore score;
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

        // Smaller anchor support, fewer live candidates, more live anchors, higher query degree, smaller vertex id
        bool isBetterFrontier(ui lhs_u, const FrontierScore &lhs, ui rhs_u, const FrontierScore &rhs) const
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

        ui selectBestFrontierVertex(const vector<ui> &frontier_candidates) const
        {
            assert(!frontier_candidates.empty());

            vector<ui> anchors_buf;
            ui best_u = frontier_candidates.front();
            FrontierScore best_score = scoreFrontier(best_u, anchors_buf);

            for (ui i = 1; i < frontier_candidates.size(); ++i) {
                ui u = frontier_candidates[i];
                FrontierScore score = scoreFrontier(u, anchors_buf);
                if (isBetterFrontier(u, score, best_u, best_score)) {
                    best_u = u;
                    best_score = score;
                }
            }

            return best_u;
        }

        void collectComponent(ui best_u, FrontierState &state)
        {
            state.component_vertices.clear();
            state.component_frontier.clear();

            if (solver.in_Mq[best_u]) {
                return;
            }

            nextToken();

            queue<ui> q;
            solver.frontier_visit[best_u] = solver.frontier_token;
            q.push(best_u);

            while (!q.empty()) {
                ui curr = q.front();
                q.pop();

                state.component_vertices.push_back(curr);
                if (solver.frontier_pos[curr] != -1) {
                    state.component_frontier.push_back(curr);
                }

                for (ui nbr : solver.q_neighbors[curr]) {
                    if (!solver.in_Mq[nbr] && solver.frontier_visit[nbr] != solver.frontier_token) {
                        solver.frontier_visit[nbr] = solver.frontier_token;
                        q.push(nbr);
                    }
                }
            }
        }

        // Collect active frontier vertices from the parent's remaining unmatched component
        void collectPreferredFrontier(const FrontierState *parent_state, vector<ui> &preferred_frontier)
        {
            preferred_frontier.clear();
            if (parent_state == nullptr || parent_state->component_vertices.empty()) return;

            nextToken();

            queue<ui> q;
            for (ui u : parent_state->component_vertices) {
                if (solver.in_Mq[u] || solver.frontier_visit[u] == solver.frontier_token) {
                    continue;
                }

                solver.frontier_visit[u] = solver.frontier_token;
                q.push(u);

                while (!q.empty()) {
                    ui curr = q.front();
                    q.pop();

                    if (solver.frontier_pos[curr] != -1) {
                        preferred_frontier.push_back(curr);
                    }

                    for (ui nbr : solver.q_neighbors[curr]) {
                        if (!solver.in_Mq[nbr] && solver.frontier_visit[nbr] != solver.frontier_token) {
                            solver.frontier_visit[nbr] = solver.frontier_token;
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
    // Procedure DFS(M_part, cost, X, u_new)
    //
    // cost:  current cost of partial match M_part
    // u_new: the most recently mapped query vertex (-1 means undefined)
    // X:     the set of excluded edges (u, ua)
    // =====================================================
    void dfs(ui cost, int u_new, const FrontierState *parent_state)
    {
        assert(part_M.size() <= qn);
        assert(cost <= threshold);

        if (part_M.size() == qn) {
            stats.recursion_calls++;
            results_ptr->push_back(part_M);
            return;
        }

        if (active_frontier.empty()) return;

        stats.recursion_calls++;

        Timer t_frontier;
        BranchSelector branch_selector(*this);
        FrontierState current_state = branch_selector.buildFrontierState(parent_state);
        const vector<ui> &U_frontier = current_state.component_frontier;
        assert(!U_frontier.empty());
        stats.frontier_time += t_frontier.elapsed();

#ifdef LOWER_BOUND
        //         Timer t_lb;
        //         if (LowerBoundPruner(*this).shouldPrune(cost, U_frontier)) {
        //             stats.lb_time += t_lb.elapsed();
        //             stats.prun_calls++;
        //             return;
        //         }
        //         stats.lb_time += t_lb.elapsed();
#endif

        Timer t_branch;
        long long child_dfs_time = 0;
        vector<pair<ui, ui>> local_X;       // Records changes to is_excluded
        vector<pair<ui, ui>> local_x_cand;  // Records changes to x_cand

        ui current_cost = cost;

        for (ui u : U_frontier) {
            vector<ui> U_anchor;
            branch_selector.collectLiveAnchors(u, U_anchor);

            bool threshold_exceeded = false;

            // --- Try matching u from each anchor's neighbors ---
            for (ui ua : U_anchor) {
                assert(!is_excluded[u][ua]);

                // branch 1: matching u by (u, ua) edge
                vector<ui> cand_v_list;
                ui deg; const ui *nbrs = data_graph->getVertexNeighbors(mapped_q[ua], deg);
                for (ui j = 0; j < deg; ++j) {
                    ui v = nbrs[j];
                    if (candidates[u].contains(v)) {
                        if (mapped_g[v] != -1) continue;
                        if (x_cand[u].contains(v)) continue;
                        cand_v_list.push_back(v);
                    }
                }

                for (ui v : cand_v_list) {
                    assert(mapped_g[v] == -1);
                    assert(!x_cand[u].contains(v));

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
                    part_M.push_back({ u, v });

                    updateFrontier(u, true);

                    Timer t_child;
                    dfs(current_cost + delta, (int)u, &current_state);
                    child_dfs_time += t_child.elapsed();

                    updateFrontier(u, false);

                    part_M.pop_back();
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

                for (ui v : cand_v_list) {
                    if (!x_cand[u].contains(v)) {
                        x_cand[u].insert(v);
                        local_x_cand.push_back({ u, v });
                    }
                }
            }

            if (threshold_exceeded) {
                break;
            }

            updateFrontierStatus(u);
        }

        stats.branch_time += t_branch.elapsed() - child_dfs_time;

        // backtracking
        for (auto &e : local_X) {
            is_excluded[e.first][e.second] = 0;
            is_excluded[e.second][e.first] = 0;

            assert(!in_Mq[e.first]);
            assert(in_Mq[e.second]);

            anchor_count[e.first]++;
            updateFrontierStatus(e.first);
        }

        for (auto &p : local_x_cand) {
            x_cand[p.first].remove(p.second);
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
