#ifndef MATCHING_ALGORITHMS_CDE_EDGE_IE_H_
#define MATCHING_ALGORITHMS_CDE_EDGE_IE_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
#include <limits>
#ifndef NDEBUG
#include <map>
#endif

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
        initEdgeOrdering();

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
#ifndef NDEBUG
        if (kDebugInitialRoot >= 0) {
            assert((ui)kDebugInitialRoot < qn);
            root = (ui)kDebugInitialRoot;
        }
#endif
        printf("Selected initial root: u=%u with %zu candidates\n", root, (size_t)candidates[root].size());

        for (ui v0 : candidates[root]) {
#ifndef NDEBUG
            recordBranchOrderDebug(root);
#endif
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
        long long frontier_time = 0; // building and ordering U_frontier
        long long frontier_preferred_time = 0;
        long long frontier_select_time = 0;
        long long frontier_component_time = 0;
        long long frontier_sort_time = 0;
        long long frontier_score_time = 0;
        long long frontier_score_live_anchor_time = 0;
        long long frontier_score_live_candidate_time = 0;
        long long frontier_score_query_degree_time = 0;
        long long frontier_score_anchor_support_time = 0;
        long long branch_time = 0;   // candidate enumeration & matching in dfs
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
        size_t result_count = 0;
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        printf("\n--- CDE-Edge-IE Time Analysis ---\n");
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("    - Preferred:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_preferred_time / 1000.0, pct(stats.frontier_preferred_time, stats.frontier_time));
        printf("    - Select Best:   %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_select_time / 1000.0, pct(stats.frontier_select_time, stats.frontier_time));
        printf("    - Component:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_component_time / 1000.0, pct(stats.frontier_component_time, stats.frontier_time));
        printf("    - Sort Hook:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_sort_time / 1000.0, pct(stats.frontier_sort_time, stats.frontier_time));
        printf("    - Score Total:   %.4lf ms (%.2f%% of Select+Sort)\n", stats.frontier_score_time / 1000.0, pct(stats.frontier_score_time, stats.frontier_select_time + stats.frontier_sort_time));
        printf("      - Live Anchors:   %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_live_anchor_time / 1000.0, pct(stats.frontier_score_live_anchor_time, stats.frontier_score_time));
        printf("      - Live Candidates:%.4lf ms (%.2f%% of Score)\n", stats.frontier_score_live_candidate_time / 1000.0, pct(stats.frontier_score_live_candidate_time, stats.frontier_score_time));
        printf("      - Query Degree:   %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_query_degree_time / 1000.0, pct(stats.frontier_score_query_degree_time, stats.frontier_score_time));
        printf("      - Anchor Support: %.4lf ms (%.2f%% of Score)\n", stats.frontier_score_anchor_support_time / 1000.0, pct(stats.frontier_score_anchor_support_time, stats.frontier_score_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
        printf("Results Found:       %zu\n", stats.result_count);
#ifndef NDEBUG
        printBranchOrderDebugStats();
#endif
        printf("-----------------------------------------------------------\n");
        fflush(stdout);
    }

private:
#ifndef NDEBUG
    // Set to a query vertex id for local root experiments; -1 keeps automatic selection.
    enum { kDebugInitialRoot = -1 };
    // Internal debug knob: set to 0 to disable branch-order prefix counting.
    enum { kDebugBranchOrderDepth = 1 };
#endif

    struct OrderedQueryEdge {
        ui u = 0;
        ui v = 0;
        ui weight = 0;
        ui rank = 0;

        bool operator<(const OrderedQueryEdge &other) const
        {
            if (weight != other.weight) return weight < other.weight;
            if (u != other.u) return u < other.u;
            return v < other.v;
        }
    };

    struct ActiveEdge {
        ui u = 0;       // unmatched endpoint
        ui anchor = 0;  // matched endpoint
        int rank = std::numeric_limits<int>::max();
    };

    const Graph *query_graph;
    const Graph *data_graph;
    vector<vector<pair<ui, ui>>> *results_ptr;
    ui threshold;
    ui qn, gn;
    ui label_count;
    ui max_g_deg;
    vector<vector<char>> q_matrix;
    vector<vector<ui>> q_neighbors;
    vector<OrderedQueryEdge> ordered_query_edges;
    vector<vector<int>> q_edge_rank;

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
#ifndef NDEBUG
    map<vector<ui>, unsigned long long> branch_order_prefix_counts;
    map<vector<ui>, unsigned long long> branch_order_prefix_reach_counts;
    map<vector<ui>, unsigned long long> branch_order_prefix_answer_counts;
#endif

    struct FrontierScore {
        explicit FrontierScore(ui frontier_u = 0) : u(frontier_u) {}

        ui u = 0;
        vector<ui> live_anchors;
        ui best_anchor_support = std::numeric_limits<ui>::max();
        ui live_candidate_count = std::numeric_limits<ui>::max();
        ui live_anchor_count = 0;
        ui query_degree = 0;

        bool live_anchors_ready = false;
        bool best_anchor_support_ready = false;
        bool live_candidate_count_ready = false;
        bool query_degree_ready = false;
    };

    // Frontier state for the selected unmatched query component in the current DFS step.
    struct FrontierState {
        // All unmatched query vertices in the selected connected component
        vector<ui> component_vertices;
        // Active frontier vertices within the selected component
        vector<ui> component_frontier;
        // Lazy score cache computed while choosing/sorting frontier vertices in this state.
        vector<FrontierScore> frontier_score_cache;
        vector<char> frontier_score_cached;
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
        ordered_query_edges.clear();
        q_edge_rank.clear();

        anchor_count.assign(qn, 0);
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
#ifndef NDEBUG
        branch_order_prefix_counts.clear();
        branch_order_prefix_reach_counts.clear();
        branch_order_prefix_answer_counts.clear();
#endif
        frontier_visit.assign(qn, 0);
        frontier_token = 1;

        lb_match_right.assign(gn, -1);
        lb_seen_right.assign(gn, 0);
        lb_seen_token = 1;
        lb_data_frontier_mark.assign(gn, 0);
        lb_data_frontier_token = 1;

        stats = TimeStats();
    }

#ifndef NDEBUG
    void recordBranchOrderDebug(ui extension_u)
    {
        if (kDebugBranchOrderDepth == 0) {
            return;
        }

        ui order_size = (ui)part_M.size() + 1;
        ui prefix_limit = std::min((ui)kDebugBranchOrderDepth, order_size);
        vector<ui> prefix;
        prefix.reserve(prefix_limit);

        for (ui i = 0; i < prefix_limit; ++i) {
            ui u = (i < part_M.size()) ? part_M[i].first : extension_u;
            prefix.push_back(u);
            branch_order_prefix_counts[prefix]++;
        }

        if (order_size <= (ui)kDebugBranchOrderDepth) {
            branch_order_prefix_reach_counts[prefix]++;
        }
    }

    void recordBranchOrderAnswerDebug()
    {
        if (kDebugBranchOrderDepth == 0) {
            return;
        }

        ui prefix_limit = std::min((ui)kDebugBranchOrderDepth, (ui)part_M.size());
        vector<ui> prefix;
        prefix.reserve(prefix_limit);

        for (ui i = 0; i < prefix_limit; ++i) {
            prefix.push_back(part_M[i].first);
            branch_order_prefix_answer_counts[prefix]++;
        }
    }

    void printOrderPrefix(const vector<ui> &prefix) const
    {
        printf("[");
        for (ui i = 0; i < prefix.size(); ++i) {
            if (i > 0) printf(",");
            printf("%u", prefix[i]);
        }
        printf("]");
    }

    void printBranchOrderDebugStats() const
    {
        if (kDebugBranchOrderDepth == 0) {
            return;
        }

        printf("\n--- Branch Order Prefix Debug ---\n");
        printf("Debug Depth:         %u\n", (ui)kDebugBranchOrderDepth);
        printf("Counted Branches:    matching branches + exclusion branches\n");
        printf("Reached:             exact branches that arrive at this prefix depth\n");
        printf("Tail Branches:       subtree branches after the prefix is reached\n");
        printf("Counted Answers:     complete matches grouped by their extension-order prefix\n");
        printf("Depth Rule:          the initial root expansion is depth 1\n");
        printf("Prefix Rule:         each branch is added to every recorded prefix of its extension order\n");

        if (branch_order_prefix_counts.empty()) {
            printf("No branch-order data recorded.\n");
            return;
        }

        for (ui depth = 1; depth <= kDebugBranchOrderDepth; ++depth) {
            unsigned long long reach_depth_total = 0;
            unsigned long long branch_depth_total = 0;
            unsigned long long tail_branch_depth_total = 0;
            unsigned long long answer_depth_total = 0;
            bool printed_header = false;

            for (map<vector<ui>, unsigned long long>::const_iterator it = branch_order_prefix_counts.begin();
                it != branch_order_prefix_counts.end(); ++it) {
                if (it->first.size() != depth) {
                    continue;
                }

                if (!printed_header) {
                    printf("Depth %u:\n", depth);
                    printed_header = true;
                }

                printf("  ");
                printOrderPrefix(it->first);
                unsigned long long reach_count = 0;
                unsigned long long answer_count = 0;
                map<vector<ui>, unsigned long long>::const_iterator reach_it =
                    branch_order_prefix_reach_counts.find(it->first);
                if (reach_it != branch_order_prefix_reach_counts.end()) {
                    reach_count = reach_it->second;
                }

                map<vector<ui>, unsigned long long>::const_iterator answer_it =
                    branch_order_prefix_answer_counts.find(it->first);
                if (answer_it != branch_order_prefix_answer_counts.end()) {
                    answer_count = answer_it->second;
                }

                unsigned long long tail_branch_count = it->second >= reach_count ? it->second - reach_count : 0;
                printf(" -> reached=%llu, tail_branches=%llu, subtree_branches=%llu, answers=%llu\n",
                    reach_count, tail_branch_count, it->second, answer_count);
                reach_depth_total += reach_count;
                branch_depth_total += it->second;
                tail_branch_depth_total += tail_branch_count;
                answer_depth_total += answer_count;
            }

            if (printed_header) {
                printf("  Total -> reached=%llu, tail_branches=%llu, subtree_branches=%llu, answers=%llu\n",
                    reach_depth_total, tail_branch_depth_total, branch_depth_total, answer_depth_total);
            }
        }
    }
#endif

    void initEdgeOrdering()
    {
        ordered_query_edges.clear();
        q_edge_rank.assign(qn, vector<int>(qn, -1));

        for (ui u = 0; u < qn; ++u) {
            for (ui v : q_neighbors[u]) {
                if (u >= v) continue;

                OrderedQueryEdge edge;
                edge.u = u;
                edge.v = v;
                edge.weight = (ui)candidates[u].size() + (ui)candidates[v].size()
                    + (ui)q_neighbors[u].size() + (ui)q_neighbors[v].size();
                ordered_query_edges.push_back(edge);
            }
        }

        sort(ordered_query_edges.begin(), ordered_query_edges.end());
        for (ui i = 0; i < ordered_query_edges.size(); ++i) {
            OrderedQueryEdge &edge = ordered_query_edges[i];
            edge.rank = i;
            q_edge_rank[edge.u][edge.v] = (int)i;
            q_edge_rank[edge.v][edge.u] = (int)i;
        }
    }

    bool findBestActiveEdge(const vector<ui> &component_frontier, ActiveEdge &best) const
    {
        bool found = false;
        for (ui u : component_frontier) {
            if (u >= qn || in_Mq[u] || frontier_pos[u] == -1) {
                continue;
            }

            for (ui anchor : q_neighbors[u]) {
                if (!in_Mq[anchor] || is_excluded[u][anchor]) {
                    continue;
                }

                int rank = q_edge_rank[u][anchor];
                if (rank < 0) {
                    continue;
                }

                if (!found || rank < best.rank ||
                    (rank == best.rank && (u < best.u || (u == best.u && anchor < best.anchor)))) {
                    best.u = u;
                    best.anchor = anchor;
                    best.rank = rank;
                    found = true;
                }
            }
        }
        return found;
    }

    void collectEdgeCandidates(ui u, ui anchor, vector<ui> &cand_v_list) const
    {
        cand_v_list.clear();
        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
        for (ui i = 0; i < deg; ++i) {
            ui v = nbrs[i];
            if (!candidates[u].contains(v)) continue;
            if (mapped_g[v] != -1) continue;
            if (x_cand[u].contains(v)) continue;
            cand_v_list.push_back(v);
        }
    }

    ui computeLiveAnchorDelta(ui u, ui v, ui selected_anchor) const
    {
        ui delta = 0;
        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (!in_Mq[anchor] || is_excluded[u][anchor]) continue;
            if (!data_graph->hasEdge(v, (ui)mapped_q[anchor])) {
                delta++;
            }
        }
        return delta;
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
#ifdef ENABLE_ADVANCED_FILTERING
            if (!filterSpoke()) return false;
            if (!filterOneHop()) return false;
#endif

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
            fflush(stdout);
        }
    };

    bool runCandidateFiltering()
    {
        return CandidateFilter(*this).run();
    }
    // ========================================================================

    struct BranchSelector {
    private:
        MatchingSolver &solver;
        // Valid while this BranchSelector builds one unchanged DFS state.
        mutable FrontierState *score_state = nullptr;
        mutable const FrontierState *parent_score_state = nullptr;
        mutable int score_newly_matched_u = -1;

    public:
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

        void sortAnchorsBySupport(ui u, vector<ui> &anchors) const
        {
            if (anchors.size() <= 1) {
                return;
            }

            vector<pair<ui, ui>> scored_anchors;
            scored_anchors.reserve(anchors.size());

            for (ui anchor : anchors) {
                scored_anchors.push_back({ countAnchorSupport(u, anchor), anchor });
            }

            sort(scored_anchors.begin(), scored_anchors.end());

            for (ui i = 0; i < scored_anchors.size(); ++i) {
                anchors[i] = scored_anchors[i].second;
            }
        }

        FrontierState buildFrontierState(const FrontierState *parent_state, int newly_matched_u)
        {
            FrontierState state;
            score_state = &state;
            parent_score_state = parent_state;
            score_newly_matched_u = newly_matched_u;

            assert(!solver.active_frontier.empty());

            Timer t_select;
            ui best_u = selectBestFrontierVertex(solver.active_frontier);
            solver.stats.frontier_select_time += t_select.elapsed();

            Timer t_component;
            collectComponent(best_u, state);
            solver.stats.frontier_component_time += t_component.elapsed();

#ifdef ENABLE_FRONTIER_ORDERING
            Timer t_sort;
            sortFrontier(state.component_frontier);
            solver.stats.frontier_sort_time += t_sort.elapsed();
#endif
            return state;
        }

    private:
        void sortFrontier(vector<ui> &frontier) const
        {
            if (frontier.size() <= 1) {
                return;
            }

            sort(frontier.begin(), frontier.end(), [&](ui lhs_u, ui rhs_u) {
                FrontierScore &lhs = scoreFor(lhs_u);
                FrontierScore &rhs = scoreFor(rhs_u);
                return isBetterFrontier(lhs, rhs);
                });
        }

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

        // Frontier scoring intentionally ignores injectivity; this keeps scores
        // stable across most child states after one data vertex becomes matched.
        ui countScoreCandidates(ui u) const
        {
            ui cnt = 0;
            for (ui v : solver.candidates[u]) {
                if (!solver.x_cand[u].contains(v)) {
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

        ui countScoreAnchorSupport(ui u, ui anchor) const
        {
            ui support = 0;
            ui deg = 0;
            const ui *nbrs = solver.data_graph->getVertexNeighbors((ui)solver.mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (solver.x_cand[u].contains(v)) {
                    continue;
                }
                if (solver.candidates[u].contains(v)) {
                    support++;
                }
            }
            return support;
        }

        void ensureScoreStorage(FrontierState &state) const
        {
            if (state.frontier_score_cache.empty()) {
                state.frontier_score_cache.resize(solver.qn);
                state.frontier_score_cached.assign(solver.qn, 0);
            }
        }

        bool canReuseParentScore(ui u) const
        {
            if (parent_score_state == nullptr) {
                return false;
            }
            if (u >= parent_score_state->frontier_score_cached.size() ||
                !parent_score_state->frontier_score_cached[u]) {
                return false;
            }
            if (score_newly_matched_u >= 0) {
                ui newly_matched_u = (ui)score_newly_matched_u;
                if (u == newly_matched_u || solver.q_matrix[u][newly_matched_u]) {
                    return false;
                }
            }
            return true;
        }

        // Returns the cached lazy frontier score for u in the current DFS state.
        FrontierScore &scoreFor(ui u) const
        {
            assert(u < solver.qn);
            assert(score_state != nullptr);
            ensureScoreStorage(*score_state);
            if (!score_state->frontier_score_cached[u]) {
                if (canReuseParentScore(u)) {
                    score_state->frontier_score_cache[u] = parent_score_state->frontier_score_cache[u];
                }
                else {
                    score_state->frontier_score_cache[u] = FrontierScore(u);
                }
                score_state->frontier_score_cached[u] = 1;
            }
            return score_state->frontier_score_cache[u];
        }

        void collectAnchors(FrontierScore &score) const
        {
            if (score.live_anchors_ready) {
                return;
            }

            Timer t;
            collectLiveAnchors(score.u, score.live_anchors);
            score.live_anchor_count = (ui)score.live_anchors.size();
            score.live_anchors_ready = true;

            long long elapsed = t.elapsed();
            solver.stats.frontier_score_live_anchor_time += elapsed;
            solver.stats.frontier_score_time += elapsed;
        }

        ui cachedBestAnchorSupportCount(FrontierScore &score) const
        {
            if (!score.best_anchor_support_ready) {
                collectAnchors(score);

                Timer t;
                for (ui anchor : score.live_anchors) {
                    score.best_anchor_support = std::min(score.best_anchor_support, countScoreAnchorSupport(score.u, anchor));
                }

                long long elapsed = t.elapsed();
                solver.stats.frontier_score_anchor_support_time += elapsed;
                solver.stats.frontier_score_time += elapsed;

                if (score.best_anchor_support == std::numeric_limits<ui>::max()) {
                    score.best_anchor_support = cachedLiveCandidateCount(score);
                }
                score.best_anchor_support_ready = true;
            }

            return score.best_anchor_support;
        }

        ui cachedLiveCandidateCount(FrontierScore &score) const
        {
            if (!score.live_candidate_count_ready) {
                Timer t;
                score.live_candidate_count = countScoreCandidates(score.u);
                score.live_candidate_count_ready = true;

                long long elapsed = t.elapsed();
                solver.stats.frontier_score_live_candidate_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
            }

            return score.live_candidate_count;
        }

        ui cachedLiveAnchorCount(FrontierScore &score) const
        {
            collectAnchors(score);
            return score.live_anchor_count;
        }

        ui cachedQueryDegreeCount(FrontierScore &score) const
        {
            if (!score.query_degree_ready) {
                Timer t;
                score.query_degree = (ui)solver.q_neighbors[score.u].size();
                score.query_degree_ready = true;

                long long elapsed = t.elapsed();
                solver.stats.frontier_score_query_degree_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
            }

            return score.query_degree;
        }

        // Smaller anchor support, fewer live candidates, more live anchors, higher query degree, smaller vertex id
        bool isBetterFrontier(FrontierScore &lhs, FrontierScore &rhs) const
        {
            ui lhs_best_anchor_support = cachedBestAnchorSupportCount(lhs);
            ui rhs_best_anchor_support = cachedBestAnchorSupportCount(rhs);
            if (lhs_best_anchor_support != rhs_best_anchor_support) {
                return lhs_best_anchor_support < rhs_best_anchor_support;
            }

            // ui lhs_live_candidate_count = cachedLiveCandidateCount(lhs);
            // ui rhs_live_candidate_count = cachedLiveCandidateCount(rhs);
            // if (lhs_live_candidate_count != rhs_live_candidate_count) {
            //     return lhs_live_candidate_count < rhs_live_candidate_count;
            // }

            ui lhs_live_anchor_count = cachedLiveAnchorCount(lhs);
            ui rhs_live_anchor_count = cachedLiveAnchorCount(rhs);
            if (lhs_live_anchor_count != rhs_live_anchor_count) {
                return lhs_live_anchor_count > rhs_live_anchor_count;
            }

            ui lhs_query_degree = cachedQueryDegreeCount(lhs);
            ui rhs_query_degree = cachedQueryDegreeCount(rhs);
            if (lhs_query_degree != rhs_query_degree) {
                return lhs_query_degree > rhs_query_degree;
            }

            return lhs.u < rhs.u;
        }

        ui selectBestFrontierVertex(const vector<ui> &frontier_candidates) const
        {
            assert(!frontier_candidates.empty());

            FrontierScore *best_score = &scoreFor(frontier_candidates.front());

            for (ui i = 1; i < frontier_candidates.size(); ++i) {
                FrontierScore &score = scoreFor(frontier_candidates[i]);
                if (isBetterFrontier(score, *best_score)) {
                    best_score = &score;
                }
            }

            return best_score->u;
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
#ifndef NDEBUG
            recordBranchOrderAnswerDebug();
#endif
#ifdef NDEBUG
            stats.result_count++;
#else
            results_ptr->push_back(part_M);
            stats.result_count++;
#endif
            return;
        }

        if (active_frontier.empty()) return;

        stats.recursion_calls++;

        (void)parent_state;
        (void)u_new;

        Timer t_frontier;
        FrontierState current_state;
        current_state.component_frontier = active_frontier;
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
        vector<pair<ui, ui>> local_X;
        vector<pair<ui, ui>> local_x_cand;
        vector<ui> cand_v_list;

        ui current_cost = cost;

        while (current_cost <= threshold) {
            ActiveEdge edge;
            if (!findBestActiveEdge(U_frontier, edge)) {
                break;
            }

            ui u = edge.u;
            ui ua = edge.anchor;
            assert(!in_Mq[u]);
            assert(in_Mq[ua]);
            assert(!is_excluded[u][ua]);

            collectEdgeCandidates(u, ua, cand_v_list);

            // Include branch: use this frontier-anchor edge to add exactly one
            // new query vertex in the recursive child state.
            for (ui v : cand_v_list) {
                assert(mapped_g[v] == -1);
                assert(!x_cand[u].contains(v));

                ui delta = computeLiveAnchorDelta(u, v, ua);
                if (current_cost + delta > threshold) continue;

#ifndef NDEBUG
                recordBranchOrderDebug(u);
#endif
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

            // Exclude branch: keep this decision in the current frame, then
            // continue to the next globally ordered active edge without a
            // recursive call. Every recursive child above still adds a vertex.
            current_cost++;
            is_excluded[u][ua] = 1;
            is_excluded[ua][u] = 1;
            anchor_count[u]--;
            local_X.push_back({ u, ua });
            updateFrontierStatus(u);

            if (current_cost > threshold) {
                break;
            }

#ifndef NDEBUG
            recordBranchOrderDebug(u);
#endif
            for (ui v : cand_v_list) {
                if (!x_cand[u].contains(v)) {
                    x_cand[u].insert(v);
                    local_x_cand.push_back({ u, v });
                }
            }
        }

        stats.branch_time += t_branch.elapsed() - child_dfs_time;

        for (auto it = local_X.rbegin(); it != local_X.rend(); ++it) {
            const pair<ui, ui> &e = *it;
            is_excluded[e.first][e.second] = 0;
            is_excluded[e.second][e.first] = 0;

            assert(!in_Mq[e.first]);
            assert(in_Mq[e.second]);

            anchor_count[e.first]++;
            updateFrontierStatus(e.first);
        }

        for (auto it = local_x_cand.rbegin(); it != local_x_cand.rend(); ++it) {
            x_cand[it->first].remove(it->second);
        }
    }
    // ========================================================================
};

// ============================================================
// Top-level function: Approximate_CDE_EdgeIE
// ============================================================
void Approximate_CDE_EdgeIE(const Graph *query_graph, const Graph *data_graph, vector<vector<pair<ui, ui> > > &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm_ged::set_reported_result_count(solver.stats.result_count);
}

#endif
