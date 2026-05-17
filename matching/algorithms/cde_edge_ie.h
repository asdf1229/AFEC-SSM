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

#ifndef CDE_LB_MWPM_MAX_ROWS
#define CDE_LB_MWPM_MAX_ROWS 0
#endif
#ifndef CDE_LB_MWPM_MAX_COLS
#define CDE_LB_MWPM_MAX_COLS 0
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

#ifdef ENABLE_BRIDGE_FILTERING
        identifyQueryBridges();
        initQueryBridgeLabelCounts();
#endif

        initGlobalLabelCounts(query_graph, Lq_counts, Lq_degrees);
        initGlobalLabelCounts(data_graph, Lg_counts, Lg_degrees);
        Lq_residual_counts = Lq_counts;
        Lg_residual_counts = Lg_counts;
        Lq_residual_degrees = Lq_degrees;
        Lg_residual_degrees = Lg_degrees;

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
#ifndef NDEBUG
        // if (kDebugInitialRoot >= 0) {
        //     assert((ui)kDebugInitialRoot < qn);
        //     root = (ui)kDebugInitialRoot;
        // }
#endif
        printf("Selected initial root: u=%u with %zu candidates\n", root, (size_t)candidates[root].size());

        for (ui v0 : candidates[root]) {
#ifndef NDEBUG
            recordBranchOrderDebug(root);
#endif
            updateResidualLabelCounts(query_graph, root, mapped_q, Lq_residual_counts, Lq_residual_degrees, false);
            updateResidualLabelCounts(data_graph, v0, mapped_g, Lg_residual_counts, Lg_residual_degrees, false);
            mapped_q[root] = (int)v0;
            mapped_g[v0] = (int)root;
            in_Mq[root] = 1;
            part_M.push_back({ root, v0 });

            updateFrontier(root, true);

            dfs(0, (int)root, nullptr, nullptr);

            updateFrontier(root, false);

            mapped_q[root] = -1;
            mapped_g[v0] = -1;
            in_Mq[root] = 0;
            part_M.pop_back();
            updateResidualLabelCounts(query_graph, root, mapped_q, Lq_residual_counts, Lq_residual_degrees, true);
            updateResidualLabelCounts(data_graph, v0, mapped_g, Lg_residual_counts, Lg_residual_degrees, true);
        }

        stats.dfs_time = t_search.elapsed();
        stats.total_time = stats.init_time + stats.dfs_time;
    }

    struct TimeStats {
        long long total_time = 0;
        // init breakdown
        long long init_time = 0;
        long long filter_time = 0;
        long long filter_nlf_time = 0;
        long long filter_bridge_time = 0;
        long long filter_spoke_time = 0;
        long long filter_onehop_time = 0;
        ui filter_candidate_count = 0;
        // search breakdown
        long long dfs_time = 0;
        long long lb_time = 0;       // computeLowerBound
        long long lb_state_time = 0;
        long long lb_sf_time = 0;
        long long lb_uu_label_time = 0;
        long long lb_uu_unsupported_time = 0;
        long long lb_light_spoke_time = 0;
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
#ifndef NDEBUG
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };
#endif

        long long lb_accounted_time = stats.lb_state_time + stats.lb_sf_time
            + stats.lb_uu_label_time + stats.lb_uu_unsupported_time
            + stats.lb_light_spoke_time;
        long long lb_other_time = stats.lb_time > lb_accounted_time
            ? stats.lb_time - lb_accounted_time : 0;
        long long search_accounted_time = stats.lb_time + stats.frontier_time + stats.branch_time;
        long long search_other_time = stats.dfs_time > search_accounted_time
            ? stats.dfs_time - search_accounted_time : 0;

        printf("\n--- CDE-Edge-IE Time Analysis ---\n");
#ifdef NDEBUG
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms\n", stats.init_time / 1000.0);
        printf("  - Filter Time:     %.4lf ms\n", stats.filter_time / 1000.0);
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms\n", stats.dfs_time / 1000.0);
        printf("  - LowerBound Time: %.4lf ms\n", stats.lb_time / 1000.0);
        printf("  - Frontier Time:   %.4lf ms\n", stats.frontier_time / 1000.0);
        printf("  - Branch Time:     %.4lf ms\n", stats.branch_time / 1000.0);
        printf("  - Search Other:    %.4lf ms\n", search_other_time / 1000.0);
#else
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("    - NLF:           %.4lf ms (%.2f%% of Filter)\n", stats.filter_nlf_time / 1000.0, pct(stats.filter_nlf_time, stats.filter_time));
        printf("    - Bridge:        %.4lf ms (%.2f%% of Filter)\n", stats.filter_bridge_time / 1000.0, pct(stats.filter_bridge_time, stats.filter_time));
        printf("    - Spoke:         %.4lf ms (%.2f%% of Filter)\n", stats.filter_spoke_time / 1000.0, pct(stats.filter_spoke_time, stats.filter_time));
        printf("    - OneHop:        %.4lf ms (%.2f%% of Filter)\n", stats.filter_onehop_time / 1000.0, pct(stats.filter_onehop_time, stats.filter_time));
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
        printf("    - State:         %.4lf ms (%.2f%% of LB)\n", stats.lb_state_time / 1000.0, pct(stats.lb_state_time, stats.lb_time));
        printf("    - S-F:           %.4lf ms (%.2f%% of LB)\n", stats.lb_sf_time / 1000.0, pct(stats.lb_sf_time, stats.lb_time));
        printf("    - U-U Label:     %.4lf ms (%.2f%% of LB)\n", stats.lb_uu_label_time / 1000.0, pct(stats.lb_uu_label_time, stats.lb_time));
        printf("    - U-U Unsup:     %.4lf ms (%.2f%% of LB)\n", stats.lb_uu_unsupported_time / 1000.0, pct(stats.lb_uu_unsupported_time, stats.lb_time));
        printf("    - Light Spoke:   %.4lf ms (%.2f%% of LB)\n", stats.lb_light_spoke_time / 1000.0, pct(stats.lb_light_spoke_time, stats.lb_time));
        printf("    - Other:         %.4lf ms (%.2f%% of LB)\n", lb_other_time / 1000.0, pct(lb_other_time, stats.lb_time));
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
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.dfs_time));
#endif
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
    enum { kDebugInitialRoot = 4 };
    // Internal debug knob: set to 0 to disable branch-order prefix counting.
    enum { kDebugBranchOrderDepth = 1 };
#endif

    struct ActiveEdge {
        ui u = 0;       // unmatched endpoint
        ui anchor = 0;  // matched endpoint
        ui anchor_support = std::numeric_limits<ui>::max();
        ui live_anchor_count = 0;
        ui query_degree = 0;
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
    vector<vector<char>> q_bridge_matrix;
    vector<pair<ui, ui>> q_bridge_edges;
    vector<vector<ui>> Lq_bridge_counts, Lq_non_bridge_counts;
    vector<ui> Lq_bridge_degrees, Lq_non_bridge_degrees;

    vector<MyBitset> candidates;
    vector<vector<ui>> Lq_counts, Lg_counts;
    vector<ui>  Lq_degrees, Lg_degrees;
    vector<vector<ui>> Lq_residual_counts, Lg_residual_counts;
    vector<ui>  Lq_residual_degrees, Lg_residual_degrees;

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
    vector<ui>  lb_component_mark;
    ui          lb_component_token;
    vector<int> lb_km_lx, lb_km_ly, lb_km_mx, lb_km_my;
    vector<int> lb_km_slack, lb_km_slackmy;
    vector<int> lb_km_prev, lb_km_queue;
    vector<char> lb_km_vis_x, lb_km_vis_y;
    vector<int> lb_v_to_col;
    vector<vector<ui>> lb_light_spoke_adj;
    vector<char> lb_light_spoke_is_bridge;
    vector<ui> lb_light_spoke_right_vertices;

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

    inline void insertXCand(ui u, ui v)
    {
        x_cand[u].insert(v);
    }

    inline void removeXCand(ui u, ui v)
    {
        x_cand[u].remove(v);
    }

    void updateResidualLabelCounts(const Graph *g, ui vertex, const vector<int> &mapped_status,
        vector<vector<ui>> &counts, vector<ui> &degrees, bool is_add)
    {
        LabelID label = g->getVertexLabel(vertex);
        if (label < 0 || (ui)label >= label_count) {
            return;
        }

        ui deg = 0;
        const ui *neighbors = g->getVertexNeighbors(vertex, deg);
        for (ui i = 0; i < deg; ++i) {
            ui neighbor = neighbors[i];
            if (mapped_status[neighbor] != -1) {
                continue;
            }

            if (is_add) {
                counts[neighbor][(ui)label]++;
                degrees[neighbor]++;
            }
            else {
                assert(counts[neighbor][(ui)label] > 0);
                assert(degrees[neighbor] > 0);
                counts[neighbor][(ui)label]--;
                degrees[neighbor]--;
            }
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
        q_bridge_matrix.clear();
        q_bridge_edges.clear();
        Lq_bridge_counts.clear();
        Lq_non_bridge_counts.clear();
        Lq_bridge_degrees.clear();
        Lq_non_bridge_degrees.clear();
        Lq_residual_counts.clear();
        Lg_residual_counts.clear();
        Lq_residual_degrees.clear();
        Lg_residual_degrees.clear();

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
        lb_component_mark.assign(qn, 0);
        lb_component_token = 1;
        lb_km_lx.clear();
        lb_km_ly.clear();
        lb_km_mx.clear();
        lb_km_my.clear();
        lb_km_slack.clear();
        lb_km_slackmy.clear();
        lb_km_prev.clear();
        lb_km_queue.clear();
        lb_km_vis_x.clear();
        lb_km_vis_y.clear();
        lb_v_to_col.assign(gn, -1);
        lb_light_spoke_adj.assign(qn, vector<ui>());
        lb_light_spoke_is_bridge.assign(qn, 0);
        lb_light_spoke_right_vertices.clear();

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

    void markQueryBridge(ui u, ui v)
    {
        if (q_bridge_matrix[u][v]) {
            return;
        }

        q_bridge_matrix[u][v] = 1;
        q_bridge_matrix[v][u] = 1;
        if (u < v) {
            q_bridge_edges.push_back({ u, v });
        }
        else {
            q_bridge_edges.push_back({ v, u });
        }
    }

    void dfsQueryBridges(ui u, ui parent, vector<int> &dfn, vector<int> &low, int &tim)
    {
        dfn[u] = low[u] = ++tim;

        for (ui v : q_neighbors[u]) {
            if (dfn[v] == 0) {
                dfsQueryBridges(v, u, dfn, low, tim);
                low[u] = std::min(low[u], low[v]);
                if (low[v] > dfn[u]) {
                    markQueryBridge(u, v);
                }
            }
            else if (v != parent) {
                low[u] = std::min(low[u], dfn[v]);
            }
        }
    }

    void identifyQueryBridges()
    {
        q_bridge_matrix.assign(qn, vector<char>(qn, 0));
        q_bridge_edges.clear();

        vector<int> dfn(qn, 0);
        vector<int> low(qn, 0);
        int tim = 0;

        for (ui u = 0; u < qn; ++u) {
            if (dfn[u] == 0) {
                dfsQueryBridges(u, qn, dfn, low, tim);
            }
        }
    }

    void initQueryBridgeLabelCounts()
    {
        Lq_bridge_counts.assign(qn, vector<ui>(label_count, 0));
        Lq_non_bridge_counts.assign(qn, vector<ui>(label_count, 0));
        Lq_bridge_degrees.assign(qn, 0);
        Lq_non_bridge_degrees.assign(qn, 0);

        for (ui u = 0; u < qn; ++u) {
            for (ui u1 : q_neighbors[u]) {
                LabelID label = query_graph->getVertexLabel(u1);
                if (label < 0 || (ui)label >= label_count) {
                    continue;
                }

                if (q_bridge_matrix[u][u1]) {
                    Lq_bridge_counts[u][(ui)label]++;
                    Lq_bridge_degrees[u]++;
                }
                else {
                    Lq_non_bridge_counts[u][(ui)label]++;
                    Lq_non_bridge_degrees[u]++;
                }
            }
        }
    }

    ui countEdgeCandidates(ui u, ui anchor) const
    {
        ui count = 0;
        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
        for (ui i = 0; i < deg; ++i) {
            ui v = nbrs[i];
            if (!candidates[u].contains(v)) continue;
            if (mapped_g[v] != -1) continue;
            if (x_cand[u].contains(v)) continue;
            count++;
        }
        return count;
    }

    ui countLiveAnchors(ui u) const
    {
        ui count = 0;
        for (ui anchor : q_neighbors[u]) {
            if (in_Mq[anchor] && !is_excluded[u][anchor]) {
                count++;
            }
        }
        return count;
    }

    // Active-edge scoring order:
    // 1. Smaller cand(u, anchor), also called anchor support.
    // 2. More live anchors.
    // 3. Higher query degree.
    // 4. Smaller id.
    bool isBetterActiveEdge(const ActiveEdge &lhs, const ActiveEdge &rhs) const
    {
        if (lhs.anchor_support != rhs.anchor_support) {
            return lhs.anchor_support < rhs.anchor_support;
        }
        if (lhs.live_anchor_count != rhs.live_anchor_count) {
            return lhs.live_anchor_count > rhs.live_anchor_count;
        }
        if (lhs.query_degree != rhs.query_degree) {
            return lhs.query_degree > rhs.query_degree;
        }
        if (lhs.u != rhs.u) {
            return lhs.u < rhs.u;
        }
        return lhs.anchor < rhs.anchor;
    }

    bool collectTopActiveEdges(const vector<ui> &component_frontier,
        ui max_count, vector<ActiveEdge> &top_edges) const
    {
        top_edges.clear();
        if (max_count == 0) {
            return false;
        }

        for (ui u : component_frontier) {
            if (u >= qn || in_Mq[u] || frontier_pos[u] == -1) {
                continue;
            }

            ui live_anchor_count = countLiveAnchors(u);
            if (live_anchor_count == 0) {
                continue;
            }

            for (ui anchor : q_neighbors[u]) {
                if (!in_Mq[anchor] || is_excluded[u][anchor]) {
                    continue;
                }

                ActiveEdge edge;
                edge.u = u;
                edge.anchor = anchor;
                edge.anchor_support = countEdgeCandidates(u, anchor);
                edge.live_anchor_count = live_anchor_count;
                edge.query_degree = (ui)q_neighbors[u].size();
                top_edges.push_back(edge);
            }
        }

        if (top_edges.empty()) {
            return false;
        }

        sort(top_edges.begin(), top_edges.end(), [&](const ActiveEdge &lhs, const ActiveEdge &rhs) {
            return isBetterActiveEdge(lhs, rhs);
            });
        if (top_edges.size() > max_count) {
            top_edges.resize(max_count);
        }
        return true;
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
        vector<int>         spoke_match_left;
        vector<bool>        spoke_vis;
        vector<char>        spoke_is_bridge;
        MyBitset            onehop_vis;

        explicit CandidateFilter(MatchingSolver &solver)
            : solver(solver),
            spoke_matrix(solver.qn, vector<ui>()),
            spoke_match(solver.max_g_deg, -1),
            spoke_match_left(solver.qn, -1),
            spoke_vis(solver.max_g_deg, false),
            spoke_is_bridge(solver.qn, 0),
            onehop_vis((int)solver.gn)
        {
        }

        bool run()
        {
            Timer t;
            // -- NLF Filtering --
            bool ok = filterNLF();
            solver.stats.filter_nlf_time += t.elapsed();
            if (!ok) return false;
#ifdef ENABLE_BRIDGE_FILTERING
            t.restart();
            ok = filterBridgeCandidates();
            solver.stats.filter_bridge_time += t.elapsed();
            if (!ok) return false;
#endif

            // --- Spoke Filtering ---
#ifdef ENABLE_SPOKE_FILTERING
            t.restart();
            ok = filterSpoke();
            solver.stats.filter_spoke_time += t.elapsed();
            if (!ok) return false;
#ifdef ENABLE_BRIDGE_FILTERING
            t.restart();
            ok = filterBridgeCandidates();
            solver.stats.filter_bridge_time += t.elapsed();
            if (!ok) return false;
#endif
#endif

            // --- One-Hop Filtering ---
#ifdef ENABLE_ONEHOP_FILTERING
            t.restart();
            ok = filterOneHop();
            solver.stats.filter_onehop_time += t.elapsed();
            if (!ok) return false;
#ifdef ENABLE_BRIDGE_FILTERING
            t.restart();
            ok = filterBridgeCandidates();
            solver.stats.filter_bridge_time += t.elapsed();
            if (!ok) return false;
#endif
#endif

            solver.stats.filter_candidate_count = 0;
            for (ui u = 0; u < solver.qn; ++u) {
                solver.stats.filter_candidate_count += (ui)solver.candidates[u].size();
            }

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
#ifdef ENABLE_BRIDGE_FILTERING
                ui bridge_need = solver.Lq_bridge_counts[u][i];
                ui data_count = solver.Lg_counts[v][i];
                if (bridge_need > data_count) {
                    return solver.threshold + 1;
                }

                ui non_bridge_need = solver.Lq_non_bridge_counts[u][i];
                ui non_bridge_available = data_count - bridge_need;
                if (non_bridge_need > non_bridge_available) {
                    diff += (non_bridge_need - non_bridge_available);
                }
#else
                if (solver.Lq_counts[u][i] > solver.Lg_counts[v][i]) {
                    diff += (solver.Lq_counts[u][i] - solver.Lg_counts[v][i]);
                }
#endif
            }
            return diff;
            }

        bool filterNLF()
        {
            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = solver.query_graph->getVertexLabel(u);
                for (ui v = 0; v < solver.gn; ++v) {
                    if (lu != solver.data_graph->getVertexLabel(v)) continue;
#ifdef ENABLE_BRIDGE_FILTERING
                    if (solver.Lq_bridge_degrees[u] > solver.Lg_degrees[v]) continue;
                    ui non_bridge_degree_capacity = solver.Lg_degrees[v] - solver.Lq_bridge_degrees[u];
                    if (solver.Lq_non_bridge_degrees[u] > non_bridge_degree_capacity + solver.threshold) continue;
#else
                    if (solver.Lq_degrees[u] > solver.Lg_degrees[v] + solver.threshold) continue;
#endif
                    if (computeNLF(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

#ifdef ENABLE_BRIDGE_FILTERING
        bool hasBridgeCandidateSupport(ui v, ui other_u) const
        {
            ui deg = 0;
            const ui *nbrs = solver.data_graph->getVertexNeighbors(v, deg);
            for (ui i = 0; i < deg; ++i) {
                if (solver.candidates[other_u].contains(nbrs[i])) {
                    return true;
                }
            }
            return false;
        }

        bool pruneBridgeSide(ui u, ui other_u, bool &changed)
        {
            vector<ui> to_remove;
            for (ui v : solver.candidates[u]) {
                if (!hasBridgeCandidateSupport(v, other_u)) {
                    to_remove.push_back(v);
                }
            }

            if (to_remove.empty()) {
                return true;
            }

            changed = true;
            for (ui v : to_remove) {
                solver.candidates[u].remove(v);
            }
            return !solver.candidates[u].empty();
        }

        bool filterBridgeCandidates()
        {
            if (solver.q_bridge_edges.empty()) {
                return true;
            }

            bool changed = true;
            while (changed) {
                changed = false;
                for (const pair<ui, ui> &edge : solver.q_bridge_edges) {
                    if (!pruneBridgeSide(edge.first, edge.second, changed)) return false;
                    if (!pruneBridgeSide(edge.second, edge.first, changed)) return false;
                }
            }
            return true;
        }
#endif

        bool canSkipApproxLocalFilter(ui u) const
        {
#ifdef ENABLE_BRIDGE_FILTERING
            if (solver.Lq_bridge_degrees[u] != 0) {
                return false;
            }
#endif
            return solver.q_neighbors[u].size() <= solver.threshold;
        }

        bool dfsMatchSpoke(ui left_idx, const vector<vector<ui>> &adj)
        {
            for (ui right_idx : adj[left_idx]) {
                if (spoke_vis[right_idx]) continue;
                spoke_vis[right_idx] = true;
                if (spoke_match[right_idx] < 0 ||
                    dfsMatchSpoke((ui)spoke_match[right_idx], adj)) {
                    spoke_match[right_idx] = (int)left_idx;
                    spoke_match_left[left_idx] = (int)right_idx;
                    return true;
                }
            }
            return false;
        }

        ui computeMaxMatchSpoke(const vector<vector<ui>> &adj, ui left_size, ui right_size)
        {
            ui mu = 0;
            std::fill(spoke_match.begin(), spoke_match.begin() + right_size, -1);
            std::fill(spoke_match_left.begin(), spoke_match_left.begin() + left_size, -1);
            for (ui i = 0; i < left_size; ++i) {
                std::fill(spoke_vis.begin(), spoke_vis.begin() + right_size, false);
                if (dfsMatchSpoke(i, adj)) mu++;
            }
            return mu;
        }

#ifdef ENABLE_BRIDGE_FILTERING
        bool computeBridgeAwareSpokeMissing(const vector<vector<ui>> &adj, ui left_size,
            ui right_size, ui &missing_non_bridge)
        {
            std::fill(spoke_match.begin(), spoke_match.begin() + right_size, -1);
            std::fill(spoke_match_left.begin(), spoke_match_left.begin() + left_size, -1);

            ui non_bridge_count = 0;
            for (ui i = 0; i < left_size; ++i) {
                if (spoke_is_bridge[i]) {
                    std::fill(spoke_vis.begin(), spoke_vis.begin() + right_size, false);
                    if (!dfsMatchSpoke(i, adj)) {
                        return false;
                    }
                }
                else {
                    non_bridge_count++;
                }
            }

            ui non_bridge_match = 0;
            for (ui i = 0; i < left_size; ++i) {
                if (spoke_is_bridge[i]) {
                    continue;
                }

                std::fill(spoke_vis.begin(), spoke_vis.begin() + right_size, false);
                if (dfsMatchSpoke(i, adj)) {
                    non_bridge_match++;
                }
            }

            missing_non_bridge = non_bridge_count - non_bridge_match;
            return true;
        }
#endif

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
#ifdef ENABLE_BRIDGE_FILTERING
                spoke_is_bridge[i] = solver.q_bridge_matrix[u][u1];
#endif
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        spoke_matrix[i].push_back(j);
                    }
                }
            }

#ifdef ENABLE_BRIDGE_FILTERING
            ui missing_non_bridge = 0;
            if (!computeBridgeAwareSpokeMissing(spoke_matrix, deg_u, deg_v, missing_non_bridge)) {
                return solver.threshold + 1;
            }
            return missing_non_bridge;
#else
            ui match_size = computeMaxMatchSpoke(spoke_matrix, deg_u, deg_v);
            return deg_u - match_size;
#endif
        }

        bool filterSpoke()
        {
            queue<ui> q;
            vector<char> in_q(solver.qn, 1);
            for (ui u = 0; u < solver.qn; ++u) q.push(u);

            while (!q.empty()) {
                ui u = q.front(); q.pop();
                in_q[u] = 0;

                if (canSkipApproxLocalFilter(u)) {
                    continue;
                }

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
        vector<ui> buildOneHopOrder(const vector<ui> &u_neighbors,
            const vector<vector<ui>> &cand, const vector<char> &spoke_bridge) const
        {
            ui deg_u = (ui)u_neighbors.size();
            vector<ui> ord(deg_u);
            iota(ord.begin(), ord.end(), 0);

            sort(ord.begin(), ord.end(), [&](ui a, ui b) {
                if (spoke_bridge[a] != spoke_bridge[b]) {
                    return spoke_bridge[a] > spoke_bridge[b];
                }
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
            const vector<vector<ui>> &cand, const vector<char> &spoke_bridge) const
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
                if (!has_free && spoke_bridge[u1]) {
                    return solver.threshold + 1;
                }
            }
            return rem;
        }

        // DFS for one-hop matching
        // Branch 1: skip u1, adding one missing spoke edge
        // Branch 2: match u1 to v1, adding newly missing inner edges
        bool dfsOneHop(ui pos, vector<int> &state, ui cost,
            const vector<ui> &ord, const vector<ui> &u_neighbors, const vector<vector<ui>> &cand,
            const vector<char> &spoke_bridge)
        {
            if (cost > solver.threshold) return false;
            if (pos == ord.size()) return true;

            ui rem_lb = computeRemainLBOneHop(pos, ord, cand, spoke_bridge);
            if (cost + rem_lb > solver.threshold) return false;

            ui i = ord[pos];
            ui u1 = u_neighbors[i];

            // branch 1: skip u1
            if (!spoke_bridge[i]) {
                state[i] = -1;
                if (dfsOneHop(pos + 1, state, cost + 1, ord, u_neighbors, cand, spoke_bridge)) {
                    state[i] = -2;
                    return true;
                }
            }

            // branch 2: match u1 to v1
            for (ui v1 : cand[i]) {
                if (onehop_vis.contains(v1)) continue;

                ui delta_inner = 0;
                bool missing_bridge = false;
                for (ui j = 0; j < u_neighbors.size(); ++j) {
                    if (state[j] < 0) continue;
                    ui u2 = u_neighbors[j];
                    if (!solver.q_matrix[u1][u2]) continue;

                    ui v2 = (ui)state[j];
                    if (!solver.data_graph->hasEdge(v1, v2)) {
#ifdef ENABLE_BRIDGE_FILTERING
                        if (solver.q_bridge_matrix[u1][u2]) {
                            missing_bridge = true;
                            break;
                        }
#endif
                        delta_inner++;
                    }
                }

                if (missing_bridge) continue;
                if (cost + delta_inner > solver.threshold) continue;

                onehop_vis.insert(v1);
                state[i] = (int)v1;

                if (dfsOneHop(pos + 1, state, cost + delta_inner, ord, u_neighbors, cand, spoke_bridge)) {
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
            vector<char> spoke_bridge(deg_u, 0);

            for (ui i = 0; i < deg_u; ++i) {
                ui u1 = u_neighbors[i];
#ifdef ENABLE_BRIDGE_FILTERING
                spoke_bridge[i] = solver.q_bridge_matrix[u][u1];
#endif
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        cand[i].push_back(v1);
                    }
                }
            }

            // matching order
            vector<ui> ord = buildOneHopOrder(u_neighbors, cand, spoke_bridge);

            // the current DFS state of u_neighbors[i]:
            // -2: unprocessed
            // -1: skipped / unmatched
            // >=0: matched to the corresponding data vertex
            vector<int> state(deg_u, -2);

            return dfsOneHop(0, state, 0, ord, u_neighbors, cand, spoke_bridge);
        }

        bool filterOneHop()
        {
            vector<pair<ui, ui>> to_remove;

            for (ui u = 0; u < solver.qn; ++u) {
                if (solver.q_neighbors[u].size() <= 1) continue;
                if (countInnerEdges(u) == 0) continue;

                for (ui v : solver.candidates[u]) {
                    ui missing_edges = computeLBSpoke(u, v);
                    if (missing_edges > solver.threshold) {
                        to_remove.push_back({ u, v });
                        continue;
                    }
                    if (solver.threshold - missing_edges > (ui)ONEHOP_FILTER_MISSING_GAP) {
                        continue;
                    }
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
            vector<ActiveEdge> top_edges;
            bool found_edge = solver.collectTopActiveEdges(solver.active_frontier, 1, top_edges);
            assert(found_edge);
            (void)found_edge;
            solver.stats.frontier_select_time += t_select.elapsed();

            Timer t_component;
            collectComponent(top_edges.front().u, state);
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

        // Count current live candidates for u without conditioning on one anchor.
        ui countScoreCandidates(ui u) const
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

        ui countScoreAnchorSupport(ui u, ui anchor) const
        {
            return countAnchorSupport(u, anchor);
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
            (void)u;
            return false;
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

        // Smaller current anchor support, more live anchors, higher query degree, smaller vertex id.
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
    struct LightweightSpokeLowerBound {
        MatchingSolver &solver;

        explicit LightweightSpokeLowerBound(MatchingSolver &solver) : solver(solver) {}

        bool shouldPrune(ui mapped_cost, ui u, ui v)
        {
            if (mapped_cost > solver.threshold) {
                return true;
            }

            ui remaining_budget = solver.threshold - mapped_cost;
            ui missing = computeUnmatchedSpokeMissing(u, v);
            return missing > remaining_budget;
        }

    private:
        ui boundInf() const
        {
            return (ui)INF;
        }

        bool isCandidateAvailable(ui u, ui v, ui chosen_v) const
        {
            return v != chosen_v && solver.mapped_g[v] == -1 &&
                !solver.x_cand[u].contains(v);
        }

        void nextDataToken()
        {
            if (++solver.lb_data_frontier_token == 0) {
                std::fill(solver.lb_data_frontier_mark.begin(),
                    solver.lb_data_frontier_mark.end(), 0);
                solver.lb_data_frontier_token = 1;
            }
        }

        void nextSeenToken()
        {
            if (++solver.lb_seen_token == 0) {
                std::fill(solver.lb_seen_right.begin(), solver.lb_seen_right.end(), 0);
                solver.lb_seen_token = 1;
            }
        }

        bool findAugment(ui left_idx)
        {
            const vector<ui> &adj = solver.lb_light_spoke_adj[left_idx];
            for (ui v : adj) {
                if (solver.lb_seen_right[v] == solver.lb_seen_token) {
                    continue;
                }
                solver.lb_seen_right[v] = solver.lb_seen_token;

                if (solver.lb_match_right[v] < 0 ||
                    findAugment((ui)solver.lb_match_right[v])) {
                    solver.lb_match_right[v] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        ui computeUnmatchedSpokeMissing(ui u, ui v)
        {
            ui deg_v = 0;
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            nextDataToken();
            solver.lb_light_spoke_right_vertices.clear();

            ui left_size = 0;
#ifdef ENABLE_BRIDGE_FILTERING
            ui non_bridge_count = 0;
#endif

            for (ui w : solver.q_neighbors[u]) {
                if (solver.in_Mq[w] || solver.is_excluded[u][w]) {
                    continue;
                }

                vector<ui> &adj = solver.lb_light_spoke_adj[left_size];
                adj.clear();

#ifdef ENABLE_BRIDGE_FILTERING
                solver.lb_light_spoke_is_bridge[left_size] = solver.q_bridge_matrix[u][w];
                if (!solver.lb_light_spoke_is_bridge[left_size]) {
                    non_bridge_count++;
                }
#endif

                for (ui i = 0; i < deg_v; ++i) {
                    ui t = v_neighbors[i];
                    if (!isCandidateAvailable(w, t, v) ||
                        !solver.candidates[w].contains(t)) {
                        continue;
                    }

                    adj.push_back(t);
                    if (solver.lb_data_frontier_mark[t] != solver.lb_data_frontier_token) {
                        solver.lb_data_frontier_mark[t] = solver.lb_data_frontier_token;
                        solver.lb_light_spoke_right_vertices.push_back(t);
                    }
                }

                left_size++;
            }

            if (left_size == 0) {
                return 0;
            }

            for (ui t : solver.lb_light_spoke_right_vertices) {
                solver.lb_match_right[t] = -1;
            }

#ifdef ENABLE_BRIDGE_FILTERING
            for (ui i = 0; i < left_size; ++i) {
                if (!solver.lb_light_spoke_is_bridge[i]) {
                    continue;
                }

                nextSeenToken();
                if (!findAugment(i)) {
                    return boundInf();
                }
            }

            ui non_bridge_match = 0;
            for (ui i = 0; i < left_size; ++i) {
                if (solver.lb_light_spoke_is_bridge[i]) {
                    continue;
                }

                nextSeenToken();
                if (findAugment(i)) {
                    non_bridge_match++;
                }
            }

            return non_bridge_count - non_bridge_match;
#else
            ui match_size = 0;
            for (ui i = 0; i < left_size; ++i) {
                nextSeenToken();
                if (findAugment(i)) {
                    match_size++;
                }
            }

            return left_size - match_size;
#endif
        }
    };

    struct LowerBoundOptions {
        bool use_sf_assignment = false;
        bool use_uu_unsupported = false;
        bool use_edge_labels = false;
        bool restrict_data_neighbors_by_candidates = true;
    };

    struct LowerBoundComponentCache {
        vector<ui> vertices;
        vector<ui> frontier;
        ui lb = 0;
    };

    struct LowerBoundState {
        vector<LowerBoundComponentCache> components;

        void clear()
        {
            components.clear();
        }
    };

    struct LowerBoundPruner {
        MatchingSolver &solver;

        explicit LowerBoundPruner(MatchingSolver &solver) : solver(solver) {}

        bool shouldPrune(ui current_miss, const vector<ui> &frontier_vertices,
            const LowerBoundState *parent_cache = nullptr, LowerBoundState *out_cache = nullptr)
        {
            (void)frontier_vertices;
            LowerBoundOptions options;
            return computeLowerBound(current_miss, options, parent_cache, out_cache) > solver.threshold;
        }

    private:
        ui boundInf() const
        {
            return (ui)INF;
        }

        ui addBound(ui lhs, ui rhs) const
        {
            if (lhs >= boundInf() || rhs >= boundInf()) {
                return boundInf();
            }
            if (lhs > boundInf() - rhs) {
                return boundInf();
            }
            return lhs + rhs;
        }

        bool isCandidateAvailable(ui u, ui v) const
        {
            return solver.mapped_g[v] == -1 && !solver.x_cand[u].contains(v);
        }

        bool hasAvailableCandidate(ui u) const
        {
            for (ui v : solver.candidates[u]) {
                if (isCandidateAvailable(u, v)) {
                    return true;
                }
            }
            return false;
        }

        void collectLiveAnchors(ui u, vector<ui> &anchors) const
        {
            anchors.clear();
            for (ui a : solver.q_neighbors[u]) {
                if (solver.in_Mq[a] && !solver.is_excluded[u][a]) {
                    anchors.push_back(a);
                }
            }
        }

        void collectActiveFrontier(vector<ui> &frontier_vertices) const
        {
            frontier_vertices.clear();
            frontier_vertices.reserve(solver.active_frontier.size());
            for (ui u : solver.active_frontier) {
                if (!solver.in_Mq[u] && solver.frontier_pos[u] != -1) {
                    frontier_vertices.push_back(u);
                }
            }
        }

        void collectUnmatched(vector<ui> &unmatched_vertices) const
        {
            unmatched_vertices.clear();
            unmatched_vertices.reserve(solver.qn - (ui)solver.part_M.size());
            for (ui u = 0; u < solver.qn; ++u) {
                if (!solver.in_Mq[u]) {
                    unmatched_vertices.push_back(u);
                }
            }
        }

        struct LBComponent {
            vector<ui> vertices;
            vector<ui> frontier;
        };

        ui doubledLimit(ui limit) const
        {
            if (limit >= boundInf() / 2) {
                return boundInf();
            }
            return limit * 2;
        }

        bool doubledExceedsLimit(ui doubled_value, ui limit) const
        {
            if (doubled_value >= boundInf()) {
                return true;
            }
            return (doubled_value + 1) / 2 > limit;
        }

        void nextComponentToken()
        {
            if (++solver.lb_component_token == 0) {
                std::fill(solver.lb_component_mark.begin(), solver.lb_component_mark.end(), 0);
                solver.lb_component_token = 1;
            }
        }

        void nextDataFrontierToken()
        {
            if (++solver.lb_data_frontier_token == 0) {
                std::fill(solver.lb_data_frontier_mark.begin(), solver.lb_data_frontier_mark.end(), 0);
                solver.lb_data_frontier_token = 1;
            }
        }

        void collectUnmatchedComponents(vector<LBComponent> &components)
        {
            components.clear();
            nextComponentToken();

            queue<ui> q;
            for (ui u = 0; u < solver.qn; ++u) {
                if (solver.in_Mq[u] || solver.lb_component_mark[u] == solver.lb_component_token) {
                    continue;
                }

                components.push_back(LBComponent());
                LBComponent &component = components.back();
                solver.lb_component_mark[u] = solver.lb_component_token;
                q.push(u);

                while (!q.empty()) {
                    ui curr = q.front();
                    q.pop();

                    component.vertices.push_back(curr);
                    if (solver.frontier_pos[curr] != -1) {
                        component.frontier.push_back(curr);
                    }

                    for (ui nbr : solver.q_neighbors[curr]) {
                        if (solver.in_Mq[nbr] || solver.is_excluded[curr][nbr] ||
                            solver.lb_component_mark[nbr] == solver.lb_component_token) {
                            continue;
                        }
                        solver.lb_component_mark[nbr] = solver.lb_component_token;
                        q.push(nbr);
                    }
                }

                sort(component.vertices.begin(), component.vertices.end());
                sort(component.frontier.begin(), component.frontier.end());
            }
        }

        const LowerBoundComponentCache *findCachedComponent(
            const LowerBoundState *parent_cache, const LBComponent &component) const
        {
            if (parent_cache == nullptr) {
                return nullptr;
            }

            for (const LowerBoundComponentCache &cached : parent_cache->components) {
                if (cached.vertices == component.vertices && cached.frontier == component.frontier) {
                    return &cached;
                }
            }
            return nullptr;
        }

        ui computeResidualLabelDeficit(ui u, ui v) const
        {
            ui deficit = 0;
            for (ui label = 0; label < solver.label_count; ++label) {
                ui q_count = solver.Lq_residual_counts[u][label];
                ui g_count = solver.Lg_residual_counts[v][label];
                if (q_count > g_count) {
                    deficit += q_count - g_count;
                }
            }
            return deficit;
        }

        ui computeRowDoubledCost(ui u, ui v, const vector<ui> &anchors) const
        {
            ui sf = computeAnchorCutDelta(v, anchors);
            ui uu = computeResidualLabelDeficit(u, v);
            ui doubled_sf = addBound(sf, sf);
            return addBound(doubled_sf, uu);
        }

        ui computeRemainDoubled(const LBComponent &component, ui limit) const
        {
            ui sum = 0;
            for (ui u : component.vertices) {
                if (binary_search(component.frontier.begin(), component.frontier.end(), u)) {
                    continue;
                }

                ui best = boundInf();
                for (ui v : solver.candidates[u]) {
                    if (!isCandidateAvailable(u, v)) {
                        continue;
                    }

                    ui deficit = computeResidualLabelDeficit(u, v);
                    if (deficit < best) {
                        best = deficit;
                        if (best == 0) {
                            break;
                        }
                    }
                }

                if (best >= boundInf()) {
                    return boundInf();
                }

                sum = addBound(sum, best);
                if (doubledExceedsLimit(sum, limit)) {
                    return boundInf();
                }
            }
            return sum;
        }

        ui computeFrontierIndependentDoubled(const vector<ui> &rows,
            const vector<vector<ui>> &anchors_by_row, ui limit) const
        {
            ui sum = 0;
            for (ui i = 0; i < rows.size(); ++i) {
                ui u = rows[i];
                const vector<ui> &anchors = anchors_by_row[i];
                ui best = boundInf();

                for (ui v : solver.candidates[u]) {
                    if (!isCandidateAvailable(u, v)) {
                        continue;
                    }

                    ui cost = computeRowDoubledCost(u, v, anchors);
                    if (cost < best) {
                        best = cost;
                        if (best == 0) {
                            break;
                        }
                    }
                }

                if (best >= boundInf()) {
                    return boundInf();
                }

                sum = addBound(sum, best);
                if (doubledExceedsLimit(sum, limit)) {
                    return boundInf();
                }
            }
            return sum;
        }

        ui solveMWPM(const vector<vector<ui>> &cost_matrix, ui limit)
        {
            if (cost_matrix.empty()) {
                return 0;
            }

            int n = (int)cost_matrix.size();
            int m = (int)cost_matrix[0].size();
            if (n > m) {
                return limit + 1;
            }

            const int INF_INT = (int)boundInf();
            if (solver.lb_km_lx.size() < (size_t)n) {
                solver.lb_km_lx.resize(n);
                solver.lb_km_mx.resize(n);
                solver.lb_km_prev.resize(n);
                solver.lb_km_queue.resize(n);
                solver.lb_km_vis_x.resize(n);
            }
            if (solver.lb_km_ly.size() < (size_t)m) {
                solver.lb_km_ly.resize(m);
                solver.lb_km_my.resize(m);
                solver.lb_km_slack.resize(m);
                solver.lb_km_slackmy.resize(m);
                solver.lb_km_vis_y.resize(m);
            }

            fill(solver.lb_km_lx.begin(), solver.lb_km_lx.begin() + n, 0);
            fill(solver.lb_km_ly.begin(), solver.lb_km_ly.begin() + m, 0);
            fill(solver.lb_km_mx.begin(), solver.lb_km_mx.begin() + n, -1);
            fill(solver.lb_km_my.begin(), solver.lb_km_my.begin() + m, -1);

            unsigned long long row_lb = 0;
            for (int i = 0; i < n; ++i) {
                int min_val = INF_INT;
                for (int j = 0; j < m; ++j) {
                    if ((int)cost_matrix[i][j] < min_val) {
                        min_val = (int)cost_matrix[i][j];
                    }
                }
                if (min_val >= INF_INT) {
                    return limit + 1;
                }
                solver.lb_km_lx[i] = min_val;
                row_lb += (unsigned int)min_val;
                if (row_lb > limit) {
                    return limit + 1;
                }
            }

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    if (solver.lb_km_my[j] == -1 &&
                        (int)cost_matrix[i][j] == solver.lb_km_lx[i]) {
                        solver.lb_km_mx[i] = j;
                        solver.lb_km_my[j] = i;
                        break;
                    }
                }
            }

            for (int root = 0; root < n; ++root) {
                if (solver.lb_km_mx[root] != -1) {
                    continue;
                }

                fill(solver.lb_km_vis_x.begin(), solver.lb_km_vis_x.begin() + n, 0);
                fill(solver.lb_km_vis_y.begin(), solver.lb_km_vis_y.begin() + m, 0);

                for (int j = 0; j < m; ++j) {
                    solver.lb_km_slack[j] =
                        (int)cost_matrix[root][j] - solver.lb_km_lx[root] - solver.lb_km_ly[j];
                    solver.lb_km_slackmy[j] = root;
                }

                solver.lb_km_vis_x[root] = 1;
                int q_size = 0;
                solver.lb_km_queue[q_size++] = root;
                int target_y = -1;
                int last_x = -1;

                while (true) {
                    int q_ptr = 0;
                    while (q_ptr < q_size) {
                        int x = solver.lb_km_queue[q_ptr++];
                        for (int y = 0; y < m; ++y) {
                            if (solver.lb_km_vis_y[y] ||
                                (int)cost_matrix[x][y] != solver.lb_km_lx[x] + solver.lb_km_ly[y]) {
                                continue;
                            }

                            if (solver.lb_km_my[y] == -1) {
                                last_x = x;
                                target_y = y;
                                goto found_path;
                            }

                            solver.lb_km_vis_y[y] = 1;
                            int next_x = solver.lb_km_my[y];
                            if (solver.lb_km_vis_x[next_x]) {
                                continue;
                            }
                            solver.lb_km_vis_x[next_x] = 1;
                            solver.lb_km_prev[next_x] = x;
                            solver.lb_km_queue[q_size++] = next_x;
                            for (int k = 0; k < m; ++k) {
                                int slack = (int)cost_matrix[next_x][k] -
                                    solver.lb_km_lx[next_x] - solver.lb_km_ly[k];
                                if (!solver.lb_km_vis_y[k] && slack < solver.lb_km_slack[k]) {
                                    solver.lb_km_slack[k] = slack;
                                    solver.lb_km_slackmy[k] = next_x;
                                }
                            }
                        }
                    }

                    int delta = INF_INT;
                    for (int y = 0; y < m; ++y) {
                        if (!solver.lb_km_vis_y[y] && solver.lb_km_slack[y] < delta) {
                            delta = solver.lb_km_slack[y];
                        }
                    }
                    if (delta >= INF_INT) {
                        return limit + 1;
                    }

                    for (int x = 0; x < n; ++x) {
                        if (solver.lb_km_vis_x[x]) {
                            solver.lb_km_lx[x] += delta;
                        }
                    }
                    for (int y = 0; y < m; ++y) {
                        if (solver.lb_km_vis_y[y]) {
                            solver.lb_km_ly[y] -= delta;
                        }
                        else {
                            solver.lb_km_slack[y] -= delta;
                        }
                    }

                    q_size = 0;
                    for (int y = 0; y < m; ++y) {
                        if (solver.lb_km_vis_y[y] || solver.lb_km_slack[y] != 0) {
                            continue;
                        }

                        if (solver.lb_km_my[y] == -1) {
                            last_x = solver.lb_km_slackmy[y];
                            target_y = y;
                            goto found_path;
                        }

                        solver.lb_km_vis_y[y] = 1;
                        int next_x = solver.lb_km_my[y];
                        if (solver.lb_km_vis_x[next_x]) {
                            continue;
                        }
                        solver.lb_km_vis_x[next_x] = 1;
                        solver.lb_km_prev[next_x] = solver.lb_km_slackmy[y];
                        solver.lb_km_queue[q_size++] = next_x;
                        for (int k = 0; k < m; ++k) {
                            int slack = (int)cost_matrix[next_x][k] -
                                solver.lb_km_lx[next_x] - solver.lb_km_ly[k];
                            if (!solver.lb_km_vis_y[k] && slack < solver.lb_km_slack[k]) {
                                solver.lb_km_slack[k] = slack;
                                solver.lb_km_slackmy[k] = next_x;
                            }
                        }
                    }
                }

            found_path:
                while (true) {
                    int prev_y = solver.lb_km_mx[last_x];
                    solver.lb_km_mx[last_x] = target_y;
                    solver.lb_km_my[target_y] = last_x;
                    if (last_x == root) {
                        break;
                    }
                    target_y = prev_y;
                    last_x = solver.lb_km_prev[last_x];
                }
            }

            unsigned long long total_cost = 0;
            for (int i = 0; i < n; ++i) {
                if (solver.lb_km_mx[i] < 0) {
                    return limit + 1;
                }
                total_cost += cost_matrix[i][solver.lb_km_mx[i]];
                if (total_cost > limit) {
                    return limit + 1;
                }
            }
            return (ui)total_cost;
        }

        void resetDataColumns(const vector<ui> &cols)
        {
            for (ui v : cols) {
                solver.lb_v_to_col[v] = -1;
            }
        }

        ui computeFrontierMWPMDoubled(const LBComponent &component, ui limit)
        {
            if (component.frontier.empty()) {
                return 0;
            }

            vector<pair<LabelID, ui>> labeled_rows;
            labeled_rows.reserve(component.frontier.size());
            for (ui u : component.frontier) {
                labeled_rows.push_back({ solver.query_graph->getVertexLabel(u), u });
            }
            sort(labeled_rows.begin(), labeled_rows.end());

            ui front_doubled = 0;
            size_t ptr = 0;
            while (ptr < labeled_rows.size()) {
                LabelID label = labeled_rows[ptr].first;
                size_t qtr = ptr;
                while (qtr < labeled_rows.size() && labeled_rows[qtr].first == label) {
                    qtr++;
                }

                vector<ui> rows;
                vector<vector<ui>> anchors_by_row;
                rows.reserve(qtr - ptr);
                anchors_by_row.reserve(qtr - ptr);
                for (size_t i = ptr; i < qtr; ++i) {
                    rows.push_back(labeled_rows[i].second);
                    anchors_by_row.push_back(vector<ui>());
                    collectLiveAnchors(rows.back(), anchors_by_row.back());
                }

                nextDataFrontierToken();
                vector<ui> v_list;
                for (ui i = 0; i < rows.size(); ++i) {
                    ui u = rows[i];
                    for (ui anchor : anchors_by_row[i]) {
                        assert(solver.mapped_q[anchor] >= 0);
                        ui deg = 0;
                        const ui *nbrs = solver.data_graph->getVertexNeighbors(
                            (ui)solver.mapped_q[anchor], deg);
                        for (ui j = 0; j < deg; ++j) {
                            ui v = nbrs[j];
                            if (solver.data_graph->getVertexLabel(v) != label ||
                                !solver.candidates[u].contains(v) || !isCandidateAvailable(u, v)) {
                                continue;
                            }
                            if (solver.lb_data_frontier_mark[v] != solver.lb_data_frontier_token) {
                                solver.lb_data_frontier_mark[v] = solver.lb_data_frontier_token;
                                solver.lb_v_to_col[v] = (int)v_list.size();
                                v_list.push_back(v);
                            }
                        }
                    }
                }

                bool over_cap = false;
#if CDE_LB_MWPM_MAX_ROWS > 0
                over_cap = over_cap || rows.size() > (size_t)CDE_LB_MWPM_MAX_ROWS;
#endif
#if CDE_LB_MWPM_MAX_COLS > 0
                over_cap = over_cap || v_list.size() + rows.size() > (size_t)CDE_LB_MWPM_MAX_COLS;
#endif
                if (over_cap) {
                    resetDataColumns(v_list);
                    ui independent = computeFrontierIndependentDoubled(rows, anchors_by_row, limit);
                    if (independent >= boundInf()) {
                        return boundInf();
                    }
                    front_doubled = addBound(front_doubled, independent);
                    if (doubledExceedsLimit(front_doubled, limit)) {
                        return boundInf();
                    }
                    ptr = qtr;
                    continue;
                }

                size_t row_n = rows.size();
                size_t col_n = v_list.size() + row_n;
                const ui INF_COST = boundInf();
                vector<vector<ui>> cost_matrix(row_n, vector<ui>(col_n, INF_COST));

                for (size_t i = 0; i < row_n; ++i) {
                    ui u = rows[i];
                    const vector<ui> &anchors = anchors_by_row[i];
                    ui virtual_cost = boundInf();

                    for (ui v : solver.candidates[u]) {
                        if (!isCandidateAvailable(u, v)) {
                            continue;
                        }

                        ui cost = computeRowDoubledCost(u, v, anchors);
                        int col = solver.lb_v_to_col[v];
                        if (col >= 0) {
                            cost_matrix[i][(size_t)col] = std::min(cost_matrix[i][(size_t)col], cost);
                        }
                        else {
                            virtual_cost = std::min(virtual_cost, cost);
                        }
                    }

                    if (virtual_cost < INF_COST) {
                        for (size_t k = 0; k < row_n; ++k) {
                            cost_matrix[i][v_list.size() + k] = virtual_cost;
                        }
                    }
                }

                ui doubled_limit = doubledLimit(limit);
                if (front_doubled > doubled_limit) {
                    resetDataColumns(v_list);
                    return boundInf();
                }

                ui local_limit = doubled_limit - front_doubled;
                ui local_cost = solveMWPM(cost_matrix, local_limit);
                resetDataColumns(v_list);
                if (local_cost > local_limit) {
                    return boundInf();
                }

                front_doubled = addBound(front_doubled, local_cost);
                if (doubledExceedsLimit(front_doubled, limit)) {
                    return boundInf();
                }

                ptr = qtr;
            }
            return front_doubled;
        }

        ui computeComponentMWPMLowerBound(const LBComponent &component, ui limit)
        {
            ui remain_doubled = computeRemainDoubled(component, limit);
            if (remain_doubled >= boundInf()) {
                return boundInf();
            }

            ui doubled_limit = doubledLimit(limit);
            if (remain_doubled > doubled_limit) {
                return boundInf();
            }

            ui front_limit = (doubled_limit - remain_doubled + 1) / 2;
            ui front_doubled = computeFrontierMWPMDoubled(component, front_limit);
            if (front_doubled >= boundInf()) {
                return boundInf();
            }

            ui total_doubled = addBound(remain_doubled, front_doubled);
            if (doubledExceedsLimit(total_doubled, limit)) {
                return boundInf();
            }
            return (total_doubled + 1) / 2;
        }

        ui computeComponentMWPMLowerBound(ui current_miss,
            const LowerBoundState *parent_cache, LowerBoundState *out_cache)
        {
            if (out_cache != nullptr) {
                out_cache->clear();
            }

            if (current_miss > solver.threshold) {
                return current_miss;
            }

            vector<LBComponent> components;
            collectUnmatchedComponents(components);

            ui residual_limit = solver.threshold - current_miss;
            ui residual_lb = 0;
            for (const LBComponent &component : components) {
                ui component_limit = residual_limit >= residual_lb ? residual_limit - residual_lb : 0;
                ui component_lb = boundInf();

#ifdef CDE_LB_COMPONENT_MWPM_CACHE
                const LowerBoundComponentCache *cached = findCachedComponent(parent_cache, component);
                if (cached != nullptr) {
                    component_lb = cached->lb;
                }
#else
                (void)parent_cache;
#endif
                if (component_lb >= boundInf()) {
                    component_lb = computeComponentMWPMLowerBound(component, component_limit);
                }

                if (out_cache != nullptr) {
                    LowerBoundComponentCache cache_entry;
                    cache_entry.vertices = component.vertices;
                    cache_entry.frontier = component.frontier;
                    cache_entry.lb = component_lb;
                    out_cache->components.push_back(cache_entry);
                }

                residual_lb = addBound(residual_lb, component_lb);
                if (residual_lb > residual_limit) {
                    return addBound(current_miss, residual_lb);
                }
                }

            return addBound(current_miss, residual_lb);
            }

        ui computeLowerBound(ui current_miss, const LowerBoundOptions &options,
            const LowerBoundState *parent_cache, LowerBoundState *out_cache)
        {
            // The DFS cost already contains fixed S-S missing edges plus
            // excluded S-F anchor edges. The lower bound below is residual:
            // only live S-F anchors and U-U edges are estimated here.
            (void)options.use_edge_labels; // Graph currently exposes vertex labels only.

#ifdef CDE_LB_COMPONENT_MWPM
            (void)options;
            return computeComponentMWPMLowerBound(current_miss, parent_cache, out_cache);
#else
            (void)parent_cache;
            (void)out_cache;

            if (current_miss > solver.threshold) {
                return current_miss;
            }

            Timer t_part;
            vector<ui> unmatched_vertices;
            collectUnmatched(unmatched_vertices);
            for (ui u : unmatched_vertices) {
                if (!hasAvailableCandidate(u)) {
                    solver.stats.lb_state_time += t_part.elapsed();
                    return boundInf();
                }
            }

            vector<ui> frontier_vertices;
            collectActiveFrontier(frontier_vertices);

            vector<vector<ui>> anchors_by_frontier_idx(frontier_vertices.size());
            for (ui i = 0; i < frontier_vertices.size(); ++i) {
                collectLiveAnchors(frontier_vertices[i], anchors_by_frontier_idx[i]);
            }
            solver.stats.lb_state_time += t_part.elapsed();

            t_part.restart();
            ui sf_lb = options.use_sf_assignment
                ? computeSFAssignmentLowerBound(frontier_vertices, anchors_by_frontier_idx)
                : computeSFMinLowerBound(frontier_vertices, anchors_by_frontier_idx);
            solver.stats.lb_sf_time += t_part.elapsed();
            if (sf_lb >= boundInf()) {
                return boundInf();
            }

            t_part.restart();
            ui uu_lb = computeUULabelLowerBound(unmatched_vertices, options);
            solver.stats.lb_uu_label_time += t_part.elapsed();
            if (uu_lb >= boundInf()) {
                return boundInf();
            }

            if (options.use_uu_unsupported) {
                t_part.restart();
                ui unsupported_lb = computeUUUnsupportedLowerBound(unmatched_vertices);
                solver.stats.lb_uu_unsupported_time += t_part.elapsed();
                if (unsupported_lb >= boundInf()) {
                    return boundInf();
                }
                uu_lb = std::max(uu_lb, unsupported_lb);
            }

            return addBound(addBound(current_miss, sf_lb), uu_lb);
#endif
        }

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

        ui computeSFMinLowerBound(const vector<ui> &frontier_vertices,
            const vector<vector<ui>> &anchors_by_frontier_idx) const
        {
            ui lb = 0;
            for (ui i = 0; i < frontier_vertices.size(); ++i) {
                ui u = frontier_vertices[i];
                const vector<ui> &anchors = anchors_by_frontier_idx[i];
                ui best_anchor_cut = boundInf();

                for (ui v : solver.candidates[u]) {
                    if (!isCandidateAvailable(u, v)) {
                        continue;
                    }

                    ui alpha = computeAnchorCutDelta(v, anchors);
                    best_anchor_cut = std::min(best_anchor_cut, alpha);
                }

                if (best_anchor_cut >= boundInf()) {
                    return boundInf();
                }

                lb = addBound(lb, best_anchor_cut);
                if (lb > solver.threshold) {
                    return lb;
                }
            }
            return lb;
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

        ui computeSFAssignmentLowerBound(const vector<ui> &frontier_vertices,
            const vector<vector<ui>> &anchors_by_frontier_idx)
        {
            ui frontier_size = (ui)frontier_vertices.size();
            if (frontier_size == 0) {
                return 0;
            }
            if (frontier_size == 1) {
                return computeSFMinLowerBound(frontier_vertices, anchors_by_frontier_idx);
            }

            vector<vector<vector<ui>>> exact_adj(frontier_size, vector<vector<ui>>(solver.qn + 1));
            vector<ui> right_vertices;
            vector<char> has_left_candidate(frontier_size, 0);

            if (++solver.lb_data_frontier_token == 0) {
                std::fill(solver.lb_data_frontier_mark.begin(), solver.lb_data_frontier_mark.end(), 0);
                solver.lb_data_frontier_token = 1;
            }

            for (ui i = 0; i < frontier_size; ++i) {
                ui u = frontier_vertices[i];
                const vector<ui> &anchors = anchors_by_frontier_idx[i];

                for (ui v : solver.candidates[u]) {
                    if (!isCandidateAvailable(u, v)) {
                        continue;
                    }

                    ui alpha = computeAnchorCutDelta(v, anchors);
                    assert(alpha <= solver.qn);
                    exact_adj[i][alpha].push_back(v);
                    has_left_candidate[i] = 1;
                    if (solver.lb_data_frontier_mark[v] != solver.lb_data_frontier_token) {
                        solver.lb_data_frontier_mark[v] = solver.lb_data_frontier_token;
                        right_vertices.push_back(v);
                    }
                }

                if (!has_left_candidate[i]) {
                    return boundInf();
                }
            }

            vector<vector<ui>> cumulative_adj(frontier_size);
            ui assignment_lb = 0;

            for (ui k = 0; k <= solver.qn; ++k) {
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
                if (k == solver.qn) {
                    return deficit == 0 ? assignment_lb : boundInf();
                }

                assignment_lb = addBound(assignment_lb, deficit);
                if (assignment_lb > solver.threshold) {
                    return assignment_lb;
                }
            }

            return assignment_lb;
        }

        void collectUUQueryCounts(ui u, vector<ui> &uu_neighbors,
            vector<ui> &query_counts) const
        {
            uu_neighbors.clear();
            std::fill(query_counts.begin(), query_counts.end(), 0);

            for (ui w : solver.q_neighbors[u]) {
                if (solver.in_Mq[w] || solver.is_excluded[u][w]) {
                    continue;
                }

                uu_neighbors.push_back(w);
                LabelID label = solver.query_graph->getVertexLabel(w);
                if (label >= 0 && (ui)label < solver.label_count) {
                    query_counts[(ui)label]++;
                }
            }
        }

        bool canSupportSomeUUNeighbor(ui data_vertex,
            const vector<ui> &uu_neighbors) const
        {
            for (ui w : uu_neighbors) {
                if (solver.candidates[w].contains(data_vertex) &&
                    isCandidateAvailable(w, data_vertex)) {
                    return true;
                }
            }
            return false;
        }

        ui computeUUPairDeficit(ui u, ui v, const vector<ui> &uu_neighbors,
            const vector<ui> &query_counts, const LowerBoundOptions &options,
            vector<ui> &data_counts) const
        {
            (void)u;
            std::fill(data_counts.begin(), data_counts.end(), 0);

            ui deg = 0;
            const ui *nbrs = solver.data_graph->getVertexNeighbors(v, deg);
            for (ui i = 0; i < deg; ++i) {
                ui t = nbrs[i];
                if (t == v || solver.mapped_g[t] != -1) {
                    continue;
                }
                if (options.restrict_data_neighbors_by_candidates &&
                    !canSupportSomeUUNeighbor(t, uu_neighbors)) {
                    continue;
                }

                LabelID label = solver.data_graph->getVertexLabel(t);
                if (label >= 0 && (ui)label < solver.label_count) {
                    data_counts[(ui)label]++;
                }
            }

            ui deficit = 0;
            for (ui label = 0; label < solver.label_count; ++label) {
                if (query_counts[label] > data_counts[label]) {
                    deficit += query_counts[label] - data_counts[label];
                }
            }
            return deficit;
        }

        ui computeUULabelLowerBound(const vector<ui> &unmatched_vertices,
            const LowerBoundOptions &options) const
        {
            vector<ui> uu_neighbors;
            vector<ui> query_counts(solver.label_count, 0);
            vector<ui> data_counts(solver.label_count, 0);
            ui sum_def = 0;

            for (ui u : unmatched_vertices) {
                collectUUQueryCounts(u, uu_neighbors, query_counts);

                ui best_deficit = boundInf();
                for (ui v : solver.candidates[u]) {
                    if (!isCandidateAvailable(u, v)) {
                        continue;
                    }

                    ui deficit = computeUUPairDeficit(
                        u, v, uu_neighbors, query_counts, options, data_counts);
                    best_deficit = std::min(best_deficit, deficit);
                }

                if (best_deficit >= boundInf()) {
                    return boundInf();
                }

                sum_def = addBound(sum_def, best_deficit);
                if ((sum_def + 1) / 2 > solver.threshold) {
                    return (sum_def + 1) / 2;
                }
            }

            return (sum_def + 1) / 2;
        }

        bool hasSupportedUUCandidatePair(ui u, ui w) const
        {
            for (ui yu : solver.candidates[u]) {
                if (!isCandidateAvailable(u, yu)) {
                    continue;
                }

                ui deg = 0;
                const ui *nbrs = solver.data_graph->getVertexNeighbors(yu, deg);
                for (ui i = 0; i < deg; ++i) {
                    ui yw = nbrs[i];
                    if (yw == yu) {
                        continue;
                    }
                    if (solver.candidates[w].contains(yw) &&
                        isCandidateAvailable(w, yw)) {
                        return true;
                    }
                }
            }
            return false;
        }

        ui computeUUUnsupportedLowerBound(const vector<ui> &unmatched_vertices) const
        {
            (void)unmatched_vertices;
            ui lb = 0;
            for (ui u = 0; u < solver.qn; ++u) {
                if (solver.in_Mq[u]) {
                    continue;
                }

                for (ui w : solver.q_neighbors[u]) {
                    if (u >= w || solver.in_Mq[w] || solver.is_excluded[u][w]) {
                        continue;
                    }

                    if (!hasSupportedUUCandidatePair(u, w)) {
                        lb++;
                        if (lb > solver.threshold) {
                            return lb;
                        }
                    }
                }
            }
            return lb;
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
    void dfs(ui cost, int u_new, const FrontierState *parent_state,
        const LowerBoundState *parent_lb_state)
    {
        (void)parent_lb_state;
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

        Timer t_frontier;
        BranchSelector branch_selector(*this);
        FrontierState current_state = branch_selector.buildFrontierState(parent_state, u_new); // sorted
        const vector<ui> &U_frontier = current_state.component_frontier;
        assert(!U_frontier.empty());
        stats.frontier_time += t_frontier.elapsed();

        const LowerBoundState *child_lb_state = nullptr;
#if defined(LOWER_BOUND) && !defined(CDE_LB_LIGHTWEIGHT_SPOKE)
        assert(threshold >= cost);
        LowerBoundState current_lb_state;
        LowerBoundState *current_lb_state_ptr = nullptr;
#ifdef CDE_LB_COMPONENT_MWPM_CACHE
        current_lb_state_ptr = &current_lb_state;
        child_lb_state = current_lb_state_ptr;
#endif
        if (threshold - cost <= (ui)LOWER_BOUND_MISSING_GAP) {
            Timer t_lb;
            if (LowerBoundPruner(*this).shouldPrune(
                cost, U_frontier, parent_lb_state, current_lb_state_ptr)) {
                stats.lb_time += t_lb.elapsed();
                stats.prun_calls++;
                return;
            }
            stats.lb_time += t_lb.elapsed();
        }
#endif

        Timer t_branch;
        long long child_dfs_time = 0;
        long long local_lb_time = 0;
        vector<pair<ui, ui>> local_X;       // Records changes to is_excluded
        vector<pair<ui, ui>> local_x_cand;  // Records changes to x_cand
        vector<ui> cand_v_list;
#if defined(LOWER_BOUND) && defined(CDE_LB_LIGHTWEIGHT_SPOKE)
        LightweightSpokeLowerBound light_lb(*this);
#endif

        ui current_cost = cost;

        while (current_cost <= threshold) {
            vector<ActiveEdge> top_edges;
            ui top_k = threshold - current_cost + 1;
            if (!collectTopActiveEdges(U_frontier, top_k, top_edges)) {
                break;
            }

            bool made_progress = false;
            for (const ActiveEdge &edge : top_edges) {
                if (current_cost > threshold) {
                    break;
                }

                ui u = edge.u;
                ui ua = edge.anchor;
                if (u >= qn || in_Mq[u] || !in_Mq[ua] ||
                    is_excluded[u][ua] || frontier_pos[u] == -1) {
                    continue;
                }

                collectEdgeCandidates(u, ua, cand_v_list);

                // Include branch: use this frontier-anchor edge to add exactly one
                // new query vertex in the recursive child state.
                for (ui v : cand_v_list) {
                    assert(mapped_g[v] == -1);
                    assert(!x_cand[u].contains(v));

                    ui delta = computeLiveAnchorDelta(u, v, ua);
                    if (current_cost + delta > threshold) continue;

#if defined(LOWER_BOUND) && defined(CDE_LB_LIGHTWEIGHT_SPOKE)
                    {
                        Timer t_lb;
                        bool prune_by_light_lb = light_lb.shouldPrune(
                            current_cost + delta, u, v);
                        long long elapsed = t_lb.elapsed();
                        stats.lb_time += elapsed;
                        stats.lb_light_spoke_time += elapsed;
                        local_lb_time += elapsed;
                        if (prune_by_light_lb) {
                            stats.prun_calls++;
                            continue;
                        }
                    }
#endif

#ifndef NDEBUG
                    recordBranchOrderDebug(u);
#endif
                    updateResidualLabelCounts(query_graph, u, mapped_q,
                        Lq_residual_counts, Lq_residual_degrees, false);
                    updateResidualLabelCounts(data_graph, v, mapped_g,
                        Lg_residual_counts, Lg_residual_degrees, false);
                    mapped_q[u] = (int)v;
                    mapped_g[v] = (int)u;
                    in_Mq[u] = 1;
                    part_M.push_back({ u, v });

                    updateFrontier(u, true);

                    Timer t_child;
                    dfs(current_cost + delta, (int)u, &current_state, child_lb_state);
                    child_dfs_time += t_child.elapsed();

                    updateFrontier(u, false);

                    part_M.pop_back();
                    in_Mq[u] = 0;
                    mapped_g[v] = -1;
                    mapped_q[u] = -1;
                    updateResidualLabelCounts(query_graph, u, mapped_q,
                        Lq_residual_counts, Lq_residual_degrees, true);
                    updateResidualLabelCounts(data_graph, v, mapped_g,
                        Lg_residual_counts, Lg_residual_degrees, true);
                }

                // Exclude branch: keep this decision in the current frame, then
                // continue to the next selected active edge without a recursive
                // call. Every recursive child above still adds a vertex.
                current_cost++;
                is_excluded[u][ua] = 1;
                is_excluded[ua][u] = 1;
                anchor_count[u]--;
                local_X.push_back({ u, ua });
                updateFrontierStatus(u);
                made_progress = true;

                if (current_cost > threshold) {
                    break;
                }

#ifndef NDEBUG
                recordBranchOrderDebug(u);
#endif
                for (ui v : cand_v_list) {
                    if (!x_cand[u].contains(v)) {
                        insertXCand(u, v);
                        local_x_cand.push_back({ u, v });
                    }
                }
            }

            if (!made_progress) {
                break;
            }
        }

        {
            long long branch_elapsed = t_branch.elapsed();
            long long excluded_time = child_dfs_time + local_lb_time;
            stats.branch_time += branch_elapsed > excluded_time
                ? branch_elapsed - excluded_time : 0;
        }

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
            removeXCand(it->first, it->second);
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
