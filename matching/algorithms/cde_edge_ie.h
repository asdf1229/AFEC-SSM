#ifndef MATCHING_ALGORITHMS_CDE_EDGE_IE_H_
#define MATCHING_ALGORITHMS_CDE_EDGE_IE_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
using namespace std;

// Anchor-support scoring defaults to the current maintained-score path.
// Define CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT to use the older on-demand path:
// every score reads the current candidates and recomputes support.
#if defined(CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT) && defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT)
#error "Choose only one cde_edge_ie anchor-support mode."
#endif
#ifndef CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT
#define CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
#endif

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
#ifndef NDEBUG
        printf("CDE-Edge-IE label_count=%u (query=%u, data=%u)\n",
            label_count, query_graph->getLabelsCount(), data_graph->getLabelsCount());
#endif

        if (qn == 0 || gn == 0) {
            stats.init_time = t_init.elapsed();
            return false;
        }

        resetState();

        initQueryGraphIndex();

        Timer t_filter;
        bool res = runCandidateFiltering();
        stats.filter_time = t_filter.elapsed();
        stats.init_time = t_init.elapsed();
        if (!res) return false;
        return true;
    }

    void match(vector<vector<pair<ui, ui>>> &results)
    {
        Timer t_search;
        t_search.restart();

        results_ptr = &results;
        results_ptr->clear();
        initDfsBuffer();

        BranchSelector branch_selector(*this);
        ui root = branch_selector.selectInitialRoot();
#ifndef NDEBUG
        // if (kDebugInitialRoot >= 0) {
        //     assert((ui)kDebugInitialRoot < qn);
        //     root = (ui)kDebugInitialRoot;
        // }
        printf("Selected initial root: u=%u with %zu candidates\n", root, (size_t)candidates[root].size());
#endif

        for (ui v0 : candidates[root]) {
#ifndef NDEBUG
            recordBranchOrderDebug(root);
#endif

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_MAPPED_VERTEX_SUPPORT)
            adjustSupportForMappedDataVertex(v0, false);
#endif
            mapped_q[root] = (int)v0;
            mapped_g[v0] = (int)root;
            part_M.push_back({ root, v0 });

            updateFrontier(root, true);

            dfs(0);

            part_M.pop_back();

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_MAPPED_VERTEX_SUPPORT)
            mapped_q[root] = -1;
            updateFrontier(root, false);
            mapped_g[v0] = -1;
            adjustSupportForMappedDataVertex(v0, true);
#else
            mapped_g[v0] = -1;
            mapped_q[root] = -1;
            updateFrontier(root, false);
#endif
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
        long long lb_light_spoke_time = 0;
        long long frontier_time = 0; // building and ordering U_frontier
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

#if defined(CDE_LB_LIGHTWEIGHT_SPOKE) && !defined(NDEBUG)
        long long lb_accounted_time = stats.lb_light_spoke_time;
        long long lb_other_time = stats.lb_time > lb_accounted_time
            ? stats.lb_time - lb_accounted_time : 0;
#endif
        long long search_accounted_time = stats.frontier_time + stats.branch_time;
#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
        search_accounted_time += stats.lb_time;
#endif
        long long search_other_time = stats.dfs_time > search_accounted_time
            ? stats.dfs_time - search_accounted_time : 0;

        printf("\n--- CDE-Edge-IE Time Analysis ---\n");
#ifdef NDEBUG
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms\n", stats.filter_time / 1000.0);
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
#endif
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("  - Branch Time:     %.4lf ms (%.2f%% of Search)\n", stats.branch_time / 1000.0, pct(stats.branch_time, stats.dfs_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.dfs_time));
#else
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
        printf("    - NLF:           %.4lf ms (%.2f%% of Filter)\n", stats.filter_nlf_time / 1000.0, pct(stats.filter_nlf_time, stats.filter_time));
        printf("    - Bridge:        %.4lf ms (%.2f%% of Filter)\n", stats.filter_bridge_time / 1000.0, pct(stats.filter_bridge_time, stats.filter_time));
#ifdef ENABLE_SPOKE_FILTERING
        printf("    - Spoke:         %.4lf ms (%.2f%% of Filter)\n", stats.filter_spoke_time / 1000.0, pct(stats.filter_spoke_time, stats.filter_time));
#endif
#ifdef ENABLE_ONEHOP_FILTERING
        printf("    - OneHop:        %.4lf ms (%.2f%% of Filter)\n", stats.filter_onehop_time / 1000.0, pct(stats.filter_onehop_time, stats.filter_time));
#endif
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
        printf("  - LowerBound Time: %.4lf ms (%.2f%% of Search)\n", stats.lb_time / 1000.0, pct(stats.lb_time, stats.dfs_time));
        printf("    - Light Spoke:   %.4lf ms (%.2f%% of LB)\n", stats.lb_light_spoke_time / 1000.0, pct(stats.lb_light_spoke_time, stats.lb_time));
        printf("    - Other:         %.4lf ms (%.2f%% of LB)\n", lb_other_time / 1000.0, pct(lb_other_time, stats.lb_time));
#endif
        printf("  - Frontier Time:   %.4lf ms (%.2f%% of Search)\n", stats.frontier_time / 1000.0, pct(stats.frontier_time, stats.dfs_time));
        printf("    - Select Best:   %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_select_time / 1000.0, pct(stats.frontier_select_time, stats.frontier_time));
        printf("    - Component:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_component_time / 1000.0, pct(stats.frontier_component_time, stats.frontier_time));
#ifdef ENABLE_FRONTIER_ORDERING
        printf("    - Sort Hook:     %.4lf ms (%.2f%% of Frontier)\n", stats.frontier_sort_time / 1000.0, pct(stats.frontier_sort_time, stats.frontier_time));
#endif
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

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
    struct SupportDelta {
        ui u = 0;
        ui anchor = 0;
        ui amount = 0;
    };
#endif

    struct EdgeScoreCache {
        vector<vector<ActiveEdge>> active_edges_by_vertex;
        vector<char> active_edges_cached;
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
    vector<ui> q_degree;
    vector<vector<char>> q_neighbor_is_bridge;

    vector<MyBitset> candidates;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<vector<char>> excluded_edges;
    vector<MyBitset> excluded_cands;
    vector<pair<ui, ui>> part_M;

    vector<ui> anchor_count;
#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
    vector<vector<ui>> support;
#endif
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
        ui live_anchor_count = 0;
        ui query_degree = 0;

        bool live_anchors_ready = false;
        bool query_degree_ready = false;
    };

    // Frontier state for the selected unmatched query component in the current DFS step.
    struct FrontierState {
        // All unmatched query vertices in the selected connected component
        vector<ui> component_vertices;
        // Active frontier vertices within the selected component
        vector<ui> component_frontier;
        ActiveEdge selected_global_best_edge;
        EdgeScoreCache edge_score_cache;
        // Lazy score cache computed while choosing/sorting frontier vertices in this state.
        vector<FrontierScore> frontier_score_cache;
        vector<char> frontier_score_cached;
    };
    vector<ui> frontier_visit;
    ui frontier_token;

#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
    vector<int> lb_match_right;
    vector<int> lb_match_left;
    vector<char> lb_light_spoke_vis;
    vector<vector<ui>> lb_light_spoke_adj;
    vector<char> lb_light_spoke_is_bridge;
#endif

    void resetState()
    {
        mapped_q.assign(qn, -1);
        mapped_g.assign(gn, -1);
        excluded_edges.assign(qn, vector<char>(qn, 0));
        part_M.clear();
        part_M.reserve(qn);

        candidates.clear();
        candidates.assign(qn, MyBitset(gn));

        excluded_cands.clear();
        excluded_cands.assign(qn, MyBitset(gn));

        q_matrix.clear();
        q_neighbors.clear();
        q_degree.clear();
        q_neighbor_is_bridge.clear();

        anchor_count.assign(qn, 0);
#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
        support.assign(qn, vector<ui>(qn, 0));
#endif
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
#ifndef NDEBUG
        branch_order_prefix_counts.clear();
        branch_order_prefix_reach_counts.clear();
        branch_order_prefix_answer_counts.clear();
#endif
        frontier_visit.assign(qn, 0);
        frontier_token = 1;

#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
        lb_match_right.assign(max_g_deg, -1);
        lb_match_left.assign(qn, -1);
        lb_light_spoke_vis.assign(max_g_deg, 0);
        lb_light_spoke_adj.assign(qn, vector<ui>());
        lb_light_spoke_is_bridge.assign(qn, 0);
#endif

        stats = TimeStats();
    }

    void initQueryGraphIndex()
    {
        q_matrix.assign(qn, vector<char>(qn, 0));
        q_neighbors.assign(qn, vector<ui>());

        q_degree.assign(qn, 0);

        for (ui u = 0; u < qn; ++u) {
            ui deg = 0;
            const ui *nbrs = query_graph->getVertexNeighbors(u, deg);

            q_neighbors[u].reserve(deg);
            for (ui i = 0; i < deg; ++i) {
                ui u1 = nbrs[i];
                q_matrix[u][u1] = 1;
                q_neighbors[u].push_back(u1);
            }
            q_degree[u] = (ui)q_neighbors[u].size();
        }
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

    // ========================================================================
    // Filtering
    // ========================================================================
    struct CandidateFilter {
        MatchingSolver &solver;
        struct LabelCount {
            LabelID label = 0;
            ui count = 0;
        };
        struct BridgeLabelCount {
            LabelID label = 0;
            ui bridge_count = 0;
            ui non_bridge_count = 0;
        };
        struct OneHopCandKey {
            ui u = 0;
            ui v = 0;

            bool operator==(const OneHopCandKey &other) const {
                return u == other.u && v == other.v;
            }
        };
        struct OneHopCandKeyHash {
            size_t operator()(const OneHopCandKey &key) const {
                return (uint64_t(key.u) << 32) | uint64_t(key.v);
            }
        };
        vector<vector<char>> q_bridge_matrix;
        vector<vector<ui>> q_bridge_adj;

        // NLF filtering
        vector<vector<LabelCount>> Lg_counts;
        vector<ui> Lg_degrees;
        vector<vector<BridgeLabelCount>> Lq_bridge_label_counts;
        vector<ui> Lq_bridge_degrees, Lq_non_bridge_degrees;

        // spoke filtering
        vector<vector<ui>>  spoke_matrix;
        vector<int>         spoke_match_right;
        vector<int>         spoke_match_left;
        vector<char>        spoke_vis;
        vector<char>        spoke_is_bridge;

        // One-hop filtering
        unordered_set<OneHopCandKey, OneHopCandKeyHash> onehop_cands;
        MyBitset            onehop_vis;
        vector<ui>          onehop_inner_edges;
        vector<vector<ui>>  onehop_inner_degree;

        explicit CandidateFilter(MatchingSolver &solver)
            : solver(solver),
            spoke_matrix(solver.qn, vector<ui>()),
            spoke_match_right(solver.max_g_deg, -1),
            spoke_match_left(solver.qn, -1),
            spoke_vis(solver.max_g_deg, 0),
            spoke_is_bridge(solver.qn, 0),
            onehop_vis((int)solver.gn)
        {}

        bool run()
        {
            Timer t;
            initBridge();
            solver.stats.filter_bridge_time += t.elapsed();

            // --- NLF Filtering ---
            t.restart();
            bool ok = filterNLF();
            solver.stats.filter_nlf_time += t.elapsed();
            if (!ok) return false;
            t.restart();
            ok = filterBridge();
            solver.stats.filter_bridge_time += t.elapsed();
            if (!ok) return false;

            // --- Spoke Filtering ---
#ifdef ENABLE_SPOKE_FILTERING
            t.restart();
            ok = filterSpoke();
            solver.stats.filter_spoke_time += t.elapsed();
            if (!ok) return false;
            t.restart();
            ok = filterBridge();
            solver.stats.filter_bridge_time += t.elapsed();
            if (!ok) return false;
#endif

            // --- One-Hop Filtering ---
#ifdef ENABLE_ONEHOP_FILTERING
            t.restart();
            ok = filterOneHop();
            solver.stats.filter_onehop_time += t.elapsed();
            if (!ok) return false;
            t.restart();
            ok = filterBridge();
            solver.stats.filter_bridge_time += t.elapsed();
            if (!ok) return false;
#endif

            solver.stats.filter_candidate_count = 0;
            for (ui u = 0; u < solver.qn; ++u) {
                solver.stats.filter_candidate_count += (ui)solver.candidates[u].size();
            }

#ifdef ENABLE_CAND_STATS
            printCandStats();
#endif
            return true;
        }

    private:
        void initBridge()
        {
            // find query bridges
            findQueryBridges();
            // build query bridge neighbor index
            buildBridgeNeighbors();
        }

        void tarjan(ui u, ui fa, vector<int> &dfn, vector<int> &low, int &tim)
        {
            dfn[u] = low[u] = ++tim;

            for (ui v : solver.q_neighbors[u]) {
                if (dfn[v] == 0) {
                    tarjan(v, u, dfn, low, tim);
                    low[u] = std::min(low[u], low[v]);
                    if (low[v] > dfn[u]) {
                        if (!q_bridge_matrix[u][v]) {
                            q_bridge_matrix[u][v] = 1;
                            q_bridge_matrix[v][u] = 1;
                            q_bridge_adj[u].push_back(v);
                            q_bridge_adj[v].push_back(u);
                        }
                    }
                }
                else if (v != fa) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
        }

        void findQueryBridges()
        {
            q_bridge_matrix.assign(solver.qn, vector<char>(solver.qn, 0));
            q_bridge_adj.assign(solver.qn, vector<ui>());

            vector<int> dfn(solver.qn, 0);
            vector<int> low(solver.qn, 0);
            int tim = 0;

            for (ui u = 0; u < solver.qn; ++u) {
                if (dfn[u] == 0) {
                    tarjan(u, solver.qn, dfn, low, tim);
                }
            }
        }

        void buildBridgeNeighbors()
        {
            solver.q_neighbor_is_bridge.assign(solver.qn, vector<char>());

            for (ui u = 0; u < solver.qn; ++u) {
                const vector<ui> &u_neighbors = solver.q_neighbors[u];
                solver.q_neighbor_is_bridge[u].assign(solver.q_degree[u], 0);
                for (ui i = 0; i < solver.q_degree[u]; ++i) {
                    solver.q_neighbor_is_bridge[u][i] = q_bridge_matrix[u][u_neighbors[i]];
                }
            }
        }

        void buildQueryLabelCounts()
        {
            Lq_bridge_label_counts.assign(solver.qn, vector<BridgeLabelCount>());
            Lq_bridge_degrees.assign(solver.qn, 0);
            Lq_non_bridge_degrees.assign(solver.qn, 0);

            for (ui u = 0; u < solver.qn; ++u) {
                vector<pair<LabelID, char>> labels;
                labels.reserve(solver.q_degree[u]);

                for (ui i = 0; i < solver.q_degree[u]; ++i) {
                    ui u1 = solver.q_neighbors[u][i];
                    LabelID label = solver.query_graph->getVertexLabel(u1);
                    assert(label >= 0 && (ui)label < solver.label_count);

                    if (solver.q_neighbor_is_bridge[u][i]) {
                        Lq_bridge_degrees[u]++;
                        labels.push_back({ label, 1 });
                    }
                    else {
                        Lq_non_bridge_degrees[u]++;
                        labels.push_back({ label, 0 });
                    }
                }

                std::sort(labels.begin(), labels.end(),
                    [](const pair<LabelID, char> &lhs, const pair<LabelID, char> &rhs) {
                        return lhs.first < rhs.first;
                    });

                for (const auto &entry : labels) {
                    if (Lq_bridge_label_counts[u].empty() ||
                        Lq_bridge_label_counts[u].back().label != entry.first) {
                        Lq_bridge_label_counts[u].push_back({ entry.first, 0, 0 });
                    }

                    BridgeLabelCount &count = Lq_bridge_label_counts[u].back();
                    if (entry.second) {
                        count.bridge_count++;
                    }
                    else {
                        count.non_bridge_count++;
                    }
                }
            }
        }

        void buildDataLabelCounts()
        {
            Lg_counts.assign(solver.gn, vector<LabelCount>());
            Lg_degrees.assign(solver.gn, 0);

            for (ui u = 0; u < solver.gn; ++u) {
                ui deg = 0;
                const ui *neighbors = solver.data_graph->getVertexNeighbors(u, deg);
                Lg_degrees[u] = deg;
                vector<LabelID> labels;
                labels.reserve(deg);

                for (ui i = 0; i < deg; ++i) {
                    LabelID label = solver.data_graph->getVertexLabel(neighbors[i]);
                    assert(label >= 0 && (ui)label < solver.label_count);
                    labels.push_back(label);
                }

                std::sort(labels.begin(), labels.end());
                for (LabelID label : labels) {
                    if (Lg_counts[u].empty() || Lg_counts[u].back().label != label) {
                        Lg_counts[u].push_back({ label, 1 });
                    }
                    else {
                        Lg_counts[u].back().count++;
                    }
                }
            }
        }

        ui computeNLF(ui u, ui v) const
        {
            ui diff = 0;
            const vector<BridgeLabelCount> &query_counts = Lq_bridge_label_counts[u];
            const vector<LabelCount> &data_counts = Lg_counts[v];
            size_t data_idx = 0;

            for (const auto &query_count : query_counts) {
                while (data_idx < data_counts.size() && data_counts[data_idx].label < query_count.label) {
                    data_idx++;
                }

                ui data_count = 0;
                if (data_idx < data_counts.size() && data_counts[data_idx].label == query_count.label) {
                    data_count = data_counts[data_idx].count;
                }

                ui bridge_need = query_count.bridge_count;
                if (bridge_need > data_count) {
                    return solver.threshold + 1;
                }

                ui non_bridge_need = query_count.non_bridge_count;
                ui non_bridge_available = data_count - bridge_need;
                if (non_bridge_need > non_bridge_available) {
                    diff += (non_bridge_need - non_bridge_available);
                }
                if (diff > solver.threshold) {
                    return diff;
                }
            }
            return diff;
        }

        bool filterNLF()
        {
            buildQueryLabelCounts();
            buildDataLabelCounts();
            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = solver.query_graph->getVertexLabel(u);
                for (ui v = 0; v < solver.gn; ++v) {
                    if (lu != solver.data_graph->getVertexLabel(v)) continue;
                    if (Lq_bridge_degrees[u] > Lg_degrees[v]) continue;
                    if (Lq_non_bridge_degrees[u] + Lq_bridge_degrees[u] > Lg_degrees[v] + solver.threshold) continue;
                    if (computeNLF(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        bool hasBridgeCandidate(ui v, ui other_u) const
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

        bool pruneBridge(ui u, ui other_u, bool &changed)
        {
            vector<ui> to_remove;
            for (ui v : solver.candidates[u]) {
                if (!hasBridgeCandidate(v, other_u)) {
                    to_remove.push_back(v);
                }
            }

            if (to_remove.empty()) return true;

            changed = true;
            for (ui v : to_remove) {
                solver.candidates[u].remove(v);
            }
            return !solver.candidates[u].empty();
        }

        bool filterBridge()
        {
            queue<ui> q;
            vector<char> in_q(solver.qn, 0);
            for (ui u = 0; u < solver.qn; ++u) {
                if (q_bridge_adj[u].empty()) continue;
                q.push(u);
                in_q[u] = 1;
            }

            while (!q.empty()) {
                ui u = q.front(); q.pop();
                in_q[u] = 0;

                bool changed = false;
                for (ui other_u : q_bridge_adj[u]) if (!pruneBridge(u, other_u, changed)) {
                    return false;
                }
                if (!changed) continue;

                for (ui nbr_u : q_bridge_adj[u]) {
                    if (!in_q[nbr_u]) {
                        q.push(nbr_u);
                        in_q[nbr_u] = 1;
                    }
                }
            }
            return true;
        }

        bool dfsMatchSpoke(ui left_idx, const vector<vector<ui>> &adj)
        {
            for (ui right_idx : adj[left_idx]) {
                if (spoke_vis[right_idx]) continue;
                spoke_vis[right_idx] = 1;
                if (spoke_match_right[right_idx] < 0 ||
                    dfsMatchSpoke((ui)spoke_match_right[right_idx], adj)) {
                    spoke_match_right[right_idx] = (int)left_idx;
                    spoke_match_left[left_idx] = (int)right_idx;
                    return true;
                }
            }
            return false;
        }

        bool maxMatchSpokeWithBridge(const vector<vector<ui>> &adj, ui left_size, ui right_size, ui &missing_non_bridge)
        {
            std::fill(spoke_match_right.begin(), spoke_match_right.begin() + right_size, -1);
            std::fill(spoke_match_left.begin(), spoke_match_left.begin() + left_size, -1);

            ui non_bridge_count = 0;
            for (ui i = 0; i < left_size; ++i) {
                if (spoke_is_bridge[i]) {
                    std::fill(spoke_vis.begin(), spoke_vis.begin() + right_size, 0);
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

                std::fill(spoke_vis.begin(), spoke_vis.begin() + right_size, 0);
                if (dfsMatchSpoke(i, adj)) {
                    non_bridge_match++;
                }
            }

            missing_non_bridge = non_bridge_count - non_bridge_match;
            return true;
        }

        // If u is matched to v, compute the minimum number of spoke edges
        // from u to its neighbors that cannot be supported by neighbors of v.
        ui computeLBSpoke(ui u, ui v)
        {
            // build bipartite graph
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            ui deg_u = solver.q_degree[u];
            ui deg_v = 0;
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            for (ui i = 0; i < deg_u; ++i) {
                spoke_matrix[i].clear();
                ui u1 = u_neighbors[i];
                spoke_is_bridge[i] = solver.q_neighbor_is_bridge[u][i];
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        spoke_matrix[i].push_back(j);
                    }
                }
            }

            // matching with all bridge spokes matched
            ui missing_non_bridge = 0;
            if (!maxMatchSpokeWithBridge(spoke_matrix, deg_u, deg_v, missing_non_bridge)) {
                return solver.threshold + 1;
            }
            return missing_non_bridge;
        }

        bool filterSpoke()
        {
            queue<ui> q;
            vector<char> in_q(solver.qn, 1);
            for (ui u = 0; u < solver.qn; ++u) q.push(u);

            while (!q.empty()) {
                ui u = q.front(); q.pop();
                in_q[u] = 0;

                if(Lq_bridge_degrees[u] == 0 && solver.q_degree[u] <= solver.threshold) continue;

                vector<ui> to_remove;
                for (ui v : solver.candidates[u]) {
                    ui missing_edges = computeLBSpoke(u, v);
                    if(missing_edges > solver.threshold) {
                        to_remove.push_back(v);
                    }
#ifdef ENABLE_ONEHOP_FILTERING
                    else if (solver.threshold - missing_edges <= (ui)ONEHOP_FILTER_MISSING_GAP) {
                        onehop_cands.insert({ u, v });
                    }
#endif
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

        void initOneHop()
        {
            onehop_inner_edges.assign(solver.qn, 0);
            onehop_inner_degree.assign(solver.qn, vector<ui>());

            for (ui u = 0; u < solver.qn; ++u) {
                ui deg_u = solver.q_degree[u];
                onehop_inner_degree[u].assign(deg_u, 0);

                const vector<ui> &u_neighbors = solver.q_neighbors[u];
                for (ui i = 0; i < deg_u; ++i) {
                    for (ui j = i + 1; j < deg_u; ++j) {
                        if (!solver.q_matrix[u_neighbors[i]][u_neighbors[j]]) continue;

                        onehop_inner_edges[u]++;
                        onehop_inner_degree[u][i]++;
                        onehop_inner_degree[u][j]++;
                    }
                }
            }
        }

        // Build the DFS order for one-hop matching.
        // Fewer candidates first.
        // More inner edges first if tied.
        vector<ui> buildOneHopOrder(ui u, const vector<vector<ui>> &cand, const vector<char> &spoke_bridge) const
        {
            ui deg_u = solver.q_degree[u];
            vector<ui> ord(deg_u);
            iota(ord.begin(), ord.end(), 0);

            sort(ord.begin(), ord.end(), [&](ui a, ui b) {
                if (spoke_bridge[a] != spoke_bridge[b]) {
                    return spoke_bridge[a] > spoke_bridge[b];
                }
                if (cand[a].size() != cand[b].size()) {
                    return cand[a].size() < cand[b].size();
                }
                return onehop_inner_degree[u][a] > onehop_inner_degree[u][b];
                });
            return ord;
        }

        ui computeRemainLBOneHop(ui pos, const vector<ui> &ord, const vector<vector<ui>> &cand, const vector<char> &spoke_bridge) const
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
            const vector<ui> &ord, const vector<ui> &u_neighbors,
            const vector<vector<ui>> &cand, const vector<char> &spoke_bridge)
        {
            if (cost > solver.threshold) return false;
            if (pos == ord.size()) return true;

            ui rem_lb = computeRemainLBOneHop(pos, ord, cand, spoke_bridge);
            if (cost + rem_lb > solver.threshold) return false;

            ui i = ord[pos];
            ui u1 = u_neighbors[i];

            // branch 1: skip u1
            bool can_skip_spoke = (spoke_bridge[i] == 0);
            if (can_skip_spoke) {
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
                        if (q_bridge_matrix[u1][u2]) {
                            missing_bridge = true;
                            break;
                        }
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
            ui deg_u = solver.q_degree[u];
            ui deg_v = 0;
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);
            vector<vector<ui>> cand(deg_u);
            const vector<char> &spoke_bridge = solver.q_neighbor_is_bridge[u];

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
            vector<ui> ord = buildOneHopOrder(u, cand, spoke_bridge);

            // the current DFS state of u_neighbors[i]:
            // -2: unprocessed
            // -1: skipped / unmatched
            // >=0: matched to the corresponding data vertex
            vector<int> state(deg_u, -2);

            return dfsOneHop(0, state, 0, ord, u_neighbors, cand, spoke_bridge);
        }

        bool filterOneHop()
        {
            initOneHop();

            vector<pair<ui, ui>> to_remove;

            for (const auto &cand : onehop_cands) {
                ui u = cand.u, v = cand.v;
                if(!solver.candidates[u].contains(v)) continue;
                if (solver.q_degree[u] <= 1) continue;
                if (onehop_inner_edges[u] == 0) continue;
                if (!checkOneHop(u, v)) {
                    to_remove.push_back({ u, v });
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

#ifdef ENABLE_CAND_STATS
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
                ui deg_u = solver.q_degree[u];
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
#endif
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

    public:
        explicit BranchSelector(MatchingSolver &solver) : solver(solver) {}

        ui selectInitialRoot() const
        {
            ui root = 0;
            for (ui u = 1; u < solver.qn; ++u) {
                if (solver.candidates[u].size() < solver.candidates[root].size() ||
                    (solver.candidates[u].size() == solver.candidates[root].size() &&
                        solver.q_degree[u] > solver.q_degree[root])) {
                    root = u;
                }
            }
            return root;
        }

        void collectLiveAnchors(ui u, vector<ui> &anchors) const
        {
            anchors.clear();
            for (ui u1 : solver.q_neighbors[u]) {
                if (solver.mapped_q[u1] != -1 && !solver.excluded_edges[u][u1]) {
                    anchors.push_back(u1);
                }
            }
        }

        void collectCandVertices(ui u, ui anchor, vector<ui> &cand_v_list) const
        {
            cand_v_list.clear();
            enumerateAnchorCandidates(u, anchor, [&](ui v) {
                cand_v_list.push_back(v);
                });
        }

        bool collectTopActiveEdges(const vector<ui> &component_frontier, ui max_count,
            vector<ActiveEdge> &top_edges, EdgeScoreCache *edge_score_cache = nullptr) const
        {
            top_edges.clear();
            if (max_count == 0) {
                return false;
            }

            vector<ActiveEdge> uncached_edges;
            for (ui u : component_frontier) {
                if (u >= solver.qn) {
                    continue;
                }
                if (edge_score_cache != nullptr) {
                    const vector<ActiveEdge> &cached_edges = cachedActiveEdgesForVertex(u, *edge_score_cache);
                    top_edges.insert(top_edges.end(), cached_edges.begin(), cached_edges.end());
                }
                else {
                    collectActiveEdgesForVertex(u, uncached_edges);
                    top_edges.insert(top_edges.end(), uncached_edges.begin(), uncached_edges.end());
                }
            }

            if (top_edges.empty()) {
                return false;
            }

            auto better_edge = [&](const ActiveEdge &lhs, const ActiveEdge &rhs) {
                return isBetterActiveEdge(lhs, rhs);
                };
            if (top_edges.size() > max_count) {
                partial_sort(top_edges.begin(), top_edges.begin() + max_count,
                    top_edges.end(), better_edge);
                top_edges.resize(max_count);
            }
            else {
                sort(top_edges.begin(), top_edges.end(), better_edge);
            }
            return true;
        }

        FrontierState buildFrontierState(vector<ActiveEdge> &top_edges)
        {
            FrontierState state;
            score_state = &state;

            assert(!solver.active_frontier.empty());

            Timer t_select;
            bool found_edge = collectTopActiveEdges(solver.active_frontier, 1,
                top_edges, &state.edge_score_cache);
            assert(found_edge);
            (void)found_edge;
            state.selected_global_best_edge = top_edges.front();
            solver.stats.frontier_select_time += t_select.elapsed();

            Timer t_component;
            collectComponent(state.selected_global_best_edge.u, state);
            solver.stats.frontier_component_time += t_component.elapsed();

#ifdef ENABLE_FRONTIER_ORDERING
            Timer t_sort;
            sortFrontier(state.component_frontier);
            solver.stats.frontier_sort_time += t_sort.elapsed();
#endif
            return state;
        }

    private:
        template <typename Visitor>
        ui enumerateAnchorCandidates(ui u, ui anchor, Visitor visit) const
        {
            ui count = 0;
            ui deg = 0;
            const ui *nbrs = solver.data_graph->getVertexNeighbors((ui)solver.mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (!solver.candidates[u].contains(v)) continue;
                if (solver.mapped_g[v] != -1) continue;
                if (solver.excluded_cands[u].contains(v)) continue;
                count++;
                visit(v);
            }
            return count;
        }

        ui countLiveAnchors(ui u) const
        {
            ui count = 0;
            for (ui anchor : solver.q_neighbors[u]) {
                if (solver.mapped_q[anchor] != -1 && !solver.excluded_edges[u][anchor]) {
                    count++;
                }
            }
            return count;
        }

        ui countEdgeSupport(ui u, ui anchor) const
        {
            if (u >= solver.qn || anchor >= solver.qn || solver.mapped_q[anchor] == -1) {
                return 0;
            }

#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
            return solver.support[u][anchor];
#else
            return solver.computeEdgeSupport(u, anchor);
#endif
        }

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

        void collectActiveEdgesForVertex(ui u, vector<ActiveEdge> &edges) const
        {
            edges.clear();
            if (u >= solver.qn || solver.mapped_q[u] != -1 || solver.frontier_pos[u] == -1) {
                return;
            }

            ui live_anchor_count = countLiveAnchors(u);
            if (live_anchor_count == 0) {
                return;
            }

            for (ui anchor : solver.q_neighbors[u]) {
                if (solver.mapped_q[anchor] == -1 || solver.excluded_edges[u][anchor]) {
                    continue;
                }

                ActiveEdge edge;
                edge.u = u;
                edge.anchor = anchor;
                edge.anchor_support = countEdgeSupport(u, anchor);
                edge.live_anchor_count = live_anchor_count;
                edge.query_degree = solver.q_degree[u];
                edges.push_back(edge);
            }
        }

        const vector<ActiveEdge> &cachedActiveEdgesForVertex(ui u, EdgeScoreCache &cache) const
        {
            if (cache.active_edges_by_vertex.empty()) {
                cache.active_edges_by_vertex.resize(solver.qn);
                cache.active_edges_cached.assign(solver.qn, 0);
            }
            if (!cache.active_edges_cached[u]) {
                collectActiveEdgesForVertex(u, cache.active_edges_by_vertex[u]);
                cache.active_edges_cached[u] = 1;
            }
            return cache.active_edges_by_vertex[u];
        }

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

        void ensureScoreStorage(FrontierState &state) const
        {
            if (state.frontier_score_cache.empty()) {
                state.frontier_score_cache.resize(solver.qn);
                state.frontier_score_cached.assign(solver.qn, 0);
            }
        }

        // Returns the cached lazy frontier score for u in the current DFS state.
        FrontierScore &scoreFor(ui u) const
        {
            assert(u < solver.qn);
            assert(score_state != nullptr);
            ensureScoreStorage(*score_state);
            if (!score_state->frontier_score_cached[u]) {
                score_state->frontier_score_cache[u] = FrontierScore(u);
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

        ui cachedLiveAnchorCount(FrontierScore &score) const
        {
            collectAnchors(score);
            return score.live_anchor_count;
        }

        ui cachedQueryDegreeCount(FrontierScore &score) const
        {
            if (!score.query_degree_ready) {
                Timer t;
                score.query_degree = solver.q_degree[score.u];
                score.query_degree_ready = true;

                long long elapsed = t.elapsed();
                solver.stats.frontier_score_query_degree_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
            }

            return score.query_degree;
        }

        // More live anchors, higher query degree, smaller vertex id.
        // Edge support itself is scored per frontier edge in collectTopActiveEdges.
        bool isBetterFrontier(FrontierScore &lhs, FrontierScore &rhs) const
        {
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

        void collectComponent(ui best_u, FrontierState &state)
        {
            state.component_vertices.clear();
            state.component_frontier.clear();

            if (solver.mapped_q[best_u] != -1) {
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
                    if (solver.mapped_q[nbr] == -1 && solver.frontier_visit[nbr] != solver.frontier_token) {
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
#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
    struct LightweightSpokeLowerBound {
        MatchingSolver &solver;

        explicit LightweightSpokeLowerBound(MatchingSolver &solver) : solver(solver) {}

        bool dfsMatchLightSpoke(ui left_idx, const vector<vector<ui>> &adj)
        {
            for (ui right_idx : adj[left_idx]) {
                if (solver.lb_light_spoke_vis[right_idx]) continue;
                solver.lb_light_spoke_vis[right_idx] = 1;
                if (solver.lb_match_right[right_idx] < 0 ||
                    dfsMatchLightSpoke((ui)solver.lb_match_right[right_idx], adj)) {
                    solver.lb_match_right[right_idx] = (int)left_idx;
                    solver.lb_match_left[left_idx] = (int)right_idx;
                    return true;
                }
            }
            return false;
        }

        bool maxMatchLightSpokeWithBridge(const vector<vector<ui>> &adj, ui left_size, ui right_size, ui &missing_non_bridge)
        {
            std::fill(solver.lb_match_right.begin(),
                solver.lb_match_right.begin() + right_size, -1);
            std::fill(solver.lb_match_left.begin(),
                solver.lb_match_left.begin() + left_size, -1);

            ui non_bridge_count = 0;
            for (ui i = 0; i < left_size; ++i) {
                if (solver.lb_light_spoke_is_bridge[i]) {
                    std::fill(solver.lb_light_spoke_vis.begin(), solver.lb_light_spoke_vis.begin() + right_size, 0);
                    if (!dfsMatchLightSpoke(i, adj)) return false;
                }
                else {
                    non_bridge_count++;
                }
            }

            ui non_bridge_match = 0;
            for (ui i = 0; i < left_size; ++i) {
                if (solver.lb_light_spoke_is_bridge[i]) {
                    continue;
                }

                std::fill(solver.lb_light_spoke_vis.begin(),
                    solver.lb_light_spoke_vis.begin() + right_size, 0);
                if (dfsMatchLightSpoke(i, adj)) {
                    non_bridge_match++;
                }
            }

            missing_non_bridge = non_bridge_count - non_bridge_match;
            return true;
        }

        ui computeLightSpokeLB(ui u, ui v)
        {
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            ui deg_v = 0;
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            ui left_size = 0;
            for (ui i = 0; i < solver.q_degree[u]; ++i) {
                ui u1 = u_neighbors[i];
                if (solver.mapped_q[u1] != -1 || solver.excluded_edges[u][u1]) {
                    continue;
                }

                vector<ui> &adj = solver.lb_light_spoke_adj[left_size];
                adj.clear();
                solver.lb_light_spoke_is_bridge[left_size] =
                    solver.q_neighbor_is_bridge[u][i];

                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (v1 == v || solver.mapped_g[v1] != -1 ||
                        solver.excluded_cands[u1].contains(v1) ||
                        !solver.candidates[u1].contains(v1)) {
                        continue;
                    }
                    adj.push_back(j);
                }

                left_size++;
            }

            if (left_size == 0) {
                return 0;
            }

            ui missing_non_bridge = 0;
            if (!maxMatchLightSpokeWithBridge(solver.lb_light_spoke_adj,
                left_size, deg_v, missing_non_bridge)) {
                return (ui)INF;
            }
            return missing_non_bridge;
        }
    };
#endif
    // ========================================================================

    inline bool isInFrontier(ui u) const
    {
        return u < qn && frontier_pos[u] != -1;
    }

    ui computeEdgeSupport(ui u, ui anchor) const
    {
        if (u >= qn || anchor >= qn || mapped_q[anchor] == -1) {
            return 0;
        }

        ui count = 0;
        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
        for (ui i = 0; i < deg; ++i) {
            ui v = nbrs[i];
            if (!candidates[u].contains(v)) continue;
#if defined(CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT) || defined(ENABLE_MAPPED_VERTEX_SUPPORT)
            if (mapped_g[v] != -1) continue;
#endif
#if defined(CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT) || defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
            if (excluded_cands[u].contains(v)) continue;
#endif
            count++;
        }
        return count;
    }

    inline void clearFrontierEdgeSupport(ui u, ui anchor)
    {
#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
        if (u < qn && anchor < qn) {
            support[u][anchor] = 0;
        }
#else
        (void)u;
        (void)anchor;
#endif
    }

    inline void refreshFrontierEdgeSupport(ui u, ui anchor)
    {
#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
        if (u >= qn || anchor >= qn) {
            return;
        }

        if (isInFrontier(u) && mapped_q[anchor] != -1 && !excluded_edges[u][anchor]) {
            support[u][anchor] = computeEdgeSupport(u, anchor);
        }
        else {
            support[u][anchor] = 0;
        }
#else
        (void)u;
        (void)anchor;
#endif
    }

    void refreshFrontierVertexSupport(ui u)
    {
#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
        if (u >= qn) {
            return;
        }

        std::fill(support[u].begin(), support[u].end(), 0);
        if (!isInFrontier(u)) {
            return;
        }

        for (ui anchor : q_neighbors[u]) {
            refreshFrontierEdgeSupport(u, anchor);
        }
#else
        (void)u;
#endif
    }

    void clearFrontierVertexSupport(ui u)
    {
#ifdef CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT
        if (u < qn) {
            std::fill(support[u].begin(), support[u].end(), 0);
        }
#else
        (void)u;
#endif
    }

    // Update active_frontier and, in maintained mode, frontier-edge scores.
    inline void updateFrontierStatus(ui u)
    {
        bool should_be = (mapped_q[u] == -1 && anchor_count[u] > 0);
        bool is_in = (frontier_pos[u] != -1);

        if (should_be && !is_in) {
            frontier_pos[u] = active_frontier.size();
            active_frontier.push_back(u);
            refreshFrontierVertexSupport(u);
        }
        else if (!should_be && is_in) {
            clearFrontierVertexSupport(u);
            ui idx = frontier_pos[u];
            ui last_u = active_frontier.back();
            active_frontier[idx] = last_u;
            frontier_pos[last_u] = idx;
            active_frontier.pop_back();
            frontier_pos[u] = -1;
        }
    }

    // Update anchor counts and active_frontier when vertex u becomes matched/unmatched.
    void updateFrontier(ui u, bool matched)
    {
        if (matched) {
            updateFrontierStatus(u);
            clearFrontierVertexSupport(u);
        }

        for (ui u1 : q_neighbors[u]) {
            if (!excluded_edges[u1][u]) {
                if (matched) {
                    anchor_count[u1]++;
                    updateFrontierStatus(u1);
                    refreshFrontierEdgeSupport(u1, u);
                }
                else {
                    clearFrontierEdgeSupport(u1, u);
                    anchor_count[u1]--;
                    updateFrontierStatus(u1);
                }
            }
        }

        if (!matched) {
            updateFrontierStatus(u);
            refreshFrontierVertexSupport(u);
        }
    }

    void excludeFrontierEdge(ui u, ui anchor)
    {
        clearFrontierEdgeSupport(u, anchor);
        excluded_edges[u][anchor] = 1;
        excluded_edges[anchor][u] = 1;
        assert(anchor_count[u] > 0);
        anchor_count[u]--;
        updateFrontierStatus(u);
    }

    void restoreFrontierEdge(ui u, ui anchor)
    {
        bool was_frontier = isInFrontier(u);
        excluded_edges[u][anchor] = 0;
        excluded_edges[anchor][u] = 0;
        anchor_count[u]++;
        if (was_frontier) {
            refreshFrontierEdgeSupport(u, anchor);
        }
        updateFrontierStatus(u);
    }

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_MAPPED_VERTEX_SUPPORT)
    void adjustSupportForMappedDataVertex(ui v, bool restore)
    {
        if (v >= gn) {
            return;
        }

        for (ui u : active_frontier) {
            if (u >= qn || mapped_q[u] != -1) continue;
            if (!candidates[u].contains(v)) continue;
#ifdef ENABLE_EXCLUDED_EDGE_SUPPORT
            if (excluded_cands[u].contains(v)) continue;
#endif

            for (ui anchor : q_neighbors[u]) {
                if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
                if (!data_graph->hasEdge(v, (ui)mapped_q[anchor])) continue;

                if (restore) {
                    support[u][anchor]++;
                }
                else {
                    assert(support[u][anchor] > 0);
                    support[u][anchor]--;
                }
            }
        }
    }
#endif

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
    void decrementSupportForExcludedCandidate(ui u, ui excluded_anchor, ui v,
        vector<SupportDelta> &support_deltas)
    {
        if (u >= qn || mapped_q[u] != -1 || !isInFrontier(u)) {
            return;
        }

        for (ui anchor : q_neighbors[u]) {
            if (anchor == excluded_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
            if (!data_graph->hasEdge(v, (ui)mapped_q[anchor])) continue;

            assert(support[u][anchor] > 0);
            support[u][anchor]--;
            support_deltas.push_back({ u, anchor, 1 });
        }
    }

    void restoreSupportDeltas(const vector<SupportDelta> &support_deltas)
    {
        for (auto it = support_deltas.rbegin(); it != support_deltas.rend(); ++it) {
            assert(it->u < qn && it->anchor < qn);
            support[it->u][it->anchor] += it->amount;
        }
    }
#endif

    // ========================================================================
    // DFS Buffer
    // ========================================================================
    struct DfsBuffer {
        vector<ActiveEdge> top_edges;
        vector<ui> cand_v_list;
        BranchSelector branch_selector;
#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
        LightweightSpokeLowerBound light_lb;
#endif
        explicit DfsBuffer(MatchingSolver &solver)
            : branch_selector(solver)
#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
            , light_lb(solver)
#endif
        {}

        void reserve(ui threshold, ui max_g_deg, ui gn)
        {
            top_edges.reserve((size_t)threshold + 1);
            cand_v_list.reserve(max_g_deg);
            local_excluded_edges.reserve((size_t)threshold + 1);
            size_t cand_size = std::min((size_t)gn, ((size_t)threshold + 1) * (size_t)max_g_deg);
            local_excluded_cands.reserve(cand_size);
#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
            local_support_deltas.reserve(cand_size);
#endif
        }

        void clearLocalLogs()
        {
            local_excluded_edges.clear();
            local_excluded_cands.clear();
#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
            local_support_deltas.clear();
#endif
        }

        void recordExcludedEdge(ui u, ui anchor)
        {
            local_excluded_edges.push_back({ u, anchor });
        }

        void recordExcludedCandidate(ui u, ui v)
        {
            local_excluded_cands.push_back({ u, v });
        }

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
        vector<SupportDelta> &supportDeltas()
        {
            return local_support_deltas;
        }
#endif

        void restoreLocalChanges(MatchingSolver &solver)
        {
            for (auto &p : local_excluded_cands) {
                solver.excluded_cands[p.first].remove(p.second);
            }

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
            solver.restoreSupportDeltas(local_support_deltas);
#endif

            for (auto it = local_excluded_edges.rbegin(); it != local_excluded_edges.rend(); ++it) {
                const auto &e = *it;
                ui u = e.first;
                ui ua = e.second;
                assert(solver.mapped_q[u] == -1);
                assert(solver.mapped_q[ua] != -1);

                solver.restoreFrontierEdge(u, ua);
            }
        }

    private:
        vector<pair<ui, ui>> local_excluded_edges;
        vector<pair<ui, ui>> local_excluded_cands;
#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
        vector<SupportDelta> local_support_deltas;
#endif
    };

    vector<DfsBuffer> dfs_buffers;

    void reserveDfsBuffer(DfsBuffer &buf) const
    {
        buf.reserve(threshold, max_g_deg, gn);
    }

    void initDfsBuffer()
    {
        dfs_buffers.clear();
        dfs_buffers.reserve((size_t)qn + 1);
    }

    DfsBuffer &dfsBufferForDepth(size_t depth)
    {
        assert(depth <= qn);
        assert(dfs_buffers.capacity() >= (size_t)qn + 1);
        while (dfs_buffers.size() <= depth) {
            dfs_buffers.emplace_back(*this);
            reserveDfsBuffer(dfs_buffers.back());
        }
        return dfs_buffers[depth];
    }
    // ========================================================================

    // =====================================================
    // Procedure DFS(M_part, cost, X)
    //
    // cost:  current cost of partial match M_part
    // X:     the set of excluded edges (u, ua)
    // =====================================================
    void dfs(ui cost)
    {
        assert(part_M.size() <= qn);
        assert(cost <= threshold);

        if (part_M.size() == qn) {
            stats.recursion_calls++;
            stats.result_count++;
#ifndef NDEBUG
            recordBranchOrderAnswerDebug();
            results_ptr->push_back(part_M);
#endif
            return;
        }

        // TODO
        assert(!active_frontier.empty());
        if (active_frontier.empty()) return;

        stats.recursion_calls++;

        DfsBuffer &buf = dfsBufferForDepth(part_M.size());

        // choose the best edge (u, ua) and the corresponding frontier state
        Timer t_frontier;
        FrontierState current_state = buf.branch_selector.buildFrontierState(buf.top_edges);
        const vector<ui> &U_frontier = current_state.component_frontier;
        assert(!U_frontier.empty());
        stats.frontier_time += t_frontier.elapsed();

        Timer t_branch;
        long long child_dfs_time = 0;
        long long local_lb_time = 0;
        buf.clearLocalLogs();
        vector<ui> &cand_v_list = buf.cand_v_list;

        vector<ActiveEdge> &top_edges = buf.top_edges;
        ui current_cost = cost;

        while (current_cost <= threshold) {
            if (!buf.branch_selector.collectTopActiveEdges(U_frontier, 1, top_edges)) {
                break;
            }

            ActiveEdge edge = top_edges.front();
            ui u = edge.u;
            ui ua = edge.anchor;
            if (u >= qn || ua >= qn || mapped_q[u] != -1 ||
                mapped_q[ua] == -1 || excluded_edges[u][ua] ||
                frontier_pos[u] == -1) {
                break;
            }
            assert(q_matrix[u][ua]);

            buf.branch_selector.collectCandVertices(u, ua, cand_v_list);

            // Include branch: use this frontier-anchor edge to add exactly one new query vertex.
            for (ui v : cand_v_list) {
                ui delta = 0;
                for (ui anchor : q_neighbors[u]) {
                    if (anchor == ua) continue;
                    if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
                    if (!data_graph->hasEdge(v, (ui)mapped_q[anchor])) delta++;
                }
                if (current_cost + delta > threshold) continue;

#ifdef CDE_LB_LIGHTWEIGHT_SPOKE
                {
                    Timer t_lb;
                    ui light_spoke_lb = buf.light_lb.computeLightSpokeLB(u, v);
                    long long elapsed = t_lb.elapsed();
                    stats.lb_time += elapsed;
                    stats.lb_light_spoke_time += elapsed;
                    local_lb_time += elapsed;
                    if (light_spoke_lb > threshold - (current_cost + delta)) {
                        stats.prun_calls++;
                        continue;
                    }
                }
#endif

#ifndef NDEBUG
                recordBranchOrderDebug(u);
#endif

#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_MAPPED_VERTEX_SUPPORT)
                adjustSupportForMappedDataVertex(v, false);
#endif
                mapped_q[u] = (int)v;
                mapped_g[v] = (int)u;
                part_M.push_back({ u, v });

                updateFrontier(u, true);

                Timer t_child;
                dfs(current_cost + delta);
                child_dfs_time += t_child.elapsed();

                part_M.pop_back();
#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_MAPPED_VERTEX_SUPPORT)
                mapped_q[u] = -1;
                updateFrontier(u, false);
                mapped_g[v] = -1;
                adjustSupportForMappedDataVertex(v, true);
#else
                mapped_g[v] = -1;
                mapped_q[u] = -1;

                updateFrontier(u, false);
#endif
            }

            // Exclude branch: keep this decision in the current frame, then
            // continue to the next selected active edge without a recursive
            // call. Every recursive child above still adds a vertex.
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);

            if (current_cost > threshold) {
                break;
            }

#ifndef NDEBUG
            recordBranchOrderDebug(u);
#endif
            for (ui v : cand_v_list) {
                if (!excluded_cands[u].contains(v)) {
                    excluded_cands[u].insert(v);
                    buf.recordExcludedCandidate(u, v);
#if defined(CDE_EDGE_IE_MAINTAIN_ANCHOR_SUPPORT) && defined(ENABLE_EXCLUDED_EDGE_SUPPORT)
                    decrementSupportForExcludedCandidate(u, ua, v, buf.supportDeltas());
#endif
                }
            }
        }

        {
            long long branch_elapsed = t_branch.elapsed();
            long long excluded_time = child_dfs_time + local_lb_time;
            stats.branch_time += branch_elapsed > excluded_time
                ? branch_elapsed - excluded_time : 0;
        }

        buf.restoreLocalChanges(*this);
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
