#ifndef MATCHING_ALGORITHMS_CDE_EDGE_IE_H_
#define MATCHING_ALGORITHMS_CDE_EDGE_IE_H_

#include "graph/graph.h"
#include "matching/run_matching.h"
#include "utility/utility.h"
#include "utility/mybitset.h"
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

        initQueryGraphIndex();

        Timer t_filter;
        bool res = runCandidateFiltering();
        stats.filter_time = t_filter.elapsed();
        if (!res) {
            stats.init_time = t_init.elapsed();
            return false;
        }

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        initDataLabelDegreeIndex();
#endif
#if CDE_EDGE_IE_FIXED_ORDER
        initFixedEdgePriorities();
#endif
        stats.init_time = t_init.elapsed();
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
        const char *forced_root_env = std::getenv("CDE_EDGE_IE_ROOT");
        if (forced_root_env != nullptr && forced_root_env[0] != '\0') {
            char *end = nullptr;
            long forced_root = std::strtol(forced_root_env, &end, 10);
            if (end != forced_root_env && *end == '\0' &&
                forced_root >= 0 && (ui)forced_root < qn) {
                root = (ui)forced_root;
            }
            else {
                fprintf(stderr, "Ignoring invalid CDE_EDGE_IE_ROOT=%s\n", forced_root_env);
            }
        }
        ui forced_root_data = 0;
        bool has_forced_root_data = parseEnvUi("CDE_EDGE_IE_ROOT_DATA", forced_root_data);
        ui forced_root_cand_index = 0;
        bool has_forced_root_cand_index =
            parseEnvUi("CDE_EDGE_IE_ROOT_CAND_INDEX", forced_root_cand_index);
        // if (kDebugInitialRoot >= 0) {
        //     assert((ui)kDebugInitialRoot < qn);
        //     root = (ui)kDebugInitialRoot;
        // }
        printf("Selected initial root: u=%u with %zu candidates\n", root, (size_t)candidates[root].size());
        printBranchProfileConfig(root, has_forced_root_data, forced_root_data,
            has_forced_root_cand_index, forced_root_cand_index);
        ui root_cand_index = 0;
#endif

        for (ui v0 : candidates[root]) {
#ifndef NDEBUG
            ui this_root_cand_index = root_cand_index++;
            if (has_forced_root_data && v0 != forced_root_data) continue;
            if (has_forced_root_cand_index && this_root_cand_index != forced_root_cand_index) continue;
            recordBranchOrderDebug(root);
#endif

#ifndef NDEBUG
            recordBranchProfileRootStart(root, v0);
#endif
            mapped_q[root] = (int)v0;
            mapped_g[v0] = (int)root;
            part_M.push_back({ root, v0 });

            updateFrontier(root, true);

            dfs(0);

            part_M.pop_back();

            mapped_g[v0] = -1;
            mapped_q[root] = -1;
            updateFrontier(root, false);
#ifndef NDEBUG
            recordBranchProfileRootFinish(root, v0);
#endif

            if (outputLimitReached()) break;
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
        ui filter_candidate_count = 0;
        // search breakdown
        long long dfs_time = 0;
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
        long long branch_cal_edge_support_time = 0;
        long long branch_count_anchors_time = 0;
        long long support_update_time = 0;
        long long candidate_loop_time = 0;
        long long exclude_update_time = 0;
        // counters
        long long recursion_calls = 0;
        long long prun_calls = 0;
        size_t result_count = 0;
        bool output_limit_reached = false;
#ifndef NDEBUG
        unsigned long long terminal_tail_calls = 0;
        unsigned long long terminal_prune_calls = 0;
        unsigned long long terminal_delayed_vertices = 0;
        unsigned long long terminal_bucket_candidate_checks = 0;
#endif
    } stats;

    void printStats() const
    {
        auto pct = [](long long part, long long whole) -> double {
            return whole > 0 ? (double)part / whole * 100.0 : 0.0;
            };

        long long search_accounted_time = stats.frontier_time + stats.branch_time;
        long long search_other_time = stats.dfs_time > search_accounted_time
            ? stats.dfs_time - search_accounted_time : 0;

        printf("\n--- CDE-Edge-IE Time Analysis ---\n");
#ifdef NDEBUG
        printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
        printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
        printf("  - Filter Time:     %.4lf ms\n", stats.filter_time / 1000.0);
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
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
        printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
        printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.dfs_time / 1000.0, pct(stats.dfs_time, stats.total_time));
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
        printf("    - branch_cal_edge_support_ms: %.4lf (%.2f%% of Branch)\n", stats.branch_cal_edge_support_time / 1000.0, pct(stats.branch_cal_edge_support_time, stats.branch_time));
        printf("    - branch_count_anchors_ms:    %.4lf (%.2f%% of Branch)\n", stats.branch_count_anchors_time / 1000.0, pct(stats.branch_count_anchors_time, stats.branch_time));
        printf("    - support_update_ms:          %.4lf (%.2f%% of Branch)\n", stats.support_update_time / 1000.0, pct(stats.support_update_time, stats.branch_time));
        printf("    - candidate_loop_ms:          %.4lf (%.2f%% of Branch)\n", stats.candidate_loop_time / 1000.0, pct(stats.candidate_loop_time, stats.branch_time));
        printf("    - exclude_update_ms:          %.4lf (%.2f%% of Branch)\n", stats.exclude_update_time / 1000.0, pct(stats.exclude_update_time, stats.branch_time));
        printf("  - Search Other:    %.4lf ms (%.2f%% of Search)\n", search_other_time / 1000.0, pct(search_other_time, stats.dfs_time));
#endif
        printf("Recursion Calls:     %lld\n", stats.recursion_calls);
        printf("Pruning Calls:       %lld\n", stats.prun_calls);
#ifndef NDEBUG
        printf("Terminal Tail Calls: %llu\n", stats.terminal_tail_calls);
        printf("Terminal Prunes:     %llu\n", stats.terminal_prune_calls);
        printf("Terminal Delayed:    %llu\n", stats.terminal_delayed_vertices);
        printf("Terminal Cand Checks:%llu\n", stats.terminal_bucket_candidate_checks);
#endif
        printf("Results Found:       %zu\n", stats.result_count);
#if MATCH_OUTPUT_LIMIT > 0
        printf("Output Limit:        %zu%s\n",
            (size_t)MATCH_OUTPUT_LIMIT,
            stats.output_limit_reached ? " (reached)" : "");
#endif
#ifndef NDEBUG
        printBranchOrderDebugStats();
        printBranchProfileStats("final", branch_profile_timer.elapsed());
#endif
        printf("-----------------------------------------------------------\n");
        fflush(stdout);
    }

private:
#ifndef NDEBUG
    // Internal debug knob: set to 0 to disable branch-order prefix counting.
    enum { kDebugBranchOrderDepth = 1 };
#endif

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
    struct DataLabelDegreeCount {
        LabelID label = 0;
        ui count = 0;
    };
#endif

    struct ActiveEdge {
        ui u = 0;       // unmatched endpoint
        ui anchor = 0;  // matched endpoint
        ui anchor_support = std::numeric_limits<ui>::max();
        double rank_support = std::numeric_limits<double>::max();
        ui live_anchor_count = 0;
        ui query_degree = 0;
    };

    struct EdgeScoreCache {
        vector<vector<ActiveEdge>> active_edges_by_vertex;
        vector<char> active_edges_cached;
    };

    struct SupportSnapshot {
        ui u = 0;
        ui anchor = 0;
        ui value = 0;
        char dirty = 0;
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
#if CDE_EDGE_IE_FIXED_ORDER
    vector<vector<ui>> fixed_edge_priority;
#endif
#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
    vector<vector<DataLabelDegreeCount>> data_label_degrees;
#endif

    vector<MyBitset> candidates;

    vector<int> mapped_q;
    vector<int> mapped_g;
    vector<vector<char>> excluded_edges;
    vector<MyBitset> excluded_cands;
    vector<pair<ui, ui>> part_M;

    vector<ui> anchor_count;
    vector<vector<ui>> support;
    vector<vector<char>> support_dirty;
    vector<SupportSnapshot> *active_support_snapshots = nullptr;
    vector<int> frontier_pos;
    vector<ui> active_frontier;
    vector<ui> data_vertex_mark;
    vector<ui> data_vertex_mark_pos;
#ifndef NDEBUG
    vector<char> delayed_query_vertices;
    vector<char> preferred_query_vertices;
#endif
    bool terminal_buckets_enabled = true;
    ui data_vertex_mark_token = 0;
    int branch_timing_depth = 0;

    struct TerminalScan {
        ui unmatched_count = 0;
        ui terminal_count = 0;
        ui terminal_frontier_count = 0;
        ui nonterminal_frontier_count = 0;

        bool allRemainingTerminal() const
        {
            return unmatched_count > 0 && unmatched_count == terminal_count;
        }

        bool hasNonterminalFrontier() const
        {
            return nonterminal_frontier_count > 0;
        }
    };

    struct TerminalTailVertex {
        ui u = 0;
        ui feasible_count = 0;
        ui min_delta = std::numeric_limits<ui>::max();
        vector<vector<ui>> buckets;
    };

#ifndef NDEBUG
    struct BranchProfileCounter {
        unsigned long long root_enters = 0;
        unsigned long long root_finishes = 0;
        unsigned long long dfs_calls = 0;
        unsigned long long dfs_returns = 0;
        unsigned long long selected_edges = 0;
        unsigned long long support_candidates = 0;
        unsigned long long include_considered = 0;
        unsigned long long include_taken = 0;
        unsigned long long include_cost_pruned = 0;
        unsigned long long exclude_taken = 0;
        unsigned long long exclude_over_budget = 0;
        unsigned long long forced_zero_excludes = 0;
        unsigned long long prunes = 0;
        unsigned long long answers = 0;
    };

    struct SecondBranchKey {
        ui root_v = 0;
        ui u = 0;
        ui v = 0;

        bool operator<(const SecondBranchKey &other) const
        {
            if (root_v != other.root_v) return root_v < other.root_v;
            if (u != other.u) return u < other.u;
            return v < other.v;
        }
    };

    struct RootEdgeBranchKey {
        ui root_v = 0;
        ui u = 0;
        ui anchor = 0;

        bool operator<(const RootEdgeBranchKey &other) const
        {
            if (root_v != other.root_v) return root_v < other.root_v;
            if (u != other.u) return u < other.u;
            return anchor < other.anchor;
        }
    };

    bool branch_profile_enabled = false;
    ui branch_profile_root = std::numeric_limits<ui>::max();
    ui branch_profile_top_n = 20;
    ui branch_profile_path_depth = 0;
    long long branch_profile_interval_us = 2000000;
    unsigned long long branch_profile_event_ticks = 0;
    mutable Timer branch_profile_timer;
    mutable long long branch_profile_last_dump_us = 0;
    ui active_branch_profile_root_v = std::numeric_limits<ui>::max();
    map<ui, BranchProfileCounter> root_branch_profiles;
    map<SecondBranchKey, BranchProfileCounter> second_branch_profiles;
    map<RootEdgeBranchKey, BranchProfileCounter> root_edge_branch_profiles;
    map<vector<pair<ui, ui>>, BranchProfileCounter> path_branch_profiles;
    map<pair<ui, ui>, BranchProfileCounter> extension_depth_profiles;
    vector<char> profile_query_vertices;
    map<vector<ui>, unsigned long long> branch_order_prefix_counts;
    map<vector<ui>, unsigned long long> branch_order_prefix_reach_counts;
    map<vector<ui>, unsigned long long> branch_order_prefix_answer_counts;
#endif

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
#if CDE_EDGE_IE_FIXED_ORDER
        fixed_edge_priority.clear();
#endif
#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        data_label_degrees.clear();
#endif

        anchor_count.assign(qn, 0);
        support.assign(qn, vector<ui>(qn, 0));
        support_dirty.assign(qn, vector<char>(qn, 1));
        active_support_snapshots = nullptr;
        frontier_pos.assign(qn, -1);
        active_frontier.clear();
        data_vertex_mark.assign(gn, 0);
        data_vertex_mark_pos.assign(gn, 0);
#ifndef NDEBUG
        delayed_query_vertices.assign(qn, 0);
        preferred_query_vertices.assign(qn, 0);
        profile_query_vertices.assign(qn, 0);
#endif
        data_vertex_mark_token = 0;
#ifndef NDEBUG
        branch_order_prefix_counts.clear();
        branch_order_prefix_reach_counts.clear();
        branch_order_prefix_answer_counts.clear();
#endif
        stats = TimeStats();
        initTerminalTailConfig();
#ifndef NDEBUG
        initBranchProfileConfig();
        initDelayedQueryVertices();
#endif
    }

#ifndef NDEBUG
    static bool parseEnvUi(const char *name, ui &value)
    {
        const char *env = std::getenv(name);
        if (env == nullptr || env[0] == '\0') {
            return false;
        }
        char *end = nullptr;
        unsigned long parsed = std::strtoul(env, &end, 10);
        if (end == env || *end != '\0') {
            fprintf(stderr, "Ignoring invalid %s=%s\n", name, env);
            return false;
        }
        value = (ui)parsed;
        return true;
    }

    bool parseQueryVertexListEnv(const char *name, vector<char> &marks)
    {
        const char *env = std::getenv(name);
        if (env == nullptr || env[0] == '\0') {
            return false;
        }

        bool parsed_any = false;
        const char *p = env;
        while (*p != '\0') {
            while (*p == ',' || *p == ' ' || *p == '\t') {
                p++;
            }
            if (*p == '\0') {
                break;
            }

            char *end = nullptr;
            unsigned long parsed = std::strtoul(p, &end, 10);
            if (end == p) {
                fprintf(stderr, "Ignoring invalid %s near: %s\n", name, p);
                break;
            }
            if (parsed < qn) {
                marks[(ui)parsed] = 1;
                parsed_any = true;
            }
            else {
                fprintf(stderr, "Ignoring out-of-range query vertex in %s: %lu\n", name, parsed);
            }
            p = end;
        }
        return parsed_any;
    }

    void initDelayedQueryVertices()
    {
        parseQueryVertexListEnv("CDE_EDGE_IE_DELAY_VERTICES", delayed_query_vertices);
        parseQueryVertexListEnv("CDE_EDGE_IE_PREFER_VERTICES", preferred_query_vertices);
    }
#endif

    void initTerminalTailConfig()
    {
        terminal_buckets_enabled = CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT != 0;
    }

#ifndef NDEBUG
    static unsigned long long branchProfileBranchTotal(const BranchProfileCounter &counter)
    {
        return counter.include_taken + counter.exclude_taken + counter.forced_zero_excludes;
    }

    void initBranchProfileConfig()
    {
        root_branch_profiles.clear();
        second_branch_profiles.clear();
        root_edge_branch_profiles.clear();
        path_branch_profiles.clear();
        extension_depth_profiles.clear();
        branch_profile_enabled = false;
        branch_profile_root = std::numeric_limits<ui>::max();
        branch_profile_top_n = 20;
        branch_profile_path_depth = 0;
        branch_profile_interval_us = 2000000;
        branch_profile_event_ticks = 0;
        active_branch_profile_root_v = std::numeric_limits<ui>::max();

        const char *profile_env = std::getenv("CDE_EDGE_IE_PROFILE");
        if (profile_env != nullptr && profile_env[0] != '\0' &&
            std::strcmp(profile_env, "0") != 0) {
            branch_profile_enabled = true;
        }

        ui parsed = 0;
        if (parseEnvUi("CDE_EDGE_IE_PROFILE_ROOT", parsed)) {
            branch_profile_enabled = true;
            branch_profile_root = parsed;
        }
        else if (branch_profile_enabled && parseEnvUi("CDE_EDGE_IE_ROOT", parsed)) {
            branch_profile_root = parsed;
        }

        if (parseEnvUi("CDE_EDGE_IE_PROFILE_TOP_N", parsed) && parsed > 0) {
            branch_profile_top_n = parsed;
        }
        if (parseEnvUi("CDE_EDGE_IE_PROFILE_DEPTH", parsed)) {
            branch_profile_enabled = true;
            branch_profile_path_depth = std::min(parsed, qn);
        }
        else if (branch_profile_enabled) {
            branch_profile_path_depth = std::min((ui)2, qn);
        }
        if (parseQueryVertexListEnv("CDE_EDGE_IE_PROFILE_VERTICES", profile_query_vertices)) {
            branch_profile_enabled = true;
        }
        if (parseEnvUi("CDE_EDGE_IE_PROFILE_INTERVAL_MS", parsed)) {
            branch_profile_interval_us = (long long)parsed * 1000;
        }

        branch_profile_timer.restart();
        branch_profile_last_dump_us = 0;
    }

    bool shouldProfileRoot(ui root) const
    {
        return branch_profile_enabled &&
            (branch_profile_root == std::numeric_limits<ui>::max() ||
                branch_profile_root == root);
    }

    bool shouldProfileCurrentBranch() const
    {
        return branch_profile_enabled && !part_M.empty() &&
            shouldProfileRoot(part_M[0].first);
    }

    BranchProfileCounter *currentRootProfile()
    {
        if (!shouldProfileCurrentBranch()) {
            return nullptr;
        }
        return &root_branch_profiles[part_M[0].second];
    }

    BranchProfileCounter *currentSecondProfile()
    {
        if (!shouldProfileCurrentBranch() || part_M.size() < 2) {
            return nullptr;
        }
        SecondBranchKey key;
        key.root_v = part_M[0].second;
        key.u = part_M[1].first;
        key.v = part_M[1].second;
        return &second_branch_profiles[key];
    }

    BranchProfileCounter *currentRootEdgeProfile(ui u, ui anchor)
    {
        if (!shouldProfileCurrentBranch() || part_M.size() != 1) {
            return nullptr;
        }
        RootEdgeBranchKey key;
        key.root_v = part_M[0].second;
        key.u = u;
        key.anchor = anchor;
        return &root_edge_branch_profiles[key];
    }

    template <typename Updater>
    void updateBranchProfilePaths(Updater updater)
    {
        if (!shouldProfileCurrentBranch() || branch_profile_path_depth == 0) {
            return;
        }

        ui limit = std::min(branch_profile_path_depth, (ui)part_M.size());
        vector<pair<ui, ui>> prefix;
        prefix.reserve(limit);
        for (ui i = 0; i < limit; ++i) {
            prefix.push_back(part_M[i]);
            updater(path_branch_profiles[prefix]);
            updater(extension_depth_profiles[{ (ui)prefix.size(), prefix.back().first }]);
        }
    }

    void maybePrintBranchProfile()
    {
        if (!branch_profile_enabled || branch_profile_interval_us <= 0) {
            return;
        }
        branch_profile_event_ticks++;
        if ((branch_profile_event_ticks & 0x3fffULL) != 0) {
            return;
        }

        long long elapsed_us = branch_profile_timer.elapsed();
        if (elapsed_us - branch_profile_last_dump_us >= branch_profile_interval_us) {
            branch_profile_last_dump_us = elapsed_us;
            printBranchProfileStats("partial", elapsed_us);
        }
    }

    void printBranchProfileConfig(ui root, bool has_forced_root_data, ui forced_root_data,
        bool has_forced_root_cand_index, ui forced_root_cand_index)
    {
        if (!shouldProfileRoot(root)) {
            return;
        }

        printf("\n--- CDE Root Branch Profile Enabled ---\n");
        printf("Profile Query Root:  u=%u\n", root);
        printf("Root Candidates:     %zu\n", (size_t)candidates[root].size());
        printf("Top Rows:            %u\n", branch_profile_top_n);
        printf("Path Depth:          %u\n", branch_profile_path_depth);
        printf("Interval:            %.3lf ms\n", branch_profile_interval_us / 1000.0);
        bool has_profile_focus = false;
        for (ui u = 0; u < profile_query_vertices.size(); ++u) {
            if (profile_query_vertices[u]) {
                has_profile_focus = true;
                break;
            }
        }
        if (has_profile_focus) {
            printf("Profile Vertices:");
            for (ui u = 0; u < profile_query_vertices.size(); ++u) {
                if (profile_query_vertices[u]) {
                    printf(" u%u", u);
                }
            }
            printf("\n");
        }
        bool has_preferred_vertices = false;
        for (ui u = 0; u < preferred_query_vertices.size(); ++u) {
            if (preferred_query_vertices[u]) {
                has_preferred_vertices = true;
                break;
            }
        }
        if (has_preferred_vertices) {
            printf("Preferred Vertices:");
            for (ui u = 0; u < preferred_query_vertices.size(); ++u) {
                if (preferred_query_vertices[u]) {
                    printf(" u%u", u);
                }
            }
            printf("\n");
        }
        if (has_forced_root_data) {
            printf("Root Data Filter:    v=%u\n", forced_root_data);
        }
        if (has_forced_root_cand_index) {
            printf("Root Index Filter:   #%u\n", forced_root_cand_index);
        }
        printf("Root Candidate Order:");
        ui idx = 0;
        for (ui v : candidates[root]) {
            if (has_forced_root_data && v != forced_root_data) {
                idx++;
                continue;
            }
            if (has_forced_root_cand_index && idx != forced_root_cand_index) {
                idx++;
                continue;
            }
            printf(" #%u=%u", idx, v);
            idx++;
        }
        printf("\n-----------------------------------------------------------\n");
        fflush(stdout);
    }

    void recordBranchProfileRootStart(ui root, ui v0)
    {
        if (!shouldProfileRoot(root)) {
            return;
        }
        active_branch_profile_root_v = v0;
        root_branch_profiles[v0].root_enters++;
        maybePrintBranchProfile();
    }

    void recordBranchProfileRootFinish(ui root, ui v0)
    {
        if (!shouldProfileRoot(root)) {
            return;
        }
        root_branch_profiles[v0].root_finishes++;
        active_branch_profile_root_v = std::numeric_limits<ui>::max();
        maybePrintBranchProfile();
    }

    void recordBranchProfileDfsEnter()
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->dfs_calls++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->dfs_calls++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.dfs_calls++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileDfsReturn()
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->dfs_returns++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->dfs_returns++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.dfs_returns++;
        });
        maybePrintBranchProfile();
    }

    struct BranchProfileDfsScope {
        MatchingSolver &solver;

        explicit BranchProfileDfsScope(MatchingSolver &solver)
            : solver(solver)
        {}

        ~BranchProfileDfsScope()
        {
            solver.recordBranchProfileDfsReturn();
        }
    };

    void recordBranchProfilePrune()
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->prunes++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->prunes++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.prunes++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileAnswer()
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->answers++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->answers++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.answers++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileSelectedEdge(ui u, ui anchor, size_t support_count)
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->selected_edges++;
        root_counter->support_candidates += support_count;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->selected_edges++;
            second_counter->support_candidates += support_count;
        }
        BranchProfileCounter *edge_counter = currentRootEdgeProfile(u, anchor);
        if (edge_counter != nullptr) {
            edge_counter->selected_edges++;
            edge_counter->support_candidates += support_count;
        }
        updateBranchProfilePaths([support_count](BranchProfileCounter &counter) {
            counter.selected_edges++;
            counter.support_candidates += support_count;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileForcedZeroExclude(ui u, ui anchor)
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->forced_zero_excludes++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->forced_zero_excludes++;
        }
        BranchProfileCounter *edge_counter = currentRootEdgeProfile(u, anchor);
        if (edge_counter != nullptr) {
            edge_counter->forced_zero_excludes++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.forced_zero_excludes++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileIncludeConsidered(ui u, ui anchor)
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->include_considered++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->include_considered++;
        }
        BranchProfileCounter *edge_counter = currentRootEdgeProfile(u, anchor);
        if (edge_counter != nullptr) {
            edge_counter->include_considered++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.include_considered++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileIncludeTaken(ui u, ui anchor)
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->include_taken++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->include_taken++;
        }
        BranchProfileCounter *edge_counter = currentRootEdgeProfile(u, anchor);
        if (edge_counter != nullptr) {
            edge_counter->include_taken++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.include_taken++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileIncludeCostPruned(ui u, ui anchor)
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->include_cost_pruned++;
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->include_cost_pruned++;
        }
        BranchProfileCounter *edge_counter = currentRootEdgeProfile(u, anchor);
        if (edge_counter != nullptr) {
            edge_counter->include_cost_pruned++;
        }
        updateBranchProfilePaths([](BranchProfileCounter &counter) {
            counter.include_cost_pruned++;
        });
        maybePrintBranchProfile();
    }

    void recordBranchProfileExclude(ui u, ui anchor, bool over_budget)
    {
        BranchProfileCounter *root_counter = currentRootProfile();
        if (root_counter == nullptr) {
            return;
        }
        root_counter->exclude_taken++;
        if (over_budget) {
            root_counter->exclude_over_budget++;
        }
        BranchProfileCounter *second_counter = currentSecondProfile();
        if (second_counter != nullptr) {
            second_counter->exclude_taken++;
            if (over_budget) {
                second_counter->exclude_over_budget++;
            }
        }
        BranchProfileCounter *edge_counter = currentRootEdgeProfile(u, anchor);
        if (edge_counter != nullptr) {
            edge_counter->exclude_taken++;
            if (over_budget) {
                edge_counter->exclude_over_budget++;
            }
        }
        updateBranchProfilePaths([over_budget](BranchProfileCounter &counter) {
            counter.exclude_taken++;
            if (over_budget) {
                counter.exclude_over_budget++;
            }
        });
        maybePrintBranchProfile();
    }

    void printBranchProfileCounter(const BranchProfileCounter &counter) const
    {
        unsigned long long open_dfs =
            counter.dfs_calls >= counter.dfs_returns ?
                counter.dfs_calls - counter.dfs_returns : 0;
        printf("branches=%llu, dfs=%llu, done=%llu, open=%llu, edge_select=%llu, support_sum=%llu, "
            "inc_considered=%llu, inc_taken=%llu, inc_cost_pruned=%llu, "
            "exclude=%llu, exclude_over=%llu, "
            "forced0=%llu, prunes=%llu, answers=%llu",
            branchProfileBranchTotal(counter), counter.dfs_calls,
            counter.dfs_returns, open_dfs,
            counter.selected_edges, counter.support_candidates,
            counter.include_considered, counter.include_taken,
            counter.include_cost_pruned, counter.exclude_taken, counter.exclude_over_budget,
            counter.forced_zero_excludes, counter.prunes, counter.answers);
    }

    void printBranchProfilePath(const vector<pair<ui, ui>> &path) const
    {
        printf("[");
        for (ui i = 0; i < path.size(); ++i) {
            if (i > 0) printf(",");
            printf("u%u->v%u", path[i].first, path[i].second);
        }
        printf("]");
    }

    void printBranchProfileStats(const char *phase, long long elapsed_us) const
    {
        if (!branch_profile_enabled) {
            return;
        }

        printf("\n--- CDE Root Branch Profile (%s, %.3lf ms) ---\n",
            phase, elapsed_us / 1000.0);
        printf("Results Found So Far:%zu\n", stats.result_count);
        printf("Recursion So Far:    %lld\n", stats.recursion_calls);
        if (active_branch_profile_root_v != std::numeric_limits<ui>::max()) {
            printf("Active Root Data:    v=%u\n", active_branch_profile_root_v);
        }

        vector<pair<ui, BranchProfileCounter>> root_rows;
        root_rows.reserve(root_branch_profiles.size());
        for (map<ui, BranchProfileCounter>::const_iterator it = root_branch_profiles.begin();
            it != root_branch_profiles.end(); ++it) {
            root_rows.push_back(*it);
        }
        std::sort(root_rows.begin(), root_rows.end(),
            [](const pair<ui, BranchProfileCounter> &a,
                const pair<ui, BranchProfileCounter> &b) {
                unsigned long long ab = branchProfileBranchTotal(a.second);
                unsigned long long bb = branchProfileBranchTotal(b.second);
                if (ab != bb) return ab > bb;
                if (a.second.dfs_calls != b.second.dfs_calls) {
                    return a.second.dfs_calls > b.second.dfs_calls;
                }
                return a.first < b.first;
            });

        printf("Root candidate rows:\n");
        ui printed = 0;
        for (vector<pair<ui, BranchProfileCounter>>::const_iterator it = root_rows.begin();
            it != root_rows.end() && printed < branch_profile_top_n; ++it, ++printed) {
            const BranchProfileCounter &counter = it->second;
            printf("  v10=%u, state=%s, enters=%llu, finishes=%llu, ",
                it->first,
                (it->first == active_branch_profile_root_v ? "active" :
                    (counter.root_finishes > 0 ? "done" : "pending")),
                counter.root_enters, counter.root_finishes);
            printBranchProfileCounter(counter);
            printf("\n");
        }
        if (root_rows.empty()) {
            printf("  <none>\n");
        }

        struct EdgeRow {
            RootEdgeBranchKey key;
            BranchProfileCounter counter;
        };
        vector<EdgeRow> edge_rows;
        edge_rows.reserve(root_edge_branch_profiles.size());
        for (map<RootEdgeBranchKey, BranchProfileCounter>::const_iterator it =
            root_edge_branch_profiles.begin(); it != root_edge_branch_profiles.end(); ++it) {
            EdgeRow row;
            row.key = it->first;
            row.counter = it->second;
            edge_rows.push_back(row);
        }
        std::sort(edge_rows.begin(), edge_rows.end(),
            [](const EdgeRow &a, const EdgeRow &b) {
                unsigned long long ab = branchProfileBranchTotal(a.counter);
                unsigned long long bb = branchProfileBranchTotal(b.counter);
                if (ab != bb) return ab > bb;
                if (a.counter.support_candidates != b.counter.support_candidates) {
                    return a.counter.support_candidates > b.counter.support_candidates;
                }
                if (a.key.root_v != b.key.root_v) return a.key.root_v < b.key.root_v;
                if (a.key.u != b.key.u) return a.key.u < b.key.u;
                return a.key.anchor < b.key.anchor;
            });

        printf("Root-depth edge rows:\n");
        printed = 0;
        for (vector<EdgeRow>::const_iterator it = edge_rows.begin();
            it != edge_rows.end() && printed < branch_profile_top_n; ++it, ++printed) {
            printf("  v10=%u, edge=(u%u-u%u), ",
                it->key.root_v, it->key.u, it->key.anchor);
            printBranchProfileCounter(it->counter);
            printf("\n");
        }
        if (edge_rows.empty()) {
            printf("  <none>\n");
        }

        struct SecondRow {
            SecondBranchKey key;
            BranchProfileCounter counter;
        };
        vector<SecondRow> second_rows;
        second_rows.reserve(second_branch_profiles.size());
        for (map<SecondBranchKey, BranchProfileCounter>::const_iterator it =
            second_branch_profiles.begin(); it != second_branch_profiles.end(); ++it) {
            SecondRow row;
            row.key = it->first;
            row.counter = it->second;
            second_rows.push_back(row);
        }
        std::sort(second_rows.begin(), second_rows.end(),
            [](const SecondRow &a, const SecondRow &b) {
                unsigned long long ab = branchProfileBranchTotal(a.counter);
                unsigned long long bb = branchProfileBranchTotal(b.counter);
                if (ab != bb) return ab > bb;
                if (a.counter.dfs_calls != b.counter.dfs_calls) {
                    return a.counter.dfs_calls > b.counter.dfs_calls;
                }
                if (a.key.root_v != b.key.root_v) return a.key.root_v < b.key.root_v;
                if (a.key.u != b.key.u) return a.key.u < b.key.u;
                return a.key.v < b.key.v;
            });

        printf("Second-extension rows:\n");
        printed = 0;
        for (vector<SecondRow>::const_iterator it = second_rows.begin();
            it != second_rows.end() && printed < branch_profile_top_n; ++it, ++printed) {
            printf("  v10=%u, second=(u%u->v%u), ",
                it->key.root_v, it->key.u, it->key.v);
            printBranchProfileCounter(it->counter);
            printf("\n");
        }
        if (second_rows.empty()) {
            printf("  <none>\n");
        }

        struct PathRow {
            vector<pair<ui, ui>> path;
            BranchProfileCounter counter;
        };
        vector<PathRow> path_rows;
        path_rows.reserve(path_branch_profiles.size());
        for (map<vector<pair<ui, ui>>, BranchProfileCounter>::const_iterator it =
            path_branch_profiles.begin(); it != path_branch_profiles.end(); ++it) {
            PathRow row;
            row.path = it->first;
            row.counter = it->second;
            path_rows.push_back(row);
        }
        std::sort(path_rows.begin(), path_rows.end(),
            [](const PathRow &a, const PathRow &b) {
                unsigned long long ab = branchProfileBranchTotal(a.counter);
                unsigned long long bb = branchProfileBranchTotal(b.counter);
                if (ab != bb) return ab > bb;
                if (a.counter.dfs_calls != b.counter.dfs_calls) {
                    return a.counter.dfs_calls > b.counter.dfs_calls;
                }
                return a.path < b.path;
            });

        printf("Path-prefix rows:\n");
        if (branch_profile_path_depth == 0) {
            printf("  <disabled>\n");
        }
        for (ui depth = 1; depth <= branch_profile_path_depth; ++depth) {
            ui depth_printed = 0;
            bool depth_header_printed = false;
            for (vector<PathRow>::const_iterator it = path_rows.begin();
                it != path_rows.end() && depth_printed < branch_profile_top_n; ++it) {
                if (it->path.size() != depth) {
                    continue;
                }
                if (!depth_header_printed) {
                    printf("  Depth %u:\n", depth);
                    depth_header_printed = true;
                }
                printf("    ");
                printBranchProfilePath(it->path);
                printf(", ");
                printBranchProfileCounter(it->counter);
                printf("\n");
                depth_printed++;
            }
            if (!depth_header_printed && !path_rows.empty()) {
                printf("  Depth %u: <none>\n", depth);
            }
        }

        bool has_profile_focus = false;
        for (ui u = 0; u < profile_query_vertices.size(); ++u) {
            if (profile_query_vertices[u]) {
                has_profile_focus = true;
                break;
            }
        }

        struct ExtensionDepthRow {
            ui depth = 0;
            ui u = 0;
            BranchProfileCounter counter;
        };
        vector<ExtensionDepthRow> extension_rows;
        extension_rows.reserve(extension_depth_profiles.size());
        for (map<pair<ui, ui>, BranchProfileCounter>::const_iterator it =
            extension_depth_profiles.begin(); it != extension_depth_profiles.end(); ++it) {
            if (has_profile_focus && !profile_query_vertices[it->first.second]) {
                continue;
            }
            ExtensionDepthRow row;
            row.depth = it->first.first;
            row.u = it->first.second;
            row.counter = it->second;
            extension_rows.push_back(row);
        }
        std::sort(extension_rows.begin(), extension_rows.end(),
            [](const ExtensionDepthRow &a, const ExtensionDepthRow &b) {
                if (a.depth != b.depth) return a.depth < b.depth;
                unsigned long long ab = branchProfileBranchTotal(a.counter);
                unsigned long long bb = branchProfileBranchTotal(b.counter);
                if (ab != bb) return ab > bb;
                if (a.counter.dfs_calls != b.counter.dfs_calls) {
                    return a.counter.dfs_calls > b.counter.dfs_calls;
                }
                return a.u < b.u;
            });

        printf("Extension-depth rows%s:\n", has_profile_focus ? " (focused)" : "");
        if (branch_profile_path_depth == 0) {
            printf("  <disabled>\n");
        }
        else if (extension_rows.empty()) {
            printf("  <none>\n");
        }
        for (ui depth = 1; depth <= branch_profile_path_depth; ++depth) {
            ui depth_printed = 0;
            bool depth_header_printed = false;
            for (vector<ExtensionDepthRow>::const_iterator it = extension_rows.begin();
                it != extension_rows.end() && depth_printed < branch_profile_top_n; ++it) {
                if (it->depth != depth) {
                    continue;
                }
                if (!depth_header_printed) {
                    printf("  Depth %u:\n", depth);
                    depth_header_printed = true;
                }
                printf("    u%u, ", it->u);
                printBranchProfileCounter(it->counter);
                printf("\n");
                depth_printed++;
            }
            if (!depth_header_printed && !extension_rows.empty()) {
                printf("  Depth %u: <none>\n", depth);
            }
        }
        printf("-----------------------------------------------------------\n");
        fflush(stdout);
    }
#else
    void recordBranchProfileRootStart(ui, ui) {}
    void recordBranchProfileRootFinish(ui, ui) {}
    void recordBranchProfileDfsEnter() {}

    struct BranchProfileDfsScope {
        explicit BranchProfileDfsScope(MatchingSolver &) {}
    };

    void recordBranchProfilePrune() {}
    void recordBranchProfileAnswer() {}
    void recordBranchProfileSelectedEdge(ui, ui, size_t) {}
    void recordBranchProfileForcedZeroExclude(ui, ui) {}
    void recordBranchProfileIncludeConsidered(ui, ui) {}
    void recordBranchProfileIncludeTaken(ui, ui) {}
    void recordBranchProfileIncludeCostPruned(ui, ui) {}
    void recordBranchProfileExclude(ui, ui, bool) {}
#endif

    // init q_matrix and q_neighbors
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
    private:
        MatchingSolver &solver;

        struct LabelFreq {
            LabelID label = 0;
            ui count = 0;
        };

        struct QueryLabelReq {
            LabelID label = 0;
            ui bridge = 0;
            ui normal = 0;
        };

        struct BridgeArc {
            ui from = 0;
            ui to = 0;
        };

        struct BridgeNbr {
            ui to = 0;
            ui support_arc = 0;
        };

        // static graph/profile cache
        vector<LabelID> query_label;
        vector<LabelID> data_label;
        vector<ui> data_degree;

        // NLF filtering
        vector<vector<LabelFreq>> data_label_freqs;
        vector<vector<QueryLabelReq>> query_label_reqs;
        vector<vector<ui>> data_by_label;
        vector<ui> query_bridge_degree;

        // Bridge filtering
        vector<BridgeArc> bridge_arcs;
        vector<vector<BridgeNbr>> bridge_nbrs;
        vector<ui> bridge_support;
        queue<pair<ui, ui>> removed;

        // spoke filtering
        vector<vector<ui>>  spoke_adj;
        vector<int>         match_right;
        vector<ui>          seen_right;
        ui                  seen_token = 1;
        vector<char>        left_is_bridge;
#ifdef ENABLE_SPOKE_FILTERING
        queue<ui>           pending_spokes;
        vector<char>        queued_spoke;
#endif

    public:
        explicit CandidateFilter(MatchingSolver &solver)
            : solver(solver),
            spoke_adj(solver.qn, vector<ui>()),
            match_right(solver.max_g_deg, -1),
            seen_right(solver.max_g_deg, 0),
            left_is_bridge(solver.qn, 0)
#ifdef ENABLE_SPOKE_FILTERING
            , queued_spoke(solver.qn, 0)
#endif
        {}

        bool run()
        {
            if (!timed(&MatchingSolver::TimeStats::filter_bridge_time, [&] {
                buildBridgeIndex();
                return true;
            })) return false;

            if (!timed(&MatchingSolver::TimeStats::filter_nlf_time, [&] {
                return filterNLF();
            })) return false;

#if CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING
            if (!timed(&MatchingSolver::TimeStats::filter_bridge_time, [&] {
                return initBridgeSupports() && propagateFilterClosure();
            })) return false;
#endif

#ifdef ENABLE_SPOKE_FILTERING
            if (!timed(&MatchingSolver::TimeStats::filter_spoke_time, [&] {
                enqueueAllSpokeVertices();
                return propagateFilterClosure();
            })) return false;
#endif

            updateCandidateCount();
            return true;
        }

    private:
        template <typename Fn>
        bool timed(long long MatchingSolver::TimeStats::*field, Fn &&fn)
        {
#ifndef NDEBUG
            Timer t;
#endif

            bool ok = fn();

#ifndef NDEBUG
            solver.stats.*field += t.elapsed();
#endif

            return ok;
        }

        void updateCandidateCount()
        {
            solver.stats.filter_candidate_count = 0;
            for (ui u = 0; u < solver.qn; ++u) {
                solver.stats.filter_candidate_count += (ui)solver.candidates[u].size();
            }
        }

        void buildBridgeIndex()
        {
            bridge_arcs.clear();
            bridge_nbrs.assign(solver.qn, vector<BridgeNbr>());
            solver.q_neighbor_is_bridge.assign(solver.qn, vector<char>());

            for (ui u = 0; u < solver.qn; ++u) {
                solver.q_neighbor_is_bridge[u].assign(solver.q_degree[u], 0);
            }

#if CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING
            vector<int> dfn(solver.qn, 0);
            vector<int> low(solver.qn, 0);
            int tim = 0;

            for (ui u = 0; u < solver.qn; ++u) {
                if (dfn[u] == 0) {
                    tarjan(u, solver.qn, dfn, low, tim);
                }
            }
#endif
        }

        void tarjan(ui u, ui parent, vector<int> &dfn, vector<int> &low, int &time)
        {
            dfn[u] = low[u] = ++time;

            for (ui v : solver.q_neighbors[u]) {
                if (dfn[v] == 0) {
                    tarjan(v, u, dfn, low, time);
                    low[u] = std::min(low[u], low[v]);
                    if (low[v] > dfn[u]) {
                        addBridge(u, v);
                    }
                }
                else if (v != parent) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
        }

        ui addBridgeArc(ui from, ui to)
        {
            bridge_arcs.push_back({ from, to });
            return (ui)bridge_arcs.size() - 1;
        }

        void addBridge(ui a, ui b)
        {
            ui ab = addBridgeArc(a, b);
            ui ba = addBridgeArc(b, a);
            bridge_nbrs[a].push_back({ b, ba });
            bridge_nbrs[b].push_back({ a, ab });
            markBridgeNeighbor(a, b);
            markBridgeNeighbor(b, a);
        }

        void markBridgeNeighbor(ui from, ui to)
        {
            bool marked = false;
            const vector<ui> &neighbors = solver.q_neighbors[from];
            for (ui i = 0; i < solver.q_degree[from]; ++i) {
                if (neighbors[i] == to) {
                    solver.q_neighbor_is_bridge[from][i] = 1;
                    marked = true;
                    break;
                }
            }
            assert(marked);
            (void)marked;
        }

        void buildVertexCache()
        {
            query_label.assign(solver.qn, 0);
            for (ui u = 0; u < solver.qn; ++u) {
                query_label[u] = solver.query_graph->getVertexLabel(u);
            }

            data_label.assign(solver.gn, 0);
            data_degree.assign(solver.gn, 0);
            for (ui v = 0; v < solver.gn; ++v) {
                data_label[v] = solver.data_graph->getVertexLabel(v);
                data_degree[v] = solver.data_graph->getVertexDegree(v);
            }
        }

        void buildQueryLabelReqs()
        {
            query_label_reqs.assign(solver.qn, vector<QueryLabelReq>());
            query_bridge_degree.assign(solver.qn, 0);
            vector<ui> bridge_counts(solver.label_count, 0);
            vector<ui> non_bridge_counts(solver.label_count, 0);
            vector<ui> touched_labels;

            for (ui u = 0; u < solver.qn; ++u) {
                query_label_reqs[u].reserve(std::min(solver.q_degree[u], solver.label_count));
                touched_labels.clear();

                for (ui i = 0; i < solver.q_degree[u]; ++i) {
                    ui u1 = solver.q_neighbors[u][i];
                    LabelID label = query_label[u1];
                    assert(label >= 0 && (ui)label < solver.label_count);
                    ui label_idx = (ui)label;

                    if (bridge_counts[label_idx] == 0 && non_bridge_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    if (solver.q_neighbor_is_bridge[u][i]) {
                        query_bridge_degree[u]++;
                        bridge_counts[label_idx]++;
                    }
                    else {
                        non_bridge_counts[label_idx]++;
                    }
                }

                sort(touched_labels.begin(), touched_labels.end());
                for (ui label : touched_labels) {
                    ui bridge_count = bridge_counts[label];
                    ui non_bridge_count = non_bridge_counts[label];
                    query_label_reqs[u].push_back({
                        (LabelID)label, bridge_count, non_bridge_count
                    });
                    bridge_counts[label] = 0;
                    non_bridge_counts[label] = 0;
                }
            }
        }

        void buildDataLabelFreqs()
        {
            data_label_freqs.assign(solver.gn, vector<LabelFreq>());
            data_by_label.assign(solver.label_count, vector<ui>());
            vector<ui> label_counts(solver.label_count, 0);
            vector<ui> touched_labels;

            for (ui v = 0; v < solver.gn; ++v) {
                ui deg = 0;
                const ui *neighbors = solver.data_graph->getVertexNeighbors(v, deg);
                data_label_freqs[v].reserve(std::min(deg, solver.label_count));
                touched_labels.clear();

                LabelID vertex_label = data_label[v];
                assert(vertex_label >= 0 && (ui)vertex_label < solver.label_count);
                if (vertex_label >= 0 && (ui)vertex_label < solver.label_count) {
                    data_by_label[(ui)vertex_label].push_back(v);
                }

                for (ui i = 0; i < deg; ++i) {
                    LabelID label = data_label[neighbors[i]];
                    assert(label >= 0 && (ui)label < solver.label_count);
                    ui label_idx = (ui)label;
                    if (label_counts[label_idx] == 0) {
                        touched_labels.push_back(label_idx);
                    }
                    label_counts[label_idx]++;
                }

                sort(touched_labels.begin(), touched_labels.end());
                for (ui label : touched_labels) {
                    data_label_freqs[v].push_back({ (LabelID)label, label_counts[label] });
                    label_counts[label] = 0;
                }
            }
        }

        ui nlfDiff(ui u, ui v) const
        {
            ui diff = 0;
            const vector<QueryLabelReq> &query_reqs = query_label_reqs[u];
            const vector<LabelFreq> &data_freqs = data_label_freqs[v];
            size_t data_idx = 0;

            for (const auto &req : query_reqs) {
                while (data_idx < data_freqs.size() && data_freqs[data_idx].label < req.label) {
                    data_idx++;
                }

                ui data_count = 0;
                if (data_idx < data_freqs.size() && data_freqs[data_idx].label == req.label) {
                    data_count = data_freqs[data_idx].count;
                }

                ui bridge_need = req.bridge;
                if (bridge_need > data_count) {
                    return solver.threshold + 1;
                }

                ui non_bridge_need = req.normal;
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
            buildVertexCache();
            buildQueryLabelReqs();
            buildDataLabelFreqs();

            for (ui u = 0; u < solver.qn; ++u) {
                LabelID lu = query_label[u];
                assert(lu >= 0 && (ui)lu < data_by_label.size());
                if (lu < 0 || (ui)lu >= data_by_label.size()) return false;

                for (ui v : data_by_label[(ui)lu]) {
                    if (query_bridge_degree[u] > data_degree[v]) continue;
                    if (solver.q_degree[u] > data_degree[v] + solver.threshold) continue;
                    if (nlfDiff(u, v) > solver.threshold) continue;
                    solver.candidates[u].insert(v);
                }
                if (solver.candidates[u].empty()) return false;
            }
            return true;
        }

        ui &support(ui arc_id, ui v)
        {
            return bridge_support[(size_t)arc_id * solver.gn + v];
        }

        bool pruneCandidate(ui u, ui v)
        {
            if (!solver.candidates[u].contains(v)) {
                return true;
            }

            solver.candidates[u].remove(v);
            if (solver.candidates[u].empty()) {
                return false;
            }
            removed.push({ u, v });
#ifdef ENABLE_SPOKE_FILTERING
            for (ui nbr_u : solver.q_neighbors[u]) {
                enqueueSpokeVertex(nbr_u);
            }
#endif
            return true;
        }

        bool initBridgeSupports()
        {
            bridge_support.assign((size_t)bridge_arcs.size() * solver.gn, 0);
            vector<pair<ui, ui>> zero_support_candidates;

            for (ui arc_id = 0; arc_id < (ui)bridge_arcs.size(); ++arc_id) {
                const BridgeArc &arc = bridge_arcs[arc_id];
                for (ui v : solver.candidates[arc.from]) {
                    ui deg = 0;
                    const ui *nbrs = solver.data_graph->getVertexNeighbors(v, deg);
                    ui support = 0;
                    for (ui i = 0; i < deg; ++i) {
                        if (solver.candidates[arc.to].contains(nbrs[i])) {
                            support++;
                        }
                    }
                    this->support(arc_id, v) = support;
                    if (support == 0) {
                        zero_support_candidates.push_back({ arc.from, v });
                    }
                }
            }

            for (const auto &candidate : zero_support_candidates) {
                if (!pruneCandidate(candidate.first, candidate.second)) {
                    return false;
                }
            }
            return true;
        }

        bool propagateBridgeRemovals()
        {
            while (!removed.empty()) {
                ui removed_u = removed.front().first;
                ui removed_v = removed.front().second;
                removed.pop();

                ui deg = 0;
                const ui *nbrs = solver.data_graph->getVertexNeighbors(removed_v, deg);

                for (const BridgeNbr &bridge_nbr : bridge_nbrs[removed_u]) {
                    ui affected_u = bridge_nbr.to;
                    ui arc_id = bridge_nbr.support_arc;
                    for (ui i = 0; i < deg; ++i) {
                        ui v = nbrs[i];
                        if (!solver.candidates[affected_u].contains(v)) {
                            continue;
                        }
                        ui &support_count = support(arc_id, v);
                        if (support_count == 0) {
                            continue;
                        }
                        support_count--;
                        if (support_count == 0 &&
                            !pruneCandidate(affected_u, v)) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        bool propagateFilterClosure()
        {
            while (true) {
                if (!propagateBridgeRemovals()) {
                    return false;
                }
#ifdef ENABLE_SPOKE_FILTERING
                if (pending_spokes.empty()) {
                    break;
                }

                ui u = pending_spokes.front();
                pending_spokes.pop();
                queued_spoke[u] = 0;
                if (!processSpokeVertex(u)) {
                    return false;
                }
#else
                break;
#endif
            }
            return true;
        }

        bool augmentSpoke(ui left_idx)
        {
            for (ui right_idx : spoke_adj[left_idx]) {
                if (seen_right[right_idx] == seen_token) continue;
                seen_right[right_idx] = seen_token;
                if (match_right[right_idx] < 0 ||
                    augmentSpoke((ui)match_right[right_idx])) {
                    match_right[right_idx] = (int)left_idx;
                    return true;
                }
            }
            return false;
        }

        bool tryAugmentSpoke(ui left_idx)
        {
            seen_token++;
            if (seen_token == 0) {
                std::fill(seen_right.begin(), seen_right.end(), 0);
                seen_token = 1;
            }
            return augmentSpoke(left_idx);
        }

        void buildSpokeAdj(ui u, ui v, ui &deg_u, ui &deg_v)
        {
            const vector<ui> &u_neighbors = solver.q_neighbors[u];
            deg_u = solver.q_degree[u];
            const ui *v_neighbors = solver.data_graph->getVertexNeighbors(v, deg_v);

            for (ui i = 0; i < deg_u; ++i) {
                spoke_adj[i].clear();
                ui u1 = u_neighbors[i];
                left_is_bridge[i] = solver.q_neighbor_is_bridge[u][i];
                for (ui j = 0; j < deg_v; ++j) {
                    ui v1 = v_neighbors[j];
                    if (solver.candidates[u1].contains(v1)) {
                        spoke_adj[i].push_back(j);
                    }
                }
            }
        }

        bool spokeFeasible(ui u, ui v, ui budget)
        {
            ui deg_u = 0;
            ui deg_v = 0;
            buildSpokeAdj(u, v, deg_u, deg_v);

            std::fill(match_right.begin(), match_right.begin() + deg_v, -1);

            ui optional_count = 0;
            for (ui i = 0; i < deg_u; ++i) {
                if (left_is_bridge[i]) {
                    if (!tryAugmentSpoke(i)) {
                        return false;
                    }
                }
                else {
                    optional_count++;
                }
            }

            if (optional_count <= budget) {
                return true;
            }

            ui required = optional_count - budget;
            ui matched = 0;
            ui processed = 0;

            for (ui i = 0; i < deg_u; ++i) {
                if (left_is_bridge[i]) {
                    continue;
                }

                processed++;
                if (tryAugmentSpoke(i)) {
                    matched++;
                    if (matched >= required) {
                        return true;
                    }
                }

                ui remaining = optional_count - processed;
                if (matched + remaining < required) {
                    return false;
                }
            }

            return matched >= required;
        }

        ui maxMissingIncidentEdges(ui u) const
        {
            if (solver.q_degree[u] == 0) return 0;
            return std::min(solver.threshold, solver.q_degree[u] - 1);
        }

#ifdef ENABLE_SPOKE_FILTERING
        void enqueueSpokeVertex(ui u)
        {
            if (queued_spoke[u]) {
                return;
            }
            pending_spokes.push(u);
            queued_spoke[u] = 1;
        }

        void enqueueAllSpokeVertices()
        {
            for (ui u = 0; u < solver.qn; ++u) {
                enqueueSpokeVertex(u);
            }
        }

        bool processSpokeVertex(ui u)
        {
            ui budget = maxMissingIncidentEdges(u);

            vector<ui> to_remove;
            for (ui v : solver.candidates[u]) {
                if (!spokeFeasible(u, v, budget)) {
                    to_remove.push_back(v);
                }
            }
            if (to_remove.empty()) return true;

            for (ui v : to_remove) {
                if (!pruneCandidate(u, v)) {
                    return false;
                }
            }
            return true;
        }
#endif
    };

    bool runCandidateFiltering()
    {
        return CandidateFilter(*this).run();
    }

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
    void initDataLabelDegreeIndex()
    {
        data_label_degrees.assign(gn, vector<DataLabelDegreeCount>());
        vector<ui> label_counts(label_count, 0);

        for (ui v = 0; v < gn; ++v) {
            ui deg = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(v, deg);
            data_label_degrees[v].reserve(std::min(deg, label_count));

            for (ui i = 0; i < deg; ++i) {
                LabelID label = data_graph->getVertexLabel(neighbors[i]);
                assert(label >= 0 && (ui)label < label_count);
                label_counts[(ui)label]++;
            }

            for (ui label = 0; label < label_count; ++label) {
                ui count = label_counts[label];
                if (count != 0) {
                    data_label_degrees[v].push_back({ (LabelID)label, count });
                    label_counts[label] = 0;
                }
            }
        }
    }

    ui dataLabelDegree(ui v, LabelID label) const
    {
        if (v >= gn || label < 0 || (ui)label >= label_count) {
            return 0;
        }

        const vector<DataLabelDegreeCount> &counts = data_label_degrees[v];
        auto it = std::lower_bound(counts.begin(), counts.end(), label,
            [](const DataLabelDegreeCount &item, LabelID target_label) {
                return item.label < target_label;
            });
        if (it == counts.end() || it->label != label) {
            return 0;
        }
        return it->count;
    }

    double initialTopkDecayRankSupport(ui u, ui anchor)
    {
        if (u >= qn || anchor >= qn || mapped_q[anchor] == -1) {
            return 0.0;
        }

        LabelID label = query_graph->getVertexLabel(u);
        ui label_degree = dataLabelDegree((ui)mapped_q[anchor], label);
        ui candidate_count = (ui)candidates[u].size();
        return (double)std::min(candidate_count, label_degree);
    }
#endif

#if CDE_EDGE_IE_FIXED_ORDER
    struct FixedEdgePriorityEntry {
        ui u = 0;
        ui anchor = 0;
        unsigned long long pair_support = 0;
        ui u_candidate_count = 0;
        ui anchor_candidate_count = 0;
    };

    unsigned long long countStaticCandidateEdgePairs(ui u, ui u1) const
    {
        unsigned long long count = 0;
        for (ui v : candidates[u]) {
            ui degree = 0;
            const ui *neighbors = data_graph->getVertexNeighbors(v, degree);
            for (ui i = 0; i < degree; ++i) {
                if (candidates[u1].contains(neighbors[i])) {
                    count++;
                }
            }
        }
        return count;
    }

    void initFixedEdgePriorities()
    {
        const ui invalid_priority = std::numeric_limits<ui>::max();
        fixed_edge_priority.assign(qn, vector<ui>(qn, invalid_priority));

        vector<FixedEdgePriorityEntry> entries;
        size_t directed_edge_count = 0;
        for (ui u = 0; u < qn; ++u) {
            directed_edge_count += q_neighbors[u].size();
        }
        entries.reserve(directed_edge_count);

        for (ui u = 0; u < qn; ++u) {
            for (ui u1 : q_neighbors[u]) {
                if (u >= u1) {
                    continue;
                }

                ui u_candidate_count = (ui)candidates[u].size();
                ui u1_candidate_count = (ui)candidates[u1].size();
                ui scan_u = u_candidate_count <= u1_candidate_count ? u : u1;
                ui target_u = scan_u == u ? u1 : u;
                unsigned long long pair_support =
                    countStaticCandidateEdgePairs(scan_u, target_u);

                FixedEdgePriorityEntry forward;
                forward.u = u;
                forward.anchor = u1;
                forward.pair_support = pair_support;
                forward.u_candidate_count = u_candidate_count;
                forward.anchor_candidate_count = u1_candidate_count;
                entries.push_back(forward);

                FixedEdgePriorityEntry reverse;
                reverse.u = u1;
                reverse.anchor = u;
                reverse.pair_support = pair_support;
                reverse.u_candidate_count = u1_candidate_count;
                reverse.anchor_candidate_count = u_candidate_count;
                entries.push_back(reverse);
            }
        }

        std::sort(entries.begin(), entries.end(),
            [&](const FixedEdgePriorityEntry &lhs,
                const FixedEdgePriorityEntry &rhs) {
                // pair_support / |Cand(anchor)| estimates the support seen after
                // fixing one candidate for the anchor. Compare the ratios
                // exactly so the priority is deterministic.
                __uint128_t lhs_scaled =
                    (__uint128_t)lhs.pair_support *
                    std::max((ui)1, rhs.anchor_candidate_count);
                __uint128_t rhs_scaled =
                    (__uint128_t)rhs.pair_support *
                    std::max((ui)1, lhs.anchor_candidate_count);
                if (lhs_scaled != rhs_scaled) {
                    return lhs_scaled < rhs_scaled;
                }
                if (lhs.u_candidate_count != rhs.u_candidate_count) {
                    return lhs.u_candidate_count < rhs.u_candidate_count;
                }
                if (lhs.pair_support != rhs.pair_support) {
                    return lhs.pair_support < rhs.pair_support;
                }
                if (q_degree[lhs.u] != q_degree[rhs.u]) {
                    return q_degree[lhs.u] > q_degree[rhs.u];
                }
                if (q_degree[lhs.anchor] != q_degree[rhs.anchor]) {
                    return q_degree[lhs.anchor] > q_degree[rhs.anchor];
                }
                if (lhs.u != rhs.u) {
                    return lhs.u < rhs.u;
                }
                return lhs.anchor < rhs.anchor;
            });

        assert(entries.size() <= (size_t)std::numeric_limits<ui>::max());
        for (size_t rank = 0; rank < entries.size(); ++rank) {
            const FixedEdgePriorityEntry &entry = entries[rank];
            fixed_edge_priority[entry.u][entry.anchor] = (ui)rank;
        }
    }
#endif
    // ========================================================================

    // ========================================================================
    // Branching
    // ========================================================================
    struct BranchSelector {
    private:
        MatchingSolver &solver;

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

        bool collectTopActiveEdges(const vector<ui> &component_frontier, ui max_count,
            vector<ActiveEdge> &top_edges, EdgeScoreCache *edge_score_cache = nullptr,
            const vector<char> *skip_query_vertices = nullptr) const
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
                if (shouldSkipQueryVertex(u, skip_query_vertices)) {
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

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
            selectTopActiveEdgesWithDecay(max_count, top_edges);
#else
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
#endif
            return true;
        }

        bool restrictTopEdgesToCoveredComponent(vector<ActiveEdge> &top_edges,
            EdgeScoreCache &edge_score_cache, vector<int> &component_id,
            vector<vector<ui>> &component_frontiers, vector<char> &component_checked,
            vector<ActiveEdge> &component_edges, vector<ActiveEdge> &best_component_edges,
            double &best_support_sum, bool &has_zero_support_component,
            const vector<char> *skip_query_vertices = nullptr) const
        {
            best_component_edges.clear();
            best_support_sum = std::numeric_limits<double>::max();
            has_zero_support_component = false;
            if (top_edges.empty()) {
                return false;
            }

            size_t component_count = labelUnmatchedComponents(component_id,
                component_frontiers, skip_query_vertices);
            component_checked.assign(component_count, 0);

            for (const ActiveEdge &edge : top_edges) {
                if (edge.u >= component_id.size() || component_id[edge.u] < 0) {
                    continue;
                }

                size_t id = (size_t)component_id[edge.u];
                if (component_checked[id]) {
                    continue;
                }
                component_checked[id] = 1;

                const vector<ui> &component_frontier = component_frontiers[id];
                if (component_frontier.empty()) {
                    continue;
                }

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
                component_edges.clear();
                bool covered = true;
                for (ui component_u : component_frontier) {
                    const vector<ActiveEdge> &cached_edges =
                        cachedActiveEdgesForVertex(component_u, edge_score_cache);
                    for (const ActiveEdge &component_edge : cached_edges) {
                        if (!containsActiveEdge(top_edges, component_edge)) {
                            covered = false;
                            break;
                        }
                    }
                    if (!covered) {
                        break;
                    }
                }

                if (!covered) {
                    continue;
                }

                double support_sum = 0.0;
                for (const ActiveEdge &top_edge : top_edges) {
                    if (top_edge.u < component_id.size() &&
                        component_id[top_edge.u] == (int)id) {
                        component_edges.push_back(top_edge);
                        support_sum += top_edge.rank_support;
                    }
                }
                if (component_edges.empty()) {
                    continue;
                }
#else
                if (!collectTopActiveEdges(component_frontier,
                    std::numeric_limits<ui>::max(), component_edges, &edge_score_cache,
                    skip_query_vertices)) {
                    continue;
                }

                double support_sum = 0.0;
                for (const ActiveEdge &component_edge : component_edges) {
                    support_sum += component_edge.anchor_support;
                }
                if (support_sum == 0) {
                    has_zero_support_component = true;
                    return false;
                }

                bool covered = true;
                for (const ActiveEdge &component_edge : component_edges) {
                    if (!containsActiveEdge(top_edges, component_edge)) {
                        covered = false;
                        break;
                    }
                }

                if (covered) {
                    if (best_component_edges.empty() || support_sum < best_support_sum) {
                        best_component_edges = component_edges;
                        best_support_sum = support_sum;
                    }
                }
#endif
#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
                if (best_component_edges.empty() || support_sum < best_support_sum) {
                    best_component_edges = component_edges;
                    best_support_sum = support_sum;
                }
#endif
            }

            if (best_component_edges.empty()) {
                return false;
            }

            top_edges = best_component_edges;
            return true;
        }

    private:
        bool shouldSkipQueryVertex(ui u, const vector<char> *skip_query_vertices) const
        {
            return skip_query_vertices != nullptr &&
                u < skip_query_vertices->size() && (*skip_query_vertices)[u];
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
            if (u >= solver.qn || anchor >= solver.qn ||
                solver.mapped_q[anchor] == -1 || solver.excluded_edges[u][anchor]) {
                return 0;
            }

            if (solver.support_dirty[u][anchor]) {
#ifndef NDEBUG
                Timer t_score;
#endif
                solver.recordSupportSnapshot(u, anchor);
                solver.support[u][anchor] = solver.calEdgeSupport(u, anchor, [](ui) {});
                solver.support_dirty[u][anchor] = 0;
#ifndef NDEBUG
                long long elapsed = t_score.elapsed();
                solver.stats.frontier_score_anchor_support_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
#endif
            }
#ifndef NDEBUG
            else {
                Timer t_score;
                ui exact_count = solver.calEdgeSupport(u, anchor, [](ui) {});
                long long elapsed = t_score.elapsed();
                solver.stats.frontier_score_anchor_support_time += elapsed;
                solver.stats.frontier_score_time += elapsed;
                assert(solver.support[u][anchor] >= exact_count);
            }
#endif
            return solver.support[u][anchor];
        }

        bool isBetterActiveEdge(const ActiveEdge &lhs, const ActiveEdge &rhs) const
        {
#if CDE_EDGE_IE_FIXED_ORDER
            ui lhs_priority = solver.fixed_edge_priority[lhs.u][lhs.anchor];
            ui rhs_priority = solver.fixed_edge_priority[rhs.u][rhs.anchor];
            if (lhs_priority != rhs_priority) {
                return lhs_priority < rhs_priority;
            }
            if (lhs.u != rhs.u) {
                return lhs.u < rhs.u;
            }
            return lhs.anchor < rhs.anchor;
#elif CDE_EDGE_IE_TOPK_SUPPORT_DECAY
#ifndef NDEBUG
            char lhs_preferred = solver.preferred_query_vertices.empty() ? 0 : solver.preferred_query_vertices[lhs.u];
            char rhs_preferred = solver.preferred_query_vertices.empty() ? 0 : solver.preferred_query_vertices[rhs.u];
            if (lhs_preferred != rhs_preferred) {
                return lhs_preferred > rhs_preferred;
            }
            char lhs_delayed = solver.delayed_query_vertices.empty() ? 0 : solver.delayed_query_vertices[lhs.u];
            char rhs_delayed = solver.delayed_query_vertices.empty() ? 0 : solver.delayed_query_vertices[rhs.u];
            if (lhs_delayed != rhs_delayed) {
                return lhs_delayed < rhs_delayed;
            }
#endif
            double lhs_scaled =
                lhs.rank_support * (double)std::max((ui)1, rhs.live_anchor_count);
            double rhs_scaled =
                rhs.rank_support * (double)std::max((ui)1, lhs.live_anchor_count);
            double scale = std::max(1.0,
                std::max(std::fabs(lhs_scaled), std::fabs(rhs_scaled)));
            if (std::fabs(lhs_scaled - rhs_scaled) > 1e-12 * scale) {
                return lhs_scaled < rhs_scaled;
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
#else
#ifndef NDEBUG
            char lhs_preferred = solver.preferred_query_vertices.empty() ? 0 : solver.preferred_query_vertices[lhs.u];
            char rhs_preferred = solver.preferred_query_vertices.empty() ? 0 : solver.preferred_query_vertices[rhs.u];
            if (lhs_preferred != rhs_preferred) {
                return lhs_preferred > rhs_preferred;
            }
            char lhs_delayed = solver.delayed_query_vertices.empty() ? 0 : solver.delayed_query_vertices[lhs.u];
            char rhs_delayed = solver.delayed_query_vertices.empty() ? 0 : solver.delayed_query_vertices[rhs.u];
            if (lhs_delayed != rhs_delayed) {
                return lhs_delayed < rhs_delayed;
            }
#endif
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
#endif
        }

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        void selectTopActiveEdgesWithDecay(ui max_count, vector<ActiveEdge> &top_edges) const
        {
            size_t selected_limit = top_edges.size();
            if ((size_t)max_count < selected_limit) {
                selected_limit = (size_t)max_count;
            }

            const double gamma = (double)CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA;
            for (size_t selected_idx = 0; selected_idx < selected_limit; ++selected_idx) {
                size_t best_idx = selected_idx;
                for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
                    if (isBetterActiveEdge(top_edges[i], top_edges[best_idx])) {
                        best_idx = i;
                    }
                }

                if (best_idx != selected_idx) {
                    std::swap(top_edges[selected_idx], top_edges[best_idx]);
                }

                const ActiveEdge &selected = top_edges[selected_idx];
                double candidate_count = (double)solver.candidates[selected.u].size();
                if (candidate_count <= 0.0) {
                    continue;
                }

                double factor = 1.0 - gamma * selected.rank_support / candidate_count;
                if (factor < 0.0) {
                    factor = 0.0;
                }
                else if (factor > 1.0) {
                    factor = 1.0;
                }

                if (factor == 1.0) {
                    continue;
                }

                for (size_t i = selected_idx + 1; i < top_edges.size(); ++i) {
                    if (top_edges[i].u == selected.u) {
                        top_edges[i].rank_support *= factor;
                    }
                }
            }

            top_edges.resize(selected_limit);
        }
#endif

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
#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
                edge.rank_support = solver.initialTopkDecayRankSupport(u, anchor);
#else
                edge.anchor_support = countEdgeSupport(u, anchor);
#endif
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

        bool containsActiveEdge(const vector<ActiveEdge> &edges, const ActiveEdge &target) const
        {
            for (const ActiveEdge &edge : edges) {
                if (edge.u == target.u && edge.anchor == target.anchor) {
                    return true;
                }
            }
            return false;
        }

        size_t labelUnmatchedComponents(vector<int> &component_id,
            vector<vector<ui>> &component_frontiers,
            const vector<char> *skip_query_vertices) const
        {
            // Label every unmatched component in one graph sweep, then reuse the
            // labels for all top edges in this DFS state.
            component_id.assign(solver.qn, -1);
            for (vector<ui> &frontier : component_frontiers) {
                frontier.clear();
            }

            size_t component_count = 0;
            queue<ui> q;
            for (ui start = 0; start < solver.qn; ++start) {
                if (solver.mapped_q[start] != -1 ||
                    shouldSkipQueryVertex(start, skip_query_vertices) ||
                    component_id[start] != -1) {
                    continue;
                }

                int id = (int)component_count++;
                if ((size_t)id == component_frontiers.size()) {
                    component_frontiers.emplace_back();
                }
                component_id[start] = id;
                q.push(start);

                while (!q.empty()) {
                    ui curr = q.front();
                    q.pop();

                    if (solver.frontier_pos[curr] != -1) {
                        component_frontiers[(size_t)id].push_back(curr);
                    }

                    for (ui nbr : solver.q_neighbors[curr]) {
                        if (solver.mapped_q[nbr] == -1 &&
                            !shouldSkipQueryVertex(nbr, skip_query_vertices) &&
                            component_id[nbr] == -1) {
                            component_id[nbr] = id;
                            q.push(nbr);
                        }
                    }
                }
            }
            return component_count;
        }

    };

    ui countAnchorsByMark(ui u, ui selected_anchor, const vector<ui> &candidate_vertices, vector<ui> &anchor_counts)
    {
        anchor_counts.assign(candidate_vertices.size(), 0);
        if(++data_vertex_mark_token == 0) {
            std::fill(data_vertex_mark.begin(), data_vertex_mark.end(), 0);
            data_vertex_mark_token = 1;
        }
        for (ui i = 0; i < candidate_vertices.size(); ++i) {
            ui v = candidate_vertices[i];
            assert(v < gn);
            data_vertex_mark[v] = data_vertex_mark_token;
            data_vertex_mark_pos[v] = i;
        }
        ui token = data_vertex_mark_token;
        ui anchor_num = 0;

        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;

            anchor_num++;
            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (data_vertex_mark[v] == token) {
                    anchor_counts[data_vertex_mark_pos[v]]++;
                }
            }
        }

        return anchor_num;
    }

    ui countAnchorsByHasEdge(ui u, ui selected_anchor, const vector<ui> &candidate_vertices, vector<ui> &anchor_counts) const
    {
        anchor_counts.assign(candidate_vertices.size(), 0);
        ui anchor_num = 0;
        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
            anchor_num++;
        }

        for (ui i = 0; i < candidate_vertices.size(); ++i) {
            ui v = candidate_vertices[i];
            for (ui anchor : q_neighbors[u]) {
                if (anchor == selected_anchor) continue;
                if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;
                if (data_graph->hasEdge(v, (ui)mapped_q[anchor])) {
                    anchor_counts[i]++;
                }
            }
        }

        return anchor_num;
    }

    // heuristic to decide whether to use mark-based counting or has-edge checking.
    bool useMarkForAnchorCount(ui u, ui selected_anchor, size_t vertex_count) const
    {
        auto binary_search_cost = [](ui degree) -> size_t {
            size_t cost = 1;
            for (; degree > 1; degree >>= 1) ++cost;
            return cost;
        };

        size_t mark_cost = vertex_count * 2;
        size_t has_edge_cost = 0;

        for (ui anchor : q_neighbors[u]) {
            if (anchor == selected_anchor) continue;
            if (mapped_q[anchor] == -1 || excluded_edges[u][anchor]) continue;

            ui mapped_anchor = (ui)mapped_q[anchor];
            ui degree = data_graph->getVertexDegree(mapped_anchor);

            mark_cost += degree;
            has_edge_cost += vertex_count * binary_search_cost(degree);
        }

        return mark_cost <= has_edge_cost;
    }

    // count u's candidate vertices that are adjacent to each anchor's mapped data vertex
    // and store the counts in anchor_counts
    // return the number of anchors in query graph excluding selected_anchor and excluded ones.
    ui countAnchors(ui u, ui ua, const vector<ui> &candidate_vertices, vector<ui> &anchor_counts)
    {
        if (candidate_vertices.empty()) return 0;
        if (useMarkForAnchorCount(u, ua, candidate_vertices.size())) {
            return countAnchorsByMark(u, ua, candidate_vertices, anchor_counts);
        }
        return countAnchorsByHasEdge(u, ua, candidate_vertices, anchor_counts);
    }

    // ========================================================================
    // Maintain Frontier and Anchor Support
    // ========================================================================
    // calculate support[u][anchor] or calculate all candidates of u by (u, anchor)
    template <typename Visitor>
    ui calEdgeSupport(ui u, ui anchor, Visitor visit) const
    {
        if (u >= qn || anchor >= qn || mapped_q[anchor] == -1 || excluded_edges[u][anchor]) {
            return 0;
        }

        ui count = 0;
        ui deg = 0;
        const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
        for (ui i = 0; i < deg; ++i) {
            ui v = nbrs[i];
            if (!candidates[u].contains(v)) continue;
            if (mapped_g[v] != -1) continue;
            if (excluded_cands[u].contains(v)) continue;
            count++;
            visit(v);
        }
        return count;
    }

    inline void recordSupportSnapshot(ui u, ui anchor)
    {
        if (active_support_snapshots == nullptr) {
            return;
        }
        SupportSnapshot snapshot;
        snapshot.u = u;
        snapshot.anchor = anchor;
        snapshot.value = support[u][anchor];
        snapshot.dirty = support_dirty[u][anchor];
        active_support_snapshots->push_back(snapshot);
    }

    inline void markSupportDirty(ui u, ui anchor)
    {
        if (u >= qn || anchor >= qn || support_dirty[u][anchor]) {
            return;
        }
        recordSupportSnapshot(u, anchor);
        support_dirty[u][anchor] = 1;
    }

    inline void markLiveAnchorSupportDirty(ui u)
    {
        if (u >= qn || mapped_q[u] != -1 || frontier_pos[u] == -1) {
            return;
        }
        for (ui anchor : q_neighbors[u]) {
            if (mapped_q[anchor] != -1 && !excluded_edges[u][anchor]) {
                markSupportDirty(u, anchor);
            }
        }
    }

    // update active_frontier, frontier_pos, anchor_support
    inline void updateFrontierStatus(ui u)
    {
        bool should_be = (mapped_q[u] == -1 && anchor_count[u] > 0);
        bool is_in = (frontier_pos[u] != -1);

        if (should_be && !is_in) {
            // add u to frontier
            frontier_pos[u] = active_frontier.size();
            active_frontier.push_back(u);
        }
        else if (!should_be && is_in) {
            // remove u from frontier
            ui idx = frontier_pos[u];
            ui last_u = active_frontier.back();
            active_frontier[idx] = last_u;
            frontier_pos[last_u] = idx;
            active_frontier.pop_back();
            frontier_pos[u] = -1;
        }
    }

    // when u becomes matched/unmatched, update anchor_count, active_frontier and anchor_support
    void updateFrontier(ui u, bool matched)
    {
        if (matched) updateFrontierStatus(u);

        // update u's neighbors in frontier
        for (ui u1 : q_neighbors[u]) {
            if(excluded_edges[u][u1]) continue;

            if (matched) {
                anchor_count[u1]++;
                updateFrontierStatus(u1);

                if (mapped_q[u1] == -1 && frontier_pos[u1] != -1) {
                    markSupportDirty(u1, u);
                }
            }
            else {
                anchor_count[u1]--;
                updateFrontierStatus(u1);
            }
        }

        if (!matched) updateFrontierStatus(u);
    }

    // mark (u, anchor) as excluded, update anchor_count[u], active_frontier and anchor_support
    void excludeFrontierEdge(ui u, ui anchor)
    {
        excluded_edges[u][anchor] = 1;
        excluded_edges[anchor][u] = 1;
        assert(anchor_count[u] > 0);
        anchor_count[u]--;
        if(anchor_count[u] == 0) updateFrontierStatus(u);
        markLiveAnchorSupportDirty(u);
    }

    // restore (u, anchor), update anchor_count[u], active_frontier and anchor_support
    void restoreFrontierEdge(ui u, ui anchor)
    {
        excluded_edges[u][anchor] = 0;
        excluded_edges[anchor][u] = 0;
        anchor_count[u]++;
        if(anchor_count[u] == 1) updateFrontierStatus(u);
    }

    struct SupportUndoScope {
        MatchingSolver &solver;
        vector<SupportSnapshot> *previous_snapshots;
        vector<SupportSnapshot> &snapshots;

        SupportUndoScope(MatchingSolver &solver, vector<SupportSnapshot> &snapshots)
            : solver(solver),
            previous_snapshots(solver.active_support_snapshots),
            snapshots(snapshots)
        {
            solver.active_support_snapshots = &snapshots;
        }

        ~SupportUndoScope()
        {
            for (auto it = snapshots.rbegin(); it != snapshots.rend(); ++it) {
                solver.support[it->u][it->anchor] = it->value;
                solver.support_dirty[it->u][it->anchor] = it->dirty;
            }
            snapshots.clear();
            solver.active_support_snapshots = previous_snapshots;
        }
    };
    // ========================================================================

    // ========================================================================
    // DFS Buffer
    // ========================================================================
    struct DfsBuffer {
        vector<ActiveEdge> top_edges;
        vector<ActiveEdge> component_edges;
        vector<ActiveEdge> best_component_edges;
        vector<ui> candidate_vertices;
        vector<ui> candidate_anchor_counts;
        vector<int> component_id;
        vector<vector<ui>> component_frontiers;
        vector<char> component_checked;
        vector<char> terminal_skip;
        vector<ui> terminal_vertices;
        vector<ui> active_terminal_vertices;
        vector<TerminalTailVertex> terminal_tail_vertices;
        EdgeScoreCache edge_score_cache;
        BranchSelector branch_selector;
        vector<SupportSnapshot> local_support_snapshots;
        explicit DfsBuffer(MatchingSolver &solver)
            : branch_selector(solver)
        {}

        // reserve space for dfs buffers
        void reserve(ui threshold, ui max_g_deg, ui qn, ui gn)
        {
            top_edges.reserve((size_t)threshold + 1);
            component_edges.reserve((size_t)threshold + 1);
            best_component_edges.reserve((size_t)threshold + 1);
            candidate_vertices.reserve(max_g_deg);
            candidate_anchor_counts.reserve(max_g_deg);
            component_id.reserve(qn);
            component_frontiers.reserve(qn);
            component_checked.reserve(qn);
            terminal_skip.reserve(qn);
            terminal_vertices.reserve(qn);
            active_terminal_vertices.reserve(qn);
            terminal_tail_vertices.reserve(qn);
            local_excluded_edges.reserve((size_t)threshold + 1);
            size_t cand_size = std::min((size_t)gn, ((size_t)threshold + 1) * (size_t)max_g_deg);
            local_excluded_cands.reserve(cand_size);
        }

        void clearLocal()
        {
            for (vector<ActiveEdge> &edges : edge_score_cache.active_edges_by_vertex) {
                edges.clear();
            }
            std::fill(edge_score_cache.active_edges_cached.begin(),
                edge_score_cache.active_edges_cached.end(), 0);
            component_id.clear();
            component_checked.clear();
            component_edges.clear();
            best_component_edges.clear();
            terminal_vertices.clear();
            active_terminal_vertices.clear();
            terminal_tail_vertices.clear();
            local_excluded_edges.clear();
            local_excluded_cands.clear();
            local_support_snapshots.clear();
        }

        void recordExcludedEdge(ui u, ui anchor)
        {
            local_excluded_edges.push_back({ u, anchor });
        }

        void recordExcludedCands(ui u, ui v)
        {
            local_excluded_cands.push_back({ u, v });
        }

        void restoreLocalChanges(MatchingSolver &solver)
        {
            for (auto &p : local_excluded_cands) {
                solver.excluded_cands[p.first].remove(p.second);
            }

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
    };

    vector<DfsBuffer> dfs_buffers;

    void reserveDfsBuffer(DfsBuffer &buf) const
    {
        buf.reserve(threshold, max_g_deg, qn, gn);
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

    // ========================================================================
    // Terminal-tail enumeration
    // ========================================================================
    bool isTerminalQueryVertex(ui u) const
    {
        if (u >= qn || mapped_q[u] != -1 || anchor_count[u] == 0) {
            return false;
        }

        for (ui nbr : q_neighbors[u]) {
            if (excluded_edges[u][nbr]) {
                continue;
            }
            if (mapped_q[nbr] == -1) {
                return false;
            }
        }
        return true;
    }

    ui terminalLiveAnchorCount(ui u) const
    {
        ui count = 0;
        for (ui anchor : q_neighbors[u]) {
            if (!excluded_edges[u][anchor] && mapped_q[anchor] != -1) {
                count++;
            }
        }
        return count;
    }

    TerminalScan markTerminalVertices(DfsBuffer &buf)
    {
        TerminalScan scan;
        if (buf.terminal_skip.size() != qn) {
            buf.terminal_skip.assign(qn, 0);
        }
        else {
            std::fill(buf.terminal_skip.begin(), buf.terminal_skip.end(), 0);
        }
        buf.terminal_vertices.clear();
        buf.active_terminal_vertices.clear();

        for (ui u = 0; u < qn; ++u) {
            if (mapped_q[u] != -1) {
                continue;
            }

            scan.unmatched_count++;
            if (isTerminalQueryVertex(u)) {
                scan.terminal_count++;
                buf.terminal_vertices.push_back(u);
                if (frontier_pos[u] != -1) {
                    buf.terminal_skip[u] = 1;
                    buf.active_terminal_vertices.push_back(u);
                    scan.terminal_frontier_count++;
                }
            }
            else if (frontier_pos[u] != -1) {
                scan.nonterminal_frontier_count++;
            }
        }
        return scan;
    }

    ui terminalMissingDelta(ui u, ui v, ui limit = std::numeric_limits<ui>::max()) const
    {
        ui delta = 0;
        for (ui anchor : q_neighbors[u]) {
            if (excluded_edges[u][anchor] || mapped_q[anchor] == -1) {
                continue;
            }
            if (!data_graph->hasEdge(v, (ui)mapped_q[anchor])) {
                delta++;
                if (delta > limit) {
                    return delta;
                }
            }
        }
        return delta;
    }

    template <typename Visitor>
    bool visitTerminalSupportedCandidates(ui u, Visitor visit)
    {
        if (++data_vertex_mark_token == 0) {
            std::fill(data_vertex_mark.begin(), data_vertex_mark.end(), 0);
            data_vertex_mark_token = 1;
        }
        ui token = data_vertex_mark_token;

        for (ui anchor : q_neighbors[u]) {
            if (excluded_edges[u][anchor] || mapped_q[anchor] == -1) {
                continue;
            }

            ui deg = 0;
            const ui *nbrs = data_graph->getVertexNeighbors((ui)mapped_q[anchor], deg);
            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (data_vertex_mark[v] == token) {
                    continue;
                }
                data_vertex_mark[v] = token;
#ifndef NDEBUG
                stats.terminal_bucket_candidate_checks++;
#endif

                if (!candidates[u].contains(v)) {
                    continue;
                }
                if (mapped_g[v] != -1) {
                    continue;
                }
                if (excluded_cands[u].contains(v)) {
                    continue;
                }
                if (!visit(v)) {
                    return false;
                }
            }
        }
        return true;
    }

    bool buildTerminalCandidateBuckets(ui u, ui cost, vector<vector<ui>> &buckets,
        ui &feasible_count, ui &min_delta)
    {
        assert(cost <= threshold);
        ui remaining_budget = threshold - cost;
        buckets.assign((size_t)remaining_budget + 1, vector<ui>());
        feasible_count = 0;
        min_delta = std::numeric_limits<ui>::max();
        ui live_anchor_count = terminalLiveAnchorCount(u);
        if (live_anchor_count == 0) {
            return false;
        }

        visitTerminalSupportedCandidates(u, [&](ui v) -> bool {
            ui missing_delta = terminalMissingDelta(u, v, remaining_budget);
            if (missing_delta > remaining_budget) {
                return true;
            }
            if (missing_delta >= live_anchor_count) {
                return true;
            }

            buckets[missing_delta].push_back(v);
            feasible_count++;
            if (missing_delta < min_delta) {
                min_delta = missing_delta;
            }
            return true;
        });

        return feasible_count > 0;
    }

    bool terminalVertexHasFeasibleCandidate(ui u, ui cost)
    {
        assert(cost <= threshold);
        ui remaining_budget = threshold - cost;
        ui live_anchor_count = terminalLiveAnchorCount(u);
        if (live_anchor_count == 0) {
            return false;
        }

        bool found = false;
        visitTerminalSupportedCandidates(u, [&](ui v) -> bool {
            ui missing_delta = terminalMissingDelta(u, v, remaining_budget);
            if (missing_delta <= remaining_budget &&
                missing_delta < live_anchor_count) {
                found = true;
                return false;
            }
            return true;
        });
        return found;
    }

    bool activeTerminalVerticesHaveCandidate(const vector<ui> &terminal_vertices,
        ui cost)
    {
        for (ui u : terminal_vertices) {
            if (!terminalVertexHasFeasibleCandidate(u, cost)) {
                return false;
            }
        }
        return true;
    }

    bool buildTerminalTailVertices(const vector<ui> &terminal_vertices, ui cost,
        vector<TerminalTailVertex> &tail_vertices)
    {
        tail_vertices.clear();
        tail_vertices.reserve(terminal_vertices.size());

        for (ui u : terminal_vertices) {
            TerminalTailVertex tail_vertex;
            tail_vertex.u = u;
            if (!buildTerminalCandidateBuckets(u, cost, tail_vertex.buckets,
                tail_vertex.feasible_count, tail_vertex.min_delta)) {
                return false;
            }
            tail_vertices.push_back(std::move(tail_vertex));
        }

        std::sort(tail_vertices.begin(), tail_vertices.end(),
            [this](const TerminalTailVertex &lhs, const TerminalTailVertex &rhs) {
                if (lhs.feasible_count != rhs.feasible_count) {
                    return lhs.feasible_count < rhs.feasible_count;
                }
                if (lhs.min_delta != rhs.min_delta) {
                    return lhs.min_delta < rhs.min_delta;
                }
                if (q_degree[lhs.u] != q_degree[rhs.u]) {
                    return q_degree[lhs.u] > q_degree[rhs.u];
                }
                return lhs.u < rhs.u;
            });
        return true;
    }

    void recordTerminalPrune()
    {
        stats.prun_calls++;
#ifndef NDEBUG
        stats.terminal_prune_calls++;
#endif
        recordBranchProfilePrune();
    }

    void enumerateTerminalTail(size_t pos, ui cost,
        vector<TerminalTailVertex> &tail_vertices)
    {
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }

#ifndef NDEBUG
        stats.terminal_tail_calls++;
#endif
        assert(cost <= threshold);

        if (pos == tail_vertices.size()) {
            assert(part_M.size() == qn);
            stats.result_count++;
            noteOutputLimitIfReached();
            recordBranchProfileAnswer();
#ifndef NDEBUG
            recordBranchOrderAnswerDebug();
            results_ptr->push_back(part_M);
#endif
            return;
        }

        TerminalTailVertex &tail_vertex = tail_vertices[pos];
        ui u = tail_vertex.u;
        assert(mapped_q[u] == -1);

        ui remaining_budget = threshold - cost;
        ui max_delta = std::min((ui)tail_vertex.buckets.size() - 1, remaining_budget);
        for (ui missing_delta = 0; missing_delta <= max_delta; ++missing_delta) {
            const vector<ui> &bucket = tail_vertex.buckets[missing_delta];
            for (ui v : bucket) {
                if (mapped_g[v] != -1) {
                    continue;
                }

#ifndef NDEBUG
                recordBranchOrderDebug(u);
#endif
                mapped_q[u] = (int)v;
                mapped_g[v] = (int)u;
                part_M.push_back({ u, v });

                enumerateTerminalTail(pos + 1, cost + missing_delta, tail_vertices);

                part_M.pop_back();
                mapped_g[v] = -1;
                mapped_q[u] = -1;

                if (outputLimitReached()) {
                    return;
                }
            }
        }
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
        // printf("part_M.size() = %zu, cost = %u\n", part_M.size(), cost);
        if (outputLimitReached()) {
            stats.output_limit_reached = true;
            return;
        }

        assert(part_M.size() <= qn);
        assert(cost <= threshold);

        recordBranchProfileDfsEnter();
        BranchProfileDfsScope branch_profile_dfs_scope(*this);
        if (part_M.size() == qn) {
            stats.recursion_calls++;
            stats.result_count++;
            noteOutputLimitIfReached();
            recordBranchProfileAnswer();
#ifndef NDEBUG
            recordBranchOrderAnswerDebug();
            results_ptr->push_back(part_M);
#endif
            return;
        }

        stats.recursion_calls++;

        DfsBuffer &buf = dfsBufferForDepth(part_M.size());
        buf.clearLocal();
        SupportUndoScope support_undo_scope(*this, buf.local_support_snapshots);
        vector<ui> &candidate_vertices = buf.candidate_vertices;
        vector<ui> &candidate_anchor_counts = buf.candidate_anchor_counts;
        vector<ActiveEdge> &top_edges = buf.top_edges;
        ui current_cost = cost;
        double selected_component_support_sum = std::numeric_limits<double>::max();
#if !CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        bool selected_covered_component = false;
#endif
        bool has_zero_support_component = false;
        const vector<char> *terminal_skip_vertices = nullptr;

        if (terminal_buckets_enabled) {
            TerminalScan terminal_scan = markTerminalVertices(buf);

            if (terminal_scan.allRemainingTerminal()) {
                if (!buildTerminalTailVertices(buf.terminal_vertices, current_cost,
                    buf.terminal_tail_vertices)) {
                    recordTerminalPrune();
                    return;
                }

                enumerateTerminalTail(0, current_cost, buf.terminal_tail_vertices);
                return;
            }

            if (terminal_scan.terminal_frontier_count > 0) {
                if (!activeTerminalVerticesHaveCandidate(buf.active_terminal_vertices,
                    current_cost)) {
                    recordTerminalPrune();
                    return;
                }

                if (terminal_scan.hasNonterminalFrontier()) {
                    terminal_skip_vertices = &buf.terminal_skip;
#ifndef NDEBUG
                    stats.terminal_delayed_vertices += terminal_scan.terminal_frontier_count;
#endif
                }
            }
        }

        Timer t_frontier;
        ui max_branch_edges = threshold - current_cost + 1;
        {
#ifndef NDEBUG
            Timer t_select;
#endif
            if (!buf.branch_selector.collectTopActiveEdges(active_frontier,
                max_branch_edges, top_edges, &buf.edge_score_cache,
                terminal_skip_vertices)) {
                stats.frontier_time += t_frontier.elapsed();
                return;
            }
#ifndef NDEBUG
            stats.frontier_select_time += t_select.elapsed();
#endif
        }

#if !CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        bool all_selected_edges_have_zero_support = std::all_of(
            top_edges.begin(), top_edges.end(),
            [](const ActiveEdge &edge) {
                return edge.anchor_support == 0;
            });
        if (all_selected_edges_have_zero_support) {
            stats.frontier_time += t_frontier.elapsed();
            stats.prun_calls++;
            recordBranchProfilePrune();
            return;
        }
#endif

        {
#ifndef NDEBUG
            Timer t_component;
#endif
#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY
            (void)buf.branch_selector.restrictTopEdgesToCoveredComponent(top_edges,
                buf.edge_score_cache, buf.component_id, buf.component_frontiers,
                buf.component_checked, buf.component_edges,
                buf.best_component_edges, selected_component_support_sum,
                has_zero_support_component, terminal_skip_vertices);
#else
            selected_covered_component = buf.branch_selector.restrictTopEdgesToCoveredComponent(top_edges,
                buf.edge_score_cache, buf.component_id, buf.component_frontiers,
                buf.component_checked, buf.component_edges,
                buf.best_component_edges, selected_component_support_sum,
                has_zero_support_component, terminal_skip_vertices);
#endif
#ifndef NDEBUG
            stats.frontier_component_time += t_component.elapsed();
#endif
        }
        stats.frontier_time += t_frontier.elapsed();

#if !CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        if (has_zero_support_component ||
            (selected_covered_component && selected_component_support_sum == 0)) {
            stats.prun_calls++;
            recordBranchProfilePrune();
            return;
        }
#endif

#ifndef NDEBUG
        branch_timing_depth++;
#endif
        Timer t_branch;
        long long child_dfs_time = 0;

        ui first_branch_edge = 0;
        bool pruned_by_forced_zero = false;
#if !CDE_EDGE_IE_TOPK_SUPPORT_DECAY
        for (; first_branch_edge < top_edges.size(); ++first_branch_edge) {
            const ActiveEdge &edge = top_edges[first_branch_edge];
            if (edge.anchor_support != 0) break;

            ui u = edge.u;
            ui ua = edge.anchor;
            assert(u < qn && ua < qn);
            assert(mapped_q[ua] != -1 && mapped_q[u] == -1);
            assert(!excluded_edges[u][ua] && frontier_pos[u] != -1);
            assert(q_matrix[u][ua]);

#ifndef NDEBUG
            Timer t_exclude_update;
#endif
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);
            recordBranchProfileForcedZeroExclude(u, ua);
#ifndef NDEBUG
            stats.exclude_update_time += t_exclude_update.elapsed();
#endif

            if (current_cost > threshold) {
                stats.prun_calls++;
                recordBranchProfilePrune();
                pruned_by_forced_zero = true;
                break;
            }
        }
#endif

        for (ui edge_idx = first_branch_edge; !pruned_by_forced_zero && edge_idx < top_edges.size(); ++edge_idx) {
            const ActiveEdge &edge = top_edges[edge_idx];
            if (current_cost > threshold) break;

            ui u = edge.u;
            ui ua = edge.anchor;
            assert(u < qn && ua < qn);
            assert(mapped_q[ua] != -1 && mapped_q[u] == -1);
            assert(!excluded_edges[u][ua] && frontier_pos[u] != -1);
            assert(q_matrix[u][ua]);

            // calculate the candidate vertices of u supported by anchor ua
            candidate_vertices.clear();
#ifndef NDEBUG
            {
                Timer t_cal_edge_support;
                calEdgeSupport(u, ua, [&](ui v) {candidate_vertices.push_back(v);});
                stats.branch_cal_edge_support_time += t_cal_edge_support.elapsed();
            }
#else
            calEdgeSupport(u, ua, [&](ui v) {candidate_vertices.push_back(v);});
#endif
            recordBranchProfileSelectedEdge(u, ua, candidate_vertices.size());
            // number of live anchors of u excluding anchor ua
            ui anchor_num = 0;
#ifndef NDEBUG
            {
                Timer t_count_anchors;
                anchor_num = countAnchors(u, ua, candidate_vertices, candidate_anchor_counts);
                stats.branch_count_anchors_time += t_count_anchors.elapsed();
            }
#else
            anchor_num = countAnchors(u, ua, candidate_vertices, candidate_anchor_counts);
#endif

            // Include branch: use this frontier-anchor edge to add exactly one new query vertex.
#ifndef NDEBUG
            Timer t_candidate_loop;
            long long candidate_child_time_before = child_dfs_time;
            long long candidate_support_update_time = 0;
#endif
            for (ui i = 0; i < candidate_vertices.size(); ++i) {
                ui v = candidate_vertices[i];
                recordBranchProfileIncludeConsidered(u, ua);
                ui delta = anchor_num - candidate_anchor_counts[i];
                ui next_cost = current_cost + delta;
                if (next_cost > threshold) {
                    recordBranchProfileIncludeCostPruned(u, ua);
                    continue;
                }

                recordBranchProfileIncludeTaken(u, ua);

#ifndef NDEBUG
                recordBranchOrderDebug(u);
#endif

                mapped_q[u] = (int)v;
                mapped_g[v] = (int)u;
                part_M.push_back({ u, v });

#ifndef NDEBUG
                long long support_update_before = stats.support_update_time;
#endif
                updateFrontier(u, true);
#ifndef NDEBUG
                candidate_support_update_time += stats.support_update_time - support_update_before;
#endif

                Timer t_child;
                dfs(next_cost);
                child_dfs_time += t_child.elapsed();

                part_M.pop_back();
                mapped_g[v] = -1;
                mapped_q[u] = -1;

#ifndef NDEBUG
                support_update_before = stats.support_update_time;
#endif
                updateFrontier(u, false);
#ifndef NDEBUG
                candidate_support_update_time += stats.support_update_time - support_update_before;
#endif

                if (outputLimitReached()) break;
            }
#ifndef NDEBUG
            {
                long long candidate_elapsed = t_candidate_loop.elapsed();
                long long candidate_excluded_time =
                    (child_dfs_time - candidate_child_time_before) +
                    candidate_support_update_time;
                stats.candidate_loop_time += candidate_elapsed > candidate_excluded_time
                    ? candidate_elapsed - candidate_excluded_time : 0;
            }
#endif

            if (outputLimitReached()) break;

            // Exclude branch: skip this edge, no recursion, just update cost and state.
#ifndef NDEBUG
            Timer t_exclude_update;
#endif
            current_cost++;
            excludeFrontierEdge(u, ua);
            buf.recordExcludedEdge(u, ua);

            if (current_cost > threshold) {
                recordBranchProfileExclude(u, ua, true);
#ifndef NDEBUG
                stats.exclude_update_time += t_exclude_update.elapsed();
#endif
                break;
            }
            recordBranchProfileExclude(u, ua, false);

#ifndef NDEBUG
            recordBranchOrderDebug(u);
#endif

            for (ui v : candidate_vertices) {
                if (!excluded_cands[u].contains(v)) {
                    excluded_cands[u].insert(v);
                    buf.recordExcludedCands(u, v);
                }
            }
#ifndef NDEBUG
            stats.exclude_update_time += t_exclude_update.elapsed();
#endif
        }

        {
            long long branch_elapsed = t_branch.elapsed();
#ifndef NDEBUG
            branch_timing_depth--;
#endif
            long long excluded_time = child_dfs_time;
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
