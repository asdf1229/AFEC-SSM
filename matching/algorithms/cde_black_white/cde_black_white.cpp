#include "matching/algorithms/cde_black_white/cde_black_white.h"

#ifndef CDE_BLACK_WHITE_ALGORITHM_KEY
#define CDE_BLACK_WHITE_ALGORITHM_KEY "cde_black_white"
#endif

#ifndef CDE_BLACK_WHITE_DISPLAY_NAME
#define CDE_BLACK_WHITE_DISPLAY_NAME "CDE-Black-White"
#endif

namespace cde_black_white {

MatchingSolver::MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

ui MatchingSolver::chooseRoot()
{
    ui root = 0;
#if CDE_BLACK_WHITE_FIXED_ORDER
    // Match TreeSpan's root rule: minimum filtered candidate count, then
    // maximum query degree, then the smallest vertex id.
    for (ui u = 1; u < qn; ++u) {
        size_t cand_u = candidates[u].size();
        size_t cand_root = candidates[root].size();
        if (cand_u < cand_root ||
            (cand_u == cand_root && q_degree[u] > q_degree[root])) {
            root = u;
        }
    }
#else
    for (ui u = 1; u < qn; ++u) {
        size_t cand_u = candidates[u].size();
        size_t cand_root = candidates[root].size();

        ui deg_u = q_degree[u];
        ui deg_root = q_degree[root];

        if (cand_u * deg_root < cand_root * deg_u) root = u;
    }
#endif
    return root;
}

void MatchingSolver::tarjan(ui u, ui parent, vector<int> &dfn, vector<int> &low, int &tim)
{
    dfn[u] = low[u] = ++tim;
    for (ui v : q_neighbors[u]) {
        if (dfn[v] == 0) {
            tarjan(v, u, dfn, low, tim);
            low[u] = std::min(low[u], low[v]);
            if (low[v] > dfn[u]) {
                q_bridge_matrix[u][v] = 1;
                q_bridge_matrix[v][u] = 1;
            }
        }
        else if (v != parent) {
            low[u] = std::min(low[u], dfn[v]);
        }
    }
}

void MatchingSolver::initQueryBridge()
{
    q_bridge_matrix.assign(qn, vector<char>(qn, 0));
    vector<int> dfn(qn, 0);
    vector<int> low(qn, 0);
    int tim = 0;
    for (ui u = 0; u < qn; ++u) {
        if (dfn[u] == 0) {
            tarjan(u, qn, dfn, low, tim);
        }
    }
}

bool MatchingSolver::isQueryBridgeEdge(ui u, ui v) const
{
    assert(u < qn && v < qn);
    return q_bridge_matrix[u][v] != 0;
}

bool MatchingSolver::init(const Graph *q, const Graph *g, ui match_threshold)
{
    Timer t_init;
    t_init.restart();

    query_graph = q;
    data_graph = g;
    threshold = match_threshold;
    qn = query_graph->getVerticesCount();
    gn = data_graph->getVerticesCount();
    if (qn == 0 || gn == 0) return false;
    max_g_deg = data_graph->getMaxDegree();
    label_count = max(query_graph->getLabelsCount(), data_graph->getLabelsCount());
    resetState();

    // init q_matrix
    q_neighbors.assign(qn, vector<ui>());
    q_degree.assign(qn, 0);
    for (ui u = 0; u < qn; ++u) {
        ui deg = 0;
        const ui *nbrs = query_graph->getVertexNeighbors(u, deg);
        q_neighbors[u].reserve(deg);
        for (ui i = 0; i < deg; ++i) q_neighbors[u].push_back(nbrs[i]);
        q_degree[u] = (ui)q_neighbors[u].size();
    }
    initQueryBridge();

    Timer t_filter;
    bool res = runCandidateFiltering();
    stats.filter_time = t_filter.elapsed();
    if (!res) {
        stats.init_time = t_init.elapsed();
        return false;
    }

    buildAdjIndex();

#if CDE_BLACK_WHITE_FIXED_ORDER
    initFixedEdgePriorities();
#elif !CDE_BLACK_WHITE_RANDOM_ORDER
    initStaticEdgeSupports();
#endif

    stats.init_time = t_init.elapsed();
    return true;
}

void MatchingSolver::match(vector<vector<pair<ui, ui>>> &results)
{
    Timer t_search;
    t_search.restart();

    results_ptr = &results;
    results_ptr->clear();

    ui root = chooseRoot();
    SearchState state;
    initState(state);
    size_t root_mark = mark();

    for (ui v : candidates[root]) {
        ui next_cost = 0;
        if (tryMapBlackWithDelta(state, 0, root, v, 0, next_cost)) search(state, next_cost);
        rollback(state, root_mark);
        if (outputLimitReached()) break;
    }

    stats.search_time = t_search.elapsed();
    stats.total_time = stats.init_time + stats.search_time;
}

void MatchingSolver::printStats() const
{
    auto pct = [](long long part, long long whole) -> double {
        return whole > 0 ? (double)part / whole * 100.0 : 0.0;
        };

    printf("\n--- CDE-Black-White Time Analysis ---\n");
#ifdef NDEBUG
    printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
    printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
    printf("  - Filter Time:     %.4lf ms\n", stats.filter_time / 1000.0);
    printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
    printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
#else
    printf("Total Time:          %.4lf ms\n", stats.total_time / 1000.0);
    printf("Init Time:           %.4lf ms (%.2f%%)\n", stats.init_time / 1000.0, pct(stats.init_time, stats.total_time));
    printf("  - Filter Time:     %.4lf ms (%.2f%% of Init)\n", stats.filter_time / 1000.0, pct(stats.filter_time, stats.init_time));
    printf("    - NLF:           %.4lf ms (%.2f%% of Filter)\n", stats.filter_nlf_time / 1000.0, pct(stats.filter_nlf_time, stats.filter_time));
    printf("    - Bridge:        %.4lf ms (%.2f%% of Filter)\n", stats.filter_bridge_time / 1000.0, pct(stats.filter_bridge_time, stats.filter_time));
#if CDE_BLACK_WHITE_ENABLE_SPOKE_FILTERING
    printf("    - Spoke:         %.4lf ms (%.2f%% of Filter)\n", stats.filter_spoke_time / 1000.0, pct(stats.filter_spoke_time, stats.filter_time));
#endif
    printf("  - Filter Candidates:%u\n", stats.filter_candidate_count);
    printf("Search Time:         %.4lf ms (%.2f%%)\n", stats.search_time / 1000.0, pct(stats.search_time, stats.total_time));
#endif
    printf("Recursion Calls:     %lld\n", stats.recursion_calls);
    printf("Pruning Calls:       %lld\n", stats.prun_calls);
    printf("Terminal Tail Calls: %lld\n", stats.terminal_tail_calls);
    printf("White Rebuilds:      %lld\n", stats.white_bucket_rebuilds);
    printf("Graph hasEdge Calls: %lld\n", stats.graph_has_edge_checks);
    printf("Range Hits/Misses:   %lld / %lld\n",
        stats.candidate_range_hits, stats.candidate_range_misses);
    printf("Range Intersections: %lld\n", stats.candidate_intersection_calls);
    printf("Range Edge Checks:   %lld\n", stats.candidate_edge_check_calls);
#if CDE_BLACK_WHITE_ENABLE_SPLIT
    printf("Split Calls:         %lld\n", stats.split_calls);
    printf("Split Branches:      %lld\n", stats.split_branches);
#endif
    printf("Results Found:       %zu\n", stats.result_count);
#if MATCH_OUTPUT_LIMIT > 0
    printf("Output Limit:        %zu%s\n",
        (size_t)MATCH_OUTPUT_LIMIT,
        stats.output_limit_reached ? " (reached)" : "");
#endif
    printf("-----------------------------------------------------------\n");
    fflush(stdout);
}

bool MatchingSolver::outputLimitReached() const
{
    return (size_t)MATCH_OUTPUT_LIMIT > 0 && stats.result_count >= (size_t)MATCH_OUTPUT_LIMIT;
}

void MatchingSolver::emitResult(const SearchState &state)
{
    assert(state.part_M.size() == qn);
    stats.result_count++;
    if (outputLimitReached()) stats.output_limit_reached = true;
#ifndef NDEBUG
    results_ptr->push_back(state.part_M);
#endif
}

void MatchingSolver::resetState()
{
    candidates.clear();
    candidates.assign(qn, MyBitset(gn));

    q_neighbors.clear();
    q_degree.clear();
    q_bridge_matrix.clear();

    resetBuffers();
    undo_stack.clear();
    candidate_adj_index.clear();
    candidate_adj_pool.clear();
    stats = TimeStats();
    static_candidate_count.clear();
    static_edge_support.clear();
#if CDE_BLACK_WHITE_FIXED_ORDER
    static_edge_priority.clear();
#endif
#if CDE_BLACK_WHITE_RANDOM_ORDER
    // A fixed seed keeps the random-order ablation reproducible.
    random_order_rng.seed(0xC0DEu);
#endif
}

} // namespace cde_black_white

void Approximate_CDE_BlackWhite(const Graph *query_graph, const Graph *data_graph,
    std::vector<std::vector<std::pair<ui, ui>>> &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    cde_black_white::MatchingSolver solver;
    if (solver.init(query_graph, data_graph, threshold)) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
    solver.printStats();
    ssm::set_reported_result_count(solver.stats.result_count);
}

namespace ssm {

    namespace {

        void run_cde_black_white(const Graph *query_graph, const Graph *data_graph, MatchResults &results, ui threshold)
        {
            Approximate_CDE_BlackWhite(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            CDE_BLACK_WHITE_ALGORITHM_KEY,
            CDE_BLACK_WHITE_DISPLAY_NAME,
            &run_cde_black_white
        };
        return definition;
    }

} // namespace ssm
