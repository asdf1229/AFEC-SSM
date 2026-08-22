#include "matching/algorithms/afec/afec.h"

#ifndef AFEC_ALGORITHM_KEY
#define AFEC_ALGORITHM_KEY "afec"
#endif

#ifndef AFEC_DISPLAY_NAME
#define AFEC_DISPLAY_NAME "AFEC"
#endif

namespace afec {

MatchingSolver::MatchingSolver() : query_graph(nullptr), data_graph(nullptr), results_ptr(nullptr) {}

ui MatchingSolver::chooseRoot()
{
    ui root = 0;
    for (ui u = 1; u < qn; ++u) {
        size_t cand_u = candidates[u].size();
        size_t cand_root = candidates[root].size();
        ui deg_u = q_degree[u];
        ui deg_root = q_degree[root];

        if (cand_u * deg_root < cand_root * deg_u) root = u;
    }
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

#if AFEC_ENABLE_ANCHOR_FRONTIER
    buildAdjIndex();
#if AFEC_ANCHOR_ORDER_FIXED
    initFixedEdgePriorities();
#elif AFEC_ANCHOR_ORDER_DYNAMIC
    initStaticEdgeSupports();
#endif
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
#if AFEC_ENABLE_ANCHOR_FRONTIER
    SearchState state;
    initState(state);
    size_t root_mark = mark();

    for (ui v : candidates[root]) {
        ui next_cost = 0;
        if (tryMapBlackWithDelta(state, 0, root, v, 0, next_cost)) search(state, next_cost);
        rollback(state, root_mark);
        if (outputLimitReached()) break;
    }
#else
    buildNoAFOrder(root);
    NoAnchorFrontierState state;
    state.mapped_q.assign(qn, -1);
    state.used_data_flag.assign(gn, 0);
    state.part_M.reserve(qn);
    searchNoAF(state, 0, 0);
#endif

    stats.search_time = t_search.elapsed();
    stats.total_time = stats.init_time + stats.search_time;
}

void MatchingSolver::printStats() const
{
    auto pct = [](long long part, long long whole) -> double {
        return whole > 0 ? (double)part / whole * 100.0 : 0.0;
        };

    printf("\n--- AFEC Time Analysis ---\n");
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
#if AFEC_ENABLE_SPOKE_FILTERING
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
#if AFEC_ENABLE_COST_SPLIT
    printf("Cost-Split Calls:    %lld\n", stats.split_calls);
    printf("Cost-Split Branches: %lld\n", stats.split_branches);
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
    ssm::report_result_progress(stats.result_count, stats.recursion_calls);
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

#if !AFEC_ENABLE_ANCHOR_FRONTIER
    no_af_order.clear();
    no_af_visited.clear();
    no_af_queue.clear();
#endif

    resetBuffers();
    undo_stack.clear();
    candidate_adj_index.clear();
    candidate_adj_pool.clear();
    stats = TimeStats();
#if AFEC_ANCHOR_ORDER_DYNAMIC
    static_candidate_count.clear();
    static_edge_support.clear();
#endif
#if AFEC_ANCHOR_ORDER_FIXED
    static_edge_priority.clear();
#endif
#if AFEC_ANCHOR_ORDER_RANDOM
    // A fixed seed keeps the random-order ablation reproducible.
    random_order_rng.seed(0xC0DEu);
#endif
}

} // namespace afec

void Approximate_AFEC(const Graph *query_graph, const Graph *data_graph,
    std::vector<std::vector<std::pair<ui, ui>>> &M_ANS, ui threshold)
{
    Timer t_total;
    t_total.restart();

    afec::MatchingSolver solver;
    const bool initialized = solver.init(query_graph, data_graph, threshold);
    if (initialized) {
        ssm::set_reported_filter_candidates(solver.stats.filter_candidate_count);
    }
    ssm::report_preprocessing_complete(solver.stats.init_time);
    if (initialized) {
        solver.match(M_ANS);
    }

    solver.stats.total_time = t_total.elapsed();
#ifndef NDEBUG
    solver.printStats();
#endif
    ssm::set_reported_phase_times(solver.stats.init_time,
        solver.stats.search_time);
    ssm::set_reported_recursion_calls(solver.stats.recursion_calls);
    ssm::set_reported_result_count(solver.stats.result_count);
}

namespace ssm {

    namespace {

        void run_afec(const Graph *query_graph, const Graph *data_graph, MatchResults &results, ui threshold)
        {
            Approximate_AFEC(query_graph, data_graph, results, threshold);
        }

    } // namespace

    const AlgorithmDefinition &create_algorithm_definition()
    {
        static const AlgorithmDefinition definition = {
            AFEC_ALGORITHM_KEY,
            AFEC_DISPLAY_NAME,
            &run_afec
        };
        return definition;
    }

} // namespace ssm
