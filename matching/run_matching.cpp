#include "matching/run_matching.h"

#include "utility/popl.hpp"

#include <array>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <limits>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#include <unistd.h>
#endif

using namespace popl;
using namespace std;

namespace ssm {
namespace detail {
    size_t next_output_report_count = numeric_limits<size_t>::max();
}
}

namespace {

    constexpr array<size_t, 4> output_checkpoints = {
        1000000ULL,
        10000000ULL,
        100000000ULL,
        1000000000ULL,
    };

    bool reported_result_count_set = false;
    size_t reported_result_count = 0;
    ssm::AlgorithmPhaseTimes reported_phase_times;
    ssm::AlgorithmKeyMetrics reported_key_metrics;
    string progress_algorithm;
    chrono::steady_clock::time_point progress_start;
    size_t next_output_checkpoint_index = output_checkpoints.size();
    bool progress_reporting_active = false;

    void begin_output_progress_reporting(const string &algorithm)
    {
        progress_algorithm = algorithm;
        progress_start = chrono::steady_clock::now();
        next_output_checkpoint_index = 0;
        progress_reporting_active = true;
        ssm::detail::next_output_report_count = output_checkpoints.front();
    }

    void end_output_progress_reporting()
    {
        progress_reporting_active = false;
        next_output_checkpoint_index = output_checkpoints.size();
        ssm::detail::next_output_report_count = numeric_limits<size_t>::max();
    }

    const char *format_optional_size(
        bool available, size_t value, char *buffer, size_t buffer_size)
    {
        if (!available) return "NA";
        snprintf(buffer, buffer_size, "%zu", value);
        return buffer;
    }

    const char *format_optional_ull(
        bool available, unsigned long long value, char *buffer, size_t buffer_size)
    {
        if (!available) return "NA";
        snprintf(buffer, buffer_size, "%llu", value);
        return buffer;
    }

    bool parse_control_fd(const char *environment_name, int &fd)
    {
#if defined(__unix__) || defined(__APPLE__)
        const char *fd_text = getenv(environment_name);
        if (fd_text == nullptr || *fd_text == '\0') return false;

        errno = 0;
        char *end = nullptr;
        long parsed_fd = strtol(fd_text, &end, 10);
        if (errno != 0 || end == fd_text || *end != '\0' ||
            parsed_fd < 0 || parsed_fd > INT_MAX) {
            return false;
        }
        fd = static_cast<int>(parsed_fd);
        return true;
#else
        (void)environment_name;
        (void)fd;
        return false;
#endif
    }

    bool notify_runner_ready(const string &algorithm, long long load_time_us,
        int requested_threshold, int effective_threshold)
    {
        printf("SSM_READY algorithm=%s load_ms=%.4lf peak_rss_kb=%zu "
               "requested_threshold=%d effective_threshold=%d\n",
            algorithm.c_str(), load_time_us / 1000.0,
            ssm::peak_resident_memory_kb(), requested_threshold,
            effective_threshold);
        fflush(stdout);

#if defined(__unix__) || defined(__APPLE__)
        int ready_fd = -1;
        int ack_fd = -1;
        const bool has_ready_fd = parse_control_fd("SSM_READY_FD", ready_fd);
        const bool has_ack_fd = parse_control_fd("SSM_ACK_FD", ack_fd);
        if (!has_ready_fd && !has_ack_fd) return true;
        if (!has_ready_fd || !has_ack_fd) return false;

        const char token[] = "READY\n";
        size_t offset = 0;
        bool ready_written = true;
        while (offset < sizeof(token) - 1) {
            ssize_t written = write(ready_fd, token + offset,
                sizeof(token) - 1 - offset);
            if (written > 0) {
                offset += static_cast<size_t>(written);
            }
            else if (written < 0 && errno == EINTR) {
                continue;
            }
            else {
                ready_written = false;
                break;
            }
        }
        close(ready_fd);
        if (!ready_written) {
            close(ack_fd);
            return false;
        }

        char ack = '\0';
        while (true) {
            ssize_t bytes_read = read(ack_fd, &ack, 1);
            if (bytes_read == 1) break;
            if (bytes_read < 0 && errno == EINTR) continue;
            close(ack_fd);
            return false;
        }
        close(ack_fd);
        return ack == 'A';
#else
        return true;
#endif
    }

    LabelID label2int(const string str, map<string, LabelID> &M)
    {
        if (M.find(str) == M.end()) M[str] = M.size();
        return M[str];
    }

    void load_graph(const string &input_graph, Graph *graph,
        map<string, LabelID> &vM, map<string, LabelID> &eM)
    {
        ifstream fin(input_graph);
        string line;

        while (getline(fin, line)) if (!line.empty() && line[0] == 't') break;

        if (!fin) {
            printf("!!! Cannot open graph file %s !!!\n", input_graph.c_str());
            assert(false);
            exit(1);
        }

        istringstream head(line);
        char tchar;
        string sharp, id;
        if (!(head >> tchar >> sharp >> id) || tchar != 't') {
            fprintf(stderr, "!!! Invalid graph header line: %s !!!\n", line.c_str());
            assert(false);
            exit(1);
        }

        vector<pair<ui, LabelID> > vertices;
        vector<pair<pair<ui, ui>, LabelID> > undirected_edges;

        while (getline(fin, line)) {
            if (line.empty()) continue;

            char type = line[0];
            if (type == 't') break;

            istringstream iss(line);

            if (type == 'v') {
                char c;
                ui vid;
                string vlab;
                if (!(iss >> c >> vid >> vlab)) {
                    fprintf(stderr, "!!! Invalid vertex line: %s !!!\n", line.c_str());
                    assert(false);
                    exit(1);
                }
                vertices.emplace_back(vid, label2int(vlab, vM));
            }
            else if (type == 'e') {
                char c;
                ui u, v;
                string elab;
                if (!(iss >> c >> u >> v)) {
                    fprintf(stderr, "!!! Invalid edge line: %s !!!\n", line.c_str());
                    assert(false);
                    exit(1);
                }
                if (!(iss >> elab)) elab = "__NO_EDGE_LABEL__";
                if (u == v) continue;
                LabelID L = label2int(elab, eM);
                undirected_edges.emplace_back(make_pair(min(u, v), max(u, v)), L);
            }
            else {
                fprintf(stderr, "!!! Unknown line type: %s !!!\n", line.c_str());
                assert(false);
                exit(1);
            }
        }

        sort(vertices.begin(), vertices.end());
        sort(undirected_edges.begin(), undirected_edges.end());
        vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());
        undirected_edges.erase(unique(undirected_edges.begin(), undirected_edges.end()), undirected_edges.end());

        vector<pair<pair<ui, ui>, LabelID> > edges;
        edges.reserve(undirected_edges.size() * 2);

        for (auto &e : undirected_edges) {
            edges.emplace_back(make_pair(e.first.first, e.first.second), e.second);
            edges.emplace_back(make_pair(e.first.second, e.first.first), e.second);
        }

        sort(edges.begin(), edges.end());

        graph->build_graph(id, vertices, edges);
    }

#ifndef NDEBUG
    static string mapping_to_key(const vector<pair<ui, ui> > &mapping)
    {
        vector<pair<ui, ui> > tmp = mapping;
        sort(tmp.begin(), tmp.end());

        ostringstream oss;
        for (const auto &p : tmp) {
            oss << p.first << "->" << p.second << ";";
        }
        return oss.str();
    }

    void check_correctness(const Graph *query_graph, const Graph *data_graph,
        const vector<vector<pair<ui, ui> > > &M_ANS, const ui threshold)
    {
        printf("--- Verifying Results ---\n");
        ui qn = query_graph->getVerticesCount();
        ui gn = data_graph->getVerticesCount();

        unordered_set<string> seen;
        vector<vector<pair<ui, ui> > > unique_results;
        size_t duplicate_count = 0;

        unique_results.reserve(M_ANS.size());
        for (const auto &mapping : M_ANS) {
            string key = mapping_to_key(mapping);
            if (seen.insert(key).second) unique_results.push_back(mapping);
            else duplicate_count++;
        }

        printf(">>> Duplicate mappings removed: %zu\n", duplicate_count);
        printf(">>> Unique mappings count   : %zu\n", unique_results.size());

        size_t invalid_count = 0;

        for (size_t i = 0; i < unique_results.size(); ++i) {
            const auto &mapping = unique_results[i];
            bool is_valid = true;

            if (mapping.size() != qn) {
                printf("Result %zu Error: Mapping size %zu != Query size %u\n", i, mapping.size(), qn);
                invalid_count++;
                continue;
            }

            vector<int> map_q(qn, -1);

            for (const auto &p : mapping) {
                ui u = p.first;
                ui v = p.second;

                assert(u < qn && v < gn);
                if (u >= qn || v >= gn) {
                    printf("Result %zu Error: Out of range mapping u%u -> v%u\n", i, u, v);
                    is_valid = false;
                    break;
                }

                map_q[u] = (int)v;

                if (query_graph->getVertexLabel(u) != data_graph->getVertexLabel(v)) {
                    printf("Result %zu Error: Label Mismatch u%u(L%u) -> v%u(L%u)\n",
                        i, u, query_graph->getVertexLabel(u),
                        v, data_graph->getVertexLabel(v));
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid) {
                invalid_count++;
                continue;
            }

            for (ui u = 0; u < qn; ++u) {
                if (map_q[u] == -1) {
                    printf("Result %zu Error: Unassigned query vertex u%u\n", i, u);
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid) {
                invalid_count++;
                continue;
            }

            for (ui u = 0; u < qn; u++) {
                assert(map_q[u] != -1);
            }
            ui start_u = 0;

            vector<bool> visited_q(qn, false);

            queue<ui> bfs_q;
            bfs_q.push(start_u);
            ui vis_count = 0;

            while (!bfs_q.empty()) {
                ui u = bfs_q.front(); bfs_q.pop();
                if (visited_q[u]) continue;
                visited_q[u] = true;
                vis_count++;
                ui deg_u; const ui *nbrs_u = query_graph->getVertexNeighbors(u, deg_u);
                for (ui j = 0; j < deg_u; ++j) {
                    ui u_nbr = nbrs_u[j];
                    if (!visited_q[u_nbr] && data_graph->hasEdge((ui)map_q[u], map_q[u_nbr])) {
                        bfs_q.push(u_nbr);
                    }
                }
            }

            if (vis_count < qn) {
                printf("Result %zu Error: Disconnected preserved-query-edge graph (visited %u / %u)\n",
                    i, vis_count, qn);
                invalid_count++;
                continue;
            }

            ui missing_edges = 0;
            for (ui u = 0; u < qn; ++u) {
                ui count;
                const ui *neighbors = query_graph->getVertexNeighbors(u, count);
                for (ui k = 0; k < count; ++k) {
                    ui u_next = neighbors[k];
                    if (u < u_next) {
                        int v = map_q[u];
                        int v_next = map_q[u_next];

                        if (!data_graph->hasEdge((ui)v, (ui)v_next)) {
                            missing_edges++;
                        }
                    }
                }
            }

            if (missing_edges > threshold) {
                printf("Result %zu Error: Missing Edges %u > Threshold %u\n",
                    i, missing_edges, threshold);
                is_valid = false;
            }

            if (!is_valid) invalid_count++;
        }

        if (invalid_count == 0) {
            printf(">>> All %lu unique results PASSED verification.\n",
                (unsigned long)unique_results.size());
        }
        else {
            printf(">>> Verification FAILED: %zu/%lu unique results are invalid.\n",
                invalid_count, (unsigned long)unique_results.size());
        }
        printf("-------------------------\n");
    }
#endif

    bool is_connected_graph(const Graph *graph)
    {
        ui n = graph->getVerticesCount();
        if (n == 0) return false;
        if (n == 1) return true;

        vector<bool> visited(n, false);
        queue<ui> q;

        q.push(0);
        visited[0] = true;

        ui visited_count = 0;

        while (!q.empty()) {
            ui u = q.front();
            q.pop();
            visited_count++;

            ui deg;
            const ui *nbrs = graph->getVertexNeighbors(u, deg);

            for (ui i = 0; i < deg; ++i) {
                ui v = nbrs[i];
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }

        return visited_count == n;
    }

    double to_ms(long long us)
    {
        return us / 1000.0;
    }

} // namespace

namespace ssm {

    size_t peak_resident_memory_kb()
    {
#if defined(__unix__) || defined(__APPLE__)
        // Process-lifetime high-water RSS.  This intentionally includes graph
        // loading as well as preprocessing and search, matching whole-process
        // peak-memory accounting in compare.sh's parent-side runner.
        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss <= 0) {
            return 0;
        }
#if defined(__APPLE__)
        // macOS reports bytes; Linux and the other supported Unix benchmark
        // environments report KiB.
        return (static_cast<size_t>(usage.ru_maxrss) + 1023U) / 1024U;
#else
        return static_cast<size_t>(usage.ru_maxrss);
#endif
#else
        return 0;
#endif
    }

    namespace detail {

        void emit_output_checkpoint(size_t count, unsigned long long recursion_calls)
        {
            if (!progress_reporting_active) {
                next_output_report_count = numeric_limits<size_t>::max();
                return;
            }

            while (next_output_checkpoint_index < output_checkpoints.size() &&
                   count >= output_checkpoints[next_output_checkpoint_index]) {
                const size_t output_limit = output_checkpoints[next_output_checkpoint_index];
                const long long elapsed_us = chrono::duration_cast<chrono::microseconds>(
                    chrono::steady_clock::now() - progress_start).count();
                const size_t peak_rss_kb = peak_resident_memory_kb();
                char recursion_value[32];
                char filter_candidates_buffer[32];
                char preprocessing_buffer[32];
                char search_buffer[32];
                const char *filter_candidates_value = format_optional_size(
                    reported_key_metrics.filter_candidates_available,
                    reported_key_metrics.filter_candidates,
                    filter_candidates_buffer, sizeof(filter_candidates_buffer));
                const char *preprocessing_value = "NA";
                const char *search_value = "NA";
                if (reported_phase_times.available) {
                    snprintf(preprocessing_buffer, sizeof(preprocessing_buffer),
                        "%.4lf", to_ms(reported_phase_times.preprocessing_us));
                    const long long search_us = std::max(
                        0LL, elapsed_us - reported_phase_times.preprocessing_us);
                    snprintf(search_buffer, sizeof(search_buffer),
                        "%.4lf", to_ms(search_us));
                    preprocessing_value = preprocessing_buffer;
                    search_value = search_buffer;
                }
                if (recursion_calls == numeric_limits<unsigned long long>::max()) {
                    snprintf(recursion_value, sizeof(recursion_value), "NA");
                }
                else {
                    snprintf(recursion_value, sizeof(recursion_value), "%llu", recursion_calls);
                }

                printf("SSM_CHECKPOINT algorithm=%s output_limit=%zu count=%zu "
                       "run_ms=%.4lf preprocessing_ms=%s search_ms=%s "
                       "peak_rss_kb=%zu filter_candidates=%s recursion_calls=%s\n",
                    progress_algorithm.c_str(), output_limit, count,
                    to_ms(elapsed_us), preprocessing_value, search_value,
                    peak_rss_kb, filter_candidates_value, recursion_value);
                // compare.sh may terminate the process later.  Make every
                // completed checkpoint visible in the redirected output now.
                fflush(stdout);
                next_output_checkpoint_index++;
            }

            next_output_report_count =
                next_output_checkpoint_index < output_checkpoints.size()
                ? output_checkpoints[next_output_checkpoint_index]
                : numeric_limits<size_t>::max();
        }

    } // namespace detail

    void clear_reported_result_count()
    {
        reported_result_count_set = false;
        reported_result_count = 0;
    }

    void set_reported_result_count(size_t count)
    {
        reported_result_count_set = true;
        reported_result_count = count;
    }

    size_t get_reported_result_count(const MatchResults &results)
    {
        return reported_result_count_set ? reported_result_count : results.size();
    }

    void clear_reported_phase_times()
    {
        reported_phase_times = AlgorithmPhaseTimes();
    }

    void report_preprocessing_complete(long long preprocessing_us)
    {
        set_reported_phase_times(preprocessing_us, 0);
        char filter_candidates_buffer[32];
        const char *filter_candidates_value = format_optional_size(
            reported_key_metrics.filter_candidates_available,
            reported_key_metrics.filter_candidates,
            filter_candidates_buffer, sizeof(filter_candidates_buffer));
        printf("SSM_PREPROCESSED algorithm=%s preprocessing_ms=%.4lf "
               "filter_candidates=%s peak_rss_kb=%zu\n",
            progress_algorithm.c_str(), to_ms(reported_phase_times.preprocessing_us),
            filter_candidates_value, peak_resident_memory_kb());
        fflush(stdout);
    }

    void set_reported_phase_times(long long preprocessing_us, long long search_us)
    {
        reported_phase_times.preprocessing_us = std::max(0LL, preprocessing_us);
        reported_phase_times.search_us = std::max(0LL, search_us);
        reported_phase_times.available = true;
    }

    AlgorithmPhaseTimes get_reported_phase_times()
    {
        return reported_phase_times;
    }

    void clear_reported_key_metrics()
    {
        reported_key_metrics = AlgorithmKeyMetrics();
    }

    void set_reported_filter_candidates(size_t count)
    {
        reported_key_metrics.filter_candidates = count;
        reported_key_metrics.filter_candidates_available = true;
    }

    void set_reported_recursion_calls(unsigned long long count)
    {
        reported_key_metrics.recursion_calls = count;
        reported_key_metrics.recursion_calls_available = true;
    }

    AlgorithmKeyMetrics get_reported_key_metrics()
    {
        return reported_key_metrics;
    }

    int run_algorithm_main(int argc, char *argv[], const AlgorithmDefinition &algorithm)
    {
        // Keep newline-terminated diagnostics visible even when stdout is
        // redirected by compare.sh. Structured milestones still call fflush
        // explicitly because the runner may terminate the process later.
        setvbuf(stdout, nullptr, _IOLBF, 0);
#ifndef NDEBUG
        printf("**** AFEC suite [%s] verbose diagnostics build at %s %s ***\n",
            algorithm.display_name.c_str(), __TIME__, __DATE__);
#endif

        int threshold = 0;

        OptionParser op("Allowed options");
        auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
        auto data_option = op.add<Value<string> >("d", "data graph", "\'data graph file name\'");
        auto query_option = op.add<Value<string> >("q", "query graph", "\'query graph file name\'");
        auto threshold_option = op.add<Value<int> >("t", "threshold", "\'threshold", 0, &threshold);
        op.parse(argc, argv);

        if (help_option->is_set() || argc == 1) {
            cout << op << endl;
            return 0;
        }
        if (!data_option->is_set() || !query_option->is_set()) {
            printf("!!! Data graph file name or query graph file name is not provided! Exit !!!\n");
            return 1;
        }
        if (!threshold_option->is_set()) {
            printf("!!! Threshold is not provided! Exit !!!\n");
            return 1;
        }
        if (threshold < 0) {
            printf("!!! Threshold must be non-negative! Exit !!!\n");
            return 1;
        }

        string data_graph_file = data_option->value();
        string query_graph_file = query_option->value();
        const int requested_threshold = threshold;

        Timer t;
        long long load_time = 0;
        long long run_time = 0;

        Graph *query_graph = new Graph();
        Graph *data_graph = new Graph();
        map<string, int> vM, eM;

        t.restart();
        load_graph(query_graph_file, query_graph, vM, eM);
        if (!is_connected_graph(query_graph)) {
            printf("!!! Query graph is disconnected! Exit !!!\n");
            delete query_graph;
            delete data_graph;
            return 1;
        }
        load_graph(data_graph_file, data_graph, vM, eM);

        assert(query_graph->getEdgesCount() % 2 == 0);
        ui query_undirected_edges = query_graph->getEdgesCount() / 2;
        assert(query_undirected_edges >= query_graph->getVerticesCount() - 1);
        threshold = min(threshold, (int)(query_undirected_edges - query_graph->getVerticesCount() + 1));
        load_time = t.elapsed();
        if (!notify_runner_ready(algorithm.key, load_time,
                requested_threshold, threshold)) {
            fprintf(stderr, "SSM handshake failed before algorithm start.\n");
            delete query_graph;
            delete data_graph;
            return 125;
        }

        MatchResults results;

        t.restart();
        clear_reported_result_count();
        clear_reported_phase_times();
        clear_reported_key_metrics();
        begin_output_progress_reporting(algorithm.key);
        algorithm.entry(query_graph, data_graph, results, threshold);
        run_time = t.elapsed();
        end_output_progress_reporting();
        size_t result_count = get_reported_result_count(results);
        size_t peak_rss_kb = peak_resident_memory_kb();
        AlgorithmPhaseTimes phase_times = get_reported_phase_times();
        if (!phase_times.available) {
            phase_times.search_us = run_time;
        }
        AlgorithmKeyMetrics key_metrics = get_reported_key_metrics();
        char filter_candidates_buffer[32];
        char recursion_calls_buffer[32];
        const char *filter_candidates_value = format_optional_size(
            key_metrics.filter_candidates_available, key_metrics.filter_candidates,
            filter_candidates_buffer, sizeof(filter_candidates_buffer));
        const char *recursion_calls_value = format_optional_ull(
            key_metrics.recursion_calls_available, key_metrics.recursion_calls,
            recursion_calls_buffer, sizeof(recursion_calls_buffer));

        printf("%s Results:\n", algorithm.display_name.c_str());
        printf("Algorithm:           %s\n", algorithm.key.c_str());
        printf("Total Time:          %.4lf ms\n", to_ms(run_time));
        printf("Init Time:           %.4lf ms\n", to_ms(phase_times.preprocessing_us));
        printf("  - Filter Candidates: %s\n", filter_candidates_value);
        printf("Search Time:         %.4lf ms\n", to_ms(phase_times.search_us));
        printf("Recursion Calls:     %s\n", recursion_calls_value);
        printf("Results Found:       %zu\n", result_count);
        printf("Output Limit:        %zu%s\n",
            static_cast<size_t>(MATCH_OUTPUT_LIMIT),
            result_count >= static_cast<size_t>(MATCH_OUTPUT_LIMIT)
                ? " (reached)" : " (not reached)");
        printf("Load Graphs Time:    %.4lf ms\n", to_ms(load_time));
        printf("Peak Process RSS:    %zu KiB\n", peak_rss_kb);
        printf("SSM_SUMMARY algorithm=%s count=%zu load_ms=%.4lf run_ms=%.4lf "
               "total_ms=%.4lf preprocessing_ms=%.4lf search_ms=%.4lf "
               "peak_rss_kb=%zu filter_candidates=%s recursion_calls=%s "
               "requested_threshold=%d effective_threshold=%d "
               "output_limit=%zu output_limit_reached=%d\n",
            algorithm.key.c_str(), result_count, to_ms(load_time), to_ms(run_time),
            to_ms(load_time + run_time), to_ms(phase_times.preprocessing_us),
            to_ms(phase_times.search_us), peak_rss_kb, filter_candidates_value,
            recursion_calls_value, requested_threshold, threshold,
            static_cast<size_t>(MATCH_OUTPUT_LIMIT),
            result_count >= static_cast<size_t>(MATCH_OUTPUT_LIMIT) ? 1 : 0);
        printf("SSM_PHASES algorithm=%s preprocessing_ms=%.4lf search_ms=%.4lf "
               "run_ms=%.4lf end_to_end_ms=%.4lf\n",
            algorithm.key.c_str(), to_ms(phase_times.preprocessing_us),
            to_ms(phase_times.search_us), to_ms(run_time),
            to_ms(load_time + run_time));
        fflush(stdout);

#ifndef NDEBUG
        if (!results.empty()) {
            check_correctness(query_graph, data_graph, results, threshold);
        }
#endif

        delete query_graph;
        delete data_graph;

        return 0;
    }

} // namespace ssm
