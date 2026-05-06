#include "matching/run_matching.h"

#include "utility/popl.hpp"

using namespace popl;
using namespace std;

namespace {

    bool reported_result_count_set = false;
    size_t reported_result_count = 0;

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
                printf("Result %zu Error: Disconnected mapping in data graph image (visited %u / %u)\n",
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

    double to_ms(long long us)
    {
        return us / 1000.0;
    }

} // namespace

namespace ssm_ged {

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

    int run_algorithm_main(int argc, char *argv[], const AlgorithmDefinition &algorithm)
    {
#ifndef NDEBUG
        printf("**** SSM-GED [%s] (Debug) build at %s %s ***\n",
            algorithm.display_name.c_str(), __TIME__, __DATE__);
#else
        printf("**** SSM-GED [%s] (Release) build at %s %s ***\n",
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

        string data_graph_file = data_option->value();
        string query_graph_file = query_option->value();

        Timer t;
        long long load_time = 0;
        long long run_time = 0;

        Graph *query_graph = new Graph();
        Graph *data_graph = new Graph();
        map<string, int> vM, eM;

        t.restart();
        load_graph(query_graph_file, query_graph, vM, eM);
        load_graph(data_graph_file, data_graph, vM, eM);

        assert(query_graph->getEdgesCount() % 2 == 0);
        ui query_undirected_edges = query_graph->getEdgesCount() / 2;
        assert(query_undirected_edges >= query_graph->getVerticesCount() - 1);
        threshold = min(threshold, (int)(query_undirected_edges - query_graph->getVerticesCount() + 1));
        load_time = t.elapsed();

        MatchResults results;

        t.restart();
        clear_reported_result_count();
        algorithm.entry(query_graph, data_graph, results, threshold);
        run_time = t.elapsed();
        size_t result_count = get_reported_result_count(results);

        printf("%s Results:\n", algorithm.display_name.c_str());
        printf("  Result count: %zu\n", result_count);
        printf("  Load Graphs Time: %.4lf ms\n", to_ms(load_time));
        printf("  Run Time: %.4lf ms\n", to_ms(run_time));
        printf("  Total Time: %.4lf ms\n", to_ms(load_time + run_time));
        printf("SSM_GED_SUMMARY algorithm=%s count=%zu load_ms=%.4lf run_ms=%.4lf total_ms=%.4lf\n",
            algorithm.key.c_str(), result_count, to_ms(load_time), to_ms(run_time), to_ms(load_time + run_time));

#ifndef NDEBUG
        if (!results.empty()) {
            check_correctness(query_graph, data_graph, results, threshold);
        }
#endif

        delete query_graph;
        delete data_graph;

        return 0;
    }

} // namespace ssm_ged
