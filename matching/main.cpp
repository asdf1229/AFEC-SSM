#include "utility/utility.h"
#include "graph/graph.h"
#include "utility/popl.hpp"
#include "matching/algo.h"

using namespace popl;
using namespace std;

LabelID label2int(const string str, map<string, LabelID> &M)
{
	if(M.find(str) == M.end()) M[str] = M.size();
	return M[str];
}

void load_graph(const string &input_graph, Graph *graph,
                map<string, LabelID> &vM, map<string, LabelID> &eM)
{
    ifstream fin(input_graph);
    string line;

    while (getline(fin, line)) if (!line.empty() && line[0] == 't') break;

    if(!fin) {
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

        if(type == 'v') {
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

    for(auto &e: undirected_edges) {
        edges.emplace_back(make_pair(e.first.first, e.first.second), e.second);
        edges.emplace_back(make_pair(e.first.second, e.first.first), e.second);
    }

    sort(edges.begin(), edges.end());

    graph->build_graph(id, vertices, edges);
}

bool calVerticesFilter(const Graph *query_graph, const Graph *data_graph, vector<vector<ui> > &candidates)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();
    candidates.clear();
    candidates.resize(qn);

    for (ui i = 0; i < qn; ++i) {
        LabelID label_i = query_graph->getVertexLabel(i);
        for (ui j = 0; j < gn; ++j) {
            LabelID label_j = data_graph->getVertexLabel(j);
            if (label_i == label_j) candidates[i].push_back(j);
        }
        if (candidates[i].empty())  return false;
    }
    return true;
}

// Usage: [0]exe [1]input_graph [2]k
int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("**** SSM-GED (Debug) build at %s %s ***\n", __TIME__, __DATE__);
#else
	printf("**** SSM-GED (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

    int threshold = 0;

    OptionParser op("Allowed options");
    auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
    auto data_option = op.add<Value<string>>("d", "data graph", "\'data graph file name\'");
	auto query_option = op.add<Value<string>>("q", "query graph", "\'query graph file name\'");
    auto threshold_option = op.add<Value<int>>("t", "threshold", "\'threshold", 0, &threshold);
    op.parse(argc, argv);

    if(help_option->is_set() || argc == 1) cout << op << endl;
	if(!data_option->is_set()||!query_option->is_set()) {
		printf("!!! Data graph file name or query graph file name is not provided! Exit !!!\n");
		return 0;
	}
    if(!threshold_option->is_set()) {
        printf("!!! Threshold is not provided! Exit !!!\n");
        return 0;
    }

    string data_graph_file = data_option->value();
	string query_graph_file = query_option->value();

    Timer t;
    long long load_time = 0;
    long long filter_candidates_time = 0;
    long long matching_time = 0;

    // query graph q, data graph G, 
    Graph *query_graph = new Graph();
    Graph *data_graph = new Graph();
    map<string, int> vM, eM;

    // load graphs
    t.restart();
    load_graph(data_graph_file, data_graph, vM, eM);
    load_graph(query_graph_file, query_graph, vM, eM);
    load_time = t.elapsed();

    // candidate filtering
    vector<vector<ui> > candidates;
    t.restart();
    bool res = calVerticesFilter(query_graph, data_graph, candidates);
    filter_candidates_time = t.elapsed();

    // print candidates size
    // printf("Number of query vertices: %u\n", query_graph->getVerticesCount());
    // for(ui i = 0; i < candidates.size(); i++) {
    //     printf("Number of candidates for query vertex %d, label = %d: %lu\n", i, query_graph->getVertexLabel(i), candidates[i].size());
    // }

    if(res == false) {
        printf("!!! No candidates for some query vertex. No possible mapping. Exit !!!\n");
        assert(false);
        return 1;
    }

    vector<vector<pair<ui, ui> > > M_ANS;
    M_ANS.clear();

    t.restart();
    Approximate_Matching(query_graph, data_graph, candidates, M_ANS, threshold);
    matching_time = t.elapsed();

    printf("count: %lu\n", M_ANS.size());
    printf("Load Graphs Time: %.4lf ms\n", load_time / 1000.0);
    printf("Filter Candidates Time: %.4lf ms\n", filter_candidates_time / 1000.0);
    printf("Matching Time: %.4lf ms\n", matching_time / 1000.0);
    printf("Total Time: %.4lf ms\n", (load_time + filter_candidates_time + matching_time) / 1000.0);

    return 0;
}