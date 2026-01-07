#include "utility/utility.h"
#include "graph/graph.h"
#include "utility/popl.hpp"
#include "matching/algo.h"

using namespace popl;
using namespace std;

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
    long long generate_order_time = 0;
    long long matching_time = 0;

    // query graph q, data graph G, 
    Graph *query_graph = new Graph();
    Graph *data_graph = new Graph();

    map<string, int> vM, eM;

    t.restart();
    load_graph(data_graph_file, data_graph, vM, eM);
    load_graph(query_graph_file, query_graph, vM, eM);
    load_time = t.elapsed();

    // print graphs
    // query_graph->print_graph();
    // data_graph->print_graph();

    // print label counts for data graph
    // map<int, int> label_counts;
    // for(ui i = 0; i < data_graph->getVerticesCount(); ++i) {
    //     label_counts[data_graph->getVertexLabel(i)]++;
    // }
    // printf("Data graph label counts:\n");
    // for(auto const& [label, count] : label_counts) {
    //     printf("Label %d: %d vertices\n", label, count);
    // }

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
    Approximate_Matching_v2(query_graph, data_graph, candidates, M_ANS, threshold);
    matching_time = t.elapsed();

    // ui count = 0;
    printf("count: %lu\n", M_ANS.size());

    printf("Load Graphs Time: %.4lf ms\n", load_time / 1000.0);
    printf("Filter Candidates Time: %.4lf ms\n", filter_candidates_time / 1000.0);
    printf("Generate Order Time: %.4lf ms\n", generate_order_time / 1000.0);
    printf("Matching Time: %.4lf ms\n", matching_time / 1000.0);
    printf("Total Time: %.4lf ms\n", (load_time + filter_candidates_time + generate_order_time + matching_time) / 1000.0);
    // for(auto &m: M_ANS) {
    //     printf("Mapping %d:\n", ++count);
    //     for(auto v: m) {
    //         printf("%d %d\n", v.first, v.second);
    //     }
    //     printf("\n");
    // }

    return 0;
}