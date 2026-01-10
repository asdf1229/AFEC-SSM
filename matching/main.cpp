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

// 将一个 mapping 规范化为 key：按 u 排序后序列化（顺序无关判重）
static std::string mapping_to_key(const std::vector<std::pair<ui, ui>> &mapping)
{
    std::vector<std::pair<ui, ui>> tmp = mapping;
    std::sort(tmp.begin(), tmp.end()); // pair 默认按 first 再 second

    std::ostringstream oss;
    for (const auto &p : tmp) {
        oss << p.first << "->" << p.second << ";";
    }
    return oss.str();
}

void check_correctness(const Graph *query_graph, const Graph *data_graph,
                       const std::vector<std::vector<std::pair<ui, ui>>> &M_ANS, ui threshold)
{
    printf("--- Verifying Results ---\n");
    ui qn = query_graph->getVerticesCount();

    // ============================================================
    // 0) 去重：映射对应完全相同（u->v 相同，pair 顺序无关）视为重复
    // ============================================================
    std::unordered_set<std::string> seen;
    std::vector<std::vector<std::pair<ui, ui>>> unique_results;
    size_t duplicate_count = 0;

    unique_results.reserve(M_ANS.size());
    for (const auto &mapping : M_ANS) {
        std::string key = mapping_to_key(mapping);
        if (seen.insert(key).second) {
            unique_results.push_back(mapping); // 第一次出现，保留
        } else {
            duplicate_count++;                 // 重复，丢弃
        }
    }

    printf(">>> Duplicate mappings removed: %zu\n", duplicate_count);
    printf(">>> Unique mappings count   : %zu\n", unique_results.size());

    // ============================================================
    // 1) 对去重后的结果做正确性验证
    // ============================================================
    size_t invalid_count = 0;

    for (size_t i = 0; i < unique_results.size(); ++i) {
        const auto &mapping = unique_results[i];
        bool is_valid = true;

        // 1. 映射大小检查
        if (mapping.size() != qn) {
            printf("Result %zu Error: Mapping size %zu != Query size %u\n",
                   i, mapping.size(), qn);
            invalid_count++;
            continue;
        }

        // 构建快速查找表: map_q[u] -> v
        std::vector<int> map_q(qn, -1);

        // 2. 检查节点标签 (Label Check)
        for (const auto &p : mapping) {
            ui u = p.first;
            ui v = p.second;

            // 基本范围保护（可选，但很建议）
            if (u >= qn || v >= data_graph->getVerticesCount()) {
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

        // （可选但推荐）检查是否有未被赋值的 u（理论上不会）
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

        // 3. 检查缺失边 (Missing Edges Check)
        ui missing_edges = 0;
        for (ui u = 0; u < qn; ++u) {
            ui count;
            const ui *neighbors = query_graph->getVertexNeighbors(u, count);
            for (ui k = 0; k < count; ++k) {
                ui u_next = neighbors[k];
                // 为了避免无向图重复计算，只在 u < u_next 时检查
                if (u < u_next) {
                    int v = map_q[u];
                    int v_next = map_q[u_next];

                    // 如果查询图有边 (u, u_next)，但数据图没有边 (v, v_next)，则记为缺失
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
    } else {
        printf(">>> Verification FAILED: %zu/%lu unique results are invalid.\n",
               invalid_count, (unsigned long)unique_results.size());
    }
    printf("-------------------------\n");
}


// void check_correctness(const Graph *query_graph, const Graph *data_graph, 
//                        const vector<vector<pair<ui, ui>>> &M_ANS, ui threshold)
// {
//     printf("--- Verifying Results ---\n");
//     ui qn = query_graph->getVerticesCount();
//     size_t invalid_count = 0;

//     for (size_t i = 0; i < M_ANS.size(); ++i) {
//         const auto &mapping = M_ANS[i];
//         bool is_valid = true;
        
//         // 1. 映射大小检查
//         if (mapping.size() != qn) {
//             printf("Result %zu Error: Mapping size %zu != Query size %u\n", i, mapping.size(), qn);
//             invalid_count++;
//             continue;
//         }

//         // 构建快速查找表: map_q[u] -> v
//         vector<int> map_q(qn, -1);

//         // 2. 检查节点标签 (Label Check)
//         for (const auto &p : mapping) {
//             ui u = p.first;
//             ui v = p.second;
//             map_q[u] = v;

//             if (query_graph->getVertexLabel(u) != data_graph->getVertexLabel(v)) {
//                 printf("Result %zu Error: Label Mismatch u%u(L%u) -> v%u(L%u)\n", 
//                        i, u, query_graph->getVertexLabel(u), v, data_graph->getVertexLabel(v));
//                 is_valid = false;
//                 break;
//             }
//         }
//         if (!is_valid) {
//             invalid_count++;
//             continue;
//         }

//         // 3. 检查缺失边 (Missing Edges Check)
//         ui missing_edges = 0;
//         for (ui u = 0; u < qn; ++u) {
//             ui count;
//             const ui* neighbors = query_graph->getVertexNeighbors(u, count);
//             for (ui k = 0; k < count; ++k) {
//                 ui u_next = neighbors[k];
//                 // 为了避免无向图重复计算，只在 u < u_next 时检查
//                 if (u < u_next) {
//                     int v = map_q[u];
//                     int v_next = map_q[u_next];
                    
//                     // 如果查询图有边 (u, u_next)，但数据图没有边 (v, v_next)，则记为缺失
//                     if (!data_graph->hasEdge(v, v_next)) {
//                         missing_edges++;
//                     }
//                 }
//             }
//         }

//         if (missing_edges > threshold) {
//             printf("Result %zu Error: Missing Edges %u > Threshold %u\n", i, missing_edges, threshold);
//             is_valid = false;
//         }

//         if (!is_valid) invalid_count++;
//     }

//     if (invalid_count == 0) {
//         printf(">>> All %lu results PASSED verification.\n", M_ANS.size());
//     } else {
//         printf(">>> Verification FAILED: %zu/%lu results are invalid.\n", invalid_count, M_ANS.size());
//     }
//     printf("-------------------------\n");
// }

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
    long long matching_time = 0;
    long long treespan_time = 0;

    // query graph q, data graph G, 
    Graph *query_graph = new Graph();
    Graph *data_graph = new Graph();
    map<string, int> vM, eM;

    // load graphs
    t.restart();
    load_graph(data_graph_file, data_graph, vM, eM);
    load_graph(query_graph_file, query_graph, vM, eM);
    load_time = t.elapsed();

    vector<vector<pair<ui, ui> > > M_ANS;
    M_ANS.clear();

    // 1. matching
    t.restart();
    Approximate_Matching(query_graph, data_graph, M_ANS, threshold);
    matching_time = t.elapsed();

    printf("Our Matching Results:\n");
    printf("  count: %lu\n", M_ANS.size());
    printf("  Load Graphs Time: %.4lf ms\n", load_time / 1000.0);
    printf("  Matching Time: %.4lf ms\n", matching_time / 1000.0);
    printf("  Total Time: %.4lf ms\n", (load_time + matching_time) / 1000.0);

    // Check Correctness for Method 1
    if (!M_ANS.empty()) {
        check_correctness(query_graph, data_graph, M_ANS, threshold);
    }

    M_ANS.clear();

    // 2. Treespan
    t.restart();
    Approximate_TreeSpan(query_graph, data_graph, M_ANS, threshold);
    treespan_time = t.elapsed();
    
    printf("TreeSpan Results:\n");
    printf("  count: %lu\n", M_ANS.size());
    printf("  Load Graphs Time: %.4lf ms\n", load_time / 1000.0);
    printf("  Matching Time: %.4lf ms\n", treespan_time / 1000.0);
    printf("  Total Time: %.4lf ms\n", (load_time + treespan_time) / 1000.0);

    // Check Correctness for Method 2
    if (!M_ANS.empty()) {
        check_correctness(query_graph, data_graph, M_ANS, threshold);
    }

    return 0;
}