#include "graph.h"

void Graph::print_graph() const {
    std::cout << "=========================================\n";
    std::cout << " Graph ID: " << graph_id << "\n";
    std::cout << "-----------------------------------------\n";

    std::cout << "Vertices (" << n << "):\n";
    // for (ui u = 0; u < n; ++u) {
    //     std::cout << "  v " << u << " " << vlabels[u] << "\n";
    // }

    std::cout << "-----------------------------------------\n";

    std::cout << "Edges (" << m << "):\n";
//     for (ui u = 0; u < n; ++u) {
//         for (ept eid = pstart[u]; eid < pstart[u + 1]; ++eid) {
//             ui v = edges[eid];
//             std::cout << "  e " << u << " " << v << "\n";
//         }
//     }

    std::cout << "=========================================\n\n";
}

void Graph::build_graph(const std::string &id,
                        std::vector<std::pair<ui, LabelID> > &vertices, 
                        std::vector<std::pair<std::pair<ui, ui>, LabelID> > &edges_list)
{
    graph_id = id;
    n = vertices.size();
    m = edges_list.size();

    pstart = new ept[n + 1];
    edges = new ui[m];
    vlabels = new LabelID[n];

    std::sort(vertices.begin(), vertices.end());
    std::sort(edges_list.begin(), edges_list.end());

#ifndef NDEBUG
    ui last_u = 0, last_v = 0;
    for(ept i = 0; i < m; i++) {
        ui u = edges_list[i].first.first;
        ui v = edges_list[i].first.second;
        assert(u >= 0 && u < n);
        assert(v >= 0 && v < n);
        assert(u > last_u || (u == last_u && v > last_v));
        last_u = u;
        last_v = v;
    }
#endif

    for(ui i = 0; i < n; i++) {
        assert(vertices[i].first == i);
        vlabels[vertices[i].first] = vertices[i].second;
        assert(vertices[i].second >= 0);
        if ((ui)vertices[i].second >= labels_count) {
            labels_count = vertices[i].second + 1;
        }
    }
    for(ept i = 0; i < m; i++) {
        edges[i] = edges_list[i].first.second;
    }

    ept idx = 0;
    pstart[0] = idx;
    for(ui i = 0; i < n; i++) {
        while(idx < m && edges_list[idx].first.first == i) idx++;
        pstart[i + 1] = idx;
    }
    assert(pstart[n] == m);

    _max_degree = 0;
    for(ui i = 0; i < n; i++) {
        ui degree = getVertexDegree(i);
        if(degree > _max_degree) {
            _max_degree = degree;
        }
    }
}
