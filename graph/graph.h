#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "utility/utility.h"
#include "configuration/types.h"
#include "configuration/config.h"

class Graph
{
private:
	std::string graph_id;	// graph id
	ui 			n;			// number of vertices
	ept 		m;			// number of edges
	ept *pstart;	// pstart[u]: starting position of vertex u in edges array
	ui *edges;		// edges array
	LabelID *vlabels;	// vertex labels array

	ui			_max_degree; // maximum degree of the graph
	// std::unordered_map<LabelID, ui> vertex_labels_frequency; // frequency of each vertex label
	ui          labels_count;
public:
	Graph()
	{
		graph_id = "";
		n = 0;
		m = 0;
		pstart = nullptr;
		edges = nullptr;
		vlabels = nullptr;
		_max_degree = 0;
		labels_count = 0;
	}

	~Graph()
	{
		if (pstart != nullptr) { delete[] pstart; pstart = nullptr; }
		if (edges != nullptr) { delete[] edges; edges = nullptr; }
		if (vlabels != nullptr) { delete[] vlabels; vlabels = nullptr; }
	}

	void build_graph(const std::string &id,
		std::vector<std::pair<ui, LabelID> > &vertices,
		std::vector<std::pair<std::pair<ui, ui>, LabelID> > &edges_list);

	void print_graph() const;

	ui getVerticesCount() const
	{
		return n;
	}

	// directed edges count
	ept getEdgesCount() const
	{
		return m;
	}

	ui getVertexDegree(ui u) const
	{
		assert(u >= 0 && u < n);
		return pstart[u + 1] - pstart[u];
	}

	ui getMaxDegree() const
	{
		return _max_degree;
	}

	LabelID getVertexLabel(ui u) const
	{
		assert(u >= 0 && u < n);
		return vlabels[u];
	}

	ui getLabelsCount() const
	{
		return labels_count;
	}

	const ui *getVertexNeighbors(const ui id, ui &count) const
	{
		count = pstart[id + 1] - pstart[id];
		return edges + pstart[id];
	}

	// Check if there is an edge u -> v
	bool hasEdge(ui u, ui v) const
	{
		assert(u < n && v < n);
		ept l = pstart[u], r = pstart[u + 1];
		// binary search in sorted adjacency list
		return std::binary_search(edges + l, edges + r, v);
	}
};

#endif