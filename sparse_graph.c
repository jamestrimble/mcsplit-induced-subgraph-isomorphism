#include "sparse_graph.h"

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_set>

constexpr int BITS_PER_UNSIGNED_INT (CHAR_BIT * sizeof(unsigned int));

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

SparseGraph::SparseGraph(unsigned int n) : n(n), adj_lists(n), in_edge_lists(n), label(n) { }

void add_edge(SparseGraph& g, int v, int w, bool directed=false, unsigned int val=1) {
    if (v != w) {
        if (directed) {
            g.adj_lists[v].push_back(w);
            g.in_edge_lists[w].push_back(v);
        } else {
            g.adj_lists[v].push_back(w);
            g.adj_lists[w].push_back(v);
        }
    } else {
        // To indicate that a vertex has a loop, we set the most
        // significant bit of its label to 1
        g.label[v] |= (1u << (BITS_PER_UNSIGNED_INT-1));
    }
}

SparseGraph induced_subgraph(struct SparseGraph& g, std::vector<int> vv) {
    std::vector<int> vv_inv(vv.size(), -1);
    for (int i=0, sz=vv.size(); i<sz; i++) {
        vv_inv[vv[i]] = i;
    }
    SparseGraph subg(vv.size());
    for (int i=0, sz=vv.size(); i<sz; i++) {
        int v = vv[i];
        for (int w : g.adj_lists[v]) {
            int w_inv = vv_inv[w];
            if (w_inv != -1) subg.adj_lists[i].push_back(vv_inv[w]);
        }
    }
    for (int i=0, sz=vv.size(); i<sz; i++) {
        int v = vv[i];
        for (int w : g.in_edge_lists[v]) {
            int w_inv = vv_inv[w];
            if (w_inv != -1) subg.in_edge_lists[i].push_back(vv_inv[w]);
        }
    }

    for (int i=0; i<subg.n; i++)
        subg.label[i] = g.label[vv[i]];
    return subg;
}

bool SparseGraph::has_edge(int v, int w) const
{
    auto & lst = this->adj_lists[v];
    return std::find(lst.begin(), lst.end(), w) != lst.end();
}

bool SparseGraph::has_filtered_edge(int v, int w) const
{
    auto & lst = this->filtered_adj_lists[v];
    return std::find(lst.begin(), lst.end(), w) != lst.end();
}

struct SparseGraph readDimacsGraph(char* filename, bool directed, bool vertex_labelled) {
    struct SparseGraph g(0);
    return g;
}

struct SparseGraph readLadGraph(char* filename, bool directed) {
    struct SparseGraph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int w;

    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = SparseGraph(nvertices);

    // If edge (v,w) has been seen, then edges_seen will contain v * nvertices + w.
    std::unordered_set<unsigned long long> edges_seen;

    for (int i=0; i<nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)
            fail("Number of edges not read correctly.\n");
        if (directed) {
            for (int j=0; j<edge_count; j++) {
                if (fscanf(f, "%d", &w) != 1)
                    fail("An edge was not read correctly.\n");
                if (edges_seen.find(i * nvertices + w) == edges_seen.end()) {
                    edges_seen.insert(i * nvertices + w);
                    add_edge(g, i, w, true);
                }
            }
        } else {
            for (int j=0; j<edge_count; j++) {
                if (fscanf(f, "%d", &w) != 1)
                    fail("An edge was not read correctly.\n");
                int lo = std::min(i, w);
                int hi = std::max(i, w);
                if (edges_seen.find(lo * nvertices + hi) == edges_seen.end()) {
                    edges_seen.insert(lo * nvertices + hi);
                    add_edge(g, i, w, false);
                }
            }
        }
    }

    fclose(f);
    return g;
}

struct SparseGraph readGfdGraph(char* filename) {
    bool directed = true;

    struct SparseGraph g(0);
    FILE* f;

    char junk[2048];
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int v, w;

    if (fscanf(f, "%s", junk) != 1)
        fail("File name not read correctly.\n");
    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = SparseGraph(nvertices);

    // If edge (v,w) has been seen, then edges_seen will contain v * nvertices + w.
    std::unordered_set<unsigned long long> edges_seen;

    for (int i=0; i<nvertices; i++) {
        int label;
        if (fscanf(f, "%d", &label) != 1)
            fail("Label not read correctly.\n");
        g.label[i] = label;
    }
    int edge_count;
    if (fscanf(f, "%d", &edge_count) != 1)
        fail("Number of edges not read correctly.\n");
    for (int j=0; j<edge_count; j++) {
        if (fscanf(f, "%d %d", &v, &w) != 2)
            fail("An edge was not read correctly.\n");
        if (edges_seen.find(v * nvertices + w) == edges_seen.end()) {
            edges_seen.insert(v * nvertices + w);
            add_edge(g, v, w, directed);
        }
    }

    fclose(f);
    return g;
}

struct SparseGraph readVfGraph(char* filename) {
    bool directed = true;

    struct SparseGraph g(0);
    FILE* f;

    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int v, w;

    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = SparseGraph(nvertices);

    // If edge (v,w) has been seen, then edges_seen will contain v * nvertices + w.
    std::unordered_set<unsigned long long> edges_seen;

    for (int i=0; i<nvertices; i++) {
        int label;
        if (fscanf(f, "%d %d", &v, &label) != 2)
            fail("Label not read correctly.\n");
        g.label[v] = label;
    }
    for (int i=0; i<nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)
            fail("Number of edges not read correctly.\n");
        for (int j=0; j<edge_count; j++) {
            if (fscanf(f, "%d %d", &v, &w) != 2)
                fail("An edge was not read correctly.\n");
            if (edges_seen.find(v * nvertices + w) == edges_seen.end()) {
                edges_seen.insert(v * nvertices + w);
                add_edge(g, v, w, directed);
            }
        }
    }

    fclose(f);
    return g;
}

struct SparseGraph readBinaryGraph(char* filename, bool directed, bool edge_labelled,
        bool vertex_labelled)
{
    struct SparseGraph g(0);
    return g;
}

struct SparseGraph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled) {
    struct SparseGraph g(0);
    if (format=='G') g = readGfdGraph(filename);
    else if (format=='V') g = readVfGraph(filename);
    else if (format=='D') g = readDimacsGraph(filename, directed, vertex_labelled);
    else if (format=='L') g = readLadGraph(filename, directed);
    else if (format=='B') g = readBinaryGraph(filename, directed, edge_labelled, vertex_labelled);
    else fail("Unknown graph format\n");
    return g;
}

void filter_adj_lists(SparseGraph & g, const std::vector<bool> & active_vertices) {
    g.filtered_adj_lists = std::vector<std::vector<int>>(g.n);
    for (int v=0; v<g.n; v++) {
        for (int w : g.adj_lists[v]) {
            if (active_vertices[w]) {
                g.filtered_adj_lists[v].push_back(w);
            }
        }
    }
}
