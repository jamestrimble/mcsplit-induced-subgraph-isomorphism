#include <limits.h>

#include <vector>

struct Graph {
    int n;
    std::vector<std::vector<int>> adj_lists;
    std::vector<std::vector<int>> filtered_adj_lists;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
    bool has_edge(int v, int w) const;
    bool has_filtered_edge(int v, int w) const;
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

void filter_adj_lists(Graph & g, const std::vector<bool> & active_vertices);
