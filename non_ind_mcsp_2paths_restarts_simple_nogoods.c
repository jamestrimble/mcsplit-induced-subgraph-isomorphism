#include "graph.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <random>

#include <argp.h>
#include <limits.h>

using std::set;
using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

std::random_device rd;
std::mt19937 mt19937(rd());

const unsigned int MAX_NOGOOD_SIZE = 10;
const unsigned int MAX_NOGOOD_LIST_LEN = 20000;

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format";
static char args_doc[] = "FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"directed", 'i', 0, 0, "Use directed graphs"},
    {"enumerate", 'e', 0, 0, "Count solutions"},
    {"labelled", 'a', 0, 0, "Use edge and vertex labels"},
    {"vertex-labelled-only", 'x', 0, 0, "Use vertex labels, but not edge labels"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    { 0 }
};

static struct {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool directed;
    bool enumerate;
    bool edge_labelled;
    bool vertex_labelled;
    char *filename1;
    char *filename2;
    int timeout;
    int arg_num;
} arguments;

static std::atomic<bool> abort_due_to_timeout;

void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.directed = false;
    arguments.enumerate = false;
    arguments.edge_labelled = false;
    arguments.vertex_labelled = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 0;
    arguments.arg_num = 0;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'i':
            arguments.directed = true;
            break;
        case 'e':
            arguments.enumerate = true;
            break;
        case 'a':
            if (arguments.vertex_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.edge_labelled = true;
            arguments.vertex_labelled = true;
            break;
        case 'x':
            if (arguments.edge_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.vertex_labelled = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 1) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                     Stats
*******************************************************************************/

unsigned long long nodes{ 0 };

/*******************************************************************************
*******************************************************************************/

#define SMALL_ARR_SZ 16
class IntVec {
    int small_arr[SMALL_ARR_SZ];
    int sz;
    int max_capacity;
    int *vals;
public:
    IntVec(int max_capacity) : sz(0), max_capacity(max_capacity), vals(small_arr) {}
    ~IntVec() { if (sz > SMALL_ARR_SZ) delete[] vals; }
    IntVec(const IntVec& a) : sz(0), max_capacity(a.max_capacity), vals(small_arr) {
        std::cerr << "IntVec copy constructor" << std::endl;
        // TODO: make this more efficient
        for (const int x : a)
            push_back(x);
    }
    IntVec& operator=(const IntVec& a) {
        sz = a.sz;
        max_capacity = a.max_capacity;
        if (sz > SMALL_ARR_SZ) {
            vals = new int[max_capacity];
        } else {
            vals = small_arr;
        }
        for (int i=0; i<sz; i++)
            vals[i] = a.vals[i];
        return *this;
    }
    IntVec(IntVec&& a) : sz(a.sz), max_capacity(a.max_capacity) {
//        std::cerr << "IntVec move constructor" << std::endl;
        if (sz > SMALL_ARR_SZ) {
            vals = a.vals;
            a.vals = nullptr;
            a.sz = 0;
        } else {
            for (int i=0; i<sz; i++)
                small_arr[i] = a.small_arr[i];
            vals = small_arr;
        }
    }
    IntVec& operator=(IntVec&& a) {
//        std::cerr << "IntVec move assignment operator" << std::endl;
        if (this != &a) {
            max_capacity = a.max_capacity;
            if (sz > SMALL_ARR_SZ)
                delete[] vals;
            sz = a.sz;
            if (sz > SMALL_ARR_SZ) {
                vals = a.vals;
                a.vals = nullptr;
                a.sz = 0;
            } else {
                for (int i=0; i<sz; i++)
                    small_arr[i] = a.small_arr[i];
                vals = small_arr;
            }
        }
        return *this;
    }
    void push_back(int x) {
        if (sz == SMALL_ARR_SZ) {
            vals = new int[max_capacity];
            for (int i=0; i<SMALL_ARR_SZ; i++) {
                vals[i] = small_arr[i];
            }
        }
        vals[sz++] = x;
    }
    int *begin() { return vals; }
    int *begin() const { return vals; }
    int *cbegin() const { return vals; }
    int *end() { return vals + sz; }
    int *end() const { return vals + sz; }
    int *cend() const { return vals + sz; }
    int& operator[](std::size_t idx) { return vals[idx]; }
    const int& operator[](std::size_t idx) const { return vals[idx]; }
    int size() const { return sz; }
    void erase_vals(vector<unsigned char>& vals_to_erase) {
        int k=0;
        for (int i=0; i<sz; i++) {
            if (!vals_to_erase[vals[i]]) {
                vals[k] = vals[i];
                ++k;
            }
        }
        if (k <= SMALL_ARR_SZ && sz > SMALL_ARR_SZ) {
            for (int i=0; i<k; i++)
                small_arr[i] = vals[i];
            delete[] vals;
            vals = small_arr;
        }
        sz = k;
    }
    int get_max_capacity() const { return max_capacity; }
};

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

inline bool operator< (const VtxPair& lhs, const VtxPair& rhs){ return lhs.v<rhs.v || (lhs.v==rhs.v && lhs.w<rhs.w); }
inline bool operator> (const VtxPair& lhs, const VtxPair& rhs){ return rhs < lhs; }
inline bool operator<=(const VtxPair& lhs, const VtxPair& rhs){ return !(lhs > rhs); }
inline bool operator>=(const VtxPair& lhs, const VtxPair& rhs){ return !(lhs < rhs); }

typedef vector<set<vector<VtxPair>>> NogoodSets;

struct Bidomain {
    IntVec left_set;
    IntVec right_set;
    bool is_adjacent;
    int left_len() const { return left_set.size(); };
    int right_len() const { return right_set.size(); };
    Bidomain(IntVec left_set, IntVec right_set, bool is_adjacent):
            left_set(std::move(left_set)),
            right_set(std::move(right_set)),
            is_adjacent(is_adjacent) { };
    Bidomain(Bidomain&& a) : left_set(std::move(a.left_set)), right_set(std::move(a.right_set)),
            is_adjacent(a.is_adjacent)
    {
//        std::cerr << "Bidomain move constructor" << std::endl;
    }
    Bidomain& operator=(Bidomain&& a) {
//        std::cerr << "Bidomain move assignment operator" << std::endl;
        if (this != &a) {
            is_adjacent = a.is_adjacent;
            left_set = std::move(a.left_set);
            right_set = std::move(a.right_set);
        }
        return *this;
    }
};

void show(const vector<VtxPair>& current, const vector<Bidomain> &domains)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i=0; i<domains.size(); i++) {
        const struct Bidomain& bd = domains[i];
        cout << "Left  ";
        for (int v : bd.left_set)
            cout << v << " ";
        cout << std::endl;
        cout << "Right  ";
        for (int w : bd.right_set)
            cout << w << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
}

bool check_sol(const Graph & g0, const Graph & g1 , const vector<VtxPair> & solution) {
    return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i=0; i<solution.size(); i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j=i+1; j<solution.size(); j++) {
            struct VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w])
                return false;
        }
    }
    return true;
}

int calc_bound(const vector<Bidomain>& domains)
{
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += std::min(bd.left_len(), bd.right_len());
    }
    return bound;
}

// Returns false if we can backtrack
bool propagate_alldiff(vector<Bidomain>& domains, const Graph& g0, const Graph& g1)
{
    // This is really a vector of boolean values
    vector<unsigned char> vv1(g1.n, 0);

    int vv0_count = 0;
    int vv1_count = 0;
    vector<unsigned char> ww_to_erase;
    for (int i=0; i<(int)domains.size(); i++) {
        auto& bd = domains[i];
        if (ww_to_erase.size()) {
            bd.right_set.erase_vals(ww_to_erase);
//            bd.right_set.erase(std::remove_if(bd.right_set.begin(), bd.right_set.end(), 
//                        [&](int w){ return ww_to_erase[w]; }),
//                    bd.right_set.end());
        }

        vv0_count += bd.left_len();

        for (int w : bd.right_set) {
            if (!vv1[w]) {
                vv1[w] = true;
                vv1_count++;
            }
        }

        if (vv0_count > vv1_count)
            return false;
        else if (vv0_count == vv1_count)
            ww_to_erase = vv1;
    }
    return true;
}

bool assignment_impossible_by_2path_count(int v, int w, const vector<VtxPair>& current,
        const vector<vector<int>>& g0_2p, const vector<vector<int>>& g1_2p
        )
{
    for (auto pair : current) {
        if (g0_2p[v][pair.v] > g1_2p[w][pair.w]) {
            return true;
        }
    }
    return false;
}

// Returns a pair, whose first value is the min domain size, and whose second value is a tie-breaker
// giving the lowest vertex index of a vertex whose domain size is minimal
std::pair<int, int> bidomain_score(
        const Bidomain &bd,
        const vector<vector<int>>& g0_2p, const vector<vector<int>>& g1_2p,
        const vector<int>& g0_deg,
        const vector<int>& g1_deg,
        const vector<VtxPair>& current,
        std::pair<int, int> incumbent)
{
    auto best = incumbent;
    for (int v : bd.left_set) {
        auto vtx_score = std::make_pair(0, v);
        for (int w : bd.right_set) {
            if (g0_deg[v] <= g1_deg[w] && !assignment_impossible_by_2path_count(v, w, current, g0_2p, g1_2p)) {
                vtx_score.first++;
                if (vtx_score > best)
                    break;
            }
        }
        if (vtx_score < best) {
            best = vtx_score;
        }
    }
    return best;
}

std::pair<int, int> select_bidomain_and_branching_var(const vector<Bidomain>& domains,
        const vector<vector<int>>& g0_2p, const vector<vector<int>>& g1_2p,
        const vector<int>& g0_deg,
        const vector<int>& g1_deg,
        const vector<VtxPair>& current)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    auto best_score = std::make_pair(INT_MAX, INT_MAX);
    int best = -1;
    for (int i=0; i<(int)domains.size(); i++) {
        const Bidomain& bd = domains[i];
        auto score = bidomain_score(bd, g0_2p, g1_2p, g0_deg, g1_deg, current, best_score);
        if (score < best_score) {
            if (score.first == 0)
                return std::make_pair(-1, -1);
            best_score = score;
            best = i;
            if (best_score.first == 1)
                break;
        }
    }
    return std::make_pair(best, best_score.second);
}

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains(const vector<Bidomain> & d,
        const Graph & g0, const Graph & g1, int v, int w)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain& old_bd : d) {
        // left_with_edge is the set of vertices (not including v) with edges to v
        IntVec left_with_edge(old_bd.left_set.get_max_capacity());
//        left_with_edge.reserve(old_bd.left_set.size());
        // left_without_edge is the set of vertices (not including v) without any edges to v
        IntVec left_without_edge(old_bd.left_set.get_max_capacity());
//        left_without_edge.reserve(old_bd.left_set.size());
        for (int u : old_bd.left_set) {
            if (u != v) {
                if (g0.adjmat[v][u]) {
                    left_with_edge.push_back(u);
                } else {
                    left_without_edge.push_back(u);
                }
            }
        }

        // right_with_edge is the set of vertices (not including w) with edges to w
        IntVec right_with_edge(old_bd.right_set.get_max_capacity());
//        right_with_edge.reserve(old_bd.right_set.size());
        // right_without_w is the right set with w removed
        IntVec right_without_w(old_bd.right_set.get_max_capacity());
//        right_without_w.reserve(old_bd.right_set.size());
        for (int u : old_bd.right_set) {
            if (u != w) {
                right_without_w.push_back(u);
                if (g1.adjmat[w][u]) {
                    right_with_edge.push_back(u);
                }
            }
        }

        if ((left_without_edge.size() > right_without_w.size()) || (left_with_edge.size() > right_with_edge.size())) {
            // Stop early if we know that there are vertices in the first graph that can't be matched
            return {};
        }
        if (left_without_edge.size() && right_without_w.size())
            new_d.push_back({std::move(left_without_edge), std::move(right_without_w), old_bd.is_adjacent});
        if (left_with_edge.size() && right_with_edge.size())
            new_d.push_back({std::move(left_with_edge), std::move(right_with_edge), true});
    }

    std::sort(new_d.begin(), new_d.end(),
            [](const Bidomain& a, const Bidomain& b) { return a.right_len() < b.right_len(); });

    return new_d;
}

#define FOUND_SOLUTION 0
#define NO_SOLUTION 1
#define TIMEOUT 2
#define REACHED_NODE_LIMIT 3

int solve(const Graph & g0, const Graph & g1,
        const vector<int>& g0_deg, const vector<int>& g1_deg,
        const vector<vector<int>> & g0_2p, const vector<vector<int>> & g1_2p,
        vector<VtxPair> & incumbent, vector<VtxPair> & current, vector<Bidomain> & domains,
        long long & solution_count, long long& node_count, long long node_limit,
        NogoodSets& nogoods)
{
//    cout << node_count << " "<< node_limit << endl;

    if (node_count >= node_limit) {
//        cout << "!" << endl;
        return REACHED_NODE_LIMIT;
    }

    if (abort_due_to_timeout)
        return TIMEOUT;

    if (arguments.verbose) show(current, domains);
    nodes++;

    if (current.size() > incumbent.size())
        incumbent = current;

    if (current.size()==(unsigned)g0.n) {
        solution_count++;
        return FOUND_SOLUTION;
    }

    if (!arguments.enumerate && incumbent.size()==(unsigned)g0.n)
        return FOUND_SOLUTION;

    if ((int)current.size() + calc_bound(domains) < g0.n)
        return NO_SOLUTION;

    if (!propagate_alldiff(domains, g0, g1))
        return NO_SOLUTION;

    auto bd_idx_and_v = select_bidomain_and_branching_var(domains,
            g0_2p, g1_2p, g0_deg, g1_deg, current);
    int bd_idx = bd_idx_and_v.first;
    int v = bd_idx_and_v.second;

    if (bd_idx == -1)   // Return if there's nothing left to branch on
        return NO_SOLUTION;
    Bidomain &bd = domains[bd_idx];

    std::shuffle(std::begin(bd.right_set), std::end(bd.right_set), mt19937);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    for (int w : bd.right_set) {
        if (g0_deg[v] <= g1_deg[w] && !assignment_impossible_by_2path_count(v, w, current, g0_2p, g1_2p)) {
            auto new_domains = filter_domains(domains, g0, g1, v, w);
            current.push_back(VtxPair(v, w));
            if (current.size() > MAX_NOGOOD_SIZE || nogoods[current.size()].count(current)==0) {
                int result = solve(g0, g1, g0_deg, g1_deg, g0_2p, g1_2p, incumbent, current, new_domains, solution_count,
                        node_count, node_limit, nogoods);
                if (result == FOUND_SOLUTION)
                    return FOUND_SOLUTION;
                if (result==NO_SOLUTION && current.size()<=MAX_NOGOOD_SIZE && nogoods[current.size()].size()<MAX_NOGOOD_LIST_LEN)
                    nogoods[current.size()].insert(current);
            }
            current.pop_back();
            if (node_count >= node_limit)
                return REACHED_NODE_LIMIT;
        }
    }
    node_count++;
    return NO_SOLUTION;
}

// TODO: change values from negative to positive
vector<vector<int>> count_2paths(const Graph & g)
{
    vector<vector<int>> num_2paths(g.n, vector<int>(g.n, 0));
    for (int i=0; i<g.n; i++) {
        for (int k=0; k<g.n; k++) {
            if (g.adjmat[i][k]) {
                for (int j=0; j<i; j++) {
                    if (g.adjmat[k][j]) {
                        num_2paths[i][j]++;
                        num_2paths[j][i]++;
                    }
                }
            }
        }
    }
    return num_2paths;
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++)
        for (int w=0; w<g.n; w++)
            if (g.adjmat[v][w]) degree[v]++;
    return degree;
}

// Returns a common subgraph and the number of induced subgraph isomorphisms found
std::pair<vector<VtxPair>, long long> mcs(const Graph & g0, const Graph & g1)
{
    vector<vector<int>> g0_2p = count_2paths(g0);
    vector<vector<int>> g1_2p = count_2paths(g1);

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);

    vector<int> g0_deg_sorted = g0_deg;
    vector<int> g1_deg_sorted = g1_deg;
    std::sort(g0_deg_sorted.begin(), g0_deg_sorted.end(), std::greater<int>());
    std::sort(g1_deg_sorted.begin(), g1_deg_sorted.end(), std::greater<int>());
    for (int i=0; i<(int)g0_deg_sorted.size(); i++) {
        if (g1_deg_sorted[i] < g0_deg_sorted[i]) {
            return {{}, 0};
        }
    }

    auto domains = vector<Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    for (int is_isolated=0; is_isolated<=1; is_isolated++) {
        // Create a bidomain for each label that appears in both graphs
        for (unsigned int label : labels) {
            IntVec left_set(g0.n);
            IntVec right_set(g1.n);

            for (int i=0; i<g0.n; i++)
                if (g0.label[i]==label && is_isolated==(g0_deg[i]==0))
                    left_set.push_back(i);
            for (int i=0; i<g1.n; i++)
                if (g1.label[i]==label && (is_isolated || g1_deg[i]>0))
                    right_set.push_back(i);

            if (left_set.size() && right_set.size())
                domains.push_back({std::move(left_set), std::move(right_set), false});
        }
    }

    vector<VtxPair> incumbent;
    vector<VtxPair> current;
    long long solution_count = 0;

    NogoodSets nogoods(MAX_NOGOOD_SIZE+1);

    long long node_count;
    long long node_limit_multiplier = 100;
    int num_attempts = 0;
    vector<long long> luby_seq = {1};
    int luby_seq_pos = 0;
    long long node_limit;
    int result;
    do {
        if (luby_seq_pos == (int) luby_seq.size()) {
            int luby_seq_len = luby_seq.size();
            for (int i=0; i<luby_seq_len; i++) {
                luby_seq.push_back(luby_seq[i]);
            }
            luby_seq.push_back(luby_seq.back() * 2);
        }
        node_limit = luby_seq[luby_seq_pos++] * node_limit_multiplier;
        node_count = 0;
//        cout << "Number of restarts: " << num_attempts << endl;
//        cout << "   Limit: " << node_limit << endl;
        num_attempts++;
//        std::cout << incumbent.size() << " " << g0.n << endl;
//        for (unsigned int i=0; i<=MAX_NOGOOD_SIZE; i++)
//            std::cout << " * " << i << " " << nogoods[i].size() << endl;
        result = solve(g0, g1, g0_deg, g1_deg, g0_2p, g1_2p, incumbent, current, domains, solution_count,
                node_count, node_limit, nogoods);
    } while (result == REACHED_NODE_LIMIT);
    
//    for (unsigned int i=0; i<=MAX_NOGOOD_SIZE; i++)
//        std::cout << " * " << i << " " << nogoods[i].size() << endl;

    cout << "Number of restarts: " << (num_attempts-1) << endl;

    return {incumbent, solution_count};
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    struct Graph g0 = readGraph(arguments.filename1, format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled);
    struct Graph g1 = readGraph(arguments.filename2, format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled);

    if (g0.n > g1.n) {
        std::cout << "Error: pattern graph has more vertices than target graph." << std::endl;
        return 1;
    }

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
                auto abort_time = std::chrono::steady_clock::now() + std::chrono::seconds(arguments.timeout);
                {
                    /* Sleep until either we've reached the time limit,
                     * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
    }

    auto start = std::chrono::steady_clock::now();

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);

    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g0_deg[a]>g0_deg[b];
    });
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
        return g1_deg[a]>g1_deg[b];
    });

    struct Graph g0_sorted = induced_subgraph(g0, vv0);
    struct Graph g1_sorted = induced_subgraph(g1, vv1);

    auto result = mcs(g0_sorted, g1_sorted);
    vector<VtxPair> solution = result.first;
    long long num_sols = result.second;

    // Convert to indices from original, unsorted graphs
    for (auto& vtx_pair : solution) {
        vtx_pair.v = vv0[vtx_pair.v];
        vtx_pair.w = vv1[vtx_pair.w];
    }

    auto stop = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }


    cout << "Nodes:                      " << nodes << endl;
    cout << "CPU time (ms):              " << time_elapsed << endl;
    if (aborted) {
        cout << "TIMEOUT" << endl;
    } else {
        if (!check_sol(g0, g1, solution))
            fail("*** Error: Invalid solution\n");

        if (arguments.enumerate) {
            std::cout << "Number of solutions: " << num_sols << std::endl;
        }
        if ((int)solution.size() == std::min(g0.n, g1.n)) {
            cout << "Solution size " << solution.size() << std::endl;
            std::cout << "SATISFIABLE" << std::endl;
            for (int i=0; i<g0.n; i++)
                for (unsigned int j=0; j<solution.size(); j++)
                    if (solution[j].v == i)
                        cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
            cout << std::endl;
        } else {
            std::cout << "UNSATISFIABLE" << std::endl;
        }
    }
}

