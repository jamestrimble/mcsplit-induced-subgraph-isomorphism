#include "graph.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>

#include <argp.h>
#include <limits.h>

using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

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
    {"opposite-value-heuristic", 'o', 0, 0, "Use opposite value heuristic"},
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
    bool opposite_value_heuristic;
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
    arguments.opposite_value_heuristic = false;
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
        case 'o':
            arguments.opposite_value_heuristic = true;
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

#define POOL_SIZE (1 << 17)

struct Pool {
    int num_chunks_outstanding;
    int *vals;
    bool active_for_allocation;
    Pool() : num_chunks_outstanding(0), active_for_allocation(true) {
        vals = new int[POOL_SIZE];
    }
    ~Pool() {
        delete[] vals;
    }
};

struct PoolChunk {
    int *vals;
    Pool *pool;
    void release() {
        pool->num_chunks_outstanding--;
        if (!pool->num_chunks_outstanding && !pool->active_for_allocation)
            delete pool;
    }
};

class PoolAllocator {
    Pool *pool;
    int num_used;
public:
    PoolAllocator() : num_used(0) {
        pool = new Pool();
    }
    PoolChunk get_chunk(int size) {
        if (num_used + size > POOL_SIZE) {
            pool->active_for_allocation = false;
            pool = new Pool();
            num_used = 0;
        }
        int *start = &pool->vals[num_used];
        pool->num_chunks_outstanding++;
        num_used += size;
        return { start, pool };
    }
};

PoolAllocator pool_allocator;

/*******************************************************************************
*******************************************************************************/

class IntVec {
    PoolChunk pool_chunk;
    int sz;
    int max_capacity;
public:
    IntVec(int max_capacity) : pool_chunk(pool_allocator.get_chunk(max_capacity)), sz(0), max_capacity(max_capacity) {}
    ~IntVec() {
        pool_chunk.release();
    }
    IntVec(const IntVec& a) : pool_chunk(pool_allocator.get_chunk(a.max_capacity)), sz(0), max_capacity(a.max_capacity) {
        std::cerr << "IntVec copy constructor" << std::endl;
        // TODO: make this more efficient
        for (const int x : a)
            push_back(x);
    }
    friend void swap(IntVec& first, IntVec& second) {
        using std::swap;
        swap(first.sz, second.sz);
        swap(first.pool_chunk, second.pool_chunk);
    }
    IntVec& operator=(IntVec a) {
        swap(*this, a);
        return *this;
    }
    IntVec(IntVec&& a) : IntVec(0) {
        swap(*this, a);
    }
    void push_back(int x) {
        pool_chunk.vals[sz++] = x;
    }
    int *begin() { return pool_chunk.vals; }
    int *begin() const { return pool_chunk.vals; }
    int *cbegin() const { return pool_chunk.vals; }
    int *end() { return pool_chunk.vals + sz; }
    int *end() const { return pool_chunk.vals + sz; }
    int *cend() const { return pool_chunk.vals + sz; }
    int& operator[](std::size_t idx) { return pool_chunk.vals[idx]; }
    const int& operator[](std::size_t idx) const { return pool_chunk.vals[idx]; }
    int size() const { return sz; }
    void erase_vals(vector<unsigned char>& vals_to_erase) {
        int k=0;
        for (int i=0; i<sz; i++) {
            if (!vals_to_erase[pool_chunk.vals[i]]) {
                pool_chunk.vals[k] = pool_chunk.vals[i];
                ++k;
            }
        }
        sz = k;
    }
};

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

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
    friend void swap(Bidomain& first, Bidomain& second) {
        using std::swap;
        swap(first.left_set, second.left_set);
        swap(first.right_set, second.right_set);
        swap(first.is_adjacent, second.is_adjacent);
    }
    Bidomain& operator=(Bidomain a) {
        swap(*this, a);
        return *this;
    }
};

struct UsefulStuff {
    vector<int> g0_deg;
    vector<int> g1_deg;
    vector<vector<int>> g0_2p;
    vector<vector<int>> g1_2p;
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
        if (ww_to_erase.size())
            bd.right_set.erase_vals(ww_to_erase);

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
        const UsefulStuff& useful_stuff,
        const vector<VtxPair>& current,
        std::pair<int, int> incumbent)
{
    auto best = incumbent;
    for (int v : bd.left_set) {
        auto vtx_score = std::make_pair(0, v);
        for (int w : bd.right_set) {
            if (useful_stuff.g0_deg[v] <= useful_stuff.g1_deg[w] &&
                    !assignment_impossible_by_2path_count(v, w, current, useful_stuff.g0_2p, useful_stuff.g1_2p)) {
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
        const UsefulStuff& useful_stuff, const vector<VtxPair>& current)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    auto best_score = std::make_pair(INT_MAX, INT_MAX);
    int best = -1;
    for (int i=0; i<(int)domains.size(); i++) {
        const Bidomain& bd = domains[i];
        auto score = bidomain_score(bd, useful_stuff, current, best_score);
        if (score < best_score) {
            if (score.first == 0)
                return std::make_pair(-1, -1);
            best_score = score;
            best = i;
        }
    }
    return std::make_pair(best, best_score.second);
}

vector<Bidomain> filter_domains(const vector<Bidomain> & d,
        const Graph & g0, const Graph & g1, int v, int w)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());

    auto& v_adjrow = g0.adjmat[v];
    auto& w_adjrow = g1.adjmat[w];

    for (const Bidomain& old_bd : d) {
        // left_with_edge is the set of vertices (not including v) with edges to v
        IntVec left_with_edge(old_bd.left_set.size());
        // left_without_edge is the set of vertices (not including v) without any edges to v
        IntVec left_without_edge(old_bd.left_set.size());
        for (int u : old_bd.left_set) {
            if (u != v) {
                if (v_adjrow[u]) {
                    left_with_edge.push_back(u);
                } else {
                    left_without_edge.push_back(u);
                }
            }
        }

        // right_with_edge is the set of vertices (not including w) with edges to w
        IntVec right_with_edge(old_bd.right_set.size());
        // right_without_w is the right set with w removed
        IntVec right_without_w(old_bd.right_set.size());

        // The next bit is could be done more simply at the cost of a little bit of efficiency;
        // I've tried to reduce the number of checks of whether u == w
        if (left_without_edge.size() && left_with_edge.size()) {
            int *p;
            int *end = old_bd.right_set.end();
            for (p=old_bd.right_set.begin(); p<end; p++) {
                int u = *p;
                if (u == w) break;
                right_without_w.push_back(u);
                if (w_adjrow[u])
                    right_with_edge.push_back(u);
            }
            p++;
            for (; p<end; p++) {
                int u = *p;
                right_without_w.push_back(u);
                if (w_adjrow[u])
                    right_with_edge.push_back(u);
            }
        } else if (left_without_edge.size()) {
            int *p;
            int *end = old_bd.right_set.end();
            for (p=old_bd.right_set.begin(); p<end; p++) {
                int u = *p;
                if (u == w) break;
                right_without_w.push_back(*p);
            }
            p++;
            for (; p<end; p++) {
                right_without_w.push_back(*p);
            }
        } else if (left_with_edge.size()) {
            for (int u : old_bd.right_set)
                if (w_adjrow[u])
                    right_with_edge.push_back(u);
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
            [](const Bidomain& a, const Bidomain& b) {
                return a.right_len() < b.right_len() ||
                        (a.right_len()==b.right_len() && a.left_set[0] < b.left_set[0]);
            });

    return new_d;
}

void solve(const Graph & g0, const Graph & g1,
        const UsefulStuff& useful_stuff,
        vector<VtxPair> & incumbent, vector<VtxPair> & current, vector<Bidomain> & domains,
        long long & solution_count)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(current, domains);
    nodes++;

    if (current.size() > incumbent.size())
        incumbent = current;

    if (current.size()==(unsigned)g0.n) {
        solution_count++;
        return;
    }

    if (!domains.size())
        return;

    if (!arguments.enumerate && incumbent.size()==(unsigned)g0.n)
        return;

    if (!propagate_alldiff(domains, g0, g1))
        return;

    auto bd_idx_and_v = select_bidomain_and_branching_var(domains, useful_stuff, current);
    int bd_idx = bd_idx_and_v.first;
    int v = bd_idx_and_v.second;

    if (bd_idx == -1)   // Return if there's nothing left to branch on
        return;
    Bidomain &bd = domains[bd_idx];

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    for (int i=0; i<bd.right_len(); i++) {
        int w = bd.right_set[i];

        if (useful_stuff.g0_deg[v] <= useful_stuff.g1_deg[w] &&
                !assignment_impossible_by_2path_count(v, w, current, useful_stuff.g0_2p, useful_stuff.g1_2p)) {
            auto new_domains = filter_domains(domains, g0, g1, v, w);
            current.push_back(VtxPair(v, w));
            solve(g0, g1, useful_stuff, incumbent, current, new_domains, solution_count);
            current.pop_back();
        }
    }
}

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

            int left_min_deg = INT_MAX;
            for (int i=0; i<g0.n; i++) {
                if (g0.label[i]==label && is_isolated==(g0_deg[i]==0)) {
                    left_set.push_back(i);
                    int deg = g0_deg[i];
                    if (deg < left_min_deg)
                        left_min_deg = deg;
                }
            }
            for (int i=0; i<g1.n; i++)
                if (g1.label[i]==label && (is_isolated || g1_deg[i]>0) && g1_deg[i]>=left_min_deg)
                    right_set.push_back(i);

            if (left_set.size() > right_set.size())
                return {{}, 0};

            if (left_set.size())
                domains.push_back({std::move(left_set), std::move(right_set), false});
        }
    }

    vector<VtxPair> incumbent;
    vector<VtxPair> current;
    long long solution_count = 0;

    UsefulStuff useful_stuff = {std::move(g0_deg), std::move(g1_deg), std::move(g0_2p), std::move(g1_2p)};
    solve(g0, g1, useful_stuff, incumbent, current, domains, solution_count);

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
    if (arguments.opposite_value_heuristic)
        std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) { return g1_deg[a]<g1_deg[b]; });
    else
        std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) { return g1_deg[a]>g1_deg[b]; });

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

