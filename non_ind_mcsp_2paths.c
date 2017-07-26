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
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

struct Bidomain {
    int l;        // start index of left set
    int left_len;
    vector<int> right_set;
    bool is_adjacent;
    int right_len() const { return right_set.size(); };
    Bidomain(int l, int left_len, vector<int> right_set, bool is_adjacent):
            l(l),
            left_len(left_len),
            right_set(right_set),
            is_adjacent(is_adjacent) { };
};

void show(const vector<VtxPair>& current, const vector<Bidomain> &domains,
        const vector<int>& left)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i=0; i<domains.size(); i++) {
        struct Bidomain bd = domains[i];
        cout << "Left  ";
        for (int j=0; j<bd.left_len; j++)
            cout << left[bd.l + j] << " ";
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
        bound += std::min(bd.left_len, bd.right_len());
    }
    return bound;
}

// Returns false if we can backtrack
bool propagate_alldiff(vector<Bidomain>& domains,
        const Graph& g0, const Graph& g1, vector<int>& left)
{
    // This is really a vector of boolean values
    vector<unsigned char> vv1(g1.n, 0);

    int vv0_count = 0;
    int vv1_count = 0;
    vector<unsigned char> ww_to_erase;
    for (int i=0; i<(int)domains.size(); i++) {
        auto& bd = domains[i];
        if (ww_to_erase.size()) {
            bd.right_set.erase(std::remove_if(bd.right_set.begin(), bd.right_set.end(), 
                        [&](int w){ return ww_to_erase[w]; }),
                    bd.right_set.end());
        }

        vv0_count += bd.left_len;

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

int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

bool assignment_impossible_by_2path_count(int v, int w, const vector<VtxPair>& current,
        const vector<vector<int>>& g0_2p, const vector<vector<int>>& g1_2p
        )
{
    for (auto pair : current) {
        if (g0_2p[v][pair.v] < g1_2p[w][pair.w]) {
        //if (g0_2p.at(v).at(pair.v) < g1_2p.at(w).at(pair.w)) {
            return true;
        }
    }
    return false;
}

// Returns a pair, whose first value is the min domain size, and whose second value is a tie-breaker
// giving the lowest vertex index of a vertex whose domain size is minimal
std::pair<int, int> bidomain_score(
        const Bidomain &bd,
        const vector<int>& left,
        const vector<vector<int>>& g0_2p, const vector<vector<int>>& g1_2p,
        const vector<int>& g0_deg,
        const vector<int>& g1_deg,
        const vector<VtxPair>& current,
        std::pair<int, int> incumbent)
{
    auto best = incumbent;
//    int best_v = INT_MAX;
//    int lowest_num_possible_assignments = INT_MAX;
    for (int i=0; i<bd.left_len; i++) {
        int v = left[bd.l + i];
        //int num_possible_assignments = 0;
        auto vtx_score = std::make_pair(0, v);
        for (int w : bd.right_set) {
//            std::cout << w << std::endl;
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
        const vector<int> & left,
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
        auto score = bidomain_score(bd, left, g0_2p, g1_2p, g0_deg, g1_deg, current, best_score);
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

// Returns length of left half of array
int partition(vector<int>& all_vv, int start, int len, const vector<unsigned int> & adjrow) {
    int i=0;
    for (int j=0; j<len; j++) {
        if (adjrow[all_vv[start+j]]) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains(const vector<Bidomain> & d,
        vector<int> & left, const Graph & g0, const Graph & g1, int v, int w)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        // After this partition, left_len is the length of the
        // array of vertices with edges from v
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);

        // right_with_edge is the set of vertices (not including w) with edges to w
        vector<int> right_with_edge;
        right_with_edge.reserve(old_bd.right_set.size());
        // right_without_w is the right set with w removed
        vector<int> right_without_w;
        right_without_w.reserve(old_bd.right_set.size());

        for (int u : old_bd.right_set) {
            if (u != w) {
                right_without_w.push_back(u);
                if (g1.adjmat[w][u]) {
                    right_with_edge.push_back(u);
                }
            }
        }

        int left_len_noedge = old_bd.left_len - left_len;
        if ((left_len_noedge > (int)right_without_w.size()) || (left_len > (int)right_with_edge.size())) {
            // Stop early if we know that there are vertices in the first graph that can't be matched
            return {};
        }
        if (left_len_noedge && right_without_w.size())
            new_d.push_back({l+left_len, left_len_noedge, std::move(right_without_w), old_bd.is_adjacent});
        if (left_len && right_with_edge.size())
            new_d.push_back({l, left_len, std::move(right_with_edge), true});
    }

    std::sort(new_d.begin(), new_d.end(),
            [](const Bidomain& a, const Bidomain& b) { return a.right_len() < b.right_len(); });

    return new_d;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(const vector<int>& arr, int start_idx, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i=0; i<len; i++) {
        if (arr[start_idx + i]>w && arr[start_idx + i]<smallest) {
            smallest = arr[start_idx + i];
            idx = i;
        }
    }
    return idx;
}

void remove_vtx_from_left_domain(vector<int>& left, Bidomain& bd, int v)
{
    int i = 0;
    while(left[bd.l + i] != v) i++;
    std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
    bd.left_len--;
}

void solve(const Graph & g0, const Graph & g1,
        const vector<int>& g0_deg, const vector<int>& g1_deg,
        const vector<vector<int>> & g0_2p, const vector<vector<int>> & g1_2p,
        vector<VtxPair> & incumbent, vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<int> & left, long long & solution_count)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(current, domains, left);
    nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        //if (!arguments.quiet) cout << "Incumbent size: " << incumbent.size() << endl;
    }

    if (current.size()==(unsigned)g0.n) {
        solution_count++;
        return;
    }

    if (!arguments.enumerate && incumbent.size()==(unsigned)g0.n)
        return;

    if ((int)current.size() + calc_bound(domains) < g0.n)
        return;

    if (!propagate_alldiff(domains, g0, g1, left))
        return;

    auto bd_idx_and_v = select_bidomain_and_branching_var(domains, left,
            g0_2p, g1_2p, g0_deg, g1_deg, current);
    int bd_idx = bd_idx_and_v.first;
    int v = bd_idx_and_v.second;

    if (bd_idx == -1)   // Return if there's nothing left to branch on
        return;
    Bidomain &bd = domains[bd_idx];

    remove_vtx_from_left_domain(left, domains[bd_idx], v);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    for (int i=0; i<bd.right_len(); i++) {
        int w = bd.right_set[i];

        if (g0_deg[v] <= g1_deg[w] && !assignment_impossible_by_2path_count(v, w, current, g0_2p, g1_2p)) {
            auto new_domains = filter_domains(domains, left, g0, g1, v, w);
            current.push_back(VtxPair(v, w));
            solve(g0, g1, g0_deg, g1_deg, g0_2p, g1_2p, incumbent, current, new_domains, left, solution_count);
            current.pop_back();
        }
    }
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
                        num_2paths[i][j]--;
                        num_2paths[j][i]--;
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

    vector<int> left;  // the big list of vertex indices for the left partitions

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
            int start_l = left.size();
            vector<int> right_set;

            for (int i=0; i<g0.n; i++)
                if (g0.label[i]==label && is_isolated==(g0_deg[i]==0))
                    left.push_back(i);
            for (int i=0; i<g1.n; i++)
                if (g1.label[i]==label && (is_isolated || g1_deg[i]>0))
                    right_set.push_back(i);

            int left_len = left.size() - start_l;

            if (left_len && right_set.size())
                domains.push_back({start_l, left_len, right_set, false});
        }
    }

    vector<VtxPair> incumbent;
    vector<VtxPair> current;
    long long solution_count = 0;
    solve(g0, g1, g0_deg, g1_deg, g0_2p, g1_2p, incumbent, current, domains, left, solution_count);

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

