#include "graph.h"

#include <algorithm>
#include <numeric>
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
#include <list>

#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE_FOR_STRONGER_BOUND 30
#define MAX_RATIO_FOR_STRONGER_BOUND 1.5

using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { min_max, min_product };

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format\vHEURISTIC can be min_max or min_product";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
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
    Heuristic heuristic;
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
                if (std::string(arg) == "min_max")
                    arguments.heuristic = min_max;
                else if (std::string(arg) == "min_product")
                    arguments.heuristic = min_product;
                else
                    fail("Unknown heuristic (try min_max or min_product)");
            } else if (arguments.arg_num == 1) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 2) {
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
    int l,        r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    bool active;
    int mod_index;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent),
            active(true),
            mod_index(INT_MAX)
            { };
    Bidomain(): Bidomain(0, 0, 0, 0, false) { };
};

using It = vector<int>::iterator;

struct NewBidomain;

using BDLL = std::list<NewBidomain>;

using BdIt = std::list<NewBidomain>::iterator;

struct NewBidomain {
    It l, r;
    It l_mid, r_mid;
    It l_end, r_end;
    bool is_adjacent;
    bool active;
    unsigned mod_index;
    BdIt reinsertion_point;
    NewBidomain(It l, It r, It l_end, It r_end, bool is_adjacent):
            l(l),
            r(r),
            l_mid(l_end),
            r_mid(r_end),
            l_end(l_end),
            r_end(r_end),
            is_adjacent (is_adjacent),
            active(true),
            mod_index(INT_MAX)
            { };
};

//// A doubly-linked list of bidomains with dummy head and tail nodes
//struct BDLL
//{
//    Bidomain head;
//    Bidomain tail;
//
//    BDLL() {
//        head.next = &tail;
//        tail.prev = &head;
//    }
//
//    Bidomain *add_bd_before(Bidomain *bd) {
//        Bidomain *new_bd = new Bidomain();
//        new_bd->next = bd;
//        new_bd->prev = bd->prev;
//        bd->prev->next = new_bd;
//        bd->prev = new_bd;
//        return new_bd;
//    }
//
//    Bidomain *add_bd_after(Bidomain *bd) {
//        Bidomain *new_bd = new Bidomain();
//        new_bd->prev = bd;
//        new_bd->next = bd->next;
//        bd->next->prev = new_bd;
//        bd->next = new_bd;
//        return new_bd;
//    }
//
//    void remove_bd(Bidomain *bd) {
//        bd->prev->next = bd->next;
//        bd->next->prev = bd->prev;
//    }
//};

struct Ptrs
{
    It vtx_it;
    BdIt bd_it;
};

void show(const vector<Ptrs> & left_ptrs, const vector<Ptrs> & right_ptrs,
        const Graph & g0, const Graph & g1, const vector<VtxPair>& current,
        const BDLL & bdll, const vector<Bidomain> &domains,
        const vector<int>& left, const vector<int>& right)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    if (!current.empty() > 0) {
        std::cout << "Adjacent to " << current.back().v << " in g0: ";
        for (int x : g0.adj_lists[current.back().v]) std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Adjacent to " << current.back().w << " in g1: ";
        for (int x : g1.adj_lists[current.back().w]) std::cout << x << " ";
        std::cout << std::endl;
    }
    cout << "---------------------" << std::endl;
    for (unsigned int i=0; i<domains.size(); i++) {
        struct Bidomain bd = domains[i];
        cout << "Left  ";
        for (int j=0; j<bd.left_len; j++)
            cout << left[bd.l + j] << " ";
        cout << std::endl;
        cout << "Right  ";
        for (int j=0; j<bd.right_len; j++)
            cout << right[bd.r + j] << " ";
        cout << std::endl;
    }
    cout << "---------------------" << std::endl;
//    for (int i=0; i<g0.n; i++) {
//        if (*left_ptrs[i].vtx_it != i) {
//            cout << "PROBLEM WITH g0 VERTEX " << i << "!!!" << endl;
//        }
//    }
//    for (int i=0; i<g1.n; i++) {
//        cout << i << " " << g1.n << std::endl;
//        cout << &*right_ptrs[i].vtx_it << std::endl;
//        if (*right_ptrs[i].vtx_it != i) {
//            cout << "PROBLEM WITH g1 VERTEX " << i << "!!!" << endl;
//        }
//    }
    for (auto & bd : bdll) {
        cout << "Left  ";
        for (It it=bd.l; it!=bd.l_end; it++)
            cout << *it << " ";
        cout << std::endl;
        cout << "Right  ";
        for (It it=bd.r; it!=bd.r_end; it++)
            cout << *it << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
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

bool can_backtrack_using_degrees_within_bidomain(const Bidomain& bd, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1)
{
    if (bd.left_len > bd.right_len)
        return true;

    // Don't bother with this procedure if it's likely to be expensive or ineffective
    if (bd.left_len == 1 || bd.left_len > MAX_SIZE_FOR_STRONGER_BOUND || bd.right_len > bd.left_len * MAX_RATIO_FOR_STRONGER_BOUND)
        return false;

    std::vector<int> left_deg(bd.left_len, 0);   // Degree in the subgraph induced by the left set of bd
    std::vector<int> right_deg(bd.right_len, 0);   // Degree in the subgraph induced by the right set of bd
    for (int i=0; i<bd.left_len; i++) {
        int v = left[bd.l + i];
        for (int j=0; j<i; j++) {
            int w = left[bd.l + j];
            if (g0.adjmat[v][w]) {
                left_deg[i]++;
                left_deg[j]++;
            }
        }
    }
    for (int i=0; i<bd.right_len; i++) {
        int v = right[bd.r + i];
        for (int j=0; j<i; j++) {
            int w = right[bd.r + j];
            if (g1.adjmat[v][w]) {
                right_deg[i]++;
                right_deg[j]++;
            }
        }
    }

    std::sort(std::begin(left_deg), std::end(left_deg));
    std::sort(std::begin(right_deg), std::end(right_deg));

    //std::cout << left_deg[0] << " " << right_deg[0] << std::endl;
    for (int i=0; i<bd.left_len; i++) {
//        std::cout << i << " " << left_deg.size() << " " << bd.left_len << std::endl;
        // Check neighbourhood degree sequence
        if (left_deg[bd.left_len-1-i] > right_deg[bd.right_len-1-i])
            return true;
        // Check neighbourhood degree sequence in complement graph
        if (bd.left_len-1-left_deg[i] > bd.right_len-1-right_deg[i])
            return true;
    }

    return false;
}

int calc_bound(const vector<Bidomain>& domains, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int target)
{
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
#ifdef TIGHTER_BOUNDING
    if (bound < target)
        return 0;   // bactrack
    for (const Bidomain &bd : domains) {
        if (can_backtrack_using_degrees_within_bidomain(bd, left, right, g0, g1)) {
            return 0;  // backtrack
        }
    }
#endif
    return bound;
}

int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

int select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,
        int current_matching_size)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        int len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
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

std::pair<vector<BdIt>, BDLL> new_filter_domains(
        BDLL & bdll, vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        bool multiway)
{
    // TODO: quit early if solution is impossible?
    vector<BdIt> split_bds;
    BDLL deleted_bds;
    for (int u : g0.adj_lists[v]) {
        auto bd_it = left_ptrs[u].bd_it;
        auto & bd = *bd_it;
        if (!bd.active) continue;
        It u_it = left_ptrs[u].vtx_it;
        if (u_it >= bd.l_end) continue;
        if (bd.mod_index >= split_bds.size() || split_bds[bd.mod_index] != bd_it) {
            bd.l_mid = bd.l_end;
            bd.r_mid = bd.r_end;
            bd.mod_index = split_bds.size();
            split_bds.push_back(bd_it);
        }
        It m = std::prev(bd.l_mid);
        std::swap(*u_it, *m);
        left_ptrs[u].vtx_it = m;
        left_ptrs[*u_it].vtx_it = u_it;
        --bd.l_mid;
    }
    for (int u : g1.adj_lists[w]) {
        auto bd_it = right_ptrs[u].bd_it;
        auto & bd = *bd_it;
        if (!bd.active) continue;
        It u_it = right_ptrs[u].vtx_it;
        if (u_it >= bd.r_end) continue;
        if (bd.mod_index >= split_bds.size() || split_bds[bd.mod_index] != bd_it) {
            bd.l_mid = bd.l_end;
            bd.r_mid = bd.r_end;
            bd.mod_index = split_bds.size();
            split_bds.push_back(bd_it);
        }
        It m = std::prev(bd.r_mid);
        std::swap(*u_it, *m);
        right_ptrs[u].vtx_it = m;
        right_ptrs[*u_it].vtx_it = u_it;
        --bd.r_mid;
    }
    for (auto bd_it : split_bds) {
        // TODO: fix connectedness
        bdll.emplace(std::next(bd_it), bd_it->l_mid, bd_it->r_mid, bd_it->l_end, bd_it->r_end, false);
        for (It it=bd_it->l_mid; it<bd_it->l_end; it++) {
            left_ptrs[*it].bd_it = std::next(bd_it);
        }
        for (It it=bd_it->r_mid; it<bd_it->r_end; it++) {
            right_ptrs[*it].bd_it = std::next(bd_it);
        }
        bd_it->l_end = bd_it->l_mid;
        bd_it->r_end = bd_it->r_mid;
        BdIt bd_its[2] = {bd_it, std::next(bd_it)};
        for (int i=0; i<2; i++) {
            BdIt bd_it = bd_its[i];
            // Delete old and new BDs if necessary
            auto & bd = *bd_it;
            if (bd.l == bd.l_end || bd.r == bd.r_end) {
                bd.active = false;
                bd.reinsertion_point = std::next(bd_it);
                deleted_bds.splice(deleted_bds.end(), bdll, bd_it);
            }
        }
    }
    return {std::move(split_bds), std::move(deleted_bds)};
}

void unfilter_domains(
        BDLL & bdll,
        vector<Ptrs> & left_ptrs,
        vector<Ptrs> & right_ptrs,
        vector<BdIt> & split_bds,
        BDLL & deleted_bds,
        const Graph & g0,
        const Graph & g1)
{
    while (!deleted_bds.empty()) {
        NewBidomain & bd = deleted_bds.back();
        bd.active = true;
        bdll.splice(bd.reinsertion_point, deleted_bds, std::prev(deleted_bds.end()));
    }

    while(!split_bds.empty()) {
        BdIt bd_it = split_bds.back();
        auto & nxt = *std::next(bd_it);
        for (It it=nxt.l; it<nxt.l_end; it++) {
            left_ptrs[*it].bd_it = bd_it;
        }
        for (It it=nxt.r; it<nxt.r_end; it++) {
            right_ptrs[*it].bd_it = bd_it;
        }
        bd_it->l_end = nxt.l_end;
        bd_it->r_end = nxt.r_end;
        bdll.erase(std::next(bd_it));
        split_bds.pop_back();
    }
}

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains(BDLL & bdll, const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        bool multiway)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
        int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if ((left_len_noedge > right_len_noedge) || (left_len > right_len)) {
            // Stop early if we know that there are vertices in the first graph that can't be matched
            // TODO: improve this for the edge-labelled case
            //return new_d;
        }
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            auto& adjrow_v = g0.adjmat[v];
            auto& adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin+left_len, [&](int a, int b)
                    { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin+right_len, [&](int a, int b)
                    { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l<l_top && r<r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
                    do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
                    new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }
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

void remove_bidomain(vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

void assign(int v, int w,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs)
{
    auto & bd = *left_ptrs[v].bd_it;
    It v_it = left_ptrs[v].vtx_it;
    It l_last = std::prev(bd.l_end);
    std::swap(*v_it, *l_last);
    left_ptrs[v].vtx_it = l_last;
    left_ptrs[*v_it].vtx_it = v_it;
    --bd.l_end;

    It w_it = right_ptrs[w].vtx_it;
    It r_last = std::prev(bd.r_end);
    std::swap(*w_it, *r_last);
    right_ptrs[w].vtx_it = r_last;
    right_ptrs[*w_it].vtx_it = w_it;
    --bd.r_end;
}

void unassign(int v,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs)
{
    auto & bd = *left_ptrs[v].bd_it;
//    cout << (bd.l_end - bd.l) << "------------" << endl;
    ++bd.l_end;
    ++bd.r_end;
}

void solve(const Graph & g0, const Graph & g1, vector<VtxPair> & incumbent,
        vector<VtxPair> & current,
        BDLL & bdll, vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, long long & solution_count)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll, domains, left, right);
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

    int bound = current.size() + calc_bound(
            domains, left, right, g0, g1, g0.n - current.size());
    if (bound < g0.n)
        return;

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1)   // Return if there's nothing left to branch on
        return;
    Bidomain &bd = domains[bd_idx];

    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    int w = -1;
    bd.right_len--;
    for (int i=0; i<=bd.right_len; i++) {
        int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
        w = right[bd.r + idx];

        // swap w to the end of its colour class
        right[bd.r + idx] = right[bd.r + bd.right_len];
        right[bd.r + bd.right_len] = w;

////        cout << "BEFORE ASSIGNMENT" << std::endl;
////        if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll, domains, left, right);
        assign(v, w, left_ptrs, right_ptrs);

        BDLL removed_bd_lst;
        if (bd.left_len == 0) {
            BdIt bd_it = left_ptrs[v].bd_it;
            bd_it->active = false;
            bd_it->reinsertion_point = std::next(bd_it);
            removed_bd_lst.splice(removed_bd_lst.end(), bdll, bd_it);
        }

////        cout << "BEFORE FILTERING" << std::endl;
////        if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll, domains, left, right);

        vector<BdIt> split_bds;
        BDLL deleted_bds;
        tie(split_bds, deleted_bds) = new_filter_domains(bdll, left_ptrs, right_ptrs, domains, left, right, g0, g1, v, w,
                arguments.directed || arguments.edge_labelled);

        auto new_domains = filter_domains(bdll, domains, left, right, g0, g1, v, w,
                arguments.directed || arguments.edge_labelled);
        current.push_back(VtxPair(v, w));
////        cout << "AFTER FILTERING" << std::endl;
////        if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll, new_domains, left, right);
        solve(g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, new_domains, left, right, solution_count);
        current.pop_back();

        unfilter_domains(bdll, left_ptrs, right_ptrs, split_bds, deleted_bds, g0, g1);

        if (!removed_bd_lst.empty()) {
            removed_bd_lst.front().active = true;
            bdll.splice(removed_bd_lst.front().reinsertion_point, removed_bd_lst, removed_bd_lst.begin());
        }
        unassign(v, left_ptrs, right_ptrs);

////        cout << "AFTER unassign()" << std::endl;
////        if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll, domains, left, right);

        if (!arguments.enumerate && incumbent.size()==(unsigned)g0.n)
            break;
    }
}

// Returns a common subgraph and the number of induced subgraph isomorphisms found
std::pair<vector<VtxPair>, long long> mcs(Graph & g0, Graph & g1)
{
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions
    vector<int> new_left;
    vector<int> new_right;
    new_left.reserve(g0.n);
    new_right.reserve(g1.n);

    vector<Ptrs> left_ptrs(g0.n);
    vector<Ptrs> right_ptrs(g1.n);

    vector<bool> g0_active_vertices(g0.n);
    vector<bool> g1_active_vertices(g1.n);

    auto domains = vector<Bidomain> {};
    BDLL bdll;

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

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (int i=0; i<g0.n; i++) {
            if (g0.label[i]==label) {
                new_left.push_back(i);
                left.push_back(i);
                left_ptrs[i].vtx_it = std::prev(new_left.end());
                g0_active_vertices[i] = true;
            }
        }
        for (int i=0; i<g1.n; i++) {
            if (g1.label[i]==label) {
                new_right.push_back(i);
                right.push_back(i);
                right_ptrs[i].vtx_it = std::prev(new_right.end());
                //cout << *right_ptrs[i].vtx_it << endl;
                g1_active_vertices[i] = true;
            }
        }

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
        bdll.emplace_back(new_left.begin() + start_l,
                          new_right.begin() + start_r,
                          new_left.begin() + start_l + left_len,
                          new_right.begin() + start_r + right_len,
                          false);

        for (int i=0; i<g0.n; i++) {
            if (g0.label[i]==label) {
                left_ptrs[i].bd_it = std::prev(bdll.end());
            }
        }
        for (int i=0; i<g1.n; i++) {
            if (g1.label[i]==label) {
                right_ptrs[i].bd_it = std::prev(bdll.end());
            }
        }
    }

    make_adj_lists(g0, g0_active_vertices);
    make_adj_lists(g1, g1_active_vertices);

    vector<VtxPair> incumbent;
    vector<VtxPair> current;
    long long solution_count = 0;
    solve(g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, domains, left, right, solution_count);

    return {incumbent, solution_count};
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++) {
        for (int w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
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

    // TODO: As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the IJCAI 2017 paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) > g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) > g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
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

