
#include "sparse_graph.h"

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

using std::vector;
using std::cout;
using std::endl;

/******************* BITSETS *******************/

void set_first_n_bits(vector<unsigned long long> & bitset, int n)
{
    int i = 0;
    while (n >= 64) {
        bitset[i] = ~0ull;
        n -= 64;
        ++i;
    }
    if (n > 0)
        bitset[i] = (1ull << n) - 1;
}

bool test_bit(const vector<unsigned long long> & bitset, int bit)
{
    return 0 != (bitset[bit/BITS_PER_WORD] & (1ull << (bit%BITS_PER_WORD)));
}

void set_bit(vector<unsigned long long> & bitset, int bit)
{
    bitset[bit/BITS_PER_WORD] |= (1ull << (bit%BITS_PER_WORD));
}

void unset_bit(vector<unsigned long long> & bitset, int bit)
{
    bitset[bit/BITS_PER_WORD] &= ~(1ull << (bit%BITS_PER_WORD));
}

void unset_bit_if(vector<unsigned long long> & bitset, int bit, bool condition)
{
    bitset[bit/BITS_PER_WORD] &= ~((unsigned long long) condition << (bit%BITS_PER_WORD));
}

int bitset_popcount(const vector<unsigned long long> & bitset)
{
    int count = 0;
    for (unsigned i=0; i<bitset.size(); i++)
        count += __builtin_popcountll(bitset[i]);
    return count;
}

int bitset_popcount_of_intersection(
        const vector<unsigned long long> & bitset1,
        const vector<unsigned long long> & bitset2)
{
    int count = 0;
    for (unsigned i=0; i<bitset1.size(); i++)
        count += __builtin_popcountll(bitset1[i] & bitset2[i]);
    return count;
}

int bitset_popcount_of_difference(
        const vector<unsigned long long> & bitset1,
        const vector<unsigned long long> & bitset2)
{
    int count = 0;
    for (unsigned i=0; i<bitset1.size(); i++)
        count += __builtin_popcountll(bitset1[i] & ~bitset2[i]);
    return count;
}

bool bitset_empty(const vector<unsigned long long> & bitset)
{
    for (unsigned i=0; i<bitset.size(); i++)
        if (bitset[i] != 0)
            return false;
    return true;
}

int first_set_bit(const vector<unsigned long long> & bitset)
{
    for (unsigned i=0; i<bitset.size(); i++)
        if (bitset[i] != 0)
            return i*BITS_PER_WORD + __builtin_ctzll(bitset[i]);
    return -1;
}

int last_set_bit(const vector<unsigned long long> & bitset)
{
    for (unsigned i=bitset.size(); i--; )
        if (bitset[i] != 0)
            return i*BITS_PER_WORD + BITS_PER_WORD - 1 - __builtin_clzll(bitset[i]);
    return -1;
}

void bitset_intersection(const vector<unsigned long long> & src1,
                                     const vector<unsigned long long> & src2,
                                     vector<unsigned long long> & dst)
{
    for (unsigned i=0; i<src1.size(); i++)
        dst[i] = src1[i] & src2[i];
}

void bitset_intersect_with(vector<unsigned long long> & bitset1,
                                     const vector<unsigned long long> & bitset2)
{
    for (unsigned i=0; i<bitset1.size(); i++)
        bitset1[i] &= bitset2[i];
}

void bitset_add_all(vector<unsigned long long> & bitset1,
                                     const vector<unsigned long long> & bitset2)
{
    for (unsigned i=0; i<bitset1.size(); i++)
        bitset1[i] |= bitset2[i];
}

void bitset_remove_all(vector<unsigned long long> & bitset,
                                     const vector<unsigned long long> & bitset2)
{
    for (unsigned i=0; i<bitset.size(); i++)
        bitset[i] &= ~bitset2[i];
}

void bitset_intersection_with_complement(const vector<unsigned long long> & src1,
                                     const vector<unsigned long long> & src2,
                                     vector<unsigned long long> & dst)
{
    for (unsigned i=0; i<src1.size(); i++)
        dst[i] = src1[i] & ~src2[i];
}

void copy_bitset(const vector<unsigned long long> & src,
                        vector<unsigned long long> & dest)
{
    for (unsigned i=0; i<src.size(); i++)
        dest[i] = src[i];
}

void clear_bitset(vector<unsigned long long> & bitset)
{
    for (unsigned i=0; i<bitset.size(); i++)
        bitset[i] = 0ull;
}

template<typename F>
void bitset_foreach(const vector<unsigned long long> & bitset, F f)
{
        for (unsigned i=0; i<bitset.size(); i++) {
            unsigned long long word = bitset[i];
            while (word) {
                int bit = __builtin_ctzll(word);
                word ^= (1ull << bit);
                int v = i*BITS_PER_WORD + bit;
                f(v);
            }
        }
}

/******************* END OF BITSETS *******************/

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { heur_A, heur_B };

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Subgraph isomorphism\vHEURISTIC can be A or B";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"enumerate", 'e', 0, 0, "Count solutions"},
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
    arguments.directed = false;   // TODO: remove (unused)
    arguments.enumerate = false;
    arguments.edge_labelled = false;  // TODO: remove (unused)
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
        case 'e':
            arguments.enumerate = true;
            break;
        case 'x':
            arguments.vertex_labelled = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "A")
                    arguments.heuristic = heur_A;
                else if (std::string(arg) == "B")
                    arguments.heuristic = heur_B;
                else
                    fail("Unknown heuristic (try A or B)");
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

//using BDLL = std::list<NewBidomain>;

using BdIt = NewBidomain *;  //std::list<NewBidomain>::iterator;

struct NewBidomain {
    It l;
    It l_mid;
    It l_end;
    vector<unsigned long long> right_set;
    union {
        NewBidomain *next_in_split_list;
        NewBidomain *next_in_free_list;
    };
    NewBidomain *next_in_deleted_list;
    NewBidomain *prev;
    NewBidomain *next;
    bool active;
    bool undergoing_split;

    void initialise(It l, It r, It l_end, It r_end)
    {
        this->l = l;
        this->l_end = l_end;
        this->active = true;
        this->undergoing_split = false;
        clear_bitset(right_set);
        for (auto it=r; it!=r_end; it++)
            set_bit(right_set, *it);
    }

    void initialise(It l, It l_end,
            vector<unsigned long long> & r_set1,
            vector<unsigned long long> & r_set2)
    {
        this->l = l;
        this->l_end = l_end;
        this->active = true;
        this->undergoing_split = false;
        bitset_intersection(r_set1, r_set2, right_set);
    }

    void insert_before(BdIt p)
    {
        BdIt new_prev = p->prev;
        this->prev = new_prev;
        this->next = p;
        p->prev = this;
        new_prev->next = this;
    }
    void insert_after(BdIt p)
    {
        BdIt new_next = p->next;
        this->next = new_next;
        this->prev = p;
        p->next = this;
        new_next->prev = this;
    }
    void move_to_before(BdIt p)
    {
        this->prev->next = this->next;
        this->next->prev = this->prev;
        this->insert_before(p);
    }
    void move_to_after(BdIt p)
    {
        this->prev->next = this->next;
        this->next->prev = this->prev;
        this->insert_after(p);
    }
    void remove()
    {
        this->prev->next = this->next;
        this->next->prev = this->prev;
    }
    void reinsert()
    {
        this->prev->next = this;
        this->next->prev = this;
    }
    int l_size() const { return l_end - l; }
    int r_size() const { return bitset_popcount(right_set); }
};

struct BDLL {
    // A circular doubly-linked list
    // (It's just circular to avoid the need for separate head and tail)
    NewBidomain head;
    BDLL() {
        this->head.next = &this->head;
        this->head.prev = &this->head;
    }

    NewBidomain *begin() { return head.next; }
    NewBidomain *end() { return &head; }
    NewBidomain & back() { return *head.prev; }
    int size()
    {
        int result = 0;
        for (BdIt bd_it=begin(); bd_it!=end(); bd_it=bd_it->next) {
            ++result;
        }
        return result;
    }
    bool empty() { return head.next == &head; }
};


// Some temporary storage space
struct Workspace {
    vector<BdIt> split_bds;
    struct Graph & g1;
    vector<vector<unsigned long long>> g1_adj_bitsets;
    NewBidomain *get_from_free_list()
    {
        if (bd_free_list == nullptr) {
            bd_memory_pools.push_back(vector<NewBidomain>(100));
            int bitset_sz = g1.bitset_word_count();
            for (NewBidomain & bd : bd_memory_pools.back()) {
                bd.next_in_free_list = bd_free_list;
                bd_free_list = &bd;
                bd.right_set.resize(bitset_sz);
            }
        }
        NewBidomain *bd = bd_free_list;
        bd_free_list = bd->next_in_free_list;
        return bd;
    }
    void add_to_free_list(NewBidomain * bd)
    {
        bd->next_in_free_list = bd_free_list;
        bd_free_list = bd;
    }
    Workspace(struct Graph & g1)
            : g1(g1), g1_adj_bitsets(g1.n, vector<unsigned long long>(g1.bitset_word_count()))
    {
        for (int v=0; v<g1.n; v++) {
            for (int w : g1.filtered_adj_lists[v]) {
                set_bit(g1_adj_bitsets[v], w);
            }
        }
    }
private:
    NewBidomain *bd_free_list = nullptr;
    std::list<vector<NewBidomain>> bd_memory_pools;
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
    BdIt bd_it;
    It vtx_it;
};

void show(const vector<Ptrs> & left_ptrs, const vector<Ptrs> & right_ptrs,
        const Graph & g0, const Graph & g1, const vector<VtxPair>& current,
        BDLL & bdll)
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
        for (int x : g0.filtered_adj_lists[current.back().v]) std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Adjacent to " << current.back().w << " in g1: ";
        for (int x : g1.filtered_adj_lists[current.back().w]) std::cout << x << " ";
        std::cout << std::endl;
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
    for (BdIt bd_it=bdll.begin(); bd_it!=bdll.end(); bd_it=bd_it->next) {
        cout << "Left  ";
        for (It it=bd_it->l; it!=bd_it->l_end; it++)
            cout << *it << " ";
        cout << std::endl;
        cout << "Right_ ";
        bitset_foreach(bd_it->right_set, [](int w) { cout << w << " "; });
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
            if (g0.has_edge(p0.v, p1.v) != g1.has_edge(p0.w, p1.w))
                return false;
        }
    }
    return true;
}

int calc_bound(Workspace & workspace, BDLL & bdll, const Graph & g0, const Graph & g1, int target)
{
    int bound = 0;
    for (BdIt bd_it=bdll.begin(); bd_it!=bdll.end(); bd_it=bd_it->next) {
        bound += std::min(bd_it->l_size(), bd_it->r_size());
    }
    return bound;
}

BdIt select_bidomain_heur_A(BDLL & domains)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int right_len = bd.r_size();
        if (right_len < min_size)
            min_size = right_len;
    }
    int min_tie_breaker = INT_MAX;
    BdIt best = domains.end();
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int right_len = bd.r_size();
        if (right_len != min_size)
            continue;
        int tie_breaker = *std::min_element(bd.l, bd.l_end);
        if (tie_breaker < min_tie_breaker) {
            min_tie_breaker = tie_breaker;
            best = bd_it;
        }
    }
    return best;
}

BdIt select_bidomain_heur_B(BDLL & domains, const Graph & g0)
{
    double best_score = INT_MIN;
    BdIt best = domains.end();
    for (BdIt bd_it=domains.begin(); bd_it!=domains.end(); bd_it=bd_it->next) {
        auto const & bd = *bd_it;
        int right_len = bd.r_size();
        if (right_len == 1) {
            // Special case where no branching is required
            return bd_it;
        }
        int deg = g0.adj_lists[*std::min_element(bd.l, bd.l_end)].size();
        double score = double(deg) / right_len;
        if (score > best_score) {
            best_score = score;
            best = bd_it;
        }
    }
    return best;
}

BdIt select_bidomain(BDLL & domains, const Graph & g0)
{
    if (arguments.heuristic == heur_A)
        return select_bidomain_heur_A(domains);
    else
        return select_bidomain_heur_B(domains, g0);
}

struct SplitAndDeletedLists {
    bool quit_early;
    NewBidomain *split_bds_list;
    NewBidomain *deleted_bds_list;
    SplitAndDeletedLists(bool quit_early, NewBidomain *split_bds_list, NewBidomain *deleted_bds_list)
        : quit_early(quit_early), split_bds_list(split_bds_list), deleted_bds_list(deleted_bds_list) {}
};

void partition_left(const vector<int> & left_vv, vector<Ptrs> & left_ptrs,
        vector<BdIt> & split_bds)
{
    for (int u : left_vv) {
        auto & u_ptrs = left_ptrs[u];
        auto bd_it = u_ptrs.bd_it;
        if (bd_it == nullptr) continue;
        auto & bd = *bd_it;
        if (!bd.active) continue;
        if (!bd.undergoing_split) {
            bd.l_mid = bd.l_end;
            bd.undergoing_split = true;
            split_bds.push_back(bd_it);
        }
        It u_it = u_ptrs.vtx_it;
        It m = std::prev(bd.l_mid);
        std::swap(*u_it, *m);
        u_ptrs.vtx_it = m;
        left_ptrs[*u_it].vtx_it = u_it;
        --bd.l_mid;
    }
}

// TODO: give this function a better name. It doesn't partition; it just
// checks whether we can backtrack early.
bool partition_right(Workspace & workspace, const vector<int> & right_vv, vector<Ptrs> & right_ptrs,
        vector<BdIt> & split_bds, BDLL & bdll, int w)
{
    for (BdIt bd_it : split_bds) {
        int l1_size = bd_it->l_end - bd_it->l_mid;
        int r1_size = bitset_popcount_of_intersection(bd_it->right_set, workspace.g1_adj_bitsets[w]);
        if (l1_size > r1_size) {
            return false;
        }
    }
    return true;
}

NewBidomain * do_splits(Workspace & workspace, vector<BdIt> & split_bds,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs, int w)
{
    NewBidomain *split_bds_list = nullptr;
    for (auto bd_it : split_bds) {
        BdIt new_elem = workspace.get_from_free_list();
        new_elem->insert_after(bd_it);
        new_elem->initialise(bd_it->l_mid, bd_it->l_end,
                bd_it->right_set, workspace.g1_adj_bitsets[w]);

        // Insert the new BD at the head of the linked list of split BDs
        new_elem->next_in_split_list = split_bds_list;
        split_bds_list = new_elem;

        for (It it=bd_it->l_mid; it<bd_it->l_end; it++) {
            left_ptrs[*it].bd_it = new_elem;
        }
        bd_it->l_end = bd_it->l_mid;
    }
    return split_bds_list;
}

NewBidomain * do_deletions(vector<BdIt> & split_bds)
{
    NewBidomain *deleted_bds_list = nullptr;
    for (auto bd_it : split_bds) {
        for (int i=0; i<2; i++) {
            // Delete old and new BDs if necessary
            auto & bd = *bd_it;
            if (bd.l == bd.l_end) {
                bd.active = false;

                // add to deleted list
                bd.remove();
                bd.next_in_deleted_list = deleted_bds_list;
                deleted_bds_list = bd_it;
            }
            bd_it = bd_it->next;
        }
    }
    return deleted_bds_list;
}

SplitAndDeletedLists filter_domains(
        Workspace & workspace,
        BDLL & bdll, vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        const Graph & g0, const Graph & g1, int v, int w)
{
    vector<BdIt> & split_bds = workspace.split_bds;
    split_bds.clear();

    partition_left(g0.filtered_adj_lists[v], left_ptrs, split_bds);

    bool partition_right_success = partition_right(
            workspace, g1.filtered_adj_lists[w], right_ptrs, split_bds, bdll, w);

    for (auto bd_it : split_bds) {
        bd_it->undergoing_split = false;
    }

    if (!partition_right_success) {
        // A solution extending the current partial solution is impossible,
        // so quit this function early.
        return { true, nullptr, nullptr };
    }

    NewBidomain *split_bds_list = do_splits(
            workspace, split_bds, left_ptrs, right_ptrs, w);

    NewBidomain *deleted_bds_list = do_deletions(split_bds);

    return {false, split_bds_list, deleted_bds_list};
}

void unfilter_domains(
        Workspace & workspace,
        BDLL & bdll,
        vector<Ptrs> & left_ptrs,
        vector<Ptrs> & right_ptrs,
        NewBidomain *split_bds_list,
        NewBidomain *deleted_bds,
        const Graph & g0,
        const Graph & g1)
{
    while (deleted_bds != nullptr) {
        NewBidomain & bd = *deleted_bds;
        bd.active = true;
        bd.reinsert();
        deleted_bds = bd.next_in_deleted_list;
        //deleted_bds.back().move_to_before(bd.reinsertion_point);
        //bdll.splice(bd.reinsertion_point, deleted_bds, std::prev(deleted_bds.end()));
    }

    for (NewBidomain *p=split_bds_list; p!=nullptr; ) {
        // TODO: better variable names
        //BdIt nxt_it = p;
        auto & nxt = *p;
        BdIt bd_it = nxt.prev;
        for (It it=nxt.l; it<nxt.l_end; it++) {
            left_ptrs[*it].bd_it = bd_it;
        }
        bd_it->l_end = nxt.l_end;

        p->remove();
        p = p->next_in_split_list;

        workspace.add_to_free_list(&nxt);
        //workspace.bd_free_list.splice(workspace.bd_free_list.end(), bdll, nxt_it);
        //bdll.erase(std::next(bd_it));
    }
}

// If any right-set becomes too large, return false
bool assign(int v, int w,
        BdIt bd_it,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        BDLL & bdll,
        vector<BdIt> & bds_with_w_removed /* filled in by this function */)
{
    auto & bd = *left_ptrs[v].bd_it;
    It v_it = left_ptrs[v].vtx_it;
    It l_last = std::prev(bd.l_end);
    std::swap(*v_it, *l_last);
    left_ptrs[v].vtx_it = l_last;
    left_ptrs[v].bd_it = nullptr;
    left_ptrs[*v_it].vtx_it = v_it;
    --bd.l_end;

    for (BdIt bd_it=bdll.begin(); bd_it!=bdll.end(); bd_it=bd_it->next) {
        if (test_bit(bd_it->right_set, w)) {
            bds_with_w_removed.push_back(bd_it);
            unset_bit(bd_it->right_set, w);
            if (bitset_popcount(bd_it->right_set) < bd_it->l_size()) {
                return false;
            }
        }
    }
    return true;
}

void unassign(int v, int w, BdIt bd_it,
        vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        vector<BdIt> & bds_with_w_removed)
{
    left_ptrs[v].bd_it = bd_it;
    right_ptrs[w].bd_it = bd_it;
//    cout << (bd.l_end - bd.l) << "------------" << endl;
    ++bd_it->l_end;

    for (BdIt bd_it : bds_with_w_removed) {
        set_bit(bd_it->right_set, w);
    }
}

void solve(Workspace & workspace, const Graph & g0, const Graph & g1, vector<VtxPair> & incumbent,
        vector<VtxPair> & current,
        BDLL & bdll, vector<Ptrs> & left_ptrs, vector<Ptrs> & right_ptrs,
        long long & solution_count)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(left_ptrs, right_ptrs, g0, g1, current, bdll);
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

    // TODO: don't calculate bound explicitly? Just ensure rightlen >= leftlen for each BD
    int bound = current.size() + calc_bound(workspace, bdll, g0, g1, g0.n - current.size());
    if (bound < g0.n)
        return;

    BdIt bd_it = select_bidomain(bdll, g0);
        
    if (bd_it == bdll.end())
        return;

    vector<unsigned long long> ww_bitset(bd_it->right_set);
    std::vector<int> ww;
    ww.reserve(bitset_popcount(ww_bitset));
    bitset_foreach(ww_bitset, [&ww](int w){ ww.push_back(w); });

    int v = *std::min_element(bd_it->l, bd_it->l_end);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    for (int w : ww) {
        // Not necessary, but sometimes helps:
        if (g0.adj_lists[v].size() > g1.adj_lists[w].size())
            continue;

        vector<BdIt> bds_with_w_removed {};
        if (assign(v, w, bd_it, left_ptrs, right_ptrs, bdll, bds_with_w_removed)) {
            bool removed_bd = bd_it->l_size() == 0;
            if (removed_bd) {
                bd_it->active = false;
                bd_it->remove();
            }
            //BDLL removed_bd_lst;
            //if (bd_it->l_size() == 0) {
            //    bd_it->active = false;
            //    bd_it->reinsertion_point = std::next(bd_it);
            //    bd_it->move_to_before(removed_bd_lst.end());
            //    //removed_bd_lst.splice(removed_bd_lst.end(), bdll, bd_it);
            //}

            auto filter_result = filter_domains(workspace,
                    bdll, left_ptrs, right_ptrs, g0, g1, v, w);

            if (!filter_result.quit_early) {
                current.push_back(VtxPair(v, w));

                solve(workspace, g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, solution_count);

                current.pop_back();

                unfilter_domains(workspace, bdll, left_ptrs, right_ptrs,
                        filter_result.split_bds_list, filter_result.deleted_bds_list, g0, g1);
            }

            if (removed_bd) {
                bd_it->active = true;
                bd_it->reinsert();
            }
        }
        unassign(v, w, bd_it, left_ptrs, right_ptrs, bds_with_w_removed);

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

    BDLL bdll;
    Workspace workspace {g1};

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

        NewBidomain *new_elem = workspace.get_from_free_list();
        new_elem->insert_before(&bdll.head);
        new_elem->initialise(new_left.begin() + start_l,
                          new_right.begin() + start_r,
                          new_left.begin() + start_l + left_len,
                          new_right.begin() + start_r + right_len);

        for (int i=0; i<g0.n; i++) {
            if (g0.label[i]==label) {
                left_ptrs[i].bd_it = &bdll.back();
            }
        }
        for (int i=0; i<g1.n; i++) {
            if (g1.label[i]==label) {
                right_ptrs[i].bd_it = &bdll.back();
            }
        }
    }

    filter_adj_lists(g0, g0_active_vertices);
    filter_adj_lists(g1, g1_active_vertices);

    vector<VtxPair> incumbent;
    vector<VtxPair> current;
    long long solution_count = 0;
    solve(workspace, g0, g1, incumbent, current, bdll, left_ptrs, right_ptrs, solution_count);

    return {incumbent, solution_count};
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree;
    degree.reserve(g.n);
    for (int v=0; v<g.n; v++) {
        degree.push_back(g.adj_lists[v].size());
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

    // TODO: improve the density-based decision?
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) * 2 > g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) * 2 > g0.n*(g0.n-1);
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

