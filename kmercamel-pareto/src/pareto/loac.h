#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <random>
#include <iomanip>
#include <string>
#include <stdexcept>

#include "../kmers.h"
#include "objective.h"
#include "index.h"

size_t RANDOM_SEED = 0;
size_t MAX_COUNT_WIDTH = 12;
size_t MAX_ITERS_WIDTH = 3;

/// Heuristic search for the pareto optimization on the Aho-Corasick automaton of the k-mers, represented only by sorted leaves.
/// A node is given by the index of its first leaf and its depth.
/// priority_t has to hold every budget used (see max_budget), otherwise searches would repeat expansions of leaves.
template <typename kmer_t, typename size_n_max, typename priority_t = uint8_t>
class LeafOnlyAC {
    static constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();
    static constexpr const size_t NO_PARENT = std::numeric_limits<size_t>::max();
    static constexpr const size_k_max PRIORITY_SATURATION = std::numeric_limits<priority_t>::max();

    struct SearchEntry {
        size_n_max node;
        size_k_max priority;                        /// Remaining budget
        size_k_max depth;
        size_t log;                                 /// Index of the last visited leaf in the search log
    };
    /// Leaf expanded by a search, parents give the path back to the origin
    struct LogEntry {
        size_n_max leaf;
        priority_t old_priority;                    /// Remaining priority before this search
        priority_t new_priority;                    /// Remaining priority written by this search
        bool restored;
        size_t parent;
    };
    /// State of a search
    struct SearchContext {
        std::vector<SearchEntry> stack;
        std::vector<LogEntry> log;

        size_n_max target;                                   /// Result: leaf to connect to
        std::vector<size_n_max> path;                        /// and reused leaves on the way
    };

    const size_k_max        K;
    const size_n_max        N;                      /// Number of leaves (k-mers)
    const bool              COMPLEMENTS;
    const ParetoObjective   OBJECTIVE;
    const size_k_max        PENALTY;
    const size_k_max        MAX_BUDGET;             /// Every connection not passing through the root costs less

    const std::vector<kmer_t>&       kMers;
    FailureIndex<kmer_t, size_n_max> failureIndex;

    std::vector<size_n_max> next;                   /// Next leaf in a chain
    std::vector<size_n_max> other_end;              /// The other end of a chain (valid for chain ends)
    std::vector<size_n_max> next_unused;            /// next_unused[i] == i iff leaf i is unused (= starts a chain)
    std::vector<priority_t> remaining_priorities;   /// Maximal remaining budget a leaf was expanded with
    std::vector<size_n_max> skip_to;                /// Skip pointers over leaves with high remaining priorities

    std::vector<size_n_max> backtracks;             /// Blocks [count, leaf_1, ..., leaf_count] of reused leaves
    std::vector<std::pair<size_n_max, size_t>> backtrack_positions; /// (leaf, position of its block)
    std::vector<uint64_t> has_backtrack;            /// Bitvector over leaves with rank, built before printing
    std::vector<size_t>   has_backtrack_rank;

    SearchContext context;

    // Chains

    inline bool is_used(size_n_max leaf) const { return next_unused[leaf] != leaf; }
    inline void mark_used(size_n_max leaf){ next_unused[leaf] = leaf + 1; }
    /// The first unused leaf at index >= leaf (N if none), with path halving
    inline size_n_max find_unused(size_n_max leaf){
        while (next_unused[leaf] != leaf){
            next_unused[leaf] = next_unused[next_unused[leaf]];
            leaf = next_unused[leaf];
        }
        return leaf;
    }
    inline kmer_t complement_kmer(size_n_max leaf) const { return ReverseComplement(kMers[leaf], K); }
    /// Index of the reverse complement, self-complementary k-mers are present twice and complement each other
    inline size_n_max complement_of(size_n_max leaf) const {
        kmer_t complement = complement_kmer(leaf);
        size_n_max index = failureIndex.find_kmer(complement);
        if (index == INVALID_NODE) throw std::logic_error("Reverse complement of a k-mer is missing");
        if (complement == kMers[leaf] && index == leaf) ++index;
        return index;
    }
    /// Whether the unused leaf can follow the chain ending with origin (not the own chain start nor its complement)
    inline bool can_connect(size_n_max origin, kmer_t origin_complement, size_n_max unused) const {
        if (unused == other_end[origin]) return false;
        if (COMPLEMENTS && kMers[unused] == origin_complement) return false;
        return true;
    }
    inline size_t find_backtrack(size_n_max leaf) const {
        if (!((has_backtrack[leaf / 64] >> (leaf % 64)) & 1)) return std::numeric_limits<size_t>::max();
        size_t rank = has_backtrack_rank[leaf / 64] + __builtin_popcountll(has_backtrack[leaf / 64] & ((uint64_t(1) << (leaf % 64)) - 1));
        return backtrack_positions[rank].second;
    }
    void connect_leaves(size_n_max uncompleted_leaf, size_n_max unused_leaf, const std::vector<size_n_max> &path);
    void store_backtrack(size_n_max leaf, const std::vector<size_n_max> &path, bool complemented);
    void build_backtrack_lookup();

    // Searching

    bool try_complete_leaf_phase_1(size_n_max leaf_to_complete);
    bool search(size_n_max origin, size_k_max budget);
    void push_failure(size_k_max priority, size_k_max node_depth, size_n_max leaf, size_t log);
    void restore_unexpanded();
    void phase_2(std::vector<size_n_max> &uncompleted_leaves);
    void print_remaining_stats(size_n_max leaves, size_k_max iterations, size_k_max phase, bool printHeader = false);

public:
    /// Budget of the last iteration -- the most expensive connection not through the root rises from a leaf to depth 1
    static size_k_max max_budget(size_k_max k, ParetoObjective objective, size_k_max penalty){
        uint64_t max_cost = (objective == ParetoObjective::RUNS) ? uint64_t(k - 1) + penalty
                                                                  : uint64_t(k - 1) + uint64_t(k > 2 ? k - 2 : 0) * penalty;
        if (max_cost + 1 >= std::numeric_limits<size_k_max>::max())
            throw std::invalid_argument("Penalty is too large for k-mer length " + std::to_string(k));
        return size_k_max(max_cost + 1);
    }

    LeafOnlyAC(const std::vector<kmer_t>& kmers, size_k_max k, bool use_complements, ParetoObjective objective, size_k_max penalty, const ParetoOptions &options) :
        K(k), N(kmers.size()), COMPLEMENTS(use_complements), OBJECTIVE(objective), PENALTY(penalty), MAX_BUDGET(max_budget(k, objective, penalty)),
        kMers(kmers), failureIndex(kMers, K) {};

    void compute_result();                  /// Runs the pareto optimization
    size_t print_result(std::ostream& of);  /// Prints the masked superstring, returns its objective value
};

template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::compute_result() {
    next.resize(N, INVALID_NODE);
    other_end.resize(N); for (size_n_max i = 0; i < N; ++i) other_end[i] = i;
    next_unused.resize(size_t(N) + 1); for (size_t i = 0; i <= size_t(N); ++i) next_unused[i] = size_n_max(i);
    remaining_priorities.resize(N, 0);
    skip_to.resize(N); for (size_n_max i = 0; i < N; ++i) skip_to[i] = i + 1;

    /// Phase 1 -- connect leaves overlapping by k - 1 characters
    print_remaining_stats(N, MAX_BUDGET, 1, true);
    for (size_n_max i = 0; i < N; ++i){
        size_n_max leaf = i;
        while (next[leaf] == INVALID_NODE && try_complete_leaf_phase_1(leaf)) leaf = next[leaf];
    }

    std::vector<size_n_max> uncompleted_leaves;
    for (size_n_max i = 0; i < N; ++i) if (next[i] == INVALID_NODE) uncompleted_leaves.push_back(i);
    std::shuffle(uncompleted_leaves.begin(), uncompleted_leaves.end(), std::default_random_engine(RANDOM_SEED));
    std::cerr << std::endl; WriteLog("Finished phase 1, " + std::to_string(uncompleted_leaves.size()) + " uncompleted leaves.");

    /// Phase 2 -- connect leaves via paths of increasing cost
    phase_2(uncompleted_leaves);
    std::cerr << std::endl; WriteLog("Finished computation of a masked superstring.");
}

template <typename kmer_t, typename size_n_max, typename priority_t>
inline bool LeafOnlyAC<kmer_t, size_n_max, priority_t>::try_complete_leaf_phase_1(size_n_max leaf_to_complete) {
    size_n_max first_failure_leaf = failureIndex.find_first_failure_leaf_by_index(leaf_to_complete, K - 1);
    if (first_failure_leaf == INVALID_NODE) return false;

    kmer_t leaf_complement = COMPLEMENTS ? complement_kmer(leaf_to_complete) : kmer_t(0);
    kmer_t failure_prefix = BitPrefix(kMers[first_failure_leaf], K, K - 1);
    for (size_n_max i = find_unused(first_failure_leaf); i < N && failure_prefix == BitPrefix(kMers[i], K, K - 1); i = find_unused(i + 1)){
        if (!can_connect(leaf_to_complete, leaf_complement, i)) continue;
        connect_leaves(leaf_to_complete, i, {});
        return true;
    }
    return false;
}

/// Iterative deepening: every iteration searches from all uncompleted leaves with the budget increased by one.
/// Budget 2 only finds overlaps of length k - 1, which phase 1 already found.
template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::phase_2(std::vector<size_n_max> &uncompleted_leaves) {
    for (size_k_max budget = 3; budget <= MAX_BUDGET; ++budget){
        size_t uncompleted_count = uncompleted_leaves.size();
        print_remaining_stats(uncompleted_count, MAX_BUDGET - budget + 1, 2);
        if (uncompleted_count == 0) break;

        for (size_n_max leaf : uncompleted_leaves){
            if (next[leaf] != INVALID_NODE) continue;
            if (search(leaf, budget)) connect_leaves(leaf, context.target, context.path);
        }

        size_t remaining = 0;
        for (size_n_max leaf : uncompleted_leaves) if (next[leaf] == INVALID_NODE) uncompleted_leaves[remaining++] = leaf;
        uncompleted_leaves.resize(remaining);
    }
}

template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::print_remaining_stats(size_n_max leaves, size_k_max iterations, size_k_max phase, bool printHeader) {
    if (printHeader){
        std::cerr << "  " << std::setw(MAX_COUNT_WIDTH) << leaves << ' ' << std::setw(MAX_ITERS_WIDTH) << iterations << ":" << std::endl;
    }
    for (size_k_max x = 0; x < 2 + MAX_COUNT_WIDTH + 1 + MAX_ITERS_WIDTH; ++x) std::cerr << '\b';
    std::cerr << phase << ' ' << std::setw(MAX_COUNT_WIDTH) << leaves << ' ' << std::setw(MAX_ITERS_WIDTH) << iterations;
    std::cerr.flush();
}

/// Search for a leaf to connect to the chain ending with origin; finds a connection iff its cost is less than budget
template <typename kmer_t, typename size_n_max, typename priority_t>
inline bool LeafOnlyAC<kmer_t, size_n_max, priority_t>::search(size_n_max origin, size_k_max budget) {
    context.stack.clear();
    context.log.clear();
    context.log.push_back({origin, 0, 0, true, NO_PARENT});
    push_failure(budget, K, origin, 0);

    kmer_t origin_complement = COMPLEMENTS ? complement_kmer(origin) : kmer_t(0);
    const bool use_skip_to = budget >= PENALTY;

    while (!context.stack.empty()){
        SearchEntry entry = context.stack.back(); context.stack.pop_back();
        size_k_max priority = entry.priority, depth = entry.depth;
        kmer_t prefix = BitPrefix(kMers[entry.node], K, depth);

        /// Unused leaf in the subtree (at most two unsuitable ones: the own chain start and the complement)
        for (size_n_max i = find_unused(entry.node); i < N && prefix == BitPrefix(kMers[i], K, depth); i = find_unused(i + 1)){
            if (!can_connect(origin, origin_complement, i)) continue;
            context.target = i;
            restore_unexpanded();
            context.path.clear();
            for (size_t log = entry.log; log != 0; log = context.log[log].parent) context.path.push_back(context.log[log].leaf);
            std::reverse(context.path.begin(), context.path.end());
            return true;
        }

        if (priority < 2) continue; /// Every further step costs at least one

        size_t log = entry.log;
        push_failure(priority, depth, context.log[log].leaf, log); /// Failure of the current node

        for (size_n_max i = entry.node; i < N && prefix == BitPrefix(kMers[i], K, depth); ++i){
            if (i == origin) continue;

            priority_t remaining = remaining_priorities[i];
            priority_t new_remaining = priority_t(std::min(priority, PRIORITY_SATURATION));
            if (remaining >= priority){
                if (use_skip_to){ /// Skip the whole run of leaves with high enough priority (with path compression)
                    size_n_max j = i;
                    while (j < N && remaining_priorities[j] >= priority && prefix == BitPrefix(kMers[j], K, depth)) j = skip_to[j];
                    for (size_n_max s = i; s < j; ){ size_n_max next_s = skip_to[s]; skip_to[s] = j; s = next_s; }
                    i = j - 1; /// j itself has to be examined
                }
                continue;
            }
            remaining_priorities[i] = new_remaining;
            if (use_skip_to && (skip_to[i] >= N || remaining_priorities[skip_to[i]] < priority)) skip_to[i] = i + 1;
            context.log.push_back({i, remaining, new_remaining, false, log});
            push_failure(priority, K, i, context.log.size() - 1); /// Failure of leaf i
        }
    }
    return false;
}

/// Push the failure node of the given node (a suffix of the leaf) if the rise is affordable; the cost is accumulated
/// level by level, so that no failure is searched for when the rise cannot be afforded
template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::push_failure(size_k_max priority, size_k_max node_depth, size_n_max leaf, size_t log) {
    auto level_cost = [&](size_k_max failure_depth) -> uint32_t {
        if (OBJECTIVE == ParetoObjective::RUNS)  /// A run is interrupted when the rise gets below depth K - 1
            return 1 + ((node_depth >= K - 1 && failure_depth == K - 2) ? PENALTY : 0);
        return 1 + ((node_depth == K && failure_depth == K - 1) ? 0 : PENALTY); /// A zero for each level but the first from a leaf
    };

    size_k_max failure_depth = node_depth - 1;
    uint32_t cost = level_cost(failure_depth);
    if (cost >= priority) return;
    size_n_max failure_index = failureIndex.find_first_failure_leaf_by_index(leaf, failure_depth);
    while (failure_index == INVALID_NODE){
        if (--failure_depth == 0) return;
        cost += level_cost(failure_depth);
        if (cost >= priority) return;
        failure_index = failureIndex.find_first_failure_leaf_by_index(leaf, failure_depth);
    }
    context.stack.push_back({failure_index, size_k_max(priority - cost), failure_depth, log});
}

/// After a successful search, leaves whose subsearch was not finished get their previous remaining priority back,
/// so that they do not prune later searches
template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::restore_unexpanded() {
    for (const auto &entry : context.stack){
        for (size_t log = entry.log; log != 0 && !context.log[log].restored; log = context.log[log].parent){
            LogEntry &leaf = context.log[log];
            leaf.restored = true;
            if (remaining_priorities[leaf.leaf] == leaf.new_priority) remaining_priorities[leaf.leaf] = leaf.old_priority;
        }
    }
}

template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::connect_leaves(size_n_max uncompleted_leaf, size_n_max unused_leaf, const std::vector<size_n_max> &path){
    size_n_max chain_start = other_end[uncompleted_leaf], chain_end = other_end[unused_leaf];
    size_n_max uncompleted_complement = INVALID_NODE, unused_complement = INVALID_NODE;
    size_n_max complement_start = INVALID_NODE, complement_end = INVALID_NODE;
    if (COMPLEMENTS){
        uncompleted_complement = complement_of(uncompleted_leaf);
        unused_complement = complement_of(unused_leaf);
        complement_start = other_end[unused_complement];    /// Complement of chain_end
        complement_end = other_end[uncompleted_complement]; /// Complement of chain_start
    }

    mark_used(unused_leaf);
    next[uncompleted_leaf] = unused_leaf;
    other_end[chain_start] = chain_end;
    other_end[chain_end] = chain_start;
    if (!path.empty()) store_backtrack(uncompleted_leaf, path, false);

    if (COMPLEMENTS){ /// Connect complements inversely
        mark_used(uncompleted_complement);
        next[unused_complement] = uncompleted_complement;
        other_end[complement_start] = complement_end;
        other_end[complement_end] = complement_start;
        if (!path.empty()) store_backtrack(unused_complement, path, true);
    }
}

template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::store_backtrack(size_n_max leaf, const std::vector<size_n_max> &path, bool complemented) {
    backtrack_positions.emplace_back(leaf, backtracks.size());
    backtracks.push_back(size_n_max(path.size()));
    if (!complemented) for (size_n_max leaf_on_path : path) backtracks.push_back(leaf_on_path);
    else for (size_t i = path.size(); i-- > 0; ) backtracks.push_back(complement_of(path[i]));
}

template <typename kmer_t, typename size_n_max, typename priority_t>
inline void LeafOnlyAC<kmer_t, size_n_max, priority_t>::build_backtrack_lookup() {
    std::sort(backtrack_positions.begin(), backtrack_positions.end());
    has_backtrack.assign(size_t(N) / 64 + 1, 0);
    for (const auto &p : backtrack_positions) has_backtrack[p.first / 64] |= uint64_t(1) << (p.first % 64);
    has_backtrack_rank.resize(has_backtrack.size());
    size_t rank = 0;
    for (size_t w = 0; w < has_backtrack.size(); ++w){ has_backtrack_rank[w] = rank; rank += __builtin_popcountll(has_backtrack[w]); }
}

template <typename kmer_t, typename size_n_max, typename priority_t>
inline size_t LeafOnlyAC<kmer_t, size_n_max, priority_t>::print_result(std::ostream& of) {
    build_backtrack_lookup();

    size_t total_objective_value = 0;
    std::string buffer;
    const size_t BUFFER_LIMIT = size_t(1) << 20;
    buffer.reserve(BUFFER_LIMIT + K);

    kmer_t last_kmer = 0;
    bool first = true;
    /// Print the part of the last k-mer not covered by the next one
    auto append_kmer = [&](size_n_max leaf){
        kmer_t actual_kmer = kMers[leaf];
        if (first){
            first = false;
            last_kmer = actual_kmer;
            return;
        }
        size_k_max ov = GetMaxOverlapLength(last_kmer, actual_kmer, K);
        buffer.push_back(NucleotideAtIndex(last_kmer, K, 0));
        for (size_k_max c = 1; c < K - ov; ++c) buffer.push_back(letter_to_lower(NucleotideAtIndex(last_kmer, K, c)));
        if (buffer.size() >= BUFFER_LIMIT){ of.write(buffer.data(), buffer.size()); buffer.clear(); }

        last_kmer = actual_kmer;
        total_objective_value += K - ov;
        if (OBJECTIVE == ParetoObjective::RUNS)  if (ov < K - 1) total_objective_value += PENALTY;
        if (OBJECTIVE == ParetoObjective::ZEROS) total_objective_value += size_t(PENALTY) * (K - 1 - ov);
    };

    for (size_n_max i = 0; i < N; ++i){
        if (is_used(i)) continue; /// Not a chain start, or the complementary chain was printed
        if (COMPLEMENTS) mark_used(complement_of(other_end[i])); /// Do not print the complementary chain

        for (size_n_max actual = i; actual != INVALID_NODE; actual = next[actual]){
            append_kmer(actual);
            size_t position = find_backtrack(actual);
            if (position != std::numeric_limits<size_t>::max())
                for (size_n_max c = 1; c <= backtracks[position]; ++c) append_kmer(backtracks[position + c]);
        }
    }
    buffer.push_back(NucleotideAtIndex(last_kmer, K, 0));
    for (size_k_max c = 1; c < K; ++c) buffer.push_back(letter_to_lower(NucleotideAtIndex(last_kmer, K, c)));
    of.write(buffer.data(), buffer.size());

    total_objective_value += K;
    if (OBJECTIVE == ParetoObjective::RUNS)  total_objective_value += PENALTY;                   /// First run of ones
    if (OBJECTIVE == ParetoObjective::ZEROS) total_objective_value += size_t(PENALTY) * (K - 1); /// Zeros at the end
    return total_objective_value;
}
