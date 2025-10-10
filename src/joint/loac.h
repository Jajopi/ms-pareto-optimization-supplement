#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <random>

#include "../kmers.h"
#include "objective.h"
#include "index.h"
#include "unionfind.h"

size_t RANDOM_SEED = 0;
size_t MAX_COUNT_WIDTH = 12;
size_t MAX_ITERS_WIDTH = 3;

/// The data structure for efficient heuristic search
template <typename kmer_t, typename size_n_max>
class LeafOnlyAC {
    using stack_t = std::vector<std::tuple<size_k_max, size_k_max, size_n_max, size_n_max>>;
    static const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();

    size_k_max K;                                   /// Length of a kmer
    size_n_max N;                                   /// Number of kmers (= number of leaves of a trie)
    bool COMPLEMENTS;                               /// Whether computing with complements (in bi-directional model)
    JointObjective OBJECTIVE;                       /// What to apply penalty for -- RUNS or ZEROS
    size_k_max PENALTY;                             /// Value of penalty (for runs / zeros) used in the computation
    
    const std::vector<kmer_t>&       kMers;         /// Sorted leaves of the trie -- kmers
    FailureIndex<kmer_t, size_n_max> failureIndex;  /// Index for fast retrieval of first failure node of given leaf and depth
    UnionFind<size_n_max>            components;    /// Unionfind, tracking ends of leaf chains
    
    std::vector<size_n_max> complements;            /// Index of complement of each kmer
    std::vector<size_n_max> backtracks;             /// Indexes of kmers appearing between leaves in a chain
    std::vector<size_n_max> backtrack_indexes;      /// Index into backtrack_indexes for each leaf -- TODO make more efficient
    std::vector<size_n_max> previous;               /// Previous visited leaf for each leaf, used in DFS
    std::vector<size_n_max> next;                   /// Index of the next leaf in a chain a leaf
    std::vector<size_k_max> remaining_priorities;   /// Speedup index preventing repetitive computation -- TODO explain more specifically
    std::vector<size_n_max> skip_to;                /// Speedup index preventing repetitive computation -- TODO explain more specifically
    std::vector<bool>       used;                   /// Used leaf has already been used as a next leaf for exactly one other leaf
    std::vector<size_k_max> no_unused;              /// Speedup index -- the smallest depth of node on failure path which has no unused leaves
    
    void construct_complements();
    bool try_complete_leaf_phase_1(size_n_max leaf_to_connect); /// First iteration, "get simplitigs"
    bool try_complete_leaf_phase_2(size_n_max leaf_to_connect, size_k_max priority_drop_limit, stack_t &stack); /// Until threshold, TODO parallel, blind DFS
    bool try_complete_leaf_phase_3(size_n_max leaf_to_connect, size_k_max priority_drop_limit, std::vector<size_n_max> &unused_leaves); /// Informed search
    void try_push_failure_node_into_stack(stack_t &stack, size_k_max priority, size_k_max node_depth, size_n_max node_index, size_n_max last_leaf);
    void squeeze_sparse_list(std::vector<size_n_max>& sparse_list);
    void set_backtrack_path_for_leaf(size_n_max origin_leaf, size_n_max next_leaf);
    void print_remaining_stats(size_n_max leaves, size_k_max iterations, size_k_max phase);

    bool RESULT_COMPUTED = false;
public:
    std::ostream& LOG_STREAM = std::cerr;

    LeafOnlyAC(const std::vector<kmer_t>& kmers, size_k_max k, bool complement, JointObjective objective, size_k_max penalty) :
        K(k), N(kmers.size()), PENALTY(penalty), kMers(kmers), failureIndex(kMers, K) { construct_complements(); };

    void compute_result();                      /// Runs joint optimization
    void optimize_result();                     /// Runs some other reoptimization -- currently none, possible TODO to implement
    size_n_max print_result(std::ostream& of);  /// Prints the masked superstring to output stream, return the value of objective function for printed result
};

// Constructing

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::construct_complements() {
    if (COMPLEMENTS){
        std::vector<std::pair<kmer_t, size_n_max>> complement_kmers(N);
        for (size_n_max i = 0; i < N; ++i){
            complement_kmers[i] = std::make_pair(ReverseComplement(kMers[i], K), N - i); /// For even k, swapping indexes for same pairs of kmers is needed
        }
        
        std::sort(complement_kmers.begin(), complement_kmers.end());

        complements.resize(N);
        for (size_n_max i = 0; i < N; ++i) complements[i] = N - complement_kmers[i].second;
    }
}

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::compute_result() {
    if (RESULT_COMPUTED){
        throw std::invalid_argument("Result has already been computed.");
    }

    components = UnionFind(N);
    backtracks.reserve(N); /// TODO more clever
    backtrack_indexes.resize(N, INVALID_NODE);
    previous.resize(N, INVALID_NODE);
    next.resize(N, INVALID_NODE);
    remaining_priorities.resize(N, 0);
    skip_to.resize(N); for (size_n_max i = 0; i < N; ++i) skip_to[i] = i + 1;
    used.resize(N, false);
    no_unused.resize(N, K);

    std::vector<size_n_max> uncompleted_leaves; uncompleted_leaves.reserve(N);
    
    size_k_max max_priority_drop = 0;
    if (OBJECTIVE == JointObjective::RUNS)  max_priority_drop = (K - 1) + PENALTY;
    if (OBJECTIVE == JointObjective::ZEROS) max_priority_drop = (K - 1) * (PENALTY + 1);

    /// Phase 1
    LOG_STREAM << 1 << ' ' << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << max_priority_drop << std::endl;
    
    for (size_n_max i = 0; i < N; ++i){
        if (next[i] != INVALID_NODE) continue;
        size_n_max leaf_index = i;
        bool result = try_complete_leaf_phase_1(leaf_index);
        while (result){
            leaf_index = next[leaf_index];
            result = try_complete_leaf_phase_1(leaf_index);
        }
        uncompleted_leaves.push_back(i);
    }

    std::shuffle(uncompleted_leaves.begin(), uncompleted_leaves.end(), std::default_random_engine(RANDOM_SEED));

    size_k_max priority_drop_limit;
    /// Phase 2
    {
        stack_t stack; // stack.reserve() /// TODO based on phase threshold
        for (priority_drop_limit = 2; priority_drop_limit <= max_priority_drop; ++priority_drop_limit){
            size_n_max uncompleted_leaf_count = uncompleted_leaves.size();

            print_remaining_stats(uncompleted_leaf_count, max_priority_drop - priority_drop_limit + 1, 2);

            for (size_n_max i = 0; i < uncompleted_leaf_count; ++i){
                if (next[i] != INVALID_NODE) continue;
                size_n_max leaf_index = uncompleted_leaves[i];
                bool result = try_complete_leaf_phase_2(leaf_index, priority_drop_limit, stack);
                while (result){
                    uncompleted_leaves[leaf_index] = INVALID_NODE;
                    leaf_index = next[leaf_index];
                    result = try_complete_leaf_phase_2(leaf_index, priority_drop_limit, stack);
                }
            }

            squeeze_sparse_list(uncompleted_leaves);
            stack.clear();
        }
    }

    /// Phase 3
    std::vector<size_n_max> unused_leaves; unused_leaves.reserve(uncompleted_leaves.size()); /// Sorted
    for (size_n_max i = 0; i < N; ++i){ if (!used[i]) unused_leaves.push_back(i); }

    for (/* intentionally not reseting priority_drop_limit */; priority_drop_limit <= max_priority_drop; ++priority_drop_limit){
        size_n_max uncompleted_leaf_count = uncompleted_leaves.size();

        print_remaining_stats(uncompleted_leaf_count, max_priority_drop - priority_drop_limit + 1, 3);

        for (size_n_max i = 0; i < uncompleted_leaf_count; ++i){
            bool result = try_complete_leaf_phase_3(uncompleted_leaves[i], priority_drop_limit, unused_leaves);
            if (result) uncompleted_leaves[i] = INVALID_NODE;
        }

        squeeze_sparse_list(uncompleted_leaves);
        squeeze_sparse_list(unused_leaves);
    }

    RESULT_COMPUTED = true;
    WriteLog("Finished computation of a masked superstring.");
}

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::optimize_result() {
    if (!RESULT_COMPUTED){
        throw std::invalid_argument("Result has not been computed yet.");
    }
    // WriteLog("Starting secondary optimization of the masked superstring.");
    WriteLog("No secondary optimization implemented yet.");
    // WriteLog("Finished secondary optimization of the masked superstring.");
}

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::print_remaining_stats(size_n_max leaves, size_k_max iterations, size_k_max phase) {
    for (size_k_max x = 0; x < 2 + MAX_COUNT_WIDTH + 1 + MAX_ITERS_WIDTH; ++x) LOG_STREAM << '\b';
    LOG_STREAM << phase << ' ' << std::setw(MAX_COUNT_WIDTH) << leaves << ' ' << std::setw(MAX_ITERS_WIDTH) << iterations--;
    LOG_STREAM.flush();
}

template <typename kmer_t, typename size_n_max>
inline bool LeafOnlyAC<kmer_t, size_n_max>::try_complete_leaf_phase_1(size_n_max leaf_to_complete) {
    size_n_max first_failure_leaf = failureIndex.find_first_failure_leaf_by_index(leaf_to_complete, K - 1);
    if (first_failure_leaf == INVALID_NODE) return false;

    size_n_max leaf_complement = COMPLEMENTS ? complements[leaf_to_complete] : INVALID_NODE;

    auto failure_prefix = BitPrefix(kMers[first_failure_leaf], K, K - 1);
    for (size_n_max i = first_failure_leaf; i < N && failure_prefix == BitPrefix(kMers[i], K, K - 1); ++i){
        if (used[i] || components.are_connected(leaf_to_complete, i) || leaf_complement == i) continue;

        /// Suitable leaf found
        used[i] = true;
        next[leaf_to_complete] = i;
        components.connect(leaf_to_complete, i); /// Second one pointing at the first one, order matters

        if (COMPLEMENTS){ /// Connect complements inversely
            used[leaf_complement] = true;
            next[complements[i]] = leaf_complement;
            components.connect(complements[i], leaf_complement);
        }
        return true;
    }
    return false;
}

template <typename kmer_t, typename size_n_max>
inline bool LeafOnlyAC<kmer_t, size_n_max>::try_complete_leaf_phase_2(size_n_max leaf_to_complete, size_k_max priority_drop_limit, stack_t &stack) {
    try_push_failure_node_into_stack(stack, priority_drop_limit, K, leaf_to_complete, leaf_to_complete);

    while (!stack.empty()){
        auto t = stack.back(); stack.pop_back();
        size_k_max priority    = std::get<0>(t);
        size_k_max chain_depth = std::get<1>(t);
        size_n_max leaf_index  = std::get<2>(t);
        size_n_max last_leaf   = std::get<3>(t);

        size_n_max leaf_complement = COMPLEMENTS ? complements[leaf_to_complete] : INVALID_NODE;

        kmer_t leaf_prefix = BitPrefix(kMers[leaf_index], K, chain_depth);

        if (no_unused[leaf_index] > chain_depth){ /// There is at least one unused leaf to be found
            bool skipped_unused = false;

            for (size_n_max i = leaf_index; i < N && leaf_prefix == BitPrefix(kMers[i], K, chain_depth); ++i){
                if (used[i]) continue;
                if (components.are_connected(leaf_to_complete, i) || leaf_complement == i){
                    skipped_unused = true;
                    continue;
                }

                /// Suitable leaf found
                used[i] = true;
                next[leaf_to_complete] = i;
                components.connect(leaf_to_complete, i); /// Second one pointing at the first one, order matters

                if (COMPLEMENTS){ /// Connect complements inversely
                    used[leaf_complement] = true;
                    next[complements[i]] = leaf_complement;
                    components.connect(complements[i], leaf_complement);
                }

                if (last_leaf != leaf_to_complete){
                    previous[i] = last_leaf;
                    set_backtrack_path_for_leaf(leaf_to_complete, i);
                }
                return true;
            }
            if (!skipped_unused) no_unused[leaf_index] = chain_depth;
        }

        /// Append into stack
        if (priority == 0) continue;
        try_push_failure_node_into_stack(stack, priority, chain_depth, last_leaf, last_leaf); /// Add failure of current node

        for (size_n_max i = leaf_index; i < N && leaf_prefix == BitPrefix(kMers[i], K, chain_depth); ++i){
            if (i == leaf_to_complete) continue;

            if (priority_drop_limit >= PENALTY){
                if (remaining_priorities[i] >= priority){
                    std::vector<size_n_max> skipped;
                    size_n_max j = i;
                    while (j < N && remaining_priorities[j] >= priority && leaf_prefix == BitPrefix(kMers[j], K, chain_depth)){
                        skipped.push_back(j);
                        j = skip_to[j];
                    }
                    for (size_n_max s : skipped){
                        skip_to[s] = j;
                    }
                    i = j;
                    continue;
                }
                remaining_priorities[i] = priority;
                if (remaining_priorities[skip_to[i]] < priority) skip_to[i] = i + 1;
            } else {
                if (remaining_priorities[i] >= priority) continue;
                remaining_priorities[i] = priority;
            }

            previous[i] = last_leaf;
            try_push_failure_node_into_stack(stack, priority, K, i, i); /// Add failure of leaf i
        }
    }
    return false;
}

template <typename kmer_t, typename size_n_max>
inline bool LeafOnlyAC<kmer_t, size_n_max>::try_complete_leaf_phase_3(size_n_max leaf_to_complete, size_k_max priority_drop_limit, std::vector<size_n_max> &unused_leaves) {
    return false;
}

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::try_push_failure_node_into_stack(stack_t &stack,
        size_k_max priority, size_k_max node_depth, size_n_max node_index, size_n_max last_leaf) {
    if (--priority == 0) return;

    size_k_max failure_depth = node_depth - 1;
    size_n_max failure_index = failureIndex.find_first_failure_leaf_by_index(node_index, failure_depth);
    while (failure_index == INVALID_NODE){
        if (--failure_depth == 0) return;
        if (--priority == 0) return;
        failure_index = failureIndex.find_first_failure_leaf_by_index(node_index, failure_depth);
    }

    if (OBJECTIVE == JointObjective::RUNS){
        if (node_depth == K - 1 || (node_depth == K && failure_depth < K - 1)){ /// Run will be interrupted
            if (priority < PENALTY) return;
            priority -= PENALTY;
        }
    }
    if (OBJECTIVE == JointObjective::ZEROS){
        size_k_max priority_decrease = (node_depth - failure_depth - (node_depth == K)) * PENALTY;
        if (priority < priority_decrease) return;
        priority -= priority_decrease;
    }

    stack.emplace_back(priority, failure_depth, failure_index, last_leaf);
}

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::squeeze_sparse_list(std::vector<size_n_max> &sparse_list) {
    size_n_max count = sparse_list.size(), shift = 0;

    for (size_n_max i = 0; i < count; ++i){
        if (sparse_list[i] == INVALID_NODE) ++shift;
        else sparse_list[i - shift] = sparse_list[i];
    }

    sparse_list.resize(count - shift);
}

template <typename kmer_t, typename size_n_max>
inline void LeafOnlyAC<kmer_t, size_n_max>::set_backtrack_path_for_leaf(size_n_max origin_leaf, size_n_max next_leaf) {
    size_n_max last_size = backtracks.size();

    size_n_max actual = next_leaf;
    while (actual != origin_leaf){
        backtracks.push_back(actual);
        actual = previous[actual];
    }
    backtrack_indexes[origin_leaf] = backtracks.size() - 1;

    if (COMPLEMENTS){
        size_n_max new_size = backtracks.size();
        size_n_max count = new_size - last_size;
        for (size_n_max i = 0; i < count; ++i) backtracks.push_back(INVALID_NODE);
        
        size_n_max index = new_size + count - 1;
        size_n_max actual = previous[next_leaf];
        while (index >= last_size + count){
            backtracks[index--] = complements[actual];
            actual = previous[actual];
        }
        backtrack_indexes[complements[next_leaf]] = backtracks.size() - 1;
    }
}

template <typename kmer_t, typename size_n_max>
inline size_n_max LeafOnlyAC<kmer_t, size_n_max>::print_result(std::ostream& of) {
    if (!RESULT_COMPUTED){
        throw std::invalid_argument("Result has not been computed yet.");
    }

    size_n_max total_objective_value = 0;

    size_n_max actual = INVALID_NODE;
    kmer_t last_kmer = 0;
    bool first = true;
    for (size_n_max i = 0; i < N; ++i){
        if (components.find(i) != i) continue; /// Skipping leaf as it is not the first leaf of a chain

        if (COMPLEMENTS) components.connect(i, components.find(complements[i])); /// Prevent complementary chain from being printed later

        actual = i;
        while (actual != INVALID_NODE){
            if (first){
                first = false;
                last_kmer = kMers[i];
                actual = next[i];
                continue;
            }

            kmer_t actual_kmer = kMers[actual];
            size_k_max ov = GetMaxOverlapLength(last_kmer, actual_kmer, K);
            PrintKmerMasked(last_kmer, K, of, size_k_max(K - ov));

            last_kmer = actual_kmer;
            total_objective_value += K - ov;
            if (OBJECTIVE == JointObjective::RUNS)  if (ov < K - 1) total_objective_value += PENALTY;
            if (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1 - ov);

            if (backtrack_indexes[actual] != INVALID_NODE){
                size_n_max backtrack_index = backtrack_indexes[actual];
                size_n_max actual_backtrack = backtracks[backtrack_index];
                size_n_max next_node = next[actual];
                while (actual_backtrack != next_node){
                    kmer_t actual_kmer = kMers[actual_backtrack];
                    size_k_max ov = GetMaxOverlapLength(last_kmer, actual_kmer, K);
                    PrintKmerMasked(last_kmer, K, of, size_k_max(K - ov));

                    last_kmer = actual_kmer;
                    total_objective_value += K - ov;
                    if (OBJECTIVE == JointObjective::RUNS)  if (ov < K - 1) total_objective_value += PENALTY;
                    if (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1 - ov);

                    --backtrack_index;
                    actual_backtrack = backtracks[backtrack_index];
                }
            }

            if (next[actual] == actual) break;
            actual = next[actual];
        }
    }
    PrintKmerMasked(last_kmer, K, of, K);
    total_objective_value += K;
    if (OBJECTIVE == JointObjective::RUNS)  total_objective_value += PENALTY;           /// First run of ones
    if (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1); /// Zeros at the end

    return total_objective_value;
}
