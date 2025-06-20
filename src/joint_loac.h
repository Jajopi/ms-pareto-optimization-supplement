#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <random>
#include <queue>

#include "kmers.h"
#include "parser.h"
#include "joint.h"

typedef uint8_t size_k_max;

constexpr size_t RANDOM_SEED = 0;
constexpr size_t MAX_COUNT_WIDTH = 12;
constexpr size_t MAX_ITERS_WIDTH = 3;

/// Union-find with non-comutative union operation
template<typename size_n_max>
class UnionFind {
    std::vector<size_n_max> roots;
    size_n_max component_count;
public:
    UnionFind(size_n_max size) : roots(size), component_count(size) {
        for (size_n_max i = 0; i < size; ++i) roots[i] = i;
    };
    UnionFind() = default;
    
    inline size_n_max find(size_n_max x) {
        size_n_max root = roots[x];
        if (roots[root] == root) return root;

        while (roots[root] != root) root = roots[root];
        while (x != root){
            size_n_max new_x = roots[x];
            roots[x] = root;
            x = new_x;
        }
        return root;
    }

    inline bool are_connected(size_n_max x, size_n_max y){
        return find(x) == find(y);
    }

    inline void connect(size_n_max to, size_n_max from){ // Second one points to the first one - points to the begining of a chain
        if (are_connected(from, to)) return;
        roots[from] = to;
        --component_count;
    }

    inline size_n_max count() const { return component_count; };
};


/// Index for fast searching of k-mers
/// Search operation returns std::numeric_limits<size_n_max>::max() when k-mer is not present in the set
template<typename kmer_t, typename size_n_max>
class FailureIndex {
    static constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();

    const std::vector<kmer_t> &kMers;
    size_n_max N;
    size_k_max K;

    size_k_max SPEEDUP_DEPTH;
    std::vector<size_n_max> search_speedup;

    std::vector<size_n_max> first_row;

    inline void construct_index(){
        N = kMers.size();
        SPEEDUP_DEPTH = log2(N) / 2;
        size_n_max speedup_size = (size_n_max(1) << 2 * SPEEDUP_DEPTH);
        search_speedup.resize(speedup_size + 1);

        size_n_max index = 0;
        for (kmer_t k = 0; k < kmer_t(speedup_size); ++k){
            search_speedup[k] = index;
            while (index < N && BitPrefix(kMers[index], K, SPEEDUP_DEPTH) == k) ++index;
        }
        search_speedup[speedup_size] = N;

        first_row.resize(N);
        for (size_n_max i = 0; i < N; ++i) first_row[i] = search(i, K - 1);
    }

    inline size_n_max search(size_n_max index, size_k_max depth){
        kmer_t searched = BitSuffix(kMers[index], depth);

        if (depth < SPEEDUP_DEPTH){
            size_n_max begin = search_speedup[searched << 2 * (SPEEDUP_DEPTH - depth)];
            return (BitPrefix(kMers[begin], K, depth) == searched) ? begin : INVALID_NODE;
        }
        
        kmer_t speedup_index = BitPrefix(searched, depth, SPEEDUP_DEPTH);
        size_n_max begin = search_speedup[speedup_index];
        size_n_max end = search_speedup[speedup_index + 1];

        /// Switch to bin-search on big intervals
        if (end - begin > 16){
            while (begin < end){
                size_n_max middle = (begin + end - 1) / 2;
                kmer_t current = BitPrefix(kMers[middle], K, depth);

                if (current < searched) begin = middle + 1;
                else end = middle;
            }
            return (BitPrefix(kMers[begin], K, depth) == searched) ? begin : INVALID_NODE;
        }

        for (size_n_max i = begin; i < end; ++i){
            kmer_t current = BitPrefix(kMers[i], K, depth);
            if (current == searched) return i;
            if (current > searched) return INVALID_NODE;
        }
        return INVALID_NODE;
    }
public:
    inline FailureIndex(const std::vector<kmer_t> &kmers, size_k_max k) :
            kMers(kmers), N(kmers.size()), K(k) {
        construct_index();
        WriteLog("Finished constructing index.");
    }

    inline size_n_max find_first_failure_leaf(size_n_max index, size_k_max depth){
        if (depth == K - 1) return first_row[index];
        return search(index, depth);
    }
};


/// The data structure for efficient heuristic search
template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
class LeafOnlyAC {
    static constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();

    size_k_max K;               /// Kmer-length
    size_n_max N;               /// Number of k-mers (number of leaves)
    size_k_max PENALTY;         /// Value of penalty (for runs / zeros) used in the computation
    
    std::vector<kmer_t> kMers;  /// Sorted leaves
    FailureIndex<kmer_t, size_n_max> failureIndex;
    
    std::vector<size_n_max> complements;
    UnionFind<size_n_max> components;
    std::vector<std::tuple<size_k_max, size_k_max, size_n_max, size_n_max>> stack;
    std::vector<size_n_max> backtracks;
    std::vector<size_n_max> backtrack_indexes;
    std::vector<size_n_max> previous;
    std::vector<size_n_max> next;
    std::vector<size_k_max> remaining_priorities;
    std::vector<size_n_max> skip_to;
    std::vector<bool> used;
    std::vector<size_k_max> no_unused;

    bool try_complete_leaf(size_n_max leaf_to_connect, size_k_max priority_drop_limit);
    void push_failure_node_into_stack(size_k_max priority, size_k_max node_depth, size_n_max node_index, size_n_max last_leaf);
    void squeeze_uncompleted_leaves(std::vector<size_n_max>& unclompleted_leaves);
    void set_backtrack_path_for_leaf(size_n_max origin_leaf, size_n_max next_leaf);

    bool RESULT_COMPUTED = false;
public:
    std::ostream& LOG_STREAM = std::cerr;

    LeafOnlyAC(const std::vector<kmer_t>& kmers, size_k_max k, size_k_max penalty) :
        K(k), N(kmers.size()), PENALTY(penalty), kMers(kmers), failureIndex(kMers, K) {
            if constexpr(COMPLEMENTS) construct_complements();
        };
    LeafOnlyAC(std::vector<kmer_t>&& kmers, size_k_max k, size_k_max penalty) :
        K(k), N(kmers.size()), PENALTY(penalty), kMers(std::move(kmers)), failureIndex(kMers, K) {
            if constexpr(COMPLEMENTS) construct_complements();
        };

    void construct_complements();
    void compute_result();
    void optimize_result();
    size_n_max print_result(std::ostream& of);
};

// Constructing

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::construct_complements() {
    if constexpr (COMPLEMENTS){
        complements.resize(N);
        bool even_k = K % 2 == 0;
    
        std::vector<std::pair<kmer_t, size_n_max>> complement_kmers(N);
        for (size_n_max i = 0; i < N; ++i){
            complement_kmers[i] = std::make_pair(ReverseComplement(kMers[i], K),
                                                 even_k ? N - i : i); /// For even k, swap indexes for same pairs of kmers
        }
    
        std::sort(complement_kmers.begin(), complement_kmers.end());
        
        for (size_n_max i = 0; i < N; ++i){
            complements[i] = even_k ? N - complement_kmers[i].second : complement_kmers[i].second;
        }
    }
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::compute_result() {
    if (RESULT_COMPUTED){
        throw std::invalid_argument("Result has already been computed.");
    }
    size_n_max SEARCH_CUTOFF = N / (1 << 10);

    components = UnionFind(N);
    backtracks.reserve(N); /// Chains for leaves where backtracking is needed
    backtrack_indexes.resize(N, INVALID_NODE); /// Indexes into backtracks for each leaf
    previous.resize(N, INVALID_NODE);
    next.resize(N, INVALID_NODE);
    remaining_priorities.resize(N, 0);
    skip_to.resize(N);
    for (size_n_max i = 0; i < N; ++i) skip_to[i] = i + 1;
    used.resize(N, false);
    no_unused.resize(N, K);

    std::vector<size_n_max> uncompleted_leaves(N);
    for (size_n_max i = 0; i < N; ++i) uncompleted_leaves[i] = i;
    std::shuffle(uncompleted_leaves.begin(), uncompleted_leaves.end(), std::default_random_engine(RANDOM_SEED));
    size_n_max next_preffered_leaf = INVALID_NODE;
    
    size_k_max max_priority_drop = (K - 1) + PENALTY;
    // size_k_max max_priority_drop = (K - 1);
    
    size_n_max remaining_iterations = max_priority_drop;
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations << std::endl;
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();

    for (size_k_max priority_drop_limit = 1;
        priority_drop_limit <= max_priority_drop;
        ++priority_drop_limit){
            size_n_max uncompleted_leaf_count = uncompleted_leaves.size();

            for (size_k_max x = 0; x < MAX_COUNT_WIDTH + 1 + MAX_ITERS_WIDTH; ++x) LOG_STREAM << '\b';
            LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << uncompleted_leaf_count
            << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations--; LOG_STREAM.flush();
            
        if (uncompleted_leaf_count < SEARCH_CUTOFF) break;
        
        for (size_n_max i = 0; i < uncompleted_leaf_count; ++i){
            size_n_max leaf_index = uncompleted_leaves[i];
            
            if (next_preffered_leaf != INVALID_NODE){
                leaf_index = next_preffered_leaf;
                --i;
                next_preffered_leaf = INVALID_NODE;

                if (next[leaf_index] != INVALID_NODE){
                    continue;
                }
    
                bool result = try_complete_leaf(leaf_index,
                    (priority_drop_limit < max_priority_drop) ? priority_drop_limit : (K - 1 + PENALTY));
                if (result){
                    next_preffered_leaf = next[leaf_index];
                }
            }
            else{
                if (leaf_index == INVALID_NODE) continue;
                if (next[leaf_index] != INVALID_NODE){
                    uncompleted_leaves[i] = INVALID_NODE;
                    continue;
                }
    
                bool result = try_complete_leaf(leaf_index,
                    (priority_drop_limit < max_priority_drop) ? priority_drop_limit : (K - 1 + PENALTY));
                if (result){
                    next_preffered_leaf = next[leaf_index];
                    uncompleted_leaves[i] = INVALID_NODE;
                }
            }
            
        }

        squeeze_uncompleted_leaves(uncompleted_leaves);
    }

    for (size_k_max x = 0; x < MAX_COUNT_WIDTH + 1 + MAX_ITERS_WIDTH; ++x) LOG_STREAM << '\b';
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << uncompleted_leaves.size()
        << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations << std::endl;

    RESULT_COMPUTED = true;
    WriteLog("Finished computation of a masked superstring.");
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::optimize_result()
{
    if (!RESULT_COMPUTED){
        throw std::invalid_argument("Result has not been computed yet.");
    }
    WriteLog("Starting optimization of the masked superstring.");

    WriteLog("Finished optimization of the masked superstring.");
}

// Internal functions

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline bool LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::try_complete_leaf(
        size_n_max leaf_to_complete, size_k_max priority_drop_limit) {

    if (priority_drop_limit == 1){
        size_n_max first_failure_leaf = failureIndex.find_first_failure_leaf(leaf_to_complete, K - 1);

        if (first_failure_leaf == INVALID_NODE){
            return false;
        }
        size_n_max leaf_complement;
        if constexpr (COMPLEMENTS) leaf_complement = complements[leaf_to_complete];

        for (size_n_max i = first_failure_leaf; i < N &&
                BitPrefix(kMers[first_failure_leaf], K, K - 1) == BitPrefix(kMers[i], K, K - 1);
                ++i){
            if (used[i] || components.are_connected(leaf_to_complete, i)) continue;
            if constexpr (COMPLEMENTS) if (leaf_complement == i) continue;

            used[i] = true;
            next[leaf_to_complete] = i;
            components.connect(leaf_to_complete, i); /// Second one pointing at the first one, order matters

            if constexpr (COMPLEMENTS){ /// Connect complements inversely
                used[leaf_complement] = true;
                next[complements[i]] = leaf_complement;
                components.connect(complements[i], leaf_complement);
            }
            
            return true;
        }

        return false;
    }

    stack.clear();
    push_failure_node_into_stack(priority_drop_limit, K, leaf_to_complete, leaf_to_complete);

    while (!stack.empty()){
        auto t = stack.back(); stack.pop_back();
        size_k_max priority = std::get<0>(t);
        size_k_max chain_depth = std::get<1>(t);
        size_n_max leaf_index = std::get<2>(t);
        size_n_max last_leaf = std::get<3>(t);

        size_n_max leaf_complement;
        if constexpr (COMPLEMENTS) leaf_complement = complements[leaf_to_complete];

        kmer_t leaf_prefix = BitPrefix(kMers[leaf_index], K, chain_depth);

        if (no_unused[leaf_index] > chain_depth){ /// There is at least one unused leaf to be found
            bool skipped_unused = false;

            for (size_n_max i = leaf_index; i < N && leaf_prefix == BitPrefix(kMers[i], K, chain_depth); ++i){
                if (used[i]) continue;
                if constexpr (COMPLEMENTS) if (leaf_complement == i) continue;
                if (components.are_connected(leaf_to_complete, i)){
                    skipped_unused = true;
                    continue;
                }

                used[i] = true;
                next[leaf_to_complete] = i;
                components.connect(leaf_to_complete, i); /// Second one pointing at the first one, order matters

                if constexpr (COMPLEMENTS){ /// Connect complements inversely
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
            }
            else {
                if (remaining_priorities[i] >= priority) continue;
                remaining_priorities[i] = priority;
            }

            previous[i] = last_leaf;
            push_failure_node_into_stack(priority, K, i, i); /// Add failure of that leaf
        }

        push_failure_node_into_stack(priority, chain_depth, last_leaf, last_leaf); /// Add failure of current node
    }

    return false;
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::push_failure_node_into_stack(
    size_k_max priority, size_k_max node_depth, size_n_max node_index, size_n_max last_leaf){

    size_k_max failure_depth = node_depth - 1;
    size_n_max failure_index = INVALID_NODE;

    while (failure_depth != 0){
        failure_index = failureIndex.find_first_failure_leaf(node_index, failure_depth);
        
        if (failure_index != INVALID_NODE) break;
        --failure_depth;
    }
    if (failure_depth == 0) return;

    if constexpr (OBJECTIVE == JointObjective::RUNS){
        if (priority < node_depth - failure_depth) return;
        priority -= (node_depth - failure_depth);

        if (node_depth == K - 1 || (node_depth == K && failure_depth < K - 1)){ /// Run will be interrupted
            if (priority < PENALTY) return;
            priority -= PENALTY;
        }
    }
    if constexpr (OBJECTIVE == JointObjective::ZEROS){
        if (priority < (node_depth - failure_depth) * 2 - 1) return;
        priority -= (node_depth - failure_depth) * 2 - 1;
    }

    stack.emplace_back(priority, failure_depth, failure_index, last_leaf);
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::squeeze_uncompleted_leaves(std::vector<size_n_max> &uncompleted_leaves) {
    size_n_max count = uncompleted_leaves.size(), shift = 0;
    
    for (size_n_max i = 0; i < count; ++i){
        if (uncompleted_leaves[i] == INVALID_NODE) ++shift;
        else uncompleted_leaves[i - shift] = uncompleted_leaves[i];
    }

    uncompleted_leaves.resize(count - shift);
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::set_backtrack_path_for_leaf(size_n_max origin_leaf, size_n_max next_leaf) {
    size_n_max last_size = backtracks.size();

    size_n_max actual = next_leaf;
    while (actual != origin_leaf){
        backtracks.push_back(actual);
        actual = previous[actual];
    }
    backtrack_indexes[origin_leaf] = backtracks.size() - 1;

    if constexpr (COMPLEMENTS){
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

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline size_n_max LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::print_result(std::ostream& of) {
    if (!RESULT_COMPUTED){
        throw std::invalid_argument("Result has not been computed yet.");
    }

    size_n_max total_objective_value = 0;

    size_n_max actual = INVALID_NODE;
    kmer_t last_kmer = 0;
    bool first = true;
    for (size_n_max i = 0; i < N; ++i){
        if (components.find(i) != i){
            continue;
        }
        if constexpr (COMPLEMENTS){
            /// Prevent complementary chain from being printed later
            components.connect(i, components.find(complements[i]));
        }

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
            if constexpr (OBJECTIVE == JointObjective::RUNS) if (ov < K - 1) total_objective_value += PENALTY;
            if constexpr (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1 - ov);

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
                    if constexpr (OBJECTIVE == JointObjective::RUNS) if (ov < K - 1) total_objective_value += PENALTY;
                    if constexpr (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1 - ov);

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
    if constexpr (OBJECTIVE == JointObjective::RUNS) total_objective_value += PENALTY; /// First run of ones
    if constexpr (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1); /// Zeros at the end

    return total_objective_value;
}
