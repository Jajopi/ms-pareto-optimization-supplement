#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <random>

#include "../kmers.h"
#include "../joint/objective.h"
#include "../joint/index.h"
#include "../joint/unionfind.h"

constexpr size_t RANDOM_SEED = 0;
constexpr size_t MAX_COUNT_WIDTH = 12;
constexpr size_t MAX_ITERS_WIDTH = 3;

/// The data structure for efficient heuristic search
template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
class LeafOnlyAC {
    static constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();

    size_k_max K;               /// Kmer-length
    size_n_max N;               /// Number of k-mers (number of leaves)
    size_k_max PENALTY;         /// Value of penalty (for runs / zeros) used in the computation
    
    std::vector<kmer_t>& kMers;  /// Sorted leaves
    FailureIndex<kmer_t, size_n_max> failureIndex;
    UnionFind<size_n_max>   components;
    
    std::vector<size_n_max> complements;
    std::vector<size_n_max> backtracks;
    std::vector<size_n_max> backtrack_indexes;
    std::vector<size_n_max> previous;
    // std::vector<size_n_max> bridging_kmers;
    std::vector<size_n_max> next;
    std::vector<size_k_max> remaining_priorities;
    std::vector<size_n_max> skip_to;
    std::vector<bool>       used;
    std::vector<size_k_max> no_unused;
    std::vector<std::tuple<size_k_max, size_k_max, size_n_max, size_n_max>> stack;

    bool try_complete_leaf(size_n_max leaf_to_connect, size_k_max priority_drop_limit);
    void push_failure_node_into_stack(size_k_max priority, size_k_max node_depth, size_n_max node_index, size_n_max last_leaf);
    void squeeze_uncompleted_leaves(std::vector<size_n_max>& unclompleted_leaves);
    void set_backtrack_path_for_leaf(size_n_max origin_leaf, size_n_max next_leaf);

    bool RESULT_COMPUTED = false;
public:
    std::ostream& LOG_STREAM = std::cerr;

    LeafOnlyAC(std::vector<kmer_t>& kmers, size_k_max k, size_k_max penalty) :
        K(k), N(kmers.size()), PENALTY(penalty), kMers(kmers), failureIndex(kMers, K) {
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
    
        std::vector<std::pair<kmer_t, size_n_max>> complement_kmers(N);
        for (size_n_max i = 0; i < N; ++i){
            complement_kmers[i] = std::make_pair(ReverseComplement(kMers[i], K),
                                                 N - i); /// For even k, swapping indexes for same pairs of kmers is needed
        }
    
        std::sort(complement_kmers.begin(), complement_kmers.end());
        
        for (size_n_max i = 0; i < N; ++i) complements[i] = N - complement_kmers[i].second;
    }
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::compute_result() {
    if (RESULT_COMPUTED){
        throw std::invalid_argument("Result has already been computed.");
    }
    size_n_max SEARCH_CUTOFF = N / (1 << 18);

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
    
    size_k_max max_priority_drop;
    if constexpr (OBJECTIVE == JointObjective::RUNS)  max_priority_drop = (K - 1) + PENALTY;
    if constexpr (OBJECTIVE == JointObjective::ZEROS) max_priority_drop = (K - 1) * (PENALTY + 1);
    
    size_n_max remaining_iterations = max_priority_drop;
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations << std::endl;
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();

    size_k_max bigger_step;
    if constexpr (OBJECTIVE == JointObjective::RUNS)  bigger_step = 1;
    if constexpr (OBJECTIVE == JointObjective::ZEROS) bigger_step = PENALTY;
    for (size_k_max priority_drop_limit = 1;
        priority_drop_limit <= max_priority_drop;
        priority_drop_limit += (priority_drop_limit <= K) ? 1 : bigger_step){
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
            } else {
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
    WriteLog("Starting secondary optimization of the masked superstring.");
    WriteLog("No optimization implemented yet.");
    WriteLog("Finished secondary optimization of the masked superstring.");
}

// Internal functions

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline bool LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::try_complete_leaf(
        size_n_max leaf_to_complete, size_k_max priority_drop_limit) {

    if (priority_drop_limit == 1){
        size_n_max first_failure_leaf = failureIndex.find_first_failure_leaf_by_index(leaf_to_complete, K - 1);

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
        size_k_max priority    = std::get<0>(t);
        size_k_max chain_depth = std::get<1>(t);
        size_n_max leaf_index  = std::get<2>(t);
        size_n_max last_leaf   = std::get<3>(t);

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

        /// Append into stack
        if (priority == 0) continue;
        push_failure_node_into_stack(priority, chain_depth, last_leaf, last_leaf); /// Add failure of current node

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
            push_failure_node_into_stack(priority, K, i, i); /// Add failure of leaf i
        }
    }

    return false;
}

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
inline void LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>::push_failure_node_into_stack(
    size_k_max priority, size_k_max node_depth, size_n_max node_index, size_n_max last_leaf){

    size_k_max failure_depth = node_depth - 1;
    size_n_max failure_index = failureIndex.find_first_failure_leaf_by_index(node_index, failure_depth);
    if (--priority == 0) return;

    while (failure_index == INVALID_NODE){
        if (--failure_depth == 0) return;
        if (--priority == 0) return;

        failure_index = failureIndex.find_first_failure_leaf_by_index(node_index, failure_depth);
    }

    if constexpr (OBJECTIVE == JointObjective::RUNS){
        if (node_depth == K - 1 || (node_depth == K && failure_depth < K - 1)){ /// Run will be interrupted
            if (priority < PENALTY) return;
            priority -= PENALTY;
        }
    }
    if constexpr (OBJECTIVE == JointObjective::ZEROS){
        if (priority < (node_depth - failure_depth - (node_depth == K)) * PENALTY) return;
        priority -= (node_depth - failure_depth - (node_depth == K)) * PENALTY;
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
            if constexpr (OBJECTIVE == JointObjective::RUNS)  if (ov < K - 1) total_objective_value += PENALTY;
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
                    if constexpr (OBJECTIVE == JointObjective::RUNS)  if (ov < K - 1) total_objective_value += PENALTY;
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
    if constexpr (OBJECTIVE == JointObjective::RUNS)  total_objective_value += PENALTY;           /// First run of ones
    if constexpr (OBJECTIVE == JointObjective::ZEROS) total_objective_value += PENALTY * (K - 1); /// Zeros at the end

    return total_objective_value;
}
