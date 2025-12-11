#pragma once

#include <vector>

#include "objective.h"

/// Index for fast retrieval of first failure node of given leaf and depth present in the set
/// Search operation returns std::numeric_limits<size_n_max>::max() when k-mer is not present in the set
template<typename kmer_t, typename size_n_max>
class FailureIndex {
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
        for (size_n_max i = 0; i < N; ++i) first_row[i] = search(kMers[i], K - 1);
    }
    
    inline size_n_max search(kmer_t searched, size_k_max depth){
        searched = BitSuffix(searched, depth);

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
    static constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();
    
    inline FailureIndex(const std::vector<kmer_t> &kmers, size_k_max k, bool print_progress = true) :
    kMers(kmers), N(kmers.size()), K(k) {
        construct_index();
        if (print_progress) WriteLog("Finished constructing index.");
    }

    inline size_n_max find_first_failure_leaf(kmer_t kmer, size_k_max depth){
        return search(kmer, depth);
    } 

    inline size_n_max find_first_failure_leaf_by_index(size_n_max index, size_k_max depth){
        if (depth == K - 1) return first_row[index];
        return search(kMers[index], depth);
    }

    inline bool failure_node_exists(size_n_max index, size_k_max depth){
        return find_first_failure_leaf_by_index(index, depth) != INVALID_NODE;
    }
};
