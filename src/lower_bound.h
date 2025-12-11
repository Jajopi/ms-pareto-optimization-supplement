#pragma once

#include <vector>
#include <tuple>

#include "global.h"
#include "kmers.h"
#include "joint/objective.h"
#include "joint/index.h"
#include "joint/unionfind.h"

/// Return the length of the cycle cover which lower bounds the superstring length.
template <typename kmer_t, typename kh_wrapper_t>
size_t LowerBoundLength(kh_wrapper_t wrapper, kmer_t kmer_type, std::vector<simplitig_t> simplitigs, const int k, bool complements) {
    auto cycle_cover = OverlapHamiltonianPath(wrapper, kmer_type, simplitigs, k, complements, true);
    size_t res = 0;
    for (auto &simplitig : simplitigs) {
        res += simplitig.size() / (2 - complements);
    }
    for (auto &overlap : cycle_cover.second) {
        res -= overlap;
    }
    return res / (1 + complements);
}


/// Kuhn's algorithm from https://cp-algorithms.com/graph/kuhn_maximum_bipartite_matching.html
template <typename size_n_max>
class MaximumBipartiteMatchingSolver {
    constexpr static size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();
    std::vector<std::vector<size_n_max>> g;
    size_n_max n, k;
    std::vector<size_n_max> mt;
    std::vector<bool> used;
    
    bool try_kuhn(size_n_max v) {
        if (used[v])
        return false;
        used[v] = true;
        for (size_n_max to : g[v]) {
            if (mt[to] == INVALID_NODE || try_kuhn(mt[to])) {
                mt[to] = v;
                return true;
            }
        }
        return false;
    }
public:
    size_n_max result;

    MaximumBipartiteMatchingSolver(const std::vector<std::vector<size_n_max>> &edges)
    : g(edges), n(edges.size()), k(edges.size()), result(0) {
        mt.assign(k, INVALID_NODE);

        /// Get arbitrary matching to start with
        std::vector<bool> used1(n, false);
        for (size_n_max v = 0; v < n; ++v) {
            if (v % 10000 == 0) WriteLog("Computing arbitrary matching: " + std::to_string(v) + "/" + std::to_string(n) + ".");
            for (size_n_max to : g[v]) {
                if (mt[to] == INVALID_NODE) {
                    mt[to] = v;
                    used1[v] = true;
                    result++;
                    break;
                }
            }
        }
        /// Improve matching using Kuhn's algorithm
        for (size_n_max v = 0; v < n; ++v) {
            if (v % 10000 == 0) WriteLog("Improving matching: " + std::to_string(v) + "/" + std::to_string(n) + ".");
            if (used1[v])
                continue;
            used.assign(n, false);
            result += try_kuhn(v);
        }
        WriteLog("Finished computing maximum bipartite matching of size " + std::to_string(result) + ".");
    }
};

template <typename kmer_t, typename size_n_max>
size_t compute_matchtig_count_lower_bound(std::vector<kmer_t> kMers, size_k_max k, bool complement_mode){
    constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();
    size_n_max N = kMers.size();
    WriteLog("Nodes initially: " + std::to_string(N) + ".");

    std::vector<size_n_max> contracted_indexes(N, INVALID_NODE);
    size_n_max CN; // Contracted node count
    
    std::vector<std::vector<size_n_max>> edges;
    size_n_max matchtig_count = 0;

    {
        /// Contract cycles and count in-edges -- Tarjan's algorithm for finding strongly connected components
        FailureIndex<kmer_t, size_n_max> failureIndex(kMers, k, false);
        UnionFind contracted_nodes(N);
        std::vector<size_k_max> in_edges(N, 0);
        {
            std::vector<std::tuple<size_n_max, size_n_max, size_n_max>> stack; stack.reserve(N); // Node index, already visited edges, parent node index
            std::vector<size_n_max> current_path; current_path.reserve(N);
            std::vector<bool> in_current_path(N, false);
            size_n_max next_index = 0;
            std::vector<size_n_max> indexes(N, INVALID_NODE);
            std::vector<size_n_max> lowlinks(N, INVALID_NODE);

            for (size_n_max base_node = 0; base_node < N; ++base_node){
                if (indexes[base_node] != INVALID_NODE) continue;
                stack.push_back(std::make_tuple(base_node, 0, INVALID_NODE));

                while (!stack.empty()){
                    auto p = stack.back(); stack.pop_back();
                    size_n_max node = std::get<0>(p), edge_index = std::get<1>(p), parent = std::get<2>(p);
                    if (edge_index == 0){ // Visiting for the first time
                        indexes[node] = lowlinks[node] = next_index++;
                        in_current_path[node] = true;
                        current_path.push_back(node);
                    }
                    stack.push_back(std::make_tuple(node, edge_index + 1, parent));

                    bool has_new_edge = true;
                    size_n_max failure = failureIndex.find_first_failure_leaf_by_index(node, k - 1);
                    if (failure == INVALID_NODE) has_new_edge = false; // No node to continue the search
                    size_n_max next_node = failure + edge_index;
                    if (next_node >= N || BitSuffix(kMers[node], k - 1) != BitPrefix(kMers[next_node], k, k - 1)) has_new_edge = false; // No node to continue the search

                    if (has_new_edge){
                        if (next_node != node){
                            ++in_edges[next_node];
                            if (indexes[next_node] == INVALID_NODE){
                                stack.push_back(std::make_tuple(next_node, 0, node)); // Not visited yet
                            } else if (in_current_path[next_node]){
                                if (lowlinks[node] > indexes[next_node]) lowlinks[node] = indexes[next_node]; // Back edge
                            }
                        }
                    } else {
                        stack.pop_back(); // Finished all edges from node
                        if (lowlinks[node] == indexes[node]){ // Root of an SCC
                            size_n_max contracted_count = 0;
                            size_n_max backtrack_node;
                            do {
                                backtrack_node = current_path.back(); current_path.pop_back();
                                in_current_path[backtrack_node] = false;
                                contracted_nodes.connect(node, backtrack_node);
                                ++contracted_count;
                            } while (backtrack_node != node);
                        }

                        if (parent != INVALID_NODE && lowlinks[parent] > lowlinks[node]) lowlinks[parent] = lowlinks[node];
                    }
                }
            }
        }
        WriteLog("Nodes after cycle contraction: " + std::to_string(contracted_nodes.count()));
        
        /// Contract paths
        for (size_n_max base_node = 0; base_node < N; ++base_node){
            size_n_max failure_index = failureIndex.find_first_failure_leaf_by_index(base_node, k - 1);
            if (failure_index == INVALID_NODE) continue; // No out-edges
            if (failure_index + 1 < N && BitSuffix(kMers[base_node], k - 1) == BitPrefix(kMers[failure_index + 1], k, k - 1)) continue; // Two or more out-edges

            if (in_edges[failure_index] > 1) continue; // Two or more in-edges
            contracted_nodes.connect(failure_index, base_node); // base-node -> failure_index
        }
        WriteLog("Nodes after path contraction: " + std::to_string(contracted_nodes.count()));

        {
            /// Remove single nodes and simple paths -- not really effective now, can be possibly done better (in later phases)
            std::vector<bool> was_removed(N, false);
            for (size_n_max base_node = 0; base_node < N; ++base_node){
                if (in_edges[base_node] != 0) continue;
                size_n_max last_node = contracted_nodes.find(base_node);
                if (failureIndex.find_first_failure_leaf_by_index(last_node, k - 1) != INVALID_NODE) continue;

                was_removed[last_node] = true;
                ++matchtig_count;
            }
            WriteLog("Nodes after non-branching path removal: " + std::to_string(contracted_nodes.count() - matchtig_count) + ".");
            
            CN = contracted_nodes.count() - matchtig_count; // matchtig_count are nodes removed in the previous step
            
            /// Link indices between contracted and base nodes
            std::vector<size_n_max> contracted_back_indices; contracted_back_indices.reserve(CN);
            contracted_back_indices.reserve(CN);
            for (size_n_max base_node = 0; base_node < N; ++base_node){
                size_n_max contracted_node = contracted_nodes.find(base_node);
                if (was_removed[contracted_node]) continue;

                if (contracted_indexes[contracted_node] == INVALID_NODE){
                    contracted_indexes[contracted_node] = contracted_back_indices.size();
                    contracted_back_indices.push_back(contracted_node);
                }
                contracted_indexes[base_node] = contracted_indexes[contracted_node];
            }
        }

        edges.resize(CN);
        for (size_n_max base_node = 0; base_node < N; ++base_node){
            size_n_max contracted_base = contracted_indexes[base_node];
            
            size_n_max failure_index = failureIndex.find_first_failure_leaf_by_index(base_node, k - 1);
            if (failure_index == INVALID_NODE) continue;
            for (size_n_max i = 0; i < 4; ++i){
                size_n_max next_index = failure_index + i;
                if (next_index >= N || BitSuffix(kMers[base_node], k - 1) != BitPrefix(kMers[next_index], k, k - 1)) break;
                
                size_n_max contracted_next = contracted_indexes[next_index];
                if (contracted_next != contracted_base) edges[contracted_base].push_back(contracted_next);
            }
        }
        for (size_n_max i = 0; i < CN; ++i){
            std::sort(edges[i].begin(), edges[i].end());
            edges[i].erase(std::unique(edges[i].begin(), edges[i].end()), edges[i].end());
        }
    }

    {
        /// Extend edges to reachability-edges
        std::vector<bool> has_children_complete(CN, false), is_in_stack(CN, false);
        std::vector<size_n_max> stack; stack.reserve(CN);
        for (size_n_max node = 0; node < CN; ++node){
            if (node % 10000 == 0) WriteLog("Extending reachability edges: " + std::to_string(node) + "/" + std::to_string(CN) + ".");

            if (has_children_complete[node]) continue;
            stack.clear(); stack.push_back(node);
            while (!stack.empty()){
                size_n_max reached_node = stack.back();
    
                if (!has_children_complete[reached_node]){
                    for (size_n_max next_node : edges[reached_node]){
                        // if (reached_node == next_node) continue;
                        if (is_in_stack[next_node]) continue;
                        is_in_stack[next_node] = true;
                        if (!has_children_complete[next_node]) stack.push_back(next_node);
                    }
                    has_children_complete[reached_node] = true;
                } else {
                    stack.pop_back();
                    size_n_max edge_count = edges[reached_node].size();
                    for (size_n_max e = 0; e < edge_count; ++e){
                        for (size_n_max further_node : edges[edges[reached_node][e]]){
                            edges[reached_node].push_back(further_node);
                        }
                    }
                    std::sort(edges[reached_node].begin(), edges[reached_node].end());
                    edges[reached_node].erase(std::unique(edges[reached_node].begin(), edges[reached_node].end()), edges[reached_node].end());
                }
            }
        }
    }

    /// Link complements and extend reachability in bi-directional mode -- I honestly have no idea what should happen with complements actually, have to think about it for a bit
    // if (complement_mode){
    //     std::vector<size_n_max> complements;
    //     std::vector<std::pair<kmer_t, size_n_max>> complement_kmers(N);
    //     for (size_n_max i = 0; i < N; ++i){
    //         complement_kmers[i] = std::make_pair(ReverseComplement(kMers[i], k), N - i); /// For even k, swapping indexes for same pairs of kmers fixes self-complement handling
    //     }

    //     std::sort(complement_kmers.begin(), complement_kmers.end());

    //     complements.resize(N);
    //     for (size_n_max i = 0; i < N; ++i) complements[i] = N - complement_kmers[i].second;

    //     std::vector<bool> complements_linked(CN, false);

    //     for (size_n_max base_node = 0; base_node < N; ++base_node){
    //         size_n_max contracted_base = contracted_indexes[base_node];
    //         if (complements_linked[contracted_base]) continue;

    //         size_n_max complement_node = complements[base_node];
    //         size_n_max contracted_complement = contracted_indexes[complement_node];
    //         if (contracted_base == contracted_complement) continue;
            
    //         edges[contracted_base][contracted_complement] = true;
    //         for (size_n_max next_node = 0; next_node < CN; ++next_node){
    //             if (edges[contracted_base][next_node]){
                    
    //             }
    //         }
    //         complements_linked[contracted_base] = true;
    //     }
    // }

    if (edges.size() <= 100){
        WriteLog("Final graph:");
        for (size_n_max i = 0; i < CN; ++i){
            std::cerr << "\t" << i << ": ";
            for (auto edge: edges[i]) std::cerr << edge << " ";
            std::cerr << std::endl;
        }
    }

    /// Compute maximum bipartite matching
    WriteLog("Starting maximum bipartite matching on " + std::to_string(CN) + " nodes.");
    
    auto mbms = MaximumBipartiteMatchingSolver(edges);

    return matchtig_count + (CN - mbms.result);
}

template <typename kmer_t>
size_t LowerBoundMatchtigCount(std::vector<kmer_t>&& kMerVec, size_k_max k, bool complements) {
    try {
        if (kMerVec.empty()) {
            throw std::invalid_argument("Empty input provided");
        }

        /// Add complements
        if (complements) AddComplements(kMerVec, k);

        /// Sort kmers
        std::sort(kMerVec.begin(), kMerVec.end());

        /// Remove k-mers present more times than they should be (2 for self complements, 1 otherwise)
        /// Skipping this step as input data are nice
        // RemoveDuplicateKmers(kMerVec, k % 2 == 0);

        size_t limit = kMerVec.size();
        if      (limit <= (size_t(1) << 15))
            return compute_matchtig_count_lower_bound<kmer_t, uint16_t>(kMerVec, k, complements);
        else if (limit <= (size_t(1) << 31))
            return compute_matchtig_count_lower_bound<kmer_t, uint32_t>(kMerVec, k, complements);
        else
            return compute_matchtig_count_lower_bound<kmer_t, uint64_t>(kMerVec, k, complements);
    }
    catch (const std::exception& e){
        WriteLog("Exception was thrown: " + std::string(e.what()) + ".");
        return 0;
    }
}

template <typename kmer_t, typename size_n_max>
size_t compute_joint_lower_bound(std::vector<kmer_t>& kMerVec, size_k_max k, bool complements, JointObjective objective, size_k_max penalty){
    size_t res = 0;
    size_t N = kMerVec.size();
    FailureIndex<kmer_t, size_t> failure_index(kMerVec, k);

    for (size_t i = 0; i < N; ++i){
        if (failure_index.failure_node_exists(i, k - 1)){
            ++res;
            continue;
        }
        for (size_k_max j = 2; j <= k; ++j){
            if (failure_index.failure_node_exists(i, k - j)){
                res += j;
                if (objective == JointObjective::RUNS)  res += penalty;
                if (objective == JointObjective::ZEROS) res += penalty * (j - 1);
                break;
            }
        }
    }
    return res / (1 + complements);
}

/// This lowerbound is quite simple and quite bad, it's not useful at all now
template <typename kmer_t>
size_t LowerBoundJoint(std::vector<kmer_t>&& kMerVec, const int k, bool complements, std::string objective_string, int penalty) {
    try {
        if (kMerVec.empty()) {
            throw std::invalid_argument("Empty input provided");
        }

        JointObjective objective = GetJointObjective(objective_string);

        /// Add complements
        if (complements) AddComplements(kMerVec, k);

        /// Sort kmers
        std::sort(kMerVec.begin(), kMerVec.end());

        /// Remove k-mers present more times than they should be (2 for self complements, 1 otherwise)
        /// Skipping this step as input data are nice
        // RemoveDuplicateKmers(kMerVec, k % 2 == 0);

        if (penalty == 0){
            if (objective == JointObjective::RUNS)  penalty = DEFAULT_PENALTY_RUNS;
            if (objective == JointObjective::ZEROS) penalty = DEFAULT_PENALTY_ZEROS;
            WriteLog("Using default penalty: " + std::to_string(penalty) + ".");
        }

        size_t limit = kMerVec.size();
        if      (limit <= (size_t(1) << 15))
            return compute_joint_lower_bound<kmer_t, uint16_t>(kMerVec, k, complements, objective, penalty);
        else if (limit <= (size_t(1) << 31))
            return compute_joint_lower_bound<kmer_t, uint32_t>(kMerVec, k, complements, objective, penalty);
        else
            return compute_joint_lower_bound<kmer_t, uint64_t>(kMerVec, k, complements, objective, penalty);
    }
    catch (const std::exception& e){
        WriteLog("Exception was thrown: " + std::string(e.what()) + ".");
        return 0;
    }
}
