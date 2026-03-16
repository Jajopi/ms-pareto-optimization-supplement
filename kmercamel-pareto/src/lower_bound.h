#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <limits>

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


/// Kuhn's algorithm from https://cp-algorithms.com/graph/kuhn_maximum_bipartite_matching.html (but we ignore different partition sizes)
/// Speedup:
/// - Uses given node order (reverse topological)
/// - Both heuristic and improvement phases use recursion to supply reachability edges
/// As a result, it runs in few seconds instead a few hours on covid dataset
class MaximumBipartiteMatchingSolver {
    const std::vector<std::vector<uint32_t>>& g;
    const std::vector<uint32_t>& node_order;
    uint32_t n;
    std::vector<uint32_t> mt;
    std::vector<bool> used;
    std::vector<bool> visited;

    inline bool try_kuhn(uint32_t from, uint32_t target) {
        if (visited[from]) return false;
        visited[from] = true;

        if (used[from]) return false;
        used[from] = true;

        for (uint32_t to : g[from]) {
            if (node_order[to] == INVALID_NODE){ // Skip marked nodes
                if (try_kuhn(to, target)) return true;
                continue;
            }

            if (mt[to] == INVALID_NODE || try_kuhn(mt[to], mt[to])) {
                mt[to] = target;
                return true;
            }
            if (try_kuhn(to, target)) return true;
        }
        return false;
    }
    inline bool try_heuristic_reachable_edges(uint32_t from, uint32_t target){
        if (visited[from]) return false;
        visited[from] = true;

        for (uint32_t to : g[from]) {
            if (node_order[to] == INVALID_NODE){ // Skip marked nodes
                if (try_heuristic_reachable_edges(to, target)) return true;
                continue;
            }

            if (mt[to] == INVALID_NODE) {
                mt[to] = target;
                return true;
            }
            if (try_heuristic_reachable_edges(to, target)) return true;
        }
        return false;
    }
public:
    constexpr static uint32_t INVALID_NODE = std::numeric_limits<uint32_t>::max();
    uint32_t result;

    MaximumBipartiteMatchingSolver(const std::vector<std::vector<uint32_t>> &edges, const std::vector<uint32_t> &node_order, uint32_t removed_nodes = 0)
    : g(edges), node_order(node_order), n(node_order.size()), mt(node_order.size(), INVALID_NODE), used(node_order.size(), false), result(0) {

        /// Get arbitrary matching to start with
        WriteLog("Computing arbitrary matching of graph with " + std::to_string(n - removed_nodes) + " nodes by heuristic...");
        std::vector<bool> used_by_heuristic(n, false);
        visited.assign(n, false);
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t v = node_order[i];
            if (v == INVALID_NODE) continue; // Node was removed
            if (try_heuristic_reachable_edges(v, v)){
                used_by_heuristic[v] = true;
                result++;
            }
            visited.assign(n, false);
        }
        /// Improve matching using Kuhn's algorithm
        WriteLog("Heuristic matched " + std::to_string(result) + ". Improving using Kuhn's algorithm...");
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t v = node_order[i];
            if (v == INVALID_NODE) continue; // Node was removed
            if (used_by_heuristic[v]) continue;

            used.assign(n, false);
            if (try_kuhn(v, v)) result++;
            visited.assign(n, false);
        }
        WriteLog("Finished computing maximum bipartite matching of size " + std::to_string(result) + ".");
    }
};

template <typename kmer_t, typename size_n_max>
size_t compute_matchtig_count_lower_bound(std::vector<kmer_t> kMers, size_k_max k, bool complement_mode){
    constexpr const size_n_max INVALID_NODE = std::numeric_limits<size_n_max>::max();
    size_n_max N = kMers.size();
    WriteLog("Nodes initially: " + std::to_string(N) + ". Contracting cycles...");

    std::vector<size_n_max> contracted_indexes(N, INVALID_NODE);
    size_n_max CN; // Contracted node count

    std::vector<std::vector<uint32_t>> edges;
    uint32_t matchtig_count = 0;

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
        WriteLog("Nodes after cycle contraction: " + std::to_string(contracted_nodes.count()) + ". Contracting paths...");

        /// Contract paths
        for (size_n_max base_node = 0; base_node < N; ++base_node){
            size_n_max failure_index = failureIndex.find_first_failure_leaf_by_index(base_node, k - 1);
            if (failure_index == INVALID_NODE) continue; // No out-edges
            if (failure_index + 1 < N && BitSuffix(kMers[base_node], k - 1) == BitPrefix(kMers[failure_index + 1], k, k - 1)) continue; // Two or more out-edges

            if (in_edges[failure_index] > 1) continue; // Two or more in-edges
            contracted_nodes.connect(failure_index, base_node); // base-node -> failure_index
        }
        WriteLog("Nodes after path contraction: " + std::to_string(contracted_nodes.count()) + ". Removing single nodes and simple paths...");

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
            WriteLog("Nodes after non-branching path removal: " + std::to_string(contracted_nodes.count() - matchtig_count) + ". Building graph...");

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

    if (edges.size() <= 100){
        WriteLog("Final graph:");
        for (size_n_max i = 0; i < CN; ++i){
            std::cerr << "\t" << i << ": ";
            for (auto edge: edges[i]) std::cerr << edge << " ";
            std::cerr << std::endl;
        }
    }

    std::vector<uint32_t> reverse_topological_order(CN, 0);
    {
        /// Mark nodes reverse topologically -- children of node have lower numbers than the node
        WriteLog("Sorting nodes and edges by reverse topological order...");
        std::vector<bool> has_children_complete(CN, false), is_in_stack(CN, false);
        std::vector<uint32_t> sorted_orders(CN, 0);
        std::vector<uint32_t> stack; stack.reserve(CN);
        for (uint32_t node = 0; node < CN; ++node){
            if (has_children_complete[node]) continue;

            stack.clear(); stack.push_back(node);
            while (!stack.empty()){
                uint32_t reached_node = stack.back();

                if (!has_children_complete[reached_node]){
                    for (uint32_t next_node : edges[reached_node]){
                        if (is_in_stack[next_node]) continue;
                        is_in_stack[next_node] = true;
                        if (!has_children_complete[next_node]) stack.push_back(next_node);
                    }
                    has_children_complete[reached_node] = true;
                } else {
                    stack.pop_back();
                    uint32_t edge_count = edges[reached_node].size();
                    for (uint32_t e = 0; e < edge_count; ++e){
                        sorted_orders[reached_node] = std::max(sorted_orders[reached_node], sorted_orders[edges[reached_node][e]] + 1);
                    }
                }
            }
        }
        /// Get node order
        std::iota(reverse_topological_order.begin(), reverse_topological_order.end(), 0);
        std::sort(reverse_topological_order.begin(), reverse_topological_order.end(), [&](uint32_t a, uint32_t b){
            return sorted_orders[a] < sorted_orders[b];
        });

        /// Sort edges
        for (uint32_t node = 0; node < CN; ++node){
            std::sort(edges[node].begin(), edges[node].end(), [&](uint32_t a, uint32_t b){
                return sorted_orders[a] < sorted_orders[b];
            });
        }
    }

    /// Compute maximum bipartite matching
    {
        auto mbms_all = MaximumBipartiteMatchingSolver(edges, reverse_topological_order);
        matchtig_count += CN - mbms_all.result; // Number of matchtigs is number of unmatched nodes
    }

    if (!complement_mode) return matchtig_count;

    /// Link complements in bi-directional mode, remove self-complement nodes
    WriteLog("Searching for complements and marking self-complement nodes...");
    size_n_max self_complement_count = 0; // Trimmed node count
    {
        std::vector<size_n_max> complements;
        std::vector<std::pair<kmer_t, size_n_max>> complement_kmers(N);
        for (size_n_max i = 0; i < N; ++i){
            complement_kmers[i] = std::make_pair(ReverseComplement(kMers[i], k), N - i); /// For even k, swapping indexes for same pairs of kmers fixes self-complement handling
        }
        std::sort(complement_kmers.begin(), complement_kmers.end());
        complements.resize(N);
        for (size_n_max i = 0; i < N; ++i) complements[i] = N - complement_kmers[i].second;

        for (size_n_max base_node = 0; base_node < N; ++base_node){
            size_n_max contracted_base = contracted_indexes[base_node];
            if (contracted_base == INVALID_NODE) continue;
            if (contracted_base == contracted_indexes[complements[base_node]]){
                if (reverse_topological_order[contracted_base] != MaximumBipartiteMatchingSolver::INVALID_NODE) self_complement_count++;
                reverse_topological_order[contracted_base] = MaximumBipartiteMatchingSolver::INVALID_NODE;
            }
        }

        WriteLog("Found " + std::to_string(self_complement_count) + " self-complement nodes.");
    }

    /// Compute maximum bipartite matching on graph without self-complement nodes to get duplicate matchtigs
    auto mbms_trimmed = MaximumBipartiteMatchingSolver(edges, reverse_topological_order, self_complement_count);
    matchtig_count -= (CN - self_complement_count - mbms_trimmed.result) / 2;  // Number of matchtigs is number of unmatched nodes

    return matchtig_count;
}

template <typename kmer_t>
size_t LowerBoundMatchtigCount(std::vector<kmer_t>&& kMerVec, size_k_max k, bool complements) {
    try {
        if (kMerVec.empty()) {
            throw std::invalid_argument("Empty input provided.");
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
    FailureIndex<kmer_t, size_n_max> failure_index(kMerVec, k);

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
