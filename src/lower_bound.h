#pragma once

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

template <typename kmer_t>
size_t LowerBoundMatchtigCount(std::vector<kmer_t>&& kMerVec, size_k_max k, bool complements) {
    typedef size_t size_n_max;

    if (complements) AddComplements(kMerVec, k);
    WriteLog("Finished adding complements.");
    
    std::sort(kMerVec.begin(), kMerVec.end());
    WriteLog("Finished sorting k-mers.");
    
    size_t N = kMerVec.size();
    
    FailureIndex<kmer_t, size_n_max> index(kMerVec, k); 
    auto INVALID_NODE = index.INVALID_NODE;
    UnionFind contracted_nodes(N);

    WriteLog("0");
    /// Contract cycles
    {
        std::vector<size_n_max> covered(N, INVALID_NODE);
        std::vector<size_n_max> previous(N, INVALID_NODE);

        for (size_n_max base_node = 0; base_node < N; ++base_node){
            if (covered[base_node] != INVALID_NODE) continue;
    
            std::vector<size_n_max> stack; stack.push_back(base_node);
            previous[base_node] = base_node;
            while (!stack.empty()){
                size_n_max node_index = stack.back(); stack.pop_back();
                if (covered[node_index] != INVALID_NODE) continue;
                covered[node_index] = base_node;
    
                size_n_max failure_index = index.find_first_failure_leaf_by_index(node_index, k - 1);
                if (failure_index == INVALID_NODE) continue;
                for (size_n_max i = 0; i < 4; ++i){
                    size_n_max next_index = failure_index + i;
                    if (BitSuffix(kMerVec[node_index], k - 1) != BitPrefix(kMerVec[next_index], k, k - 1)) break;
                    
                    if (covered[next_index] == INVALID_NODE){
                        stack.push_back(next_index);
                        previous[next_index] = node_index;
                    }
                    else if (covered[next_index] == base_node){
                        size_n_max backtrack = node_index;
                        while (!contracted_nodes.are_connected(backtrack, next_index)){
                            contracted_nodes.connect(next_index, backtrack);
                            // if (previous[backtrack] == backtrack) break;
                            backtrack = previous[backtrack];
                        }
                    }
                }
            }
        }
    }
    WriteLog("1");
    /// Count in-edges and construct out-edges of all contracted nodes
    std::vector<size_n_max> in_edge_counts(N, 0);
    std::vector<std::vector<size_n_max>> out_edges(N);

    for (size_n_max node_index = 0; node_index < N; ++node_index){
        size_n_max failure_index = index.find_first_failure_leaf_by_index(node_index, k - 1);
        if (failure_index == INVALID_NODE) continue;

        size_n_max contracted_index = contracted_nodes.find(node_index);
        for (size_n_max i = 0; i < 4; ++i){
            size_n_max next_index = failure_index + i;
            if (BitSuffix(kMerVec[node_index], k - 1) != BitPrefix(kMerVec[next_index], k, k - 1)) break;

            size_n_max contracted_next = contracted_nodes.find(next_index);
            if (contracted_next != contracted_index){
                ++in_edge_counts[contracted_next];
                out_edges[contracted_index].push_back(contracted_next);
            }
        }
    }
    WriteLog("2");
    /// Find all beginnings
    std::vector<size_n_max> beginnings;

    for (size_n_max node_index = 0; node_index < N; ++node_index){
        if (in_edge_counts[node_index] == 0 && contracted_nodes.find(node_index) == node_index){ /// Only one node per cycle
            beginnings.push_back(contracted_nodes.find(node_index));
        }
    }
    WriteLog("3");
    /// Count matchtigs by...
    size_n_max matchtigs = 0;
    for (size_n_max i = 0; i < beginnings.size(); ++i){
        size_n_max base_node = beginnings[i];
        if (in_edge_counts[base_node] == INVALID_NODE) continue; /// Node has been marked as processed
        
        ++matchtigs;

        std::vector<size_n_max> path; path.push_back(base_node);
        for (size_n_max p = 0; p < path.size(); ++p){
            size_n_max node_index = path[p];
            in_edge_counts[node_index] = INVALID_NODE;
            
            while (!out_edges[node_index].empty()){
                size_n_max next_index = out_edges[node_index].back(); out_edges.pop_back();
                if (in_edge_counts[next_index] != INVALID_NODE){
                    --in_edge_counts[next_index];
                    path.push_back(next_index);
                    break;
                }
            }
        }
        for (size_n_max node_index : path){
            while (!out_edges[node_index].empty() && in_edge_counts[out_edges[node_index].back()] == INVALID_NODE){
                out_edges[node_index].pop_back();
            }
            for (size_n_max j = 0; j < out_edges[node_index].size(); ++j){
                auto current_out_edges = out_edges[node_index];
                if (in_edge_counts[current_out_edges[j]] == INVALID_NODE){
                    std::swap(current_out_edges[j], current_out_edges.back());
                    current_out_edges.pop_back();
                } else {
                    --in_edge_counts[current_out_edges[j]];
                    if (in_edge_counts[current_out_edges[j]] == 0) beginnings.push_back(current_out_edges[j]);
                }
            }
        }
    }

    return matchtigs;
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
                if (objective == JointObjective::RUNS)       res += penalty;
                else if (objective == JointObjective::ZEROS) res += penalty * (j - 1);
                break;
            }
        }
    }
    return res / (1 + complements);
}

template <typename kmer_t>
size_t LowerBoundJoint(std::vector<kmer_t>&& kMerVec, const int k, bool complements, std::string objective_string, int penalty) {
    try {
        if (kMerVec.empty()) {
            throw std::invalid_argument("Empty input provided");
        }

        JointObjective objective = GetJointObjective(objective_string);

        if (complements) AddComplements(kMerVec, k);
        WriteLog("Finished adding complements.");

        std::sort(kMerVec.begin(), kMerVec.end());
        WriteLog("Finished sorting k-mers.");

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
