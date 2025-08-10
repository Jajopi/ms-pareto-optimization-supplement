#pragma once

#include "global.h"
#include "kmers.h"

#include "joint_objective.h"
#include "joint_index.h"

/// Return the length of the cycle cover which lower bounds the superstring length.
template <typename kmer_t, typename kh_wrapper_t>
size_t LowerBoundLength(kh_wrapper_t wrapper, kmer_t kmer_type, std::vector<simplitig_t> simplitigs, int k, bool complements) {
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
size_t LowerBoundMatchtigsCount(std::vector<kmer_t>&& kMerVec, int k, bool complements) {
    
    if (complements){
        /// Add complementary k-mers.
        size_t n = kMerVec.size();
        kMerVec.resize(n * 2);
        for (size_t i = 0; i < n; ++i) kMerVec[i + n] = ReverseComplement(kMerVec[i], k);
        WriteLog("Finshed adding complements.");
    }
    std::sort(kMerVec.begin(), kMerVec.end());
    
    FailureIndex<kmer_t, size_t> index(kMerVec, k);
    return 0;
}

template <typename kmer_t, typename size_n_max>
size_t compute_joint_lower_bound(std::vector<kmer_t>& kMers, size_k_max k, bool complements, JointObjective objective, size_k_max penalty){
    size_t res = 0;
    size_t N = kMers.size();
    FailureIndex<kmer_t, size_t> failure_index(kMers, k);

    for (size_t i = 0; i < N; ++i){
        if (failure_index.exists(i, k - 1)){
            ++res;
            continue;
        }
        for (size_k_max j = 2; j <= k; ++j){
            if (failure_index.exists(i, k - j)){
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
size_t LowerBoundJoint(std::vector<kmer_t>&& kMerVec, int k, bool complements, std::string objective_string) {
    try {
        if (kMerVec.empty()) {
            throw std::invalid_argument("Empty input provided");
        }

        JointObjective objective = GetJointObjective(objective_string);

        if (complements){
            /// Add complementary k-mers.
            size_t n = kMerVec.size();
            kMerVec.resize(n * 2);
            for (size_t i = 0; i < n; ++i) kMerVec[i + n] = ReverseComplement(kMerVec[i], k);
            WriteLog("Finished adding complements.");
        }

        /// Sort k-mers
        std::sort(kMerVec.begin(), kMerVec.end());
        WriteLog("Finished sorting k-mers.");

        /// Skipping this step as input data are nice
        // RemoveDuplicateKmers(kMerVec, k % 2 == 0);

        size_k_max penalty;
        if (objective == JointObjective::RUNS)       penalty = 7;
        else if (objective == JointObjective::ZEROS) penalty = 2;

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
