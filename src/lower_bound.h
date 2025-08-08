#pragma once

#include "global.h"
#include "kmers.h"

#include "joint.h"

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
size_t LowerBoundMatchtigsCount(std::vector<kmer_t>& kMerVec, int k, bool complements) {
    
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

template <typename kmer_t>
size_t LowerBoundJoint(std::vector<kmer_t>& kMerVec, int k, bool complements, std::string objective_string) {
    // JointObjective objective;
    // if (objective_string == "zeros") objective = JointObjective::ZEROS;
    // else if (objective_string == "runs") objective = JointObjective::RUNS;
    // else throw std::invalid_argument("Wrong objective provided");
    return 0;
}
