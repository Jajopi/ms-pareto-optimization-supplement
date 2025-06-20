#pragma once

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <cstdint>
#include <limits>

#include "kmers.h"
#include "parser.h"

enum JointObjective {
    RUNS, ZEROS
};

#include "joint_loac.h"

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
void compute_with_loac(LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS> loac, std::ostream& of){
    loac.compute_result();
    loac.optimize_result(); // TODO add parameters
    size_t objective = loac.print_result(of);

    std::cerr << objective << std::endl;
}

/// The joint optimization process
template <typename kmer_t, typename size_n_max>
void compute_joint_optimization(std::vector<kmer_t>& kMers, std::ostream& of, size_t k,
        bool complements, JointObjective objective){

    if (objective == JointObjective::RUNS){
        if (complements) compute_with_loac(LeafOnlyAC<kmer_t, size_n_max,  JointObjective::RUNS,  true>(kMers, size_n_max(k), 7), of);
        else             compute_with_loac(LeafOnlyAC<kmer_t, size_n_max,  JointObjective::RUNS, false>(kMers, size_n_max(k), 7), of);
    }
    else if (objective == JointObjective::ZEROS){
        if (complements) compute_with_loac(LeafOnlyAC<kmer_t, size_n_max, JointObjective::ZEROS,  true>(kMers, size_n_max(k), 7), of);
        else             compute_with_loac(LeafOnlyAC<kmer_t, size_n_max, JointObjective::ZEROS, false>(kMers, size_n_max(k), 7), of);
    }
    else {
        throw std::invalid_argument("Invalid objective function");
    }    
}

// Removes duplicate k-mers from sorted vector
template <typename kmer_t>
size_t RemoveDuplicateKmers(std::vector<kmer_t> & kMerVec, bool even_k){
    size_t size_limit = kMerVec.size() - (even_k ? 2 : 1);
    size_t removed = 0;
    for (size_t i = 0; i < size_limit; ++i){
        if (kMerVec[i] == kMerVec[i + 1 + even_k]) ++removed;
        if (removed != 0) kMerVec[i] = kMerVec[i + removed];
    }
    kMerVec.resize(kMerVec.size() - removed);
    return removed;
}

/// Get the masked superstring of the given k-mers using the joint optimization heuristic.
/// The objective is in the form |S| + c * X, where X is either:
/// - number of runs of ones in the mask    (JointObjective.RUNS)
/// - number of zeroes if the mask          (JointObjective.ZEROS)
///
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
template <typename kmer_t>
void JointOptimization(std::vector<kmer_t>& kMerVec, std::ostream& of, int k, bool complements, std::string objective_string){
    try {
        if (kMerVec.empty()) {
            throw std::invalid_argument("Empty input provided");
        }

        JointObjective objective;
        if (objective_string == "zeros") objective = JointObjective::ZEROS;
        else if (objective_string == "runs") objective = JointObjective::RUNS;
        else throw std::invalid_argument("Wrong objective provided");

        if (complements){
            /// Add complementary k-mers.
            size_t n = kMerVec.size();
            kMerVec.resize(n * 2);
            for (size_t i = 0; i < n; ++i) kMerVec[i + n] = ReverseComplement(kMerVec[i], k);
            WriteLog("Finshed adding complements.");
        }

        /// Sort k-mers
        std::sort(kMerVec.begin(), kMerVec.end());
        WriteLog("Finished sorting k-mers.");

        /// Remove k-mers present more times than the should be (2 for self complements, 1 otherwise)
        // RemoveDuplicateKmers(kMerVec, k % 2 == 0);

        size_t limit = kMerVec.size();
        if      (limit <= (size_t(1) << 15))
            compute_joint_optimization<kmer_t, uint16_t>(kMerVec, of, k, complements, objective);
        else if (limit <= (size_t(1) << 31))
            compute_joint_optimization<kmer_t, uint32_t>(kMerVec, of, k, complements, objective);
        else
            compute_joint_optimization<kmer_t, uint64_t>(kMerVec, of, k, complements, objective);
    }
    catch (const std::exception& e){
        WriteLog("Exception was thrown: " + std::string(e.what()));
    }
}
