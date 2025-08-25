#pragma once

#include <vector>
#include <algorithm>

#include "kmers.h"
#include "parser.h"
#include "joint/objective.h"
#include "joint/loac.h"

template <typename kmer_t, typename size_n_max, JointObjective OBJECTIVE, bool COMPLEMENTS>
void compute_with_loac(LeafOnlyAC<kmer_t, size_n_max, OBJECTIVE, COMPLEMENTS>&& loac, std::ostream& of){
    loac.compute_result();
    loac.optimize_result(); // TODO add parameters
    size_t objective = loac.print_result(of);

    WriteLog("Finished joint optimization of masked superstring, resulting objective: " + std::to_string(objective) + ".");
}

/// The joint optimization process
template <typename kmer_t, typename size_n_max>
void compute_joint_optimization(std::vector<kmer_t>& kMers, std::ostream& of, size_k_max k,
        bool complements, JointObjective objective){

    if (objective == JointObjective::RUNS){
        if (complements) compute_with_loac(LeafOnlyAC<kmer_t, size_n_max,  JointObjective::RUNS,  true>(kMers, size_n_max(k), DEFAULT_PENALTY_RUNS),  of);
        else             compute_with_loac(LeafOnlyAC<kmer_t, size_n_max,  JointObjective::RUNS, false>(kMers, size_n_max(k), DEFAULT_PENALTY_RUNS),  of);
    } else if (objective == JointObjective::ZEROS){
        if (complements) compute_with_loac(LeafOnlyAC<kmer_t, size_n_max, JointObjective::ZEROS,  true>(kMers, size_n_max(k), DEFAULT_PENALTY_ZEROS), of);
        else             compute_with_loac(LeafOnlyAC<kmer_t, size_n_max, JointObjective::ZEROS, false>(kMers, size_n_max(k), DEFAULT_PENALTY_ZEROS), of);
    } else {
        throw std::invalid_argument("Invalid objective function specified");
    }
}

/// Remove k-mers present more times than they should be (2 for self complements, 1 otherwise)
template <typename kmer_t>
size_t RemoveDuplicateKmers(std::vector<kmer_t>& sorted_kMerVec, bool even_k){
    size_t size_limit = sorted_kMerVec.size() - (even_k ? 2 : 1);
    size_t removed = 0;
    for (size_t i = 0; i < size_limit; ++i){
        if (sorted_kMerVec[i] == sorted_kMerVec[i + 1 + even_k]) ++removed;
        if (removed != 0) sorted_kMerVec[i] = sorted_kMerVec[i + removed];
    }
    sorted_kMerVec.resize(sorted_kMerVec.size() - removed);
    return removed;
}

/// Get the masked superstring of the given k-mers using the joint optimization heuristic.
/// The objective is in the form |S| + c * X, where X is either:
/// - number of runs of ones in the mask    (JointObjective.RUNS)
/// - number of zeroes if the mask          (JointObjective.ZEROS)
///
/// This runs in O(n^2 k^2), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
template <typename kmer_t>
void JointOptimization(std::vector<kmer_t>&& kMerVec, std::ostream& of, size_k_max k, bool complements, std::string objective_string){
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

        size_t limit = kMerVec.size();
        if      (limit <= (size_t(1) << 15))
            compute_joint_optimization<kmer_t, uint16_t>(kMerVec, of, k, complements, objective);
        else if (limit <= (size_t(1) << 31))
            compute_joint_optimization<kmer_t, uint32_t>(kMerVec, of, k, complements, objective);
        else
            compute_joint_optimization<kmer_t, uint64_t>(kMerVec, of, k, complements, objective);
    }
    catch (const std::exception& e){
        WriteLog("Exception was thrown: " + std::string(e.what()) + ".");
    }
}
