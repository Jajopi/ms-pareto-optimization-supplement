#pragma once

#include <vector>
#include <algorithm>
#include <limits>

#include "kmers.h"
#include "parser.h"
#include "pareto/objective.h"
#include "pareto/loac.h"

/// Run the optimization with the given type of leaf indexes, return the objective value
template <typename kmer_t, typename size_n_max>
size_t RunParetoOptimization(const std::vector<kmer_t>& kMers, std::ostream& of, size_k_max k, bool complements,
                             ParetoObjective objective, size_k_max penalty){
    LeafOnlyAC<kmer_t, size_n_max> loac(kMers, size_n_max(k), complements, objective, penalty);
    loac.compute_result();
    return loac.print_result(of);
}

/// Get the masked superstring of the given k-mers using the pareto optimization heuristic.
/// The objective is |S| + penalty * X, where X is the number of runs of ones (runs) or of zeros (zeros) in the mask.
/// The k-mers have to be distinct (canonical if complements are used). Returns false if the computation failed.
template <typename kmer_t>
bool ParetoOptimization(std::vector<kmer_t>&& kMerVec, std::ostream& of, size_k_max k, bool complements, ParetoOptions options){
    try {
        if (kMerVec.empty()) throw std::invalid_argument("Empty input provided");
        ParetoObjective objective = GetParetoObjective(options.objective);
        if (options.penalty < 0 || options.penalty > std::numeric_limits<size_k_max>::max())
            throw std::invalid_argument("Penalty has to be between 0 and " + std::to_string(std::numeric_limits<size_k_max>::max()));
        size_k_max penalty = size_k_max(options.penalty);
        if (penalty == 0){
            penalty = (objective == ParetoObjective::RUNS) ? DEFAULT_PENALTY_RUNS : DEFAULT_PENALTY_ZEROS;
            WriteLog("Using default penalty: " + std::to_string(penalty) + ".");
        }

        if (complements) AddComplements(kMerVec, k);
        std::sort(kMerVec.begin(), kMerVec.end());

        size_t n = kMerVec.size(), value;
        if (n <= (size_t(1) << 15))      value = RunParetoOptimization<kmer_t, uint16_t>(kMerVec, of, k, complements, objective, penalty);
        else if (n <= (size_t(1) << 31)) value = RunParetoOptimization<kmer_t, uint32_t>(kMerVec, of, k, complements, objective, penalty);
        else                             value = RunParetoOptimization<kmer_t, uint64_t>(kMerVec, of, k, complements, objective, penalty);
        WriteLog("Finished pareto optimization of masked superstring, resulting objective value: " + std::to_string(value) + ".");
    }
    catch (const std::exception& e){
        WriteLog("Exception was thrown: " + std::string(e.what()) + ".");
        return false;
    }
    return true;
}
