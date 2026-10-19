#pragma once

#include <vector>
#include <algorithm>
#include <limits>

#include "kmers.h"
#include "parser.h"
#include "pareto/objective.h"
#include "pareto/loac.h"

/// Run the optimization with the given types of leaf indexes and stored remaining priorities, return the objective value
template <typename kmer_t, typename size_n_max, typename priority_t>
size_t RunParetoOptimization(const std::vector<kmer_t>& kMers, std::ostream& of, size_k_max k, bool complements,
                             ParetoObjective objective, size_k_max penalty, const ParetoOptions &options){
    LeafOnlyAC<kmer_t, size_n_max, priority_t> loac(kMers, k, complements, objective, penalty, options);
    loac.compute_result();
    return loac.print_result(of);
}

template <typename kmer_t, typename size_n_max>
size_t RunParetoOptimization(const std::vector<kmer_t>& kMers, std::ostream& of, size_k_max k, bool complements,
                             ParetoObjective objective, size_k_max penalty, const ParetoOptions &options){
    /// Remaining priorities fit in a byte for usual penalties, larger budgets need two bytes
    if (LeafOnlyAC<kmer_t, size_n_max>::max_budget(k, objective, penalty) <= std::numeric_limits<uint8_t>::max())
        return RunParetoOptimization<kmer_t, size_n_max, uint8_t>(kMers, of, k, complements, objective, penalty, options);
    return RunParetoOptimization<kmer_t, size_n_max, uint16_t>(kMers, of, k, complements, objective, penalty, options);
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
        if (n <= (size_t(1) << 15))      value = RunParetoOptimization<kmer_t, uint16_t>(kMerVec, of, k, complements, objective, penalty, options);
        else if (n <= (size_t(1) << 31)) value = RunParetoOptimization<kmer_t, uint32_t>(kMerVec, of, k, complements, objective, penalty, options);
        else                             value = RunParetoOptimization<kmer_t, uint64_t>(kMerVec, of, k, complements, objective, penalty, options);
        WriteLog("Finished pareto optimization of masked superstring, resulting objective value: " + std::to_string(value) + ".");
    }
    catch (const std::exception& e){
        WriteLog("Exception was thrown: " + std::string(e.what()) + ".");
        return false;
    }
    return true;
}
