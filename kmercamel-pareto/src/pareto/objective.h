#pragma once

typedef uint16_t size_k_max;

enum ParetoObjective {
    RUNS, ZEROS
};

ParetoObjective GetParetoObjective(std::string(objective_string)){
    if (objective_string == "runs")       return ParetoObjective::RUNS;
    else if (objective_string == "zeros") return ParetoObjective::ZEROS;
    else throw std::invalid_argument("Wrong objective provided");
}

size_k_max DEFAULT_PENALTY_RUNS = 7;
size_k_max DEFAULT_PENALTY_ZEROS = 2;
