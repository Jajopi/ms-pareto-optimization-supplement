#pragma once

typedef uint16_t size_k_max;

enum JointObjective {
    RUNS, ZEROS
};

JointObjective GetJointObjective(std::string(objective_string)){
    if (objective_string == "runs")       return JointObjective::RUNS;
    else if (objective_string == "zeros") return JointObjective::ZEROS;
    else throw std::invalid_argument("Wrong objective provided");
}

size_k_max DEFAULT_PENALTY_RUNS = 7;
size_k_max DEFAULT_PENALTY_ZEROS = 2;
