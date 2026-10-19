// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include <iostream>
#include <iomanip>
#include "utils.h"

#include <unordered_map>

// std::unordered_map <std::string, LoopData> autoLoopData;

void read_env_string(const char *var_name, std::string &target) {
    if (std::getenv(var_name) != nullptr) {
        target = std::string(std::getenv(var_name));
    } else {
        std::cout << "[read_env_string] Couldn't read envvar " << var_name << ". Using default: " << target
                  << std::endl;
    }
}

void read_env_int(const char *var_name, int &target) {
    if (std::getenv(var_name) != nullptr) {
        target = std::stoi(std::getenv(var_name));
    } else {
        std::cout << "[read_env_int] Couldn't read envvar " << var_name << ". Using default: " << target << std::endl;
    }
}

void read_env_double(const char *var_name, double &target) {
    if (std::getenv(var_name) != nullptr) {
        target = std::stod(std::getenv(var_name));
    } else {
        std::cout << "[read_env_double] Couldn't read envvar " << var_name << ". Using default: " << target
                  << std::endl;
    }
}

void read_env_enum(const char *var_name, RewardType &target) {
    if (std::getenv(var_name) != nullptr) {
        std::string str = std::string(std::getenv(var_name));
        if (RewardTable.find(str) == RewardTable.end()) {
            std::cout << "[read_env_enum] " << var_name << " cannot use value " << str << ". Possible values: ";
            for (auto const &x: RewardTable) {
                std::cout << x.first << ", ";
            }
            std::cout << " using default: " << RewardLookup.at(int(target)) << std::endl;
        } else {
            target = RewardTable.at(str);
        }
    } else {
        std::cout << "[read_env_enum] Couldn't read envvar " << var_name << ". Using default: "
                  << RewardLookup.at(int(target)) << std::endl;
    }
}

void read_env_enum(const char *var_name, InitType &target) {
    if (std::getenv(var_name) != nullptr) {
        std::string str = std::string(std::getenv(var_name));
        if (InitTable.find(str) == InitTable.end()) {
            std::cout << "[read_env_enum] " << var_name << " cannot use value " << str << ". Possible values: ";
            for (auto const &x: InitTable) {
                std::cout << x.first << ", ";
            }
            std::cout << " using default: " << InitLookup.at(int(target)) << std::endl;
        } else {
            target = InitTable.at(str);
        }
    } else {
        std::cout << "[read_env_enum] Couldn't read envvar " << var_name << ". Using default: "
                  << InitLookup.at(int(target)) << std::endl;
    }
}

void read_env_enum(const char *var_name, PolicyType &target) {
    if (std::getenv(var_name) != nullptr) {
        std::string str = std::string(std::getenv(var_name));
        if (PolicyTable.find(str) == PolicyTable.end()) {
            std::cout << "[read_env_enum] " << var_name << " cannot use value " << str << ". Possible values: ";
            for (auto const &x: PolicyTable) {
                std::cout << x.first << ", ";
            }
            std::cout << " using default: " << PolicyLookup.at(int(target)) << std::endl;
        } else {
            target = PolicyTable.at(str);
        }
    } else {
        std::cout << "[read_env_enum] Couldn't read envvar " << var_name << ". Using default: "
                  << PolicyLookup.at(int(target)) << std::endl;
    }
}

void read_env_enum(const char *var_name, DecayType &target) {
    if (std::getenv(var_name) != nullptr) {
        std::string str = std::string(std::getenv(var_name));
        if (DecayTable.find(str) == DecayTable.end()) {
            std::cout << "[read_env_enum] " << var_name << " cannot use value " << str << ". Possible values: ";
            for (auto const &x: DecayTable) {
                std::cout << x.first << ", ";
            }
            std::cout << " using default: " << DecayLookup.at(int(target)) << std::endl;
        } else {
            target = DecayTable.at(str);
        }
    } else {
        std::cout << "[read_env_enum] Couldn't read envvar " << var_name << ". Using default: "
                  << DecayLookup.at(int(target)) << std::endl;
    }
}

void set_env_var(const char *var_name, int value) {
    setenv(var_name, std::to_string(value).c_str(), 1);
}

double sum(const double *array, int length) {
    double sum = 0.0f;

    for (int i = 0; i < length; i++)
        sum += array[i];

    return sum;
}

/*
 * Searches for the best action (index) for a given state (using the Q-Values).
 * */
int argmax(double **ref_table, int next_state, int size) {
    double actions_avg[size] = {0};
    double best_current = -9999;
    int best_index = 0;

    // Calculate the average reward for an action independent of the state
    for (int i = 0; i < size; i++) {
        double rowSum = 0.0;
        int taken_actions = 0;
        for (int j = 0; j < size; j++) {
            if (ref_table[i][j] != 0.0) {
                rowSum += ref_table[i][j];
                taken_actions++;
            }
        }
        actions_avg[i] = rowSum / taken_actions;
    }
    // select the argmax of the best-rewarded action
    for (int i = 0; i < size; i++) {
//        if (actions_avg[i] > best_current && i != 6) {
        if (actions_avg[i] > best_current) {
            best_current = actions_avg[i];
            best_index = i;
        }
    }
    return best_index;
}

void dumpArray2D(double **ref_table, int row, int col) {
    int i, j;
    for (i = 0; i < row; i++) {
        std::cout << "[ ";
        for (j = 0; j < col; j++) {
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            std::cout << ref_table[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
}

// void set_loop_data(const char *loop_id, LoopData *loop_data) {
//     autoLoopData.insert(std::pair<std::string, LoopData>(loop_id, *loop_data));
// }

// LoopData *get_loop_data(const char *loop_id) {
//     return &autoLoopData.at(loop_id);
// }