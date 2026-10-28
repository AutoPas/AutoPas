// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


#include <unordered_map>


enum PolicyType {
    EXPLORE_FIRST = 0,
    EPSILON_GREEDY = 1,
    SOFTMAX = 2
};

static std::unordered_map <std::string, PolicyType> PolicyTable = {
        {"explore-first",  PolicyType::EXPLORE_FIRST},
        {"epsilon-greedy", PolicyType::EPSILON_GREEDY},
        {"softmax",        PolicyType::SOFTMAX}
};

static std::unordered_map<int, std::string> PolicyLookup = {
        {PolicyType::EXPLORE_FIRST,  "explore-first"},
        {PolicyType::EPSILON_GREEDY, "epsilon-greedy"},
        {PolicyType::SOFTMAX,        "softmax"}
};
