// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#include <string>
#include <unordered_map>

#pragma once


enum DecayType {
    EXPONENTIAL = 0,
    STEP = 1
};

static std::unordered_map <std::string, DecayType> DecayTable = {
        {"exponential", DecayType::EXPONENTIAL},
        {"step",        DecayType::STEP}
};

static std::unordered_map<int, std::string> DecayLookup = {
        {DecayType::EXPONENTIAL, "exponential"},
        {DecayType::STEP,        "step"}
};
