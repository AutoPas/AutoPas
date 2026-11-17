// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


#include <unordered_map>


enum RewardType {
    LOOPTIME = 0,
    LOOPTIME_AVERAGE = 1,
    LOOPTIME_ROLLING_AVERAGE = 2,
    LOOPTIME_INVERSE = 3,
    LOADIMBALANCE = 4,
    ROBUSTNESS = 5,
    STDDEV = 6,
    SKEWNESS = 7,
    KURTOSIS = 8,
    COV = 9
};

static std::unordered_map <std::string, RewardType> RewardTable = {
        {"looptime",                 RewardType::LOOPTIME},
        {"looptime-average",         RewardType::LOOPTIME_AVERAGE},
        {"looptime-rolling-average", RewardType::LOOPTIME_ROLLING_AVERAGE},
        {"looptime-inverse",         RewardType::LOOPTIME_INVERSE},
        {"loadimbalance",            RewardType::LOADIMBALANCE},
        {"robustness",               RewardType::ROBUSTNESS},
        {"stddev",                   RewardType::STDDEV},
        {"skewness",                 RewardType::SKEWNESS},
        {"kurtosis",                 RewardType::KURTOSIS},
        {"cov",                      RewardType::COV}
};

static std::unordered_map<int, std::string> RewardLookup = {
        {RewardType::LOOPTIME,                 "looptime"},
        {RewardType::LOOPTIME_AVERAGE,         "looptime-average"},
        {RewardType::LOOPTIME_ROLLING_AVERAGE, "looptime-rolling-average"},
        {RewardType::LOOPTIME_INVERSE,         "looptime-inverse"},
        {RewardType::LOADIMBALANCE,            "loadimbalance"},
        {RewardType::ROBUSTNESS,               "robustness",},
        {RewardType::STDDEV,                   "stddev"},
        {RewardType::SKEWNESS,                 "skewness"},
        {RewardType::KURTOSIS,                 "kurtosis"},
        {RewardType::COV,                      "cov"}
};
