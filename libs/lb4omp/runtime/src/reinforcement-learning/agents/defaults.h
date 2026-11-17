// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#include "../rewards/kmp_reward_type.h"
#include "../initializers/kmp_init_type.h"
#include "../policies/kmp_policy_type.h"
#include "../decays/kmp_decay_type.h"

#pragma once


namespace defaults {
    const double SEED                    = 420.69f;

    const double ALPHA                   = 0.85f;
    const double ALPHA_MIN               = 0.10f;
    const double ALPHA_DECAY_FACTOR      = 0.01f;
    const double GAMMA                   = 0.95f;
    const double EPSILON                 = 0.90f;
    const double EPSILON_MIN             = 0.10f;
    const double EPSILON_DECAY_FACTOR    = 0.01f;
    const double TAU                     = 1.50f;

    const int CHUNK_TYPE                 = 8;

    const int ROLLING_AVG_WINDOW         = 10;
    const int INVERSE_REWARD_MULT        = 10;

    const RewardType REWARD_TYPE         = RewardType::LOOPTIME;
    const InitType INIT_TYPE             = InitType::ZERO;
    const PolicyType POLICY_TYPE         = PolicyType::EXPLORE_FIRST;
    const DecayType DECAY_TYPE           = DecayType::EXPONENTIAL;

    const std::string KMP_RL_AGENT_STATS = "notOutputtingStats";

    const std::string REWARD_STRING      = "0.0,-2.0,-4.0"; // Place good reward on the left (beginning) of array and bad reward on the right (end)
}
