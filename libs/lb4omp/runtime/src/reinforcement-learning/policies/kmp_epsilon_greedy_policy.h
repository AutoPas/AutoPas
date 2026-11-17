// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


//  INFORMATION:
//  ------------
//  The "Epsilon Greedy" policy has a chance (probability = 1 - ε) to explore a random action
//  from the action space. If this chance isn't met, the agent using this policy uses the greedy
//  action from the learned experience (exploit). After every learning step, ε gets decreased by
//  a constant factor (epsilon decays) to favor exploration at the start of the learning phase and
//  exploitation in the end.
//
//  Environment variables:
//      KMP_RL_EPSILON
//      KMP_RL_EPS_DECAY
//      KMP_RL_EPS_MIN


#pragma once


#include "kmp_base_policy.h"


class EpsilonGreedyPolicy : public BasePolicy {
protected:
    int policy(int episode, int timestep, Agent *agent) override;
};
