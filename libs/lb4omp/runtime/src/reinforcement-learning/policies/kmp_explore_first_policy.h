// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//
//  --------------------------------------------------------------------------------------------//
//
//  INFORMATION:
//  ------------
//  The "Explore First" policy explores all actions for (state_space * action_space) timesteps
//  (pure exploration) and then only exploits the agent's knowledge (pure exploitation).
//
//  Environment variables:
//      (none)


#pragma once

#include "kmp_base_policy.h"


class ExploreFirstPolicy : public BasePolicy {
public:
    int policy(int episode, int timestep, Agent *agent) override;
};
