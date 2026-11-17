// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


#include "../kmp_loopdata.h"
#include "kmp_base_reward.h"


class Agent;


class CovReward : public BaseReward {
public:
    double reward(LoopData *stats, Agent *agent) override;
};
