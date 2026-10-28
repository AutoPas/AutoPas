// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "../kmp_loopdata.h"
#include "kmp_looptime_inverse_reward.h"
#include "../agents/kmp_agent.h"

double LooptimeInverseReward::reward(LoopData *stats, Agent *agent) {
    double reward_signal = stats->cTime;
    int multiplier = agent->get_inverse_multiplier();

    return (1 / reward_signal) * multiplier;
};
