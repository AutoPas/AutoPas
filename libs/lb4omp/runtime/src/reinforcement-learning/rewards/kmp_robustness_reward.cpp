// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "kmp_robustness_reward.h"
#include "../kmp_loopdata.h"
#include "../agents/kmp_agent.h"


double RobustnessReward::reward(LoopData *stats, Agent *agent) {
    // TODO: andrei: investigate why this is not working
    if (stats->cTime <= agent->get_tau() * agent->get_low()) {
        double reward = agent->get_tau() * agent->get_low() - stats->cTime;
        agent->set_reward_val(reward);
        return reward;
    } else {
        agent->set_reward_val(-1.0f);
        return -1.0f;
    }
}
