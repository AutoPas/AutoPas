// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "../kmp_loopdata.h"
#include "kmp_loadimbalance_reward.h"
#include "../agents/kmp_agent.h"

double LoadimbalanceReward::reward(LoopData *stats, Agent *agent) {
#if (_RL_DEBUG > 1)
    std::cout << "[LoadimbalanceReward::reward] Getting reward ..." << std::endl;
#endif
    double reward_signal = stats->cLB;
    double low = agent->get_low();
    double high = agent->get_high();
    agent->set_reward_val(reward_signal);

#if (_RL_DEBUG > 0)
    std::cout << "[LoadimbalanceReward::reward] High: " << high << ", Low: " << low << ", Reward: " << reward_signal
              << ", cDLS/cChunk index: " << agent->get_current_action() << std::endl << std::endl;
#endif
    // Good case
    if (reward_signal < 1.05 * low) {
        if (reward_signal < low) agent->set_low(reward_signal);
        return agent->get_reward_num()[0];
    }

    // Bad case
    else if (reward_signal > 0.9 * high) {
        if (reward_signal > high) {
            double new_high = (2 * reward_signal + high) / 3.0;
            agent->set_high(new_high);
            agent->set_high(reward_signal);
        }
        return agent->get_reward_num()[2];
    }

    // Neutral case
    else return agent->get_reward_num()[1];
}
