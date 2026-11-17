// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "../kmp_loopdata.h"
#include "kmp_skewness_reward.h"
#include "../agents/kmp_agent.h"

double SkewnessReward::reward(LoopData *stats, Agent *agent) {
#if (_RL_DEBUG > 1)
    std::cout << "[SkewnessReward::reward] Getting reward ..." << std::endl;
#endif
    double reward_signal = fabs(stats->cSkew);
    double low = agent->get_low();
    double high = agent->get_high();
    agent->set_reward_val(reward_signal);

#if (_RL_DEBUG > 0)
    std::cout << "[SkewnessReward::reward] High: " << high << ", Low: " << low << ", Reward: " << reward_signal
              << ", cDLS/cChunk index: " << agent->get_current_action() << std::endl << std::endl;
#endif
    // Good case
    if ((reward_signal) < low) {
#if (_RL_DEBUG > 1)
        std::cout << "[SkewnessReward::reward] Good!" << std::endl;
#endif
        agent->set_low(reward_signal);
        return agent->get_reward_num()[0]; // by default: 0.0
    }
    // Neutral case
    if ((reward_signal > low) && (reward_signal < high)) {
#if (_RL_DEBUG > 1)
        std::cout << "[SkewnessReward::reward] Neutral." << std::endl;
#endif
        return agent->get_reward_num()[1]; // by default: -2.0
    }
    // Bad case
    if (reward_signal > high) {
#if (_RL_DEBUG > 1)
        std::cout << "[SkewnessReward::reward] Bad!" << std::endl;
#endif
        agent->set_high(reward_signal);
        return agent->get_reward_num()[2]; // by default: -4.0
    }
    return 0.0;
}
