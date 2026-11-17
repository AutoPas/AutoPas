// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "../kmp_loopdata.h"
#include "kmp_looptime_average_reward.h"
#include "../agents/kmp_agent.h"

double LooptimeAverageReward::reward(LoopData *stats, Agent *agent) {
#if (_RL_DEBUG > 1)
    std::cout << "[LooptimeAverageReward::reward] Getting reward ..." << std::endl;
#endif
    double reward_signal = stats->cTime;
    agent->set_reward_val(reward_signal);
    static int count{0};
    double avg = agent->get_high();
    double reward{0};
    count++;

#if (_RL_DEBUG > 0)
    std::cout << "[LooptimeAverageReward::reward] Count: " << count << ", Average: " << avg << ", cDLS/cChunk index: "
              << agent->get_current_action() << std::endl << std::endl;
#endif
    // Good case
    if (reward_signal <= avg) // We choose <= as comparator to be optimistic (equally fast is good also)
    {
        reward = agent->get_reward_num()[0]; // by default: 0.0
    } else { // Bad case
        reward = agent->get_reward_num()[1]; // by default: -2.0
    }

    avg = ((count - 1) * avg + reward_signal) / count;
    agent->set_high(avg);
    return reward;
};
