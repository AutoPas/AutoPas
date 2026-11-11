// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include <numeric>
#include "../kmp_loopdata.h"
#include "kmp_looptime_rolling_average_reward.h"
#include "../agents/kmp_agent.h"

double LooptimeRollingAverageReward::reward(LoopData *stats, Agent *agent) {
#if (_RL_DEBUG > 1)
    std::cout << "[LooptimeRollingAverageReward::reward] Getting reward ..." << std::endl;
#endif
    double reward_signal = stats->cTime;
    static std::vector<double> avgs;
    double reward{0};
    float reward_mult;

    if (avgs.size() >= agent->get_window()) {
        avgs.pop_back();
    }

    avgs.insert(avgs.begin(), stats->cTime);
    double sum = std::accumulate(avgs.begin(), avgs.end(), 0.0);
    double avg = sum / avgs.size();

#if (_RL_DEBUG > 0)
    std::cout << "[LooptimeRollingAverageReward::reward] Size: " << avgs.size() << ", Average: " << avg << ", Method: "
              << agent->get_current_action() << std::endl << std::endl;
#endif
    // Good case
    if (reward_signal <= avg) // We choose <= as comparator to be optimistic (equally fast is good also)
    {
#if (_RL_DEBUG > 1)
        std::cout << "[LooptimeRollingAverageReward::reward] Good!" << std::endl;
#endif
        if (avg == 0 || avg == 999) reward_mult = 1;
        else reward_mult = std::abs((avg - reward_signal) / reward_signal);
        reward = reward_mult * agent->get_reward_num()[0]; // by default: 0.0
    } else
        // Bad case
    {
#if (_RL_DEBUG > 1)
        std::cout << "[LooptimeRollingAverageReward::reward] Bad!" << std::endl;
#endif
        if (avg == 0 || avg == 999) reward_mult = 1;
        else reward_mult = std::abs((reward_signal - avg) / avg);
        reward = reward_mult * agent->get_reward_num()[1]; // by default: -2.0
    }

    agent->set_reward_val(reward_signal);
    agent->set_high(avg);

    return reward;
};
