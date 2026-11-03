// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include <random>

#include "kmp_epsilon_greedy_policy.h"
#include "../agents/defaults.h"
#include "../agents/kmp_agent.h"

int EpsilonGreedyPolicy::policy(int episode, int timestep, Agent *agent) {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<double> uniform(0, 1);

    double random_unif = uniform(re);
    int next_action;

    // Switches between exploration and exploitation with the probability of epsilon (or 1-epsilon)
    if (random_unif < agent->get_epsilon()) {
        // Explore (random action)
#if (_RL_DEBUG > 0)
        std::cout << "[EpsilonGreedyPolicy::policy] Exploring!" << std::endl;
#endif
        next_action = agent->sample_action(); // Chooses action (which is equal to the next state)
    } else {
        // Exploit (previous knowledge)
#if (_RL_DEBUG > 0)
        std::cout << "[EpsilonGreedyPolicy::policy] Exploiting!" << std::endl;
#endif
        next_action = argmax(*agent->get_table(), agent->get_current_state(), agent->get_action_space());
    }
    return next_action;
}
