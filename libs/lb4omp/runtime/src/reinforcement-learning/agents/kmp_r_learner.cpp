// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#include <random>

#include "kmp_agent.h"
#include "kmp_r_learner.h"


// public
RLearner::RLearner(int num_states, int num_actions) :
        Agent(num_states, num_actions, "R-Learner")
{
    pInit->init(&q_table, state_space, action_space);
    set_table(&q_table);
}

// private
void RLearner::update(int next_state, int next_action, int actions, double reward_value)
{
    // TODO@kurluc00: Implement
}
