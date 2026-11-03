// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#include "kmp_agent.h"
#include "kmp_sarsa_learner.h"

// public
SARSALearner::SARSALearner(int num_states, int num_actions) : Agent(num_states, num_actions, "SARSA-Learner") {
    pInit->init(&q_table, state_space, action_space);
    set_table(&q_table);
}

// private
void SARSALearner::update(int next_state, int next_action, int actions, double reward_value) {
    int next2_action = this->pPolicy->policy(this->episode, this->timestep, this); // (T+1) action from the policy
    double right_term = alpha * (reward_value + gamma * Q(next_state, next2_action) - Q(current_state, next_action));

#if (_RL_DEBUG > 0) // remember the old value for debugging
    double old = q_table[current_state][next_action];
#endif

    q_table[current_state][next_action] += right_term;

#if (_RL_DEBUG > 0) // Print the update
    std::cout << "[SARSALearner::update] Q[" << current_state << "][" << next_action << "]: " << old << " (old) + " << right_term << " (right) = " << q_table[current_state][next_action] << std::endl;
#endif

#if (_RL_DEBUG > 1) // Print the whole Q-table
    std::cout << "[SARSALearner::update] Q-Table:" << std::endl;
    for (int i = 0; i < actions; i++) {
        for (int j = 0; j < actions; j++) {
            std::cout << "Q[" << i << "][" << j << "]: " << q_table[i][j] << " ";
        }
        std::cout << std::endl;
    }
#endif
}
