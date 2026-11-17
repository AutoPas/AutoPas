// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "../utils/utils.h"
#include "../agents/kmp_agent.h"
#include "../agents/kmp_q_learner.h"

// public
QLearner::QLearner(int num_states, int num_actions) : Agent(num_states, num_actions, "Q-Learner") {
    pInit->init(&q_table, state_space, action_space);
    set_table(&q_table);
}

// private
void QLearner::update(int next_state, int next_action, int actions, double reward_value) {
    int best_action = argmax(q_table, next_state, actions);
    double right_term = alpha * (reward_value + gamma * Q(next_state, best_action) - Q(current_state, next_action));
#if (_RL_DEBUG > 0) // remember the old value for debugging
    double old = q_table[current_state][next_action];
#endif
    q_table[current_state][next_action] += right_term;


#if (_RL_DEBUG > 0) // Print the update
    std::cout << "[QLearner::update] Q[" << current_state << "][" << next_action << "]: " << old << " (old) + "
              << right_term << " (right) = " << q_table[current_state][next_action] << std::endl;
#endif

#if (_RL_DEBUG > 1) // Print the whole Q-table
    std::cout << "[QLearner::update] Q-Table:" << std::endl;
    for (int i = 0; i < actions; i++) {
        for (int j = 0; j < actions; j++) {
            std::cout << "Q[" << i << "][" << j << "]: " << q_table[i][j] << " ";
        }
        std::cout << std::endl;
    }
#endif
}
