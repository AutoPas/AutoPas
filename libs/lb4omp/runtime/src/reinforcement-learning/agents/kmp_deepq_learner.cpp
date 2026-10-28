// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "kmp_deepq_learner.h"


// public
DeepQLearner::DeepQLearner(int num_states, int num_actions)
    : Agent(num_states, num_actions, "DQN-Learner") {
  // TODO: Implementation
}

// private
void DeepQLearner::update(int next_state, int next_action,
                          int actions, double reward_value) {
  // TODO@kurluc00: Implement
}
