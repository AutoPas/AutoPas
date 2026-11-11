// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "kmp_agent.h"
#include "kmp_expected_sarsa_learner.h"


// public
ExpectedSARSALearner::ExpectedSARSALearner(int numStates, int numActions) :
        Agent(numStates, numActions, "ExpectedSARSA-Learner")
{
    pInit->init(&q_table, state_space, action_space);
    set_table(&q_table);
}

// private
void ExpectedSARSALearner::update(int next_state, int next_action, int actions, double reward_value)
{
    int best_action = argmax(q_table, next_state, actions);
    double expected_return = (1 - epsilon) * Q(next_state, best_action) + (epsilon / action_space) * sum(Q(next_state), action_space);
    q_table[current_state][next_action] += alpha * (reward_value + gamma * expected_return - Q(current_state, next_action));
}
