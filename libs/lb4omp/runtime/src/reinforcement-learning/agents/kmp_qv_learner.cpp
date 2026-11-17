// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//



#include "kmp_agent.h"
#include "kmp_qv_learner.h"


// public
QVLearner::QVLearner(int num_states, int num_actions) :
        Agent(num_states, num_actions, "QV-Learner")
{
    pInit->init(&v_table, state_space);
    pInit->init(&q_table, state_space, action_space);
    set_table(&q_table);
}

// private
void QVLearner::update(int next_state, int next_action, int actions, double reward_value)
{
    v_table[current_state] += alpha * (reward_value + (gamma * V(next_state) - V(current_state)));
    q_table[current_state][next_action] += alpha * (reward_value + gamma * V(next_state) - Q(current_state, next_action));
}

double QVLearner::V(int state)
{
    return v_table[state];
}
