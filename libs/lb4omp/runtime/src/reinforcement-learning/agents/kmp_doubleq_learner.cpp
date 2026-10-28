// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include <random>
#include "kmp_doubleq_learner.h"


// public
DoubleQLearner::DoubleQLearner(int num_states, int num_actions) :
                Agent(num_states, num_actions, "DoubleQ-Learner")
{
    pInit->init(&q_table_a, state_space, action_space);
    pInit->init(&q_table_b, state_space, action_space);
    pInit->init(&q_table_sum, state_space, action_space);
    set_table(&q_table_sum);
}

// private
void DoubleQLearner::update(int next_state, int next_action, int actions, double reward_value)
{
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<double> uniform(0, 1);
    double updateA = uniform(re);

    if (updateA > 0.5) // Update table A or B with probability of 0.5
    {
        int bestIndex = argmax(q_table_a, next_state, actions);
        q_table_a[current_state][next_state] += alpha * (reward_value + gamma * q_table_b[next_state][bestIndex] - q_table_a[current_state][next_state]);
    }
    else
    {
        int bestIndex = argmax(q_table_a, next_state, actions);
        q_table_b[current_state][next_state] += alpha * (reward_value + gamma * q_table_a[next_state][bestIndex] - q_table_b[current_state][next_state]);
    }

    q_table_sum[current_state][next_state] = q_table_a[current_state][next_state] + q_table_b[current_state][next_state]; // Update the combined table (sum) to find best action in next step
}
