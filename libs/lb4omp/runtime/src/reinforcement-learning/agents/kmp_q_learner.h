// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#pragma once

#include "kmp_agent.h"


class QLearner : public Agent {
public:
    explicit QLearner(int num_states, int num_actions);

    ~QLearner() = default;

protected:
    double **q_table{nullptr};

    /* Updates the internal values of the agent. */
    void update(int next_state, int next_action, int actions, double reward_value) override;
};
