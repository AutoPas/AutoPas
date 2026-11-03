// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once

#include <string>

#include "kmp_agent.h"


class RLearner : public Agent {
public:
    explicit RLearner(int num_states, int num_actions);

    ~RLearner() = default;

private:
    double** q_table{nullptr};

    /* Updates the internal values of the agent. */
    void update(int next_state, int next_action, int actions, double reward_value) override;
};
