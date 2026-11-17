// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


#include "kmp_agent.h"


class ExpectedSARSALearner : public Agent {
public:
    explicit ExpectedSARSALearner(int numStates, int numActions);

    ~ExpectedSARSALearner() = default;

private:
    double** q_table{nullptr};

    /* Updates the internal values of the agent. */
    void update(int next_state, int next_action, int actions, double reward_value) override;

};
