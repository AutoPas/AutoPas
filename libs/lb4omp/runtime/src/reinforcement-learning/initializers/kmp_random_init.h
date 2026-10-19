// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


#include <random>
#include "kmp_base_init.h"
#include "../agents/defaults.h"


/*
 * This class implements an initializer for the tabular data structures of the reinforcement learning agents.
 *
 * */
class RandomInit : public BaseInit {
public:
    RandomInit(double upper_bound, double lower_bound);

    void init(int **data, int size) override;

    void init(double **data, int size) override;

    void init(double ***data, int num_actions, int num_states) override;

    double get_double() const;

protected:
    double upper{0.0};
    double lower{-4.0};
};
