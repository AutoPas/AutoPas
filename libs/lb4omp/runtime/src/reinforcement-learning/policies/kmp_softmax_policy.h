//  ************************* Reinforcement Learning Extension ***********************************
//  * June 2022                                                                                  *
//  * Master Thesis                                                                              *
//  * Luc Kury, <luc.kury@unibas.ch>                                                             *
//  * University of Basel, Switzerland                                                           *
//  **********************************************************************************************
//  * INFORMATION:                                                                               *
//  * ------------                                                                               *
//  * The "Softmax" policy assigns each action a probability and selects the next actions        *
//  * according to this distribution.                                                            *
//  *                                                                                            *
//  * Environment variables:                                                                     *
//  *     TAU (Temperature)                                                                      *
//  **********************************************************************************************


#pragma once


#include <vector>
#include <random>
#include "kmp_base_policy.h"


class SoftmaxPolicy : public BasePolicy {
public:
    int policy(int episode, int timestep, Agent *agent) override;

private:
    /*
     * The temperature parameter here might be 1/temperature seen elsewhere.
     * Here, lower temperatures move the highest-weighted output
     * toward a probability of 1.0.
     * And higher temperatures tend to even out all the probabilities,
     * toward 1/<entry count>.
     * temperature's range is between 0 and +Infinity (excluding these
     * two extremes).
    **/
    std::vector<double> Softmax(const std::vector<double> &weights, double temperature);

    /*
     * Selects one index out of a vector of probabilities, "probabilities"
     * The sum of all elements in "probabilities" must be 1.
     **/
    int StochasticSelection(const std::vector<double> &probabilities);
};

/*
 * Rng class encapsulates random number generation
 * of double values uniformly distributed between 0 and 1,
 * in case you need to replace std's <random> with something else.
 **/
struct Rng {
    std::mt19937 engine;
    std::uniform_real_distribution<double> distribution;

    Rng() : distribution(0, 1) {
        std::random_device rd;
        engine.seed(rd());
    }

    double operator()() {
        return distribution(engine);
    }
};
