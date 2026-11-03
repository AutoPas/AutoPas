// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include <iostream>
#include "kmp_softmax_policy.h"
#include "../agents/kmp_agent.h"

int SoftmaxPolicy::policy(int episode, int timestep, Agent *agent) {
    int curr_state = agent->get_current_state();
    int no_actions = agent->get_action_space();

    std::vector<double> weights((*agent->get_table())[curr_state], (*agent->get_table())[curr_state] + no_actions);
    auto probabilities = Softmax(weights, agent->get_tau());
    auto action = StochasticSelection(probabilities);
//    while ((int) action == 6){
//    action = StochasticSelection(probabilities);
//    }

    return (int) action;
}

std::vector<double> SoftmaxPolicy::Softmax(const std::vector<double> &weights, double temperature) {
    std::vector<double> probabilities;
    double sum = 0;

    for (auto weight: weights) {
        double pr = std::exp(weight / temperature);
        sum += pr;
        probabilities.push_back(pr);
    }

    for (auto &pr: probabilities) {
        pr /= sum;
    }

    return probabilities;
}

int SoftmaxPolicy::StochasticSelection(const std::vector<double> &probabilities) {
    /*
     * The unit interval is divided into sub-intervals, one for each
     * entry in "probabilities".  Each sub-interval's size is proportional
     * to its corresponding probability.

     * You can imagine a roulette wheel divided into differently-sized
     * slots for each entry.  An entry's slot size is proportional to
     * its probability and all the entries' slots combine to fill
     * the entire roulette wheel.

     * The roulette "ball"'s final location on the wheel is determined
     * by generating a (pseudo)random value between 0 and 1.
     * Then a linear search finds the entry whose sub-interval contains
     * this value.  Finally, the selected entry's index is returned.
    **/

    static Rng rng;
    const double point = rng();
    double cur_cutoff = 0;

    for (std::vector<double>::size_type i = 0; i < probabilities.size() - 1; ++i) {
        cur_cutoff += probabilities[i];
        if (point < cur_cutoff) return i;
    }

    return probabilities.size() - 1;
}
