// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include <random>
#include <iostream>
#include "kmp_random_init.h"


RandomInit::RandomInit(double upper_bound, double lower_bound) : upper{upper_bound}, lower{lower_bound} {
}

double RandomInit::get_double() const {
    std::random_device rd;
    std::default_random_engine re(rd());
    static std::uniform_real_distribution<double> unif(lower, upper);
    double the_num = unif(re);
    return the_num;
}

void RandomInit::init(int **data, int size) {
    *data = new int[size];

    int i;
    for (i = 0; i < size; i++) {
        (*data)[i] = 0;
    }
}

void RandomInit::init(double **data, int size) {
    *data = new double[size];

    int i;
    for (i = 0; i < size; i++) {
        (*data)[i] = get_double();
    }
}

void RandomInit::init(double ***data, int num_states, int num_actions) {
    *data = new double *[num_states];

    int i;
    for (i = 0; i < num_states; i++) {
        (*data)[i] = new double[num_actions];
    }

    int s, a;
    for (s = 0; s < num_states; s++) {
        for (a = 0; a < num_actions; a++) {
            (*data)[s][a] = get_double();
        }
    }
}
