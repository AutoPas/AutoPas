// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#include "kmp_optimistic_init.h"


void OptimisticInit::init(int **data, int size) {
    *data = new int[size];

    int i;
    for (i = 0; i < size; i++) {
        (*data)[i] = 0;
    }
}

void OptimisticInit::init(double **data, int size) {
    *data = new double[size];

    int i;
    for (i = 0; i < size; i++) {
        (*data)[i] = 0.0f;
    }
}

void OptimisticInit::init(double ***data, int num_states, int num_actions) {
    *data = new double *[num_states];

    int i;
    for (i = 0; i < num_states; i++) {
        (*data)[i] = new double[num_actions];
    }

    int s, a;
    for (s = 0; s < num_states; s++) {
        for (a = 0; a < num_actions; a++)
            (*data)[s][a] = 0.0f;
    }
}
