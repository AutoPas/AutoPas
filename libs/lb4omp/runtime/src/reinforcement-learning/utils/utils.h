// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#pragma once

#include <string>
#include "../rewards/kmp_reward_type.h"
#include "../policies/kmp_policy_type.h"
#include "../initializers/kmp_init_type.h"
#include "../decays/kmp_decay_type.h"

#include "../kmp_loopdata.h"

/*
 * Reads the environment variable with the name 'var_name' and parses it as a string.
 * */
void read_env_string(const char *var_name, std::string &target);

void read_env_int(const char *var_name, int &target);
//int read_env_int(const char* var_name, int rdefault);

/*
 * Reads the environment variable with the name 'var_name' and parses it as a double.
 * */
void read_env_double(const char *var_name, double &target);

void read_env_enum(const char *var_name, RewardType &target);

void read_env_enum(const char *var_name, InitType &target);

void read_env_enum(const char *var_name, PolicyType &target);

void read_env_enum(const char *var_name, DecayType &target);

/*
 * Set the environment variable with the name 'var_name' to the integer 'value'.
 * */
void set_env_var(const char *var_name, int value);

/*
 * Returns the sum of values (double) in an array.
 * */
double sum(const double *array, int length);

/*
 * Searches for the best action (index) for a given state (using the Q-Values).
 * */
int argmax(double **ref_table, int next_state, int size);

void dumpArray2D(double **ref_table, int row, int col);

void set_loop_data(const char *loop_id, LoopData *loop_data);

LoopData *get_loop_data(const char *loop_id);
