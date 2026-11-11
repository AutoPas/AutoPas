// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


class BaseInit {
public:
  virtual void init(int** data, int size) = 0;
  virtual void init(double** data, int size) = 0;
  virtual void init(double*** data, int num_actions, int num_states) = 0;
};
