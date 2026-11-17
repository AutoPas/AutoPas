// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


class Agent;


class BasePolicy {
public:
  virtual int policy(int episode, int timestep, Agent* agent) = 0;
};
