// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//


#pragma once


class Agent;


class BaseReward {
public:
  virtual double reward(LoopData* stats, Agent* agent) = 0;
};
