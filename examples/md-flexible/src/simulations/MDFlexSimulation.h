/**
 * @file MDFlexSimulation.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "../Simulation.h"

/**
 * Provides base functionality for MD-Flexible simulations.
 * Maybe add unique identifier for each class or an MDFlex factory of some kind
 */
class MDFlexSimulation {
 public:
  MDFlexSimulation(int argc, char **argv);
  ~MDFlexSimulation();

  virtual void run() = 0;

 protected:
  MDFlexConfig *_configuration;
  Simulation *_simulation;
};
