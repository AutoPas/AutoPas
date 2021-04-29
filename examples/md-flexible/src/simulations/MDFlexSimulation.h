/**
 * @file MDFlexSimulation.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "../Simulation.h"

/**
 * Handles minimal initialization requriements for MD-Flexible simulations.
 * Derivce this class to create custom simulations.
 */
class MDFlexSimulation {
 public:
  /**
   * Initializes the simulation according to the arguments passed to the main function.
   * @param argc The number of arguments passed in argv.
   * @param argv The arguments passed to the program.
   */
  MDFlexSimulation(int argc, char **argv);

  ~MDFlexSimulation();

  /**
   * Runs the simulation
   */
  virtual void run() = 0;

 protected:
  /**
   * Stores the configuration used for the simulation.
   * The configuration is defined by the .yaml file passed to the application  with the '--yaml-file' argument.
   */
  MDFlexConfig *_configuration;

  /**
   * The simulation, which will be started using the run function.
   * @todo: Create a simulation interface to allow users to provide their own simulation class
   */
  Simulation *_simulation;
};
