/**
 * @file MDFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

/**
 * Runs the MD-Flex simulation on a single node.
 * This is the default demonstration of AutoPas.
 */
class MDFlexSingleNode : MDFlexSimulation {
 public:
  MDFlexSingleNode(int argc, char **argv);
  ~MDFlexSingleNode();

  /**
   * Runs the simulation
   */
  void run() override;

 private:
  /**
   * Stores the argument count passed to the constructor for later reuse.
   */
  int _argc;

  /**
   * Stores the arguments passed to the constructor for later reuse.
   */
  char **_argv;
};
