/**
 * @file MDFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

class MDFlexSingleNode : MDFlexSimulation {
 public:
  MDFlexSingleNode(int argc, char **argv);
  ~MDFlexSingleNode();

  void run() override;

 private:
  int _argc;
  char **_argv;
};
