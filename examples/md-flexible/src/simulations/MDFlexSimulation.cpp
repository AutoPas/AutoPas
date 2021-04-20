/**
 * @file MDFlexSimulation.cpp
 * @author J. Körner
 * @date 07.04.2021
 */
#include "MDFlexSimulation.h"

#include "../parsing/MDFlexParser.h"

MDFlexSimulation::MDFlexSimulation(int argc, char **argv) {
  _simulation = new Simulation();
  _configuration = new MDFlexConfig();

  // parse input and only continue of parsing went without hickups
  if (auto parserExitCode = MDFlexParser::parseInput(argc, argv, *_configuration);
      parserExitCode != MDFlexParser::exitCodes::success) {
    if (parserExitCode == MDFlexParser::exitCodes::parsingError) {
      std::cout << "Error when parsing configuration file." << std::endl;
      exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
  }

  // make sure sim box is big enough
  _configuration->calcSimulationBox();
}

MDFlexSimulation::~MDFlexSimulation() {
  delete _configuration;
  delete _simulation;
}
