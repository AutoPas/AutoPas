/**
 * @file MDFlexSimulation.cpp
 * @author J. Körner
 * @date 07.04.2021
 */
#include "MDFlexSimulation.h"

#include "../parsing/MDFlexParser.h"

MDFlexSimulation::MDFlexSimulation(int argc, char **argv) {
  _simulation = std::shared_ptr<Simulation>(new Simulation());
  _configuration = std::shared_ptr<MDFlexConfig>(new MDFlexConfig());

  // parse input and only continue if parsing went without hiccups
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
