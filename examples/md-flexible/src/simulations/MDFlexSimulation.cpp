/**
 * @file MDFlexSimulation.cpp
 * @author J. Körner
 * @date 07.04.2021
 */
#include "MDFlexSimulation.h"

#include "../parsing/MDFlexParser.h"

#include <iostream>

MDFlexSimulation::MDFlexSimulation(){
	_configuration = new MDFlexConfig();
	_simulation = new Simulation();
}

MDFlexSimulation::~MDFlexSimulation(){
	delete _configuration;
	delete _simulation;
}

void MDFlexSimulation::initialize(int argc, char** argv){
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

 	// print config to console
 	std::cout << *_configuration;

	initializeAutoPas();
}

void MDFlexSimulation::run(){
 	_simulation->simulate(*_autopas);
}

void MDFlexSimulation::finalize(int argc, char** argv){}

