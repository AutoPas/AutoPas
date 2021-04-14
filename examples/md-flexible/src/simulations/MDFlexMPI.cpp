/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

void MDFlexMPI::run(){
	autopas::AutoPas<Simulation::ParticleType> autopas();	
	_simulation->initialize(*_configuration, autopas);

	_simulation->simulate(autopas)

	_simulation->printStatistics(autopas)
}

